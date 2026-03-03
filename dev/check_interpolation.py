from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from datapackage import Package
import sparse

from trails.datapackage import load_matrices_from_package, interpolate_to_annual


@dataclass
class SliceLookup:
    """Sorted key/value lookup for a sparse 2D slice."""

    keys: np.ndarray
    vals: np.ndarray


def _build_lookup(matrix_2d: sparse.COO) -> SliceLookup:
    ncols = int(matrix_2d.shape[1])
    i = matrix_2d.coords[0].astype(np.int64, copy=False)
    j = matrix_2d.coords[1].astype(np.int64, copy=False)
    v = matrix_2d.data.astype(np.float64, copy=False)
    k = i * ncols + j
    order = np.argsort(k)
    return SliceLookup(keys=k[order], vals=v[order])


def _lookup_values(lookup: SliceLookup, keys: np.ndarray) -> np.ndarray:
    pos = np.searchsorted(lookup.keys, keys)
    out = np.zeros(keys.shape[0], dtype=np.float64)
    mask = (pos < lookup.keys.size) & (lookup.keys[pos] == keys)
    out[mask] = lookup.vals[pos[mask]]
    return out


@dataclass
class CheckStats:
    n_checked: int = 0
    max_abs_err: float = 0.0
    max_rel_err: float = 0.0
    max_norm_err: float = 0.0
    worst_year: int | None = None
    worst_key: int | None = None
    worst_expected: float = 0.0
    worst_actual: float = 0.0

    def update(
        self,
        *,
        year: int,
        keys: np.ndarray,
        expected: np.ndarray,
        actual: np.ndarray,
        abs_tol: float,
        rel_tol: float,
    ) -> None:
        if expected.size == 0:
            return
        abs_err = np.abs(actual - expected)
        denom = np.maximum(np.abs(expected), 1e-15)
        rel_err = abs_err / denom
        tol = float(abs_tol) + float(rel_tol) * np.abs(expected)
        norm_err = abs_err / np.maximum(tol, 1e-30)

        local_idx = int(np.argmax(norm_err))
        local_abs = float(abs_err[local_idx])
        local_rel = float(rel_err[local_idx])
        local_norm = float(norm_err[local_idx])

        self.n_checked += int(expected.size)
        self.max_abs_err = max(self.max_abs_err, float(abs_err.max(initial=0.0)))
        self.max_rel_err = max(self.max_rel_err, float(rel_err.max(initial=0.0)))
        if local_norm > self.max_norm_err:
            self.max_norm_err = local_norm
            self.max_abs_err = local_abs
            self.max_rel_err = local_rel
            self.worst_year = int(year)
            self.worst_key = int(keys[local_idx])
            self.worst_expected = float(expected[local_idx])
            self.worst_actual = float(actual[local_idx])


def _decode_key(key: int, ncols: int) -> tuple[int, int]:
    i = int(key // ncols)
    j = int(key % ncols)
    return i, j


def _pick_keys_for_interval(
    k0: np.ndarray,
    k1: np.ndarray,
    sample_size: int,
    rng: np.random.Generator,
) -> np.ndarray:
    union = np.union1d(k0, k1)
    if sample_size <= 0 or union.size <= sample_size:
        return union
    idx = rng.choice(union.size, size=sample_size, replace=False)
    return union[np.sort(idx)]


def _check_matrix_interpolation(
    *,
    matrix_name: str,
    base: sparse.COO,
    annual: sparse.COO,
    base_years: list[int],
    annual_index: dict[str, int],
    sample_size: int,
    seed: int,
    abs_tol: float,
    rel_tol: float,
) -> CheckStats:
    rng = np.random.default_rng(seed)
    stats = CheckStats()

    for interval_idx in range(len(base_years) - 1):
        y0 = int(base_years[interval_idx])
        y1 = int(base_years[interval_idx + 1])
        dt = y1 - y0
        if dt <= 0:
            continue

        m0 = base[interval_idx]
        m1 = base[interval_idx + 1]
        l0 = _build_lookup(m0)
        l1 = _build_lookup(m1)

        keys = _pick_keys_for_interval(l0.keys, l1.keys, sample_size, rng)
        if keys.size == 0:
            continue

        v0 = _lookup_values(l0, keys)
        v1 = _lookup_values(l1, keys)

        for y in range(y0 + 1, y1 + 1):
            w = float(y - y0) / float(dt)
            expected = (1.0 - w) * v0 + w * v1
            lookup_y = _build_lookup(annual[int(annual_index[str(y)])])
            actual = _lookup_values(lookup_y, keys)
            stats.update(
                year=y,
                keys=keys,
                expected=expected,
                actual=actual,
                abs_tol=abs_tol,
                rel_tol=rel_tol,
            )

    if stats.worst_key is not None:
        ncols = int(base.shape[2])
        i, j = _decode_key(stats.worst_key, ncols)
        print(
            f"[{matrix_name}] worst mismatch at year={stats.worst_year}, "
            f"coord=({i}, {j}), expected={stats.worst_expected:.12g}, "
            f"actual={stats.worst_actual:.12g}, abs_err={stats.max_abs_err:.3e}, "
            f"norm_err={stats.max_norm_err:.3f}"
        )
    else:
        print(f"[{matrix_name}] no comparable entries found.")

    print(
        f"[{matrix_name}] checked={stats.n_checked} entries, "
        f"max_abs_err={stats.max_abs_err:.3e}, max_rel_err={stats.max_rel_err:.3e}, "
        f"max_norm_err={stats.max_norm_err:.3f}"
    )
    return stats


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Diagnose annual interpolation quality for a TRAILS datapackage. "
            "Compares annual values to expected linear interpolation from "
            "anchor scenario years."
        )
    )
    parser.add_argument(
        "--dataset",
        type=Path,
        required=True,
        help="Path to datapackage.json",
    )
    parser.add_argument(
        "--sample-size",
        type=int,
        default=100_000,
        help=(
            "Max keys sampled per anchor interval (union of nonzero coords). "
            "Use <=0 to check all keys."
        ),
    )
    parser.add_argument(
        "--abs-tol",
        type=float,
        default=1e-8,
        help="Absolute tolerance threshold for pass/fail.",
    )
    parser.add_argument(
        "--rel-tol",
        type=float,
        default=1e-6,
        help="Relative tolerance threshold for pass/fail.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=0,
        help="Random seed for sampled checks.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    dataset = args.dataset.expanduser().resolve()
    if not dataset.exists():
        print(f"Dataset not found: {dataset}")
        return 2

    print(f"Loading datapackage: {dataset}")
    package = Package(str(dataset))
    A, B, labels, _, _, _ = load_matrices_from_package(package)
    years = sorted(int(lbl) for lbl in labels)
    print(
        f"Loaded anchors: years={years[0]}..{years[-1]} "
        f"({len(years)} slices), A shape={A.shape}, B shape={B.shape}"
    )

    A_i, B_i, annual_labels, annual_index = interpolate_to_annual(A, B, labels)
    annual_years = [int(y) for y in annual_labels]
    print(
        f"Interpolated annual slices: years={annual_years[0]}..{annual_years[-1]} "
        f"({len(annual_years)} slices)"
    )

    # Reorder anchors by year so interval indices map to sorted years.
    order = np.argsort(np.array([int(lbl) for lbl in labels], dtype=int))
    A_sorted = A[order]
    B_sorted = B[order]

    print("\nRunning interpolation checks...")
    stats_a = _check_matrix_interpolation(
        matrix_name="A",
        base=A_sorted,
        annual=A_i,
        base_years=years,
        annual_index=annual_index,
        sample_size=int(args.sample_size),
        seed=int(args.seed),
        abs_tol=float(args.abs_tol),
        rel_tol=float(args.rel_tol),
    )
    stats_b = _check_matrix_interpolation(
        matrix_name="B",
        base=B_sorted,
        annual=B_i,
        base_years=years,
        annual_index=annual_index,
        sample_size=int(args.sample_size),
        seed=int(args.seed) + 1,
        abs_tol=float(args.abs_tol),
        rel_tol=float(args.rel_tol),
    )

    passed = max(stats_a.max_norm_err, stats_b.max_norm_err) <= 1.0

    print("\nSummary")
    print(f"- abs_tol={float(args.abs_tol):.3e}, rel_tol={float(args.rel_tol):.3e}")
    print(f"- max_abs_err={max(stats_a.max_abs_err, stats_b.max_abs_err):.3e}")
    print(f"- max_rel_err={max(stats_a.max_rel_err, stats_b.max_rel_err):.3e}")
    print(f"- max_norm_err={max(stats_a.max_norm_err, stats_b.max_norm_err):.3f}")
    print(f"- result={'PASS' if passed else 'FAIL'}")

    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
