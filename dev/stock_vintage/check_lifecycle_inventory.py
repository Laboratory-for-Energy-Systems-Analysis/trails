"""Check scoped lifecycle rewrites against a complete local Brightway matrix.

No database is written and no licensed inventory is exported. The JSON report
contains identities, rewritten capital coefficients and aggregate error metrics.
The pure producer modules are loaded directly so this check can use the existing
Brightway environment independently of premise's application Python version.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
import sys
import types

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu


def load_rewrite(root):
    namespace = types.ModuleType("stock_pilot")
    namespace.__path__ = [str(root / "premise")]
    sys.modules["stock_pilot"] = namespace
    for name in ("stock_vintage", "stock_lifecycle"):
        spec = importlib.util.spec_from_file_location(
            f"stock_pilot.{name}", root / "premise" / f"{name}.py"
        )
        module = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)
    return module.split_embedded_lifecycle


def identity(record, exchange=False):
    return (
        record["name"],
        record["product"] if exchange else record["reference product"],
        record["unit"],
        record["location"],
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--premise-root", type=Path, required=True)
    parser.add_argument("--project", default="ecoinvent-3.12-cutoff")
    parser.add_argument("--database", default="ecoinvent-3.12-cutoff")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    split = load_rewrite(args.premise_root)
    import bw2data as bd
    import bw2calc as bc

    pilots = {
        "passenger_bev_chassis": [
            "a02324de7a448dd7b3a111f44bbe43f7",
            "95cb62e15c32b1b3ab69f1baf8ab581a",
            "7da6cdfd32c431ed2ba85292f4245601",
            "b0fd1b4c1f6ae4e070b3aab88f515752",
        ],
        "heavy_truck": [
            "14ef57aa246952105aa2e6c5c5e32f36",
            "6ef32683d04ca65bb1345468dc8dcba5",
            "d989e6f9515d9e2ca291436755880c08",
            "feaa6a90af6c4675516a9b263df9d18e",
        ],
    }
    if args.project not in {p.name for p in bd.projects}:
        raise ValueError("An existing Brightway project is required")
    previous = bd.projects.current
    reports = []
    try:
        bd.projects.set_current(args.project)
        if args.database not in bd.databases:
            raise ValueError("The requested database is unavailable")
        for pilot, codes in pilots.items():
            activities = [bd.get_activity((args.database, code)) for code in codes]
            lca = bc.LCA({activities[0]: 1.0})
            lca.lci()
            a, b = lca.technosphere_matrix.tocsc(), lca.biosphere_matrix.tocsc()
            precision_differences = []
            # Read only the four reviewed activities; all other dependencies
            # remain represented by the complete Brightway matrices below.
            datasets = []
            supplier_nodes = {}
            for activity in activities:
                row = {
                    k: activity[k]
                    for k in (
                        "name",
                        "reference product",
                        "unit",
                        "location",
                        "database",
                        "code",
                    )
                }
                row["exchanges"] = []
                seen_pairs = set()
                for exchange in activity.exchanges():
                    e = exchange.as_dict()
                    supplier = exchange.input
                    e.update(name=supplier["name"], unit=supplier["unit"])
                    if e["type"] in ("production", "technosphere"):
                        e.update(
                            product=supplier["reference product"],
                            location=supplier["location"],
                        )
                        supplier_nodes[identity(e, exchange=True)] = supplier
                    pair = (e["type"] == "biosphere", supplier.id)
                    if pair in seen_pairs:
                        raise ValueError(
                            "Aggregate duplicate source matrix coordinates before this check"
                        )
                    seen_pairs.add(pair)
                    column = lca.dicts.activity[activity.id]
                    stored = (
                        float(b[lca.dicts.biosphere[supplier.id], column])
                        if e["type"] == "biosphere"
                        else float(a[lca.dicts.product[supplier.id], column])
                        * (1 if e["type"] == "production" else -1)
                    )
                    relative = abs(stored - float(e["amount"])) / max(
                        abs(float(e["amount"])), 1e-30
                    )
                    if relative > np.finfo(np.float32).eps:
                        raise ValueError(
                            "Stored matrix and source coefficient differ beyond float32 precision"
                        )
                    precision_differences.append(relative)
                    # Prove equality in the complete processed Brightway system,
                    # without mixing its rounded coefficients and raw JSON values.
                    e["amount"] = stored
                    row["exchanges"].append(e)
                datasets.append(row)
            updated, audit = split(
                datasets,
                caller=datasets[0],
                paths=[datasets[1:]],
                context_id=pilot,
                uncertainty_mode="deterministic",
            )
            new = updated[len(datasets) :]
            new_index = {identity(d): i for i, d in enumerate(new)}
            n = a.shape[0]
            if a.shape[1] != n:
                raise ValueError("This check requires a square allocated source matrix")
            ar, ac, av, br, bc_, bv = [], [], [], [], [], []
            for j, dataset in enumerate(new):
                for e in dataset["exchanges"]:
                    if e["type"] == "biosphere":
                        br.append(lca.dicts.biosphere[bd.get_activity(e["input"]).id])
                        bc_.append(j)
                        bv.append(float(e["amount"]))
                    elif e["type"] in ("production", "technosphere"):
                        key = identity(e, exchange=True)
                        row = (
                            n + new_index[key]
                            if key in new_index
                            else lca.dicts.product[supplier_nodes[key].id]
                        )
                        ar.append(row)
                        ac.append(j)
                        av.append(
                            float(e["amount"])
                            * (1 if e["type"] == "production" else -1)
                        )
                    else:
                        raise ValueError(
                            f"Unsupported source exchange type: {e['type']}"
                        )
            extension = sparse.coo_matrix(
                (av, (ar, ac)), shape=(n + len(new), len(new))
            ).tocsc()
            augmented_a = sparse.hstack(
                [sparse.vstack([a, sparse.csc_matrix((len(new), n))]), extension]
            ).tocsc()
            augmented_b = sparse.hstack(
                [b, sparse.coo_matrix((bv, (br, bc_)), shape=(b.shape[0], len(new)))]
            ).tocsc()
            rhs = np.zeros(n + len(new))
            rhs[n + new_index[identity(audit["caller"])]] = 1.0
            # Independently solve the small copied subnetwork, then prove its
            # operator identity against ALL rows of the full source matrices.
            # This avoids different pivots of an ill-conditioned background
            # being mistaken for a modelling difference.
            local_supply = splu(extension[n:, :]).solve(rhs[n:])
            removed_supply = np.zeros(n)
            original_nodes = {
                identity(d): node for d, node in zip(datasets, activities)
            }
            for j, d in enumerate(new):
                if d["stock_vintage_role"] in ("capital", "service"):
                    node = original_nodes[
                        identity(d["stock_vintage_original_identity"])
                    ]
                    removed_supply[lca.dicts.activity[node.id]] += local_supply[j]
            a_identity_residual = (
                a @ removed_supply - extension[:n, :] @ local_supply - lca.demand_array
            )
            b_identity_residual = b @ removed_supply - augmented_b[:, n:] @ local_supply
            if (
                np.max(np.abs(a_identity_residual)) > 1e-12
                or np.max(np.abs(b_identity_residual)) > 1e-12
            ):
                raise AssertionError(
                    "Rewritten subnetwork fails the full matrix operator identity"
                )
            supply = np.concatenate([lca.supply_array - removed_supply, local_supply])
            actual = np.asarray(augmented_b @ supply).ravel()
            expected = np.asarray(lca.inventory.sum(axis=1)).ravel()
            error = np.abs(actual - expected)
            tolerance = 1e-11 + 1e-9 * np.abs(expected)
            if np.any(error > tolerance):
                raise AssertionError(f"{pilot}: rewritten biosphere inventory differs")
            residual = float(np.max(np.abs(augmented_a @ supply - rhs)))
            reports.append(
                {
                    "pilot": pilot,
                    "functional_unit": dict(audit["original_caller"], amount=1.0),
                    "activity_count": n,
                    "added_activity_count": len(new),
                    "biosphere_flow_count": len(expected),
                    "max_absolute_flow_error": float(error.max()),
                    "relative_l1_flow_error": float(
                        error.sum() / max(np.abs(expected).sum(), 1e-20)
                    ),
                    "linear_system_residual": residual,
                    "technosphere_operator_identity_residual": float(
                        np.max(np.abs(a_identity_residual))
                    ),
                    "biosphere_operator_identity_residual": float(
                        np.max(np.abs(b_identity_residual))
                    ),
                    "max_raw_to_processed_coefficient_relative_difference": max(
                        precision_differences
                    ),
                    "coefficient_basis": "processed Brightway matrix; raw-to-processed differences checked against float32 precision",
                    "rewrite": audit,
                }
            )
            print(
                f"{pilot}: {len(expected)} biosphere flows agree; max absolute difference {error.max():.3g}"
            )
    finally:
        bd.projects.set_current(previous)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(
            {
                "project": args.project,
                "database": args.database,
                "validation_scope": "static deterministic LCI equivalence; timing and pilot approval not established",
                "pilots": reports,
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
