"""Validate real premise CCGT packages through TRAILS routing and static LCI.

The reference is an independent sparse solve of the exported constant-technology
matrices. Report every biosphere-flow residual, capital pulse dates, the scoped
wrapper, signed quantities and the legacy calendar behaviour. Restricted pilot
packages and reports remain outside version control.
"""

import argparse
from collections import defaultdict
import gc
import hashlib
import importlib
import json
import math
from pathlib import Path

from datapackage import Package
import numpy as np
from scipy.sparse.linalg import splu

from trails import Trails


def identity(metadata):
    return tuple(metadata[k] for k in ("name", "reference product", "unit", "location"))


def dense(values):
    return np.asarray(values.todense() if hasattr(values, "todense") else values)


def independent_graph_reference(model, reference_a, caller, amount):
    """Prove quantity closure, then independently accumulate per-year static LCI.

    The full database is ill-conditioned. Comparing one SuperLU root solve with
    many UMFPACK frontier solves mixes routing errors and numerical differences.
    Use the same numerical backend for an independently assembled reference,
    and retain the single-root comparison separately. No TRAILS solve or
    biosphere-accumulation helper is used to construct this reference.
    """
    lca_module = importlib.import_module("trails.lca")
    n = reference_a.shape[0]
    expanded = defaultdict(lambda: np.zeros(n))
    demands = defaultdict(lambda: np.zeros(n))
    frontier = np.zeros(n)
    factors, matrices = {}, {}
    for node, data in model.graph.nodes(data=True):
        year, act = node[0], node[-1]
        boundary = float(data.get("frontier_amount") or 0)
        expanded[year][act] += float(data.get("amount") or 0) - boundary
        frontier[act] += boundary * reference_a[act, act]
        if boundary:
            roots = data.get("frontier_roots") or {act: boundary}
            if not math.isclose(
                math.fsum(roots.values()), boundary, rel_tol=1e-10, abs_tol=1e-30
            ):
                raise AssertionError("Root attribution does not sum to frontier demand")
            t = model.scenario_index[str(model._map_year_to_scenario_year(year))]
            for root, value in roots.items():
                demands[(year, root)][act] += value * float(model.A[t, act, act])
    x = sum(expanded.values(), start=np.zeros(n))
    rhs = np.zeros(n)
    rhs[caller] = amount
    residual = reference_a @ x + frontier - rhs
    magnitude = abs(reference_a) @ abs(x) + abs(frontier) + abs(rhs)
    scaled_residual = float(np.max(abs(residual) / np.maximum(magnitude, 1e-30)))
    result = np.zeros(model.B.shape[2])
    backend = "umfpack" if lca_module.UmfpackContext is not None else "superlu"
    for year in sorted(expanded):
        t = model.scenario_index[str(model._map_year_to_scenario_year(year))]
        if t not in matrices:
            a = model.A[t].to_scipy_sparse().T.tocsc()
            a.sort_indices()
            b = model.B[t].to_scipy_sparse().T.tocsc()
            matrices[t] = (a, b)
            if backend == "umfpack":
                from scikits.umfpack import UmfpackContext, UMFPACK_A

                context = UmfpackContext()
                context.numeric(a)
                factors[t] = context
            else:
                factors[t] = splu(a)
        a, b = matrices[t]
        supply = expanded[year].copy()
        for (demand_year, _), demand in demands.items():
            if demand_year != year:
                continue
            if backend == "umfpack":
                supply += factors[t].solve(UMFPACK_A, a, demand)
            else:
                supply += factors[t].solve(demand)
        result += b @ supply
    return result, {
        "backend": backend,
        "operator_max_absolute_residual": float(abs(residual).max()),
        "operator_max_componentwise_relative_residual": scaled_residual,
        "operator_conservation_passed": scaled_residual <= 1e-10,
    }


def validate(
    directory,
    *,
    interpolate_annual=True,
    solver_mode="direct",
    cache=False,
    variant="both",
    lci_years=(2023,),
):
    evidence = json.loads((directory / "export-evidence.json").read_text())
    annual = {r["service_year"]: r for r in evidence["annual"]}
    if not lci_years or set(lci_years) - set(annual):
        raise ValueError("Every requested LCI year must have an annual cohort record")
    audit = evidence["rewrite"]
    q = audit["preserved_capital_amount"]
    reports = {}
    for package_variant, resource in evidence["packages"].items():
        if variant != "both" and package_variant != variant:
            continue
        print(f"Checking {package_variant}", flush=True)
        path = Path(resource["path"])
        if hashlib.sha256(path.read_bytes()).hexdigest() != resource["sha256"]:
            raise ValueError("Package differs from the exporter evidence")
        corrected = package_variant == "corrected"
        with Trails(
            Package(str(path)),
            stock_vintage=corrected,
            cache_interpolation=cache,
            interpolate_annual=interpolate_annual,
            value_dtype=np.float64,
        ) as model:
            indices = model.activity_indices[model.template_labels[0]]
            lookup = {identity(d): i for i, d in indices.items()}
            caller = lookup[
                identity(audit["caller" if corrected else "original_caller"])
            ]
            capital = lookup[
                identity(audit["capital_root" if corrected else "original_root"])
            ]
            constructor_identity = dict(audit["zero_shift_bindings"][0]["supplier"])
            if not corrected:
                suffix = f" [stock context {audit['context_id']}: capital]"
                if not constructor_identity["name"].endswith(suffix):
                    raise ValueError("Unrecognised scoped constructor identity")
                constructor_identity["name"] = constructor_identity[
                    "name"
                ].removesuffix(suffix)
            construction = lookup[identity(constructor_identity)]
            a = model.A[0].to_scipy_sparse().T.tocsc()
            b = model.B[0].to_scipy_sparse().T.tocsc()
            rhs = np.zeros(a.shape[0])
            rhs[caller] = 1.0
            factor = splu(a)
            reference = b @ factor.solve(rhs)
            pulse_checks, lci_checks = [], []
            for year in (2030, 2022, 2023, 2025, 2022):
                for amount in (1.0, -2.0):
                    direct = model.expand_temporal_exchanges(
                        year=year, act_idx=caller, amount=amount
                    )
                    pulses = {
                        y: row[capital]
                        for y, row in direct.items()
                        if capital in row and row[capital] != 0
                    }
                    total = math.fsum(pulses.values())
                    if not math.isclose(total, amount * q, rel_tol=1e-10, abs_tol=0):
                        raise AssertionError("Capital exchange total changed")
                    raw_construction = defaultdict(float)
                    for event, value in pulses.items():
                        inner = model.expand_temporal_exchanges(
                            year=event, act_idx=capital, amount=value
                        )
                        for y, row in inner.items():
                            if construction in row:
                                raw_construction[y] += row[construction]
                    if not math.isclose(
                        math.fsum(raw_construction.values()),
                        total,
                        rel_tol=1e-10,
                        abs_tol=0,
                    ):
                        raise AssertionError(
                            "Market-to-constructor exchange total changed"
                        )
                    if corrected:
                        expected = dict(
                            zip(
                                annual[year]["event_years"],
                                [amount * q * w for w in annual[year]["weights"]],
                            )
                        )
                        if set(pulses) != set(expected) or not all(
                            math.isclose(
                                pulses[y], expected[y], rel_tol=1e-10, abs_tol=0
                            )
                            for y in expected
                        ):
                            raise AssertionError(
                                "Actual annual capital calendar differs from the curated cohorts"
                            )
                        # Internal calls use the event year, not the service
                        # year or nearest matrix anchor. Their profile is zero.
                        for event, value in pulses.items():
                            inner = model.expand_temporal_exchanges(
                                year=event, act_idx=capital, amount=value
                            )
                            constructor_pulses = {
                                y: row[construction]
                                for y, row in inner.items()
                                if construction in row and row[construction] != 0
                            }
                            if set(constructor_pulses) != {event} or not math.isclose(
                                constructor_pulses[event],
                                value,
                                rel_tol=1e-10,
                                abs_tol=0,
                            ):
                                raise AssertionError(
                                    "Capital wrapper duplicated or changed the vintage shift"
                                )
                    pulse_checks.append(
                        {
                            "year": year,
                            "amount": amount,
                            "capital_total": total,
                            "capital_mean_age_years": math.fsum(
                                (year - y) * v for y, v in pulses.items()
                            )
                            / total,
                            "capital_event_years": sorted(pulses),
                            "capital_pulses": pulses,
                            "raw_construction_pulses": dict(raw_construction),
                            "raw_construction_mean_age_years": math.fsum(
                                (year - y) * v for y, v in raw_construction.items()
                            )
                            / total,
                        }
                    )
            print(f"Capital pulses passed: {package_variant}", flush=True)
            # A non-anchor service year exercises annual lookup and interpolation.
            # Positive and negative FUs require separate solves to test signs.
            for year, amount in (
                (year, amount) for year in lci_years for amount in (1.0, -2.0)
            ):
                model.temporal_routing(
                    start_year=year,
                    start_act_idx=caller,
                    amount=amount,
                    max_depth=2,
                    adaptive_relative_score_cutoff=None,
                    min_amount=0,
                    show_progress=False,
                    attribute_to_roots=True,
                )
                capital_years = sorted(
                    {node[0] for node in model.graph if node[-1] == capital}
                )
                construction_years = sorted(
                    {node[0] for node in model.graph if node[-1] == construction}
                )
                if corrected and (
                    capital_years != sorted(annual[year]["event_years"])
                    or construction_years != capital_years
                ):
                    raise AssertionError(
                        "Router changed the physical construction calendar"
                    )
                inventory = model.lci(
                    show_progress=False,
                    solver_mode=solver_mode,
                    inventory_backend=(
                        "factorized" if solver_mode == "direct" else "chunked"
                    ),
                )
                print(
                    f"LCI returned: {package_variant}, {year}, amount={amount}",
                    flush=True,
                )
                # Use the same bounded native reduction used by FaIR. Xarray's
                # generic skip-NaN reduction can broadcast an enormous sparse
                # mask over the activity/flow/year/root tensor.
                flows = np.asarray(inventory.coords["flow"].values, dtype=int)
                reduced = model._inventory_builder.reduce_activity_for_flows(
                    flows, show_progress=False
                )
                actual, single_root_target = (
                    dense(reduced.sum(axis=(1, 2)))[flows],
                    reference[flows] * amount,
                )
                independently_routed, operator = independent_graph_reference(
                    model, a, caller, amount
                )
                target = independently_routed[flows]
                biosphere = model.biosphere_indices[model.template_labels[0]]
                fossil = next(
                    i
                    for i, row in biosphere.items()
                    if row["name"] == "Carbon dioxide, fossil"
                    and row["compartment"] == "air"
                    and row["subcompartment"] == "unspecified"
                )
                capital_co2 = dense(reduced[fossil, :, capital])
                capital_co2_by_year = {
                    int(y): float(v)
                    for y, v in zip(inventory.coords["year"].values, capital_co2)
                    if v != 0
                }
                difference = np.abs(actual - target)
                passed = bool(np.allclose(actual, target, rtol=1e-7, atol=1e-12))
                expanded, frontier = np.zeros(a.shape[0]), np.zeros(a.shape[0])
                for node, data in model.graph.nodes(data=True):
                    index = node[-1]
                    boundary = float(data.get("frontier_amount") or 0)
                    expanded[index] += float(data.get("amount") or 0) - boundary
                    frontier[index] += boundary * a[index, index]
                closure = a @ expanded + frontier - rhs * amount
                static_graph = b @ (expanded + factor.solve(frontier))
                worst = np.argsort(difference / (1e-12 + 1e-7 * np.abs(target)))[-10:][
                    ::-1
                ]
                lci_checks.append(
                    {
                        "year": year,
                        "amount": amount,
                        "static_conservation_passed": passed
                        and operator["operator_conservation_passed"],
                        "flow_accumulation_passed": passed,
                        "operator_conservation": operator,
                        "single_root_superlu_max_flow_difference": float(
                            np.max(abs(actual - single_root_target))
                        ),
                        "graph_operator_max_residual": float(np.abs(closure).max()),
                        "graph_static_max_flow_error": float(
                            np.max(np.abs(static_graph - reference * amount))
                        ),
                        "largest_scaled_errors": [
                            {
                                "flow": int(flows[i]),
                                "actual": float(actual[i]),
                                "expected": float(target[i]),
                            }
                            for i in worst
                        ],
                        "biosphere_flows_checked": len(flows),
                        "capital_fossil_co2_kg_by_year": capital_co2_by_year,
                        "max_absolute_flow_error": float(difference.max()),
                        "l1_error_relative_to_reference": float(
                            difference.sum() / np.abs(target).sum()
                        ),
                        "capital_graph_years": capital_years,
                        "construction_graph_years": construction_years,
                        "inventory_year_bounds": [
                            int(inventory.coords["year"].min()),
                            int(inventory.coords["year"].max()),
                        ],
                    }
                )
            reports[package_variant] = {
                "activity_count": a.shape[0],
                "biosphere_flow_count": b.shape[0],
                "pulse_checks": pulse_checks,
                "lci_checks": lci_checks,
            }
        del model
        gc.collect()
    if not reports:
        raise ValueError(
            "The requested package variant is absent from the export evidence"
        )
    return {
        "scope": "real ecoinvent inventory with constant technology at every anchor",
        "case": evidence["case"],
        "interpolate_annual": interpolate_annual,
        "solver_mode": solver_mode,
        "cache_interpolation": cache,
        "lci_service_years": list(lci_years),
        "reports": reports,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--no-interpolation", action="store_true")
    parser.add_argument("--cache", action="store_true")
    parser.add_argument(
        "--variant", choices=["both", "legacy", "corrected"], default="both"
    )
    parser.add_argument("--solver", choices=["direct", "bw2calc"], default="direct")
    parser.add_argument("--lci-years", type=int, nargs="+", default=[2023])
    args = parser.parse_args()
    report = validate(
        args.directory,
        interpolate_annual=not args.no_interpolation,
        solver_mode=args.solver,
        cache=args.cache,
        variant=args.variant,
        lci_years=args.lci_years,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    if not all(
        check["static_conservation_passed"]
        and check["operator_conservation"]["operator_conservation_passed"]
        for check in report["reports"].get("corrected", {}).get("lci_checks", [])
    ):
        raise AssertionError(
            f"Static conservation check failed; diagnostic report: {args.output}"
        )
    if "corrected" in report["reports"]:
        print(f"Verified real CCGT capital calendars and signed LCI: {args.output}")
    else:
        print(
            f"Recorded legacy CCGT comparison, including conservation diagnostics: {args.output}"
        )


if __name__ == "__main__":
    main()
