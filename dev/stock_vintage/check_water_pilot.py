"""Check the real water rewrite, annual component calendars and signed LCI.

Uses native factorized inventory reduction. Restricted package inputs and
derived inventory time series remain local; no Brightway database is modified.
"""

import argparse
from contextlib import redirect_stdout
import gc
import hashlib
import io
import json
import math
from pathlib import Path

from datapackage import Package
import numpy as np
from scipy.sparse.linalg import splu

from trails import Trails

try:
    from .check_ccgt_pilot import identity, independent_graph_reference
except ImportError:
    from check_ccgt_pilot import identity, independent_graph_reference


def prove_scoped_rewrite(a, b, indices, audit):
    """Prove the full operator equality without comparing sparse pivots."""
    lookup = {identity(d): i for i, d in indices.items()}
    prefix = f" [stock context {audit['context_id']}: "
    scoped = sorted(i for i, d in indices.items() if prefix in d["name"])
    caller = lookup[identity(audit["caller"])]
    original = lookup[identity(audit["original_caller"])]
    demand = np.zeros(len(scoped))
    demand[scoped.index(caller)] = 1
    local = splu(a[scoped, :][:, scoped].tocsc()).solve(demand)
    delta = np.zeros(a.shape[0])
    delta[scoped] = -local
    for i, value in zip(scoped, local):
        record = dict(indices[i])
        if not record["name"].endswith((prefix + "capital]", prefix + "service]")):
            continue
        record["name"] = record["name"].split(prefix)[0]
        delta[lookup[identity(record)]] += value
    rhs = np.zeros(a.shape[0])
    rhs[original] = 1
    rhs[caller] = -1
    residual = a @ delta - rhs
    scale = abs(a) @ abs(delta) + abs(rhs)
    relative = float(np.max(abs(residual) / np.maximum(scale, 1e-30)))
    bio = b @ delta
    bio_scale = abs(b) @ abs(delta)
    bio_relative = float(np.max(abs(bio) / np.maximum(bio_scale, 1e-30)))
    if relative > 1e-10 or bio_relative > 1e-10:
        raise AssertionError("Scoped water rewrite changed the source inventory")
    return {
        "scoped_activities": len(scoped),
        "technosphere_max_relative_residual": relative,
        "biosphere_max_relative_residual": bio_relative,
        "net_caller_production": float(a[caller, caller]),
        "static_equivalence_passed": True,
    }


def stream_inventory_totals(builder, flows, fossil, root_roles, n_years):
    """Bound memory while checking all flows and retaining four dated summaries.

    An activity-reduced sparse cube can still require billions of intermediate
    entries for many roots and retirement years. Consume native bounded windows
    and discard each one after updating the required aggregates.
    """
    n_flows = max(flows) + 1
    totals = np.zeros(n_flows)
    labels = sorted(set(root_roles.values()))
    root_lookup = np.full(builder.n_activities, -1, dtype=int)
    for root, role in root_roles.items():
        root_lookup[root] = labels.index(role)
    summaries = np.zeros(len(labels) * n_years)
    count = 0
    for _, _, batch_flows, years, roots, values in builder.iter_entries_for_flows(
        flows, window_entries=100_000, show_progress=False
    ):
        totals += np.bincount(batch_flows, weights=values, minlength=n_flows)
        selected = (batch_flows == fossil) & (root_lookup[roots] >= 0)
        if np.any(selected):
            keys = root_lookup[roots[selected]] * n_years + years[selected]
            summaries += np.bincount(
                keys, weights=values[selected], minlength=len(summaries)
            )
        count += len(values)
    return totals, dict(zip(labels, summaries.reshape(len(labels), n_years))), count


def check_direct_biosphere_dates(model, inventory, b, lookup, audit, year, amount, net):
    """Check pure biosphere adapters, including lifetime occupation at service."""
    result = []
    years = np.asarray(inventory.coords["year"].values, dtype=int)
    for row in audit.get("lifted_biosphere", []):
        supplier = lookup[identity(row["supplier"])]
        column = b[:, supplier].tocoo()
        if len(column.data) != 1 or not math.isclose(float(column.data[0]), 1.0):
            raise AssertionError("Direct biosphere adapter changed its unit output")
        flow = int(column.row[0])
        actual = np.zeros(len(years))
        for (
            _,
            _,
            _,
            dates,
            roots,
            values,
        ) in model._inventory_builder.iter_entries_for_flows(
            [flow], window_entries=100_000, show_progress=False
        ):
            selected = roots == supplier
            actual += np.bincount(
                dates[selected], weights=values[selected], minlength=len(years)
            )
        expected = np.zeros(len(years))
        if np.count_nonzero(years == year) != 1:
            raise AssertionError(
                "Requested service year is absent from biosphere calendar"
            )
        expected[years == year] = amount * row["amount_per_calling_dataset"] / net
        if not np.allclose(actual, expected, rtol=1e-10, atol=1e-13):
            raise AssertionError(
                "Service-year biosphere burden was dated outside service"
            )
        result.append(
            {
                "supplier": row["supplier"],
                "service_year": year,
                "amount": amount,
                "flow_index": flow,
                "unit": row["supplier"]["unit"],
                "events": {int(y): float(v) for y, v in zip(years, actual) if v != 0},
                "passed": True,
            }
        )
    return result


def validate(
    directory,
    *,
    variant="corrected",
    interpolate=True,
    cache=False,
    lci_years=(2023,),
    solver="direct",
    pulses_only=False,
):
    evidence = json.loads((directory / "export-evidence.json").read_text())
    resource = evidence["packages"][variant]
    path = Path(resource["path"])
    if hashlib.sha256(path.read_bytes()).hexdigest() != resource["sha256"]:
        raise ValueError("Exported package changed")
    corrected = variant == "corrected"
    audit = evidence["rewrite"]
    loading_log = io.StringIO()
    with redirect_stdout(loading_log):
        model = Trails(
            Package(str(path)),
            stock_vintage=corrected,
            cache_interpolation=cache,
            interpolate_annual=interpolate,
            value_dtype=np.float64,
        )
    print(loading_log.getvalue(), end="", flush=True)
    with model:
        indices = model.activity_indices[model.template_labels[0]]
        lookup = {identity(d): i for i, d in indices.items()}
        caller = lookup[identity(audit["caller" if corrected else "original_caller"])]
        capital = lookup[
            identity(audit["capital_root" if corrected else "original_root"])
        ]
        constructor = dict(
            next(
                row["supplier"]
                for row in audit["zero_shift_bindings"]
                if identity(row["caller"]) == identity(audit["capital_root"])
            )
        )
        if not corrected:
            constructor["name"] = constructor["name"].removesuffix(
                f" [stock context {audit['context_id']}: capital]"
            )
        construction = lookup[identity(constructor)]
        a = model.A[0].to_scipy_sparse().T.tocsc()
        b = model.B[0].to_scipy_sparse().T.tocsc()
        net = float(a[caller, caller])
        if not math.isclose(
            net, audit["caller_net_reference_production"], rel_tol=1e-12
        ):
            raise AssertionError("Water distribution losses changed")
        report = {
            "package_sha256": resource["sha256"],
            "variant": variant,
            "case": evidence["case"],
            "interpolate_annual": interpolate,
            "cache_interpolation": cache,
            "cache_loaded": "Loaded interpolated matrices from cache:"
            in loading_log.getvalue(),
            "solver": solver,
            "activity_count": a.shape[0],
            "biosphere_flow_count": b.shape[0],
            "pulse_checks": [],
            "lci_checks": [],
        }
        if corrected:
            report["rewrite_equivalence"] = prove_scoped_rewrite(a, b, indices, audit)
        expected = evidence["expected_exchanges"]
        for year in (2030, 2022, 2023, 2025, 2022):
            for amount in (1.0, -2.0):
                direct = model.expand_temporal_exchanges(
                    year=year, act_idx=caller, amount=amount / net
                )
                rows = expected if corrected else [expected[0]]
                for row in rows:
                    supplier = (
                        lookup[identity(row["supplier"])] if corrected else capital
                    )
                    pulses = {
                        y: r[supplier]
                        for y, r in direct.items()
                        if r.get(supplier, 0) != 0
                    }
                    target = amount * row["amount_per_calling_dataset"] / net
                    if not math.isclose(
                        math.fsum(pulses.values()), target, rel_tol=1e-10, abs_tol=1e-30
                    ):
                        raise AssertionError(
                            "Signed water/component coefficient changed"
                        )
                    if corrected:
                        annual = next(
                            r for r in row["annual"] if r["service_year"] == year
                        )
                        reference = {
                            y: target * w
                            for y, w in zip(annual["event_years"], annual["weights"])
                            if w > 0
                        }
                        if pulses.keys() != reference.keys() or any(
                            not math.isclose(pulses[y], v, rel_tol=1e-10, abs_tol=1e-30)
                            for y, v in reference.items()
                        ):
                            raise AssertionError("Annual water event calendar changed")
                    report["pulse_checks"].append(
                        {
                            "service_year": year,
                            "amount": amount,
                            "role": row["physical_role"],
                            "profile_id": row["profile_id"],
                            "total": math.fsum(pulses.values()),
                            "event_years": sorted(pulses),
                            "pulses": pulses,
                        }
                    )
                    if row["physical_role"] in {
                        "pipe_construction",
                        "parent_construction",
                    }:
                        raw = {}
                        for event, value in pulses.items():
                            nested = model.expand_temporal_exchanges(
                                year=event, act_idx=capital, amount=value
                            )
                            for y, values in nested.items():
                                if values.get(construction, 0) != 0:
                                    raw[y] = raw.get(y, 0) + values[construction]
                        if not math.isclose(
                            math.fsum(raw.values()), target, rel_tol=1e-10
                        ):
                            raise AssertionError(
                                "Raw network construction amount changed"
                            )
                        report["pulse_checks"][-1]["raw_construction_pulses"] = raw
                        report["pulse_checks"][-1]["raw_construction_mean_age"] = (
                            math.fsum((year - y) * value for y, value in raw.items())
                            / target
                        )
        if corrected:
            # Every wrapper is instantaneous over the entire relevant support,
            # including retirement outside inventory anchors. Negative-reference
            # suppliers require signed production conversion, not absolute values.
            profile_support = {
                y for r in report["pulse_checks"] for y in r["event_years"]
            }
            for binding in audit["zero_shift_bindings"]:
                parent, supplier = (
                    lookup[identity(binding["caller"])],
                    lookup[identity(binding["supplier"])],
                )
                for year in sorted(profile_support):
                    result = model.expand_temporal_exchanges(
                        year=year, act_idx=parent, amount=1
                    )
                    pulses = {
                        y: r[supplier]
                        for y, r in result.items()
                        if r.get(supplier, 0) != 0
                    }
                    target = -float(a[supplier, parent]) / float(a[supplier, supplier])
                    if set(pulses) != {year} or not math.isclose(
                        pulses[year], target, rel_tol=1e-10, abs_tol=1e-30
                    ):
                        raise AssertionError(
                            "Wrapper duplicated a component or retirement shift"
                        )
            report["zero_shift_bindings_checked"] = len(audit["zero_shift_bindings"])
        print(f"Water signed pulses and wrappers passed: {variant}", flush=True)
        for year in (() if pulses_only else lci_years):
            for amount in (1.0, -2.0):
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
                if corrected:
                    for row in expected:
                        supplier = lookup[identity(row["supplier"])]
                        record = next(
                            r for r in row["annual"] if r["service_year"] == year
                        )
                        years = {n[0] for n in model.graph if n[-1] == supplier}
                        if years != set(record["event_years"]):
                            raise AssertionError("Router lost physical event dates")
                inventory = model.lci(
                    show_progress=False,
                    solver_mode=solver,
                    inventory_backend="factorized" if solver == "direct" else "chunked",
                )
                print(f"Water LCI returned: {variant}, {year}, {amount}", flush=True)
                if corrected and audit.get("lifted_biosphere"):
                    report.setdefault("direct_biosphere_checks", []).extend(
                        check_direct_biosphere_dates(
                            model, inventory, b, lookup, audit, year, amount, net
                        )
                    )
                flows = np.asarray(inventory.coords["flow"].values, dtype=int)
                biosphere = model.biosphere_indices[model.template_labels[0]]
                fossil = next(
                    i
                    for i, r in biosphere.items()
                    if r["name"] == "Carbon dioxide, fossil"
                    and r["compartment"] == "air"
                    and r["subcompartment"] == "unspecified"
                )
                root_roles = (
                    {
                        lookup[identity(row["supplier"])]: row["physical_role"]
                        for row in expected
                    }
                    if corrected
                    else {capital: "embedded_network_capital"}
                )
                totals, categories, streamed_entries = stream_inventory_totals(
                    model._inventory_builder,
                    flows,
                    fossil,
                    root_roles,
                    len(inventory.coords["year"]),
                )
                actual = totals[flows]
                print(
                    f"Streamed {streamed_entries} inventory entries: {year}, {amount}",
                    flush=True,
                )
                target, operator = independent_graph_reference(model, a, caller, amount)
                error = abs(actual - target[flows])
                passed = bool(np.allclose(actual, target[flows], rtol=1e-7, atol=1e-12))
                time_series = {
                    k: {
                        int(y): float(v)
                        for y, v in zip(inventory.coords["year"].values, series)
                        if v != 0
                    }
                    for k, series in categories.items()
                }
                report["lci_checks"].append(
                    {
                        "service_year": year,
                        "amount": amount,
                        "operator_conservation": operator,
                        "flow_accumulation_passed": passed,
                        "static_conservation_passed": passed
                        and operator["operator_conservation_passed"],
                        "max_absolute_flow_error": float(error.max()),
                        "biosphere_flows_checked": len(flows),
                        "streamed_inventory_entries": streamed_entries,
                        "fossil_co2_kg_by_role_and_year": time_series,
                        "calendar_bounds": [
                            int(inventory.coords["year"].min()),
                            int(inventory.coords["year"].max()),
                        ],
                    }
                )
    gc.collect()
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--variant", choices=["legacy", "corrected"], default="corrected"
    )
    parser.add_argument("--no-interpolation", action="store_true")
    parser.add_argument("--cache", action="store_true")
    parser.add_argument("--pulses-only", action="store_true")
    parser.add_argument("--lci-years", type=int, nargs="+", default=[2023])
    parser.add_argument("--solver", choices=["direct", "bw2calc"], default="direct")
    args = parser.parse_args()
    report = validate(
        args.directory,
        variant=args.variant,
        interpolate=not args.no_interpolation,
        cache=args.cache,
        lci_years=args.lci_years,
        solver=args.solver,
        pulses_only=args.pulses_only,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    if args.variant == "corrected" and not all(
        r["static_conservation_passed"] for r in report["lci_checks"]
    ):
        raise AssertionError(
            f"Water conservation failed; diagnostic report: {args.output}"
        )
    print(f"Recorded water pilot evidence: {args.output}")


if __name__ == "__main__":
    main()
