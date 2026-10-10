"""Archive completed pilot evidence locally, with verified hashes and source code.

This copies restricted inputs and outputs; the destination must be outside a Git
checkout. It creates a new directory and never replaces an existing archive.
Legacy failures are retained as diagnostics, not treated as corrected evidence.
"""

import argparse
from copy import deepcopy
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import xml.etree.ElementTree as ET

REQUIRED = {
    "ccgt-real-npi": [
        ("checks-corrected-audited.json", 2023),
        ("checks-2030-no-interpolation.json", 2030),
    ],
    "pv-real-primary": [("checks-warm-2023.json", 2023), ("checks-2030.json", 2030)],
    "water-real-primary": [
        ("checks-streamed-2023.json", 2023),
        ("checks-streamed-2030.json", 2030),
    ],
    "trucks-real-primary": [
        ("checks-warm-2023.json", 2023),
        ("checks-2030.json", 2030),
    ],
    "bev-real-primary": [("checks-warm-2023.json", 2023), ("checks-2030.json", 2030)],
    "pv-polymer-manufacture": [("checks-2030.json", 2030)],
    "pv-mounting-eol": [("checks-2030.json", 2030)],
}
COHORTS = [
    "ccgt-npi-cohorts.json",
    "ccgt-pk650-cohorts.json",
    "pv-npi2025-cohorts.json",
    "pv-pkbudg650-cohorts.json",
    "water-cohorts.json",
    "trucks-npi2025-cohorts.json",
    "trucks-pkbudg650-cohorts.json",
    "bev-npi2025-cohorts.json",
    "bev-pkbudg650-cohorts.json",
]


def sha256(path):
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def inside(root, path):
    path = path.resolve()
    if not path.is_relative_to(root.resolve()):
        raise ValueError(f"Input escapes the declared source directory: {path}")
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--research-dir", type=Path, required=True)
    parser.add_argument("--public-sources-dir", type=Path, required=True)
    parser.add_argument("--iam-root", type=Path, required=True)
    parser.add_argument("--premise-repo", type=Path, required=True)
    parser.add_argument("--premise-python", type=Path, required=True)
    parser.add_argument("--trails-python", type=Path, required=True)
    parser.add_argument("--destination", type=Path, required=True)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    root, destination = args.research_dir.resolve(), args.destination.resolve()
    if destination.exists():
        raise ValueError("Use a fresh archive destination")
    if any((p / ".git").exists() for p in (destination, *destination.parents)):
        raise ValueError("Restricted evidence must stay outside a Git checkout")
    selected, exports, evidence_checks = {}, {}, []
    test_suites = []
    for name in ("trails-routing.xml", "trails-integration.xml", "premise-pilots.xml"):
        document = ET.parse(root / "validation" / name)
        suites = list(document.iter("testsuite"))
        tests = sum(int(s.get("tests", 0)) for s in suites)
        failed = sum(
            int(s.get("failures", 0)) + int(s.get("errors", 0)) for s in suites
        )
        if not tests or failed:
            raise ValueError(f"Required regression suite incomplete or failed: {name}")
        test_suites.append(
            {"file": name, "tests": tests, "failures_and_errors": failed}
        )

    def register(source, relative, expected=None):
        relative = Path(relative)
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError("Archive paths must be safe relative names")
        actual = sha256(source)
        if expected is not None and actual != expected:
            raise ValueError(f"Input changed from recorded SHA-256: {source}")
        record = {
            "source": str(source.resolve()),
            "path": str(relative),
            "sha256": actual,
            "bytes": source.stat().st_size,
        }
        old = selected.get(str(relative))
        if old is not None and old["sha256"] != actual:
            raise ValueError("Conflicting archive input")
        selected[str(relative)] = record

    for pilot, required in REQUIRED.items():
        directory = root / pilot
        export = json.loads((directory / "export-evidence.json").read_text())
        exports[pilot] = export
        for filename, year in required:
            report = json.loads((directory / filename).read_text())
            corrected = report.get("reports", {}).get("corrected", report)
            rows = corrected.get("lci_checks", [])
            if {(r.get("service_year", r.get("year")), r["amount"]) for r in rows} != {
                (year, 1.0),
                (year, -2.0),
            }:
                raise ValueError(
                    f"Required signed LCI evidence is incomplete: {pilot}/{filename}"
                )
            if not all(
                r["static_conservation_passed"]
                and r["flow_accumulation_passed"]
                and r["operator_conservation"]["operator_conservation_passed"]
                for r in rows
            ):
                raise ValueError(f"Corrected conservation failed: {pilot}/{filename}")
            if (
                report.get("package_sha256", export["packages"]["corrected"]["sha256"])
                != export["packages"]["corrected"]["sha256"]
            ):
                raise ValueError("Validation report refers to a different package")
            evidence_checks.append(
                {
                    "pilot": pilot,
                    "file": filename,
                    "year": year,
                    "signed_full_lci_passed": True,
                }
            )
        for variant, package in export["packages"].items():
            path = inside(root, Path(package["path"]))
            if path.suffix != ".zip":
                raise ValueError("Expected an exported ZIP package")
            register(
                path, Path("pilots") / pilot / variant / path.name, package["sha256"]
            )
        for path in sorted(directory.iterdir()):
            if path.is_file() and path.suffix in {".json", ".log"}:
                name = (
                    "export-evidence.source.json"
                    if path.name == "export-evidence.json"
                    else path.name
                )
                register(path, Path("pilots") / pilot / name)

    inventory = root / "ecoinvent-3.12-cutoff-local.json.gz"
    source_manifest = inventory.with_suffix(inventory.suffix + ".manifest.json")
    source_hash = json.loads(source_manifest.read_text())["sha256"]
    if any(e["inventory_sha256"] != source_hash for e in exports.values()):
        raise ValueError("Pilots refer to different source inventories")
    register(inventory, Path("inputs") / inventory.name, source_hash)
    register(source_manifest, Path("inputs") / source_manifest.name)
    for filename in COHORTS:
        path = root / filename
        register(path, Path("cohorts") / filename)
        report = json.loads(path.read_text())
        for source in report.get("public_sources", report.get("sources", [])):
            name = source.get("local_research_filename", source.get("filename"))
            if not name or Path(name).name != name:
                raise ValueError("Public source requires a safe basename")
            source_path = inside(
                args.public_sources_dir, args.public_sources_dir / name
            )
            register(source_path, Path("public-sources") / name, source["sha256"])
        if "iam_source" in report:
            source = report["iam_source"]
            path = inside(args.iam_root, args.iam_root / source["file"])
            register(path, Path("iam-inputs") / source["file"], source["sha256"])
    cohort_hashes = {
        r["sha256"] for name, r in selected.items() if name.startswith("cohorts/")
    }
    if any(e["cohorts_sha256"] not in cohort_hashes for e in exports.values()):
        raise ValueError("Exported cohort input is missing or has changed")
    for name in (
        "ccgt-sensitivity",
        "pv-sensitivity",
        "water-sensitivity",
        "vehicle-sensitivity",
        "pv-endpoints",
        "parent-calendars",
        "performance",
        "compatibility",
        "baseline-audit",
        "validation",
    ):
        directory = root / name
        for path in sorted(directory.iterdir()):
            if path.is_file() and path.suffix in {
                ".csv",
                ".json",
                ".png",
                ".svg",
                ".pdf",
                ".xml",
            }:
                register(path, Path("figures") / directory.name / path.name)

    for name in (
        "trails-routing-tests.log",
        "trails-integration-tests.log",
        "premise-pilot-tests.log",
    ):
        register(root / name, Path("validation") / name)
    for name in ("synthetic-roundtrip", "lifecycle-roundtrip"):
        for path in sorted((root / name).iterdir()):
            if path.is_file() and path.suffix in {".json", ".zip"}:
                register(path, Path("synthetic-fixtures") / name / path.name)

    repositories = {
        "premise": args.premise_repo.resolve(),
        "trails": Path(__file__).resolve().parents[2],
        "trails_baseline": Path(__file__).resolve().parents[2],
    }
    revisions = {}
    for name, repository in repositories.items():
        status = subprocess.check_output(
            ["git", "-C", str(repository), "status", "--porcelain"], text=True
        )
        if status.strip():
            raise ValueError(f"Commit reviewed changes before archiving {name}")
        revision = (
            "b7d1bb78f0388ba3916270e6654353e01624807e"
            if name == "trails_baseline"
            else "HEAD"
        )
        revisions[name] = subprocess.check_output(
            ["git", "-C", str(repository), "rev-parse", revision], text=True
        ).strip()
    environment_probe = (
        "import json,sys; from importlib.metadata import distributions; "
        "print(json.dumps({'python':sys.version,'executable':sys.executable,"
        "'packages':sorted([{'name':d.metadata['Name'],'version':d.version} "
        "for d in distributions()],key=lambda r:r['name'].lower())}))"
    )
    environments = {
        name: json.loads(
            subprocess.check_output(
                [str(executable), "-c", environment_probe], text=True
            )
        )
        for name, executable in (
            ("premise", args.premise_python),
            ("trails", args.trails_python),
        )
    }
    if args.dry_run:
        print(
            json.dumps(
                {
                    "files": len(selected),
                    "input_bytes": sum(r["bytes"] for r in selected.values()),
                    "revisions": revisions,
                    "checked_full_lci_reports": len(evidence_checks),
                    "destination": str(destination),
                }
            )
        )
        return
    destination.mkdir(parents=True, mode=0o700)
    os.chmod(destination, 0o700)
    records = []
    for record in selected.values():
        source, target = Path(record["source"]), destination / record["path"]
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        os.chmod(target, 0o600)
        if sha256(target) != record["sha256"] or sha256(source) != record["sha256"]:
            raise ValueError(f"Input changed during archive copy: {source}")
        records.append(record)
    for pilot, export in exports.items():
        relocated = deepcopy(export)
        for variant, package in relocated["packages"].items():
            package["path"] = str(
                destination / "pilots" / pilot / variant / Path(package["path"]).name
            )
        if "cohorts_path" in relocated:
            relocated["cohorts_path"] = str(
                destination / "cohorts" / Path(relocated["cohorts_path"]).name
            )
        target = destination / "pilots" / pilot / "export-evidence.json"
        target.write_text(json.dumps(relocated, indent=2, allow_nan=False) + "\n")
        records.append(
            {
                "path": str(target.relative_to(destination)),
                "sha256": sha256(target),
                "bytes": target.stat().st_size,
                "derivation": "Archive-local package/cohort paths; numerical evidence and source copy unchanged",
            }
        )
    (destination / "code").mkdir()
    for name, repository in repositories.items():
        target = destination / "code" / f"{name}-{revisions[name][:12]}.tar.gz"
        subprocess.run(
            [
                "git",
                "-C",
                str(repository),
                "archive",
                "--format=tar.gz",
                f"--output={target}",
                revisions[name],
            ],
            check=True,
        )
        records.append(
            {
                "path": str(target.relative_to(destination)),
                "sha256": sha256(target),
                "bytes": target.stat().st_size,
                "git_revision": revisions[name],
            }
        )
    (destination / "manifest.json").write_text(
        json.dumps(
            {
                "schema_version": 1,
                "rights": "Private restricted ecoinvent/IAM research evidence; do not publish or redistribute",
                "revisions": revisions,
                "environments": environments,
                "files": records,
                "corrected_lci_checks": evidence_checks,
                "regression_suites": test_suites,
                "limitations": "Archive checks conservation evidence, not empirical validity or completion of all methodological release requirements. Legacy failures remain diagnostic records.",
            },
            indent=2,
        )
        + "\n"
    )
    (destination / "README.md").write_text("""# Private stock-vintage pilot evidence

This directory contains restricted ecoinvent and IAM inputs and derived results.
Keep it local. The manifest pins every copied artifact and both code revisions.
`pilots/*/export-evidence.source.json` preserves the original export record;
`export-evidence.json` changes only paths to use packages inside this archive.
Corrected signed full-LCI checks passed before copying. Legacy conservation
failures are retained as diagnostics and are not corrected validation results.

Use the premise code/environment for acquisition, cohort generation and export;
use the TRAILS code/environment for `dev/stock_vintage/check_*_pilot.py` and
sensitivity reporting. The archived Markdown methods describe exact commands,
assumptions and source rights. Save rerun reports under fresh filenames.
This archive is not itself a public release or a claim of native IAM vintages.
""")
    print(
        json.dumps(
            {
                "destination": str(destination),
                "files": len(records),
                "corrected_full_lci_reports": len(evidence_checks),
            }
        )
    )


if __name__ == "__main__":
    main()
