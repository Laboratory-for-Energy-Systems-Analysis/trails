"""Build a local PDF review from the verified stock-vintage evidence archive.

Uses Matplotlib only; no LCI is recomputed. Every consumed archive file is
checked against its manifest. PDF, figures and plotted values are restricted
derived research outputs and belong outside version control (e.g. build/).
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import sys

import numpy as np


EXPECTED_MANIFEST = (
    "c15467bdb209ece5c556248c16edc6257dae034dd121b9a6c935b4f5fe726231"
)
DEFAULT_ARCHIVE = Path.home() / (
    "Library/Application Support/trails/stock-vintage-evidence/"
    "2026-10-10-pilots"
)
ROOT = Path(__file__).resolve().parents[2]
BLUE, TEAL, ORANGE = "#1769A6", "#008876", "#CF6A16"
INK, MUTED, LIGHT = "#183345", "#536977", "#EDF3F6"
GREY, PURPLE = "#929DA5", "#8456A6"
YEAR_COLORS = {2022: BLUE, 2025: TEAL, 2030: ORANGE}
PILOTS = {
    "bev": ("Passenger BEV", "bev-npi2025", "bev-pkbudg650", "bev-real-primary"),
    "truck": (
        "Heavy truck", "trucks-npi2025", "trucks-pkbudg650", "trucks-real-primary"
    ),
    "pv": ("Solar PV", "pv-npi2025", "pv-pkbudg650", "pv-real-primary"),
    "ccgt": ("Gas CCGT", "ccgt-npi", "ccgt-pk650", "ccgt-real-npi"),
    "water": ("Water pipes", "water", None, "water-real-primary"),
}


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class Evidence:
    def __init__(self, root):
        self.root = root
        manifest = root / "manifest.json"
        if sha256(manifest) != EXPECTED_MANIFEST:
            raise ValueError("This report expects the signed-off 10 October archive")
        self.manifest = json.loads(manifest.read_text())
        self.files = {r["path"]: r for r in self.manifest["files"]}
        self.used = {"manifest.json": EXPECTED_MANIFEST}
        self.data = {}
        self.calendars = {}
        self.validation = []
        for key, (_, primary, alternative, directory) in PILOTS.items():
            self.data[key] = [self.read(f"cohorts/{primary}-cohorts.json")]
            if alternative:
                self.data[key].append(
                    self.read(f"cohorts/{alternative}-cohorts.json")
                )
            self.calendars[key] = self.read(
                f"pilots/{directory}/checks-parent-calendar-2030.json"
            )
            export = self.read(f"pilots/{directory}/export-evidence.json")
            source_name = Path(export["cohorts_path"]).name
            if source_name != f"{primary}-cohorts.json" or export["case"] != "primary":
                raise ValueError("Pilot export is not linked to the plotted primary")
            if export["cohorts_sha256"] != self.used[f"cohorts/{source_name}"]:
                raise ValueError("Export and cohort hashes disagree")
        self.check_profiles()

    def read(self, relative):
        path = self.root / relative
        record = self.files[relative]
        digest = sha256(path)
        if digest != record["sha256"] or path.stat().st_size != record["bytes"]:
            raise ValueError(f"Archive file differs: {relative}")
        self.used[relative] = digest
        return json.loads(path.read_text())

    def case(self, key, scenario=0, name="primary", asset="pipes"):
        case = self.data[key][scenario]["cases"][name]
        return case["assets"][asset] if key == "water" else case

    def annual(self, key, year, scenario=0, name="primary", asset="pipes"):
        return next(
            r for r in self.case(key, scenario, name, asset)["annual"]
            if r["service_year"] == year
        )

    def check_profiles(self):
        self.profile_rows = []
        self.summaries = []
        for key, datasets in self.data.items():
            for dataset in datasets:
                scenario = dataset.get("iam_source", {}).get(
                    "scenario", "Constant stock renewal"
                )
                for case_name, outer in dataset["cases"].items():
                    assets = outer["assets"] if key == "water" else {"parent": outer}
                    for asset, case in assets.items():
                        for role, records in case.items():
                            if not isinstance(records, list):
                                continue
                            for record in records:
                                if (
                                    not isinstance(record, dict)
                                    or "weights" not in record
                                ):
                                    continue
                                years = np.array(record["event_years"], dtype=int)
                                weights = np.array(record["weights"], dtype=float)
                                if (
                                    len(years) != len(weights)
                                    or len(set(years)) != len(years)
                                    or not np.all(np.isfinite(weights))
                                    or np.any(weights < 0)
                                    or not math.isclose(
                                        math.fsum(weights), 1, abs_tol=1e-10
                                    )
                                ):
                                    raise ValueError("Invalid annual probability mass")
                                service = record["service_year"]
                                manufacture = role in (
                                    "annual", "battery_manufacture",
                                    "inverter_manufacture",
                                )
                                if manufacture and np.any(years[weights > 0] > service):
                                    raise ValueError("Manufacture after service")
                                if "retirement" in role and np.any(
                                    years[weights > 0] < service
                                ):
                                    raise ValueError("Retirement before service")
                                mean = float(np.dot(years, weights))
                                if "mean_age_years" in record and not math.isclose(
                                    service - mean, record["mean_age_years"],
                                    abs_tol=1e-9,
                                ):
                                    raise ValueError(
                                        "Recorded age disagrees with weights"
                                    )
                                base = dict(
                                    pilot=key, scenario=scenario, case=case_name,
                                    asset=asset, role=role, service_year=service,
                                )
                                self.summaries.append({
                                    **base, "mean_event_year": mean,
                                    "mean_offset_from_service": mean - service,
                                    "mass": math.fsum(weights),
                                })
                                self.profile_rows.extend(
                                    {**base, "event_year": int(y), "weight": float(w)}
                                    for y, w in zip(years, weights)
                                )
        for key, calendar in self.calendars.items():
            weights = dict(zip(
                self.annual(key, 2030)["event_years"],
                self.annual(key, 2030)["weights"],
            ))
            corrected = calendar["reports"]["corrected"]
            pulses = corrected["routed_construction_pulses"]
            total = math.fsum(pulses.values())
            if not math.isclose(
                total, corrected["quantity_per_net_service_unit"], rel_tol=1e-10
            ):
                raise ValueError("Routed parent quantity is not conserved")
            actual = {int(y): v / total for y, v in pulses.items()}
            if any(
                not math.isclose(actual.get(y, 0), weights.get(y, 0), abs_tol=1e-10)
                for y in set(actual) | set(weights)
            ):
                raise ValueError("Plotted cohorts do not match actual routed dates")
        self.validation.append(
            "Every normalized profile checked for finite nonnegative mass, unit sum, "
            "chronology and stored mean age; five primary routed calendars match."
        )


class Report:
    """Small fixed-page layout with measured text wrapping and vector plots."""

    W, H = 595.276, 841.89
    LEFT, RIGHT = 51, 544

    def __init__(self, output):
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
        from matplotlib.font_manager import FontProperties

        self.plt, self.FontProperties = plt, FontProperties
        self.output = output
        (output / "figures").mkdir(parents=True, exist_ok=True)
        self.pdf_path = output / "stock-vintage-pilots-report.pdf"
        self.pdf = PdfPages(self.pdf_path, metadata={
            "Title": "Stock vintages for temporal LCA — five bounded pilots",
            "Author": "Prepared for Romain | premise and TRAILS",
            "Subject": "Methods, generated annual distributions and validation",
            "Keywords": "stock vintage, temporal LCA, premise, TRAILS, cohorts",
        })
        self.page_number = 0
        self.page_titles = []
        self.text_log = []
        self.figure_number = 0

    def start(self, title, section, subtitle=None):
        self.fig = self.plt.figure(figsize=(self.W / 72, self.H / 72), dpi=72)
        self.fig.patch.set_facecolor("white")
        self.page_number += 1
        self.page_titles.append(title)
        self.label(51, 31, "PREMISE  /  TRAILS", size=8, color=TEAL, weight="bold")
        self.label(544, 31, section.upper(), size=8, color=MUTED, ha="right")
        self.rule(51, 47, 544)
        y = self.para(title, 66, size=23, leading=28, weight="bold", color=INK)
        if subtitle:
            y = self.para(subtitle, y + 9, size=10, color=MUTED, leading=14)
        self.y = y + 16
        return self.y

    def label(self, x, y, text, size=10, color=INK, weight="normal", **kwargs):
        self.fig.text(
            x / self.W, 1 - y / self.H, text,
            fontsize=size, color=color, weight=weight,
            va="top", family="DejaVu Sans", **kwargs,
        )
        self.text_log.append(text)

    def rule(self, x1, y, x2, color="#D9E2E7"):
        from matplotlib.lines import Line2D
        self.fig.add_artist(Line2D(
            [x1 / self.W, x2 / self.W], [1 - y / self.H] * 2,
            transform=self.fig.transFigure, color=color, linewidth=0.8,
        ))

    def box(self, x, top, width, height, color=LIGHT):
        from matplotlib.patches import Rectangle
        self.fig.add_artist(Rectangle(
            (x / self.W, 1 - (top + height) / self.H),
            width / self.W, height / self.H,
            transform=self.fig.transFigure, facecolor=color,
            edgecolor="none", zorder=-2,
        ))

    def wrap(self, text, width, size, weight):
        renderer = self.fig.canvas.get_renderer()
        font = self.FontProperties(family="DejaVu Sans", size=size, weight=weight)
        lines = []
        for paragraph in text.split("\n"):
            line = ""
            for word in paragraph.split():
                candidate = f"{line} {word}" if line else word
                measured = renderer.get_text_width_height_descent(
                    candidate, font, ismath=False
                )[0]
                if measured > width and line:
                    lines.append(line)
                    line = word
                else:
                    line = candidate
            lines.append(line)
        return lines

    def para(self, text, y, x=51, width=493, size=10, leading=14.7,
             weight="normal", color=INK):
        lines = self.wrap(text, width, size, weight)
        end = y + leading * len(lines)
        if end > 791:
            raise ValueError(
                f"Page {self.page_number} text ends at {end:.1f}: {text[:60]}"
            )
        self.label(x, y, "\n".join(lines), size=size, color=color,
                   weight=weight, linespacing=leading / size)
        return end

    def heading(self, text, y):
        return self.para(text, y, size=12, leading=17, weight="bold", color=TEAL) + 5

    def note(self, text, y, label=None):
        if label:
            y = self.para(label, y, size=10, weight="bold", leading=15)
        return self.para(text, y, size=9, leading=13.2, color=MUTED) + 8

    def table(self, headers, rows, y, widths, size=9, row_pad=9):
        x0 = self.LEFT
        all_rows = [headers, *rows]
        for i, row in enumerate(all_rows):
            weight = "bold" if i == 0 else "normal"
            wrapped = [self.wrap(str(cell), w - 14, size, weight)
                       for cell, w in zip(row, widths)]
            height = max(map(len, wrapped)) * (size * 1.42) + row_pad * 2
            if y + height > 782:
                raise ValueError(f"Table overflows page {self.page_number}")
            self.box(x0, y, sum(widths), height, LIGHT if i % 2 == 0 else "#FAFCFD")
            x = x0
            for cell, width in zip(wrapped, widths):
                self.label(x + 7, y + row_pad, "\n".join(cell), size=size,
                           weight=weight, linespacing=1.42)
                x += width
            y += height
        return y + 10

    def axes(self, x, top, width, height):
        ax = self.fig.add_axes([
            x / self.W, 1 - (top + height) / self.H,
            width / self.W, height / self.H,
        ])
        ax.spines[["top", "right"]].set_visible(False)
        ax.spines[["left", "bottom"]].set_color("#ADBDC6")
        ax.tick_params(colors=MUTED, labelsize=8)
        ax.grid(axis="y", color="#E4EBEF", linewidth=0.6)
        ax.set_axisbelow(True)
        return ax

    def caption(self, text, y):
        self.figure_number += 1
        return self.note(f"Figure {self.figure_number}. {text}", y)

    def finish(self, slug):
        self.rule(51, 806, 544)
        self.label(51, 816, "10 OCT 2026  •  LOCAL RESEARCH REPORT",
                   size=7, color=MUTED)
        self.label(544, 816, f"{self.page_number:02d}", size=8, color=MUTED, ha="right")
        self.pdf.savefig(self.fig)
        self.fig.savefig(self.output / "figures" / f"{self.page_number:02d}-{slug}.pdf")
        self.plt.close(self.fig)

    def close(self):
        self.pdf.close()


def pmf(ax, record, label, color, linestyle="-", width=1.8):
    years = np.array(record["event_years"], dtype=int)
    lo, hi = years.min(), years.max()
    values = np.zeros(hi - lo + 1)
    values[years - lo] = np.array(record["weights"]) * 100
    ax.stairs(values, np.arange(lo - 0.5, hi + 1.5), label=label,
              color=color, linestyle=linestyle, linewidth=width)


def mean_age(record):
    return record["service_year"] - math.fsum(
        y * w for y, w in zip(record["event_years"], record["weights"])
    )


def get_record(case, role, year=2030):
    return next(r for r in case[role] if r["service_year"] == year)


def style_distribution(ax, ylabel="Share in annual cohort (%)"):
    ax.set_ylabel(ylabel, fontsize=9, color=INK)
    ax.set_xlabel("Manufacture / first-use / commissioning year", fontsize=9, color=INK)
    ax.set_ylim(bottom=0)
    ax.ticklabel_format(axis="x", style="plain", useOffset=False)


def save_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def opening_pages(r, e):
    r.start("Stock vintages\nfor temporal LCA", "Five bounded pilots",
            "Methods, generated distributions and validation • 10 October 2026")
    y = r.para(
        "The opt-in revision replaces generic stock-age curves with explicit annual "
        "cohort profiles. Five bounded pilots satisfy the three release requirements "
        "under declared common-amortisation and proxy assumptions. Original exchange "
        "totals are retained; their attribution to calendar years changes.", r.y + 2,
        size=11, leading=16,
    )
    for x, number, label in [(51, "5", "validated pilots"),
                              (220, "399", "passing tests"),
                              (389, "12", "full signed LCI reports")]:
        r.box(x, y + 16, 155, 70)
        r.label(x + 12, y + 23, number, size=27, color=TEAL, weight="bold")
        r.label(x + 12, y + 60, label, size=8.5, color=MUTED)
    top = y + 134
    ax = r.axes(159, top, 373, 198)
    entries = [(k, "pipes", v[0]) for k, v in PILOTS.items()]
    entries.append(("water", "storage", "Water storage"))
    for i, (key, asset, label) in enumerate(entries):
        a = mean_age(e.annual(key, 2022, asset=asset))
        b = mean_age(e.annual(key, 2030, asset=asset))
        ax.plot([a, b], [i, i], color="#C9D8E2", linewidth=3)
        ax.scatter(a, i, color=BLUE, s=30, zorder=3, label="2022" if i == 0 else None)
        ax.scatter(b, i, color=ORANGE, s=30, zorder=3, label="2030" if i == 0 else None)
        ax.annotate(f"{b:.1f}", (b, i), xytext=(8, 0), textcoords="offset points",
                    va="center", fontsize=8, color=ORANGE)
    ax.set_yticks(range(6), [x[2] for x in entries])
    ax.invert_yaxis()
    ax.set_xlim(0, 43)
    ax.grid(axis="x", color="#E4EBEF")
    ax.set_xlabel("Service-weighted mean parent age (years)", fontsize=9)
    ax.legend(frameon=False, loc="lower right", fontsize=9)
    y = r.caption(
        "Generated mean ages: observed 2022 reference stocks and the 2030 primary "
        "reconstruction. REMIND SSP2-NPi2025 supplies relative vehicle/power paths; "
        "water uses constant-stock renewal. Ages describe timing, not extra burdens.",
        top + 239,
    )
    y = r.heading("How to read this report", y + 8)
    y = r.para(
        "Pages 2–3 explain the method and data. Pages 4–8 show the annual vintage "
        "distributions. Pages 9–13 compare service weights, lifecycle events, legacy "
        "timing and sensitivities. Pages 14–17 cover validation, implementation, "
        "remaining gaps and reproducibility.", y,
    )
    r.note(
        "The PDF contains local IAM- and inventory-derived research results. "
        "The source archive remains private. Public defaults are unchanged; "
        "the feature branches remain unmerged.", y + 12,
    )
    r.finish("overview")

    r.start("What changes methodologically", "Method",
            "Stock composition, service contribution and amortisation "
            "have separate roles.")
    y = r.heading("1. Recover the serving stock", r.y)
    y = r.para(
        "A survival curve describes the probability that an asset remains in use "
        "at a given age. The stock-age distribution also depends on past additions, "
        "retirements and transfers. A rapidly growing BEV fleet therefore need not "
        "have the same age distribution as a mature fleet with the same lifetimes.", y,
    )
    y = r.heading("2. Weight each cohort by the service it provides", y + 10)
    r.label(72, y + 5,
            r"$w_{c,t}=\frac{N_{c,t}\,u_{c,t}}{\sum_j N_{j,t}\,u_{j,t}}$",
            size=20, color=TEAL)
    y = r.para(
        "N is the surviving stock in cohort c and service year t; u is its service "
        "per stock unit. For power, observed generation provides the service weights. "
        "Equal service per vehicle, pipe kilometre or tank is an explicit pilot "
        "assumption, not a property of stock data.", y + 52,
    )
    y = r.heading("3. Distribute the original amortised coefficient", y + 10)
    r.label(72, y + 4,
            r"$a_{c,t}=a_t\,w_{c,t}\qquad\sum_c a_{c,t}=a_t$", size=19, color=TEAL)
    y = r.para(
        "The original signed capital coefficient a remains the amount attributed "
        "to one unit of service. No second division by lifetime and no additional "
        "battery or inverter multiplier is applied. If current use is equal and "
        "the common lifetime-service denominator is retained, stock shares are "
        "the timing weights. Equal calendar lifetimes alone do not imply equal "
        "annual use or equal total lifetime service.", y + 47,
    )
    y = r.heading("4. Evolve observed survivors conditionally", y + 10)
    r.label(64, y + 5,
            r"$N_{c,t}=N_{c,t_0}\,\frac{S_c(t-c)}{S_c(t_0-c)}$", size=19, color=TEAL)
    y = r.para(
        "Observed initial assets have already survived to the reference date. "
        "Applying unconditional survival again would remove them twice. New "
        "cohorts enter separately; additions, retirements and exits reconcile to "
        "the annual stock target. Observed assets outside a finite survival prior "
        "receive a documented residual-life treatment.", y + 52,
    )
    y = r.heading("5. Date components and lifecycle events separately", y + 10)
    y = r.para(
        "A parent asset, its current replacement component, maintenance during "
        "service and eventual disposal have different event dates. Scoped "
        "caller/supplier bindings stop the same vintage delay from being applied "
        "again inside a production chain. Physical dates can precede or follow "
        "the available inventory anchors.", y,
    )
    r.note(
        "Scope: timing under common amortisation. Heterogeneous asset-specific "
        "lifetime-service allocation and physical material-demand forecasting "
        "would require additional methods and data.", y + 12,
    )
    r.finish("method")

    r.start("Data used and what they support", "Evidence",
            "All pilots start from a bounded 2022 reference population "
            "and evolve annually to 2030.")
    y = r.table(
        ["Pilot", "Observed reference", "Future / service basis"],
        [
            ["UK passenger BEV",
             "628,318 cars with known 2010–2022 first use. DfT [S1].",
             "EUR relative BEV stocks and service; equal use per vehicle."],
            ["UK heavy truck", "14,231 articulated diesel vehicles, 32–40 t. DfT [S2].",
             "EUR relative 40-t liquid-fuel truck path; equal use per vehicle."],
            ["US-WECC PV", "363 fixed crystalline plants; "
             "1,294.8 MW AC / 1,629.2 MW DC [S3–S4].",
             "USA relative PV capacity/output; observed output weights."],
            ["US-WECC CCGT", "91 complete non-CHP blocks; 41,476.3 MW AC [S3–S4].",
             "USA CCGT capacity/output; observed output weights."],
            ["Quebec water", "45,114 km pipes and 974 storage assets; "
             "construction bands [S5].",
             "Constant-stock renewal; equal use per km / storage asset."],
        ], r.y, [100, 190, 203], size=8.8,
    )
    y = r.heading("REMIND is the usable scenario bridge", y + 8)
    y = r.para(
        "The screened local REMIND 3.5.2 exports include transport Stock, Sales "
        "and ES variables, plus power Cap, New Cap, generation and technology "
        "Lifetime variables. Two pinned files, SSP2-NPi2025 and SSP2-PkBudg650, "
        "support the vehicle and power reconstructions. Exact native vintage "
        "arrays and exact-run transport survival parameters were unavailable.", y,
    )
    y = r.para(
        "Pinned EDGE/reporting code clarifies annual vehicle sales and the "
        "survival convention. Vehicle support parameters (20 years for cars, "
        "15 for trucks) are not mean lifetimes. Power uses a different survival "
        "formulation; reported additions require explicit period handling.", y + 9,
    )
    y = r.heading("IMAGE and other local exports", y + 14)
    y = r.para(
        "The earlier local audit screened eight IMAGE 3.4 workbooks and found "
        "no physical stock, additions, lifetime or cohort variable families. "
        "Screened GCAM, TIAM-UCL, WITCH and MESSAGE exports likewise did not "
        "supply usable stock-building families. This describes the available "
        "exports, not the full capabilities of those models. No IAM was rerun.", y,
    )
    y = r.heading("The legacy audit supports a scoped revision", y + 14)
    y = r.para(
        "The audit covered 987 stock rows (964 unique identities; 23 "
        "duplicates). All 317 legacy lognormal rows produced uniform "
        "fallback weights. The 670 triangular rows received independent "
        "CDF diagnostics. Original table and runtime semantics remain "
        "unchanged outside the exact opt-in pilot bindings.", y,
        size=9.5, leading=14,
    )
    r.note(
        "Initial observations and scenario constraints remain distinct. Country "
        "stocks combined with larger IAM regions are labelled geographic proxies. "
        "Source identifiers refer to the references on page 17.", y + 12,
    )
    r.finish("data-sources")


PILOT_TEXT = {
    "bev": {
        "title": "Passenger cars: UK battery electric fleet",
        "subtitle": "Observed first-use cohorts • chassis timing • "
                    "equal service per vehicle",
        "intro": "The reference includes 628,318 cars with known first-use dates "
                 "from 2010 onwards (99.8941% of the reported BEV stock). "
                 "The young initial fleet evolves using conditional survival "
                 "and the relative REMIND EUR stock path.",
        "boundary": "The source retains 918.22 kg chassis per 150,000 vehicle-km. "
                    "Battery manufacture/disposal, maintenance and 17 parent "
                    "disposal ports have separate timing. The source already "
                    "contains 1.5 equivalent battery packs; the revision adds "
                    "no replacement multiplier.",
        "limits": "First use proxies manufacture. The primary excludes 309 "
                  "pre-2010 cars and 357 unknown dates. Mileage by age and battery "
                  "histories are unobserved. The primary battery interval is "
                  "13.33 years; 10- and 20-year alternatives are sensitivity cases.",
    },
    "truck": {
        "title": "Heavy trucks: UK articulated diesel fleet",
        "subtitle": "32–40 tonne gross weight • first-use cohorts • "
                    "equal service per vehicle",
        "intro": "The 2022 stock contains 14,231 vehicles. Unknown first-use "
                 "dates for 370 trucks are assigned proportionally; one pre-1980 "
                 "vehicle is spread across a declared 1950–1979 tail. The "
                 "primary future follows relative EUR liquid-fuel 40-t stocks.",
        "boundary": "The source EURO6 lorry coefficient remains 9.65 × 10⁻⁸ "
                    "lorry units per tonne-km. Manufacture, used-lorry disposal "
                    "and service-year maintenance are separated. Manufacturing "
                    "wastewater stays with manufacture.",
        "limits": "The observed fleet has no matching EURO-class breakdown. "
                  "A EURO6 inventory is therefore a technology proxy. "
                  "Age-specific mileage, loading and territorial exports are "
                  "not observed. Reported sales and stock targets do not "
                  "uniquely determine the same additions history.",
    },
    "pv": {
        "title": "Solar PV: US-WECC fixed crystalline plants",
        "subtitle": "Observed commissioning dates • generation-weighted "
                    "service • AC stock balance",
        "intro": "EIA data identify 363 complete plants with 1,294.8 MW AC "
                 "and 1,629.2 MW DC. Observed net output weights the stock "
                 "serving electricity demand; zero-output capacity remains "
                 "in the capacity balance. Future relative paths use USA PV.",
        "boundary": "The source 570 kWp installation coefficient stays "
                    "4.42343 × 10⁻⁸ units per kWh. Installation, panels, active "
                    "inverter, losses and disposal have reviewed roles. "
                    "Land occupation occurs during service; land transformation "
                    "stays at construction.",
        "limits": "Fixed crystalline EIA classification does not establish "
                  "multi-Si or ground mounting. The source is a technology "
                  "proxy. Primary parent life is fixed at 30 years and inverter "
                  "renewal at 15 years. These are timing assumptions; "
                  "replacement histories are unavailable.",
    },
    "ccgt": {
        "title": "Gas CCGT: US-WECC non-CHP blocks",
        "subtitle": "Observed plant cohorts • generation-weighted service • "
                    "USA capacity/output proxy",
        "intro": "The initial fleet comprises 91 complete blocks and "
                 "41,476.3 MW AC. Output-weighted construction cohorts show "
                 "a pronounced early-2000s concentration. A total of 410.5 MW "
                 "lies beyond the survival prior's support and is retained "
                 "with an explicit residual-life extension.",
        "boundary": "The source 400 MW plant coefficient remains "
                    "1.38889 × 10⁻¹¹ plant units per kWh. Plant construction "
                    "uses the serving stock's dates. Independent factory "
                    "capital keeps its own role. No absent disposal inventory "
                    "is invented.",
        "limits": "A quartic continuous-capacity survival analogue uses a "
                  "35-year parameter and a five-year residual-life tail for "
                  "observed overage assets. It is not an exact reproduction "
                  "of REMIND's native discrete vintage logic. Commissioning "
                  "does not identify refurbishment dates.",
    },
}


def pilot_pages(r, e):
    for key, content in PILOT_TEXT.items():
        r.start(content["title"], "Annual distributions", content["subtitle"])
        r.para(content["intro"], r.y)
        ax = r.axes(85, 254, 447, 202)
        for year in (2022, 2025, 2030):
            pmf(ax, e.annual(key, year), f"Service in {year}", YEAR_COLORS[year],
                linestyle="--" if year == 2025 else "-")
        ax.axvspan(2022.5, 2030.5, color="#EDF3F6", zorder=0)
        ax.axvline(2022.5, color="#9BAEB9", linewidth=0.8, linestyle=":")
        all_years = e.annual(key, 2030)["event_years"]
        ax.set_xlim(min(all_years) - 1, 2031)
        style_distribution(ax)
        ax.legend(frameon=False, fontsize=8, loc="upper left")
        r.caption(
            "Explicit annual probability masses from the saved primary "
            "cohort files, without smoothing. Each coloured profile sums to "
            "100%. Grey shading marks post-2022 additions reconstructed "
            "using SSP2-NPi2025. These are shares of service-attributed "
            "parent production timing.", 493,
        )
        r.table(
            ["Service year", "2022", "2025", "2030"],
            [["Mean parent age", *[
                f"{mean_age(e.annual(key, yr)):.2f} years" for yr in (2022, 2025, 2030)
            ]]], 551, [163, 110, 110, 110], size=8.8, row_pad=7,
        )
        y = r.heading("Inventory boundary", 620)
        y = r.para(content["boundary"], y, size=9.4, leading=13.2)
        y = r.heading("Interpretation and limits", y + 10)
        r.para(content["limits"], y, size=9.4, leading=13.2)
        r.finish(f"{key}-distributions")

    r.start("Water infrastructure: Quebec pipes and storage", "Annual distributions",
            "Construction bands annualised explicitly • "
            "constant-stock replacement scenario")
    r.para(
        "Statistics Canada reports 45,114 km of potable-water pipes and 974 "
        "storage assets. Construction dates are unknown for 5.16% of pipe "
        "length and 25.56% of storage counts. The primary distributes known "
        "bands uniformly by year and assigns unknown dates proportionally.", r.y,
    )
    for asset, top, label in [("pipes", 250, "Pipe / network capital"),
                              ("storage", 457, "Storage / tank capital")]:
        ax = r.axes(85, top, 447, 140)
        for year in (2022, 2025, 2030):
            pmf(ax, e.annual("water", year, asset=asset), str(year), YEAR_COLORS[year],
                linestyle="--" if year == 2025 else "-")
        ax.set_title(label, loc="left", fontsize=10, color=INK, pad=10)
        ax.set_xlim(1848, 2032)
        ax.axvspan(2022.5, 2030.5, color=LIGHT, zorder=0)
        style_distribution(ax, "Annual share (%)")
        ax.set_xlabel("Construction-completion year (annualised)", fontsize=9)
        ax.legend(frameon=False, fontsize=8, loc="upper left", ncol=3)
    y = r.caption(
        "Primary pipe and storage manufacture profiles. Each line sums to 100%; "
        "small early tails are retained. Annual precision within historical "
        "survey bands is imputed, not measured. Future additions replace "
        "survival losses while keeping each stock total constant.", 636,
    )
    y = r.para(
        "Mean age changes from 37.12 to 36.78 years for pipes and from 28.01 "
        "to 31.46 years for storage between 2022 and 2030. A Weibull prior "
        "with mean 70 years and shape 3 controls future survival. Equal use "
        "per kilometre / storage asset is assumed; storage counts proxy "
        "tank volume and type.", y + 2, size=9.4, leading=13.2,
    )
    r.note(
        "The original water losses and net production (0.191 / 0.809) are "
        "preserved. Pipe and tank manufacture/disposal remain separate. "
        "No IAM water-demand or network-expansion forecast is claimed.", y + 9,
    )
    r.finish("water-distributions")


def analytical_pages(r, e):
    r.start("Stock shares and service shares", "Weighting",
            "The 2022 power observations illustrate why physical stock "
            "is not always the timing weight.")
    r.para(
        "Installed MW describe the stock. Generated MWh describe its "
        "contribution to electricity service. The primary power profiles "
        "weight construction years by observed generation, retaining "
        "zero-output capacity in the separate stock balance.", r.y,
    )
    stats = []
    for key, top, title in [("pv", 233, "PV: initial reference fleet"),
                            ("ccgt", 466, "CCGT: initial reference fleet")]:
        record = e.annual(key, 2022)
        capacities = record["cohort_capacity_mw_ac"]
        total = math.fsum(capacities.values())
        stock = {"event_years": [int(y) for y in capacities],
                 "weights": [v / total for v in capacities.values()]}
        ax = r.axes(85, top, 447, 153)
        pmf(ax, stock, "Capacity share (MW AC)", GREY, "--")
        pmf(ax, record, "Service share (net MWh)", BLUE)
        ax.set_title(title, loc="left", fontsize=10, color=INK, pad=10)
        style_distribution(ax, "Annual share (%)")
        ax.legend(frameon=False, fontsize=8)
        stock_age = math.fsum((2022 - y) * w
                             for y, w in zip(stock["event_years"], stock["weights"]))
        stats.append((key, stock_age, mean_age(record)))
    y = r.caption(
        "Stock and service distributions normalised separately to 100%. "
        "Both come from the same selected 2022 fleet. Differences reflect "
        "cohort-specific utilisation; they do not represent different "
        "numbers of plants or changed capital quantities.", 662,
    )
    values = "; ".join(
        f"{key.upper()}: {stock:.2f} → {service:.2f} years"
        for key, stock, service in stats
    )
    y = r.para(f"Mean age, capacity-weighted → service-weighted: {values}.",
               y + 2, size=9.5, leading=14)
    r.note(
        "For vehicles and water, equal service per stock unit makes the "
        "primary shares coincide. Declining-use sensitivities expose the "
        "consequence of relaxing that assumption.", y + 9,
    )
    r.finish("stock-versus-service")

    r.start("Components have their own event calendars", "Lifecycle distributions",
            "Service in 2030 • primary SSP2-NPi2025 • "
            "each profile normalised independently")
    r.para(
        "The active battery or inverter can be younger than the asset "
        "containing it. Its retirement is the earlier of replacement and "
        "conditional parent retirement. These calendars distribute the "
        "original role coefficients; they do not predict pack or inverter counts.", r.y,
    )
    panels = [
        ("bev", "annual", "battery_manufacture", "BEV manufacture", 85, 249),
        ("bev", "parent_retirement", "battery_retirement", "BEV disposal", 341, 249),
        ("pv", "annual", "inverter_manufacture", "PV manufacture", 85, 475),
        ("pv", "parent_retirement", "inverter_retirement", "PV disposal", 341, 475),
    ]
    for key, parent, component, title, x, top in panels:
        case = e.case(key)
        ax = r.axes(x, top, 190, 144)
        comp = "Battery" if key == "bev" else "Inverter"
        for role, label, color in [(parent, "Parent", BLUE), (component, comp, TEAL)]:
            record = get_record(case, role)
            pmf(ax, record, label, color, "--" if role == component else "-")
        style_distribution(ax, "Annual share (%)")
        ax.set_xlabel("Event year", fontsize=9)
        ax.set_title(title, loc="left", fontsize=10, color=INK, pad=10)
        ax.legend(frameon=False, fontsize=7.8)
    y = r.caption(
        "Construction and disposal dates from the saved annual profiles. "
        "Each curve has unit mass within its own role; curves must not be "
        "added as physical component counts. The BEV interval is 13.33 years "
        "and the PV inverter interval 15 years in the primary cases.", 663,
    )
    y = r.para(
        "For 2030 service, mean BEV chassis and battery manufacture occur "
        "in 2026.39 and 2026.45; mean disposal occurs in 2043.46 and 2039.99. "
        "PV parent and inverter manufacture occur in 2024.30 and 2025.91; "
        "mean disposal occurs in 2054.30 and 2040.91.", y + 3,
        size=9.5, leading=14,
    )
    r.note("Replacement histories are assumed and tested through sensitivity cases.",
           y + 8)
    r.finish("component-calendars")

    r.start("Conditional retirement extends beyond 2030", "Lifecycle distributions",
            "Disposal attributed to service in 2030 • future event dates "
            "are not clamped to inventory years")
    r.para(
        "A service-year profile can attribute disposal to later calendar "
        "years. These are conditional event distributions for assets serving "
        "in 2030, not a forecast of all future scrap flows. Physical disposal "
        "requires assumptions distinct from territorial or operating-stock exits.", r.y,
    )
    ax = r.axes(85, 240, 447, 139)
    record = get_record(e.case("truck"), "parent_retirement")
    pmf(ax, record, "Truck parent disposal", BLUE)
    ax.set_title("Heavy truck", loc="left", fontsize=10, color=INK, pad=10)
    style_distribution(ax, "Annual share (%)")
    ax.set_xlabel("Disposal year", fontsize=9)
    ax.legend(frameon=False, fontsize=8)
    ax = r.axes(85, 468, 447, 143)
    for asset, label, color in [("pipes", "Pipe / network", BLUE),
                                 ("storage", "Storage / tank", TEAL)]:
        record = get_record(e.case("water", asset=asset), "retirement")
        pmf(ax, record, label, color, "--" if asset == "storage" else "-")
    ax.set_title("Water infrastructure", loc="left", fontsize=10, color=INK, pad=10)
    style_distribution(ax, "Annual share (%)")
    ax.set_xlabel("Disposal year", fontsize=9)
    ax.legend(frameon=False, fontsize=8)
    y = r.caption(
        "Full saved support is shown, including very small terminal tails. "
        "Mean time after 2030 is 6.28 years for trucks, 41.96 for pipes "
        "and 44.89 for storage. The water prior uses Weibull mean 70 years "
        "and shape 3, conditioned on the observed surviving cohorts.", 653,
    )
    y = r.para(
        "Water disposal is assumed to coincide with retirement. Abandonment "
        "in place and demolition delay are unobserved. Vehicle territorial "
        "exits are accounted for separately; an export is not automatically "
        "physical destruction.", y + 3, size=9.5, leading=14,
    )
    r.note(
        "No CCGT disposal curve is shown because the reviewed source lacks "
        "that end-of-life inventory. A missing burden is not created by "
        "assigning a retirement date.", y + 8,
    )
    r.finish("retirement-calendars")

    r.start("Legacy timing and corrected timing", "Direct comparison",
            "Cumulative parent manufacture distributions for one unit of 2030 service")
    r.para(
        "The legacy raw distribution and its actually routed calendar differ "
        "because dates were clamped to inventory support. The corrected "
        "route preserves the physical years. The curves below use the "
        "archived real-export routing evidence, normalised by parent quantity.", r.y,
    )
    variants = [
        ("legacy", "raw_construction_pulses", "Legacy raw", GREY, "--"),
        ("legacy", "routed_construction_pulses", "Legacy routed", ORANGE, ":"),
        ("corrected", "routed_construction_pulses", "Corrected", BLUE, "-"),
    ]
    for i, (key, meta) in enumerate(PILOTS.items()):
        top = 244 + (i // 2) * 174
        ax = r.axes(85 + (i % 2) * 255, top, 191, 109)
        for variant, field, label, color, style in variants:
            pulses = e.calendars[key]["reports"][variant][field]
            xs = sorted(map(int, pulses))
            ws = [pulses[str(x)] for x in xs]
            mass = math.fsum(ws)
            ax.step([xs[0] - 1, *xs, 2031],
                    [0, *list(np.cumsum(ws) * 100 / mass), 100],
                    where="post", color=color, linestyle=style, linewidth=1.7,
                    label=label)
        ax.set_title(meta[0], loc="left", fontsize=10, color=INK, pad=10)
        ax.set_ylim(-2, 105)
        ax.set_xlabel("Manufacture year", fontsize=8)
        ax.set_ylabel("Cumulative share (%)", fontsize=8)
        ax.set_yticks([0, 50, 100])
        ax.ticklabel_format(axis="x", style="plain", useOffset=False)
    ax = r.axes(340, 592, 191, 109)
    ax.axis("off")
    for _, _, label, color, style in variants:
        ax.plot([], [], color=color, linestyle=style, label=label, linewidth=2)
    ax.legend(frameon=False, loc="upper left", fontsize=9)
    r.caption(
        "A steep jump means many attributed events occur in one year. "
        "The selected legacy truck production row has no stock-age profile "
        "and occurs in the service year. Parent-only comparisons exclude "
        "components, maintenance, disposal and background factory capital.", 748,
    )
    r.finish("legacy-calendar-comparison")

    r.start("Sensitivity and scenario dependence", "Robustness",
            "Service-weighted parent mean ages, 2022–2030 • scenario ranges "
            "are not confidence intervals")
    r.para(
        "Solid lines are primary cases. Shaded bands span every documented "
        "sensitivity for that group and scenario. They combine alternative "
        "assumptions about utilisation, survival, missing dates and "
        "stock/additions constraints; they have no probability interpretation.", r.y,
    )
    entries = [(k, "pipes", v[0]) for k, v in PILOTS.items()]
    entries.append(("water", "storage", "Water storage"))
    for i, (key, asset, title) in enumerate(entries):
        top = 247 + (i // 2) * 167
        ax = r.axes(85 + (i % 2) * 255, top, 191, 106)
        for j, dataset in enumerate(e.data[key]):
            years = list(range(2022, 2031))
            profiles = [
                [mean_age(e.annual(key, yr, j, name, asset)) for yr in years]
                for name in dataset["cases"]
            ]
            values = np.array(profiles)
            primary = [mean_age(e.annual(key, yr, j, asset=asset)) for yr in years]
            color = BLUE if j == 0 else ORANGE
            label = "NPi2025" if j == 0 else "PkBudg650"
            if key == "water":
                label = "Constant stock"
            ax.fill_between(years, values.min(axis=0), values.max(axis=0),
                            color=color, alpha=0.14, linewidth=0)
            ax.plot(years, primary, color=color, label=label, linewidth=1.8,
                    linestyle="-" if j == 0 else "--")
        ax.set_title(title, loc="left", fontsize=10, color=INK, pad=10)
        ax.set_ylabel("Mean age (years)", fontsize=8)
        ax.set_xlabel("Service year", fontsize=8)
        ax.set_xticks([2022, 2026, 2030])
        ax.ticklabel_format(axis="x", style="plain", useOffset=False)
        ax.legend(frameon=False, fontsize=7, loc="best")
    r.caption(
        "All saved cases: 14 for BEVs, 12 for trucks, eight for PV, seven "
        "for CCGT (each in two REMIND scenarios), and 13 per water asset. "
        "Intervals span assumptions rather than statistical uncertainty. "
        "Component-only sensitivities can leave the parent mean unchanged.", 746,
    )
    r.finish("scenario-and-sensitivity")


def closing_pages(r, e):
    r.start("Three requirements, five validated pilots", "Validation",
            "Evidence from the actual premise exporter and TRAILS loader/router")
    y = r.table(
        ["Pilot", "1 · Boundaries", "2 · Cohorts", "3 · Integration"],
        [[label, "Pass", "Pass", "Pass"] for label in
         ["Passenger BEV", "Heavy truck", "Solar PV", "Gas CCGT", "Water network"]],
        r.y, [157, 112, 112, 112], size=9, row_pad=5,
    )
    y = r.para(
        "1. Caller/supplier roles, capital/component units and embedded "
        "maintenance/disposal are explicit, preventing repeated vintage shifts.\n"
        "2. Initial cohorts, service weights, conditional survival, future "
        "stock balance and retained amortisation reconcile; proxies are declared.\n"
        "3. Real exports reproduce annual event dates and conserved signed "
        "amounts, with interpolation, cache and lifecycle checks.", y + 4,
        size=9.5, leading=14,
    )
    y = r.heading("Full inventory comparisons", y + 12)
    y = r.para(
        "Twelve corrected full-LCI reports cover all five primary pilots in "
        "2023 with annual interpolation and in 2030 without interpolation, "
        "plus two PV waste-timing endpoints. Both +1 and −2 demands are "
        "checked across all 9,847 biosphere rows against an independently "
        "assembled year-wise calculation and a separate graph-operator identity.", y,
        size=9.5, leading=14,
    )
    y = r.table(
        ["Pilot", "Max |flow difference|*", "Max relative operator residual"],
        [["BEV", "2.85 × 10⁻¹⁴", "1.83 × 10⁻¹⁶"],
         ["Truck", "4.00 × 10⁻¹⁵", "3.04 × 10⁻¹⁶"],
         ["PV + endpoints", "2.67 × 10⁻¹⁴", "2.57 × 10⁻¹⁶"],
         ["CCGT", "2.22 × 10⁻¹⁶", "< 1.2 × 10⁻¹⁶"],
         ["Water", "7.60 × 10⁻¹⁶", "6.45 × 10⁻¹⁶"]],
        y + 12, [120, 175, 198], size=8.7, row_pad=4,
    )
    y = r.note("* Absolute differences are in each biosphere flow's native unit; "
               "these are numerical conservation diagnostics, not impact scores.", y)
    y = r.heading("Regression and independent archive checks", y + 1)
    y = r.para(
        "The signed-off regression suites contain 276 TRAILS and 123 premise "
        "passing tests. The archive has 163 independently verified artifacts. "
        "Relocated smoke checks pass 210 BEV and 340 water signed annual "
        "cases, plus eight full synthetic lifecycle cases.", y,
        size=9.5, leading=14,
    )
    r.note(
        "These are archived implementation results, not freshly rerun LCIs. "
        "Report generation separately verifies consumed source hashes, "
        "profile normalisation/chronology and plotted-vs-routed primary dates.", y + 9,
    )
    r.finish("release-evidence")

    r.start("Implementation and operational limits", "Delivery",
            "Both repositories use feat/stock-vintage-pilots; "
            "the feature remains opt-in.")
    y = r.heading("The producer/consumer contract", r.y)
    y = r.para(
        "Premise binds annual profiles to exact caller and supplier "
        "identities, roles and units. A hashed annual resource carries "
        "consecutive service years and explicit event years. TRAILS uses "
        "stock_vintage=True, checks resource hashes even on warm matrix "
        "caches and rejects missing annual service years. Foreground "
        "overrides take precedence without corrupting the base package cache.", y,
    )
    y = r.para(
        "Physical event year and inventory scenario year are separate. "
        "Manufacture may occur before 2022 and disposal after 2030. "
        "Inventory or adaptive screening-score lookup selects an available "
        "inventory year while the graph retains the physical date.", y + 9,
    )
    y = r.heading("Use the paired implementation revisions", y + 14)
    y = r.para(
        "The preserved old TRAILS reader can load the corrected package "
        "but ignores its annual resource. In the compatibility fixture it "
        "routes 0.05 to 2009 and 0.20 to 2020 instead of the required "
        "0.10 to 2010 and 0.15 to 2021. Unknown metadata alone does not "
        "protect older readers from using incompatible timing.", y,
    )
    y = r.para(
        "The archive pins executable code at premise da7edd7ab87e and "
        "TRAILS 62eacddb764d. Later sign-off commits only change Markdown "
        "documentation. The report script does not modify runtime behaviour.", y + 9,
    )
    y = r.heading("Measured loading and routing", y + 14)
    y = r.table(
        ["Measure", "Corrected packages on this machine"],
        [["Cold matrix-cache load", "5.91–6.54 seconds"],
         ["Confirmed warm-cache load", "1.12–1.29 seconds"],
         ["Depth-two routing", "0.60–1.12 seconds"],
         ["Peak RSS through load/route",
          "Cold: 2,089–2,246 MiB; warm: 1,064–1,160 MiB"],
         ["Package / annual resource",
          "66.3–67.4 MB ZIP; 33 kB–7.4 MB annual resource"]],
        y, [183, 310], size=8.8, row_pad=7,
    )
    y = r.para(
        "Twenty subprocess cases cover legacy/corrected and cold/warm "
        "loading/routing. Full signed verification took minutes to about "
        "152 minutes for water and includes the independent calculation "
        "and concurrent load. Those times are not isolated native LCI "
        "solve benchmarks. Water streams roughly 9.6 billion entries "
        "per signed demand.", y + 4, size=9.5, leading=14,
    )
    r.note(
        "The full 2023 water LCI loaded cold. A separate confirmed warm-cache "
        "pulse check supplies its cache evidence. Full-LCI peak memory was "
        "not instrumented; further performance work is needed for global use.", y + 9,
    )
    r.finish("implementation-and-performance")

    r.start("What remains uncertain", "Scope and next steps",
            "No remaining blocker for these bounded pilots "
            "under their stated assumptions.")
    y = r.table(
        ["Gap", "Consequence and next data priority"],
        [["Exact native IAM vintages",
          "Reconstructions cannot claim exact model reproduction. "
          "Obtain exact-run cohort arrays, reporting versions "
          "and survival parameters."],
         ["Age-specific vehicle service",
          "Equal mileage/loading is assumed. Link cohorts to "
          "mileage, occupancy and payload before broader fleet attribution."],
         ["Component histories", "Battery/inverter renewal and some PV waste dates are "
          "assumptions. Seek installation, replacement and retirement records."],
         ["Water storage detail",
          "25.56% unknown dates; tank counts proxy volume/type. "
          "Prioritise capacity, material, refurbishment and abandonment data."],
         ["Geographic / technology proxies",
          "UK → EUR and WECC → USA paths are bounded proxies. "
          "Truck EURO class and PV multi-Si/ground mounting remain unmatched."],
         ["Heterogeneous amortisation",
          "Individual lifetime-service denominators and Monte "
          "Carlo correlation preservation are outside the validated claim."]],
        r.y, [157, 336], size=9, row_pad=9,
    )
    y = r.heading("Broader release still needs evidence", y + 8)
    y = r.para(
        "The 23-group plan remains the expansion roadmap. Buildings, "
        "bridges and some buried networks have plausible survey/register "
        "routes. Grid equipment, rail components, harbours and industrial "
        "equipment are more dependent on operator records and explicit "
        "component priors. These groups are not calibrated by the five pilots.", y,
    )
    y = r.para(
        "Before changing defaults: extend exact exchange-role audits, "
        "validate geography/technology mappings, resolve compatibility "
        "behaviour, benchmark larger runs and reassess study-level results. "
        "Historical inventory technology also needs its own reconstruction.", y + 9,
    )
    y = r.heading("Interpret the comparisons within their design", y + 14)
    r.para(
        "All real packages use the same ecoinvent 3.12 cut-off technology "
        "at the 2022, 2025 and 2030 anchors to isolate timing and conservation. "
        "They do not demonstrate a fully IAM-transformed background, "
        "historical manufacturing technology, GWP or dynamic climate effects. "
        "The sensitivities are alternative assumptions, not confidence bounds.", y,
        size=9.5, leading=14,
    )
    r.finish("gaps-and-next-steps")

    r.start("Sources and reproducibility", "References",
            "This PDF summarises the signed-off local evidence; "
            "no model or LCI was rerun to draw the figures.")
    references = [
        ("[S1] UK Department for Transport. VEH1111 first-use tables.",
         "https://assets.publishing.service.gov.uk/media/"
         "69ef3555606c20d41216341a/veh1111.ods"),
        ("[S2] UK Department for Transport. VEH0520 goods-vehicle cohorts.",
         "https://assets.publishing.service.gov.uk/media/"
         "69ef3f04ed93f72cf816340e/df_VEH0520.csv"),
        ("[S3] US EIA. Form EIA-860, 2022 generator/plant records.",
         "https://www.eia.gov/electricity/data/eia860/archive/xls/eia8602022.zip"),
        ("[S4] US EIA. Form EIA-923, 2022 net electricity generation.",
         "https://www.eia.gov/electricity/data/eia923/archive/xls/f923_2022.zip"),
        ("[S5] Statistics Canada. Table 34-10-0289, infrastructure construction age.",
         "https://www150.statcan.gc.ca/n1/tbl/csv/34100289-eng.zip"),
        ("[S6] EDGE-Transport. Pinned fleet/depreciation source, 1124c9967e96.",
         "https://github.com/pik-piam/edgeTransport/tree/"
         "1124c9967e9644b41977c4750640629d51b5494b"),
        ("[S7] reporttransport. Pinned reporting source, eb4046a859b2.",
         "https://github.com/pik-piam/reporttransport/tree/"
         "eb4046a859b25363eb8460ced06561010cb4653d"),
    ]
    y = r.y
    for label, url in references:
        r.label(51, y, label, size=8.7, color=BLUE, url=url)
        y += 23
    y = r.note(
        "Clickable source labels refer to the files/code used in the local "
        "campaign (retrieved 9–10 October 2026). Full source hashes, units, "
        "selection audits and rights are recorded in the cohort JSON files "
        "and premise acquisition manifest. REMIND input files and ecoinvent "
        "inventories are restricted local sources.", y + 3,
    )
    y = r.heading("Authoritative local documentation", y + 7)
    y = r.para(
        "TRAILS: docs/stock-vintage-pilot-review.md and "
        "docs/stock-vintage-development.md.\n"
        "Premise implementation worktree: docs/development/stock-asset-temporal/ "
        "contains implementation-status.md, the vehicle/PV/CCGT/water methods, "
        "lifecycle contract, IAM assessment and 23-group search strategy.", y,
        size=9, leading=13.5,
    )
    y = r.heading("Pinned private archive", y + 12)
    y = r.para(
        "~/Library/Application Support/trails/stock-vintage-evidence/\n"
        "2026-10-10-pilots\n"
        "Manifest SHA-256:\n"
        "c15467bdb209ece5c556248c16edc6257dae\n"
        "034dd121b9a6c935b4f5fe726231",
        y, size=8.8, leading=13,
    )
    y = r.note(
        "The archive comprises 163 verified artifacts (1,128,568,509 bytes). "
        "Its sibling verification directory records independent hash checks "
        "and relocated smoke results. The archive's 12 full-LCI reports and "
        "399 tests establish implementation evidence, not empirical certainty.", y + 8,
    )
    y = r.heading("Rebuild this report", y + 7)
    y = r.para(
        "From the TRAILS checkout, using the existing Python environment:\n"
        "python dev/stock_vintage/build_review_pdf.py\n"
        "Optional arguments: --archive PATH --output-dir PATH",
        y, size=9, leading=14,
    )
    r.note(
        "Outputs include the PDF, individual vector pages, a distribution "
        "overview (PDF/PNG), annual profile values, profile statistics and "
        "report-provenance.json. The provenance file lists every consumed "
        "archive hash and the checks performed. Keep derived outputs local.", y + 8,
    )
    r.finish("references-and-reproduction")


def standalone_figure(e, output):
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(3, 2, figsize=(11, 11.5))
    entries = [(k, "pipes", v[0]) for k, v in PILOTS.items()]
    entries.append(("water", "storage", "Water storage"))
    for ax, (key, asset, title) in zip(axes.flat, entries):
        for year in (2022, 2025, 2030):
            pmf(ax, e.annual(key, year, asset=asset), str(year), YEAR_COLORS[year],
                "--" if year == 2025 else "-")
        style_distribution(ax)
        ax.set_xlabel("Cohort year", fontsize=10)
        ax.set_title(title, loc="left", fontsize=12)
        ax.grid(axis="y", color="#E4EBEF", linewidth=0.7)
        ax.spines[["top", "right"]].set_visible(False)
        ax.axvspan(2022.5, 2030.5, color=LIGHT, zorder=0)
        ax.legend(frameon=False, fontsize=9)
    fig.suptitle(
        "Generated annual stock/service vintage distributions\n"
        "Primary cases • SSP2-NPi2025 for vehicles/power; "
        "constant-stock renewal for water",
        fontsize=13, color=INK, y=0.98,
    )
    fig.subplots_adjust(left=0.075, right=0.985, top=0.90, bottom=0.08,
                        hspace=0.35, wspace=0.18)
    fig.text(0.5, 0.016,
             "Every profile sums to 100%. Historical water bands are annualised; "
             "post-2022 additions are reconstructed. Local research evidence.",
             ha="center", fontsize=9, color=MUTED)
    for extension in ("pdf", "png"):
        fig.savefig(output / f"vintage-distributions-overview.{extension}", dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, default=DEFAULT_ARCHIVE)
    parser.add_argument("--output-dir", type=Path,
                        default=ROOT / "build/reports/stock-vintage-pilots")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(args.output_dir / ".mplconfig"))
    import matplotlib

    matplotlib.use("Agg")
    matplotlib.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9,
        "axes.labelcolor": INK, "text.color": INK,
        "axes.titlecolor": INK, "pdf.fonttype": 42,
        "savefig.facecolor": "white",
    })
    evidence = Evidence(args.archive)
    for check in evidence.manifest["corrected_lci_checks"]:
        if not check["signed_full_lci_passed"]:
            raise ValueError("An archived full signed LCI check did not pass")
        evidence.read(f"pilots/{check['pilot']}/{check['file']}")
    if sum(x["tests"] for x in evidence.manifest["regression_suites"]) != 399:
        raise ValueError("Unexpected archived regression count")
    report = Report(args.output_dir)
    try:
        opening_pages(report, evidence)
        pilot_pages(report, evidence)
        analytical_pages(report, evidence)
        closing_pages(report, evidence)
    finally:
        report.close()
    if report.page_number != 17:
        raise ValueError("Unexpected report pagination")
    standalone_figure(evidence, args.output_dir)
    save_csv(args.output_dir / "annual-profile-values.csv", evidence.profile_rows)
    save_csv(args.output_dir / "profile-statistics.csv", evidence.summaries)
    (args.output_dir / "report-text.txt").write_text("\n\n".join(report.text_log))
    provenance = {
        "report": report.pdf_path.name,
        "report_sha256": sha256(report.pdf_path),
        "script_sha256": sha256(Path(__file__)),
        "archive": str(args.archive.resolve()),
        "archive_manifest_sha256": EXPECTED_MANIFEST,
        "consumed_archive_files": evidence.used,
        "documentation_sources": {
            str(path): sha256(path)
            for path in [
                ROOT / "docs/stock-vintage-pilot-review.md",
                ROOT / "docs/stock-vintage-development.md",
                *[
                    ROOT.parent / "premise-stock-vintage-pilots"
                    / "docs/development/stock-asset-temporal" / name
                    for name in (
                        "implementation-status.md", "iam-scenario-assessment.md",
                        "group-solutions.md", "vehicle-pilot-method.md",
                        "pv-pilot-method.md", "ccgt-pilot-method.md",
                        "water-pilot-method.md",
                    )
                ],
            ]
            if path.exists()
        },
        "checks": evidence.validation,
        "profile_count": len(evidence.summaries),
        "annual_probability_mass_rows": len(evidence.profile_rows),
        "pages": report.page_titles,
        "figure_count_in_report": report.figure_number,
        "python": sys.version.split()[0],
        "matplotlib": matplotlib.__version__,
        "numpy": np.__version__,
        "scope": "Report rendering from existing evidence; no LCI rerun",
        "rights": "Restricted local derived evidence; do not publish or redistribute",
    }
    (args.output_dir / "report-provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n"
    )
    (args.output_dir / "README.md").write_text(
        "# Stock-vintage pilot report\n\n"
        "Open `stock-vintage-pilots-report.pdf` for the 17-page review.\n"
        "`vintage-distributions-overview.pdf` and `.png` collect the six "
        "parent profiles (five pilots, with separate water pipes/storage).\n\n"
        "All plotted values come from hash-checked archived evidence. "
        "`annual-profile-values.csv` includes every saved scenario, case, "
        "service year and lifecycle role; `profile-statistics.csv` contains "
        "their mass and mean event year. `report-provenance.json` pins inputs "
        "and describes the fresh report checks. Individual report pages are "
        "available as vector PDFs in `figures/`.\n\n"
        "Keep these restricted IAM-/inventory-derived outputs local and "
        "outside version control. Rebuild from the TRAILS checkout with "
        "`python dev/stock_vintage/build_review_pdf.py`.\n"
    )
    print(json.dumps({
        "pdf": str(report.pdf_path), "pages": report.page_number,
        "figures": report.figure_number,
        "profiles_checked": len(evidence.summaries),
        "input_files_verified": len(evidence.used),
    }))


if __name__ == "__main__":
    main()
