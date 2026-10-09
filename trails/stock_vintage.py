"""Strict, opt-in annual stock-vintage profiles independent of matrix anchors.

Profiles store calendar event years. Bindings use activity identities, never
maintained matrix indices. A missing service year is an error, not permission to
reuse an endpoint profile. The legacy temporal exchange parser is unchanged.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
from numbers import Real
from pathlib import Path
from typing import Any

from .temporal_distributions import TemporalExchange

IDENTITY_FIELDS = ("name", "reference product", "unit", "location")
ROLES = {
    "existing_asset_service",
    "new_asset_component",
    "replacement_component",
    "lifecycle_service",
}


def _identity(value: dict) -> tuple[str, ...]:
    if not isinstance(value, dict):
        raise ValueError("A stock-vintage binding needs an activity identity")
    fields = tuple(value.get(key) for key in IDENTITY_FIELDS)
    if any(not isinstance(item, str) or not item.strip() for item in fields):
        raise ValueError(f"Incomplete stock-vintage activity identity: {fields!r}")
    return fields


def _year(value: Any) -> int:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"Expected an integer calendar year, got {value!r}")
    if not math.isfinite(value) or int(value) != value:
        raise ValueError(f"Expected an integer calendar year, got {value!r}")
    return int(value)


@dataclass(frozen=True)
class AnnualProfile:
    """Validated manufacture/event timing under common amortisation."""

    profile_id: str
    event_role: str
    asset_unit: str
    service_unit: str
    years: dict[int, TemporalExchange]
    provenance: dict


class StockVintageProfiles:
    """Validated annual profiles bound to exact package activity identities."""

    def __init__(self, payload: dict, activity_indices: dict, *, sha256: str):
        if not isinstance(payload, dict) or payload.get("schema_version") != 1:
            raise ValueError("Unsupported stock-vintage schema_version")
        self.sha256 = sha256
        self.profiles: dict[str, AnnualProfile] = {}
        self.bindings: dict[tuple[str, int, int], str] = {}
        raw_profiles = payload.get("profiles")
        if not isinstance(raw_profiles, list) or not raw_profiles:
            raise ValueError("Stock-vintage profiles must be a nonempty list")
        for profile in raw_profiles:
            if not isinstance(profile, dict):
                raise ValueError("A stock-vintage profile must be an object")
            profile_id = profile.get("id")
            if not isinstance(profile_id, str) or not profile_id.strip():
                raise ValueError("A stock-vintage profile needs a nonempty id")
            if profile_id in self.profiles:
                raise ValueError(f"Duplicate stock-vintage profile: {profile_id}")
            role = profile.get("event_role")
            if role not in ROLES:
                raise ValueError(f"Unsupported stock-vintage event_role: {role!r}")
            if profile.get("allocation_basis") != "common_amortisation":
                raise ValueError("Stock-vintage v1 requires common_amortisation")
            for field in ("asset_unit", "service_unit"):
                if (
                    not isinstance(profile.get(field), str)
                    or not profile[field].strip()
                ):
                    raise ValueError(f"Profile {profile_id} requires {field}")
            provenance = profile.get("provenance")
            if not isinstance(provenance, dict) or not provenance:
                raise ValueError(f"Profile {profile_id} requires provenance")
            records = profile.get("years")
            if not isinstance(records, list) or not records:
                raise ValueError(f"Profile {profile_id} needs annual records")
            annual = {}
            for record in records:
                if not isinstance(record, dict):
                    raise ValueError("A stock-vintage annual record must be an object")
                year = _year(record.get("service_year"))
                if year in annual:
                    raise ValueError(f"Duplicate service year {year} in {profile_id}")
                events, weights = record.get("event_years"), record.get("weights")
                if not isinstance(events, list) or not isinstance(weights, list):
                    raise ValueError("Calendar event years and weights must be lists")
                if not events or len(events) != len(weights):
                    raise ValueError(
                        "Event years and weights need equal nonzero length"
                    )
                event_years = [_year(value) for value in events]
                if len(set(event_years)) != len(event_years):
                    raise ValueError("Aggregate duplicate calendar event years first")
                if role != "lifecycle_service" and max(event_years) > year:
                    raise ValueError(
                        "An asset cannot be manufactured after its service"
                    )
                if any(
                    isinstance(value, bool)
                    or not isinstance(value, Real)
                    or not math.isfinite(value)
                    or value < 0
                    for value in weights
                ):
                    raise ValueError(
                        "Stock-vintage weights must be finite and nonnegative"
                    )
                if not math.isclose(math.fsum(weights), 1.0, rel_tol=0, abs_tol=1e-10):
                    raise ValueError("Stock-vintage weights must sum to one")
                pairs = sorted(zip(event_years, weights))
                offsets = tuple(event - year for event, _ in pairs)
                annual[year] = TemporalExchange(
                    distribution=6,
                    loc=None,
                    scale=None,
                    offset_min=min(offsets),
                    offset_max=max(offsets),
                    amount_source="port",
                    offsets=offsets,
                    weights=tuple(float(weight) for _, weight in pairs),
                )
            self.profiles[profile_id] = AnnualProfile(
                profile_id,
                role,
                profile["asset_unit"],
                profile["service_unit"],
                annual,
                provenance,
            )

        if not activity_indices:
            raise ValueError("Stock-vintage bindings need package activity indices")
        indices = {}
        for label, activities in activity_indices.items():
            lookup = {}
            for index, metadata in activities.items():
                key = _identity(metadata)
                if key in lookup:
                    raise ValueError(f"Ambiguous activity identity in {label}: {key}")
                lookup[key] = int(index)
            indices[str(label)] = lookup
        raw_bindings = payload.get("bindings")
        if not isinstance(raw_bindings, list) or not raw_bindings:
            raise ValueError("Stock-vintage bindings must be a nonempty list")
        for binding in raw_bindings:
            if not isinstance(binding, dict):
                raise ValueError("A stock-vintage binding must be an object")
            profile_id = binding.get("profile_id")
            if profile_id not in self.profiles:
                raise ValueError(f"Unknown stock-vintage profile: {profile_id!r}")
            caller, supplier = _identity(binding.get("caller")), _identity(
                binding.get("supplier")
            )
            if caller == supplier:
                raise ValueError("A stock-vintage binding cannot target production")
            profile = self.profiles[profile_id]
            if caller[2] != profile.service_unit or supplier[2] != profile.asset_unit:
                raise ValueError("Profile asset/service units do not match the binding")
            for label, lookup in indices.items():
                if caller not in lookup or supplier not in lookup:
                    raise ValueError(
                        f"Unmatched stock-vintage binding in {label}: {binding}"
                    )
                key = (label, lookup[caller], lookup[supplier])
                if key in self.bindings:
                    raise ValueError(f"Conflicting stock-vintage bindings: {key}")
                self.bindings[key] = profile_id

    def get(
        self, service_year: int, template_year: int, caller: int, supplier: int
    ) -> TemporalExchange | None:
        """Resolve an exact service year, using the template only for identities."""
        profile_id = self.bindings.get((str(template_year), int(caller), int(supplier)))
        if profile_id is None:
            return None
        try:
            return self.profiles[profile_id].years[int(service_year)]
        except KeyError as error:
            raise ValueError(
                f"Stock-vintage profile {profile_id!r} has no service year "
                f"{service_year}; "
                "generate an explicit annual profile or documented extrapolation"
            ) from error

    def offset_bounds(self) -> tuple[int, int]:
        """Return support across annual records for inventory horizon accounting."""
        exchanges = [
            tex for profile in self.profiles.values() for tex in profile.years.values()
        ]
        return (
            min(tex.offset_min for tex in exchanges),
            max(tex.offset_max for tex in exchanges),
        )

    def remove_callers(self, labels: list[str], callers: set[int]) -> None:
        """Let explicit foreground replacements supersede package bindings."""
        removed = [
            key for key in self.bindings if key[0] in labels and key[1] in callers
        ]
        for key in removed:
            del self.bindings[key]


def load_stock_vintage_profiles(
    package: Any, activity_indices: dict, *, enabled: bool
) -> StockVintageProfiles | None:
    """Read the declared, hashed resource; require explicit consumer opt-in."""
    spec = getattr(package, "descriptor", {}).get("stock_vintage")
    if spec is None:
        return None
    if not enabled:
        raise ValueError(
            "This package requires annual stock-vintage support; "
            "pass stock_vintage=True"
        )
    if not isinstance(spec, dict) or spec.get("schema_version") != 1:
        raise ValueError("Unsupported stock_vintage package metadata")
    resources = [
        resource
        for resource in package.resources
        if resource.name == spec.get("resource")
    ]
    if len(resources) != 1:
        raise ValueError("Expected exactly one declared stock-vintage resource")
    raw = Path(resources[0].source).read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    if digest != spec.get("sha256"):
        raise ValueError(
            "Stock-vintage resource hash mismatch; regenerate package metadata"
        )
    return StockVintageProfiles(json.loads(raw), activity_indices, sha256=digest)
