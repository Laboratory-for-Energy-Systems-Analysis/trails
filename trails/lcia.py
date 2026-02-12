from __future__ import annotations

import json

import numpy as np
from scipy import sparse
from scipy.sparse import csr_matrix

from .filesystem_constants import DATA_DIR

import logging

logger = logging.getLogger(__name__)

LCIA_METHODS_EI310 = DATA_DIR / "lcia_ei310.json"
LCIA_METHODS_EI311 = DATA_DIR / "lcia_ei311.json"


_LCIA_METHODS_CACHE = {}


def get_lcia_method_names(ei_version: str = "3.11") -> list[str]:
    """Get lcia method names.

    :param ei_version: Value for `ei_version`.
    :type ei_version: str
    :returns: Return value.
    :rtype: list[str]"""

    if ei_version == "3.11":
        filepath = LCIA_METHODS_EI311
    else:
        filepath = LCIA_METHODS_EI310

    with open(filepath, "r") as f:
        data = json.load(f)

    return [" - ".join(x["name"]) for x in data]


def format_lcia_method_exchanges(
    method: dict,
) -> dict[tuple[str, str, str], float]:
    """Format lcia method exchanges.

    :param method: Value for `method`.
    :type method: dict
    :returns: Return value.
    :rtype: dict[tuple[str, str, str], float]"""

    return {
        (
            x["name"],
            x["categories"][0],
            x["categories"][1] if len(x["categories"]) > 1 else "unspecified",
        ): x["amount"]
        for x in method["exchanges"]
    }


def get_lcia_methods(
    methods: list[str] | None = None, ei_version: str = "3.11"
) -> dict[str, dict[tuple[str, str, str], float]]:
    """Get lcia methods.

    :param methods: Value for `methods`.
    :type methods: list[str] | None
    :param ei_version: Value for `ei_version`.
    :type ei_version: str
    :returns: Return value.
    :rtype: dict[str, dict[tuple[str, str, str], float]]"""
    key = (ei_version, tuple(methods) if methods else None)
    if key in _LCIA_METHODS_CACHE:
        return _LCIA_METHODS_CACHE[key]

    filepath = LCIA_METHODS_EI311 if ei_version == "3.11" else LCIA_METHODS_EI310

    with open(filepath, "r") as f:
        data = json.load(f)

    if methods:
        data = [x for x in data if " - ".join(x["name"]) in methods]

    out = {" - ".join(x["name"]): format_lcia_method_exchanges(x) for x in data}
    _LCIA_METHODS_CACHE[key] = out
    return out


def fill_characterization_factors_matrices(
    methods: list[str],
    biosphere_matrix_dict: dict[int, int],
    biosphere_dict: dict[tuple[str, str, str], int],
    debug: bool = False,
) -> csr_matrix:
    """Fill characterization factors matrices.

    :param methods: Value for `methods`.
    :type methods: list[str]
    :param biosphere_matrix_dict: Value for `biosphere_matrix_dict`.
    :type biosphere_matrix_dict: dict[int, int]
    :param biosphere_dict: Value for `biosphere_dict`.
    :type biosphere_dict: dict[tuple[str, str, str], int]
    :param debug: Value for `debug`.
    :type debug: bool
    :returns: Return value.
    :rtype: csr_matrix"""

    if debug:
        logger.info(
            "LCIA: building CF matrix for methods=%d biosphere_flows=%d",
            len(methods),
            len(biosphere_matrix_dict),
        )

    lcia_data = get_lcia_methods(methods=methods)

    # Prepare data for efficient creation of the sparse matrix
    data = []
    rows = []
    cols = []
    cfs = []

    matched, unmatched = 0, 0

    for m, method in enumerate(methods):
        method_data = lcia_data[method]

        for flow_name in method_data:
            if flow_name in biosphere_dict:
                idx = biosphere_dict[flow_name]
                if idx in biosphere_matrix_dict:
                    data.append(method_data[flow_name])
                    rows.append(biosphere_matrix_dict[idx])
                    cols.append(m)
                    cfs.append((method, flow_name, idx, method_data[flow_name]))
                    matched += 1
                else:
                    unmatched += 1
            else:
                unmatched += 1

    if debug:
        logger.info("LCIA: CF nonzeros=%d", len(data))
        if len(data) == 0:
            logger.warning(
                "LCIA: CF matrix has zero entries -> all scores will be zero."
            )

        logger.info("LCIA: matched_flows=%d unmatched_flows=%d", matched, unmatched)

    # Efficiently create the sparse matrix
    matrix = sparse.csr_matrix(
        (data, (cols, rows)),
        shape=(len(methods), len(biosphere_matrix_dict)),
        dtype=np.float64,
    )

    if debug:
        # sort l by method and flow
        cfs = sorted(cfs, key=lambda x: (x[0], x[1]))
        for x in cfs:
            method, flow, f, value = x
            logger.info(
                f"LCIA method: {method}, Flow: {flow}, Index: {f}, Value: {value}"
            )

    return matrix
