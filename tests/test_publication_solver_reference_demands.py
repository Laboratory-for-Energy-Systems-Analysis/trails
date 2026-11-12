"""Check the benchmark's activity-to-reference-function solve boundary."""

from pathlib import Path
import importlib.util
import sys

import networkx as nx
import numpy as np
import pytest
import scipy.sparse as sp
import sparse

from trails import Trails


def benchmark_module():
    path = Path(__file__).resolve().parents[1] / "dev/publication/solver_stress_test.py"
    spec = importlib.util.spec_from_file_location("a18_solver_benchmark", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("scaling", [1.0, -1.0])
def test_benchmark_frontier_reference_demands(example_package, monkeypatch, scaling):
    module = benchmark_module()
    with Trails(example_package, interpolate_annual=False) as model:
        a = np.eye(model.A.shape[1])
        a[1, 1] = 5.0
        a[2, 2] = -4.0
        model.A = sparse.COO.from_numpy(
            np.repeat(a[None, :, :], model.A.shape[0], axis=0)
        )
        model.graph = nx.DiGraph()
        model.graph.add_node(1, year=2005, act_idx=1, frontier_amount=0.6 * scaling)
        model.graph.add_node(2, year=2005, act_idx=2, frontier_amount=-0.25 * scaling)
        monkeypatch.setattr(model, "temporal_routing", lambda **kwargs: None)
        cases, *_ = module._cases_from_routing_with_stats(
            trails=model,
            start_year=2005,
            activity_index=0,
            amount=1.0,
            max_depth=1,
            adaptive_relative_score_cutoff=None,
            adaptive_methods=None,
            selected_years={2005},
        )
        assert len(cases) == 1
        reference = cases[0].activity_demand
        assert reference[1] == pytest.approx(3.0 * scaling)
        assert reference[2] == pytest.approx(-1.0 * scaling)
        rhs = module._activity_demand_to_product_rhs(
            sp.csc_matrix(a.T), reference, activity_offset=0, product_offset=0
        )
        assert rhs == pytest.approx({1: 3.0 * scaling, 2: 1.0 * scaling})
