# `TRAILS`: Temporal Routing And Aggregation of Impacts across Life-cycle Systems

<p align="center">
  <img src="https://github.com/Laboratory-for-Energy-Systems-Analysis/trails/blob/main/assets/permanent/trails_logo_light_gray_bg_dark_frame.png" height="110"/>
</p>

[![PyPI version](https://badge.fury.io/py/trails.svg)](https://badge.fury.io/py/trails)

`TRAILS` is a Python library for **temporal and prospective Life Cycle Assessment (LCA)**, 
designed to analyze how environmental impacts propagate through **exchanges with temporal distributions**.

It provides a formal framework for **temporal graph traversal** of supply chains, enabling 
the routing, aggregation, and attribution of impacts across **multiple time horizons, scenarios, 
and technological transitions**.

`TRAILS` integrates data packages produced by `premise`.

---

## Usage

Below is a minimal example that loads a Frictionless data package, runs a
temporal LCA, and plots the resulting impact time series.

```python
from datapackage import Package

from trails import Trails, lca, get_lcia_method_names, plot_temporal_scores

# Load a Frictionless data package exported by premise (or compatible tooling)
package = Package("path/to/datapackage.json")

# Initialize the TRAILS wrapper (with optional annual interpolation)
trails = Trails(package, interpolate_annual=True)

# Pick an activity index from the metadata
activity_indices = next(iter(trails.activity_indices.values()))
start_act_idx = next(iter(activity_indices.keys()))

# Choose an LCIA method bundled with TRAILS
method = get_lcia_method_names(ei_version="3.11")[0]

# Run a temporal LCA
results = lca(
    trails=trails,
    start_year=2030,
    start_act_idx=start_act_idx,
    methods=[method],
    max_depth=2,
)

# Plot temporal impact scores
fig = plot_temporal_scores(results, trails, method_label=method)
fig.show()
```

---

## Motivation

Conventional LCA frameworks treat time implicitly or exogenously. Impacts are typically 
computed for a single static system, even when future scenarios or dynamic technologies 
are considered.

`TRAILS` addresses this limitation by introducing:

* Handling of **temporal dimensions** in technosphere and biosphere matrices  
* **Time-aware routing of exchanges** across supply chains  
* **Scenario-dependent inventories and impacts**  

Instead of asking *“What is the impact of this system?”*, `TRAILS` allows you to ask:

> *When, where, and through which pathways do impacts occur across the life cycle?*

---

## Core Concepts

### 1. Temporal graph traversal
Life-cycle systems are represented as **time-indexed graphs**, where exchanges may occur at 
different points in time relative to the functional unit.

### 2. Routing of impacts
Impacts are **routed along supply-chain paths**, allowing attribution to:
* specific suppliers,
* specific time periods,
* specific traversal depths.

### 3. Aggregation across scenarios and horizons
Impacts can be aggregated or compared across:
* years (e.g., 2020 → 2050 → 2100),
* scenarios (e.g., SSPs, decarbonization pathways),
* temporal horizons (short-term vs long-term effects).

---

## Key Features

* Temporal LCA engine with explicit time handling  
* Deep supply-chain traversal  
* Scenario-aware computation

---

## Installation

```bash
pip install trails
```

---

## Documentation

https://trails.readthedocs.io/en/latest/index.html

---

## License

MIT License.
