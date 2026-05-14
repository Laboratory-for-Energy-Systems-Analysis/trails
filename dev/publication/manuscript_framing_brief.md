# Manuscript Framing Brief

This brief distills the completed questionnaire into a concrete plan for
rewriting the manuscript. It is intended to guide a full rewrite of
`manuscript.md`, not to preserve the current section text.

## 1. Working Title

Preferred title:

> A graph-matrix hybrid approach for deep temporalization in time-explicit LCA

Possible variants:

- Deep temporalization in time-explicit life cycle assessment with TRAILS
- TRAILS: graph-matrix hybrid modeling for deep temporalization in prospective
  life cycle assessment
- Connecting prospective databases through time for deep temporalization in LCA

Recommendation: keep "deep temporalization" in the title or subtitle. This is
the manuscript's differentiating concept.

## 2. Target Paper Type and Venue

Target venue: Journal of Industrial Ecology.

Paper type: combined methods-and-software paper, supported by application
demonstrations.

Expected tone:

- accessible to applied LCA practitioners and prospective LCA users;
- enough mathematical detail to make the method credible;
- not a pure software note;
- not a case-study paper with substantive conclusions about specific products;
- moderate length, around 5000-7000 words.

The manuscript should foreground reproducibility and open-source
implementation. Jupyter notebooks should reproduce all main result figures and
be provided as electronic supplementary information.

## 3. Central Claim

Draft central claim:

> TRAILS enables easy and computationally efficient time-explicit LCA by
> combining prospective LCA databases, dynamic LCA, and deep temporalization,
> allowing long-lived systems and systems embedded in fast-evolving background
> systems to be assessed with temporal dynamics throughout the supply chain.

Sharper version for the abstract:

> TRAILS extends time-explicit LCA from foreground temporalization to deep
> temporalization by representing year-specific foreground and background
> inventories as a connected three-dimensional technosphere and combining
> temporal graph routing with year-wise matrix solving.

One-week reader memory:

> TRAILS is an easy-to-use, fast, graph-based, time-explicit LCA tool that
> combines prospective and dynamic LCA, supports deep temporalization, and can
> couple inventories to a climate emulator.

## 4. Contribution Hierarchy

Primary contribution:

1. Deep temporalization: temporal distributions can occur throughout the supply
   chain, including background processes, not only in foreground inventories.

Secondary contributions:

2. Graph-matrix hybrid algorithm: temporal routing builds time-indexed demands;
   year-wise matrix solves avoid one large time-expanded technosphere solve.
3. Two temporal amount semantics: after routing an exchange over time, TRAILS
   can use either the anchor-year amount (`port`) or matrix-sourced target-year
   amounts (`matrix`).

Demonstration-level contribution:

4. Coupling to FaIR climate response modeling, with background climate
   scenarios that can be made consistent with premise scenario assumptions.

The FaIR coupling should not be presented as the core methodological novelty.
It should be used as evidence that TRAILS can feed dynamic impact models once a
time-resolved inventory has been produced.

## 5. Claims to Avoid

The manuscript should not claim:

- that TRAILS is better than `bw_timex`;
- that TRAILS is the first tool combining prospective and dynamic LCA;
- that TRAILS replaces Brightway;
- that TRAILS performs uncertainty sampling;
- that TRAILS captures sub-annual dynamics;
- that case-study results are definitive product conclusions unless the case
  studies are designed and validated for that purpose.

The safer claim is:

> Existing time-explicit tools have already combined temporal distributions with
> time-specific background databases. TRAILS addresses a complementary gap:
> representing foreground and background inventories as a connected
> three-dimensional technosphere so that temporalized exchanges can be
> propagated deeper into the supply chain.

## 6. Audience Stance

Assume readers know LCA and may know prospective LCA, but do not assume they are
fluent in Brightway, premise, or matrix notation.

Main text should briefly explain:

- static matrix LCA at the level of `A`, `B`, demand vector, and inventory;
- prospective LCA databases as year- and scenario-specific snapshots;
- dynamic LCA as temporal differentiation of demands and emissions;
- time-explicit LCA as the combination of temporal distributions and
  time-specific database selection.

Detailed implementation, data schema, and premise temporal parameterization
should go to ESI.

## 7. Literature Positioning

The introduction should build toward the gap of deep temporal traversal.

Recommended positioning:

- Beloin-Saint-Pierre et al. and related ESPA work: establish process-relative
  temporal distributions and traversal/convolution ideas.
- Tiruta-Barna, Pigne, and DyPLCA: establish operational full dynamic LCI and
  dynamic LCA tooling.
- Cardellini et al. and Temporalis/bw_temporalis: establish Brightway-based
  graph traversal and dynamic characterization.
- Muller et al. and Diepers et al. / `bw_timex`: establish modern
  time-explicit LCA combining temporal distributions with time-specific
  background databases.
- premise: establish prospective LCA database generation as the natural source
  of year- and scenario-specific databases.
- TRAILS: connects those year-specific databases into a single
  three-dimensional technosphere that can be traversed deeply.

Comparator framing:

> `bw_timex` is the closest comparator and should be treated respectfully as a
> conceptual predecessor and adjacent solution. The difference to emphasize is
> not "TRAILS is better", but "TRAILS addresses deep temporalization of
> connected year-specific databases."

## 8. Proposed Manuscript Structure

### Abstract

Structure:

1. Problem: dynamic and prospective LCA are both needed for long-lived systems,
   but most workflows temporalize mainly the foreground or connect to
   time-specific backgrounds without deep temporal propagation.
2. Method: TRAILS represents inventories as a three-dimensional
   scenario-year technosphere and combines temporal graph routing with
   year-wise matrix solving.
3. Demonstration: 4-5 case studies compare static, foreground-only
   temporalized, and deeply temporalized results.
4. Finding: deep temporalization can materially change time-resolved LCA scores
   in sensitive systems.
5. Availability: open-source Python implementation and reproducible notebooks.

### 1. Introduction

Purpose:

- motivate time-explicit LCA for long-lived systems and systems embedded in
  rapidly changing background systems;
- explain the distinction between prospective database evolution and dynamic
  temporal distributions;
- position existing dynamic and time-explicit tools;
- define the remaining gap as deep temporalization through connected
  year-specific databases.

End with contributions:

1. connected 3D technosphere representation for deep temporalization;
2. graph-matrix hybrid routing and solving;
3. port and matrix temporal amount semantics;
4. open-source software with reproducible demonstrations and optional FaIR
   coupling.

### 2. Methods

Recommended subsections:

2.1 Static and prospective matrix LCA notation

- Use `A_y` for technosphere matrix in year `y`.
- Use `B_y` for intervention or biosphere matrix in year `y`.
- Use `s_y` for the demand or scaling vector.
- Use `f_y = A_y^{-1} s_y` for the supply/scaling solution.
- Use `G_y = B_y diag(f_y)` or the equivalent inventory expression.
- Use `q` or `Q` for characterization factors and `H_y` for characterized
  results.

Note: check notation carefully during drafting. In standard matrix LCA, the
usual relation is supply vector `x = A^{-1} f`, where `f` is final demand. The
questionnaire uses `f = A^-1 * s`. We should choose one convention and define
it explicitly to avoid reviewer confusion.

2.2 Three-dimensional scenario-year matrices

- Explain `A_y` and `B_y` as slices of a sparse 3D matrix tensor.
- Explain that years are typically interpolated between scenario database years,
  for example 2005-2100.
- Keep interpolation details brief in the main text.

2.3 Temporal exchange distributions

- Define relative temporal distributions on exchanges.
- Define annual offsets.
- Define `port` and `matrix` semantics as a methodological feature.

2.4 Temporal routing

- Define nodes as `(process, year)`.
- Define edges as technosphere exchanges with optional temporal distributions.
- Define edge weights as exchange amounts propagated over time.
- Define depth, cutoff, roots, and frontier.
- Present this section with an algorithm box.

2.5 Year-wise matrix solving

- Explain that routing produces time-indexed frontier demand vectors.
- For each year, TRAILS solves the corresponding matrix system.
- Explain why this avoids a single massive time-expanded solve.
- Mention the iterative Krylov solver with Jacobi preconditioning in one
  sentence; move solver details to ESI.

2.6 Inventory, characterization, and optional climate response

- Explain accumulation by emission year.
- Explain LCIA application.
- Present FaIR coupling as optional demonstration, not as the core method.

### 3. Implementation

Keep concise in the main text:

- open-source Python package;
- sparse matrices;
- Frictionless datapackage input;
- premise as one supported input pathway;
- Jupyter notebooks for figure reproduction;
- Brightway compatibility through data conversion, but TRAILS is not a
  replacement for Brightway.

Move to ESI:

- full datapackage schema;
- how premise adds temporal distributions across future ecoinvent exchanges;
- solver implementation details;
- cache behavior;
- additional software API details.

### 4. Demonstration Case Studies

Main purpose:

- show that deep temporalization can matter;
- show how to use TRAILS;
- avoid overclaiming substantive conclusions about specific products.

Recommended design:

- include 4-5 systems chosen because they are likely sensitive to upstream
  temporalization;
- compare at least three modeling modes:
  1. static LCA;
  2. foreground-only temporalization or time-specific background linking;
  3. deep temporalization with TRAILS;
- vary routing depth to show convergence or sensitivity;
- include temporal distribution sensitivity for at least one case.

Current uncertainty:

- passenger cars and carbon removal versus carbon avoidance may or may not be
  retained;
- case selection should be based on screening for sensitivity to deep
  temporalization.

### 5. Results

Organize results by the manuscript's core claim, not by software feature.

Recommended result subsections:

5.1 Deep temporalization changes time-resolved scores in selected systems

- main multi-case figure;
- show static, foreground-only, and deep temporalization results.

5.2 Routing depth reveals where temporal dynamics enter the supply chain

- show score evolution as routing depth increases;
- identify cases where foreground temporalization is sufficient and cases where
  it is not.

5.3 Temporal amount semantics can affect conclusions

- show `port` versus `matrix` for one or two examples where the distinction
  matters.

5.4 Optional climate response demonstration

- use FaIR only if it reinforces timing relevance and scenario consistency.

### 6. Discussion

Recommended discussion claims:

- TRAILS is most useful for long-lived systems, infrastructure-dominated supply
  chains, carbon storage/removal systems, and prospective LCA studies where
  upstream systems evolve quickly.
- TRAILS may not be worth the extra complexity for short-lived systems or
  systems with impacts concentrated in the foreground.
- TRAILS complements Brightway and `bw_timex`; it is not a replacement.
- Deep temporalization should be applied where it changes decision-relevant
  results, not automatically for every study.

Prominent limitations:

- deterministic calculations; no uncertainty sampling;
- annual time steps; no sub-annual dynamics;
- case studies demonstrate capability rather than definitive product rankings;
- metadata and temporal distribution assumptions matter;
- some Brightway context is abstracted away when operating through TRAILS.

Future work:

- dynamic assessment for other impact areas beyond climate;
- time-dependent air pollution and ecotoxicity pathways;
- pollutant migration through media and delayed exposure;
- uncertainty propagation;
- sub-annual dynamics where relevant;
- improved database linking and richer temporal metadata.

## 9. Essential Figures

Figure 1: Conceptual method figure

- Optional but recommended if space allows.
- Show datapackage input, connected year-specific matrices, temporal routing,
  frontier demand vectors, year-wise solves, inventory aggregation, LCIA, and
  optional FaIR.

Figure 2: Main result figure across 4-5 cases

- Required.
- Compare static, foreground-only temporalization, and deep temporalization.
- Could use panels, one per case.
- This is the main evidence that TRAILS adds something distinct.

Figure 3: Routing depth sensitivity

- Required or strongly recommended.
- Show how results evolve as routing depth increases.
- Supports the "deep" part of deep temporalization.

Figure 4: Temporal amount semantics

- Recommended.
- Show `port` versus `matrix` behavior in a selected case.
- Makes this methodological contribution concrete.

Figure 5: FaIR climate response demonstration

- Optional.
- Include only if it reinforces the core story without diluting the paper.

No performance figure is needed unless reviewers or results make scalability a
central concern.

## 10. Case-Study Selection Criteria

Because current case studies are not fixed, choose cases by screening for:

- long service lives;
- infrastructure-heavy supply chains;
- relevant upstream capital goods;
- dependence on energy systems that evolve strongly over time;
- biomass growth, fleet age, lifetime, or storage dynamics;
- carbon storage, carbon removal, or delayed release/removal;
- clear contrast between foreground-only and deep temporalization.

A publishable case result is not "this product is better." It is:

> A clear demonstration that deeper temporalization changes the time-resolved
> score or interpretation relative to static or foreground-only temporalization.

Recommended next analytical step:

- run a screening notebook over candidate systems;
- calculate scores at multiple routing depths;
- rank systems by the difference between foreground-only and deeper
  temporalization;
- choose 4-5 cases with diverse mechanisms.

## 11. Main Text Versus Supplement

Main text should include:

- motivation and gap;
- concise method with equations and algorithm box;
- deep temporalization definition;
- port/matrix semantics;
- main figures for selected cases;
- limitations and scope;
- software availability and reproducibility.

Supplement should include:

- full Frictionless datapackage schema;
- premise temporal exchange parameterization;
- solver and cache details;
- extended sensitivity analyses;
- additional case results;
- notebooks and scripts.

## 12. Draft Highlight Statements

Potential highlights:

- TRAILS enables deep temporalization across foreground and background supply
  chains in time-explicit LCA.
- A graph-matrix hybrid algorithm combines temporal routing with year-wise
  matrix solving.
- Connected prospective databases allow temporal exchanges to propagate beyond
  foreground inventories.
- Routing depth and temporal amount semantics can change time-resolved LCA
  results.
- Open-source notebooks reproduce the demonstration figures.

## 13. Draft Keywords

- Time-explicit life cycle assessment
- Dynamic life cycle assessment
- Prospective life cycle assessment
- Deep temporalization
- Graph traversal
- Matrix LCA
- Brightway
- premise
- Climate response modeling

## 14. Abstract Skeleton

Draft skeleton:

> Time-explicit life cycle assessment is increasingly needed for long-lived
> systems and systems embedded in rapidly changing background economies.
> Existing dynamic and prospective LCA tools can represent process timing and
> time-specific databases, but temporalization often remains concentrated in the
> foreground or in links to background snapshots. Here we present TRAILS, an
> open-source Python tool for deep temporalization in time-explicit LCA. TRAILS
> represents foreground and background inventories as connected
> scenario-year-specific technosphere and biosphere matrices, applies temporal
> graph routing to exchanges throughout the supply chain, and solves the
> resulting frontier demands with year-wise matrix systems. The method supports
> two temporal amount semantics, allowing users either to distribute
> anchor-year exchange amounts or to draw target-year amounts from the
> scenario-year matrices. We demonstrate the approach on selected systems for
> which upstream temporalization changes time-resolved LCA results compared with
> static and foreground-only temporalized calculations. The implementation is
> distributed as open-source software, with reproducible notebooks for all main
> figures. TRAILS complements existing time-explicit LCA tools by making deep
> temporalization of connected prospective databases accessible to applied LCA
> practitioners.

This is a skeleton, not final prose. It needs concrete case-study results before
submission.

## 15. Immediate Rewrite Decisions

Use these decisions when rewriting `manuscript.md`:

- Replace the current placeholder introduction with a gap-driven introduction
  centered on deep temporalization.
- Keep "graph-matrix hybrid" as the method label, but make "deep
  temporalization" the conceptual hook.
- Treat `bw_timex` and DyPLCA respectfully as predecessors and adjacent tools.
- Present premise as a supported and especially important input pathway, not as
  the only way to build a TRAILS datapackage.
- Do not over-explain the datapackage schema in the main text.
- Do include equations and an algorithm box.
- Do not include a runtime/performance figure unless later results demand it.
- Use the case studies to demonstrate sensitivity to deep temporalization, not
  to claim product rankings.

