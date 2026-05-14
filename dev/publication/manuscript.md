# A graph-matrix hybrid approach for deep temporalization in time-explicit LCA

Romain Sacchi<sup>1,2</sup>\*, Tom M. Terlouw<sup>1,2</sup>, Arthus Jakobs<sup>1</sup>, Karin Treyer<sup>1</sup>, Alvaro Hahn-Menacho<sup>1,2</sup>, Christian Bauer<sup>1</sup>

<sup>1</sup> PSI Centers for Nuclear Engineering and Sciences and for Energy and Environmental Sciences, Paul Scherrer Institute, Villigen, Switzerland

<sup>2</sup> Chair of Energy Systems Analysis, Institute of Energy and Process Engineering, ETH Zürich, Zürich, Switzerland

\* Corresponding author: romain.sacchi@psi.ch

## Abstract

Time-explicit life cycle assessment (LCA) is increasingly needed for long-lived systems and systems embedded in rapidly changing background economies. Existing dynamic and prospective LCA approaches can represent process timing, emission timing, and time-specific background databases, but practical workflows often concentrate temporalization in the foreground or connect foreground events to background snapshots without propagating temporalized exchanges deeper into the supply chain. Here we present TRAILS, an open-source implementation of a graph-matrix hybrid method for deep temporalization in time-explicit LCA. TRAILS represents foreground and background inventories as connected scenario-year-specific technosphere and biosphere matrices, applies temporal graph routing to exchanges throughout the supply chain, and solves the resulting frontier demands with year-wise matrix systems. The method supports two interpretations of temporally distributed exchange amounts, allowing users either to distribute anchor-year exchange amounts or to re-evaluate exchange amounts in the target years where they occur. We demonstrate the approach using selected systems for which upstream temporalization is expected to affect time-resolved LCA results compared with static and foreground-only temporalized calculations. TRAILS complements existing time-explicit LCA tools by making deep temporalization of connected prospective databases accessible to applied LCA practitioners, while preserving the need to justify added temporal complexity by the decision context and system understanding.

## Keywords

Time-explicit life cycle assessment; dynamic life cycle assessment; prospective life cycle assessment; deep temporalization; graph traversal; matrix LCA; Brightway; premise; climate response modeling

## Highlights

- TRAILS enables deep temporalization across foreground and background supply chains.
- A graph-matrix hybrid algorithm combines temporal routing with year-wise LCA solves.
- Connected prospective databases allow temporal exchanges to propagate beyond foreground inventories.
- Alternative amount interpretations distinguish fixed anchor-year exchange coefficients from target-year re-evaluated coefficients.
- Open-source notebooks will reproduce the demonstration figures and sensitivity analyses.

## 1. Introduction

Life cycle assessment (LCA) is commonly applied as if all processes in a product system occur at a single point in time and as if the background economy remains fixed during the life cycle of the assessed system (Beloin-Saint-Pierre et al. 2020). These assumptions are often acceptable for short-lived products or for questions where temporal dynamics do not affect the decision. They become less satisfactory for long-lived infrastructure, energy systems, carbon dioxide removal, biogenic carbon storage, systems with delayed end-of-life treatment, and technologies embedded in rapidly changing supply chains (Kendall 2012; Brandão et al. 2013; Lueddeckens et al. 2020). The timing of material production, operation, replacement, disposal, emissions, and removals occuring within the studied system, but also in the rest of the modelled economy, can affect both the inventory and the impact assessment.

Dynamic LCA has addressed the temporal distribution of processes and emissions by representing when exchanges occur and, in some cases, by applying time-dependent characterization. Early work established that temporally differentiated inventories and process-relative temporal information can be propagated through product systems (Levasseur et al. 2010; Beloin-Saint-Pierre et al. 2014). Subsequent methods and tools, including ESPA, DyPLCA, and Temporalis, developed graph traversal, convolution, and operational approaches for dynamic life cycle inventory and dynamic impact assessment (Beloin-Saint-Pierre et al. 2014; Tiruta-Barna et al. 2016; Cardellini et al. 2018; Pigné et al. 2020). Reviews have shown that temporal considerations in LCA span both life cycle inventory and life cycle impact assessment, and that their relevance depends on the system, the temporal profile of exchanges, and the impact category considered (Beloin-Saint-Pierre et al. 2020; Lueddeckens et al. 2020; Sohn et al. 2020).

In parallel, prospective LCA has developed workflows for representing how background systems evolve under future scenarios (Mendoza Beltran et al. 2020; Sacchi et al. 2022). Tools such as Premise generate year- and scenario-specific LCA databases from integrated assessment model outputs and other scenario data, allowing practitioners to assess systems under changing electricity mixes, fuel supply chains, technology efficiencies, and production pathways (Sacchi et al. 2022). Prospective LCA databases are temporal in the sense that each database represents a modeled future year. However, conventional prospective LCA typically evaluates one selected year at a time and does not by itself represent the timing of processes across a product life cycle.

Time-explicit LCA brings these two perspectives together: processes occur at different times, and the background systems supplying those processes also evolve over time. Recent frameworks such as bw_timex explicitly combine temporal distributions with time-specific background databases and dynamic characterization (Müller et al. 2025; Diepers et al. 2026). Arblaster et al. (2026) further show that different temporal perspectives can change the design insights obtained from LCA, while also emphasizing that the value of added temporal complexity depends on how the system and decision context are understood. They distinguish *imminent-future static*, *advanced-future static*, *mosaic time-explicit*, and *metabolic time-explicit* perspectives. This is an important framing for TRAILS: the process of considering temporal distribution in the entire LCA system, not only in the foreground model under study, is not an end in itself, but a way to test whether upstream temporal dynamics are decision-relevant within the temporal perspective selected for the study. There remains a practical need for workflows in which temporal distributions can be represented throughout the supply chain, including in background exchanges, and not only in the foreground model or in links from foreground processes to time-specific backgrounds.

We refer to this capability as *deep temporalization*, that is the ability to represent and propagate temporal distributions on exchanges throughout the supply chain, including exchanges far upstream from the functional unit. This is particularly relevant when prospective databases are not independent snapshots but connected time slices of an evolving technosphere. For example, infrastructure exchanges, fleet turnover, biomass growth, delayed end-of-life treatment, and technology lifetimes may occur in many background processes (Pinsonnault et al. 2014; Pigné et al. 2020). If these temporalized exchanges are ignored beyond the foreground, time-resolved results may miss dynamics that are relevant for long-lived systems and systems dependent on fast-evolving supply chains. And if these temporal distirbution are absent, we cannot realistically expect the practitionner to manually include such disitributions, which are often outside of the expertise area.

We use TRAILS to implement *deep temporalization* in time-explicit LCA. TRAILS represents scenario-specific foreground and background inventories as a connected three-dimensional technosphere and biosphere system indexed by year. It applies graph-based temporal routing to build time-indexed demands and then uses year-wise matrix solves to calculate inventories and impacts while avoiding the computing struggle associated with resolving a single large time-expanded technosphere matrix (we demonstrate why this can be an issue in the ESI). This graph-matrix hybrid approach is intended to be practical for applied LCA users while preserving a clear connection to matrix-based LCA.

The contribution of this paper is fourfold. First, we formalize *deep temporalization* as temporal routing over connected year-specific foreground and background matrices. Second, we describe the graph-matrix hybrid algorithm implemented in TRAILS, where explicit temporal routing is combined with conventional year-wise LCA solves. Third, we make explicit a modeling choice that is often implicit in time-explicit calculations: whether a temporally distributed exchange preserves the coefficient from the process year or is re-evaluated using the coefficient represented in the target year where the exchange occurs. Fourth, we demonstrate the method with case studies designed to test when deeper temporalization changes time-resolved LCA results.

The paper is structured as follows. Section 2 defines the scope and concepts used in the manuscript. Section 3 describes the method. Section 4 presents TRAILS as the implementation used in the case studies. Section 5 describes the demonstration design and case-study selection strategy. Section 6 provides the planned result structure for the selected case studies. Section 7 discusses interpretation, applicability, and limitations.

## 2. Concepts and scope

### 2.1 Static, dynamic, prospective, and time-explicit LCA

In conventional matrix-based LCA, a final demand vector is used to solve for the activity levels required to supply the functional unit (Heijungs and Suh 2002). The associated biosphere exchanges are then characterized using life cycle impact assessment factors. In a static calculation, the product system and its background economy are treated as fixed.

Dynamic LCA adds temporal differentiation. It can represent when processes, product flows, elementary flows, or impacts occur, and it can apply characterization models that depend on emission timing (Beloin-Saint-Pierre et al. 2020; Sohn et al. 2020). Prospective LCA, by contrast, represents how the product system changes over future years or scenarios, typically by modifying background databases or foreground parameters for a target year (Arvidsson et al. 2024). Time-explicit LCA combines both aspects: processes and emissions occur over time, and the system supplying them may also change over time (Müller et al. 2025). In this manuscript, the term refers to a computational approach, but we follow Arblaster et al. (2026) in treating temporal modeling choices as part of system framing rather than as a purely technical upgrade.

TRAILS is designed for time-explicit LCA in which the temporal resolution is annual. Temporal resolution is itself a modeling choice in dynamic LCA (Lueddeckens et al. 2020). TRAILS is not intended to model sub-annual variation such as seasonal electricity profiles, hourly emissions, or day-night differences in exposure. It is also deterministic: it does not sample uncertainty distributions for exchange amounts or model parameters.

### 2.2 Temporal perspectives and TRAILS

Arblaster et al. (2026) distinguish four temporal perspectives for prospective and time-explicit LCA. The first two are static perspectives. In imminent-future static LCA, the product system is modeled as a steady state close to current conditions, even if some life cycle stages are conceptually located in the future. In advanced-future static LCA, the product system is also modeled as a steady state, but the database and foreground assumptions represent a future that differs substantially from current conditions.

The other two perspectives are time-explicit. Mosaic time-explicit LCA follows a single object or cohort through time. Different life cycle stages are placed in different calendar years, so a product manufactured in one year may use energy, require maintenance, or reach end-of-life in later years. Metabolic time-explicit LCA instead represents a function provided by many objects or cohorts over a time span. Its demand vector is itself time-explicit, often informed by stock-flow modeling, market introduction, or transition dynamics.

TRAILS fits most directly with mosaic time-explicit LCA. A functional unit is anchored in a start year, temporal distributions place foreground and upstream exchanges in calendar time, and year-specific matrices represent the background system encountered at those times. The distinctive contribution of TRAILS is to make this mosaic perspective deep: temporalized exchanges can occur beyond the foreground, within connected prospective background systems.

TRAILS can also support metabolic time-explicit LCA when the functional unit is formulated as a time series of demands. In that case, a stock-flow model or scenario model can provide annual demands, and TRAILS can route and solve each demand through the connected scenario-year system. TRAILS does not decide whether a study should be imminent-future static, advanced-future static, mosaic, or metabolic; that remains a goal-and-scope choice. Its methodological role is to make deep temporalization available once a time-explicit perspective is justified.

### 2.3 Deep temporalization

We distinguish foreground-only temporalization from deep temporalization. Foreground-only temporalization assigns temporal distributions to exchanges or processes close to the functional unit and links those events to the appropriate time-specific background. Deep temporalization extends this principle to exchanges throughout the supply chain. Temporal distributions can therefore occur not only in the foreground inventory but also in background processes.

This distinction matters when background processes include time-dependent mechanisms. Examples include infrastructure construction before operation, delayed replacement of capital goods, age distributions of vehicle fleets, growth and harvest cycles for biomass, storage and release of biogenic carbon, delayed waste treatment, and changing production technologies across scenarios. In such cases, the temporal profile of upstream exchanges can affect time-resolved inventories and impacts.

### 2.4 Connected prospective databases

Prospective LCA databases are usually created as scenario-year snapshots, often to represent how background systems change under integrated assessment model scenarios (Mendoza Beltran et al. 2020; Sacchi et al. 2022). A scenario may contain databases for 2020, 2030, 2050, and 2100, with changes in technology efficiencies, energy mixes, and process availability. TRAILS treats these snapshots as slices of a connected three-dimensional matrix system. Temporal routing can then move from an activity in one calendar year to an upstream activity in another calendar year according to the temporal distribution on the exchange.

The approach does not require premise, but premise-generated prospective databases are an important supported pathway. In this workflow, temporal distributions can be added to exchanges across future ecoinvent-derived databases, for example to represent lifetimes, stock turnover, biomass-related delays, and infrastructure timing. The detailed parameterization of these temporal exchanges belongs in the electronic supplementary information.

## 3. Method

### 3.1 Static and prospective matrix notation

For each calendar year or scenario year \(y\), TRAILS uses technosphere and biosphere matrices in the same practical matrix-based data model used by Brightway workflows (Mutel 2017). Let \(A_y\) denote the technosphere matrix, \(B_y\) the biosphere or intervention matrix, \(d_y\) a demand vector in year \(y\), and \(x_y\) the corresponding supply or activity scaling vector. The year-specific matrix solve can be written as:

\[
x_y = A_y^{-1} d_y
\]

The year-specific inventory is then:

\[
g_y = B_y x_y
\]

or, when activity-resolved inventories are retained, as an activity-flow inventory matrix derived from \(B_y\) and \(x_y\). Characterization factors are represented by a vector \(q\), or by a diagonal matrix \(Q\), giving characterized scores:

\[
h_y = Q g_y
\]

In prospective LCA, \(A_y\) and \(B_y\) vary by year and scenario. TRAILS builds on this formulation but separates two tasks: explicit temporal routing of selected exchanges and year-specific matrix solving for unresolved frontier demands.

### 3.2 Scenario-year matrix tensors

TRAILS stores the technosphere and biosphere systems as sparse three-dimensional arrays:

\[
A[t, i, j]
\]

and

\[
B[t, i, k]
\]

where \(t\) indexes the scenario year, \(i\) indexes the activity, \(j\) indexes the technosphere product or linked activity, and \(k\) indexes the biosphere flow. A slice \(A_y\) or \(B_y\) is selected for the year \(y\), after mapping the requested calendar year to the available scenario-year grid.

Scenario-year matrices can be loaded directly or interpolated to an annual grid. In typical prospective applications, input databases may exist for selected years such as 2005, 2020, 2050, and 2100, reflecting the scenario-year structure common in prospective LCA workflows (Mendoza Beltran et al. 2020; Sacchi et al. 2022). TRAILS can linearly interpolate matrix values to annual resolution so that routing and solving can proceed on a yearly calendar. The main text treats interpolation as an implementation detail; full data schema and interpolation behavior are documented separately.

### 3.3 Temporal exchange distributions

Each technosphere or biosphere exchange can optionally carry a temporal distribution, following the process-relative temporal descriptions used in dynamic LCA (Beloin-Saint-Pierre et al. 2014; Cardellini et al. 2018). A temporal distribution defines a set of integer year offsets and weights relative to the year of the consuming or emitting process. If a process occurs in year \(y\), and an exchange has offset \(o\) with weight \(w_o\), then the corresponding exchange pulse is assigned to year \(y + o\).

TRAILS supports several distribution families, including discrete, lognormal, normal, uniform, triangular, and empirical distributions with explicit offsets and weights. Weights are normalized over the integer offset support. The annual time step is a deliberate scope choice: TRAILS does not represent events within a year.

### 3.4 Interpreting Temporally Distributed Exchange Amounts

TRAILS makes explicit a modeling choice about the amount attached to a temporally distributed exchange.

In the first interpretation, the exchange coefficient is read from the anchor-year matrix and then distributed over target years according to the temporal weights. This is appropriate when the temporal distribution describes the timing of an amount determined by the process as represented in the anchor year.

In the second interpretation, the temporal distribution determines the target years, but the exchange coefficient is read from the matrix corresponding to each target year. This is appropriate when the occurrence of the exchange is temporally distributed and the amount should reflect the technology or supply-chain conditions of the year in which the exchange occurs.

This distinction is important in prospective databases because matrix values can change over time. A delayed exchange routed to 2050 may represent a different technology mix or process efficiency than the same exchange in 2020. The target-year interpretation allows this information to affect the routed amount directly from the connected scenario-year matrices, while the anchor-year interpretation preserves the coefficient from the process year.

### 3.5 Temporal routing graph

TRAILS builds a temporal routing graph from the functional unit. Nodes represent process-year combinations:

\[
v = (i, y)
\]

where \(i\) is an activity and \(y\) is a calendar year. Edges represent technosphere exchanges from one process-year node to another. If an exchange has no temporal distribution, the child demand remains in the same year. If an exchange has a temporal distribution, one parent node may create several child nodes in different years.

Routing is controlled by a maximum depth and a minimum amount threshold. The maximum depth limits how far explicit temporal traversal proceeds into the supply chain. The minimum amount threshold avoids excessive graph expansion from negligible flows. Nodes at the boundary of explicit traversal form the frontier. Frontier demands are then solved using conventional year-wise matrix systems.

Algorithm 1 summarizes the routing and solving logic.

```text
Algorithm 1: Graph-matrix hybrid calculation in TRAILS

Input:
  functional unit activity i0
  start year y0
  amount a0
  scenario-year matrices A_y and B_y
  temporal exchange metadata
  maximum routing depth D
  cutoff epsilon

1. Initialize graph with root node (i0, y0, a0).
2. For each expandable node (i, y, a, depth):
     a. Read technosphere exchanges from A_y for activity i.
     b. For each linked product/activity j:
          i. If no temporal distribution exists, create child (j, y).
         ii. If a temporal distribution exists, create child nodes
             (j, y + offset) for all non-zero temporal pulses.
        iii. Compute child amounts using either the anchor-year or
             target-year amount interpretation.
     c. If depth reaches D or child amount is below epsilon,
        record the demand as frontier demand.
     d. Record direct biosphere exchanges from expanded nodes.
3. Aggregate frontier demands by year.
4. For each year with frontier demand:
     a. Select A_y and B_y.
     b. Solve the year-wise matrix system.
     c. Accumulate biosphere flows by emission year.
5. Characterize inventories or store time-resolved inventory arrays.
6. Optionally pass time-resolved inventories to dynamic impact models.
```

The graph can also retain attribution to first-tier root activities. Root attribution is useful for interpreting which immediate branch of the functional unit is responsible for impacts that occur later in time, but it is not required for the core calculation.

### 3.6 Year-wise solving and inventory accumulation

After temporal routing, TRAILS aggregates unresolved frontier demands by year. For each year, it solves the corresponding static LCA system using the year-specific technosphere matrix. The main methodological point is the separation between explicit temporal routing and year-wise solving; solver options and computational details are provided in the supplementary material.

Biosphere exchanges are accumulated on an inventory year axis. Biosphere exchanges can also have temporal distributions, so emissions may occur before or after the activity year. The resulting inventory can be stored as a sparse time-resolved array with dimensions for activity, biosphere flow, year, and optionally root activity.

Characterization can be applied during score accumulation or after storing the inventory. TRAILS includes characterization data for selected ecoinvent versions and can calculate time series of characterized scores. Because characterization matching depends on biosphere flow metadata, consistent flow names and compartments are important.

### 3.7 Optional dynamic impact assessment

The time-resolved inventory produced by TRAILS can be passed to dynamic impact models. This connection is important because fixed GWP-style metrics can misrepresent emission timing, and dynamic LCA has proposed time-dependent characterization and response calculations for temporally differentiated inventories (Levasseur et al. 2010; Kendall 2012; Shimako et al. 2018). Recent work also extends dynamic climate impact modeling to short-lived climate forcers, aviation and shipping emissions, and carbon-cycle climate feedbacks (Tiruta-Barna 2026). In the current implementation, TRAILS includes an optional interface to the FaIR climate emulator. This interface maps greenhouse gas flows to FaIR species, applies inventory perturbations to a baseline emissions scenario, and calculates delta radiative forcing and delta temperature time series.

The FaIR coupling is used here as a demonstration of how time-resolved inventories can feed dynamic impact assessment. It is not the primary methodological contribution of TRAILS. The central contribution remains deep temporalization of the life cycle inventory calculation.

## 4. TRAILS implementation for the case studies

TRAILS is the software implementation used to operationalize the method described above. In the context of this manuscript, its role is to make deep temporalization calculable for realistic prospective LCA systems: it connects year-specific matrices, applies temporal routing, performs year-wise LCA solves, and returns time-resolved inventories or impact scores.

The manuscript focuses on the methodological behavior of TRAILS and on the case-study results. Technical details of the implementation are provided in the electronic supplementary information and package documentation.

For the case studies, TRAILS will be used with prospective databases that represent selected years and scenarios. TRAILS interprets these databases as connected slices of an evolving technosphere, so that temporal distributions can route demands across years and across the supply chain. This allows the case studies to compare static, foreground-temporalized, and deeply temporalized perspectives using the same underlying scenario assumptions.

Reproducibility is a central design goal. The demonstration figures in this manuscript will be accompanied by Jupyter notebooks in the electronic supplementary information. The notebooks will document the modeled systems, temporalization settings, routing depths, impact indicators, and result processing steps needed to reproduce the figures.

## 5. Demonstration design

The purpose of the case studies is not to produce definitive comparative conclusions about specific technologies. Their purpose is to demonstrate when deep temporalization matters and how TRAILS can be used to identify such cases. The final manuscript will include four to five case studies selected through a screening procedure. The case-study structure will use the temporal perspectives of Arblaster et al. (2026) to clarify what each comparison represents.

Candidate systems should have one or more of the following characteristics, reflecting classes of temporal issues discussed in dynamic LCA and prospective background-scenario literature (Brandão et al. 2013; Beloin-Saint-Pierre et al. 2020; Lueddeckens et al. 2020; Mendoza Beltran et al. 2020; Sacchi et al. 2022):

- long service lives;
- infrastructure-heavy supply chains;
- important upstream capital goods;
- dependence on energy systems that evolve strongly over time;
- biomass growth, storage, or delayed release dynamics;
- carbon removal, avoided emissions, or temporary carbon storage;
- delayed end-of-life treatment or replacement cycles;
- large differences between foreground-only and upstream temporalization.

For each selected case, the manuscript will compare temporal perspectives selected from the following set:

1. imminent-future static LCA, where the system is assessed as a near-current steady state;
2. advanced-future static LCA, where the system is assessed as a future steady state;
3. foreground-only or shallow mosaic time-explicit LCA, where the functional unit follows a single object or cohort through time but upstream temporalization remains limited;
4. deep mosaic time-explicit LCA with TRAILS, where the same object or cohort is followed through time and temporal distributions can propagate into background processes;
5. metabolic time-explicit LCA with TRAILS, where relevant, for cases where the functional unit represents many objects or cohorts produced, used, or retired over a time span.

Not every case study needs all five perspectives. For most TRAILS demonstrations, the central comparison is expected to be between foreground-only mosaic time-explicit LCA and deep mosaic time-explicit LCA. Metabolic time-explicit LCA should be used only where the decision concerns deployment, market introduction, fleet turnover, circularity transitions, or other dynamics that cannot be represented by a single object or cohort.

Routing depth will be varied to show how results evolve as temporalization reaches further into the supply chain. This sensitivity directly tests the central claim of the paper: that foreground temporalization can be insufficient for some systems, and that upstream temporalized exchanges can change time-resolved scores within a mosaic or metabolic time-explicit perspective.

At least one case should also compare the two temporal amount interpretations. This comparison is needed because the distinction is methodological, not merely computational. If the coefficient attached to a delayed exchange is re-evaluated in the target year, prospective changes in technology can affect the routed exchange in a way that differs from simply distributing the anchor-year coefficient.

## 6. Results

This section will be completed after the case-study screening and calculations. It is currently structured around the evidence needed to support the manuscript's main claim.

### 6.1 Deep temporalization can change time-resolved scores

The main result figure should compare the relevant temporal perspectives for the selected four to five cases. For most cases, this means comparing imminent-future or advanced-future static LCA with shallow mosaic time-explicit LCA and deep mosaic time-explicit LCA. Where the decision concerns deployment over time, circularity transitions, or changing stocks, a metabolic time-explicit perspective should be added. The figure should show whether deeper temporalization changes:

- the magnitude of the time-resolved score;
- the timing of impacts;
- cumulative impacts over the assessment horizon;
- the ranking or interpretation of alternatives, where alternatives are included.

The strongest case-study result would be a clear divergence between shallow mosaic time-explicit LCA and deep mosaic time-explicit LCA, or between mosaic and metabolic perspectives where the functional unit spans multiple cohorts. Cases where these approaches agree are also useful, because they help define when the added complexity of TRAILS is unnecessary.

### 6.2 Routing depth reveals where temporal dynamics enter the supply chain

A second result should show how scores change as routing depth increases. If results stabilize quickly, foreground or shallow temporalization may be sufficient. If results continue changing at greater depths, this indicates that temporalized background exchanges are relevant.

The figure should report both the impact score and a measure of traversal behavior, such as frontier size, routed amount, or contribution by depth, if this can be shown clearly without turning the paper into a performance study.

### 6.3 Delayed exchange amounts can be interpreted differently

The two amount interpretations should be illustrated in a case where the choice changes the result. The text should explain why. For example, if the target-year matrix represents a cleaner electricity mix, lower process intensity, or different technology composition, then re-evaluating the exchange coefficient in the target year can produce a different routed exchange amount than distributing the anchor-year coefficient.

If no selected case shows a meaningful difference, this section should be moved to the supplementary material and described as a modeling option rather than a result.

### 6.4 Optional climate response demonstration

If retained, the FaIR demonstration should show how a time-resolved inventory from TRAILS can be converted into delta radiative forcing or delta temperature under a consistent climate scenario. This result should support the relevance of emission timing and scenario consistency, but it should not distract from the main result on deep temporalization.

## 7. Discussion

### 7.1 Interpretation of deep temporalization

TRAILS addresses a specific gap in time-explicit LCA: the ability to propagate temporalized exchanges beyond the foreground into connected prospective background databases. In the terminology of Arblaster et al. (2026), the typical TRAILS application is a deep mosaic time-explicit assessment: one object or cohort is followed through calendar time, but temporal routing is allowed to continue into upstream background processes. When the demand itself is formulated over a time span, TRAILS can contribute to a metabolic time-explicit assessment by routing each time-specific demand through the connected background system. The added value is not that TRAILS is universally more accurate than other tools, but that it exposes a dimension of temporal dynamics that can otherwise remain hidden.

The case studies should therefore be interpreted as demonstrations of sensitivity to deep temporalization within explicitly stated temporal perspectives. If deep temporalization changes a result, the analyst should investigate which upstream exchanges are responsible and whether their temporal distributions are well supported. If it does not change a result, a simpler imminent-future static, advanced-future static, or shallow mosaic perspective may be adequate for that system.

### 7.2 Relationship to existing tools

TRAILS complements existing dynamic and time-explicit LCA tools. Temporalis and related work established graph-based propagation of temporal distributions and dynamic characterization in the Brightway ecosystem (Cardellini et al. 2018). DyPLCA and related tools operationalized dynamic LCI with time differentiation on the background database (Tiruta-Barna et al. 2016; Pigné et al. 2020). bw_timex provides a modern open-source framework for time-explicit LCA that links temporally differentiated foreground events to time-specific background databases and embeds time into matrix-based calculations (Müller et al. 2025; Diepers et al. 2026).

The treatment of delayed exchange amounts is related to, but distinct from, the temporal evolution functionality in bw_timex. In bw_timex, users can specify time-dependent scaling factors or absolute amounts for foreground exchanges, and those user-defined modifiers adjust the base exchange amount at the time of the process. TRAILS instead exposes the coefficient choice as part of temporal routing over connected scenario-year matrices: a delayed exchange can preserve the coefficient from the process year or be re-evaluated from the matrix slice corresponding to the target year. The latter does not require the user to manually provide exchange-specific scaling factors when the relevant coefficient changes are already represented in the prospective database sequence.

TRAILS should be understood as an adjacent approach focused on connected three-dimensional technospheres and deep temporalization. It is especially aligned with prospective LCA workflows in which multiple year-specific databases already exist and can be interpreted as a connected system. This positioning avoids the claim that TRAILS supersedes existing tools; rather, it clarifies the class of temporal questions for which TRAILS is designed.

### 7.3 When TRAILS is useful

TRAILS is most useful for systems where temporal dynamics are expected to propagate beyond the foreground. In Arblaster et al.'s terms, this most often means moving from a shallow mosaic time-explicit perspective to a deep mosaic time-explicit perspective. Examples include long-lived infrastructure, transport systems with evolving fuel or electricity supply, carbon removal and storage systems, bio-based products, systems with delayed end-of-life treatment, and prospective studies where upstream technology mixes change rapidly.

TRAILS may not be worth the added complexity for short-lived products, systems dominated by foreground emissions, or studies where temporal distributions are poorly known and unlikely to affect the conclusion. In such cases, imminent-future static LCA, advanced-future static LCA, or shallow mosaic time-explicit LCA may be sufficient. Metabolic time-explicit use of TRAILS should be reserved for questions where the functional unit itself represents provision of a function by many objects or cohorts over time, such as deployment pathways, fleet transitions, or circular economy transitions.

### 7.4 Limitations

TRAILS is deterministic. It does not perform uncertainty sampling from exchange uncertainty fields, temporal distribution parameters, scenario assumptions, or characterization factors. Uncertainty analysis remains an important future development.

TRAILS uses annual time steps. It cannot represent seasonal, monthly, hourly, or event-scale dynamics. This makes it unsuitable for questions where sub-annual timing is central, such as hourly electricity balancing, seasonal air pollution exposure, or day-night noise impacts (Lueddeckens et al. 2020).

The method depends on temporal distribution assumptions. Deep temporalization is only as meaningful as the temporal metadata attached to exchanges. Poorly parameterized temporal distributions can create a false sense of precision.

TRAILS abstracts the calculation away from some native Brightway context. This makes the temporal calculation more streamlined, but it can reduce direct interaction with Brightway data structures during analysis. The tool is intended to complement Brightway workflows, not replace them.

The current software focuses on inventory temporalization and selected impact pathways. Climate response can be explored through FaIR, but TRAILS does not yet implement the broader dynamic climate impact modeling now available in specialized tools, such as explicit short-lived climate forcer modeling and carbon-cycle climate feedbacks (Tiruta-Barna 2026). Other time-dependent impacts, such as toxicity, air pollution exposure, and pollutant migration through environmental media, also require further development.

### 7.5 Future work

Future work should extend dynamic impact assessment beyond the current climate response interface. For climate change, this includes better coverage of short-lived climate forcers, aviation and shipping effects, and carbon-cycle climate feedback mechanisms highlighted by Tiruta-Barna (2026). Beyond climate change, relevant directions include time-dependent air pollution impacts, ecotoxicity, pollutant transport across media, delayed exposure of human populations and ecosystems, and other impact categories where timing affects fate, exposure, or effect (Shah and Ries 2009; Lebailly et al. 2014; Shimako et al. 2017).

Further methodological development should address uncertainty propagation, sub-annual temporal resolution where needed, richer temporal metadata generation, improved database linking, and additional diagnostics to help users identify which upstream temporalized exchanges drive results.

## 8. Conclusions

TRAILS provides a graph-matrix hybrid approach for deep temporalization in time-explicit LCA. By representing foreground and background inventories as connected scenario-year matrices, it allows temporal distributions to propagate throughout the supply chain rather than remaining limited to the foreground. The method combines explicit temporal routing with year-wise matrix solving, preserving compatibility with conventional matrix LCA while avoiding a single large time-expanded solve.

The main value of TRAILS is diagnostic and methodological: it helps practitioners test whether upstream temporalization changes time-resolved LCA results. For systems where deep temporalization matters, TRAILS can reveal dynamics that are missed by static or foreground-only approaches. For systems where it does not matter, TRAILS can help justify simpler modeling choices. In this sense, TRAILS complements existing dynamic and time-explicit LCA tools by making connected prospective databases available for deep temporal analysis.

## Data Availability

The distributable form of the Python library used to produce the results presented in this article is freely available from the Python Package Index (PyPI). Its source code is available at https://github.com/Laboratory-for-Energy-Systems-Analysis/trails. Documentation is available at https://trails.readthedocs.io/en/latest/. Jupyter notebooks to reproduce the figures in the results section will be provided as electronic supplementary information. A version of the code will be archived on Zenodo before publication.

## References

Arblaster T, Guinée J, Blanco Rocha CF, Burzic I, Pretschuh C, Pérez Sánchez F, Thonemann N (2026) System understanding shapes insights for eco-design: a comparison of four temporal perspectives. Research Square. https://doi.org/10.21203/rs.3.rs-9554797/v1

Arvidsson R, Svanström M, Sandén BA, Thonemann N, Steubing B, Cucurachi S (2024) Terminology for future-oriented life cycle assessment: Review and recommendations. Int J Life Cycle Assess 29:607-613. https://doi.org/10.1007/s11367-023-02265-8

Beloin-Saint-Pierre D, Albers A, Hélias A, Tiruta-Barna L, Fantke P, Levasseur A, Benetto E, Benoist A, Collet P (2020) Addressing temporal considerations in life cycle assessment. Sci Total Environ 743:140700. https://doi.org/10.1016/j.scitotenv.2020.140700

Beloin-Saint-Pierre D, Heijungs R, Blanc I (2014) The ESPA (Enhanced Structural Path Analysis) method: A solution to an implementation challenge for dynamic life cycle assessment studies. Int J Life Cycle Assess 19:861-871. https://doi.org/10.1007/s11367-014-0710-9

Brandão M, Levasseur A, Kirschbaum MUF, Weidema BP, Cowie AL, Jørgensen SV, Hauschild MZ, Pennington DW, Chomkhamsri K (2013) Key issues and options in accounting for carbon sequestration and temporary storage in life cycle assessment and carbon footprinting. Int J Life Cycle Assess 18:230-240. https://doi.org/10.1007/s11367-012-0451-6

Cardellini G, Mutel CL, Vial E, Muys B (2018) Temporalis, a generic method and tool for dynamic Life Cycle Assessment. Sci Total Environ 645:585-595. https://doi.org/10.1016/j.scitotenv.2018.07.044

Diepers T, Müller A, Jakobs A (2026) bw_timex: A Python Package for Time-Explicit Life Cycle Assessment. J Open Source Softw 11(120):9621. https://doi.org/10.21105/joss.09621

Heijungs R, Suh S (2002) The Computational Structure of Life Cycle Assessment. Kluwer Academic Publishers, Dordrecht. https://doi.org/10.1007/978-94-015-9900-9

Kendall A (2012) Time-adjusted global warming potentials for LCA and carbon footprints. Int J Life Cycle Assess 17:1042-1049. https://doi.org/10.1007/s11367-012-0436-5

Lang-Quantzendorff L, Beermann M (2025) Prosperdyn - a tool to describe dynamic transitions in prospective life cycle assessment. Int J Life Cycle Assess 30:3214-3229. https://doi.org/10.1007/s11367-025-02515-x

Lebailly F, Levasseur A, Samson R, Deschênes L (2014) Development of a dynamic LCA approach for the freshwater ecotoxicity impact of metals and application to a case study regarding zinc fertilization. Int J Life Cycle Assess 19:1745-1754. https://doi.org/10.1007/s11367-014-0779-1

Levasseur A, Lesage P, Margni M, Deschênes L, Samson R (2010) Considering time in LCA: dynamic LCA and its application to global warming impact assessments. Environ Sci Technol 44:3169-3174. https://doi.org/10.1021/es9030003

Lueddeckens S, Saling P, Guenther E (2020) Temporal issues in life cycle assessment - a systematic review. Int J Life Cycle Assess 25:1385-1401. https://doi.org/10.1007/s11367-020-01757-1

Mendoza Beltran A, Cox B, Mutel C, van Vuuren DP, Font Vivanco D, Deetman S, Edelenbosch OY, Guinée J, Tukker A (2020) When the background matters: using scenarios from integrated assessment models in prospective life cycle assessment. J Ind Ecol 24:64-79. https://doi.org/10.1111/jiec.12825

Müller A, Diepers T, Jakobs A, Cardellini G, von der Assen N, Guinée J, Steubing B (2025) Time-explicit life cycle assessment: a flexible framework for coherent consideration of temporal dynamics. Int J Life Cycle Assess 30:3052-3071. https://doi.org/10.1007/s11367-025-02539-3

Mutel C (2017) Brightway: An open source framework for Life Cycle Assessment. J Open Source Softw 2:236. https://doi.org/10.21105/joss.00236

Pigné Y, Gutiérrez TN, Gibon T, Schaubroeck T, Popovici E, Shimako AH, Benetto E, Tiruta-Barna L (2020) A tool to operationalize dynamic LCA, including time differentiation on the complete background database. Int J Life Cycle Assess 25:267-279. https://doi.org/10.1007/s11367-019-01696-6

Pinsonnault A, Lesage P, Levasseur A, Samson R (2014) Temporal differentiation of background systems in LCA: relevance of adding temporal information in LCI databases. Int J Life Cycle Assess 19:1843-1853. https://doi.org/10.1007/s11367-014-0783-5

Sacchi R, Terlouw T, Siala K, Dirnaichner A, Bauer C, Cox B, Mutel C, Daioglou V, Luderer G (2022) PRospective environMental impact asSEment (premise): A streamlined approach to producing databases for prospective life cycle assessment using integrated assessment models. Renew Sustain Energy Rev 160:112311. https://doi.org/10.1016/j.rser.2022.112311

Shah VP, Ries RJ (2009) A characterization model with spatial and temporal resolution for life cycle impact assessment of photochemical precursors in the United States. Int J Life Cycle Assess 14:313-327. https://doi.org/10.1007/s11367-009-0084-6

Shimako AH, Tiruta-Barna L, Ahmadi A (2017) Operational integration of time dependent toxicity impact category in dynamic LCA. Sci Total Environ 599-600:806-819. https://doi.org/10.1016/j.scitotenv.2017.04.211

Shimako AH, Tiruta-Barna L, Bisinella de Faria AB, Ahmadi A, Spérandio M (2018) Sensitivity analysis of temporal parameters in a dynamic LCA framework. Sci Total Environ 624:1250-1262. https://doi.org/10.1016/j.scitotenv.2017.12.220

Sohn J, Kalbar P, Goldstein B, Birkved M (2020) Defining temporally dynamic life cycle assessment: a review. Integr Environ Assess Manag 16:314-323. https://doi.org/10.1002/ieam.4235

Tiruta-Barna L, Pigné Y, Navarrete Gutiérrez T, Benetto E (2016) Framework and computational tool for the consideration of time dependency in Life Cycle Inventory: Proof of concept. J Clean Prod 116:198-206. https://doi.org/10.1016/j.jclepro.2015.12.049

Tiruta-Barna L (2026) Expanding the dynamic climate change impact model for dynamic LCA based on the climate science. Int J Life Cycle Assess 31:27. https://doi.org/10.1007/s11367-026-02583-7

## Acknowledgments

The authors acknowledge the financial support of the Europe Horizon project UPTAKE (project ID 101081521) via the contribution from Switzerland's State Secretariat for Education, Research and Innovation SERI.

## Competing Interests

**The authors declare no conflict of interest.**

## Supplementary Information

The electronic supplementary information will include Jupyter notebooks to reproduce the figures, extended case-study results, technical details of the TRAILS library and dependencies, the input schema, and documentation of the temporal exchanges added to premise-generated prospective databases for TRAILS.

## Contributions

Author contributions will be completed after finalizing the case studies and writing responsibilities.
