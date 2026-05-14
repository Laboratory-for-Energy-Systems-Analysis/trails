# Manuscript Framing Questions

Use this file as an answer sheet. Short answers are enough where the direction is
obvious; longer answers are useful where the manuscript needs a clear decision.

## A. Core Message

1. What is the one-sentence central claim of the paper?

   Answer: this new tool allows easy but powerful time-explicit LCA, combining prospective LCA and dynamic LCA and deep temporalization, allowing to assessment the performance of long-lived systems, or systems anchored in fast-evolving larger systems.

2. What should a reader remember about TRAILS one week after reading the paper?

   Answer: easy-to-use, fast, graph-based, time-explicit LCA tool, combining prospective and dynamic LCA, with coupling to climate emulator.

3. Is this paper primarily a methods paper, a software paper, an application
   paper, or a combined methods-and-software paper?

   Answer: combined methods-and-software paper, with applications as demonstrations. We will present the method, how it's contained in the tool, position teh method with respect to existing  methods/tools, and the edges of the tool, which is the effortless integration of deep temporalization (meaning temporal links between time-specific databases), by showing a few examples.

4. What is the specific problem in dynamic/prospective LCA that TRAILS solves
   better than existing approaches?

   Answer: state-of-the-art tools like bw_times do combine prospective LCA with dynamic LCA, but they do not do deep temporalization. That is, they temporalize the foregound inventories which call on time-specific LCA databases, but those LCA databases do not connect to one another. Trails consideres the foreground and nacground inventories as a single three-dimensional technosphere matrix, allowing temporalization of exchange upstream teh supply chain, and potentially far from the foreground.

5. What should the manuscript explicitly avoid claiming?

   Answer: it should avoid claiming that it is better than other tools like bw_times, and it should avoid claiming that it is first tool to do dynamic + prospective LCA (i.e., time-explicit LCA).

## B. Audience and Venue

6. Who is the primary audience: LCA method developers, applied LCA
   practitioners, prospective LCA users, energy systems researchers, climate
   impact modelers, or another group?

   Answer: applied LCA practitioners, prospective LCA users interested in long-lived systems, or the effect of biogenic CO2 storage, CDR, etc.

7. What journal or venue are we aiming for, and what does that imply for tone,
   length, methods detail, and software emphasis?

   Answer: Journal of Industrial Ecology, which implies a methods-and-application paper with moderate length (e.g., 5000-7000 words), enough methods detail to understand the approach but not so much as to be overwhelming, and a strong emphasis on the software implementation and reproducibility.

8. Should the manuscript assume readers know Brightway, premise, and matrix LCA,
   or should these be explained from first principles?

   Answer: We do not want to assume that all readers are familiar with Brightway, premise, and matrix LCA, but we also do not want to spend too much time explaining these concepts in detail. A brief overview of these concepts should be provided in the introduction or methods section, with references to more detailed explanations for readers who are interested.

9. Should the manuscript foreground open-source implementation and
   reproducibility, or keep those as supporting material?

   Answer: The manuscript should foreground open-source implementation and reproducibility as key strengths of the work. Jupyter notebook to reproduce teh figures in teh results sections (where we will demonstrate how the LCA scores of some systems evolve as you temporalize deeper in teh supply chain) will be provided as ESI.

10. What kind of reviewer criticism do you most want the manuscript to preempt?

    Answer: That the tool do not bring anything new. That is why, the main claim will be tha the tool brings the benefits of considering deep temporalization, whihc other tools do not. And exmaples will show that, for some systems, their LCA score will change as we temporal deeper int eh graph system. That is, the temporalization of the foreground inventories is not enough to capture the full effect of temporal dynamics on the LCA results, and that temporalization of upstream inventories can lead to different conclusions about the performance of a system over time.

## C. Positioning in the Literature

11. Which existing dynamic LCA methods should TRAILS be compared against most
    directly?

    Answer: those described in literature.md. The closest being bw_timex, and DyLCA. 

12. Which existing tools or papers should be treated as conceptual predecessors
    rather than competitors?

    Answer: bw_timex, and DyLCA. But hte real inspiration is Premise, which produces prospective LCA databases (for a given year and scenario). So, Trails should eb framed as a natural extension of the idea of prospective LCA databases, which are already temporalized in the sense that they are specific to a given year, but which do not connect to one another. Trails allows to connect those databases together, and to temporalize exchanges across them, allowing to capture the effect of temporal dynamics on the LCA results.

13. How should the manuscript distinguish TRAILS from Temporalis-style dynamic
    LCA?

    Answer: Temporalis-style dynamic LCA is a convolutional-family approach that represents temporal distributions as process-relative profiles and propagates timing through supply chains using graph traversal and convolution-based products of temporal information. TRAILS, on the other hand, is a graph-matrix hybrid approach that uses temporal routing to build time-indexed demands and year-wise solves to avoid a single massive solve. While both approaches aim to capture temporal dynamics in LCA, TRAILS allows for deep temporalization, meaning that temporal distributions can occur throughout the supply chain, not just in the foreground.

14. How should the manuscript distinguish TRAILS from time-explicit background
    linking frameworks such as the Muller et al. approach?

    Answer: The Muller et al. approach (i.e., bw_timex) is a time-explicit foreground-linking framework that connects temporally differentiated foreground demands to time-specific background databases (via constructs such as “temporal markets”) and interpolation/mapping across database years. TRAILS, on the other hand, is a graph-matrix hybrid approach that uses temporal routing to build time-indexed demands and year-wise solves to avoid a single massive solve. While both approaches aim to capture temporal dynamics in LCA, TRAILS allows for deep temporalization, meaning that temporal distributions can occur throughout the supply chain, not just in the foreground. Additionally, TRAILS is designed to ingest Premise-exported packages where many exchanges are temporally distributed across scenario years, effectively forming a connected 3D technosphere that TRAILS can traverse and solve.

15. Which literature gap should the introduction build toward: deep temporal
    traversal, prospective scenario consistency, computational scalability,
    reproducibility, or something else?

    Answer: The introduction should build toward the gap of deep temporal traversal, which is the ability to capture temporal dynamics throughout the supply chain, not just in the foreground. This is a key contribution of TRAILS that distinguishes it from existing approaches and allows for more accurate and nuanced assessments of long-lived systems and systems anchored in fast-evolving larger systems.

## D. Novel Contributions

16. What are the top three original contributions of the paper, in priority
    order?

    Answer: 1. Deep temporalization: temporal distributions can occur throughout the supply chain, not only in the foreground. 2. Hybrid algorithm: temporal routing to build time-indexed demands + year-wise solves to avoid a single massive solve. 3. Two temporal-amount semantics: after distributing an exchange over time (i.e., after routing), TRAILS can apply ported amounts or matrix-sourced amounts to each target year. 4. Coupling to FaIR climate emulator, allowing to produce scenario-consistent time series of climate impacts over time.

17. Is "graph-matrix hybrid" the best label for the method, or should we use a
    different framing?

    Answer: "Graph-matrix hybrid" is a descriptive label that captures the essence of the method, which combines graph-based temporal routing with matrix-based year-wise solving. However, it may not be the most intuitive or catchy label for a broader audience. We could consider "deep temporalization" as a more intuitive label that emphasizes the key contribution of the method, which is the ability to capture temporal dynamics throughout the supply chain. Alternatively, we could use a more descriptive label such as "temporal routing and year-wise solving" to convey the specific approach used in TRAILS.

18. How important is "deep temporalization" as a contribution, and how should it
    be defined precisely?

    Answer: "Deep temporalization" is a crucial contribution of the paper, as it allows for the capture of temporal dynamics throughout the supply chain, not just in the foreground. It should be defined precisely as the ability to represent and propagate temporal distributions on exchanges throughout the entire supply chain, including background processes, rather than being limited to the foreground inventories. This means that temporal dynamics can affect not only the immediate processes called by the foreground but also upstream processes that may have significant impacts on the overall LCA results.

19. Should the port-vs-matrix temporal amount semantics be presented as a major
    methodological contribution or as an implementation detail?

    Answer: The port-vs-matrix temporal amount semantics should be presented as a major methodological contribution, as it provides users with flexibility in how they interpret and apply temporal distributions in their LCA analyses. This distinction allows for different approaches to handling temporal dynamics, which can have significant implications for the results and conclusions of an LCA study. By highlighting this as a key contribution, we can emphasize the versatility and adaptability of TRAILS in accommodating various user preferences and analytical needs.

20. Should the FaIR climate response integration be part of the core
    contribution, a demonstration, or a secondary extension?

    Answer: The FaIR climate response integration should be presented as a demonstration of the capabilities of TRAILS, rather than a core contribution. While it is an important feature that allows for scenario-consistent time series of climate impacts over time, it is not the central methodological innovation of the paper. Instead, it serves as an example of how TRAILS can be applied to produce meaningful insights in the context of climate change impacts, showcasing the practical utility of the tool in real-world applications. Also, the integration of background concentration scenarios in FAiR from scenarios that are also available in Premise allows consistency between the prospective LCA databases, and teh climate response, which is a key strength of the tool, but it is not the main methodological contribution, which is the deep temporalization and the graph-matrix hybrid approach.

## E. Scope and Boundaries

21. What should be inside the main manuscript, and what should be moved to
    supplementary information?

    Answer: The main manuscript should include a clear and concise introduction to the problem, a detailed description of the TRAILS method, key results from the case studies, and a discussion of the implications and limitations of the work. The supplementary information can include additional technical details about the implementation, extended results from sensitivity analyses, and any additional case studies that support the main findings but are not essential for understanding the core contributions of the paper. The ESI should also include how Premise adds temporal exchanges in teh data package it exports for Trails.

22. Which capabilities of TRAILS should not be discussed because they distract
    from the paper's main story?

    Answer: n/a

23. Which limitations should be stated prominently rather than hidden in the
    discussion?

    Answer:
        - TRAILS is deterministic and does not perform uncertainty sampling, which means that it does not capture the full range of possible outcomes or the uncertainty associated with input data and model parameters. This should be stated prominently to clarify the scope of the tool and manage user expectations.
        - TRAILS works at an annual resolution and does not represent sub-annual dynamics (unlike bw_timex), which means that it may not capture important temporal variations that occur within a year, such as seasonal effects or short-term fluctuations in technology performance or emissions. This should also be stated prominently to clarify the limitations of the tool and guide users in selecting appropriate applications.

24. Should the paper emphasize that TRAILS is deterministic and does not perform
    uncertainty sampling?

    Answer: yes

25. Should the paper emphasize that TRAILS uses annual time steps and does not
    represent sub-annual dynamics?

    Answer: yes

## F. Method Description

26. How mathematical should the method section be: equations and notation,
    algorithm boxes, descriptive text, or a combination?

    Answer: a combination of equations for the LCA matrix part, and algorithm boxes for the routing part, with descriptive text to explain the intuition and the steps of the method. The equations should be clear and concise, and the algorithm boxes should be well-structured and easy to follow. The descriptive text should provide context and explanations for the equations and algorithms, making it accessible to readers who may not be familiar with the mathematical details.

27. What notation should we use for the scenario-year technosphere and
    biosphere matrices?

    Answer: A for the technosphere matrix, B for the intervention matrix, f for teh demand vector, and s for the scaling vector. f=A^-1 * s, where A is the technosphere matrix, s is the scaling vector, and f is the demand vector. s is diagonalized and multipled by B to get G, the inventory matrix. q, a vector contianing characterization factors is diagonalized into Q and multiplied by G to get the characterization matrix H. For the scenario-year technosphere and biosphere matrices, we can use A_y and B_y, where y indicates the year. The demand vector for year y can be denoted as f_y, and the scaling vector for year y can be denoted as s_y. The inventory matrix for year y can be denoted as G_y, and the characterization matrix for year y can be denoted as H_y.

28. How should the temporal routing graph be defined: nodes, edges, weights,
    years, roots, and frontier?

    Answer: The temporal routing graph can be defined as follows:
- Nodes: Each node represents a process-year combination, denoted as (process, year).
- Edges: Each edge represents a technosphere exchange between two nodes, with an optional temporal distribution that indicates how the exchange is distributed over time.
- Weights: The weights on the edges can represent the amount of the exchange, which can be distributed over time according to the temporal distribution.

29. How should we explain the relationship between temporal routing and
    year-wise matrix solving?

    Answer: Temporal routing is the process of traversing the graph of process-year nodes and distributing the demand from the foreground processes to the background processes according to the temporal distributions on the edges. This results in time-indexed demand vectors for each year. Year-wise matrix solving is then performed for each year using the corresponding technosphere and biosphere matrices for that year, allowing us to compute the inventory and impact results for each year separately. The combination of temporal routing and year-wise solving allows us to capture the temporal dynamics of the system while avoiding the computational complexity of solving a single massive system that includes all years at once.

30. Should root attribution be a central concept in the paper, or only an
    optional analysis feature?

    Answer: Root attribution should be presented as an optional analysis feature rather than a central concept in the paper. While it can provide valuable insights into the contributions of different processes to the overall impacts, it is not essential for understanding the core method of TRAILS. By framing it as an optional feature, we can highlight its usefulness without making it a requirement for using the tool or understanding the main contributions of the paper.

## G. Data Model and Implementation

31. How much detail should the manuscript give about the Frictionless
    datapackage schema?

    Answer: little detail. The manuscript should provide a high-level overview of the Frictionless datapackage schema, emphasizing its role in ensuring reproducibility and standardization of inputs for TRAILS. However, the specific details of the schema, such as the required A/B columns and the scenario-year metadata structure, can be provided in supplementary material or documentation for users who want to implement their own data packages.

32. Should premise integration be described as the standard input pathway, or
    only as one supported pathway?

    Answer: it's a supported pathway. Data package cna be hand-made or produced by any other tool. However, when produced by Premise, Premise will add temporal disitrbution throughout each year-specific LCA database teh data package contains. Much of the work in developing Trails was spent on parametrizing temporal disitrbutions for Premise to add. Premise add temporal distirbution to exchanges through the futurozed ecoinvent, definintion fleet age, lifetime, etc. for throughput processes, biomass growth-related processes, infrastructure processes, etc. The ESI will document these temporal exchanges.

33. Should the paper explain annual interpolation of scenario matrices in the
    main text, and if yes, at what level of detail?

    Answer: not necessarily. we will just mention that the technosphere matrix is three dimensions, with teh additional dimension being the years, interpolated, typically from 2005 to 2100 (but cna be expanded further).

34. Should solver choices and computational performance be described in the main
    method section or mainly in supplementary material?

    Answer: not necessarily. we use a specific iteration-based solver using a Jacobi preconditioner, followed by a Krylov solver, which is implemented in SciPy. We can mention this in the main text, but the details of the solver implementation and performance can be provided in supplementary material for readers who are interested in the computational aspects of the tool.

35. What reproducibility artifacts should be promised: code archive, example
    datapackage, notebooks, scripts, generated results, or all of these?

    Answer: Jupyter notebooks to reproduce the figures.

## H. Case Studies

36. How many case studies should the paper include?

    Answer: Five case studies highlighting the importance of deep temporalization.

37. What is the main question answered by the passenger car case study?

    Answer: I do not know yet whether we will keep the passenger car study. we need to figure otu the cases that are most sensitive to deep temporalization. we have not done this yet.

38. What is the main question answered by the carbon removal versus carbon
    avoidance case study?

    Answer: not sure if we will keep this case.

39. Are the case studies meant to validate TRAILS, demonstrate capabilities, or
    produce substantive LCA conclusions?

    Answer: they are ment to show that deep temporalization can matter for the LCA results of some systems, and to show how to use the tool to capture those effects. They are not meant to produce substantive LCA conclusions about specific systems, but rather to demonstrate the capabilities of the tool and the importance of deep temporalization in certain contexts.

40. What would count as a strong, publishable result from each case study?

    Answer: A strong, publishable result from each case study would be a clear demonstration that deep temporalization can lead to different LCA results compared to a non-temporalized or foreground-only temporalized approach. For example, in the passenger car case study, we might find that the LCA score of the car changes significantly over time when we consider temporal distributions throughout the supply chain, which could have implications for how we assess the environmental performance of the car over its lifetime. In the carbon removal versus carbon avoidance case study, we might find that the timing of emissions and removals leads to different conclusions about the climate impacts of these two approaches when we consider deep temporalization. Overall, a strong result would be one that clearly illustrates the added value of TRAILS in capturing temporal dynamics that other approaches might miss.

## I. Results and Figures

41. What are the essential figures the manuscript must contain?

    Answer: figures that show, for 4-5 cases, how the LCA score evolves as we go deeper in the deep temporalization (routing), versus only foreground temporalization (like bw_times does), vs. a static LCA score.

42. Should there be a conceptual workflow figure showing datapackage input,
    routing, frontier solving, inventory aggregation, LCIA, and FaIR?

    Answer: not necessarily.

43. Should there be a figure comparing static, prospective-static, and dynamic
    TRAILS results?

    Answer: yes, this is the main figure that will show the importance of deep temporalization. We can show how the LCA score of a system changes over time as we go from a static approach (no temporalization), to a prospective-static approach (temporalization only in the foreground), to a dynamic approach (deep temporalization throughout the supply chain). This figure will illustrate the added value of TRAILS in capturing temporal dynamics that other approaches might miss.

44. Should there be a performance or scaling figure showing runtime, depth,
    frontier size, or solver behavior?

    Answer: no.

45. Which sensitivity analyses are necessary: routing depth, temporal
    distributions, interpolation, solver tolerance, LCIA method, or scenario?

    Answer: routing depth and temporal distributions are the most important sensitivity analyses to include, as they directly relate to the core contributions of the paper. We can show how the LCA results change as we vary the routing depth, demonstrating the importance of deep temporalization. We can also show how different assumptions about temporal distributions affect the results, highlighting the flexibility of TRAILS in accommodating different user preferences and analytical needs. Other sensitivity analyses, such as interpolation, solver tolerance, LCIA method, or scenario, can be included in supplementary material if they provide additional insights but are not essential for understanding the main contributions of the paper.

## J. Discussion and Framing

46. When is TRAILS the right tool, and when is it not worth the added
    complexity?

    Answer: TRAILS is the right tool for assessing long-lived systems, infrastructure-dominated supply chains, and scenario-driven prospective LCA studies where temporal dynamics are expected to play a significant role in the results. It is particularly useful when the timing of emissions and resource use can affect the overall environmental performance of a system, such as in the case of carbon removal technologies or long-lived infrastructure projects. However, TRAILS may not be worth the added complexity for short-lived systems or systems with concentrated impacts in the foreground, where temporal dynamics may have a limited effect on the overall results. In such cases, a simpler static or foreground-temporalized approach may be sufficient to capture the key insights without the need for deep temporalization.

47. What should the discussion say about the loss or preservation of Brightway
    context?

    Answer: The discussion should acknowledge that while TRAILS is built on top of the Brightway framework and can ingest data from Brightway databases, it operates in a way that abstracts away some of the specific features and context of Brightway. This means that users may lose some of the direct interaction with Brightway's data structures and functionalities when using TRAILS. However, this abstraction allows TRAILS to focus on its core contribution of deep temporalization and time-explicit LCA, providing a more streamlined and user-friendly experience for those who may not be familiar with the intricacies of Brightway. The discussion should emphasize that TRAILS is designed to complement Brightway rather than replace it, and that users can still leverage their existing Brightway knowledge and resources when using TRAILS.

48. What future work should be described: uncertainty, sub-annual dynamics,
    better database linking, more impact categories, faster solvers, or other
    directions?

    Answer: Future developments will be to dynamically assessment other areas of environmental concerns, besides teh current coupling with FAiR, like air pollution, ecotoxicity matters, where time plays a role (i.e., migration of pollutants through different mediums to reach/impact human populations and ecosystems). 

49. What title would best match the rewritten manuscript?

    Answer: "A graph-matrix hybrid approach for deep temporalization in time-explicit LCA"

50. Which papers, concepts, or software projects absolutely must be cited for
    the manuscript to feel properly grounded?

    Answer: those contained in literature.md, with a focus on those that are most closely related to the core contributions of the paper, such as bw_timex, DyLCA, and Premise. Additionally, we should cite foundational works in LCA and dynamic LCA to provide context for our contributions. We should also cite any relevant works on climate emulators, such as FAiR, to contextualize our integration of climate response modeling into TRAILS. Finally, we should ensure that we cite any relevant software projects that we build upon or integrate with, such as Brightway and SciPy, to acknowledge the tools that enable our work and to guide readers who may want to explore those resources further.

