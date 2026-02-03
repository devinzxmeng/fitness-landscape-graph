## Our responsibility
- Figure 2
- Figure S4
- Body text section: Fitness graphs reveal the impact of selection type and strength on the evolution of ESBLs. 
- Methods section: Fitness Landscape Graph Construction
## Reviewer #1
### Minor comments:
Figure 2B,C and S4: It would be helpful to highlight the wild-type containing node to judge accessible peaks from this node.

Response: 
When I replot the landscape figures in Gephi, I need to use highlight color for the node containing the wildtype.
## Reviewer #3
### Major

2a) How do trajectories/predictability change from one MIC step to the next i.e. can you compare your graphs between concentration steps? What is the path for a sequence to each higher resistance step, not optimising within a given concentration? This would be a very interesting problem to be able to address.

Interpretation: Instead of looking at evolution *within* a single concentration, can you trace how a genotype would evolve *across* increasing concentrations? E.g., a genotype optimized at 0.5 µg/mL — where does it sit in the 1 µg/mL landscape, and what path does it take to the new peak? Essentially, they want a cross-concentration evolutionary trajectory analysis.

## Reviewer #4
iii) The graph analysis reports a striking non-monotonic dependence of the aztreonam landscape topology on drug concentration and states that at intermediate concentrations a global peak “disappears as it becomes absorbed into a large connection node,” thereby losing influence on evolutionary trajectories. This is counter-intuitive and currently under-explained. The authors should clarify the biological meaning of this “absorption.” In particular, it is important to disambiguate whether the effect reflects fitness differences between the former peak and surrounding genotypes falling below the neutrality threshold used in coarse-graining, whether it reflects genuine biological flattening due to concentration-dependent assay dynamics, or whether it is sensitive to the choice of neutrality definition (which depends on a threshold set from the no-drug condition). A useful addition would be a quantitative diagnostic showing how the peak’s fitness advantage over its neighborhood changes with concentration and a robustness check demonstrating whether the peak’s disappearance persists under reasonable perturbations of the neutrality cutoff.

Interpretation: The reviewer wants two things:
1. A quantitative plot showing how the peak's fitness advantage over its neighbors changes with concentration (to distinguish real biological flattening from a threshold artifact).
2. A robustness/sensitivity analysis showing the peak disappearance holds under different neutrality cutoff values.

Response:
To make the fitness landscape graph more robust to the choice of neutral threshold, we introduce the concept of forbidden pairs. This concept come from the observation of gradually increasing the neutral threshold: initially there are some big supernodes (which represent the macroscopic structure of the landscape) but there are too many nodes to be readable, then as we increase the neutral threshold, these macroscopic landscape structure start to be merged into each other and the number of nodes decrease. We noticed that it's really hard to pick a neutral threshold parameter such that it keeps the macroscopic structure but still have readable number of nodes, and we want the algorithm to be simple, interpretable, and fast. Therefore, we design this feature of forbidden pair, which use a "tiny" neutral merge to find the macroscopic landscape and then froze it from being merged. Please check the methods section or our code for more details.