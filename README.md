# Identify single-cell variance eQTL in the OneK1K cohort

This repository contains the analysis code pipeline to identify single-cell variance eQTL (sc-veQTL) as part of the manuscript "**Genetic variants associated with cell-type-specific intra-individual gene expression variability reveal new mechanisms of genome regulation**"

In this study, we identified a novel genetic regulation mechanism, where genetic variants affect the variance and dispersion of the intra-individual gene expression independent from mean effects. This list of genes (dGene, i.e., gene with dispersion-eQTL) is enriched in the immune response and interspecies interaction and their cellular dispersion levels are associated with auto-immune disease risk.

Scripts are listed by the order in the methods section of the manuscript:
1. The collection and QC of the OneK1K cohort
2. Generate intra-individual mean, variance, and dispersion matrix for the gene expression
3. The sc- eQTL/veQTL/deQTL association
4. Simulation to evaluate the accuracy of dispersion estimation
5. Pseudotime inference
6. The G x G and G x E interaction test
7. Tras- eQTL and veQTL mapping
8. Replication analysis in non-EUR (East Asian) cohort

The repository will be updated soon after peer review.

# Citation

Angli Xue, Seyhan Yazar, José Alquicira-Hernández, Anna S E Cuomo, Anne Senabouth, Gracie Gordon, Pooja Kathail, Chun Jimme Ye, Alex W. Hewitt, Joseph E. Powell. Genetic variants associated with cell-type-specific intra-individual gene expression variability reveal new mechanisms of genome regulation. _Under Review_. 2024. ([Preprint](https://www.biorxiv.org/content/10.1101/2024.05.05.592598v1))

For questions, please email us at Angli Xue (a.xue@garvan.org.au) or Joseph E. Powell (j.powell@garvan.org.au)
