# Identify single-cell variance eQTL in the OneK1K cohort

This repository contains the analysis code pipeline to identify single-cell variance eQTL (sc-veQTL) as part of the manuscript "**Genetic variants associated with within-individual gene expression dispersion at single-cell resolution reveal new mechanisms of genome regulation**"

In this study, we identified a novel genetic regulation mechanism, where genetic variants affect the variance and dispersion of the intra-individual gene expression independent from mean effects. This list of genes (dGene, i.e., gene with dispersion-eQTL) is enriched in the immune response and interspecies interaction and their cellular dispersion levels are associated with auto-immune disease risk.

Scripts are listed by the order in the methods section of the manuscript:
1. The collection and QC of the OneK1K cohort
2. Generate intra-individual mean, variance, and dispersion matrix for the gene expression
3. The sc- eQTL/veQTL/deQTL association
4. Functional annotation and gene sets enrichment
5. The G x G and G x E interaction test
6. Integration of dGene with auto-immune disease risk

The repository will be updated soon upon manuscript submission later this year.

# Citation

Angli Xue, Seyhan Yazar, José Alquicira-Hernández, Anna S E Cuomo, Anne Senabouth, Gracie Gordon, Pooja Kathail, Chun Jimme Ye, Alex W. Hewitt, Joseph E. Powell. Genetic variants associated with within-individual gene expression dispersion at single-cell resolution reveal new mechanisms of genome regulation. _Preprint coming soon_. 2024.

For questions, please email us at Angli Xue (a.xue@garvan.org.au) or Joseph E. Powell (j.powell@garvan.org.au)
