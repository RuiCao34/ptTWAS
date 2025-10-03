# pt-TWAS
![](/Fig/FlowChart.png)

pt-TWAS is a novel pseudotime-dependent TWAS framework that incorporates single-cell transcriptomic data and disease GWAS data to infer causal genes at a cell-stage resolution.  Statistically, pt-TWAS is a two-stage functional regression framework. In the first stage, we use a function-on-scalar regression model to predict the single-cell gene expression trajectory using genetic variants, with a group-lasso penalty to impose SNP-wise sparsity. In the second stage, we test for an association between this imputed functional gene expression and a disease outcome. 

For inference, we are interested in both a global null hypothesis and a secondary hypothesis, correponding to a global null effect and a pseudotime-invariant effect respectively. In addition, to infer the causal cell stage for disease, we also visualize simultaneous confidence band for the gene effect curve. As demonstrated in both our simulation and real data results, by modeling gene expression as a function of pseudotime, our method identifies dynamic and cell-stage-specific genetic effects that can be missed by bulk-tissue or pseudobulk approaches. This functional approach not only provides a finer-resolution map of when a gene's activity contributes to disease but can also increases statistical power by borrowing information across the entire cellular continuum. 

# Installation
devtools::install_github("RuiCao34/ptTWAS")

# Vignette
To learn more about the package details, please check the vignette by calling the following code

vignette("vignette", package = "ptTWAS")
