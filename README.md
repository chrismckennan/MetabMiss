# MetabMiss
Estimate and account for non-ignorable missingness mechanisms and latent confounders in metabolomic and proteomic data

## Installation
First, make sure the packages "devtools", "parallel", "quantreg", "sn", "pbivnorm" and "mvtnorm" are installed. If not, they can be installed from CRAN. To install MetabMiss, type the following into R:

devtools::install_github("chrismckennan/MetabMiss"); library(MetabMiss)

## Package functions
The package has two functions: "EstimateMissing" and "CC.Missing". "EstimateMissing" estimates the metabolite-dependent missingness mechanisms. This function only needs to be run once per dataset, and the output should be stored with the data. The function "CC.Missing" estimates latent confounding and does inference on the coefficients of interest in a multivariate linear model. For example, "CC.Missing" will estimate the effect of disease on each metabolite's expression, as well as compute p-values for the null hypothesis that there is no effect of disease on a particular metabolite.  

## Description and usage
We have noticed that sometimes R will not load the .Rd files using the "?" command after installing and then loading the package. We found that a simple fix for this is to install the package, restart your R-session and then load the package.
