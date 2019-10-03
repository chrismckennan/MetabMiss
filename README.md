# MetabMiss
Estimate and account for non-ignorable missingness mechanisms and latent confounders in metabolomic and proteomic data. The corresponding paper can be found here: https://arxiv.org/abs/1909.02644

## Installation
First, make sure the packages "devtools", "parallel", "qvalue" and "sva" are installed. If not, the first two can be installed from CRAN and the second two can be installed from Bioconductor. To install MetabMiss, type the following into R:

devtools::install_github("chrismckennan/MetabMiss"); library(MetabMiss)

## Package functions
The package has three functions: "Num.Instruments", "EstimateMissing" and  "CC.Missing".


## Description and usage
The following instructions describe how to use these three functions to analyze metabolomic data. We recommend using this package on relatively large datasets (i.e. n > 250 independent samples).

1.) Log2 (log base 2) transform your raw MS integrated intensity data. Leave all missing values as NA and DO NOT USE IMPUTED VALUES. Call this datamatrix Y, which is a #metabolite x #sample datamatrix.
2.) Run Num.Instruments(Y = Y), which will estimate the number of potential instruments, K, to use in downstream functions. Experienced users choose the number of potential instruments using Frac1 and Frac2. It is recommended that all other uses set K = $K.recommended.
3.) Estimate each metabolite-dependent missingness mechanism using Miss.Mech <- EstimateMissing(Y = Y, K = K). This is the most time consuming function, and may take a few hours on a laptop for large datasets. Luckily, this only has to be done once per datamatrix Y. SAVE THE OUTPUT AS SOON AS THIS FINISHES.
4.) Suppose Disease\_1, ..., Disease\_d are the d>=1 variables that you are interested in. For example, these might be disease status, environmental phenotypes, etc. Let Cov = model.matrix(~Disease\_1 + ... + Disease\_d + NUISANCE COVARIATES, data=data.frame.of.covariates). NUISANCE COVARIATES are covariates that are observed but whose coefficients you are not interested in. These might include covariates like sex, age, observed diet variables, etc.
5.) Run CC.Missing as out.miss <- CC.Missing(Y = Y, X = Cov[,2:(d+1)], Z = Cov[,-(2:(d+1))], Miss.Mech = Miss.Mech). out.miss$p.titer is a #metabolites x d matrix that whose (g,j)th element is the p-value for the null hypothesis H_0: Disease\_j does no influence the concentration of metabolite g. out.miss$Beta.iter is a #metabolites x d matrix that whose (g,j)th element is the effect of Disease\_j on the concentration of metabolite g.

We have noticed that sometimes R will not load the .Rd files using the "?" command after installing and then loading the package. We found that a simple fix for this is to install the package, restart your R-session and then load the package.
