# Rhcpp

A fast R implementation using Rcpp based on the original matlab HCP
method for Normalizing RNA-Sequencing Data by Modeling Hidden
Covariates with Prior Knowledge by Mostafavi et al. PLOSone (2013).
 
# Install

From unix and alike use either:

```r
library(devtools)
install_github("mvaniterson/Rhcpp")
```
Or using the `biocLite` function from BioConductor:

```r
source("http://bioconductor.org/biocLite.R")
biocLite("mvaniterson/Rhcpp")
```

# Usage

see ?hcp or the unitTest directory for examples. The package contains
example data from the original paper. 


