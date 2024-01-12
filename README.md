TAPhelpR is a companion R packge for the [TAP pipeline](https://github.com/FrietzeLabUVM/TAP). 

# Installation

```{r}
if(!require("devtools")){
  install.packages("devtools")
}
# This optional dependency has fallen out of Bioconductor and must be installed manually. 
# It is only required to download data from ENCODE.
devtools::install_github("CharlesJB/ENCODExplorer")
devtools::install_github("FrietzeLabUVM/TAPhelpR")
```

# Docs

For documentation specific to TAPhelpR see [this vignette](https://frietzelabuvm.github.io/TAPhelpR/) and ?function help within R.

See here for a [complete workflow use case](https://frietzelabuvm.github.io/TAP/).

