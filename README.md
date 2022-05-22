
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bladdr: A Helping Hand in McConkey Land <img src='man/figures/logo.png' align="right" height="139" />

bladdr is a ‘miscellaneous’ toolbox made for the McConkey lab, but
useful to many who do routine lab tasks. Its vision is to be a package
that can automate everyday tasks in small functions that don’t try to do
too much. Larger, more full fledged ideas are pulled out into their own
packages, like `amplify` and `gp`.

# Overview

Currently, bladdr has two ‘groups’ of functions

-   EDA functions, which take a `SummarizedExperiment` object to
    visualize expression of gene(s) across optional strata

-   Helpers, like `dilute`, `make_pipette_vol`, `theme_tufte`, and
    others.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("McConkeyLab/bladdr")
```
