# bladdr: A Helping Hand in McConkey Land <img src='man/figures/logo.png' align="right" height="139" />

bladdr's vision is to be a package that can automate everyday days in small functions that don't try to do too much. 

# Guiding Principles

## Reproducibility

Any output of these functions should produce a report with all relevant information.

## Clarity

Functions should strive to be simple over elegant when one must be decided over the other.
Functions should be allowed to be opinionated about the input of data, as the input format is often not optimal for working with.
* However, input data should not be removed without precedent, only reshaped or refactored. The goal of these functions is to have the be able to be peppered in to every day workflows. Therefore they must be flexible and humble in their function, lest someone want to use them in a different way.
* For now, until it causes problems, we will assume that the tidyverse exists and use of it in our package will be encouraged for the sake of piping and clarity.

## Modularization

Functions should only do one thing: creating two small functions, one which feeds into the other, is preferred over one larger function
