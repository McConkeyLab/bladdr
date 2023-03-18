## code to prepare `mtt` dataset goes here
library(bladdr)
library(mop)

mtt <- get_gbci_file("Raw Data/SPECTRAmax/aragaki-kai/mtt/2022-09-20_upfl1-erda-pemi-gefi.txt") |>
  read_spectramax()

usethis::use_data(mtt, overwrite = TRUE)
