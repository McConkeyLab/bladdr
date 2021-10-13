df <- expand.grid(time = c("24hr", "48hr"), drug_conc = c("DMSO", "1uM"), rep = c("1", "2"))
df$sample <- paste(df$time, df$drug_conc, df$rep, sep = "_")

set.seed(808)
conc <-
  rnorm(nrow(df),
        mean = 100,
        sd = 20) |>
  abs() |>
  round(2)

dummy_rna_conc <- data.frame(sample = df$sample,
                             conc = conc)

usethis::use_data(dummy_rna_conc, overwrite = TRUE)
