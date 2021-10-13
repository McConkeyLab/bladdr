#' (Re)calculate Delta Ct Mean
#'
#' @param data A dataset output from pcr_tidy
#' @param control_probe A probe to be used as an endogenous control (eg GAPDH)
#'
#' @return A tibble
#' @export
#'
#' @importFrom rlang .data
#'
#' @examples
#' dat_path <- system.file("extdata", "untidy-pcr-example.xls", package = "bladdr")
#'
#' pcr_tidy(dat_path) |>
#' pcr_control("GAPDH")
pcr_control <- function(data, control_probe) {
  data |>
    dplyr::group_by(.data$target_name, .data$sample_name) |>
    dplyr::mutate(ct_mean = mean(.data$ct),
                  ct_sd   = stats::sd(.data$ct),
                  rep     = dplyr::n()) |>
    dplyr::ungroup() |>
    tidyr::nest(sample_nest = c(.data$ct, .data$well, .data$well_row,
                                .data$well_col, .data$well_position,
                                .data$baseline_start, .data$baseline_end)) |>
    dplyr::group_by(.data$sample_name) |>
    dplyr::mutate(delta_ct     = .data$ct_mean - .data$ct_mean[.data$target_name == control_probe],
                  delta_ct_sd  = sqrt(.data$ct_sd^2 + .data$ct_sd[.data$target_name == control_probe]^2),
                  delta_ct_se  = .data$delta_ct_sd/sqrt(.data$rep),  # Rep might not be the correct metric here
                  df           = max(1, .data$rep + .data$rep[.data$target_name == control_probe] - 2),
                  t            = stats::qt(.05/2, .data$df, lower.tail = F)) |>
    tidyr::unnest(.data$sample_nest)
}


#' Calculate library PCR concentrations
#'
#' @param tidy_pcr a dataset run that has been run through `pcr_tidy()`
#' @param dil_factor the factor to which the libraries were diluted for pcr
#'
#' @return a list, containing the original dataframe, data on the standards, and data on the samples.
#' @export
#'
#' @importFrom rlang .data
#'
#' @examples
#'
#' dat_path <- system.file("extdata", "untidy-standard-curve.xlsx", package = "bladdr")
#'
#' pcr_tidy(dat_path) |>
#' pcr_lib_calc()

pcr_lib_calc <- function(tidy_pcr, dil_factor = 1000) {

  tidy_pcr |>
    tidyr::nest(replicates = c(.data$well, .data$well_position, .data$ct,
                               .data$quantity, .data$well_row, .data$well_col,
                               dplyr::starts_with("prfdrop"),
                               dplyr::starts_with("baxrox"))) |>
    dplyr::group_by(.data$task) |>
    dplyr::arrange(.data$ct_mean) |>
    dplyr::mutate(standard_diff = .data$ct_mean - dplyr::lag(.data$ct_mean, default = .data$ct_mean[1]),
                  dil = 2^.data$standard_diff,
                  quant_actual = 6.8/cumprod(.data$dil),
                  dil = dplyr::if_else(.data$dil == 1, 0, .data$dil)) |>
    tidyr::unnest(cols = .data$replicates) |>
    dplyr::mutate(dil = dplyr::if_else(.data$task == "STANDARD", .data$dil, NA_real_),
                  standard_diff = dplyr::if_else(.data$task == "STANDARD", .data$standard_diff, NA_real_),
                  quant_actual = dplyr::if_else(.data$task == "STANDARD", .data$quant_actual, .data$quantity),
                  concentration = .data$quantity_mean * dil_factor)
}

#' Generate visual library prep pcr quality control report
#'
#' @param lib_calc_pcr an output from `pcr_lib_calc`
#' @param report_name filename of report to be generated (do not include file extension). Defaults to date-time of report generation.
#' @return a ggplot
#' @export
#'
#' @importFrom ggplot2 aes element_blank
#'
#' @examples
#'
#' dat_path <- system.file("extdata", "untidy-standard-curve.xlsx", package = "bladdr")
#'
#' pcr_tidy(dat_path, pad_zero = TRUE) |>
#'         pcr_lib_calc() |>
#'         pcr_lib_qc()
#'

pcr_lib_qc <- function(lib_calc_pcr, report_name = NULL) {
  dat <- lib_calc_pcr |>
    dplyr::select(c("task", "sample_name", "quantity_mean", "concentration", "quantity",
                    "quant_actual", "dil", "slope", "efficiency", "r2", "ct"))

  standards <- dat |>
    dplyr::filter(.data$task == "STANDARD")

  samples <- dat |>
    dplyr::filter(.data$task == "UNKNOWN")

  sample_summary <- samples |>
    dplyr::group_by(.data$sample_name) |>
    dplyr::summarize(quantity_mean = mean(.data$quantity_mean),
                     concentration_mean = mean(.data$concentration))

  standard_summary <- standards |>
    dplyr::group_by(.data$quantity) |>
    dplyr::summarize(quantity_mean = mean(.data$quantity),
                     quant_actual = mean(.data$quant_actual),
                     dil = mean(.data$dil)) |>
    tidyr::pivot_longer(cols = c(.data$quantity_mean, .data$quant_actual))

  dilution_lines <- standard_summary |>
    dplyr::filter(.data$name == "quant_actual") |>
    dplyr::mutate(line_start = 1/.data$value,
                  line_end = dplyr::lag(.data$line_start),
                  dil = dplyr::lag(.data$dil),
                  y = rep_len(c(1.1, 0.9), 5),
                  y_text = rep_len(c(1.15, 0.85), 5)) |>
    dplyr::filter(!is.na(.data$line_end)) |>
    dplyr::rowwise() |>
    dplyr::mutate(mid = sqrt(.data$line_start * .data$line_end))

  vert_lines <-
    tibble::tibble(x = c(dilution_lines$line_start, dilution_lines$line_end)) |>
    dplyr::arrange(.data$x) |>
    dplyr::mutate(y = rep(c(1.1, 0.9, 1.1, 0.9), each = 2),
                  yend = rep(c(1.05, 0.95, 1.05, 0.95), each = 2))

  standards_plot <-
    ggplot2::ggplot(standard_summary, aes(x = 1/.data$value, y = 1, color = .data$name)) +
    ggplot2::geom_point(size = 10, alpha = 0.7) +
    ggplot2::scale_color_manual(values = c("#00AAAA", "#222222")) +
    ggplot2::geom_segment(data = dilution_lines, aes(x = .data$line_start, xend = .data$line_end, y = .data$y, yend = .data$y), size = 1) +
    ggplot2::geom_segment(data = vert_lines, aes(x = .data$x, xend = .data$x, y = .data$y, yend = .data$yend), color = "#00AAAA", size = 1) +
    ggplot2::geom_text(data = dilution_lines, aes(label = round(.data$dil, 1), y = .data$y_text, x = .data$mid), size = 8) +
    ggplot2::scale_x_log10() +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid = element_blank(), axis.title = element_blank(), axis.text = element_blank(), legend.position = "none") +
    ggplot2::coord_cartesian(ylim = c(0.8, 1.2)) +
    ggplot2::geom_point(data = sample_summary, aes(x = 1/.data$quantity_mean, y = 1), color = "red")

  slope_text <-
    dplyr::tibble(x = -2.5,
                  y = 20,
                  label = paste(paste("Slope:", round(standards$slope[1], 2)),
                                paste("R\u00B2:", round(standards$r2[1], 2)),
                                paste0("Efficiency: ", round(standards$efficiency[1], 1), "%"),
                                sep = "\n"))

  slope_plot <-
    ggplot2::ggplot(standards, aes(x = log10(.data$quantity), y = .data$ct)) +
    ggplot2::geom_point(color = "#00AAAA", size = 5) +
    ggplot2::geom_smooth(method = "lm", se = F, color = "#666666") +
    ggplot2::geom_text(data = slope_text, aes(x = .data$x,  y = .data$y, label = .data$label), size = 8) +
    ggplot2::xlab("Log10(Quantity)") +
    ggplot2::ylab("CT")

  outliers_plot_dat <- dat |>
    tidyr::nest(reps = c(.data$ct, .data$quantity, .data$quant_actual, .data$concentration)) |>
    dplyr::mutate(mean_wo_outlier = purrr::map(.data$reps, find_mean)) |>
    tidyr::unnest_wider(.data$mean_wo_outlier) |>
    tidyr::unnest(cols = c(.data$no_po_mean, .data$no_po_sd, .data$reps)) |>
    dplyr::mutate(keep = dplyr::case_when(.data$no_po_mean - (3*.data$no_po_sd) < .data$ct & .data$no_po_mean + (3 * .data$no_po_sd) > .data$ct ~ TRUE,
                                          is.na(.data$no_po_sd) ~ TRUE,
                                          TRUE ~ NA),
                  sample_name = dplyr::case_when(is.na(.data$sample_name) & .data$quantity > 6.0000 ~ "1",
                                                 is.na(.data$sample_name) & .data$quantity > 0.6000 ~ "1:10",
                                                 is.na(.data$sample_name) & .data$quantity > 0.0600 ~ "1:100",
                                                 is.na(.data$sample_name) & .data$quantity > 0.0060 ~ "1:1000",
                                                 is.na(.data$sample_name) & .data$quantity > 0.0006 ~ "1:10000",
                                                 TRUE ~   .data$sample_name)) |>
    dplyr::group_by(.data$sample_name) |>
    dplyr::mutate(adj_mean = mean(.data$keep * .data$ct, na.rm = TRUE),
                  adj_sd   = stats::sd(.data$keep*.data$ct, na.rm = TRUE),
                  keep_logi = !is.na(.data$keep),
                  z = (.data$ct-.data$adj_mean)/.data$adj_sd,
                  overflow = abs(.data$z) > 10,
                  z_plot = dplyr::case_when(.data$z > 10 ~ 10,
                                            .data$z < -10 ~ -10,
                                            TRUE ~ .data$z),
                  label = dplyr::case_when(.data$z > 10 ~ paste(as.character(round(.data$z, 0)),">>>"),
                                           .data$z < -10 ~ paste("<<<", as.character(round(.data$z, 0))),
                                           TRUE ~ NA_character_)) |>
    dplyr::filter(!is.na(.data$sample_name))

  # Plot z-score (y) of samples (x) with a shading window of -3 < z < 3.
  # Z-scores are centered using the outlier-removed mean
  outliers_plot <-
    ggplot2::ggplot(outliers_plot_dat, aes(x = .data$z_plot, y = .data$sample_name, color = .data$keep_logi, shape = .data$overflow)) +
    ggplot2::scale_color_manual(values = c("#666666", "#00AAAA")) +
    ggplot2::scale_shape_manual(values = c(16, NA)) +
    ggplot2::geom_text(aes(label = .data$label), size = 5) +
    ggplot2::geom_rect(xmin = -3, xmax = 3, ymin = 0, ymax = nrow(outliers_plot_dat), fill = "#AACCCC",  alpha = 0.1,  color = NA) +
    ggplot2::geom_vline(xintercept = 0, color = "#555555") +
    ggplot2::geom_point(size = 5) +
    ggplot2::xlab("Z-score") +
    ggplot2::coord_cartesian(xlim = c(-11.5, 11.5), ylim = c(0.5, length(unique(outliers_plot_dat$sample_name)) + 0.5), expand = F) +
    ggplot2::theme(panel.grid.major.x = element_blank(), axis.title.y = element_blank(), panel.grid.minor.x = element_blank(), legend.position = "none")


  table_samples <- sample_summary |>
    dplyr::select(-"quantity_mean") |>
    dplyr::rename("Sample Name" = .data$sample_name,
                  "Concentration" = .data$concentration_mean)

  # Plot library concentration (y) of samples (x)
  conc_plot <-
    ggplot2::ggplot(table_samples, aes(x = .data$`Sample Name`, y = .data$`Concentration`)) +
    ggplot2::geom_point(color = "#00AAAA", size = 5) +
    ggplot2::theme(axis.title = element_blank(), legend.position = "none", axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5))

  report <- system.file("lib-qc.Rmd", package = "bladdr")

  # If there's a user supplied filename, use that. If not, use current date and time.
  if (is.null(report_name)) {
    report_name <- paste0(gsub(":", "-", gsub("\\s", "_", Sys.time())),"_report.html")
  } else {
    report_name <- paste0(report_name, ".html")
  }

  # Generate report
  rmarkdown::render(report,
                    output_dir = "./pcr_qc_reports",
                    output_file = report_name,
                    params = list(stan = standards_plot, out = outliers_plot, samp = table_samples, conc = conc_plot, slope = slope_plot))
}

find_mean <- function(df){
  if(nrow(df) >= 3 & !all(is.na(df$ct))) {
    hc <- df$ct |>
      stats::dist() |>
      stats::hclust()
    possible_outlier <- hc$merge[nrow(df) - 1, 1] |>
      abs()
    no_po <- df$ct[-possible_outlier]
    no_po_mean <- mean(no_po)
    no_po_sd <- stats::sd(no_po)
    list(no_po_mean = no_po_mean, no_po_sd = no_po_sd)
  } else {
    no_po_mean <- mean(df$ct)
    no_po_sd <- stats::sd(df$ct, na.rm = T)
    list(no_po_mean = no_po_mean, no_po_sd = no_po_sd)
  }
}

#' View sample plating layout
#'
#' @param tidy_pcr an output from the `pcr_tidy` function
#' @param fill the variable to use to fill the geom_tiles
#'
#' @return a ggplot
#' @export
#'
#' @importFrom ggplot2 aes
#'
#' @examples
#'
#' dat_path <- system.file("extdata", "untidy-pcr-example.xls", package = "bladdr") |>
#' pcr_tidy() |>
#' pcr_plate_view()

pcr_plate_view <- function(tidy_pcr, fill = .data$target_name) {
  usr_fill <- substitute(fill)
  if (sum(grepl("384", tidy_pcr$plate_type), na.rm = TRUE) > 0) {
    x <- 2
  } else if (sum(grepl("96", tidy_pcr$plate_type), na.rm = TRUE) > 0) {
    x <- 1
  }
  ggplot2::ggplot(tidy_pcr, aes(x = .data$well_col, y = .data$well_row, fill = eval(usr_fill))) +
    ggplot2::geom_tile(aes(size = 2)) +
    ggplot2::coord_cartesian(xlim = c(1,(x*12)), ylim = c((x*8), 1)) +
    ggplot2::scale_y_continuous(breaks = 1:(x*8), labels = LETTERS[1:(x*8)]) +
    ggplot2::scale_x_continuous(breaks = 1:(x*12), minor_breaks = NULL) +
    ggplot2::labs(fill = deparse(usr_fill)) +
    ggplot2::guides(size = FALSE)
}


#' Plot qPCR results
#'
#' @param tidy_pcr an output from the `pcr_tidy` function, or some derivative thereof
#'
#' @return a ggplot
#' @export
#'
#' @importFrom rlang .data
#'
#' @examples
#'
#' dat_path <- system.file("extdata", "untidy-pcr-example.xls", package = "bladdr")
#' pcr_tidy(dat_path) |>
#' pcr_rq("RD1") |>
#' pcr_plot()
#'

pcr_plot <- function(tidy_pcr) {
  tidy_pcr |>
    dplyr::filter(!is.na(.data$sample_name)) |>
    dplyr::distinct(.data$target_name, .data$sample_name, .keep_all = T) |>
    ggplot2::ggplot(ggplot2::aes(x = .data$sample_name, y = .data$rq, fill = .data$target_name)) +
    ggplot2::geom_col() +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$rq_min, ymax = .data$rq_max)) +
    ggplot2::facet_wrap(~.data$target_name, scales = "free") +
    ggplot2::scale_fill_viridis_d(begin = 0.1, end = 0.9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
                   axis.title.x = ggplot2::element_blank())
}

#' (Re)calculate rq for a given sample
#'
#' @param data A dataset output from pcr_tidy/pcr_control
#' @param relative_sample A sample to set others relative to (eg my_dmso_sample)
#'
#' @return A tibble
#' @export
#'
#' @importFrom rlang .data
#'
#' @examples
#' dat_path <- system.file("extdata", "untidy-pcr-example.xls", package = "bladdr")
#'
#' pcr_tidy(dat_path) |>
#' pcr_rq("U6D1")
#'
#' # Can also be run after using pcr_control:
#' pcr_tidy(dat_path) |>
#' pcr_control("GAPDH") |>
#' pcr_rq("U6D1")
pcr_rq <- function(data, relative_sample) {

  control_probe <- unique(data$control)

  data |>
    dplyr::group_by(.data$target_name, .data$sample_name) |>
    dplyr::mutate(ct_mean = mean(.data$ct),
                  ct_sd   = stats::sd(.data$ct),
                  rep     = dplyr::n()) |>
    dplyr::ungroup() |>
    tidyr::nest(sample_nest = c("ct", "well", "well_row", "well_col", "well_position",
                                "baseline_start", "baseline_end")) |>
    dplyr::group_by(.data$sample_name) |>
    dplyr::mutate(delta_ct     = .data$ct_mean - .data$ct_mean[.data$target_name == control_probe],
                  delta_ct_sd  = sqrt(.data$ct_sd^2 + .data$ct_sd[.data$target_name == control_probe]^2),
                  delta_ct_se  = .data$delta_ct_sd/sqrt(.data$rep),  # Rep might not be the correct metric here
                  df           = max(1, .data$rep + .data$rep[.data$target_name == control_probe] - 2),
                  t            = stats::qt(.05/2, .data$df, lower.tail = FALSE)) |>
    dplyr::group_by(.data$target_name) |>
    dplyr::filter(!is.na(.data$sample_name)) |>
    dplyr::mutate(delta_delta_ct = .data$delta_ct - .data$delta_ct[.data$sample_name == relative_sample],
                  rq     = 2^-(.data$delta_delta_ct),
                  rq_min = 2^-(.data$delta_delta_ct + .data$t * .data$delta_ct_se),
                  rq_max = 2^-(.data$delta_delta_ct - .data$t * .data$delta_ct_se)) |>
    tidyr::unnest(.data$sample_nest)
}

#' Tidy a PCR Excel File
#'
#' Takes in a fresh results output from qPCR and converts it into a tidy
#' format. Useful for downstream analyses.
#'
#' @param file_path Path to an excel file containing the results of a qPCR
#'   run. If left blank, will open up interactive file chooser.
#' @param pad_zero Should a leading zero be added to "Sample 1" etc?
#'
#' @return Tidy results dataframe
#' @export
#'
#' @importFrom rlang .data
#'
#' @examples
#' dat_path <- system.file("extdata", "untidy-pcr-example.xls", package = "bladdr")
#'
#' # Before tidying
#' dat_dirty <- readxl::read_excel(dat_path, sheet = "Results")
#' dat_dirty[1:10]
#'
#' # After tidying
#' dat_clean <- pcr_tidy(dat_path)
#' dat_clean[1:10]

pcr_tidy <- function(file_path = NULL, pad_zero = FALSE) {

  if (is.null(file_path)) {
    file_path <- file.choose()
  }

  dat_og <- readxl::read_excel(file_path, sheet = "Results")

  exp_type <- substr(dat_og[which(dat_og$`Block Type`=="Experiment Type"),2], 0, 4)

  ind_start <- which(dat_og[,1] == "Well")

  if (exp_type == "Stan") {
    ind_end <- nrow(dat_og)
    dat <- dat_og[ind_start:ind_end,]
  }
  if (exp_type == "Comp") {
    ind_end <- which(dat_og[,1] == "Analysis Type")
    exp_dat <- dat_og[ind_end:nrow(dat_og), 1:2] |>
      t()
    colnames(exp_dat) <- make.names(exp_dat[1,])
    exp_dat <- exp_dat[-1,]
    exp_dat <- t(exp_dat) |>
      as.data.frame()
    dat <- dat_og[-c(1:(ind_start-1), (ind_end-1):nrow(dat_og)),]
  }

  names <- gsub(" ", "_", dat[1,])
  names <- gsub("-", "_", names)
  names <- gsub("\\(superscript_2\\)", "2", names)
  names <- tolower(names)
  colnames(dat) <- names

  dat <- dat[-1,] |>
    dplyr::mutate(dplyr::across(dplyr::matches("^(delta )*ct.*|^rq|quantity|^baseline|y_intercept|r2|slope|efficiency"), as.numeric),
                  well_row = stringr::str_extract(.data$well_position, "^.{1}"),
                  well_col = as.numeric(stringr::str_extract(.data$well_position, "[:digit:]{1,2}$")),
                  well_row = as.numeric(factor(.data$well_row, levels = LETTERS)))

  if (exp_type == "Comp") {
    dat$analysis_type <- exp_dat$Analysis.Type
    dat$control <- exp_dat$Endogenous.Control
    dat$conf_int <- exp_dat$RQ.Min.Max.Confidence.Level
    dat$ref_samp <- exp_dat$Reference.Sample
  }

  if (pad_zero) {
    dat$sample_name <- pad_zero(dat$sample_name)
  }

  dat$plate_type <- colnames(dat_og)[2]
  dat$exp_type <- tolower(exp_type[[1]])

  dat
}

#' Plan PCR experiment
#'
#' @param data a data.frame, with samples as the first column (if `has_names = TRUE`) and RNA concentrations as the second (or first, if `has_names = FALSE`)
#' @param n_primers integer. Number of primers to be used in the experiment.
#' @param format integer. 96 or 384 - the number of wells of the plate planned to be used
#' @param exclude_border logical. Should the border be excluded to avoid edge effects? Default is TRUE.
#' @param primer_names character vector. Names of primers.
#' @param headless logical. If FALSE, return invisible and redirect to shiny application.
#' @param has_names logical. Is the first column the names of the samples?
#' @param make_report logical. Should an HTML report be written?
#' @param file_path Optional path to where the report should be written, as well as the file name. Defaults to temp file.
#'
#' @return a list
#' @export
pcr_plan <- function(data, n_primers, format = 384, exclude_border = TRUE,
                     primer_names = NULL, headless = TRUE, has_names = TRUE,
                     make_report = FALSE, file_path = NULL) {

  # TODO
  # Add documentation (previously comments) to their respective functions

  # Checks ---------------------------------------------------------------------
  if (!headless) {
    utils::browseURL("https://kai-a.shinyapps.io/plan-pcr/")
    return(invisible())
  }

  format <- match.arg(as.character(format), c("96", "384"))

  data <- tibble::as_tibble(data) # Allows for vector input

  if (ncol(data) == 1 & has_names) {
    stop("Data only has one column - did you mean `has_names = FALSE`?")
  }

  # Experimental constants -----------------------------------------------------
  ntc <- 1 # Non-targeting control - add one to sample number
  final_rna_conc <- 5#ng/uL Final [RNA], determined by protocol
  reps <- 3 # Perform in triplicate
  safety_reps <- 6 # Extra, since nothing is ever perfect
  rna_per_well <- 2#uL Vol RNA/well

  # Runtime constants ----------------------------------------------------------
  sample_names <- data |>
    get_sample_names(has_names)

  if (!has_names) {
    data <- cbind(sample_names, data)
  }

  final_vol <- ((n_primers * reps) + safety_reps) * rna_per_well |> as.integer()
  n_samples <- nrow(data)
  section_area <- reps * (n_samples + ntc)
  full_plate <- get_plate(format, no_border = FALSE)
  this_plate <- get_plate(format, no_border = exclude_border)
  plate_dims <- useable_plate_dims(this_plate)
  max_sections_before_flow <- (plate_dims$cols %/% reps) * (plate_dims$rows %/% (n_samples + ntc)) # max_wide * max_tall
  max_sections_after_flow <-((plate_dims$cols %/% reps) * plate_dims$rows) %/% (n_samples + ntc) # (max_wide * n_row) %/% tot_samples
  max_sections_theoretical <- (plate_dims$rows * plate_dims$cols) %/% section_area

  # Checks ---------------------------------------------------------------------
  if(max_sections_theoretical < n_primers) {
    stop("This experiment requires too many wells.")
  }

  # Primer names ---------------------------------------------------------------
  pn <- paste("Primer", 1:n_primers)

  if (!missing(primer_names)) {
    pn[1:length(primer_names)] <- primer_names
  }

  # Make sections --------------------------------------------------------------
  if (n_primers > max_sections_before_flow) {
    plate <- this_plate |>
      denote_vlane(plate_dims, reps) |>
      flow_lanes(n_primers, n_samples, ntc, reps)
  } else {
    plate <- this_plate |>
      denote_vlane(plate_dims, reps) |>
      denote_hlane(plate_dims, n_samples, ntc) |>
      section_lanes(n_primers)
  }

  plate <- plate |>
    add_sample_names(reps, sample_names) |>
    add_primer_names(pn) |>
    dplyr::right_join(full_plate, by = c("col", "row"))

  # Sample preparation  --------------------------------------------------------
  sample_prep <- data |>
    dplyr::mutate(vol_to_add = final_rna_conc * final_vol / data[[2]]) |>
    dplyr::rowwise() |>
    dplyr::mutate(dilution_factor = get_best_factor(.data$vol_to_add)) |>
    dplyr::ungroup() |>
    dplyr::mutate(diluted_concentration = data[[2]] / .data$dilution_factor,
                  final_vol = final_vol,
                  diluted_rna_to_add = final_rna_conc * .data$final_vol / .data$diluted_concentration,
                  water_to_add = .data$final_vol - .data$diluted_rna_to_add) |>
    dplyr::select(-.data$vol_to_add) |>
    dplyr::relocate(.data$final_vol, .after = dplyr::last_col())

  # Mastermix Preparation ------------------------------------------------------
  mm <- tibble::tibble(
    reagent = c("2X RT-PCR Buffer", "Primer", "25X RT-PCR Enzyme", "Nuclease Free H2O"),
    vol = c(6.25, .625, .5, 3.125) * (n_samples + ntc + 2) * reps
  )

  if (make_report) {
    if (missing(file_path)) {
      file_path <- tempfile(pattern = paste0(Sys.Date(), "_pcr-report_"),
                            fileext = ".html")
    }

    params <- list(sample_prep = sample_prep,
                   mm_prep = mm,
                   plate = plate,
                   n_primers = n_primers,
                   primer_names = primer_names,
                   format = format,
                   exclude_border = exclude_border)

    rmarkdown::render(system.file("rmd", "pcr_report-template.Rmd", package = "bladdr"), output_file = file_path,
                      params = params, envir = new.env(parent = globalenv()))

    return(list(master_mix_prep = mm, sample_prep = sample_prep,
                plate = plate, report_path = file_path))
  }
  list(master_mix_prep = mm, sample_prep = sample_prep, plate = plate)
}

#' Get or make sample names
#'
#' If the user did not supply sample names (`has_names = FALSE`), sample names
#' will be generated in the form of "Sample_n".
#'
#' @param data A data.frame
#' @param has_names logical. Should we expect sample names to be in the first
#'   column?
#'
#' @return A vector of sample names
#'
#' @keywords internal
get_sample_names <- function(data, has_names) {
  if (has_names) data[[1]] else paste("Sample", 1:nrow(data))
}

#' Create a plate with given well format
#'
#' @param format integer. How many wells should this plate have?
#' @param no_border logical. Will the user omit filling the edge wells of the
#'   plate (to avoid edge effects)?
#'
#' @return A tibble with `col` (number, factor, denotes plate column from
#'   left to right) and `row` (character, factor, denotes plate row from top to
#'   bottom). If `no_border = TRUE`, the edges of the plate will not be included.
#'
#' @keywords internal
get_plate <- function(format, no_border) {
  plate <- if (format == 96) plate_template(8, 12) else plate_template(16, 24)
  if (no_border & format == 96) {
    plate <- plate |>
      dplyr::filter(!(.data$col %in% c(1, 12) | .data$row %in% c("a", "h")))
  } else if (no_border & format == 384) {
    plate <- plate |>
      dplyr::filter(!(.data$col %in% c(1, 24) | .data$row %in% c("a", "p")))
  } else {
    plate
  }
  plate
}

#' Create a plate of given dimensions
#'
#' @param n_row integer. Number of rows plate should have.
#' @param n_col integer. Number of columns plate should have.
#'
#' @return A tibble with `col` (a factor of numbers with the lowest being the
#'   first and the highest being the last) and `row` (a factor of characters in
#'   reverse alphabetical order, to mimic ordering in real world plates when
#'   plotting)
#'
#' @keywords internal
plate_template <- function(n_row, n_col) {
  tidyr::expand_grid(col = 1:n_col,
                     row = letters[1:n_row]) |>
    dplyr::mutate(row = as.factor(.data$row) |> forcats::fct_rev(),
                  col = as.factor(.data$col))
}


#' Get dimensions of plate that are allowed to be filled
#'
#' The 'useable' dimensions of the plate depends if the user has set
#' `exclude_border` to be `TRUE` or `FALSE`. This calculates plate dimensions
#' dependent on that.
#'
#' @param plate
#'
#' @return a named list, with integer values referring to the number of rows and
#'   columns, respectively.
#'
#' @keywords internal
useable_plate_dims <- function(plate) {
  list(rows = length(unique(plate$row)),
       cols = length(unique(plate$col)))
}


#' Calculate a sensible dilution factor
#'
#' If the volume of RNA to add is < 1uL, it must be diluted. Dilutions are easy
#' to calculate in one's head only if they are integers (divisible by 5
#' preferred). Further, this dilution should be as small as reasonably possible,
#' otherwise it will become too dilute.
#'
#' @param vol_to_add numeric. 'Naive' volume to add, before dilution.
#'
#' @return an integer, either 1 (no dilution) or something divisible by 5.
#'
#' @keywords internal
get_best_factor <- function(vol_to_add) {
  if (vol_to_add < 1) {
    exact_factor <- 1 / vol_to_add
    best_factor <- ceiling(exact_factor/5) * 5 # Give something divisible by 5
  } else {
    best_factor <- 1
  }
  as.integer(best_factor)
}

denote_vlane <- function(plate, plate_dims, reps) {
  n_lanes <- plate_dims$cols %/% reps
  wells_per_lane <- plate_dims$rows * reps
  lanes <- rep(1:n_lanes, each = wells_per_lane)
  length(lanes) <- nrow(plate) # Includes 'unusable' dims (ie borders of borderless)
  dplyr::arrange(plate, .data$col, dplyr::desc(.data$row)) |>
    dplyr::mutate(lane_v = lanes)
}

denote_hlane <- function(plate, plate_dims, n_samples, ntc) {
  total_samples <- n_samples + ntc
  num_lanes <-  plate_dims$rows %/% total_samples
  wells_per_lane <- plate_dims$cols * total_samples
  lanes <- rep(1:num_lanes, each = wells_per_lane)
  length(lanes) <- nrow(plate)
  dplyr::arrange(plate, dplyr::desc(.data$row), .data$col) |>
    dplyr::mutate(lane_h = lanes)
}

flow_lanes <- function(plate, n_primers, n_samples, ntc, reps) {
  total_samples <- n_samples + ntc
  primers <- rep(1:n_primers, each = total_samples * reps)
  length(primers) <- nrow(plate)
  dplyr::arrange(plate, .data$lane_v, dplyr::desc(.data$row)) |>
    dplyr::mutate(primer = primers,
                  primer = dplyr::if_else(.data$primer <= n_primers,
                                          .data$primer,
                                          NA_integer_),
                  available_well = TRUE)
}

section_lanes <- function(laned_plate, n_primers) {
  laned_plate |>
    dplyr::mutate(primer = .data$lane_h * 100 + .data$lane_v) |> # Numerical hierarchy for ordering
    dplyr::arrange(.data$primer) |>
    dplyr::group_by(.data$primer) |>
    dplyr::mutate(primer = dplyr::if_else(dplyr::cur_group_id() <= n_primers,
                                          dplyr::cur_group_id(),
                                          NA_integer_),
                  available_well = TRUE) |>
    dplyr::ungroup()
}

add_sample_names <- function(plate, reps, sample_names) {
  plate <- plate |>
    dplyr::group_by(.data$primer) |>
    dplyr::mutate(sample = (dplyr::row_number() + 2) %/% reps,
                  sample = dplyr::if_else(!is.na(.data$primer), .data$sample, NA_real_)) |>
    dplyr::group_by(.data$sample) |>
    tidyr::nest() |>
    dplyr::ungroup()

  sample_names <- c(sample_names, "NTC")
  plotting_names <- stringr::str_trunc(sample_names, width = 15, side = "center", ellipsis = "\U2026")
  length(sample_names) <- nrow(plate)
  length(plotting_names) <- nrow(plate)

  plate |>
    dplyr::mutate(sample_name = sample_names,
                  plotting_name = plotting_names) |>
    tidyr::unnest(cols = .data$data) |>
    dplyr::group_by(sample, .data$primer) |>
    dplyr::arrange(dplyr::desc(.data$row), .data$col) |>
    dplyr::mutate(sample_name = dplyr::if_else(!is.na(.data$primer) & dplyr::row_number() == 2, .data$sample_name, NA_character_),
                  plotting_name = dplyr::if_else(!is.na(.data$primer) & dplyr::row_number() == 2, .data$plotting_name, NA_character_),
                  sample = as.factor(.data$sample))
}

add_primer_names <- function(plate, primer_names) {
  plate <- plate |>
    dplyr::group_by(.data$primer) |>
    tidyr::nest() |>
    dplyr::arrange(.data$primer) |>
    dplyr::ungroup()

  length(primer_names) <- nrow(plate) # Account for any NAs (that may or may not be present)

  plate |>
    dplyr::mutate(primer_name = primer_names) |>
    tidyr::unnest(cols = .data$data) |>
    dplyr::group_by(.data$primer) |>
    dplyr::arrange(.data$lane_v, dplyr::desc(.data$row), .data$col) |>
    dplyr::mutate(primer_name = dplyr::if_else(dplyr::row_number() == 2, .data$primer_name, NA_character_),
                  primer = as.factor(.data$primer))
}
