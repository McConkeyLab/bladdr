#' Plan PCR experiment
#'
#' @param data a data.frame, with samples as the first column (if  `has_names = TRUE`) and RNA concentrations as the second (or first, if `has_names = FALSE`)
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
plan_pcr <- function(data, n_primers, format = 384, exclude_border = TRUE,
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
    reagent = c("2X RT-PCR Buffer", "Primer", "25X Enzyme Mix", "Nuclease Free H2O"),
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

    rmarkdown::render("./inst/pcr_report-template.Rmd", output_file = file_path,
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
    dplyr::mutate(primer = .data$primers,
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
