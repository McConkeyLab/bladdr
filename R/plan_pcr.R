plan_pcr <- function(data, n_primers, format = 384, exclude_border = TRUE,
                     primer_names = NULL, headless = TRUE, has_names = TRUE,
                     make_report = FALSE, file_path = NULL) {

  data <- tibble::as_tibble(data)

  format <- match.arg(as.character(format), c("96", "384"))

  # TODO
  # Should probably include ... in middle of long string in both here and plan-pcr for display purposes.
  # Add documentation (previously comments) to their respective functions

  if (!headless) {
    browseURL("https://kai-a.shinyapps.io/plan-pcr/")
    return(invisible())
  }

  # Experimental constants -----------------------------------------------------
  ntc <- 1 # Non-targeting control - add one to sample number
  final_rna_conc <- 5#ng/uL Final [RNA], determined by protocol
  reps <- 3 # Perform in triplicate
  safety_reps <- 6 # Extra, since nothing is ever perfect
  rna_per_well <- 2#uL Vol RNA/well

  # Runtime constants ----------------------------------------------------------
  sample_names <- get_sample_names(data, has_names)

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

  ## Make sections -------------------------------------------------------------
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
    dplyr::mutate(dilution_factor = get_best_factor(vol_to_add)) |>
    dplyr::ungroup() |>
    dplyr::mutate(diluted_concentration = data[[2]] / dilution_factor,
           final_vol = final_vol,
           diluted_rna_to_add = final_rna_conc * final_vol / diluted_concentration,
           water_to_add = final_vol - diluted_rna_to_add) |>
    dplyr::select(-vol_to_add) |>
    dplyr::relocate(final_vol, .after = dplyr::last_col())

  # Mastermix Preparation ------------------------------------------------------
  mm <- tibble::tibble(
    reagent = c("2X RT-PCR Buffer", "Primer", "25X Enzyme Mix", "Nuclease Free H2O"),
    vol = c(6.25, .625, .5, 3.125) * (n_samples + ntc + 2) * reps
  )

  if (make_report) {
    if (missing(file_path)) {
      file_path <- tempfile(pattern = paste0(Sys.Date(), "_pcr-report_"), fileext = ".html")
    }

    params <- list(sample_prep = sample_prep,
                   mm_prep = mm,
                   plate = plate,
                   n_primers = n_primers,
                   primer_names = primer_names,
                   format = format,
                   exclude_border = exclude_border)

    rmarkdown::render("./R/pcr_report-template.Rmd", output_file = file_path,
                      params = params, envir = new.env(parent = globalenv()))

    return(list(master_mix_prep = mm, sample_prep = sample_prep, plate = plate, report_path = file_path))
  }
  list(master_mix_prep = mm, sample_prep = sample_prep, plate = plate)
}

useable_plate_dims <- function(plate) {
  list(rows = length(unique(plate$row)),
       cols = length(unique(plate$col)))
}

plate_template <- function(n_row, n_col) {
  tidyr::expand_grid(col = 1:n_col,
                     row = letters[1:n_row]) |>
    dplyr::mutate(row = as.factor(row) |> forcats::fct_rev(),
                  col = as.factor(col))
}

get_plate <- function(format, no_border) {
  plate <- if (format == 96) plate_template(8, 12) else plate_template(16, 24)
  if (no_border & format == 96) {
    plate <- plate |>
      dplyr::filter(!(col %in% c(1, 12) | row %in% c("a", "h")))
  } else if (no_border & format == 384) {
    plate <- plate |>
      dplyr::filter(!(col %in% c(1, 24) | row %in% c("a", "p")))
  } else {
    plate
  }
  plate
}

get_best_factor <- function(vol_to_add) {
  if (vol_to_add < 1) {
    exact_factor <- 1 / vol_to_add
    best_factor <- ceiling(exact_factor/5) * 5 # Give something divisible by 5
  } else {
    best_factor <- 1
  }
  as.integer(best_factor)
}

get_sample_names <- function(data, has_names) {
  if (ncol(data) == 1 & has_names) {
    stop("Data only has one column - did you mean `has_names = FALSE`?")
  }
  if (has_names) data[[1]] else paste("Sample", 1:nrow(data))
}

denote_vlane <- function(plate, plate_dims, reps) {
  n_lanes <- plate_dims$cols %/% reps
  wells_per_lane <- plate_dims$rows * reps
  lanes <- rep(1:n_lanes, each = wells_per_lane)
  length(lanes) <- nrow(plate) # Includes 'unusable' dims (ie borders of borderless)
  dplyr::arrange(plate, col, dplyr::desc(row)) |>
    dplyr::mutate(lane_v = lanes)
}

denote_hlane <- function(plate, plate_dims, n_samples, ntc) {
  total_samples <- n_samples + ntc
  num_lanes <-  plate_dims$rows %/% total_samples
  wells_per_lane <- plate_dims$cols * total_samples
  lanes <- rep(1:num_lanes, each = wells_per_lane)
  length(lanes) <- nrow(plate)
  dplyr::arrange(plate, dplyr::desc(row), col) |>
    dplyr::mutate(lane_h = lanes)
}

flow_lanes <- function(plate, n_primers, n_samples, ntc, reps) {
  total_samples <- n_samples + ntc
  primers <- rep(1:n_primers, each = total_samples * reps)
  length(primers) <- nrow(plate)
  dplyr::arrange(plate, lane_v, dplyr::desc(row)) |>
    dplyr::mutate(primer = primers,
           primer = dplyr::if_else(primer <= n_primers,
                            primer,
                            NA_integer_),
           available_well = TRUE)
}

section_lanes <- function(laned_plate, n_primers) {
  laned_plate |>
    dplyr::mutate(primer = lane_h * 100 + lane_v) |> # Numerical hierarchy for ordering
    dplyr::arrange(primer) |>
    dplyr::group_by(primer) |>
    dplyr::mutate(primer = dplyr::if_else(dplyr::cur_group_id() <= n_primers,
                                          dplyr::cur_group_id(),
                            NA_integer_),
                  available_well = TRUE) |>
    dplyr::ungroup()
}

add_sample_names <- function(plate, reps, sample_names) {
  plate <- plate |>
    dplyr::group_by(primer) |>
    dplyr::mutate(sample = (dplyr::row_number() + 2) %/% reps,
           sample = dplyr::if_else(!is.na(primer), sample, NA_real_)) |>
    dplyr::group_by(sample) |>
    tidyr::nest() |>
    dplyr::ungroup()

  sample_names <- c(sample_names, "NTC")
  length(sample_names) <- nrow(plate)

  plate |>
    dplyr::mutate(sample_name = sample_names) |>
    tidyr::unnest(cols = data) |>
    dplyr::group_by(sample, primer) |>
    dplyr::arrange(dplyr::desc(row), col) |>
    dplyr::mutate(sample_name = dplyr::if_else(!is.na(primer) & dplyr::row_number() == 2, sample_name, NA_character_),
           sample = as.factor(sample))
}

add_primer_names <- function(plate, primer_names) {
  plate <- plate |>
    dplyr::group_by(primer) |>
    tidyr::nest() |>
    dplyr::arrange(primer) |>
    dplyr::ungroup()

  length(primer_names) <- nrow(plate) # Account for any NAs (that may or may not be present)

  plate |>
    dplyr::mutate(primer_name = primer_names) |>
    tidyr::unnest(cols = data) |>
    dplyr::group_by(primer) |>
    dplyr::arrange(lane_v, dplyr::desc(row), col) |>
    dplyr::mutate(primer_name = dplyr::if_else(dplyr::row_number() == 2, primer_name, NA_character_),
           primer = as.factor(primer))
}
