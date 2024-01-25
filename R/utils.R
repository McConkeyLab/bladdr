get_username_from_env <- function() {
  usr <- Sys.getenv("GBCI_USERNAME")
  ifelse(usr == "", username_env_helper(), usr)
}

username_env_helper <- function() {
  rlang::abort(c(
    "The environmental username GBCI_USERNAME is unset.",
    "To set for future R sessions, add the following line to your .Rprofile:",
    "Sys.setenv(GBCI_USERNAME = \"YOUR-GBCI-USERNAME\")",
    "Open your .Rprofile with usethis::edit_r_profile"
  ))
}

#' Get indices of genes by name in a dds
#' @noRd
get_gene_index <- function(dds, genes) {
  inds <- which(rownames(dds) %in% genes)
  if (identical(inds, integer(0))) {
    stop("None of the supplied genes exist in this dataset.")
  } else if (length(genes) != length(inds)) {
    found <- rownames(dds)[inds]
    missing <- setdiff(genes, found)
    warning(
      "Dataset does not contain gene: ", paste(missing, collapse = ", "), "\n"
    )
  }
  inds
}

#' Get expression of gene by indices in a dds
#'
#' @param dds A DESeqDataSet
#' @param gene_indices integer vector that refer to the row number of genes to
#'   retrieve
#' @param assay integer, assay index that refers to the desired expression data
#'   slot
#'
#' @return A dataframe where row names are sample names, and column names are
#'   gene names
#'
#' @details Importantly, returns data as a `data.frame` without dropping
#'   dimensions.
get_gene_expression <- function(dds, gene_indices, assay) {
  SummarizedExperiment::assay(dds, assay)[gene_indices, , drop = FALSE] |>
    t()
}

#' ggplot theme in the style of Edward Tufte
#'
#' @param font_size numeric. Size of `element_text` font.
#' @param use_gillsans logical. Should Gill Sans be used for the plot?
#'
#' @return A ggplot theme
#' @export
#'
#' @examples
#' library(ggplot2)
#'
#' ggplot(dummy_rna_conc, aes(x = sample, y = conc)) +
#'   geom_point() +
#'   theme_tufte(10, use_gillsans = FALSE)
theme_tufte <- function(font_size = 30, use_gillsans = TRUE) {
  font <- "sans"

  if (use_gillsans) {
    if (Sys.info()[["sysname"]] == "Windows") {
      font <- "Gill Sans MT"
    } else if (Sys.info()[["sysname"]] %in% c("Darwin", "Linux")) {
      font <- "Gill Sans"
    }
  }

  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    panel.background = ggplot2::element_rect(
      fill = "#FFFFF8", color = "#CCCCCC"
    ),
    plot.background = ggplot2::element_rect(fill = "#FFFFF8"),
    strip.background = ggplot2::element_rect(fill = "#BBBBB0"),
    legend.background = ggplot2::element_rect(fill = "#FFFFF8"),
    legend.position = "top",
    legend.key = ggplot2::element_blank(),
    text = ggplot2::element_text(family = font, size = font_size)
  )
}

#' Calculate dilution from known concentrations
#'
#' @param c1 Numeric. Initial concentration of sample.
#' @param c2 Numeric. Target concentration of sample.
#' @param v2 Numeric. Target final volume of sample. If `round_for_pipettes =
#'   TRUE`, assumes volume is mL.
#' @param round_for_pipettes Logical. If TRUE, rounds values to the accuracy of
#'   standard pipettes using `make_pipette_vol`.
#' @param quiet Logical. If FALSE, will warn when dilution is impossible to do
#'   without concentrating sample.
#'
#' @return a named list, with `sample_to_add` as the volume of sample to add,
#'   and `add_to` as the volume to dilute the sample into.
#' @export
#'
#' @examples
#' dilute(203, 70, 10)
#' dilute(203, 70, 10, round_for_pipettes = FALSE)
#'
dilute <- function(c1,
                   c2 = min(c1),
                   v2,
                   round_for_pipettes = TRUE,
                   quiet = FALSE) {
  if (any(c2 > c1) & !quiet) {
    warning("This dilution is impossible without concentrating the sample.")
  }

  v1 <- sapply(c1, \(x) c2 * v2 / x)
  add_to <- v2 - v1

  if (round_for_pipettes == TRUE) {
    v1 <- sapply(v1, make_pipette_vol)
    add_to <- sapply(add_to, make_pipette_vol)
  }

  list(sample_to_add = v1, add_to = add_to)
}

#' Round volume to be pipette-compatible
#'
#' @param vol Numeric. Volume to be rounded
#'
#' @return Numeric. Rounded volume.
#' @export
#'
#' @examples
#' make_pipette_vol(104.13398)
#' make_pipette_vol(15.3331)
#' make_pipette_vol(9.9211)
#'
make_pipette_vol <- function(vol) {
  if (is.na(vol)) {
    return(vol)
  } else if (abs(vol) > 200) {
    vol <- round(vol)
  } else if (abs(vol) > 20) {
    vol <- round(vol / 2, 1) * 2
  } else if (abs(vol) > 10) {
    vol <- round(vol / 2, 2) * 2
  } else {
    vol <- round(vol, 2)
  }
  vol
}


#' Get a directory or file from the GBCI SharePoint
#'
#' This function will call either `get_gbci_file` or `get_gbci_dir` downstream
#' depending on the path argument.
#'
#' @param path Path to file or directory on SharePoint
#' @param dest Where to put the file, and what to name it. Defaults to
#'   tempfile/dir
#' @param overwrite Logical. Should files be overwritten if they aready exist?
#' @param drive Optional. A `Microsoft365R::ms_drive`. Can be passed from parent
#'   functions to avoid multiple calls, which can be faster.
#' @param create_dir Logical. Should the destination directory be created if it
#'   doesn't exist?
#'
#' @return The destination of the file or directory (`dest`)
#' @export
get_gbci <- function(path, dest = NULL, overwrite = FALSE, create_dir = TRUE,
                     drive = NULL) {
  if (is.null(drive)) drive <- get_gbci_drive_connection()

  if (is.null(drive$get_item_properties(path)$folder)) {
    get_gbci_file(
      path = path, dest = dest, overwrite = overwrite, drive = drive
    )
  } else {
    get_gbci_dir(
      path = path, dest = dest, overwrite = overwrite, drive = drive,
      create_dir = create_dir
    )
  }
}

#' @rdname get_gbci
#' @export
get_gbci_file <- function(path, dest = NULL, overwrite = FALSE,
                          create_dir = TRUE, drive = NULL) {
  ext <- fs::path_ext(path)
  if (is.null(drive)) drive <- get_gbci_drive_connection()
  if (is.null(dest)) dest <- fs::file_temp(ext = ext)
  if (create_dir) {
    dir.create(fs::path_dir(dest), recursive = TRUE, showWarnings = FALSE)
  }
  drive$download_file(path, dest = dest, overwrite = overwrite)
  dest
}

#' @rdname get_gbci
#' @examples
#' \dontrun{
#' get_gbci_dir("Raw Data/SPECTRAmax/aragaki-kai/", "path/to/my/dir")
#' }
#' @export
get_gbci_dir <- function(path, dest = NULL, overwrite = FALSE,
                         create_dir = TRUE, drive = NULL) {
  if (is.null(drive)) drive <- get_gbci_drive_connection()

  if (is.null(drive$get_item_properties(path)$folder)) {
    stop("Specified path is not a directory")
  }

  if (is.null(dest)) dest <- tempdir()
  items <- drive$list_items(path, full_names = TRUE)
  if (create_dir) dir.create(dest, recursive = TRUE, showWarnings = FALSE)

  apply(
    items, 1, get_recursive,
    drive = drive, og_path = path, dest = dest, overwrite = overwrite,
    simplify = FALSE
  )

  fs::path(dest, stringr::str_extract(path, "[^/]*$"))
}
# FIXME return last dir path not parent to dir path
get_recursive <- function(item, drive, og_path, dest, overwrite) {
  top_dir <- stringr::str_extract(og_path, "[^/]*$")
  fs::dir_create(fs::path(dest, top_dir))
  save_path <- stringr::str_remove(item[["name"]], paste0(og_path, "?/"))
  save_path <- fs::path(dest, top_dir, save_path)
  if (as.logical(item[["isdir"]])) {
    # FIXME redundant to above?
    dir.create(save_path, recursive = TRUE)
    items <- drive$list_items(item[["name"]], full_names = TRUE)
    apply(items, 1, get_recursive, drive, og_path, dest, simplify = FALSE)
  } else {
    drive$download_file(item[["name"]], dest = save_path, overwrite = overwrite)
  }
}

#' @export
list_gbci_dir <- function(path, recursive = FALSE) {
  sp <- Microsoft365R::get_sharepoint_site(
    site_url = "https://livejohnshopkins.sharepoint.com/sites/GBCIStorage"
  )
  drive <- sp$get_drive()
  if (is.null(drive$get_item_properties(path)$folder)) {
    stop("Specified path is not a directory")
  }
  items <- drive$list_items(path, full_names = TRUE)
  if (recursive) {
    items <- apply(items, 1, list_recursive, drive = drive, simplify = FALSE) |>
      lapply(dplyr::bind_rows) |>
      dplyr::bind_rows() # This feels...fragile
  }
  items
}

list_recursive <- function(item, drive) {
  if (as.logical(item[["isdir"]])) {
    items <- drive$list_items(item[["name"]], full_names = TRUE)
    apply(items, 1, list_recursive, drive, simplify = FALSE)
  } else {
    tibble::as_tibble(t(item))
  }
}

#' @export
get_gbci_drive_connection <- function() {
  sp <- Microsoft365R::get_sharepoint_site(
    site_url = "https://livejohnshopkins.sharepoint.com/sites/GBCIStorage"
  )
  sp$get_drive()
}
