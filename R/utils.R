#' Get indices of genes by name in a dds
#' @noRd
get_gene_index <- function(dds, genes) {
  inds <- which(rownames(dds) %in% genes)
  if(identical(inds, integer(0))) {
    stop("None of the supplied genes exist in this dataset.")
  } else if (length(genes) != length(inds)) {
    found <- rownames(dds)[inds]
    missing <- setdiff(genes, found)
    warning("Dataset does not contain gene: ", paste(missing, collapse = ", "), "\n")
  }
  inds
}

#' Get expression of gene by indices in a dds
#'
#' @param dds A DESeqDataSet
#' @param gene_indices integer vector that refer to the row number of genes to retrieve
#' @param assay integer, assay index that refers to the desired expression data slot
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
    panel.background = ggplot2::element_rect(fill = "#FFFFF8", color = "#CCCCCC"),
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
#' @param c1 numeric. Initial concentration of sample.
#' @param c2 numeric. Target concentration of sample.
#' @param v2 numeric. Target final volume of sample. If `round_for_pipettes = TRUE`, assumes volume is mL.
#' @param round_for_pipettes logical. If TRUE, rounds values to the accuracy of standard pipettes using `make_pipette_vol`.
#' @param quiet logical. If FALSE, will warn when dilution is impossible to do without concentrating sample.
#'
#' @return a named list, with `sample_to_add` as the volume of sample to add, and `add_to` as the volume to dilute the sample into.
#' @export
#'
#' @examples
#' dilute(203, 70, 10)
#' dilute(203, 70, 10, round_for_pipettes = FALSE)
#'
dilute <- function(c1, c2 = min(c1), v2, round_for_pipettes = TRUE, quiet = FALSE) {
  if (any(c2 > c1) & !quiet) {
    warning("This dilution is impossible without concentrating the sample.")
  }

  v1 <- sapply(c1, \(x) c2 * v2 / x)
  add_to <- v2-v1

  if (round_for_pipettes == TRUE) {
    v1 <- sapply(v1, make_pipette_vol)
    add_to <- sapply(add_to, make_pipette_vol)
  }

  list(sample_to_add = v1, add_to = add_to)
}

#' Round volume to be pipette-compatible
#'
#' @param vol numeric. Volume to be rounded
#'
#' @return numeric. Rounded volume.
#' @export
#'
#' @examples
#' make_pipette_vol(104.13398)
#' make_pipette_vol(15.3331)
#' make_pipette_vol(9.9211)
#'
make_pipette_vol <- function(vol) {
  if (is.na(vol)) return(vol)
  if (abs(vol) > 200) {
    vol <- round(vol)
  } else if (abs(vol) > 20) {
    vol <- round(vol/2, 1) * 2
  } else if (abs(vol) > 10) {
    vol <- round(vol/2, 2) * 2
  } else {
    vol <- round(vol, 2)
  }
  vol
}

#' Get a file from the GBCI SharePoint
#'
#' @details This function requires access to the SharePoint in the first place.
#'
#' @param path Path to file on SharePoint
#' @param dest Where to put the file (and what to name it). Defaults to a temp file.
#'
#' @return Character. The local path to the downloaded file.
#' @export
get_gbci_file <- function(path, dest = NULL) {
  ext <- fs::path_ext(path)
  sp <- Microsoft365R::get_sharepoint_site(site_url = "https://livejohnshopkins.sharepoint.com/sites/GBCIStorage")
  drive <- sp$get_drive()

  if (is.null(dest)) {
    dest <- fs::file_temp(ext = ext)
  }

  drive$download_file(path, dest = dest)
  dest
}

#' Recursively download a directory from SharePoint
#'
#' @param path Path to the directory on SharePoint
#' @param dest Path to where the file should be downloaded
#'
#' @return Returns `dest`
#' @export
#'
#' @examples
#' \dontrun{
#' get_gbci_dir("Raw Data/SPECTRAmax/aragaki-kai/", "path/to/my/dir")
#' }
get_gbci_dir <- function(path, dest) {
  sp <- Microsoft365R::get_sharepoint_site(site_url = "https://livejohnshopkins.sharepoint.com/sites/GBCIStorage")
  drive <- sp$get_drive()
  if (is.null(drive$get_item_properties(path)$folder)) {
    stop("Specified path is not a directory")
  }
  items <- drive$list_items(path, full_names = TRUE)
  dir.create(dest, recursive = TRUE)
  apply(items, 1, get_recursive, drive = drive, og_path = path, dest = dest, simplify = FALSE)
  dest
}

get_recursive <- function(item, drive, og_path, dest) {
  top_dir <- stringr::str_extract(og_path, "[^/]*$")
  save_path <- stringr::str_remove(item[["name"]], paste0(og_path, "?/"))
  save_path <- fs::path(dest, top_dir, save_path)
  if (as.logical(item[["isdir"]])) {
    dir.create(save_path, recursive = TRUE)
    items <- drive$list_items(item[["name"]], full_names = TRUE)
    apply(items, 1, get_recursive, drive, og_path, dest, simplify = FALSE)
  } else {
    drive$download_file(item[["name"]], dest = save_path)
  }
}

#' @export
list_gbci_dir <- function(path, dest) {
  sp <- Microsoft365R::get_sharepoint_site(site_url = "https://livejohnshopkins.sharepoint.com/sites/GBCIStorage")
  drive <- sp$get_drive()
  if (is.null(drive$get_item_properties(path)$folder)) {
    stop("Specified path is not a directory")
  }
  items <- drive$list_items(path, full_names = TRUE)
  apply(items, 1, list_recursive, drive = drive, simplify = FALSE) |>
    dplyr::bind_rows()
}

list_recursive <- function(item, drive) {
  if (as.logical(item[["isdir"]])) {
    items <- drive$list_items(item[["name"]], full_names = TRUE)
    apply(items, 1, list_recursive, drive, simplify = FALSE)
  } else {
    item |> t() |> tibble::as_tibble()
  }
}

