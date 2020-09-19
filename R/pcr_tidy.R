#' Tidy a PCR Excel File
#'
#' Takes in a fresh results output from qPCR and converts it into a tidy
#' format. Useful for downstream analyses.
#'
#' @param file_path Path to an excel file containing the results of a qPCR
#'   run. If left blank, will open up interactive file chooser.
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

pcr_tidy <- function(file_path = NULL) {

        if (is.null(file_path)) {
                file_path <- file.choose()
        }

        dat_og <- readxl::read_excel(file_path, sheet = "Results")

        exp_type <- dat_og[which(dat_og$`Block Type`=="Experiment Type"),2]


        ind_start <- which(dat_og[,1] == "Well")

        if (exp_type == "Standard Curve") {
                ind_end <- nrow(dat_og)
                dat <- dat_og[ind_start:ind_end,]
        }
        if (exp_type == "Comparative Cт (ΔΔCт) ") {
                ind_end <- which(dat_og[,1] == "Analysis Type")
                exp_dat <- dat_og[ind_end:nrow(dat_og), 1:2] %>%
                        t()
                colnames(exp_dat) <- make.names(exp_dat[1,])
                exp_dat <- exp_dat[-1,]
                exp_dat <- t(exp_dat) %>%
                        as.data.frame()
                dat <- dat_og[-c(1:(ind_start-1), (ind_end-1):nrow(dat_og)),]
        }

        dat <- dat_og[-c(1:(ind_start-1), (ind_end-1):nrow(dat_og)),]
        names <- gsub(" ", "_", dat[1,])
        names <- gsub("-", "_", names)
        names <- gsub("\\(superscript_2\\)", "2", names)
        names <- tolower(names)
        colnames(dat) <- names

        dat <- dat[-1,] %>%
                dplyr::mutate(dplyr::across(dplyr::matches("^(delta )*ct.*|^rq|quantity|^baseline|y_intercept|r2|slope|efficiency"), as.numeric),
                              well_row = stringr::str_extract(.data$well_position, "^.{1}"),
                              well_col = as.numeric(stringr::str_extract(.data$well_position, "[:digit:]{1,2}$")),
                              well_row = as.numeric(factor(.data$well_row, levels = LETTERS)))

        if (exp_type == "Comparative Cт (ΔΔCт) ") {
                dat$analysis_type <- exp_dat$Analysis.Type
                dat$control <- exp_dat$Endogenous.Control
                dat$conf_int <- exp_dat$RQ.Min.Max.Confidence.Level
                dat$ref_samp <- exp_dat$Reference.Sample
        }

        dat$plate_type <- colnames(dat_og)[2]
        dat$exp_type <- exp_type[[1]]

        dat
}
