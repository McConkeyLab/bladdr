spectra_tidy <- function(file_path = NULL) {

        if (is.null(file_path)) {
                file_path <- file.choose()
        }

        dat <- readxl::read_excel(file_path)
}
