pcr_tidy <- function(file_path = NULL) {

        if (is.null(file_path)) {
                file_path <- file.choose()
        }

        dat <- readxl::read_excel(file_path, sheet = "Results")

        ind_start <- which(dat[,1] == "Well")
        ind_end   <- which(dat[,1] == "Analysis Type")

        dat <- dat[-c(1:(ind_start-1), (ind_end-1):nrow(dat)),]

        colnames(dat) <- dat[1,]

        dat <- dat[-1,]
}
