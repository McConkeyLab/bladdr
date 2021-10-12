#' Add leading 0 to sample names
#'
#' @param sample_names a vector of sample names
#' @return A sample name with up to one zero padded
#' @keywords internal
#' @export

pad_zero <- function(sample_names) {
  new_names <- vector(length = length(sample_names))
  for(i in 1:length(sample_names)){
    if (grepl("^Sample \\d{1,2}$", sample_names[i])){
      new_names[i] <- paste("Sample", formatC(as.integer(regmatches(sample_names[i], regexpr("(?<=Sample )\\d{1,2}", sample_names[i], perl = T))), width = 2, format = "d", flag = "0"))
    }
    else {
      new_names[i] <- sample_names[i]
    }
  }
  new_names
}
