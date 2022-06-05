#' Make a volcano plot from a dds results
#'
#' @param results a `DESeqResults` object
#' @param pval the pvalue cutoff to consider significant
#' @param lfc the log fold change cutoff to consider significant
#'
#' @return a volcano plot
#' @export
#'
#' @importFrom rlang .data
#'
#' @examples
#'
#' DESeq2::makeExampleDESeqDataSet(n = 100, m = 5, betaSD = 1) |>
#' DESeq2::DESeq() |>
#' DESeq2::results() |>
#' volcano_plot(lfc = 1)
#'
volcano_plot <- function(results, pval = 0.05, lfc = 0) {

        volcano_data <- data.frame(Gene = rownames(results),
                                   LFC = as.vector(results$log2FoldChange),
                                   p_val = results$padj) |>
                dplyr::filter(stats::complete.cases(.)) |>
                dplyr::mutate(neg_log_p = -log10(.data$p_val),
                              color = dplyr::case_when(
                                      .data$p_val < pval & abs(LFC) > lfc ~ "both",
                                      .data$p_val < pval ~ "pval",
                                      abs(.data$LFC) > lfc ~ "lfc",
                                      TRUE ~ "neither"),
                              is_sig = dplyr::if_else(.data$p_val < pval, "Significant", "Non-Significant"),
                              lfc_cutoff = dplyr::if_else(abs(.data$LFC) > lfc, "Met", "Unmet"))

        ggplot2::ggplot(volcano_data, ggplot2::aes(x = .data$LFC, y = .data$neg_log_p, color = .data$color)) +
                ggplot2::geom_point(size = 2, alpha = 1) +
                ggplot2::ylab("-log10(P-adj)") +
                ggplot2::xlab("Log2 Fold Change") +
                ggplot2::geom_hline(yintercept = -log10(pval), linetype = "dashed") +
                ggplot2::geom_vline(xintercept =  lfc, linetype = "dashed") +
                ggplot2::geom_vline(xintercept = -lfc, linetype = "dashed") +
                ggplot2::theme(legend.position = "none")

}
