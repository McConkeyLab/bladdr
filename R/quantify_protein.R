library(tidyverse)
library(parsnip)
library(readxl)
library(shiny)
library(broom)

#' Quantify protein concentration from a MicroBCA assay
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
qp <- function(x, ...) {

}

#' @export
#' @includeRmd quantify_protein
qp.character <- function(x,
                         replicate_orientation = c("h", "v"),
                         target_conc  = 1.0,
                         target_vol   = 15,
                         sample_names = NULL,
                         remove_empty = TRUE) {

  replicate_orientation <- rlang::arg_match(replicate_orientation)

  # Check if SpectraMAX File

  abs <- read_spectramax_excel(x)

  # Theoretically I'd like to be able to work with it as a gp here

  abs_tidy <- tidy_absorbances(absrb, replicate_orientation)

  abs_annot <- annotate_absorbances(tidy_abs)



  list(absorbances = abs_annot)

}

qp_report <- function(qp) {
  # Create report
}



find_mean <- function(df){
  sample_hclust <- df$absrb |>
    dist() |>
    hclust()
  suspect_index <- sample_hclust$merge[nrow(df) - 1, 1] |>
    abs()
  df$is_suspect <- FALSE
  df$is_suspect[suspect_index] <- TRUE
  tibble(absrb = df$absrb, is_suspect = df$is_suspect)
}

read_spectramax_excel <- function(x) {
  readxl::read_excel(x) |>
    dplyr::select(-1) |>
    rlang::set_names(NULL)
}

tidy_absorbances <- function(absorbances, replicate_orientation) {
  if (replicate_orientation == "v") {
    tidy_abs <- absrb |>
      t()
    tidy_abs <- tidy_abs[,1:6]
    tidy_abs <- tidy_abs |>
      tibble::as_tibble(.name_repair = "minimal") |>
      rlang::set_names(c("1", "2", "3", "1", "2", "3"))
    tidy_abs <- dplyr::bind_rows(tidy_abs[,1:3], tidy_abs[,4:6])
  } else {
    abs <- rlang::set_names(absrb, rep(1:3, times = 4))
    tidy_abs <- dplyr::bind_rows(abs[, 1:3], abs[, 4:6], abs[, 7:9], abs[, 10:12])
  }
}

annotate_absorbances <- function(tidy_absorbances) {
  tidy_abs <- tidy_abs |>
    mutate(sample_type = c(rep("standard", times = 7), rep("sample", times = nrow(tidy_abs) - 7)),
           id = row_number()) |>
    pivot_longer(cols = -c(sample_type, id), names_to = "replicate", values_to = "absrb")

  tidy_abs <- tidy_abs |>
    mutate(conc = c(rep(c(0, 0.125, 0.250, 0.500, 1.000, 2.000, 4.000), each = 3),
                    rep(NA, times = nrow(tidy_abs) - 21))) |>
    group_by(id, sample_type, conc) |>
    group_nest() |>
    mutate(mean_no_outlier = map(data, find_mean)) |>
    select(-data) |>
    unnest(cols = c(mean_no_outlier)) |>
    group_by(id) |>
    mutate(no_out_mean = mean(absrb[!.data$is_suspect]),
           no_out_sd = sd(absrb[!.data$is_suspect]),
           keep = if_else(abs(absrb - no_out_mean) > 3 * no_out_sd, FALSE, TRUE),
           log_conc = log2(conc + .5),
           log_abs = log2(absrb))
}
}

make_qp_heatmap <- function(qp) {
  qp$absorbances |>
    set_names(as.character(1:12)) |>
    as_tibble() |>
    mutate(row = LETTERS[1:8]) |>
    pivot_longer(cols = -row, names_to = "column", values_to = "absorbance")
}

server <- function(input, output) {



  dat_fit <- eventReactive({
    input$file
    input$replicate_orientation
  },{
    standards <- filter(dat_tidy(), sample_type == "standard", keep)

    lm_fit <- lm(log_conc ~ log_abs, data = standards)

    tidy_abs <- dat_tidy() |>
      bind_cols(.pred = predict(lm_fit, dat_tidy())) |>
      group_by(id) |>
      mutate(true_mean = mean(.pred[.data$keep])) |>
      mutate(pred_conc = (2^true_mean) - .5)
  })

  tidy_fit <- eventReactive({
    input$file
    input$replicate_orientation
  },{
    standards <- filter(dat_tidy(), sample_type == "standard", keep)

    tidy_fit <- lm(log_conc ~ log_abs, data = standards) |>
      tidy()
  })

  fit_data_filt <- eventReactive({
    input$file
    input$remove_zero
    input$replicate_orientation
  },{
    if(input$remove_zero) {
      filter(dat_fit(), pred_conc > 0 | sample_type == "standard")
    } else{
      dat_fit()
    }
  })

  fit_data_named <- eventReactive({
    input$file
    input$remove_zero
    input$replicate_orientation
    input$sample_name
  },{
    name_list <- c("stnd_0.000", "stnd_0.125", "stnd_0.250", "stnd_0.500",
                   "stnd_1.000", "stnd_2.000", "stnd_4.000",
                   unlist(strsplit(input$sample_name, split = "; ", )))
    with_names <- fit_data_filt() |>
      mutate(sample_name = name_list[id],
             sample_name = if_else(is.na(sample_name), as.character(cur_group_id()), sample_name))
  })

  standards <- eventReactive({
    input$file
    input$remove_zero
    input$replicate_orientation
  },{
    fit_data_named() |>
      filter(sample_type == "standard")
  })

  samples <- eventReactive({fit_data_named()},{
    fit_data_named() |>
      filter(sample_type == "sample")
  })

  sample_summary <- eventReactive({
    input$file
    input$remove_zero
    input$replicate_orientation
    input$target_conc
    input$target_vol
    input$sample_name
  },{
    samples() |>
      ungroup() |>
      select(sample_name, pred_conc, id) |>
      group_by(sample_name, id) |>
      summarize(`Concentration (mg/mL)` = mean(pred_conc)) |>
      arrange(id) |>
      mutate(`[Target] (mg/mL)` = input$target_conc,
             `Target Volume (uL)` = input$target_vol,
             `Lysate to Add (uL)` = round(`[Target] (mg/mL)` * `Target Volume (uL)`/`Concentration (mg/mL)`, 1),
             `RIPA to Add (uL)` = as.character(round(`Target Volume (uL)` - `Lysate to Add (uL)`, 1)),
             `Lysate to Add (uL)` = as.character(`Lysate to Add (uL)`))
  }) # The 'as.character' conversion keeps formatter from adding padding 0s on the end


  output$plate_heat <- renderPlot({
    plate_heatmap()
  }, width = 600)

  std_plot <- reactive({
    ggplot(samples(), aes(x = log_abs,
                          y = true_mean,
                          color = sample_type,
                          shape = keep)) +
      scale_color_viridis_d(option = "viridis", end = 0.8) +
      geom_abline(intercept = tidy_fit()$estimate[1],
                  slope = tidy_fit()$estimate[2],
                  size = 2,
                  alpha = 0.2) +
      geom_point(size = 3, alpha = 0.7) +
      geom_point(data = standards(),
                 aes(y = log_conc),
                 size = 3, alpha = 0.7) +
      scale_shape_manual(values = c(4, 16)) +
      labs(x = "Log2(Absorbance)", y = "Log2(Concentration + 0.5)") +
      theme(legend.position = "bottom",
            panel.background = element_blank(),
            axis.ticks = element_blank())
  })


  plate_heatmap <- reactive({
    temp <- absrb()
    colnames(temp) <- sprintf('%0.2d', 1:ncol(temp))
    cbind(`Row` = LETTERS[1:8], temp) |>
      pivot_longer(cols = -`Row`, names_to = "Column", values_to = "abs") |>
      mutate(Row = fct_rev(Row)) |>
      ggplot(aes(Column, Row)) +
      geom_point(aes(color = abs), size = 20) +
      geom_text(aes(label = round(abs, 2))) +
      scale_color_gradient(low = "darkseagreen1", high = "mediumpurple3") +
      theme(legend.position = "none",
            axis.title = element_blank(),
            panel.background = element_blank())
  })

  make_filename <- eventReactive({input$file},{
    file_name <- input$file$name |>
      str_replace( "\\..*$", "_report.html")
  })

  output$get_report <- downloadHandler(
    filename = function(){
      make_filename()
    },
    content = function(file) {
      temp_report = file.path("report.Rmd")
      params = list(file = input$file,
                    sample_orientation = input$replicate_orientation,
                    remove_zero = input$remove_zero,
                    target_vol = input$target_vol,
                    target_conc = input$target_conc,
                    samples = samples(),
                    sample_summary = sample_summary(),
                    standards = standards(),
                    plate_heatmap = plate_heatmap(),
                    std_plot = std_plot())
      rmarkdown::render(temp_report, output_file = file,
                        params = params, envir = new.env(parent = globalenv()))
    }
  )
}

