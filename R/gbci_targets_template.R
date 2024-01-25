make_make <- function(project) {
  if (!file.exists("make.R")) {
    sprintf(
      c("library(targets)",
        "",
        "Sys.setenv(TAR_PROJECT = \"%s\")",
        "tar_make()"),
      project
    ) |>
      write("make.R")
  } else {
    sprintf(
      c("Sys.setenv(TAR_PROJECT = \"%s\")",
        "tar_make()"),
      project
    ) |>
      write("make.R", append = TRUE)
  }
}

make_skeleton <- function() {
  fs::dir_create(fs::path("R", "functions"))
  fs::dir_create(fs::path("R", "targets"))
  fs::dir_create("stores")
}

#' Add a new targets 'project'
#'
#' @param project Character. The name of the project. Do not include .R.
#'
#' @details
#' This function will produce the following files:
#'
#' ./R/targets/project.R
#' ./R/functions/project.R
#'
#' It will also update () the make.R script to include this new project, and add
#' it to the `_targets.yaml` file via tar_config_set. Note that if not already
#' created, running this script will create another directory, `stores`, which
#' will contain another subdirectory of the name `project`.
#'
#' @export
add_project <- function(project) {
  make_skeleton()
  make_make(project)
  targets::tar_config_set(
    script = fs::path("R", "targets", project, ext = "R"),
    store = fs::path("stores", project),
    project = project
  )
  fs::file_create(fs::path(".", "R", "functions", project, ext = "R"))
  fs::file_create(fs::path(".", "R", "targets", project, ext = "R"))

  sprintf(
    c("library(targets)",
      "library(tarchetypes)",
      "",
      "tar_option_set(",
      "  packages = c(",
      "    \"conflicted\" # MUST be loaded last",
      "  )",
      ")",
      "",
      "source(\"./R/functions/%s.R\")",
      "",
      "list(",
      "",
      ")"),
    project
  ) |>
    write(fs::path("R", "targets", project, ext = "R"))
}
