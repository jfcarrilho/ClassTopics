# =============================================================================
# Stan model resolution and compilation
# =============================================================================
#
# Design:
#   - The .stan SOURCE files ship inside the package at inst/stan/*.stan.
#     Once installed, they live at system.file("stan", ..., package = "ClassTopics").
#     inst/ is read-only after installation, so we never write compiled
#     binaries there.
#   - The compiled .exe is built ONCE per user machine into a writable cache
#     directory (tools::R_user_dir("ClassTopics", "cache")), and reused on
#     every subsequent call. This is what "compile on first use" means here.

# Names of the Stan models shipped with ClassTopics
#' @keywords internal
.stan_model_names <- c("model_betadir", "model_pureNMF", "model_test")

# Locate the installed .stan source file for a given model
#
#' @param model_name One of "model_betadir", "model_pureNMF", "model_test"
#' @return Absolute path to the installed .stan file
#' @keywords internal
.stan_source_path <- function(model_name) {

  model_name <- match.arg(model_name, .stan_model_names)

  path <- system.file(
    "stan", paste0(model_name, ".stan"),
    package = "ClassTopics",
    mustWork = FALSE
  )

  if (!nzchar(path) || !file.exists(path)) {
    stop(
      sprintf(
        "Could not find '%s.stan' inside the installed package. ",
        model_name
      ),
      "This usually means ClassTopics was not installed correctly ",
      "(inst/stan files missing). Try reinstalling the package.",
      call. = FALSE
    )
  }

  path
}
