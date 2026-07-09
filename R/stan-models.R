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

#' Names of the Stan models shipped with ClassTopics
#' @keywords internal
.stan_model_names <- c("model_betadir", "model_pureNMF", "model_test")

#' Locate the installed .stan source file for a given model
#'
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

#' Directory where compiled Stan executables are cached for this user
#'
#' Uses tools::R_user_dir(), the standard writable, per-user, per-package
#' cache location (respects XDG conventions on Linux/macOS and
#' %LOCALAPPDATA% on Windows). Falls back to tempdir() if unavailable.
#'
#' @return Absolute path to the cache directory (created if needed)
#' @keywords internal
.stan_cache_dir <- function() {

  cache_dir <- tryCatch(
    tools::R_user_dir("ClassTopics", which = "cache"),
    error = function(e) file.path(tempdir(), "ClassTopics_stan_cache")
  )

  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }

  cache_dir
}

#' Get a compiled cmdstanr model, compiling (and caching) on first use
#'
#' This is the single entry point every fitting function should call instead
#' of hardcoding "stan_files/model_xxx.exe" / ".stan" paths. It:
#'   1. Finds the shipped .stan source via system.file()
#'   2. Checks the user's cache dir for an up-to-date compiled executable
#'   3. Compiles (once) if missing or stale, writing the exe into the cache
#'   4. Returns a ready-to-use CmdStanModel object
#'
#' @param model_name One of "model_betadir", "model_pureNMF", "model_test"
#' @param quiet Logical, passed to cmdstanr's compile() to suppress
#'   compiler output (default TRUE)
#' @param force_recompile Logical, force recompilation even if a cached
#'   executable already exists (default FALSE)
#'
#' @return A cmdstanr::CmdStanModel object
#' @keywords internal
.get_stan_model <- function(model_name,
                             quiet = TRUE,
                             force_recompile = FALSE) {

  model_name <- match.arg(model_name, .stan_model_names)

  stan_file <- .stan_source_path(model_name)
  cache_dir <- .stan_cache_dir()

  # Platform-appropriate executable extension: ".exe" on Windows, "" elsewhere.
  # (cmdstanr does not export a helper for this, so we replicate its own
  # internal convention using base R.)
  exe_ext  <- if (.Platform$OS.type == "windows") ".exe" else ""
  exe_file <- file.path(cache_dir, paste0(model_name, exe_ext))

  needs_compile <- force_recompile ||
    !file.exists(exe_file) ||
    file.mtime(stan_file) > file.mtime(exe_file)

  if (needs_compile) {
    message(sprintf(
      "Compiling Stan model '%s' (first use on this machine)...",
      model_name
    ))
    mod <- cmdstanr::cmdstan_model(
      stan_file = stan_file,
      dir       = cache_dir,      # write the compiled binary into the cache
      exe_file  = exe_file,
      quiet     = quiet
    )
  } else {
    mod <- cmdstanr::cmdstan_model(
      stan_file = stan_file,
      exe_file  = exe_file,
      compile   = FALSE
    )
  }

  mod
}
