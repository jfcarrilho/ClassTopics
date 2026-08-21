#' ClassTopics: Supervised Topic Models for Binary or Categorical Response Prediction via Stan
#'
#' Fits supervised non-negative matrix factorization / topic models to count
#' data jointly with a binary or categorical response, using Stan (via cmdstanr) for
#' full Bayesian inference, with an EM warm-start for initialization,
#' cross-validation utilities, prediction on new samples, and diagnostic
#' plots for topic-response associations.
#'
#' @section Stan models:
#' The package ships three Stan model sources under \code{inst/stan/}:
#' \itemize{
#'   \item \code{model_betadir.stan} -- supervised topic modeling with a Beta-Dirichlet prior
#'   \item \code{model_pureNMF.stan} -- supervised topic modeling without that prior
#'   \item \code{model_test.stan} -- model used to fit topic proportions on held-out data
#' }
#' These are compiled on first use into a per-user cache directory
#' (see \code{tools::R_user_dir("ClassTopics", "cache")}) and reused on
#' subsequent calls. See \code{\link{ClassTopics}} for the main entry point.
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom methods new setClass setGeneric setMethod setValidity representation show
#' @importFrom stats quantile rgamma sd setNames
## usethis namespace: end
NULL
