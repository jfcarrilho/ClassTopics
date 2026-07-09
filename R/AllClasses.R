# =============================================================================
# ClassTopics S4 class definitions
#
# This file defines the S4 classes used to represent the output of
# ClassTopics's model-fitting, prediction, and cross-validation functions.
#
#   CTparameter   - a single model parameter (theta, beta, or eta) reported
#                    as a posterior mean matrix with a credible interval
#   CTprediction  - output of predict_ClassTopics_stan(): predicted class
#                    labels and predicted class probabilities for new data
#   cvCTprediction- output of cv_ClassTopics(): k-fold cross-validation
#                    accuracy/prediction summary (excluding the final model,
#                    which lives in the "fit" slot of the containing
#                    cvCTresults object)
#   CTresults     - top-level container summarizing a cmdStan fit produced
#                    by ClassTopics_results(): parameter estimates, topic
#                    correlations, interpretable summaries, in-sample
#                    predictions, and metadata
#   cvCTresults   - as CTresults, but for cv_ClassTopics_results(); the
#                    only difference is that "predictions" is replaced by
#                    "cv_predictions", a cvCTprediction object
# =============================================================================


# -----------------------------------------------------------------------
# CTparameter
# -----------------------------------------------------------------------

#' Class CTparameter
#'
#' An S4 class representing a single model parameter (e.g. theta, beta, or
#' eta) summarized as a posterior mean matrix together with lower and upper
#' credible interval bounds. All three slots are always 2-dimensional and
#' share identical dimensions.
#'
#' @slot mean Numeric matrix of posterior means.
#' @slot lower Numeric matrix of lower credible interval bounds.
#' @slot upper Numeric matrix of upper credible interval bounds.
#'
#' @name CTparameter-class
#' @rdname CTparameter-class
#' @exportClass CTparameter
setClass("CTparameter",
  representation(
    mean  = "matrix",
    lower = "matrix",
    upper = "matrix"
  )
)

setValidity("CTparameter", function(object) {
  errs <- character()

  dims <- list(mean = dim(object@mean),
               lower = dim(object@lower),
               upper = dim(object@upper))

  if (!identical(dims$mean, dims$lower) || !identical(dims$mean, dims$upper)) {
    errs <- c(errs, "'mean', 'lower', and 'upper' must have identical dimensions")
  }

  if (length(errs) == 0 &&
      any(object@lower > object@upper, na.rm = TRUE)) {
    errs <- c(errs, "'lower' bound exceeds 'upper' bound for at least one element")
  }

  if (length(errs) == 0) TRUE else errs
})

#' Constructor for CTparameter
#'
#' @param mean Numeric matrix of posterior means.
#' @param lower Numeric matrix of lower credible interval bounds.
#' @param upper Numeric matrix of upper credible interval bounds.
#'
#' @return An object of class \code{CTparameter}.
#' @export
CTparameter <- function(mean, lower, upper){
  new("CTparameter", mean = mean, lower = lower, upper = upper)
}

#' @param x An CTparameter object.
#' @param object An CTparameter object.
#' @rdname CTparameter-class
#' @export
setGeneric("paramMean", function(x) standardGeneric("paramMean"))

#' @rdname CTparameter-class
#' @export
setMethod("paramMean", "CTparameter", function(x) x@mean)

#' @rdname CTparameter-class
#' @export
setGeneric("paramCI", function(x) standardGeneric("paramCI"))

#' @rdname CTparameter-class
#' @export
setMethod("paramCI", "CTparameter", function(x) {
  list(lower = x@lower, upper = x@upper)
})

#' @rdname CTparameter-class
#' @export
setMethod("show", "CTparameter", function(object) {
  cat("CTparameter object\n")
  cat("  Dimensions:", paste(dim(object@mean), collapse = " x "), "\n")
})

# -----------------------------------------------------------------------
# CTtrainpred
# -----------------------------------------------------------------------

#' Class CTtrainpred
#'
#' An S4 class representing the summary of the predictive performance of
#' the fitted multiclass topic model (without cross-validation) on the
#' respective training data.
#'
#' @slot predicted_categories Factor vector of posterior predicted categories.
#' @slot true_categories Factor vector of true categories.
#' @slot prediction_probabilities Numeric matrix whose rows are the posterior
#'   softmax vectors for each observation.
#' @slot overall_accuracy Numeric scalar, accuracy of the model.
#' @slot category_accuracy Numeric vector, accuracy of the model stratified by
#'   class.
#' @slot confusion_matrix Integer matrix, confusion matrix of the counts at
#'   each \code{predicted_categories}/\code{true_categories} combination.
#'
#' @name CTtrainpred-class
#' @rdname CTtrainpred-class
#' @exportClass CTtrainpred
setClass("CTtrainpred",
         representation(
           predicted_categories     = "factor",
           true_categories          = "factor",
           prediction_probabilities = "matrix",
           overall_accuracy         = "numeric",
           category_accuracy        = "numeric",
           confusion_matrix         = "table"
         )
)

setValidity("CTtrainpred", function(object) {
  errs <- character()
  if (length(object@predicted_categories) != length(object@true_categories)) {
    errs <- c(errs, "predicted_categories and true_categories must be the same length")
  }
  if (object@overall_accuracy < 0 || object@overall_accuracy > 1) {
    errs <- c(errs, "overall_accuracy must be in [0, 1]")
  }
  if (length(errs)) errs else TRUE
})

#' Constructor for CTtrainpred
#'
#' @param predicted_categories Factor vector of posterior predicted categories.
#' @param true_categories Factor vector of true categories.
#' @param prediction_probabilities Numeric matrix whose rows are the posterior
#'   softmax vectors for each observation.
#' @param overall_accuracy Numeric scalar, accuracy of the model.
#' @param category_accuracy Numeric vector, accuracy of the model stratified by
#'   class.
#' @param confusion_matrix Integer matrix, confusion matrix of the counts at
#'   each \code{predicted_categories}/\code{true_categories} combination.
#'
#' @return An object of class \code{CTtrainpred}.
#' @export
CTtrainpred <- function(predicted_categories, true_categories,
                         prediction_probabilities, overall_accuracy,
                         category_accuracy, confusion_matrix){
  new("CTtrainpred",
      predicted_categories     = predicted_categories,
      true_categories          = true_categories,
      prediction_probabilities = prediction_probabilities,
      overall_accuracy         = overall_accuracy,
      category_accuracy        = category_accuracy,
      confusion_matrix         = confusion_matrix)
}


# -----------------------------------------------------------------------
# CTprediction  (output of predict_ClassTopics_stan/predict_multiclass_EM)
# -----------------------------------------------------------------------

#' Class CTprediction
#'
#' An S4 class representing the output of \code{\link{predict_ClassTopics_stan}}:
#' predicted class labels and predicted class probabilities for new
#' (out-of-sample) observations, given a fitted model's topic-gene loadings
#' and regression coefficients.
#'
#' @slot predicted_class Character vector of predicted category labels
#'   (one per observation).
#' @slot predicted_probs Numeric matrix of predicted class probabilities
#'   (observations x categories).
#'
#' @name CTprediction-class
#' @rdname CTprediction-class
#' @exportClass CTprediction
setClass("CTprediction",
  representation(
    predicted_class = "character",
    predicted_probs  = "matrix"
  )
)

setValidity("CTprediction", function(object) {
  errs <- character()

  if (length(object@predicted_class) != nrow(object@predicted_probs)) {
    errs <- c(errs,
      "length of 'predicted_class' must equal the number of rows of 'predicted_probs'")
  }

  if (length(errs) == 0) TRUE else errs
})

#' Constructor for CTprediction
#'
#' @param predicted_class Character vector of predicted category labels.
#' @param predicted_probs Numeric matrix of predicted class probabilities.
#'
#' @return An object of class \code{CTprediction}.
#' @export
CTprediction <- function(predicted_class, predicted_probs){
  new("CTprediction",
      predicted_class = predicted_class,
      predicted_probs = predicted_probs)
}

#' @param x A CTprediction object.
#' @param object A CTprediction object.
#' @rdname CTprediction-class
#' @export
setGeneric("predictedClass", function(x) standardGeneric("predictedClass"))

#' @rdname CTprediction-class
#' @export
setMethod("predictedClass", "CTprediction", function(x) x@predicted_class)

#' @rdname CTprediction-class
#' @export
setGeneric("predictedProbs", function(x) standardGeneric("predictedProbs"))

#' @rdname CTprediction-class
#' @export
setMethod("predictedProbs", "CTprediction", function(x) x@predicted_probs)

#' @rdname CTprediction-class
#' @export
setMethod("show", "CTprediction", function(object) {
  cat("CTprediction object\n")
  cat("  N observations:", length(object@predicted_class), "\n")
  cat("  N categories:  ", ncol(object@predicted_probs), "\n")
})

# -----------------------------------------------------------------------
# cvCTprediction  (output of cv_ClassTopics - minus final_model)
# -----------------------------------------------------------------------

#' Class cvCTprediction
#'
#' An S4 class summarizing k-fold cross-validation results produced by
#' \code{\link{cv_ClassTopics}}: per-fold and overall train/test accuracy,
#' per-category accuracy, the test confusion matrix, an overfitting
#' diagnostic, and pooled out-of-fold predictions. This class does not
#' include the final model fitted on the full dataset; that object is
#' stored separately as the \code{fit} slot of the containing
#' \code{\link{cvCTresults-class}} object.
#'
#' @slot cv_train_accuracy_mean Numeric scalar, mean training accuracy
#'   across folds.
#' @slot cv_train_accuracy_sd Numeric scalar, SD of training accuracy
#'   across folds.
#' @slot fold_train_accuracies Numeric vector of per-fold training accuracy.
#' @slot cv_test_accuracy Numeric scalar, overall (pooled) test accuracy.
#' @slot cv_test_accuracy_mean Numeric scalar, mean test accuracy across
#'   folds.
#' @slot cv_test_accuracy_sd Numeric scalar, SD of test accuracy across
#'   folds.
#' @slot fold_test_accuracies Numeric vector of per-fold test accuracy.
#' @slot cv_test_confusion_matrix A \code{table}, pooled test confusion
#'   matrix across all folds.
#' @slot cv_train_class_accuracy_mean Named numeric vector, mean per-category
#'   training accuracy across folds.
#' @slot cv_train_class_accuracy_sd Named numeric vector, SD of per-category
#'   training accuracy across folds.
#' @slot fold_train_class_accuracies Numeric matrix (folds x categories) of
#'   per-fold, per-category training accuracy.
#' @slot cv_test_class_accuracy_mean Named numeric vector, mean per-category
#'   test accuracy across folds.
#' @slot cv_test_class_accuracy_sd Named numeric vector, SD of per-category
#'   test accuracy across folds.
#' @slot fold_test_class_accuracies Numeric matrix (folds x categories) of
#'   per-fold, per-category test accuracy.
#' @slot overfitting_gap Numeric scalar, mean training accuracy minus mean
#'   test accuracy.
#' @slot all_test_predictions Factor of pooled out-of-fold predicted
#'   categories (one per observation, in original data order).
#' @slot all_test_pred_probs Numeric matrix of pooled out-of-fold predicted
#'   class probabilities.
#'
#' @name cvCTprediction-class
#' @rdname cvCTprediction-class
#' @exportClass cvCTprediction
setClass("cvCTprediction",
         representation(
           cv_train_accuracy_mean       = "numeric",
           cv_train_accuracy_sd         = "numeric",
           fold_train_accuracies        = "numeric",
           
           cv_test_accuracy             = "numeric",
           cv_test_accuracy_mean        = "numeric",
           cv_test_accuracy_sd          = "numeric",
           fold_test_accuracies         = "numeric",
           cv_test_confusion_matrix     = "table",
           
           cv_train_class_accuracy_mean = "numeric",
           cv_train_class_accuracy_sd   = "numeric",
           fold_train_class_accuracies  = "matrix",
           
           cv_test_class_accuracy_mean  = "numeric",
           cv_test_class_accuracy_sd    = "numeric",
           fold_test_class_accuracies   = "matrix",
           
           overfitting_gap              = "numeric",
           
           all_test_predictions         = "factor",
           all_test_pred_probs          = "matrix"
         )
)

setValidity("cvCTprediction", function(object) {
  errs <- character()
  
  if (length(object@cv_train_accuracy_mean) != 1 ||
      object@cv_train_accuracy_mean < 0 || object@cv_train_accuracy_mean > 1) {
    errs <- c(errs, "'cv_train_accuracy_mean' must be a single value in [0, 1]")
  }
  
  if (length(object@cv_test_accuracy) != 1 ||
      object@cv_test_accuracy < 0 || object@cv_test_accuracy > 1) {
    errs <- c(errs, "'cv_test_accuracy' must be a single value in [0, 1]")
  }
  
  if (length(object@all_test_predictions) != nrow(object@all_test_pred_probs)) {
    errs <- c(errs,
              "length of 'all_test_predictions' must equal the number of rows of 'all_test_pred_probs'")
  }
  
  if (nrow(object@fold_train_class_accuracies) != length(object@fold_train_accuracies) ||
      nrow(object@fold_test_class_accuracies) != length(object@fold_test_accuracies)) {
    errs <- c(errs,
              "number of folds implied by 'fold_train_class_accuracies'/'fold_test_class_accuracies' must match 'fold_train_accuracies'/'fold_test_accuracies'")
  }
  
  if (length(errs) == 0) TRUE else errs
})

#' Constructor for cvCTprediction
#'
#' @param cv_train_accuracy_mean Numeric scalar, mean training accuracy
#'   across folds.
#' @param cv_train_accuracy_sd Numeric scalar, SD of training accuracy
#'   across folds.
#' @param fold_train_accuracies Numeric vector of per-fold training accuracy.
#' @param cv_test_accuracy Numeric scalar, overall (pooled) test accuracy.
#' @param cv_test_accuracy_mean Numeric scalar, mean test accuracy across
#'   folds.
#' @param cv_test_accuracy_sd Numeric scalar, SD of test accuracy across
#'   folds.
#' @param fold_test_accuracies Numeric vector of per-fold test accuracy.
#' @param cv_test_confusion_matrix A \code{table}, pooled test confusion
#'   matrix.
#' @param cv_train_class_accuracy_mean Named numeric vector, mean
#'   per-category training accuracy.
#' @param cv_train_class_accuracy_sd Named numeric vector, SD of
#'   per-category training accuracy.
#' @param fold_train_class_accuracies Numeric matrix (folds x categories).
#' @param cv_test_class_accuracy_mean Named numeric vector, mean
#'   per-category test accuracy.
#' @param cv_test_class_accuracy_sd Named numeric vector, SD of per-category
#'   test accuracy.
#' @param fold_test_class_accuracies Numeric matrix (folds x categories).
#' @param overfitting_gap Numeric scalar.
#' @param all_test_predictions Factor of pooled out-of-fold predictions.
#' @param all_test_pred_probs Numeric matrix of pooled out-of-fold predicted
#'   probabilities.
#'
#' @return An object of class \code{cvCTprediction}.
#' @export
cvCTprediction <- function(cv_train_accuracy_mean,
                            cv_train_accuracy_sd,
                            fold_train_accuracies,
                            cv_test_accuracy,
                            cv_test_accuracy_mean,
                            cv_test_accuracy_sd,
                            fold_test_accuracies,
                            cv_test_confusion_matrix,
                            cv_train_class_accuracy_mean,
                            cv_train_class_accuracy_sd,
                            fold_train_class_accuracies,
                            cv_test_class_accuracy_mean,
                            cv_test_class_accuracy_sd,
                            fold_test_class_accuracies,
                            overfitting_gap,
                            all_test_predictions,
                            all_test_pred_probs){
  new("cvCTprediction",
      cv_train_accuracy_mean       = cv_train_accuracy_mean,
      cv_train_accuracy_sd         = cv_train_accuracy_sd,
      fold_train_accuracies        = fold_train_accuracies,
      cv_test_accuracy             = cv_test_accuracy,
      cv_test_accuracy_mean        = cv_test_accuracy_mean,
      cv_test_accuracy_sd          = cv_test_accuracy_sd,
      fold_test_accuracies         = fold_test_accuracies,
      cv_test_confusion_matrix     = cv_test_confusion_matrix,
      cv_train_class_accuracy_mean = cv_train_class_accuracy_mean,
      cv_train_class_accuracy_sd   = cv_train_class_accuracy_sd,
      fold_train_class_accuracies  = fold_train_class_accuracies,
      cv_test_class_accuracy_mean  = cv_test_class_accuracy_mean,
      cv_test_class_accuracy_sd    = cv_test_class_accuracy_sd,
      fold_test_class_accuracies   = fold_test_class_accuracies,
      overfitting_gap              = overfitting_gap,
      all_test_predictions         = all_test_predictions,
      all_test_pred_probs          = all_test_pred_probs)
}

#' @param x A cvCTprediction object.
#' @param object A cvCTprediction object.
#' @rdname cvCTprediction-class
#' @export
setGeneric("cvTestAccuracy", function(x) standardGeneric("cvTestAccuracy"))

#' @rdname cvCTprediction-class
#' @export
setMethod("cvTestAccuracy", "cvCTprediction", function(x) x@cv_test_accuracy)

#' @rdname cvCTprediction-class
#' @export
setGeneric("overfittingGap", function(x) standardGeneric("overfittingGap"))

#' @rdname cvCTprediction-class
#' @export
setMethod("overfittingGap", "cvCTprediction", function(x) x@overfitting_gap)

#' @rdname cvCTprediction-class
#' @export
setMethod("show", "cvCTprediction", function(object) {
  cat("cvCTprediction object\n")
  cat("  Folds:               ", length(object@fold_test_accuracies), "\n")
  cat("  Mean train accuracy: ", round(object@cv_train_accuracy_mean, 3),
      " (SD", round(object@cv_train_accuracy_sd, 3), ")\n")
  cat("  Mean test accuracy:  ", round(object@cv_test_accuracy_mean, 3),
      " (SD", round(object@cv_test_accuracy_sd, 3), ")\n")
  cat("  Pooled test accuracy:", round(object@cv_test_accuracy, 3), "\n")
  cat("  Overfitting gap:     ", round(object@overfitting_gap, 3), "\n")
})

# -----------------------------------------------------------------------
# cvCTprediction_plusfm  (output of cv_ClassTopics)
# -----------------------------------------------------------------------

#' Class cvCTprediction_plusfm
#'
#' An S4 class summarizing k-fold cross-validation results produced by
#' \code{\link{cv_ClassTopics}}: per-fold and overall train/test accuracy,
#' per-category accuracy, the test confusion matrix, an overfitting
#' diagnostic, and pooled out-of-fold predictions.
#'
#' @slot folds List of numeric vectors, test folds created (observations are
#'   identified by indexes)
#' @slot cv_train_accuracy_mean Numeric scalar, mean training accuracy
#'   across folds.
#' @slot cv_train_accuracy_sd Numeric scalar, SD of training accuracy
#'   across folds.
#' @slot fold_train_accuracies Numeric vector of per-fold training accuracy.
#' @slot cv_test_accuracy Numeric scalar, overall (pooled) test accuracy.
#' @slot cv_test_accuracy_mean Numeric scalar, mean test accuracy across
#'   folds.
#' @slot cv_test_accuracy_sd Numeric scalar, SD of test accuracy across
#'   folds.
#' @slot fold_test_accuracies Numeric vector of per-fold test accuracy.
#' @slot cv_test_confusion_matrix A \code{table}, pooled test confusion
#'   matrix across all folds.
#' @slot cv_train_class_accuracy_mean Named numeric vector, mean per-category
#'   training accuracy across folds.
#' @slot cv_train_class_accuracy_sd Named numeric vector, SD of per-category
#'   training accuracy across folds.
#' @slot fold_train_class_accuracies Numeric matrix (folds x categories) of
#'   per-fold, per-category training accuracy.
#' @slot cv_test_class_accuracy_mean Named numeric vector, mean per-category
#'   test accuracy across folds.
#' @slot cv_test_class_accuracy_sd Named numeric vector, SD of per-category
#'   test accuracy across folds.
#' @slot fold_test_class_accuracies Numeric matrix (folds x categories) of
#'   per-fold, per-category test accuracy.
#' @slot overfitting_gap Numeric scalar, mean training accuracy minus mean
#'   test accuracy.
#' @slot all_test_predictions Factor of pooled out-of-fold predicted
#'   categories (one per observation, in original data order).
#' @slot all_test_pred_probs Numeric matrix of pooled out-of-fold predicted
#'   class probabilities.
#' @slot final_model The raw model object fitted on the full training dataset
#'   and returned by cmdStan (e.g. a \code{CmdStanMCMC} object).
#'
#' @name cvCTprediction_plusfm-class
#' @rdname cvCTprediction_plusfm-class
#' @exportClass cvCTprediction_plusfm
setClass("cvCTprediction_plusfm",
  representation(
    folds                        = "list", 
    cv_train_accuracy_mean       = "numeric",
    cv_train_accuracy_sd         = "numeric",
    fold_train_accuracies        = "numeric",

    cv_test_accuracy             = "numeric",
    cv_test_accuracy_mean        = "numeric",
    cv_test_accuracy_sd          = "numeric",
    fold_test_accuracies         = "numeric",
    cv_test_confusion_matrix     = "table",

    cv_train_class_accuracy_mean = "numeric",
    cv_train_class_accuracy_sd   = "numeric",
    fold_train_class_accuracies  = "matrix",

    cv_test_class_accuracy_mean  = "numeric",
    cv_test_class_accuracy_sd    = "numeric",
    fold_test_class_accuracies   = "matrix",

    overfitting_gap              = "numeric",

    all_test_predictions         = "factor",
    all_test_pred_probs          = "matrix",
    
    final_model                  = "ANY"
  )
)

setValidity("cvCTprediction_plusfm", function(object) {
  errs <- character()

  if (length(object@cv_train_accuracy_mean) != 1 ||
      object@cv_train_accuracy_mean < 0 || object@cv_train_accuracy_mean > 1) {
    errs <- c(errs, "'cv_train_accuracy_mean' must be a single value in [0, 1]")
  }

  if (length(object@cv_test_accuracy) != 1 ||
      object@cv_test_accuracy < 0 || object@cv_test_accuracy > 1) {
    errs <- c(errs, "'cv_test_accuracy' must be a single value in [0, 1]")
  }

  if (length(object@all_test_predictions) != nrow(object@all_test_pred_probs)) {
    errs <- c(errs,
      "length of 'all_test_predictions' must equal the number of rows of 'all_test_pred_probs'")
  }

  if (nrow(object@fold_train_class_accuracies) != length(object@fold_train_accuracies) ||
      nrow(object@fold_test_class_accuracies) != length(object@fold_test_accuracies)) {
    errs <- c(errs,
      "number of folds implied by 'fold_train_class_accuracies'/'fold_test_class_accuracies' must match 'fold_train_accuracies'/'fold_test_accuracies'")
  }

  if (length(errs) == 0) TRUE else errs
})

#' Constructor for cvCTprediction_plusfm
#'
#' @param folds List of numeric vectors, test folds created
#' @param cv_train_accuracy_mean Numeric scalar, mean training accuracy
#'   across folds.
#' @param cv_train_accuracy_sd Numeric scalar, SD of training accuracy
#'   across folds.
#' @param fold_train_accuracies Numeric vector of per-fold training accuracy.
#' @param cv_test_accuracy Numeric scalar, overall (pooled) test accuracy.
#' @param cv_test_accuracy_mean Numeric scalar, mean test accuracy across
#'   folds.
#' @param cv_test_accuracy_sd Numeric scalar, SD of test accuracy across
#'   folds.
#' @param fold_test_accuracies Numeric vector of per-fold test accuracy.
#' @param cv_test_confusion_matrix A \code{table}, pooled test confusion
#'   matrix.
#' @param cv_train_class_accuracy_mean Named numeric vector, mean
#'   per-category training accuracy.
#' @param cv_train_class_accuracy_sd Named numeric vector, SD of
#'   per-category training accuracy.
#' @param fold_train_class_accuracies Numeric matrix (folds x categories).
#' @param cv_test_class_accuracy_mean Named numeric vector, mean
#'   per-category test accuracy.
#' @param cv_test_class_accuracy_sd Named numeric vector, SD of per-category
#'   test accuracy.
#' @param fold_test_class_accuracies Numeric matrix (folds x categories).
#' @param overfitting_gap Numeric scalar.
#' @param all_test_predictions Factor of pooled out-of-fold predictions.
#' @param all_test_pred_probs Numeric matrix of pooled out-of-fold predicted
#'   probabilities.
#' @param final_model The raw model object fitted on the full training dataset
#'   and returned by cmdStan (e.g. a \code{CmdStanMCMC} object).
#'
#' @return An object of class \code{cvCTprediction_plusfm}.
#' @export
cvCTprediction_plusfm <- function(folds,
                            cv_train_accuracy_mean,
                            cv_train_accuracy_sd,
                            fold_train_accuracies,
                            cv_test_accuracy,
                            cv_test_accuracy_mean,
                            cv_test_accuracy_sd,
                            fold_test_accuracies,
                            cv_test_confusion_matrix,
                            cv_train_class_accuracy_mean,
                            cv_train_class_accuracy_sd,
                            fold_train_class_accuracies,
                            cv_test_class_accuracy_mean,
                            cv_test_class_accuracy_sd,
                            fold_test_class_accuracies,
                            overfitting_gap,
                            all_test_predictions,
                            all_test_pred_probs,
                            final_model){
  new("cvCTprediction_plusfm",
      folds                        = folds,
      cv_train_accuracy_mean       = cv_train_accuracy_mean,
      cv_train_accuracy_sd         = cv_train_accuracy_sd,
      fold_train_accuracies        = fold_train_accuracies,
      cv_test_accuracy             = cv_test_accuracy,
      cv_test_accuracy_mean        = cv_test_accuracy_mean,
      cv_test_accuracy_sd          = cv_test_accuracy_sd,
      fold_test_accuracies         = fold_test_accuracies,
      cv_test_confusion_matrix     = cv_test_confusion_matrix,
      cv_train_class_accuracy_mean = cv_train_class_accuracy_mean,
      cv_train_class_accuracy_sd   = cv_train_class_accuracy_sd,
      fold_train_class_accuracies  = fold_train_class_accuracies,
      cv_test_class_accuracy_mean  = cv_test_class_accuracy_mean,
      cv_test_class_accuracy_sd    = cv_test_class_accuracy_sd,
      fold_test_class_accuracies   = fold_test_class_accuracies,
      overfitting_gap              = overfitting_gap,
      all_test_predictions         = all_test_predictions,
      all_test_pred_probs          = all_test_pred_probs,
      final_model                  = final_model)
}

#' @param x A cvCTprediction_plusfm object.
#' @param object A cvCTprediction_plusfm object.
#' @rdname cvCTprediction_plusfm-class
#' @export
setMethod("cvTestAccuracy", "cvCTprediction_plusfm", function(x) x@cv_test_accuracy)

#' @rdname cvCTprediction_plusfm-class
#' @export
setMethod("overfittingGap", "cvCTprediction_plusfm", function(x) x@overfitting_gap)

#' @rdname cvCTprediction_plusfm-class
#' @export
setMethod("show", "cvCTprediction_plusfm", function(object) {
  cat("cvCTprediction_plusfm object\n")
  cat("  Folds:               ", length(object@fold_test_accuracies), "\n")
  cat("  Mean train accuracy: ", round(object@cv_train_accuracy_mean, 3),
      " (SD", round(object@cv_train_accuracy_sd, 3), ")\n")
  cat("  Mean test accuracy:  ", round(object@cv_test_accuracy_mean, 3),
      " (SD", round(object@cv_test_accuracy_sd, 3), ")\n")
  cat("  Pooled test accuracy:", round(object@cv_test_accuracy, 3), "\n")
  cat("  Overfitting gap:     ", round(object@overfitting_gap, 3), "\n")
})


# -----------------------------------------------------------------------
# CTresults  (output of ClassTopics_results)
# -----------------------------------------------------------------------

#' Class CTresults
#'
#' Top-level S4 class summarizing the output of
#' \code{\link{ClassTopics_results}}. It bundles parameter estimates
#' (theta, beta, eta), topic correlations, interpretable summaries (top
#' genes, regression coefficients), in-sample prediction results, and
#' metadata, alongside the raw cmdStan fit object from which these
#' summaries were derived.
#'
#' @slot theta An \code{CTparameter} object: document-topic proportions.
#' @slot beta An \code{CTparameter} object: topic-gene loadings.
#' @slot eta An \code{CTparameter} object: topic-category regression
#'   parameters.
#' @slot topic_correlations Numeric matrix of correlations between topics.
#' @slot top_genes A named list of \code{data.frame}s, one per topic, each
#'   giving the top genes and their loading estimates (mean/lower/upper).
#' @slot regression_coefficients A \code{data.frame} summary of \code{eta},
#'   melted by (category, topic) pairs, with columns for the posterior mean
#'   and credible interval bounds.
#' @slot predictions An \code{CTprediction} object: in-sample predictions
#'   made by the fitted model.
#' @slot response_levels Character vector of category labels.
#' @slot gene_names Character vector of gene names.
#' @slot credible_interval Numeric scalar giving the credible interval
#'   width used (e.g. \code{0.95}).
#' @slot fit The raw fitted model object returned by cmdStan (e.g. a
#'   \code{CmdStanMCMC} object), of which this object is a summary.
#'
#' @name CTresults-class
#' @rdname CTresults-class
#' @exportClass CTresults
setClass("CTresults",
  representation(
    theta                   = "CTparameter",
    beta                    = "CTparameter",
    eta                     = "CTparameter",
    topic_correlations      = "matrix",
    top_genes               = "list",
    regression_coefficients = "data.frame",
    predictions             = "CTprediction",
    response_levels         = "character",
    gene_names              = "character",
    credible_interval       = "numeric",
    fit                     = "ANY"
  )
)

setValidity("CTresults", function(object) {
  errs <- character()

  if (length(object@credible_interval) != 1 ||
      object@credible_interval <= 0 || object@credible_interval >= 1) {
    errs <- c(errs, "'credible_interval' must be a single value in (0, 1)")
  }

  expected_reg_cols <- c("Category", "Topic", "Coefficient_Mean",
                          "Coefficient_Lower", "Coefficient_Upper")
  if (!all(expected_reg_cols %in% colnames(object@regression_coefficients))) {
    errs <- c(errs,
      paste0("'regression_coefficients' must contain columns: ",
             paste(expected_reg_cols, collapse = ", ")))
  }

  if (length(errs) == 0) TRUE else errs
})

#' Constructor for CTresults
#'
#' @param theta An \code{CTparameter} object.
#' @param beta An \code{CTparameter} object.
#' @param eta An \code{CTparameter} object.
#' @param topic_correlations Numeric matrix of topic correlations.
#' @param top_genes A named list of \code{data.frame}s, one per topic.
#' @param regression_coefficients A \code{data.frame} summary of \code{eta}
#'   melted by (category, topic) pairs.
#' @param predictions An \code{CTprediction} object.
#' @param response_levels Character vector of category labels.
#' @param gene_names Character vector of gene names.
#' @param credible_interval Numeric scalar, credible interval width.
#' @param fit The raw fitted model object from cmdStan.
#'
#' @return An object of class \code{CTresults}.
#' @export
CTresults <- function(theta,
                       beta,
                       eta,
                       topic_correlations,
                       top_genes,
                       regression_coefficients,
                       predictions,
                       response_levels,
                       gene_names,
                       credible_interval,
                       fit) {
  new("CTresults",
      theta                   = theta,
      beta                    = beta,
      eta                     = eta,
      topic_correlations      = topic_correlations,
      top_genes               = top_genes,
      regression_coefficients = regression_coefficients,
      predictions             = predictions,
      response_levels         = response_levels,
      gene_names              = gene_names,
      credible_interval       = credible_interval,
      fit                     = fit)
}

#' @param x A CTresults object.
#' @param ... Additional arguments (currently unused; present for generic
#'   consistency with methods that may accept extra arguments).
#' @param object A CTresults object.
#' @rdname CTresults-class
#' @export
setGeneric("getPredictions", function(x) standardGeneric("getPredictions"))

#' @rdname CTresults-class
#' @export
setMethod("getPredictions", "CTresults", function(x) x@predictions)

#' @rdname CTresults-class
#' @export
setGeneric("topGenes", function(x, ...) standardGeneric("topGenes"))

#' @rdname CTresults-class
#' @export
setMethod("topGenes", "CTresults", function(x, ...) x@top_genes)

#' @rdname CTresults-class
#' @export
setGeneric("regressionCoefficients", function(x) standardGeneric("regressionCoefficients"))

#' @rdname CTresults-class
#' @export
setMethod("regressionCoefficients", "CTresults", function(x) x@regression_coefficients)

#' @rdname CTresults-class
#' @export
setGeneric("getFit", function(x) standardGeneric("getFit"))

#' @rdname CTresults-class
#' @export
setMethod("getFit", "CTresults", function(x) x@fit)

#' @rdname CTresults-class
#' @export
setMethod("show", "CTresults", function(object) {
  cat("CTresults object (ClassTopics)\n")
  cat("  Topics:           ", ncol(object@theta@mean), "\n")
  cat("  Documents:        ", nrow(object@theta@mean), "\n")
  cat("  Genes:            ", length(object@gene_names), "\n")
  cat("  Categories:       ", length(object@response_levels), "\n")
  cat("  Credible interval:", object@credible_interval, "\n")
})


# -----------------------------------------------------------------------
# cvCTresults  (output of cv_ClassTopics_results)
# -----------------------------------------------------------------------

#' Class cvCTresults
#'
#' Top-level S4 class summarizing the output of
#' \code{\link{cv_ClassTopics_results}}. It extends \code{CTresults},
#' inheriting all of its slots (theta, beta, eta, topic correlations, top
#' genes, regression coefficients, metadata, and the final-model \code{fit}
#' object trained on the full dataset), but replaces the \code{predictions}
#' slot with \code{cv_predictions}, a \code{\link{cvCTprediction-class}}
#' object summarizing cross-validated (out-of-fold) performance rather than
#' in-sample predictions.
#'
#' @slot folds List of numeric vectors, test folds created (observations are
#'   identified by indexes)
#' @slot theta An \code{CTparameter} object: document-topic proportions.
#' @slot beta An \code{CTparameter} object: topic-gene loadings.
#' @slot eta An \code{CTparameter} object: topic-category regression
#'   parameters.
#' @slot topic_correlations Numeric matrix of correlations between topics.
#' @slot top_genes A named list of \code{data.frame}s, one per topic, each
#'   giving the top genes and their loading estimates (mean/lower/upper).
#' @slot regression_coefficients A \code{data.frame} summary of \code{eta},
#'   melted by (category, topic) pairs, with columns for the posterior mean
#'   and credible interval bounds.
#' @slot cv_predictions A \code{cvCTprediction} object summarizing k-fold
#'   cross-validation accuracy and pooled out-of-fold predictions.
#' @slot response_levels Character vector of category labels.
#' @slot gene_names Character vector of gene names.
#' @slot credible_interval Numeric scalar giving the credible interval
#'   width used (e.g. \code{0.95}).
#' @slot fit The raw fitted model object returned by cmdStan (e.g. a
#'   \code{CmdStanMCMC} object), of which this object is a summary.
#'   
#' @name cvCTresults-class
#' @rdname cvCTresults-class
#' @exportClass cvCTresults
setClass("cvCTresults",
  representation(
   folds                   = "list", 
   theta                   = "CTparameter",
   beta                    = "CTparameter",
   eta                     = "CTparameter",
   topic_correlations      = "matrix",
   top_genes               = "list",
   regression_coefficients = "data.frame",
   cv_predictions          = "cvCTprediction",
   response_levels         = "character",
   gene_names              = "character",
   credible_interval       = "numeric",
   fit                     = "ANY"
  )
)

setValidity("cvCTresults", function(object) {
  errs <- character()

  # No additional cross-slot invariants beyond what CTresults and
  # cvCTprediction already enforce individually.

  if (length(errs) == 0) TRUE else errs
})

#' Constructor for cvCTresults
#'
#' @param folds List of numeric vectors, test folds created
#' @param theta An \code{CTparameter} object.
#' @param beta An \code{CTparameter} object.
#' @param eta An \code{CTparameter} object.
#' @param topic_correlations Numeric matrix of topic correlations.
#' @param top_genes A named list of \code{data.frame}s, one per topic.
#' @param regression_coefficients A \code{data.frame} summary of \code{eta}
#'   melted by (category, topic) pairs.
#' @param cv_predictions A \code{cvCTprediction} object.
#' @param response_levels Character vector of category labels.
#' @param gene_names Character vector of gene names.
#' @param credible_interval Numeric scalar, credible interval width.
#' @param fit The final model (an \code{CTresults}-producing cmdStan fit,
#'   or the raw cmdStan fit object) trained on the full dataset.
#'
#' @return An object of class \code{cvCTresults}.
#' @export
cvCTresults <- function(folds,
                         theta,
                         beta,
                         eta,
                         topic_correlations,
                         top_genes,
                         regression_coefficients,
                         cv_predictions,
                         response_levels,
                         gene_names,
                         credible_interval,
                         fit) {
  new("cvCTresults",
      folds                   = folds,
      theta                   = theta,
      beta                    = beta,
      eta                     = eta,
      topic_correlations      = topic_correlations,
      top_genes               = top_genes,
      regression_coefficients = regression_coefficients,
      cv_predictions          = cv_predictions,
      response_levels         = response_levels,
      gene_names              = gene_names,
      credible_interval       = credible_interval,
      fit                     = fit)
}

#' @param x A cvCTresults object.
#' @param object A cvCTresults object.
#' @rdname cvCTresults-class
#' @export
setGeneric("getCVPredictions", function(x) standardGeneric("getCVPredictions"))

#' @rdname cvCTresults-class
#' @export
setMethod("getCVPredictions", "cvCTresults", function(x) x@cv_predictions)

#' @rdname cvCTresults-class
#' @export
setMethod("show", "cvCTresults", function(object) {
  cat("cvCTresults object (ClassTopics, cross-validated)\n")
  cat("  Folds:            ", length(object@folds), "\n")
  cat("  Topics:           ", ncol(object@theta@mean), "\n")
  cat("  Documents:        ", nrow(object@theta@mean), "\n")
  cat("  Genes:            ", length(object@gene_names), "\n")
  cat("  Categories:       ", length(object@response_levels), "\n")
  cat("  Credible interval:", object@credible_interval, "\n")
  cat("  CV test accuracy: ", round(object@cv_predictions@cv_test_accuracy, 3), "\n")
})
