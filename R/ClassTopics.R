# =============================================================================
# Binary or Categorical Topic Modeling -- EM warm-up for Stan initialization
# =============================================================================

#' @param H D x K  initial topic loadings (strictly positive)
#' @param W K x V  initial topic distribution estimates (strictly positive)
#' @param X D x V  count matrix
#' @param n_iter_EM number of outer EM iterations
#' @param H_only logical that controls if EM is executed for training
#'                  (default) or \code{123})
#' @param verbose logical that controls whether the function prints information
#'                regarding the algorithm steps (default = \code{TRUE})

.em_nmf <- function(H, W, X,
                   n_iter_EM = 20,
                   verbose = TRUE,
                   H_only = FALSE){
  
  D <- nrow(X)
  V <- ncol(X)
  K <- ncol(H)
  
  stopifnot(
    nrow(H) == D, ncol(H) == K,
    nrow(W) == K, ncol(W) == V
  )
  
  FLOOR <- 1e-10
  
  for (iter in seq_len(n_iter_EM)) {
    
    # E-step
    lambda <- H %*% W
    ratio_obs_est <- X / lambda
    
    # M-step: theta -- uses H and correction
    nmf_num <- ratio_obs_est %*% t(W)
    nmf_denom <- matrix(rowSums(W), nrow = D, ncol = K, byrow = TRUE)
    
    H_new <- H * (nmf_num / nmf_denom)
    
    if (!all(is.finite(H_new))) {
      if (verbose) cat(sprintf("Iter %3d | H non-finite -- stopping early\n", iter))
      break
    }
    H <- pmax(H_new, FLOOR)
    
    if(!H_only){
      # M-step: W -- no correction needed
      W_num <- t(H) %*% ratio_obs_est
      W_denom <- matrix(colSums(H), nrow = K, ncol = V, byrow = FALSE)
      W_new <- W * (W_num / W_denom)
      
      if (!all(is.finite(W_new))) {
        if (verbose) cat(sprintf("Iter %3d | W non-finite -- stopping early\n", iter))
        break
      }
      W <- pmax(W_new, FLOOR)
      
      # Rescale to prevent scale drift
      col_scale <- sqrt(colMeans(H) / colMeans(t(W)))
      H <- pmax(sweep(H, 2, col_scale, "/"), FLOOR)
      W <- pmax(sweep(W,  1, col_scale, "*"), FLOOR)
    }
    
    # Diagnostics
    if(verbose){
      lambda_new <- H %*% W
      ll_nmf <- sum(X * log(lambda_new) - lambda_new, na.rm = TRUE)
      
      cat(sprintf("Iter %3d | NMF ll: %12.2f\n", iter, ll_nmf))
    }
  }
  
  return(list(H = H))
}

# =============================================================================
# Binary or Categorical NMF-LDA -- Stan model fitting
# =============================================================================

#' Fits a supervised topic model using integrated LDA and NMF to data with
#' categorical response variable
#' 
#' @param counts count matrix
#' @param response character or factor vector on the class labels
#' @param K number of topics to be modeled (default = \code{6})
#' @param shape scalar prior shape of the NMF coefficients (default = 4)
#' @param sigma_eta scalar prior sd of the coefficients of the eta matrix
#'                  (default = 0.5)
#' @param betadir logical, controls whether beta is estimated with Dirichlet
#'                prior (default = \code{TRUE})
#' @param alpha_beta numeric hyperparameter of the symmetric Dirichlet
#'                   distribution if \code{betadir = TRUE} (default = 0.05)
#' @param lambda_ridge numeric ridge-penalty weight applied to the
#'                     regression coefficients \code{eta} (default = \code{0},
#'                     i.e. no ridge penalty)
#' @param n_iter_EM number of iterations of the EM pre-processing (default = \code{10})
#' @param verbose_EM logical that controls whether EM likelihood is printed
#'                   (default = \code{TRUE})
#' @param iter_warmup number of warm-up iterations for each chain
#'                    (default = \code{1000})
#' @param iter_sampling number of sampling iterations for each chain
#'                      (default = \code{1000})
#' @param cores number of CPU cores to use (default = \code{3})
#' @param chains number of chains to run (default = \code{3})
#' @param seed numeric seed for reproducibility (default = \code{123})
#' @param control control parameters for the Stan algorithm
#' @param verbose logical that controls whether the function prints information
#'                regarding the algorithm steps (default = \code{TRUE})
#' 
#' @return cmdstan model to interpret with the function ClassTopics_results
#' 
#' @details
#'   `counts` ~ Poisson(`lambda`), `lambda = H %*% W`    (NMF component)
#'   `response` ~ Categorical(softmax(`theta[i,] %*% eta`))  (supervised component)
#'
#' Parameters (matching Stan file names exactly):
#'   \code{H} : D x K  -- observation topic loadings
#'   \code{W} : K x V  -- topic-variable profiles
#'
#' Derived quantities (computed as transformed parameters in Stan,
#' NOT passed as initial values):
#'   \code{theta} : D x K  -- topic proportions per observation
#'   \code{beta}  : K x V  -- variable distributions per topic
#'
#' @export
ClassTopics <- function(counts,
                        response,
                        K = 6,
                        shape = 4,
                        sigma_eta = 0.5,
                        betadir = TRUE,
                        alpha_beta = 0.05,
                        lambda_ridge = 0,
                        n_iter_EM = 10,
                        verbose_EM = TRUE,
                        iter_warmup = 1000,
                        iter_sampling = 1000,
                        chains = 3,
                        cores = 3,
                        seed = 123,
                        control = list(adapt_delta = 0.90,
                                       max_treedepth = 15),
                        verbose = TRUE){
  
  # Input validation
  if(!is.matrix(counts) && !is.data.frame(counts)){
    stop("'counts' must be a matrix or data frame")
  }
  
  if(is.data.frame(counts)){
    counts <- as.matrix(counts)
  }
  
  if(!is.numeric(counts)){
    stop("'counts' must contain numeric values only")
  }
  
  if(any(counts < 0)){
    stop("'counts' must contain non-negative values only")
  }
  
  if(any(counts != floor(counts))){
    warning("Non-integer counts detected. Rounding to nearest integer.")
    counts <- round(counts)
  }
  
  counts <- matrix(as.integer(counts), nrow = nrow(counts))
  
  # Filter zero-count variables
  totals <- colSums(counts)
  nonzero_vars <- totals > 0
  n_zero_vars <- sum(!nonzero_vars)
  
  if(n_zero_vars > 0){
    stop("Remove columns with zero counts before fitting the model")
  }
  
  if(nrow(counts) != length(response)){
    stop("Number of patients must match length of response")
  }
  
  # Get dimensions
  D <- nrow(counts)
  V <- ncol(counts)
  
  # Process response
  if(!is.factor(response)){
    cat("\n'response' converted to factor")
    response <- factor(response)
  }
  response_levels <- levels(response)
  C <- length(response_levels)
  response_int <- as.integer(response)
  
  cat(sprintf("\n=== Dataset Summary ===\n"))
  cat(sprintf("Observations: %d\n", D))
  cat(sprintf("Variables: %d\n", V))
  cat(sprintf("Response categories: %d (%s)\n", C, 
              paste(response_levels, collapse = ", ")))
  
  # Calculate totals
  total_counts <- sum(counts)
  patient_totals <- rowSums(counts)
  
  cat(sprintf("Total counts: %.0f\n", total_counts))
  cat(sprintf("Average counts per observation: %.1f\n", total_counts / D))
  cat(sprintf("Counts range: [%d, %d]\n", min(counts),
              max(counts)))
  cat(sprintf("Count totals range: [%d, %d]\n", min(totals),
              max(totals)))
  
  if(any(totals == 0)){
    zero_observations <- sum(totals == 0)
    stop(sprintf("%d observation%s zero total counts. Remove before fitting", 
                 zero_observations, ifelse(zero_observations == 1, " has", "s have")))
  }
  
  cat(sprintf("\n=== Preparing ClassTopics Model with K=%d Topics ===\n", K))
  
  #### INITIALIZATION WITH EM ####
  
  mu_target <- sqrt(mean(counts) / K)
  
  rate <- shape / mu_target
  
  stan_init <- lapply(1:chains, function(chain){
    set.seed(chain)
    
    # Fresh random initialisation for each chain
    H0 <- matrix(rgamma(D * K, shape = shape, rate = rate), D, K)
    W0 <- matrix(rgamma(K * V, shape = shape, rate = rate), K, V)
    
    init <- .em_nmf(
      H = H0,
      W = W0,
      X = counts,
      n_iter_EM = n_iter_EM,
      verbose = verbose_EM,
      H_only = FALSE)
    
    return(init)
  })
  
  stan_data <- list(
    K = K,
    V = V,
    D = D,
    counts = counts,
    C = C,
    y = response_int,
    shape = shape,
    rate = rate,
    sigma_eta = sigma_eta,
    lambda_ridge_eta = lambda_ridge
  )
  
  if(betadir){
    stan_data$alpha_beta <- alpha_beta
    model_name <- "model_betadir"
  }
  else{
    model_name <- "model_pureNMF"
  }
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  if(verbose){
    cat("\n=== Stan Configuration ===\n")
    
    cat(sprintf("\nHyperparameters:\n"))
    cat(sprintf("  sigma_eta: %.3f\n", sigma_eta))
    cat(sprintf("  shape: %.3f\n", shape))
    cat(sprintf("  rate: %.3f\n", rate))
    if(betadir){
      cat(sprintf("  alpha_beta: %.3f\n", alpha_beta))
    }
  }
  
  mod <- .get_stan_model(model_name)

  # Fit model
  cat(sprintf("\n=== Fitting Model with K=%d ===\n", K))
  cat("This may take a while...\n")
  
  fit <- mod$sample(
    data = stan_data,
    chains = chains,
    parallel_chains = cores,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    init = stan_init,
    seed = seed
  )
  
  cat("\n=== Model fitted! ===\n")
  
  # Quick diagnostic
  if(verbose){
    mean_corr <- posterior::E(posterior::as_draws_rvars(
      fit$draws("topic_correlations"))$topic_correlations)
    cat("\nTopic correlations (off-diagonal should be low for good separation):\n")
    print(round(mean_corr, 3))
  }
  
  return(fit)
}

# =============================================================================
# K-FOLD CROSS-VALIDATION FOR ClassTopics
# =============================================================================

# Secondary function that prepares a list to accommodate the results from
# each fold in cross-validation
# 
#' @param D number of observations in the full dataset
#' @param response_levels character, class labels
#' @param k_folds number of folds for cross-validation
#' 
#' @return list with elements \code{all_test_predictions}, \code{all_test_pred_probs},
#'   \code{fold_train_accuracies}, \code{fold_test_accuracies},
#'   \code{fold_train_class_accuracies} and \code{fold_test_class_accuracies}
#'   
.build_cv_list <- function(D, response_levels, k_folds){
  C <- length(response_levels)
  
  pred_n_acc <- list(
    # Storage for train and test predictions
    all_test_predictions = rep(NA_character_, D), # test prediction for each obs
    all_test_pred_probs = matrix(NA, nrow = D,
                                 ncol = ifelse(C == 2, 1, C)),
    
    fold_train_accuracies = numeric(k_folds),
    fold_test_accuracies = numeric(k_folds),
    
    # Storage for per-category accuracies across folds
    fold_train_class_accuracies = matrix(NA, nrow = k_folds, ncol = C),
    fold_test_class_accuracies = matrix(NA, nrow = k_folds, ncol = C))
  
  
  colnames(pred_n_acc$fold_test_class_accuracies) = response_levels
  colnames(pred_n_acc$fold_train_class_accuracies) = response_levels
  
  return(pred_n_acc)
}

# Create Stratified Folds
#
#' @param response Factor vector of class labels to stratify by.
#' @param k_folds Integer, number of folds to create.
#'
#' @return A list of length `k_folds`, each element an integer vector of
#'   row indices assigned to that fold's test set.
#' @keywords internal
.create_stratified_folds <- function(response, k_folds){
  
  folds <- vector("list", k_folds)
  
  for(level in levels(factor(response))){
    level_indices <- which(response == level)
    n_level <- length(level_indices)
    
    level_indices_shuffled <- sample(level_indices)
    fold_assignment <- rep(1:k_folds, length.out = n_level)
    
    for(k in 1:k_folds){
      folds[[k]] <- c(folds[[k]], level_indices_shuffled[fold_assignment == k])
    }
  }
  
  return(folds)
}

# ==============================================================================
# PREDICTION FUNCTION
# ==============================================================================

#' Predict Mixed Memberships and Response for New Observations Given Learned Topics via Stan
#'
#' @param counts_test Matrix of counts
#' @param W Learned topic-variable loadings \code{K, V}
#' @param eta Learned regression coefficients \code{C, K}
#' @param response Character vector of response categories
#' @param shape scalar prior shape of H (default = 1)
#' @param cores number of CPU cores to use (default = \code{3})
#' @param chains number of chains to run (default = \code{3})
#' @param control control parameters for the Stan algorithm
#' @param seed numeric seed for reproducibility (default = \code{123})
#' @param iter_warmup number of warm-up iterations for each chain
#'                    (default = \code{1000})
#' @param iter_sampling number of sampling iterations for each chain
#'                      (default = \code{1000})
#' @param n_iter_EM Number of required EM iterations
#' @param verbose logical that controls whether EM likelihood is printed
#'                (default = \code{TRUE})
#'
#' @return List with predicted classes and probabilities
#' @export
predict_ClassTopics_stan <- function(
    counts_test,
    W,
    eta,
    response,
    shape = 1,
    cores = 3,
    chains = 3,
    control = list(adapt_delta = 0.90,
                   max_treedepth = 15),
    seed = 123,
    iter_warmup = 1000,
    iter_sampling = 1000,
    n_iter_EM = 4,
    verbose = TRUE){
  
  D <- nrow(counts_test)
  V <- ncol(W)
  K <- nrow(W)
  C <- nrow(eta)
  response_levels <- levels(response)
  
  if(C == 1){
    C <- 2
    eta <- rbind(eta, -eta)
  }
  
  mu_target <- sqrt(mean(counts_test) / K)
  
  rate <- shape / mu_target
  
  stan_init <- lapply(1:chains, function(chain){
    set.seed(chain)
    
    # Fresh random initialisation for each chain
    H0 <- matrix(rgamma(D * K, shape = shape, rate = rate), D, K)
    
    return(.em_nmf(
      H = H0,
      W = W,
      X = counts_test,
      H_only = TRUE,
      n_iter_EM = n_iter_EM,
      verbose = verbose
      ))
  })
  
  stan_data <- list(
    K = K,
    V = V,
    D = D,
    counts = counts_test,
    C = C,
    y = as.integer(response),
    shape = shape,
    rate = rate,
    W = W,
    eta = eta
  )
  
  cat("\n\nFitting theta to the test set using Stan...\n\n")
  
  mod <- .get_stan_model("model_test")
  
  fit <- mod$sample(
    data = stan_data,
    chains = chains,
    parallel_chains = cores,   # runs chains truly in parallel
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    init = stan_init,
    seed = seed
  )
  
  draws <- posterior::as_draws_rvars(fit$draws())
  
  predicted_class <- posterior::modal_category(draws$y_pred)
  predicted_probs <- posterior::E(draws$response_probs)
  
  if(C == 2){
    predicted_probs <- matrix(predicted_probs[, 1], ncol = 1)
  }
  
  return(CTprediction(
    predicted_class = response_levels[predicted_class],
    predicted_probs = predicted_probs
  ))
}

#' Predict Mixed Memberships and Response for New Patients Given Learned Topics via EM only
#'
#' @param counts_test Matrix of counts
#' @param W Learned topic-variable loadings \code{K, V}
#' @param eta Learned regression coefficients \code{C, K}
#' @param response Character vector of response categories
#' @param shape scalar prior shape of H (default = 1)
#' @param n_iter_EM Number of required EM iterations
#' @param verbose logical that controls whether EM likelihood is printed
#'                (default = \code{TRUE})
#'
#' @return List with predicted classes and probabilities
#' @export
predict_ClassTopics_EM <- function(
    counts_test,
    W,
    eta,
    response,
    shape = 1,
    n_iter_EM = 4,
    verbose = TRUE
  ){
  
  D <- nrow(counts_test)
  V <- ncol(W)
  K <- nrow(W)
  C <- nrow(eta)
  response_levels <- levels(response)
  
  if(C == 1){
    C <- 2
    eta <- rbind(eta, -eta)
  }
  
  mu_target <- sqrt(mean(counts_test) / K)
  
  rate <- shape / mu_target
  
  H0 <- matrix(rgamma(D * K, shape = shape, rate = rate), D, K)
  
  H <- .em_nmf(
    H = H0,
    W = W,
    X = counts_test,
    n_iter_EM = n_iter_EM,
    H_only = TRUE,
    verbose = verbose
  )$H
  
  u <- rowSums(W)
  
  Hu <- lapply(1:D, function(d) H[d, ] * u)
  
  sum_Hu <- sapply(Hu, sum)
  
  theta <- t(sapply(1:D, function(d) Hu[[d]] / sum_Hu[d]))
  
  theta_eta <- lapply(1:D, function(d) as.numeric(exp(eta %*% theta[d, ])))
  
  sum_theta_eta <- sapply(theta_eta, sum)
  
  predicted_probs <- t(sapply(1:D, function(d) theta_eta[[d]]/sum_theta_eta[d]))
  predicted_class <- sapply(1:D, function(d) which.max(predicted_probs[d,]))
  
  if(C == 2){
    predicted_probs <- matrix(predicted_probs[, 1], ncol = 1)
  }
  
  return(CTprediction(
    predicted_class = response_levels[predicted_class],
    predicted_probs = predicted_probs
  ))
}

# ============================================================================
# K-FOLD CROSS-VALIDATION FOR ClassTopics
# ============================================================================

# K-Fold Cross-Validation for ClassTopics (For unbiased estimation)
#
#' @param counts count matrix
#' @param response character or factor vector on the class labels
#' @param folds list whose elements are train-test partitions
#' @param fold index of the fold to be fitted a ClassTopics model
#' @param K_topics number of topics to be modeled (default = \code{3})
#' @param shape scalar prior shape of the NMF coefficients (default = 4)
#' @param shape_test scalar prior shape of the H matrix on test sets (default = 4)
#' @param sigma_eta scalar prior sd of the coefficients of the eta matrix
#'                  (default = 0.5)
#' @param betadir logical, controls whether beta is estimated with Dirichlet
#'                prior (default = \code{TRUE}) or, if \code{FALSE}, with Gamma prior
#' @param alpha_beta numeric hyperparameter of the symmetric Dirichlet
#'                   distribution if \code{betadir} = TRUE] (default = \code{0.05})
#' @param seed numeric seed for reproducibility (default = \code{123})
#' @param verbose_per_fold Logical, print details for each fold (default: \code{TRUE})
#' @param iter_warmup number of warm-up iterations for each chain
#'                    (default = \code{1000})
#' @param iter_sampling number of sampling iterations for each chain
#'                      (default = \code{1000})
#' @param n_iter_EM number of iterations of the EM pre-processing (default = \code{10})
#' @param n_iter_EM_test number of EM iterations used when fitting theta on
#'                       this fold's held-out (test) data (default = \code{4})
#' @param test_stan logical, if \code{TRUE} (default) test-fold topic
#'                  proportions are fit via Stan (see
#'                  [predict_ClassTopics_stan()]); if \code{FALSE}, via EM
#'                  only (see [predict_ClassTopics_EM()])
#' @param lambda_ridge numeric ridge-penalty weight applied to the
#'                     regression coefficients \code{eta} (default = \code{0},
#'                     i.e. no ridge penalty)
#' @param cores number of CPU cores to use (default = \code{3})
#' @param chains number of chains to run (default = \code{3})
#' @param control control parameters for the Stan algorithm
#' @param ... Additional arguments passed to [ClassTopics()] and [predict_ClassTopics_stan()]
#'
#' @return List with accuracy estimates and predictions for one of the folds
#' 
#' @details
#' This performs k-fold cross-validation:
#' - Each patient appears in test set exactly ONCE
#' - Training sets are DISJOINT (no overlap)
#' - Accuracy estimate is UNBIASED
#' 
#' Returns:
#' - cv_accuracy: Unbiased estimate of model performance
#' - final_model: Model trained on ALL data (for interpretation)
#' @keywords internal

.cv_ClassTopics_fold_by_fold <- function(
    counts,
    response,
    folds,
    fold, 
    K_topics = 3,
    shape = 4,
    shape_test = 1,
    sigma_eta = 0.5,
    betadir = TRUE,
    alpha_beta = 0.05,
    seed = 123,
    verbose_per_fold = TRUE,
    iter_warmup = 1000,
    iter_sampling = 1000,
    n_iter_EM = 4,
    n_iter_EM_test = 4,
    test_stan = TRUE,
    lambda_ridge = 0,
    chains = 3,
    cores = 3,
    control = list(adapt_delta = 0.90,
                   max_treedepth = 15),
    ...){
  
  cat("Worker PID:", Sys.getpid(), "\n")
  
  set.seed(seed)
  
  response <- as.factor(response)
  
  D <- nrow(counts)
  response_levels <- levels(response)
  C <- length(response_levels)
  
  k_folds <- length(folds)
  
  cat(sprintf("\n=== %d-Fold Cross-Validation ===\n", k_folds))
    
    if (!verbose_per_fold) {
      cat(sprintf("Fold %d/%d... ", fold, k_folds))
    } else {
      cat(sprintf("\n========== FOLD %d/%d ==========\n", fold, k_folds))
    }
    
    test_indices <- folds[[fold]]
    train_indices <- setdiff(1:D, test_indices)
    
    if (verbose_per_fold) {
      cat(sprintf("Train: %d patients, Test: %d patients\n", 
                  length(train_indices), length(test_indices)))
    }
    
    nonzero_train_cols <- (1:ncol(counts[train_indices,]))[
      colSums(counts[train_indices,]) != 0]
    
    train_counts <- counts[train_indices,
                           nonzero_train_cols,
                           drop = FALSE]

      # Fit on training fold
      fit_fold <- ClassTopics(
        counts = train_counts,
        response = response[train_indices],
        K = K_topics,
        sigma_eta = sigma_eta,
        shape = shape,
        betadir = betadir,
        alpha_beta = alpha_beta,
        seed = seed,
        verbose = FALSE,
        iter_warmup = iter_warmup,
        iter_sampling = iter_sampling,
        n_iter_EM = n_iter_EM,
        lambda_ridge = lambda_ridge,
        chains = chains,
        cores = cores,
        control = control,
        ...
      )
    
    if(verbose_per_fold){
      cat(sprintf("Model %d/%d fitted!\n", fold, k_folds))
      cat("Extracting parameters and predictions...\n")
    }
    
    draws <- posterior::as_draws_rvars(fit_fold$draws())
    
    eta_means_mat <- matrix(colMeans(posterior::E(draws$eta)),
                            nrow = C, ncol = K_topics, byrow = TRUE)
    
    beta_fold <- posterior::E(draws$beta)
    W_pred <- posterior::E(draws$W)
    eta_fold <- posterior::E(draws$eta) - eta_means_mat
    
    y_pred_mode <- posterior::modal_category(draws$y_pred)
    response_probs_mode <- posterior::E(draws$response_probs)
    
    # Calculate prediction accuracy
    fold_tr_acc <- mean(response_levels[y_pred_mode] == response[train_indices])
    
    pred_n_acc <- .build_cv_list(
      D = D,
      response_levels = response_levels,
      k_folds = k_folds)
    
    pred_n_acc$fold_train_accuracies[fold] <- fold_tr_acc
    
    # ADDED: Calculate per-category training accuracy
    train_confusion <- table(
      Predicted = response_levels[y_pred_mode],
      Truth = response[train_indices]
    )
    train_class_acc <- diag(train_confusion) / table(response[train_indices])
    pred_n_acc$fold_train_class_accuracies[fold, names(train_class_acc)] <-
      train_class_acc
    
    if(verbose_per_fold){
      cat("Extracted for training set!\n")
      cat("Predicting for test set...\n")
    }
    
    # Predict on test fold
    
    mod <- .get_stan_model("model_test")
    
    if(test_stan){
      test_pred <- predict_ClassTopics_stan(
        counts_test = counts[test_indices,
                             nonzero_train_cols,
                             drop = FALSE],
        W = W_pred,
        eta = eta_fold,
        response = response[test_indices],
        shape = shape_test,
        seed = seed,
        iter_warmup = iter_warmup,
        iter_sampling = iter_sampling,
        n_iter_EM = n_iter_EM_test,
        chains = chains,
        cores = cores,
        control = control
      )
    }
    else{
      test_pred <- predict_ClassTopics_EM(
        counts_test = counts[test_indices,
                             nonzero_train_cols,
                             drop = FALSE],
        W = W_pred,
        eta = eta_fold,
        shape = shape_test,
        response = response[test_indices],
        n_iter_EM = n_iter_EM_test
      )
    }
    
    if(verbose_per_fold){
      cat("Predicted!\n")
    }
    
    # Store test predictions
    pred_n_acc$all_test_predictions[test_indices] <-
      test_pred@predicted_class
    
    pred_n_acc$all_test_pred_probs[test_indices, ] <-
      test_pred@predicted_probs
    
    # Calculate test accuracy
    test_response_factor <- factor(response[test_indices], levels = response_levels)
    pred_n_acc$fold_test_accuracies[fold] <- mean(test_pred@predicted_class ==
                                         test_response_factor)
    
    # Calculate per-category test accuracy
    test_confusion <- table(
      Predicted = test_pred@predicted_class,
      Truth = test_response_factor
    )
    test_class_acc <- diag(test_confusion) / table(test_response_factor)
    pred_n_acc$fold_test_class_accuracies[fold, names(test_class_acc)] <- test_class_acc
    
    if(verbose_per_fold){
      cat("Fold results extracted!\n\n")
      
      cat(sprintf("Train: %.1f%%, Test: %.1f%%\n", 
                  pred_n_acc$fold_train_accuracies[fold] * 100, 
                  pred_n_acc$fold_test_accuracies[fold] * 100))
    } else {
      cat(sprintf("Fold %d Train Accuracy: %.2f%%\n", fold,
                  pred_n_acc$fold_train_accuracies[fold] * 100))
      cat(sprintf("Fold %d Test Accuracy: %.2f%%\n", fold,
                  pred_n_acc$fold_test_accuracies[fold] * 100))
    }
    
  return(pred_n_acc)
}

# K-Fold Cross-Validation for ClassTopics (For unbiased estimation)
#
#' @param counts count matrix
#' @param response character or factor vector on the class labels
#' @param folds list whose elements are train-test partitions (as produced
#'              by [.create_stratified_folds])
#' @param pred_n_acc list containing the results of the previously obtained
#'                   classification results with .cv_ClassTopics_fold_by_fold
#' @param final_model argument that may receive the estimated full model
#'                    obtained with ClassTopics. If \code{NULL} (default), said model is
#'                    fitted internally
#' @param cores number of CPU cores (default = \code{3})
#' @param chains number of chains to run (default = \code{3})
#' @param seed numeric seed for reproducibility (default = \code{123})
#' @param control control parameters for the Stan algorithm
#' @param ... Additional arguments passed to ClassTopics and predict_ClassTopics
#'            (both _stan and _EM)
#'
#' @return List with summarized results of the topic model and cross-validated
#'         classification results
#' 
#' @details
#' This performs k-fold cross-validation:
#' - Each patient appears in test set exactly ONCE
#' - Training sets are DISJOINT (no overlap)
#' - Accuracy estimate is UNBIASED
#' 
#' Returns:
#' - \code{cv_accuracy}: Unbiased estimate of model performance
#' - \code{final_model}: Model trained on ALL data (for interpretation)
#' @keywords internal

.cv_ClassTopics_overall <- function(
    counts,
    response,
    folds,
    pred_n_acc,
    final_model = NULL,
    chains = 3,
    cores = 3,
    seed = 123,
    control = list(adapt_delta = 0.90,
                   max_treedepth = 15),
    ...){
  
  
  response <- as.factor(response)
  
  D <- nrow(counts)
  response_levels <- levels(response)
  C <- length(response_levels)
  # Calculate overall metrics
  all_test_predictions_factor <- factor(pred_n_acc$all_test_predictions,
                                        levels = response_levels)
  overall_test_accuracy <- mean(all_test_predictions_factor == response,
                                na.rm = TRUE)
  
  # ADDED: Calculate mean and SD for per-category accuracies
  mean_test_class_acc <- colMeans(pred_n_acc$fold_test_class_accuracies, na.rm = TRUE)
  sd_test_class_acc <- apply(pred_n_acc$fold_test_class_accuracies, 2, sd, na.rm = TRUE)
  
  mean_train_class_acc <- colMeans(pred_n_acc$fold_train_class_accuracies, na.rm = TRUE)
  sd_train_class_acc <- apply(pred_n_acc$fold_train_class_accuracies, 2, sd, na.rm = TRUE)
  
  cat(sprintf("\n=== Cross-Validation Results ===\n"))
  
  # Print training accuracy
  cat(sprintf("\nTraining Accuracy (per fold):\n"))
  with(pred_n_acc, cat(sprintf("  Mean: %.2f%% with SD of %.2f%%\n", 
                               mean(fold_train_accuracies) * 100, sd(fold_train_accuracies) * 100)))
  with(pred_n_acc, cat(sprintf("  Range: %.2f%% - %.2f%%\n", 
                               min(fold_train_accuracies) * 100, max(fold_train_accuracies) * 100)))
  
  cat(sprintf("\nTest (Validation) Accuracy (per fold):\n"))
  with(pred_n_acc, cat(sprintf("  Mean: %.2f%% with SD of %.2f%%\n", 
                               mean(fold_test_accuracies) * 100, sd(fold_test_accuracies) * 100)))
  with(pred_n_acc, cat(sprintf("  Range: %.2f%% - %.2f%%\n", 
                               min(fold_test_accuracies) * 100, max(fold_test_accuracies) * 100)))
  
  # Check for overfitting
  overfit_gap <- with(pred_n_acc,
                      mean(fold_train_accuracies) - mean(fold_test_accuracies))
  cat(sprintf("\nOverfitting Gap: %.2f%%\n", overfit_gap * 100))
  if(overfit_gap > 0.10){
    cat("  Large gap suggests overfitting - consider regularization\n")
  } else if(overfit_gap > 0.05){
    cat("  Moderate gap - model may be slightly overfitting\n")
  } else {
    cat("  Small gap - model generalizes well\n")
  }
  
  # Confusion matrix (test)
  confusion_cv <- table(
    Predicted = all_test_predictions_factor,
    Truth = response
  )
  
  cat("\nTest Confusion Matrix:\n")
  print(confusion_cv)
  
  # ADDED: Per-category accuracies with standard deviations
  cat("\n=== Per-Category Accuracies (Mean and SD across folds) ===\n")
  cat("\nTraining:\n")
  for(i in 1:C){
    cat(sprintf("  %s: %.1f%% with SD of %.1f%%\n", 
                response_levels[i],
                mean_train_class_acc[i] * 100,
                sd_train_class_acc[i] * 100))
  }
  
  cat("\nTest (Validation):\n")
  for(i in 1:C){
    cat(sprintf("  %s: %.1f%% and SD of %.1f%%\n", 
                response_levels[i],
                mean_test_class_acc[i] * 100,
                sd_test_class_acc[i] * 100))
  }
  
  if(is.null(final_model)){
    # Fit final model on ALL data for interpretation
    cat("\n=== Fitting Final Model on Full Dataset ===\n")
    cat("(This model is for interpretation; CV accuracy reported above)\n")
    
    final_model <- ClassTopics(
      counts = counts,
      response = response,
      chains = chains,
      cores = cores,
      control = control,
      seed = seed,
      ...
    )
  }
  
  with(pred_n_acc, return(cvCTprediction_plusfm(
    # Folds considered in CV
    folds = folds,
    
    # Training metrics
    cv_train_accuracy_mean = mean(fold_train_accuracies),
    cv_train_accuracy_sd = sd(fold_train_accuracies),
    fold_train_accuracies = fold_train_accuracies,
    
    # Test (validation) metrics
    cv_test_accuracy = overall_test_accuracy,
    cv_test_accuracy_mean = mean(fold_test_accuracies),
    cv_test_accuracy_sd = sd(fold_test_accuracies),
    fold_test_accuracies = fold_test_accuracies,
    cv_test_confusion_matrix = confusion_cv,
    
    # Per-category metrics (train)
    cv_train_class_accuracy_mean = mean_train_class_acc,
    cv_train_class_accuracy_sd = sd_train_class_acc,
    fold_train_class_accuracies = fold_train_class_accuracies,
    
    # Per-category metrics (test)
    cv_test_class_accuracy_mean = mean_test_class_acc,
    cv_test_class_accuracy_sd = sd_test_class_acc,
    fold_test_class_accuracies = fold_test_class_accuracies,
    
    # ADDED: Overfitting assessment
    overfitting_gap = overfit_gap,
    
    # Predictions
    all_test_predictions = all_test_predictions_factor,
    all_test_pred_probs = all_test_pred_probs,
    
    # Final model
    final_model = final_model  # Trained on ALL data for interpretation
  )))
}

# Aggregate Per-Fold Predictions and Accuracies
#
# Combines the per-fold prediction/accuracy lists produced by
# `.cv_ClassTopics_fold_by_fold()` into pooled cross-validation summary
# statistics.
#
#' @param pred_n_acc_list A list of per-fold results, each element as
#'   returned by `.cv_ClassTopics_fold_by_fold()`.
#'
#' @return A list of pooled/aggregated cross-validation statistics (mean and
#'   SD of train/test accuracy, per-category accuracy, pooled confusion
#'   matrix, overfitting gap, and pooled out-of-fold predictions).
#' @keywords internal
.aggregate_folds <- function(pred_n_acc_list){
  
  nfolds <- length(pred_n_acc_list)
  pred_n_acc <- pred_n_acc_list[[1]]
  
  
  for(i in 2:nfolds){
    fold <- pred_n_acc_list[[i]]
    
    pred_n_acc$all_test_predictions[!is.na(fold$all_test_predictions)] <-
      fold$all_test_predictions[!is.na(fold$all_test_predictions)]
    
    pred_n_acc$all_test_pred_probs[!is.na(fold$all_test_pred_probs[,1]),] <-
      #pred_n_acc$all_test_pred_probs +
      fold$all_test_pred_probs[!is.na(fold$all_test_pred_probs[,1]),]
    
    pred_n_acc$fold_test_accuracies <- 
      pred_n_acc$fold_test_accuracies +
      pred_n_acc_list[[i]]$fold_test_accuracies
    
    pred_n_acc$fold_train_accuracies <- 
      pred_n_acc$fold_train_accuracies +
      pred_n_acc_list[[i]]$fold_train_accuracies
    
    pred_n_acc$fold_train_class_accuracies[i, ] <-
      pred_n_acc_list[[i]]$fold_train_class_accuracies[i, ]
    
    pred_n_acc$fold_test_class_accuracies[i, ] <-
      pred_n_acc_list[[i]]$fold_test_class_accuracies[i, ]
  }
  
  return(pred_n_acc)
}

# K-Fold Cross-Validation for ClassTopics (For unbiased estimation)
#
#' @param counts count matrix
#' @param response character or factor vector on the class labels
#' @param k_folds number of folds to create for cross-validation (default = 5)
#' @param K_topics number of topics to be modeled (default = 3)
#' @param shape scalar prior shape of the NMF coefficients (default = 4)
#' @param shape_test scalar prior shape of the H matrix on test sets (default = 4)
#' @param sigma_eta scalar prior sd of the coefficients of the eta matrix
#'                  (default = 0.5)
#' @param betadir logical, controls whether beta is estimated with Dirichlet
#'                prior (default = `TRUE`) or, if `FALSE`, with Gamma prior
#' @param alpha_beta numeric hyperparameter of the symmetric Dirichlet
#'                   distribution if `betadir = TRUE` (default = 0.05)
#' @param seed numeric seed for reproducibility (default = 123)
#' @param verbose_per_fold Logical, print details for each fold (default: TRUE)
#' @param iter_warmup number of warm-up iterations for each chain
#'                    (default = 1000)
#' @param iter_sampling number of sampling iterations for each chain
#'                      (default = 1000)
#' @param n_iter_EM number of iterations of the EM pre-processing (default = 10)
#' @param n_iter_EM_test number of EM iterations used when fitting theta on
#'                       each fold's held-out (test) data (default = 4)
#' @param test_stan logical, if `TRUE` (default) test-fold topic proportions
#'                  are fit via Stan (see [predict_ClassTopics_stan()]);
#'                  if `FALSE`, via EM only (see [predict_ClassTopics_EM()])
#' @param lambda_ridge numeric ridge-penalty weight applied to the
#'                     regression coefficients `eta` (default = 0, i.e. no
#'                     ridge penalty)
#' @param ... Additional arguments passed to ClassTopics and predict_ClassTopics (both _stan and _EM)
#'
#' @return List with accuracy estimates and predictions along with the
#'         full topic model
#' 
#' @details
#' This performs k-fold cross-validation:
#' - Each patient appears in test set exactly ONCE
#' - Training sets are DISJOINT (no overlap)
#' - Accuracy estimate is UNBIASED
#' 
#' Returns:
#' - cv_accuracy: Unbiased estimate of model performance
#' - final_model: Model trained on ALL data (for interpretation)
#' @keywords internal

.sequential_cv_ClassTopics <- function(
    counts,
    response,
    k_folds = 5,
    K_topics = 3,
    shape = 4,
    shape_test = 1,
    sigma_eta = 0.5,
    betadir = TRUE,
    alpha_beta = 0.05,
    seed = 123,
    verbose_per_fold = TRUE,
    iter_warmup = 1000,
    iter_sampling = 1000,
    n_iter_EM = 4,
    n_iter_EM_test = 4,
    test_stan = TRUE,
    lambda_ridge = 0,
    ...){
  
  set.seed(seed)
  
  response <- as.factor(response)
  
  D <- nrow(counts)
  response_levels <- levels(response)
  C <- length(response_levels)
  
  cat(sprintf("\n=== %d-Fold Cross-Validation ===\n", k_folds))
  
  # Create k-fold splits
  folds <- .create_stratified_folds(response, k_folds)
  
  pred_n_acc <- .build_cv_list(
    D = D,
    response_levels = response_levels,
    k_folds = k_folds)
  
  # Fit model for each fold
  for(fold in 1:k_folds){
    
    if(!verbose_per_fold){
      cat(sprintf("Fold %d/%d... ", fold, k_folds))
    } else {
      cat(sprintf("\n========== FOLD %d/%d ==========\n", fold, k_folds))
    }
    
    test_indices <- folds[[fold]]
    train_indices <- setdiff(1:D, test_indices)
    
    if(verbose_per_fold){
      cat(sprintf("Train: %d patients, Test: %d patients\n", 
                  length(train_indices), length(test_indices)))
    }
    
    nonzero_train_cols <- (1:ncol(counts[train_indices,]))[
      colSums(counts[train_indices,]) != 0]
    
    train_counts <- counts[train_indices,
                           nonzero_train_cols,
                           drop = FALSE]
  
    # Fit on training fold
    fit_fold <- ClassTopics(
      counts = train_counts,
      response = response[train_indices],
      K = K_topics,
      sigma_eta = sigma_eta,
      betadir = betadir,
      alpha_beta = alpha_beta,
      seed = seed,
      verbose = FALSE,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      n_iter_EM = n_iter_EM,
      lambda_ridge = lambda_ridge,
      ...
    )
    
    if(verbose_per_fold){
      cat(sprintf("Model %d/%d fitted!\n", fold, k_folds))
      cat("Extracting parameters and predictions...\n")
    }
    
    draws <- posterior::as_draws_rvars(fit_fold$draws())
    
    eta_means_mat <- matrix(colMeans(posterior::E(draws$eta)),
                            nrow = C, ncol = K_topics, byrow = TRUE)
    
    beta_fold <- posterior::E(draws$beta)
    W_pred <- posterior::E(draws$W)
    eta_fold <- posterior::E(draws$eta) - eta_means_mat
    
    y_pred_mode <- posterior::modal_category(draws$y_pred)
    response_probs_mode <- posterior::E(draws$response_probs)
    
    # Calculate prediction accuracy
    fold_tr_acc <- mean(response_levels[y_pred_mode] == response[train_indices])
    pred_n_acc$fold_train_accuracies[fold] <- fold_tr_acc
    
    # ADDED: Calculate per-category training accuracy
    train_confusion <- table(
      Predicted = response_levels[y_pred_mode],
      Truth = response[train_indices]
    )
    train_class_acc <- diag(train_confusion) / table(response[train_indices])
    pred_n_acc$fold_train_class_accuracies[fold, names(train_class_acc)] <-
      train_class_acc
    
    if(verbose_per_fold){
      cat("Extracted for training set!\n")
      cat("Predicting for test set...\n")
    }
    
    # Predict on test fold
    
    mod <- .get_stan_model("model_test")
    
    if(test_stan){
      test_pred <- predict_ClassTopics_stan(
        counts_test = counts[test_indices,
                             nonzero_train_cols,
                             drop = FALSE],
        W = W_pred,
        eta = eta_fold,
        response = response[test_indices],
        shape = shape_test,
        seed = seed,
        iter_warmup = iter_warmup,
        iter_sampling = iter_sampling,
        n_iter_EM = n_iter_EM,
        ...
      )
    }
    else{
      test_pred <- predict_ClassTopics_EM(
        counts_test = counts[test_indices,
                             nonzero_train_cols,
                             drop = FALSE],
        W = W_pred,
        eta = eta_fold,
        shape = shape_test,
        response = response[test_indices],
        n_iter_EM = n_iter_EM_test
      )
    }
    
    if(verbose_per_fold){
      cat("Predicted!\n")
    }
    
    # Store test predictions
    pred_n_acc$all_test_predictions[test_indices] <-
      test_pred@predicted_class
    
    pred_n_acc$all_test_pred_probs[test_indices, ] <-
      test_pred@predicted_probs
    
    # Calculate test accuracy
    test_response_factor <- factor(response[test_indices],
                                   levels = response_levels)
    pred_n_acc$fold_test_accuracies[fold] <- mean(test_pred@predicted_class ==
                                                    test_response_factor)
    
    # Calculate per-category test accuracy
    test_confusion <- table(
      Predicted = test_pred@predicted_class,
      Truth = test_response_factor
    )
    test_class_acc <- diag(test_confusion) / table(test_response_factor)
    pred_n_acc$fold_test_class_accuracies[fold, names(test_class_acc)] <-
      test_class_acc
    
    
    if(verbose_per_fold){
      cat("Fold results extracted!\n\n")
      
      cat(sprintf("Train: %.1f%%, Test: %.1f%%\n", 
                  pred_n_acc$fold_train_accuracies[fold] * 100, 
                  pred_n_acc$fold_test_accuracies[fold] * 100))
    } else {
      cat(sprintf("Fold %d Train Accuracy: %.2f%%\n", fold,
                  pred_n_acc$fold_train_accuracies[fold] * 100))
      cat(sprintf("Fold %d Test Accuracy: %.2f%%\n", fold,
                  pred_n_acc$fold_test_accuracies[fold] * 100))
    }
  }
  
  # Calculate overall metrics
  all_test_predictions_factor <- factor(pred_n_acc$all_test_predictions,
                                        levels = response_levels)
  overall_test_accuracy <- mean(all_test_predictions_factor == response,
                                na.rm = TRUE)
  
  # ADDED: Calculate mean and SD for per-category accuracies
  mean_test_class_acc <- colMeans(pred_n_acc$fold_test_class_accuracies,
                                  na.rm = TRUE)
  sd_test_class_acc <- apply(pred_n_acc$fold_test_class_accuracies, 2, sd,
                             na.rm = TRUE)
  
  mean_train_class_acc <- colMeans(pred_n_acc$fold_train_class_accuracies,
                                   na.rm = TRUE)
  sd_train_class_acc <- apply(pred_n_acc$fold_train_class_accuracies, 2, sd,
                              na.rm = TRUE)
  
  cat(sprintf("\n=== Cross-Validation Results ===\n"))
  
  # Print training accuracy
  with(pred_n_acc, {
    cat(sprintf("\nTraining Accuracy (per fold):\n"))
    cat(sprintf("  Mean: %.2f%% with SD of %.2f%%\n", 
                mean(fold_train_accuracies) * 100,
                sd(fold_train_accuracies) * 100))
    cat(sprintf("  Range: %.2f%% - %.2f%%\n", 
                min(fold_train_accuracies) * 100,
                max(fold_train_accuracies) * 100))
    
    cat(sprintf("\nTest (Validation) Accuracy (per fold):\n"))
    cat(sprintf("  Mean: %.2f%% with SD of %.2f%%\n", 
                mean(fold_test_accuracies) * 100,
                sd(fold_test_accuracies) * 100))
    cat(sprintf("  Range: %.2f%% - %.2f%%\n", 
                min(fold_test_accuracies) * 100,
                max(fold_test_accuracies) * 100))
  
    # Check for overfitting
    overfit_gap <- mean(fold_train_accuracies) - mean(fold_test_accuracies)
    cat(sprintf("\nOverfitting Gap: %.2f%%\n", overfit_gap * 100))
    if(overfit_gap > 0.10){
      cat("  Large gap suggests overfitting - consider regularization\n")
    } else if(overfit_gap > 0.05){
      cat("  Moderate gap - model may be slightly overfitting\n")
    } else {
      cat("  Small gap - model generalizes well\n")
    }
  })
  
  # Confusion matrix (test)
  confusion_cv <- table(
    Predicted = all_test_predictions_factor,
    Truth = response
  )
  
  cat("\nTest Confusion Matrix:\n")
  print(confusion_cv)
  
  # ADDED: Per-category accuracies with standard deviations
  cat("\n=== Per-Category Accuracies (Mean and SD across folds) ===\n")
  cat("\nTraining:\n")
  for(i in 1:C){
    cat(sprintf("  %s: %.1f%% with SD of %.1f%%\n", 
                response_levels[i],
                mean_train_class_acc[i] * 100,
                sd_train_class_acc[i] * 100))
  }
  
  cat("\nTest (Validation):\n")
  for(i in 1:C){
    cat(sprintf("  %s: %.1f%% and SD of %.1f%%\n", 
                response_levels[i],
                mean_test_class_acc[i] * 100,
                sd_test_class_acc[i] * 100))
  }
  
  # Fit final model on ALL data for interpretation
  cat("\n=== Fitting Final Model on Full Dataset ===\n")
  cat("(This model is for interpretation; CV accuracy reported above)\n")
  
  final_model <- ClassTopics(
    counts = counts,
    response = response,
    K = K_topics,
    seed = seed,
    sigma_eta = sigma_eta,
    alpha_beta = alpha_beta,
    betadir = betadir,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    n_iter_EM = n_iter_EM,
    ...
  )
  
  return(with(pred_n_acc, cvCTprediction_plusfm(
    # Folds considered in CV
    folds = folds,
    
    # Training metrics
    cv_train_accuracy_mean = mean(fold_train_accuracies),
    cv_train_accuracy_sd = sd(fold_train_accuracies),
    fold_train_accuracies = fold_train_accuracies,
    
    # Test (validation) metrics
    cv_test_accuracy = overall_test_accuracy,
    cv_test_accuracy_mean = mean(fold_test_accuracies),
    cv_test_accuracy_sd = sd(fold_test_accuracies),
    fold_test_accuracies = fold_test_accuracies,
    cv_test_confusion_matrix = confusion_cv,
    
    # Per-category metrics (train)
    cv_train_class_accuracy_mean = mean_train_class_acc,
    cv_train_class_accuracy_sd = sd_train_class_acc,
    fold_train_class_accuracies = fold_train_class_accuracies,
    
    # Per-category metrics (test)
    cv_test_class_accuracy_mean = mean_test_class_acc,
    cv_test_class_accuracy_sd = sd_test_class_acc,
    fold_test_class_accuracies = fold_test_class_accuracies,
    
    # ADDED: Overfitting assessment
    overfitting_gap = overfit_gap,
    
    # Predictions
    all_test_predictions = all_test_predictions_factor,
    all_test_pred_probs = all_test_pred_probs,
    
    # Final model
    final_model = final_model  # Trained on ALL data for interpretation
  )))
}

#' K-Fold Cross-Validation for ClassTopics Models (Parallel)
#'
#' Runs stratified k-fold cross-validation for the supervised topic model
#' fit by [ClassTopics()], dispatching each fold's fit -- plus one final
#' fit on the full dataset -- as concurrent background jobs via the
#' `future` package.
#'
#' @param counts count matrix
#' @param response character or factor vector on the class labels
#' @param k_folds number of folds to create for cross-validation (default = 5)
#' @param K_topics number of topics to be modeled (default = 3)
#' @param shape scalar prior shape of the NMF coefficients (default = 4)
#' @param sigma_eta scalar prior sd of the coefficients of the eta matrix
#' @param shape_test scalar prior shape of the H matrix on test sets (default = 4)
#' @param betadir logical, controls whether beta is estimated with Dirichlet
#'                prior (default = `TRUE`) or, if `FALSE`, with Gamma prior
#' @param alpha_beta numeric hyperparameter of the symmetric Dirichlet
#'                   distribution if `betadir = TRUE` (default = 0.05)
#' @param seed numeric seed for reproducibility (default = 123)
#' @param verbose_per_fold Logical, print details for each fold (default: TRUE)
#' @param iter_warmup number of warm-up iterations for each chain
#'                    (default = 1000)
#' @param iter_sampling number of sampling iterations for each chain
#'                      (default = 1000)
#' @param n_iter_EM number of iterations of the EM pre-processing (default = 4)
#' @param n_iter_EM_test number of EM iterations used when fitting theta on
#'                       each fold's held-out (test) data (default = 4)
#' @param test_stan logical, if `TRUE` (default) test-fold topic proportions
#'                  are fit via Stan (see [predict_ClassTopics_stan()]);
#'                  if `FALSE`, via EM only (see [predict_ClassTopics_EM()])
#' @param lambda_ridge numeric ridge-penalty weight applied to the
#'                     regression coefficients `eta` (default = 0, i.e. no
#'                     ridge penalty)
#' @param cores number of CPU cores to use per fold and within the final model (default = \code{3})
#' @param chains number of chains to run per fold and within the final model (default = \code{3})
#' @param control control parameters for the Stan algorithm
#' @param ... Additional arguments passed to [ClassTopics()] and
#'   [predict_ClassTopics_stan()] / [predict_ClassTopics_EM()]
#'
#' @return A [cvCTprediction_plusfm-class] object with per-fold and pooled
#'   train/test accuracy, per-category accuracy, the pooled test confusion
#'   matrix, an overfitting diagnostic, pooled out-of-fold predictions, and
#'   the final model fit on the full dataset.
#'
#' @details
#' This performs k-fold cross-validation:
#' - Each patient appears in test set exactly ONCE
#' - Training sets are DISJOINT (no overlap)
#' - Accuracy estimate is UNBIASED
#'
#' Folds (and the final full-data fit) are run concurrently via
#' `future::multisession`; the previous `future::plan()` is restored on exit.
#'
#' @export
cv_ClassTopics <- function(
    counts,
    response,
    k_folds = 5,
    K_topics = 3,
    shape = 4,
    shape_test = 1,
    sigma_eta = 0.5,
    betadir = TRUE,
    alpha_beta = 0.05,
    seed = 123,
    verbose_per_fold = TRUE,
    iter_warmup = 1000,
    iter_sampling = 1000,
    n_iter_EM = 4,
    n_iter_EM_test = 4,
    test_stan = TRUE,
    lambda_ridge = 0,
    chains = 3,
    cores = 3,
    control = list(adapt_delta = 0.90,
                   max_treedepth = 15),
    ...){
  
  set.seed(seed)
  
  response <- as.factor(response)
  
  D <- nrow(counts)
  response_levels <- levels(response)
  C <- length(response_levels)
  
  cat(sprintf("\n=== %d-Fold Cross-Validation ===\n", k_folds))
  
  # Create k-fold splits
  folds <- .create_stratified_folds(response, k_folds)
  
  # ---------------------------------------------------------------------
  # PARALLEL EXECUTION SETUP
  # ---------------------------------------------------------------------
  # Each fold is an independent job handled by
  # .cv_ClassTopics_fold_by_fold(), and the final_model fit is another
  # independent job -- these k_folds + 1 jobs are dispatched concurrently
  # using the 'future' package. Whatever cmdstan
  # chains / parallel_chains settings come from '...' or the defaults
  # inside ClassTopics()/predict_ClassTopics_stan() are passed through
  # completely untouched -- only the OUTER dispatch (fold jobs + final
  # model job) is parallelized here. This intentionally allows CPU
  # oversubscription when (chains per job) x (concurrent jobs) exceeds
  # the number of available cores.
  if(!requireNamespace("future", quietly = TRUE)){
    stop("Package 'future' is required to run cv_ClassTopics in parallel. ",
         "Please install it with install.packages('future').")
  }
  
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  
  # multisession spins up k_folds + 1 background R processes (one per
  # fold, plus one for the final_model fit) regardless of how many
  # physical cores are actually available -- this is deliberate.
  future::plan(future::multisession, workers = k_folds + 1)
  
  # ---------------------------------------------------------------------
  # Live, fold-tagged progress logging
  # ---------------------------------------------------------------------
  # Every job (each fold + the final model) writes its console output to
  # its own temp log file AS IT RUNS (via sink() inside the worker). The
  # main process polls those files in a loop and streams new lines to
  # the console immediately, prefixed with a job tag, so you get live
  # interleaved progress instead of one big buffered dump per job at the
  # end.
  log_dir <- tempfile("cv_ClassTopics_logs_")
  dir.create(log_dir)
  on.exit(unlink(log_dir, recursive = TRUE, force = TRUE), add = TRUE)
  
  fold_log_files <- file.path(log_dir, sprintf("fold_%d.log", 1:k_folds))
  final_log_file <- file.path(log_dir, "final_model.log")
  
  # Track how many lines of each log file we've already printed
  lines_read <- setNames(rep(0L, k_folds + 1),
                         c(sprintf("Fold %d", 1:k_folds), "Final Model"))
  log_files_by_tag <- setNames(c(fold_log_files, final_log_file),
                               names(lines_read))
  
  .flush_new_log_lines <- function(){
    for(tag in names(log_files_by_tag)){
      lf <- log_files_by_tag[[tag]]
      if(!file.exists(lf)) next
      all_lines <- suppressWarnings(readLines(lf, warn = FALSE))
      n_have <- lines_read[[tag]]
      if(length(all_lines) > n_have){
        new_lines <- all_lines[(n_have + 1):length(all_lines)]
        new_lines <- new_lines[nzchar(trimws(new_lines))]
        if(length(new_lines) > 0){
          cat(paste0("[", tag, "] ", new_lines, "\n"), sep = "")
        }
        lines_read[[tag]] <<- length(all_lines)
      }
    }
  }
  
  # ---------------------------------------------------------------------
  # Launch the final_model fit asynchronously (parallel to the folds),
  # redirecting its console output to its own log file as it runs
  # ---------------------------------------------------------------------
  final_model_future <- future::future({
    log_con <- file(final_log_file, open = "wt")
    sink(log_con, type = "output")
    sink(log_con, type = "message")
    on.exit({
      sink(type = "message")
      sink(type = "output")
      close(log_con)
    }, add = TRUE)
    
    ClassTopics(
      counts = counts,
      response = response,
      K = K_topics,
      seed = seed,
      sigma_eta = sigma_eta,
      shape = shape,
      alpha_beta = alpha_beta,
      betadir = betadir,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      n_iter_EM = n_iter_EM,
      chains = chains,
      cores = cores,
      control = control,
      ...
    )
  }, seed = TRUE)
  
  # ---------------------------------------------------------------------
  # Launch .cv_ClassTopics_fold_by_fold() for each fold in parallel,
  # each redirecting its console output to its own log file as it runs
  # ---------------------------------------------------------------------
  fold_futures <- vector("list", k_folds)
  for(fold in 1:k_folds){
    local({
      fold_local <- fold
      log_file_local <- fold_log_files[fold_local]
      
      fold_futures[[fold_local]] <<- future::future({
        log_con <- file(log_file_local, open = "wt")
        sink(log_con, type = "output")
        sink(log_con, type = "message")
        on.exit({
          sink(type = "message")
          sink(type = "output")
          close(log_con)
        }, add = TRUE)
        
        .cv_ClassTopics_fold_by_fold(
          counts = counts,
          response = response,
          folds = folds,
          fold = fold_local,
          K_topics = K_topics,
          sigma_eta = sigma_eta,
          shape = shape,
          shape_test = shape_test,
          betadir = betadir,
          alpha_beta = alpha_beta,
          seed = seed,
          verbose_per_fold = verbose_per_fold,
          iter_warmup = iter_warmup,
          iter_sampling = iter_sampling,
          n_iter_EM = n_iter_EM,
          n_iter_EM_test = n_iter_EM_test,
          test_stan = test_stan,
          lambda_ridge = lambda_ridge,
          chains = chains,
          cores = cores,
          control = control,
          ...
        )
      }, seed = TRUE)
    })
  }
  
  # ---------------------------------------------------------------------
  # Poll all jobs (folds + final model) until they finish, streaming
  # new log lines to the console live as they appear
  # ---------------------------------------------------------------------
  all_futures <- c(fold_futures, list(final_model_future))
  
  repeat{
    .flush_new_log_lines()
    
    if(all(future::resolved(all_futures))) break
    
    Sys.sleep(0.5)
  }
  
  # One last flush in case any lines were written between the final
  # resolved() check and job completion
  .flush_new_log_lines()
  
  # ---------------------------------------------------------------------
  # Collect fold results (values() re-raises any worker-side errors)
  # ---------------------------------------------------------------------
  pred_n_acc_list <- lapply(fold_futures, future::value)
  
  # ---------------------------------------------------------------------
  # Combine the per-fold outputs using the existing aggregation helper
  # ---------------------------------------------------------------------
  pred_n_acc <- .aggregate_folds(pred_n_acc_list)
  
  # ---------------------------------------------------------------------
  # Resolve the final_model fit that was started in parallel at the top
  # (fitted concurrently with the fold jobs)
  # ---------------------------------------------------------------------
  cat("\n=== Final Model on Full Dataset fitted (in parallel with folds) ===\n")
  cat("(This model is for interpretation; CV accuracy reported above)\n")
  final_model <- future::value(final_model_future)
  
  # ---------------------------------------------------------------------
  # Summarize everything using the existing overall-results helper
  # ---------------------------------------------------------------------
  all_results <- .cv_ClassTopics_overall(
    counts = counts,
    response = response,
    folds = folds,
    pred_n_acc = pred_n_acc,
    final_model = final_model,
    ...
  )
  
  return(all_results)
}

# ============================================================================
# RESULTS EXTRACTION AND INTERPRETATION
# ============================================================================

#' Extract and Interpret Supervised Topic Modeling Results
#'
#' @param fit Stan fit object from \code{ClassTopics}
#' @param true_response Character or factor vector on the class labels
#' @param vars_names Character vector of variable names (optional)
#' @param top_vars Integer, number of top variables to extract per topic
#' @param credible_interval Numeric, credible interval width (default: \code{0.95})
#'
#' @return List containing all results and interpretations
#' @export
ClassTopics_results <- function(fit, true_response,
                                 vars_names = NULL, top_vars = 10, 
                                 credible_interval = 0.95){
  
  # Calculate posterior summaries
  alpha_level <- (1 - credible_interval) / 2
  
  draws <- posterior::as_draws_rvars(fit$draws())
  
  response_levels <- levels(as.factor(true_response))
  
  D <- nrow(draws$theta)
  K <- ncol(draws$theta)
  V <- ncol(draws$beta)
  C <- length(response_levels)
  
  if(is.null(vars_names)){
    vars_names <- paste0("var_", 1:V)
  }
  
  eta_means_mat <- matrix(colMeans(posterior::E(draws$eta)),
                          nrow = C, ncol = K, byrow = TRUE)
  
  theta_mean <- posterior::E(draws$theta)
  theta_lower <- matrix(as.numeric(
    quantile(draws$theta, probs = alpha_level)),
    nrow = D, ncol = K)
  theta_upper <- matrix(as.numeric(
    quantile(draws$theta, probs = 1 - alpha_level)),
    nrow = D, ncol = K)
  
  theta <- CTparameter(mean = theta_mean,
                       lower = theta_lower,
                       upper = theta_upper)
  
  beta_mean <- posterior::E(draws$beta)
  beta_lower <- matrix(as.numeric(
    quantile(draws$beta, probs = alpha_level)),
    nrow = K, ncol = V)
  beta_upper <- matrix(as.numeric(
    quantile(draws$beta, probs = 1 - alpha_level)),
    nrow = K, ncol = V)
  
  beta <- CTparameter(mean = beta_mean,
                      lower = beta_lower,
                      upper = beta_upper)
  
  # Regression coefficients
  eta_mean <- posterior::E(draws$eta) - eta_means_mat
  eta_lower <- matrix(as.numeric(
    quantile(draws$eta, probs = alpha_level)),
    nrow = C, ncol = K) - eta_means_mat
  eta_upper <- matrix(as.numeric(
    quantile(draws$eta, probs = 1 - alpha_level)),
    nrow = C, ncol = K) - eta_means_mat
  
  eta <- CTparameter(mean = eta_mean,
                     lower = eta_lower,
                     upper = eta_upper)
  
  y_pred_mode <- posterior::modal_category(draws$y_pred)
  
  response_probs_mean <- posterior::E(draws$response_probs)
  
  if(C == 2){
    eta_mean <- matrix(eta_mean[1,], nrow = 1)
    eta_lower <- matrix(eta_lower[1,], nrow = 1)
    eta_upper <- matrix(eta_upper[1,], nrow = 1)
    
    response_probs_mean <- matrix(response_probs_mean[, 1], ncol = 1)
  }
  
  # Extract topic correlations
  topic_cors <- posterior::E(draws$topic_correlations)
  
  # Extract top variables for each topic
  top_vars_list <- list()
  for(k in 1:K){
    top_indices <- order(beta_mean[k, ], decreasing = TRUE)[1:top_vars]
    top_vars_list[[k]] <- data.frame(
      var = vars_names[top_indices],
      probability_mean = beta_mean[k, top_indices],
      probability_lower = beta_lower[k, top_indices],
      probability_upper = beta_upper[k, top_indices],
      topic = k,
      stringsAsFactors = FALSE
    )
  }
  names(top_vars_list) <- paste0("Topic_", 1:K)
  
  # Format regression coefficients for interpretation
  if(C == 2){
    Category <- paste0(response_levels[1], " vs ", response_levels[2])
  }
  else{
    Category <- response_levels
  }
  
    regression_results <- expand.grid(
      Category = Category,
      Topic = paste0("T", 1:K),
      stringsAsFactors = FALSE
    )
  
  regression_results$Coefficient_Mean <- as.vector(eta_mean)
  regression_results$Coefficient_Lower <- as.vector(eta_lower)
  regression_results$Coefficient_Upper <- as.vector(eta_upper)
  regression_results$Significant <- sign(regression_results$Coefficient_Lower) ==
    sign(regression_results$Coefficient_Upper)
  
  # Calculate prediction accuracy
  accuracy <- mean(response_levels[y_pred_mode] == true_response)
  
  # Create confusion matrix
  confusion_matrix <- table(
    Predicted = response_levels[y_pred_mode],
    # Actual = response_levels[actual_responses]
    Truth = true_response
  )
  
  # Per-category accuracy
  category_accuracy <- diag(confusion_matrix) / colSums(confusion_matrix)
  
  # Predictions and accuracy
  predictions <- CTtrainpred(
    predicted_categories = response_levels[y_pred_mode],
    true_categories = true_response,
    prediction_probabilities = response_probs_mean,
    overall_accuracy = accuracy,
    category_accuracy = category_accuracy,
    confusion_matrix = confusion_matrix
  )
  
  # Create comprehensive results object
  results <- CTresults(
    # Core parameters
    theta = theta,
    beta = beta,
    eta = eta,
    
    # Topic correlations
    topic_correlations = topic_cors,
    
    # Interpretable summaries
    top_vars = top_vars_list,
    regression_coefficients = regression_results,
    
    # Predictions and accuracy
    predictions = predictions,
    
    # Metadata
    response_levels = response_levels,
    vars_names = vars_names,
    credible_interval = credible_interval,
    fit = fit
  )
  
  return(results)
}

#' Extract and Interpret ClassTopics Model Results after cross-validation
#'
#' @param fit Stan fit object from cv_ClassTopics or .cv_ClassTopics_overall
#' @param true_response Character or factor vector on the class labels
#' @param vars_names Character vector of variable names (optional)
#' @param top_vars Integer, number of top variables to extract per topic
#' @param credible_interval Numeric, credible interval width (default: \code{0.95})
#'
#' @return List containing all results and interpretations
#' @export

cv_ClassTopics_results <- function(fit, true_response,
                                    vars_names = NULL, top_vars = 10, 
                                    credible_interval = 0.95){
  
  final_model_results <- ClassTopics_results(
    fit = fit@final_model,
    true_response = true_response,
    vars_names = vars_names,
    top_vars = top_vars,
    credible_interval = credible_interval)
  
  cv_predictions <- cvCTprediction(
    cv_train_accuracy_mean = fit@cv_train_accuracy_mean,
    cv_train_accuracy_sd = fit@cv_train_accuracy_sd,
    fold_train_accuracies = fit@fold_train_accuracies,
    cv_test_accuracy = fit@cv_test_accuracy,
    cv_test_accuracy_mean = fit@cv_test_accuracy_mean,
    cv_test_accuracy_sd = fit@cv_test_accuracy_sd,
    fold_test_accuracies = fit@fold_test_accuracies,
    cv_test_confusion_matrix = fit@cv_test_confusion_matrix,
    cv_train_class_accuracy_mean = fit@cv_train_class_accuracy_mean,
    cv_train_class_accuracy_sd = fit@cv_train_class_accuracy_sd,
    fold_train_class_accuracies = fit@fold_train_class_accuracies,
    cv_test_class_accuracy_mean = fit@cv_test_class_accuracy_mean,
    cv_test_class_accuracy_sd = fit@cv_test_class_accuracy_sd,
    fold_test_class_accuracies = fit@fold_test_class_accuracies,
    overfitting_gap = fit@overfitting_gap,
    all_test_predictions = fit@all_test_predictions,
    all_test_pred_probs = fit@all_test_pred_probs
  )
  
  results <- cvCTresults(
    folds = fit@folds,
    theta = final_model_results@theta,
    beta = final_model_results@beta,
    eta = final_model_results@eta,
    topic_correlations = final_model_results@topic_correlations,
    top_vars = final_model_results@top_vars,
    regression_coefficients = final_model_results@regression_coefficients,
    cv_predictions = cv_predictions,
    response_levels = final_model_results@response_levels,
    vars_names = final_model_results@vars_names,
    credible_interval = final_model_results@credible_interval,
    fit = fit@final_model
  )
  
  return(results)
}

# ============================================================================
# MODEL DIAGNOSTICS
# ============================================================================

#' Plot Topic Correlations
#'
#' @param results Results from [ClassTopics_results()] or [cv_ClassTopics_results()]
#'
#' @return A \code{ggplot} object showing a heatmap of topic correlations
#' @importFrom ggplot2 aes
#' @importFrom rlang .data
#' @export
plot_topic_correlations <- function(results){
  
  thres <- 0.15
  
  cor_data <- reshape2::melt(results@topic_correlations)
  cor_data <- dplyr::mutate(
    cor_data,
    Var1 = paste0("Topic ", .data$Var1),
    Var2 = paste0("Topic ", .data$Var2),
    value_grid = ifelse(.data$Var1 == .data$Var2, 0, .data$value),
    text = ifelse(.data$Var1 != .data$Var2, sprintf("%.3f", .data$value), "")
  )
  
  p <- ggplot2::ggplot(cor_data, aes(x = .data$Var1, y = .data$Var2, fill = .data$value_grid)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = "#FF4500", mid = "white", high = "#3A5FCD",
                                  name = "Correlation") +
    ggplot2::geom_text(aes(label = .data$text),
                       color = ifelse(abs(cor_data$value_grid) < thres,
                                      "black", "white"),
                       size = 3) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Topic Correlations")
  
  return(p)
}

#' Plot Topic-Response Associations
#'
#' @param results Results from [ClassTopics_results()] or [cv_ClassTopics_results()]
#' @param significance_only Logical, show only significant associations
#'
#' @return A \code{ggplot} object showing a heatmap of topic-response coefficients
#' @importFrom ggplot2 aes
#' @importFrom rlang .data
#' @export
plot_topic_response_heatmap <- function(results, significance_only = FALSE){
  
  coeff_data <- results@regression_coefficients
  
  if(significance_only){
    coeff_data <- coeff_data[coeff_data$Significant, ]
  }
  
  thres <- max(abs(coeff_data$Coefficient_Mean)) / 3.5
  
  p <- ggplot2::ggplot(coeff_data, aes(x = .data$Topic, y = .data$Category, fill = .data$Coefficient_Mean)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = "#FF4500", mid = "white", high = "#3A5FCD",
                                  name = "Coefficient") +
    ggplot2::geom_text(aes(label = ifelse(.data$Significant, sprintf("%.2f*", .data$Coefficient_Mean),
                                          sprintf("%.2f", .data$Coefficient_Mean))),
                       color = ifelse(abs(coeff_data$Coefficient_Mean) < thres,
                                      "black", "white"),
                       size = 3) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Topic-Response Associations",
                  subtitle = "* indicates significant associations",
                  x = "Topics", y = "Response Categories") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  return(p)
}

#' Plot Top Variables per Topic
#'
#' @param results Results from [ClassTopics_results()] or [cv_ClassTopics_results()]
#' @param n_vars Integer, number of top vars to show per topic
#'
#' @return A \code{ggplot} object showing the top vars for each topic
#' @importFrom ggplot2 aes
#' @importFrom rlang .data
#' @export
plot_top_vars <- function(results, n_vars = 5){
  
  # Combine top variables data
  top_vars_df <- do.call(rbind, lapply(names(results@top_vars), function(topic){
    df <- results@top_vars[[topic]][1:n_vars, ]
    df$Topic <- topic
    return(df)
  }))
  
  p <- ggplot2::ggplot(top_vars_df,
                       aes(x = tidytext::reorder_within(.data$var, .data$probability_mean, .data$Topic),
                           y = .data$probability_mean)) +
    ggplot2::geom_col(fill = "steelblue", alpha = 0.7) +
    ggplot2::geom_errorbar(aes(ymin = .data$probability_lower, ymax = .data$probability_upper),
                           width = 0.2, alpha = 0.6) +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~ Topic, scales = "free_y", ncol = 2) +
    tidytext::scale_x_reordered() +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Top Variables per Topic",
                  x = "Variables", y = "Probability") +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold"))
  
  return(p)
}