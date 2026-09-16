#' Class CTsynthetic
#'
#' An S4 class representing a count dataset simulation following a supervised
#' non-negative matrix factorization (NMF) approach.
#'
#' @slot X Integer matrix of data counts.
#' @slot H Numeric matrix of observation-topic loadings.
#' @slot W Numeric matrix of topic-variable loadings.
#' @slot Y Factor of observed labels
#'
#' @name CTsynthetic-class
#' @rdname CTsynthetic-class
#' @exportClass CTsynthetic
setClass("CTsynthetic",
         representation(
           X = "matrix",
           H = "matrix",
           W = "matrix",
           Y = "factor"
         )
)

setValidity("CTsynthetic", function(object) {
  errs <- character()
  
  dims <- list(X = dim(object@X),
               H = dim(object@H),
               W = dim(object@W),
               Y = length(object@Y))
  
  if (!identical(dims$X[1], dims$H[1]) ||
      !identical(dims$X[2], dims$W[2]) ||
      !identical(dims$H[2], dims$W[1]) ||
      !identical(dims$X[1], dims$Y)) {
    errs <- c(errs, "Mismatching dimensions")
  }
  
  if (length(errs) == 0 &&
      (any(object@X < 0, na.rm = TRUE) ||
       any(object@H < 0, na.rm = TRUE) ||
       any(object@W < 0, na.rm = TRUE))) {
    errs <-
      c(errs,
        "One of the matrices X, H or W is not non-negative for at least one element")
  }
  
  if (length(errs) == 0) TRUE else errs
})

#' Constructor for CTsynthetic
#'
#' @param X Integer matrix of data counts.
#' @param H Numeric matrix of observation-topic loadings.
#' @param W Numeric matrix of topic-variable loadings.
#' @param Y Factor of observed labels
#'
#' @return An object of class \code{CTsynthetic}.
#' @export
CTsynthetic <- function(X, H, W, Y){
  new("CTsynthetic", X = X, H = H, W = W, Y = Y)
}


# =============================================================================
# Count data generative model - Softmax function
# =============================================================================

#' @param vec numerical vector for softmax function application

.softmax <- function(vec){exp(vec) / sum(exp(vec))}

# =============================================================================
# Count data generative model - Test for scalars
# =============================================================================

#' @param x evaluated object

.is.integer.scalar <- function(x){
  is.numeric(x) && length(x) == 1 && x == round(x)
}

# =============================================================================
# Count data generative model - NMF-based algorithm
# =============================================================================

#' Generates a class-balanced synthetic count dataset following a supervised
#' Poisson non-negative matrix factorization approach
#' 
#' @param D_per_class integer, number of observations to draw per class
#'                    (default = 50)
#' @param V integer, number of synthetic variables (default = 100)
#' @param eta numeric or matrix with the topic-response coefficients.
#'            If `NULL` (default), it gets defined as
#'            `eta <- matrix(rnorm(15), 3, 5)`
#' @param alpha_w numeric, controls data sparsity (default = 0.05)
#' @param lambda numeric, corresponds to the magnitude of the dataset
#'               (default = 1e6)
#' @param class_spread numeric that controls level of class separation
#'                     (small `class_spread` favors class overlapping,
#'                     large `class_spread` favors class separation).
#'                     Default = 2
#' @param seed numeric seed for reproducibility (default = 123)
#' 
#' @return an object of class CTsynthetic
#' 
#' @details When 'eta' is given as a 'matrix' object, its rows correspond
#'          to classes and columns to topics. If 'eta' is a 'numeric' array,
#'          the function assumes the symmetric binary class construction and
#'          executes `eta <- rbind(eta, -eta)` before simulation so that
#'          softmax is applied smoothly.
#'          
#' @importFrom stats rnorm
#' @importFrom stats rpois
#' 
#' @examples
#' 
#' ## Generate with the default parameters
#' 
#' df <- generate_CTdata()
#' 
#' ## Specific multinomial case with 3 classes
#' 
#' eta <- matrix(c( 0.4,  0.4,  0.0,  2.0, -1.1, -0.7,
#'                 -0.1,  0.5, -0.3, -1.3,  2.2, -0.1,
#'                 -0.3, -0.9,  0.3, -0.7, -1.1,  0.8),
#'                 nrow = 3, byrow = TRUE)
#' 
#' D_per_class <- 200
#' V <- 105
#' eta <- eta
#' alpha_w <- 0.05
#' lambda <- 1e6
#' class_spread <- 2
#' seed <- 123
#' 
#' df <- generate_CTdata(D_per_class = D_per_class,
#'                       V = V,
#'                       eta = eta,
#'                       alpha_w = alpha_w,
#'                       lambda = lambda, 
#'                       class_spread = class_spread,
#'                       seed = seed)
#'                       
#' ## Symmetric binary case, reusing all parameters except for 'eta'
#' 
#' eta <- c(-1.1, 0.0, 0.5, -0.9, 0.8)
#' 
#' # Also an option that produces the same result:
#' # eta <- matrix(c(-1.1, 0.0, 0.5, -0.9, 0.8), nrow = 1)
#'                 
#' df <- generate_CTdata(D_per_class = D_per_class,
#'                       V = V,
#'                       eta = eta,
#'                       alpha_w = alpha_w,
#'                       lambda = lambda, 
#'                       class_spread = class_spread,
#'                       seed = seed)
#'                       
#' ## Asymmetric binary case, reusing all parameters except for 'eta'
#' 
#' eta <- matrix(c( 0.4,  0.4,  0.0,  2.0, -1.1, -0.7,
#'                 -0.3, -0.9,  0.3, -0.7, -1.1,  0.8),
#'                 nrow = 2, byrow = TRUE)
#'                 
#' df <- generate_CTdata(D_per_class = D_per_class,
#'                       V = V,
#'                       eta = eta,
#'                       alpha_w = alpha_w,
#'                       lambda = lambda, 
#'                       class_spread = class_spread,
#'                       seed = seed)
#' 
#' @export
generate_CTdata <- function(D_per_class = 50,
                            V = 100,
                            eta = NULL,
                            alpha_w = 0.05,
                            lambda = 1e6,
                            class_spread = 2,
                            seed = 123){
  
  if(is.null(eta)){
    set.seed(seed)
    eta <- matrix(rnorm(15), 3, 5)
  }
  
  
  stopifnot(
    .is.integer.scalar(D_per_class), D_per_class > 0,
    .is.integer.scalar(V), V > 0,
    ((class(eta)[1] == "matrix" && nrow(eta) >= 1 && ncol(eta) >= 1) ||
       class(eta)[1] != "matrix" && is.numeric(eta)),
    is.numeric(alpha_w), alpha_w > 0, length(alpha_w) == 1,
    .is.integer.scalar(lambda), lambda > 0,
    is.numeric(class_spread), class_spread > 0, length(class_spread) == 1
  )
  
  if(class(eta)[1] != "matrix"){
    eta <- rbind(eta, -eta)
  }
  else{
    if(nrow(eta) == 1){
      eta <- rbind(eta, -eta)
    }
  }
  
  C <- nrow(eta)
  K <- ncol(eta)
  D <- C * D_per_class
  Y <- rep(1:C, each = D_per_class)
  
  eta_theta <- t(sapply(1:C, function(r) .softmax(eta[r, ])))
  
  alpha_theta <- class_spread * eta_theta[Y,]
  
  set.seed(seed)
  W <- matrix(rgamma(K * V, shape = alpha_w,
                     rate = 1/sqrt(lambda)),
              K, V)
  H <- matrix(rgamma(D * K, shape = as.numeric(alpha_theta),
                     rate = class_spread/sqrt(lambda)),
              D, K)
  
  X <- matrix(rpois(D * V, as.numeric(H %*% W)),
              nrow = D, ncol = V)
  
  return(CTsynthetic(X = X,
                     H = H,
                     W = W,
                     Y = as.factor(Y)))
}
