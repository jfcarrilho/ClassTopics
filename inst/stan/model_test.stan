// =============================================================================
// Supervised NMF / sLDA Model for test sets
// =============================================================================
//
// Generative model:
//   H[d, k] ~ Gamma(shape, rate)           // topic loadings per patient
//   W and eta come from the respective training posterior means
//
//   lambda[d, v]  = dot_product(H[d,:], W[:,v]) // Poisson rate
//   counts[d, v] ~ Poisson(lambda[d, v])          // NMF likelihood
//
//   linear_pred[d, c] = eta[c,:] * theta[d,:]'
//   y[d] ~ Categorical(softmax(linear_pred[d]))        // supervised likelihood
//
// Identifiability:
//   - H and W are non-negative; their product defines the Poisson rate
//   - eta is fully free (all C rows estimated); normal(0, sigma_eta) prior
//     resolves the softmax non-identifiability, mirroring glmnet's approach
//   - theta and beta are derived as transformed parameters and are
//     the interpretable topic representations
// =============================================================================

data{
  int<lower=2> K;                              // number of topics
  int<lower=2> V;                              // number of variables
  int<lower=1> D;                              // number of patients
  array[D, V] int<lower=0> counts;             // count matrix: D x V

  int<lower=2> C;                              // number of response categories
  array[D] int<lower=1, upper=C> y;            // class label for each patient

  // Hyperparameters
  real<lower=0> shape;                         // Gamma prior shape for theta
  real<lower=0> rate;                          // Gamma prior rate for theta
  
  matrix<lower=0>[K, V] W;                     // variable-topic weights:  K x V
  matrix[C, K]          eta;                   // class-topic weights: C x K
}

transformed data{
  // ------------------------------------------------------------------
  // u : K  topic scales
  //   u[k] = sum_v W[k, v]
  //   The total variable weight attributed to topic k. Topics with larger u
  //   contribute more to the overall Poisson reconstruction, so raw
  //   loadings H[d, k] must be rescaled by u[k] before comparing
  //   across topics within a patient.
  // ------------------------------------------------------------------
  vector[K] u;
  for(k in 1:K){
    u[k] = sum(W[k, :]);
  }

}

parameters{
  matrix<lower=0>[D, K] H;                     // topic loadings:      D x K
}

transformed parameters{
  // ------------------------------------------------------------------
  // theta : D x K  topic proportions per patient
  //   Obtained according to the Poisson Non-negative Matrix Factorization to 
  //   latent Dirichlet allocation reparameterization
  //   (see Carbonetto et al. 2021).
  //
  //   Procedure:
  //     HU[d, k] = H[d, k] * u[k]     // scale by topic weight
  //     theta[d, k] = HU[d, k] / sum_k HU[d, k']  // row-normalise
  //
  //   Rows of theta sum to 1 and are interpretable as the fraction
  //   of each patient's variable expression explained by each topic.
  // ------------------------------------------------------------------
  matrix[D, K] theta;
  for(d in 1:D){
    vector[K] HU_d = H[d, :]' .* u;             // elementwise: K-vector
    theta[d, :] = HU_d' / sum(HU_d);
  }
}

model{
  // ------------------------------------------------------------------
  // Priors
  // ------------------------------------------------------------------
  for(d in 1:D){
    H[d, :] ~ gamma(shape, rate);
  }

  // ------------------------------------------------------------------
  // NMF likelihood: Poisson with rate = H * W
  // Zeros are skipped — they contribute 0 to the Poisson log-pmf
  // only when the -lambda term is accounted for separately, so we
  // use target += and handle the full expression explicitly.
  // ------------------------------------------------------------------
  for(d in 1:D){
    target += -dot_product(H[d, :], W * rep_vector(1.0, V));
    for(v in 1:V){
      if (counts[d, v] > 0){
        real lambda_dv = dot_product(H[d, :], W[:, v]);
        target += counts[d, v] * log(lambda_dv);
      }
    }
  }
}

generated quantities{
  // ------------------------------------------------------------------
  // Log-likelihoods (for model comparison, LOO-CV, etc.)
  // ------------------------------------------------------------------
  real var_log_lik     = 0;
  real response_log_lik = 0;
  real total_log_lik;

  // ------------------------------------------------------------------
  // Posterior predictive
  // ------------------------------------------------------------------
  array[D] int<lower=1, upper=C> y_pred;
  array[D] vector[C]             response_probs;

  // NMF log-likelihood (sparse: skip zero counts still accounting for -lambda)
  for(d in 1:D){
    var_log_lik += -dot_product(H[d, :], W * rep_vector(1.0, V));
    for(v in 1:V){
      if (counts[d, v] > 0){
        real lambda_dv = dot_product(H[d, :], W[:, v]);
        var_log_lik += counts[d, v] * log(lambda_dv);
      }
    }
  }

  // Categorical log-likelihood and predictions
  for(d in 1:D){
    vector[C] linear_pred;
    for(c in 1:C){
      linear_pred[c] = dot_product(eta[c, :], theta[d, :]);
    }
    response_probs[d]  = softmax(linear_pred);
    y_pred[d]          = categorical_logit_rng(linear_pred);
    response_log_lik  += categorical_logit_lpmf(y[d] | linear_pred);
  }

  total_log_lik = var_log_lik + response_log_lik;
}

