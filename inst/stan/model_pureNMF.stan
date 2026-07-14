// =============================================================================
// Supervised NMF / sLDA Model
// =============================================================================
//
// Generative model:
//   H[d, k] ~ Gamma(shape, rate)           // topic loadings per observation
//   W[k, v] ~ Gamma(shape, rate)           // variable weights per topic
//   eta[c, k]   ~ Normal(0, sigma_eta)         // class-topic regression weights
//
//   lambda[d, v]  = dot_product(H[d,:], W[:,v])   // Poisson rate
//   counts[d, v] ~ Poisson(lambda[d, v])            // NMF likelihood
//
//   linear_pred[d, c] = eta[c,:] * theta[d,:]'
//   y[d] ~ Categorical(softmax(linear_pred[d]))          // supervised likelihood
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
  real<lower=0> shape;                         // Gamma prior shape for H, W
  real<lower=0> rate;                          // Gamma prior rate  for H, W
  real<lower=0> sigma_eta;                     // Normal prior SD for eta
  
  real<lower=0> lambda_cat;                    // likelihood weight for
                                               // categorical component
}

parameters{
  matrix<lower=0>[D, K] H;                 // topic loadings:      D x K
  matrix<lower=0>[K, V] W;                 // variable-topic weights:  K x V
  matrix[C, K]          eta_raw;               // class-topic weights: C x K
}

transformed parameters{
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

  // ------------------------------------------------------------------
  // theta : D x K  topic proportions per patient
  //   Obtained according to Carbonetto et al. (2021).
  //
  //   Procedure:
  //     HU[d, k]    = H[d, k] * u[k]              // scale by topic weight
  //     theta[d, k] = HU[d, k] / sum_k HU[d, k']  // row-normalize
  //
  //   Rows of theta sum to 1 and are interpretable as the fraction
  //   of each patient's variable expression explained by eachtopic.
  // ------------------------------------------------------------------
  matrix[D, K] theta;
  for(d in 1:D){
    vector[K] HU_d = H[d, :]' .* u;             // elementwise: K-vector
    theta[d, :] = HU_d' / sum(HU_d);
  }

  // ------------------------------------------------------------------
  //   Obtained according to the Poisson Non-negative Matrix Factorization to 
  //   Multinomial Topic Model reparameterization
  //   (see Carbonetto et al. 2021).
  //   Each row sums to 1: beta[k, v] = W[k, v] / u[k]
  // ------------------------------------------------------------------
  matrix[K, V] beta;
  for(k in 1:K){
    beta[k, :] = W[k, :] / u[k];
  }
  
  matrix[C, K] eta = sigma_eta * eta_raw;   // scaled matrix used in likelihood
}

model{
  // ------------------------------------------------------------------
  // Priors
  // ------------------------------------------------------------------
  for(d in 1:D){
    H[d, :] ~ gamma(shape, rate);
  }
  
  for(k in 1:K){
    W[k, :]  ~ gamma(shape, rate);
  }
  
  for(c in 1:C){
    // Prior on eta_raw is standard normal — well-conditioned geometry
    eta_raw[c, :] ~ normal(0, 1);
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

  // ------------------------------------------------------------------
  // Supervised likelihood: categorical with softmax linear predictor
  // linear_pred[c] = eta[c, :] * theta[d, :]'
  // ------------------------------------------------------------------
  for(d in 1:D){
    vector[C] linear_pred;
    for (c in 1:C){
      linear_pred[c] = dot_product(eta[c, :], theta[d, :]);
    }
    target += lambda_cat * categorical_logit_lpmf(y[d] | linear_pred);
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

  // ------------------------------------------------------------------
  // Topic correlations (diagnostic: are topics distinguishable?)
  // Uses beta_norm so scale differences don't dominate correlation.
  // ------------------------------------------------------------------
  matrix[K, K] topic_correlations;
  for(k1 in 1:K){
    for(k2 in 1:K){
      if(k1 == k2){
        topic_correlations[k1, k2] = 1.0;
      } else {
        real mean1 = mean(beta[k1, :]);
        real mean2 = mean(beta[k2, :]);
        vector[V] dev1 = beta[k1, :]' - mean1;
        vector[V] dev2 = beta[k2, :]' - mean2;
        real cov12     = dot_product(dev1, dev2);
        real var1      = dot_self(dev1);
        real var2      = dot_self(dev2);
        topic_correlations[k1, k2] = cov12 / (sqrt(var1) * sqrt(var2));
      }
    }
  }
}

