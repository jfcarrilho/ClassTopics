// =============================================================================
// Supervised NMF / sLDA Model
// =============================================================================
//
// Generative model:
//   H[d, k] ~ Gamma(shape, shape * rate)  // topic loadings per observation
//   W[k, v] ~ Gamma(alpha_beta, rate)     // variable weights per topic
//   eta[c, k]   ~ Normal(0, sigma_eta)    // class-topic regression weights
//
//   lambda = H * W                        // Poisson rate
//   counts[d, v] ~ Poisson(lambda[d, v])  // NMF likelihood
//
//   linear_pred_mat = eta * theta'
//   y[d] ~ Categorical(softmax(linear_pred_mat[:, d])) // supervised likelihood
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
  int<lower=1> D;                              // number of observations
  array[D, V] int<lower=0> counts;             // count matrix: D x V

  int<lower=2> C;                              // number of response categories
  array[D] int<lower=1, upper=C> y;            // class label for each observation

  // Hyperparameters
  real<lower=0> shape;                         // Gamma prior shape for H
  real<lower=0> rate;                          // Gamma prior rate  for H, W
  real<lower=0> alpha_beta;                    // Gamma prior shape for W
  real<lower=0> sigma_eta;                     // Normal prior SD for eta
  
  real<lower=0> lambda_ridge_eta;              // L2 penalty on eta
  
  real<lower=0> nmf_weight;                    // multiplicative weight on the
                                               // unsupervised component of the
                                               // log-likelihood
                                               
  real<lower=0> sup_weight;                    // multiplicative weight on the
                                               // supervised component of the
                                               // log-likelihood
}

transformed data{
                                                          
  int<lower=0> NZ = 0;                  // number of non-zero elements in counts
  for(d in 1:D){
    for(v in 1:V){
      if(counts[d, v] > 0){
        NZ += 1;                        // If non-zero, increment NZ
      }
    }
  }
  
  array[NZ] int<lower=1, upper=D> nz_d; // Row indexes of non-zero entries
  array[NZ] int<lower=1, upper=V> nz_v; // Column indexes of non-zero entries
  vector[NZ] nz_counts;                 // non-zero elements of counts
  {
    int idx = 1;
    for(d in 1:D){
      for(v in 1:V){
        if(counts[d, v] > 0){
          nz_d[idx] = d;
          nz_v[idx] = v;
          nz_counts[idx] = counts[d, v];
        idx += 1;
        }
      }
    }
  }
}

parameters{
  matrix<lower=0>[D, K] H;                 // topic loadings: D x K
  matrix<lower=0>[K, V] W;                 // variable-topic loadings: K x V
  matrix[C, K] eta_raw;                    // class-topic weights: C x K
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
  //     HU[d, k] = H[d, k] * u[k]                 // scale by topic weight
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
  array[K] vector[V] beta;                  // not defined as an array of
                                            // simplices for computational
                                            // reasons
  for(k in 1:K){
    beta[k] = W[k, :]' / u[k];
  }
  
  matrix[C, K] eta = sigma_eta * eta_raw;   // scaled version used in likelihood
  
  matrix[D, V] lambda = H * W;              // Poisson rate
  
  matrix[C, D] linear_pred_mat = eta * theta'; // linear predictions for
                                               // supervised likelihood
}

model{
  // ------------------------------------------------------------------
  // Priors
  // ------------------------------------------------------------------
  for(d in 1:D){
    H[d, :] ~ gamma(shape, shape * rate);
  }
  
  for(k in 1:K){
    W[k, :] ~ gamma(alpha_beta, rate);
  }
  
  for(c in 1:C){
    // Prior on eta_raw is standard normal — well-conditioned geometry
    eta_raw[c, :] ~ normal(0, 1);
  }
  
  // Ridge on scaled eta (additional shrinkage beyond the sigma_eta scaling)
  target += -0.5 * lambda_ridge_eta * sum(eta .* eta);

  // ------------------------------------------------------------------
  // NMF likelihood: Poisson with rate = H * W
  // Zeros are skipped — they contribute 0 to the Poisson log-pmf
  // only when the -lambda term is accounted for separately, so we
  // use target += and handle the full expression explicitly.
  // ------------------------------------------------------------------
  
  vector[NZ] lambda_nz;               // Poisson rate values for non-zero counts
  for(i in 1:NZ){
    lambda_nz[i] = lambda[nz_d[i], nz_v[i]];
  }
  
  target += nmf_weight * (dot_product(nz_counts, log(lambda_nz)) - sum(lambda));

  // ------------------------------------------------------------------
  // Supervised likelihood: categorical with softmax linear predictor
  // linear_pred[c] = eta[c, :] * theta[d, :]'
  // ------------------------------------------------------------------
  
  for(d in 1:D){
    target += sup_weight * categorical_logit_lpmf(y[d] | linear_pred_mat[:, d]);
  }
}

generated quantities{
  // ------------------------------------------------------------------
  // Log-likelihoods (for model comparison, LOO-CV, etc.)
  // ------------------------------------------------------------------
  real var_log_lik = 0;
  real response_log_lik = 0;
  real total_log_lik;

  // ------------------------------------------------------------------
  // Posterior predictive
  // ------------------------------------------------------------------
  array[D] int<lower=1, upper=C> y_pred;
  array[D] vector[C] response_probs;
  
  vector[NZ] lambda_nz;
  for(i in 1:NZ){
    lambda_nz[i] = lambda[nz_d[i], nz_v[i]];
  }

  // NMF log-likelihood (sparse: skip zero counts still accounting for -lambda)
  var_log_lik += dot_product(nz_counts, log(lambda_nz)) - sum(lambda);
  
  // Categorical log-likelihood and predictions
  for(d in 1:D){
    response_probs[d] = softmax(linear_pred_mat[:, d]);
    y_pred[d] = categorical_logit_rng(linear_pred_mat[:, d]);
    response_log_lik += categorical_logit_lpmf(y[d] | linear_pred_mat[:, d]);
  }

  total_log_lik = var_log_lik + response_log_lik;

  // ------------------------------------------------------------------
  // Topic correlations (diagnostic: are topics distinguishable?)
  // Uses beta so scale differences don't dominate correlation.
  // ------------------------------------------------------------------
  matrix[K, K] topic_correlations;
  for(k1 in 1:K){
    for(k2 in 1:K){
      if(k1 == k2){
        topic_correlations[k1, k2] = 1.0;
      } else {
        vector[V] dev1 = beta[k1] - 1.0 / V;
        vector[V] dev2 = beta[k2] - 1.0 / V;
        real cov12 = dot_product(dev1, dev2);
        real var1 = dot_self(dev1);
        real var2 = dot_self(dev2);
        topic_correlations[k1, k2] = cov12 / (sqrt(var1) * sqrt(var2));
      }
    }
  }
}

