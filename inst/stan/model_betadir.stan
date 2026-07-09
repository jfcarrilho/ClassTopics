// =============================================================================
// Supervised NMF / sLDA Model
// =============================================================================
//
// Generative model:
//   H[d, k] ~ Gamma(shape, rate)          // topic loadings per patient
//   beta[k, v]  ~ Dirichlet(alpha_beta)   // gene weights per topic
//   eta[c, k]   ~ Normal(0, sigma_eta)    // class-topic regression weights
//
//   lambda[d, v]  = dot_product(H[d,:], W[:,v])   // Poisson rate
//   gene_counts[d, v] ~ Poisson(lambda[d, v])     // NMF likelihood
//
//   linear_pred[d, c] = eta[c,:] * theta[d,:]'
//   y[d] ~ Categorical(softmax(linear_pred[d]))   // supervised likelihood
//
// Identifiability:
//   - H and W are non-negative; their product defines the Poisson rate
//   - eta is fully free (all C rows estimated); normal(0, sigma_eta) prior
//     resolves the softmax non-identifiability, mirroring glmnet's approach
//   - theta is derived as a transformed parameters and is the interpretable
//     topic representation along with beta
//   - W is estimated through beta to control sparsity within topics
// =============================================================================

data{
  int<lower=2> K;                              // number of topics
  int<lower=2> V;                              // number of genes
  int<lower=1> D;                              // number of patients
  array[D, V] int<lower=0> gene_counts;        // gene count matrix: D x V

  int<lower=2> C;                              // number of response categories
  array[D] int<lower=1, upper=C> y;            // class label for each patient

  // Hyperparameters
  real<lower=0> shape;                         // Gamma prior shape for H
  real<lower=0> rate;                          // Gamma prior rate for H
  real<lower=0> alpha_beta;                    // Dirichlet prior concentration
                                               // for beta
  real<lower=0> sigma_eta;                     // Normal prior SD for eta
  
  real<lower=0> lambda_cat;                    // likelihood weight for
                                               // categorical component
  real<lower=0> lambda_ridge_eta;              // L2 penalty on eta
}

parameters{
  matrix<lower=0>[D, K] H;            // topic loadings:      D x K
  array[K] simplex[V] beta;           // gene-topic weights:  K x V
  vector<lower=0>[K] u;               // overall magnitude of each topic
  matrix[C, K] eta_raw;               // class-topic weights: C x K
}

transformed parameters{
  matrix<lower=0>[K, V] W;            // row-wise scaled beta matrix used in NMF
  for(k in 1:K){
    W[k, :] = u[k] * beta[k]';
  }

  // ------------------------------------------------------------------
  // theta_norm : D x K  topic proportions per patient
  //   Obtained according to the Poisson Non-negative Matrix Factorization to 
  //   Multinomial Topic Model reparameterization
  //   (see Carbonetto et al. 2021).
  //
  //   Procedure:
  //     HU[d, k]    = H[d, k] * u[k]              // scale by topic weight
  //     theta[d, k] = HU[d, k] / sum_k HU[d, k']  // row-normalise
  //
  //   Rows of theta sum to 1 and are interpretable as the fraction
  //   of each patient's gene expression explained by each topic.
  // ------------------------------------------------------------------
  matrix[D, K] theta;
  for(d in 1:D){
    vector[K] HU_d = H[d, :]' .* u;     // elementwise: K-vector
    theta[d, :] = HU_d' / sum(HU_d);
  }
  
  matrix[C, K] eta = sigma_eta * eta_raw;   // scaled version used in likelihood
}

model{
  // ------------------------------------------------------------------
  // Priors
  // ------------------------------------------------------------------
  for(d in 1:D){
    H[d, :] ~ gamma(shape, rate);
  }
  
  for(k in 1:K){
    beta[k] ~ dirichlet(rep_vector(alpha_beta, V));
  }
  
  u ~ gamma(shape, rate);
  
  for(c in 1:C){
    // Prior on eta_raw is standard normal — well-conditioned geometry
    eta_raw[c, :] ~ normal(0, 1);
  }
  
  // Ridge on scaled eta (additional shrinkage beyond the sigma_eta scaling)
  target += -0.5 * lambda_ridge_eta * sum(eta .* eta);

  // ------------------------------------------------------------------
  // NMF likelihood: Poisson with rate = theta * beta
  // Zeros are skipped — they contribute 0 to the Poisson log-pmf
  // only when the -lambda term is accounted for separately, so we
  // use target += and handle the full expression explicitly.
  // ------------------------------------------------------------------
  for(d in 1:D){
    target += -dot_product(H[d, :], W * rep_vector(1.0, V));
    for(v in 1:V){
      if (gene_counts[d, v] > 0){
        real lambda_dv = dot_product(H[d, :], W[:, v]);
        target += gene_counts[d, v] * log(lambda_dv);
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
  real gene_log_lik     = 0;
  real response_log_lik = 0;
  real total_log_lik;

  // ------------------------------------------------------------------
  // Posterior predictive
  // ------------------------------------------------------------------
  array[D] int<lower=1, upper=C> y_pred;
  array[D] vector[C]             response_probs;

  // NMF log-likelihood (sparse: skip zero counts still accounting for -lambda)
  for(d in 1:D){
    gene_log_lik += -dot_product(H[d, :], W * rep_vector(1.0, V));
    for(v in 1:V){
      if (gene_counts[d, v] > 0){
        real lambda_dv = dot_product(H[d, :], W[:, v]);
        gene_log_lik += gene_counts[d, v] * log(lambda_dv);
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

  total_log_lik = gene_log_lik + response_log_lik;

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
        real mean1 = mean(beta[k1, :]);
        real mean2 = mean(beta[k2, :]);
        vector[V] dev1 = beta[k1] - mean1;
        vector[V] dev2 = beta[k2] - mean2;
        real cov12     = dot_product(dev1, dev2);
        real var1      = dot_self(dev1);
        real var2      = dot_self(dev2);
        topic_correlations[k1, k2] = cov12 / (sqrt(var1) * sqrt(var2));
      }
    }
  }
}

