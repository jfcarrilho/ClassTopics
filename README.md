<<<<<<< HEAD
# ClassTopics

<!-- badges: start -->
[![R-CMD-check](https://github.com/jfcarrilho/ClassTopics/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jfcarrilho/ClassTopics/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`ClassTopics` fits supervised topic models to count data (e.g. gene
expression) jointly with a binary or categorical response. It combines
non-negative matrix factorization (NMF) with a supervised latent Dirichlet
allocation (LDA)-style reparameterization, estimated via full Bayesian
inference in [Stan](https://mc-stan.org/) (through
[`cmdstanr`](https://mc-stan.org/cmdstanr/)), with an EM warm-start used to
initialize the sampler.

The model jointly learns:

- **Topics**: interpretable, non-negative variable signatures shared
  across observations (`beta`), and each observation's topic proportions
  (`theta`).
- **A classifier**: regression coefficients (`eta`) linking topic proportions
  to a binary or multiclass outcome, fit simultaneously with the topics so
  that the discovered topics are informed by the response they need to predict.

## Installation

`ClassTopics` is not on CRAN. Install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("jfcarrilho/ClassTopics")
```

### Stan / cmdstanr setup

`ClassTopics` fits its models with [`cmdstanr`](https://mc-stan.org/cmdstanr/),
which is not on CRAN and must be installed separately:

```r
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
cmdstanr::install_cmdstan()
```

You only need to do this once per machine. The first time you call any
`ClassTopics` fitting function, the package's bundled Stan models are
compiled automatically and the resulting executables are cached for reuse
on subsequent calls -- you do not need to compile anything by hand.

## Quick example

```r
library(ClassTopics)

# counts:   a D x V matrix of non-negative counts (observations x variables)
# response: a length-D factor/character vector with 2+ categories
fit <- ClassTopics(
  counts   = counts,
  response = response,
  K        = 5,     # number of topics
  seed     = 123
)

results <- ClassTopics_results(fit, true_response = response)

# Diagnostic plots
plot_topic_correlations(results)
plot_topic_response_heatmap(results)
plot_top_vars(results, n_vars = 5)
```

## Cross-validation

```r
cv_fit <- cv_ClassTopics(
  counts   = counts,
  response = response,
  k_folds  = 5,
  K_topics = 5,
  seed     = 123
)

cvTestAccuracy(cv_fit)
overfittingGap(cv_fit)
```

See `vignette("ClassTopics-intro")` for a complete walkthrough, including
prediction on new samples and interpreting model output.

## Learning more

- `vignette("ClassTopics-intro")` -- full workflow: fitting, cross-validation,
  prediction, and interpreting results
- `?ClassTopics` -- main fitting function
- `?cv_ClassTopics` -- k-fold cross-validation
- `?ClassTopics_results` -- extract and summarize a fitted model

## Citation

If you use `ClassTopics` in your research, please cite the associated
publication (details to be added).

## License

MIT © João F. Carrilho, Marta B. Lopes, Susan P. Holmes
=======
# ClassTopics

<!-- badges: start -->
<!--[![R-CMD-check](https://github.com/jfcarrilho/ClassTopics/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jfcarrilho/ClassTopics/actions/workflows/R-CMD-check.yaml)-->
<!-- badges: end -->
 
`ClassTopics` fits supervised topic models to count data (e.g. gene
expression) jointly with a binary or categorical response. It combines
non-negative matrix factorization (NMF) with a supervised latent Dirichlet
allocation (LDA)-style reparameterization, estimated via full Bayesian
inference in [Stan](https://mc-stan.org/) (through
[`cmdstanr`](https://mc-stan.org/cmdstanr/)), with an EM warm-start used to
initialize the sampler.
 
The model jointly learns:
 
- **Topics**: interpretable, non-negative gene expression signatures shared
  across patients (`beta`), and each patient's topic proportions (`theta`).
- **A classifier**: regression coefficients (`eta`) linking topic proportions
  to a binary or multiclass outcome, fit simultaneously with the topics so
  that the discovered topics are informed by the response they need to predict.
## Installation
 
`ClassTopics` is not on CRAN. Install the development version from GitHub:
 
```r
# install.packages("remotes")
remotes::install_github("jfcarrilho/ClassTopics")
```
 
### Stan / cmdstanr setup
 
`ClassTopics` fits its models with [`cmdstanr`](https://mc-stan.org/cmdstanr/),
which is not on CRAN and must be installed separately:
 
```r
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
cmdstanr::install_cmdstan()
```
 
You only need to do this once per machine. The first time you call any
`ClassTopics` fitting function, the package's bundled Stan models are
compiled automatically and the resulting executables are cached for reuse
on subsequent calls -- you do not need to compile anything by hand.
 
## Quick example
 
```r
library(ClassTopics)
 
# gene_counts: a D x V matrix of non-negative counts (patients x genes)
# response:    a length-D factor/character vector with 2+ categories
fit <- ClassTopics(
  gene_counts = gene_counts,
  response    = response,
  K           = 5,     # number of topics
  seed        = 123
)
 
results <- ClassTopics_results(fit, true_response = response)
 
# Diagnostic plots
plot_topic_correlations(results)
plot_topic_response_heatmap(results)
plot_top_genes(results, n_genes = 5)
```
 
## Cross-validation
 
```r
cv_fit <- cv_ClassTopics(
  gene_counts = gene_counts,
  response    = response,
  k_folds     = 5,
  K_topics    = 5,
  seed        = 123
)
 
cvTestAccuracy(cv_fit@cv_predictions)
overfittingGap(cv_fit@cv_predictions)
```
 
See `vignette("ClassTopics-intro")` for a complete walkthrough, including
prediction on new samples and interpreting model output.
 
## Learning more
 
- `vignette("ClassTopics-intro")` -- full workflow: fitting, cross-validation,
  prediction, and interpreting results
- `?ClassTopics` -- main fitting function
- `?cv_ClassTopics` -- k-fold cross-validation
- `?ClassTopics_results` -- extract and summarize a fitted model
## Citation
 
If you use `ClassTopics` in your research, please cite the associated
publication (details to be added).
 
## License
 
MIT © João F. Carrilho, Marta B. Lopes, Susan P. Holmes
>>>>>>> 81565ab30b70382ad2a484a8f5c90dd9a8ba4e20
