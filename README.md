# BIDA
Code for running linear causal effect estimation using Bayesian IDA (BIDA).

## APS solver

## Running BIDA from R

After compiling the APS solver, BIDA can be run from R. As an example case, we will use the included data, which contains 200 observations from a Gaussian DAG model over 10 variables. The data has been zero-centered and standardized. 

Assuming the BIDA folder is set as the current working directory in R, we read in the data and true causal effect matrix. In the true effect matrix, a value at row i and column j will correspond to the causal effect of i on j.

``` r
data <- read.csv(file = 'example_data/data_d10.txt', header = FALSE)
true_effects <- as.matrix(read.csv(file = 'example_data/tce_d10.txt', header = FALSE))
```

We then source the R-files.

``` r
file_paths <- list.files(pattern = "[.]R$", path = "R/", full.names = TRUE)
invisible(sapply(file_paths, source))
```

``` r
We compute the BIDA posteriors, setting the maximum parent set size to 5.
bida_post <- bida(data, max_parent_size = 5)
```

Finally, we evaluate our estimates by calculating the mean squared error (MSE) between the posterior means and true causal effects. 

``` r
# Calculate mean of BIDA posteriors
bida_mean <- calc_bida_post_mean(bida_post)

# Calculate mean squared error for the mean posterior point estimates 
mse <- mean((bida_mean-true_effects)[diag(ncol(data)) == 0]^2)
```

This should give an MSE around 0.0655. Alternatively, if one is interested in ranking the causal effect by their magnitude, one can use calc_bida_post_mean_abs to calculate (numerically) the posterior absolute mean.