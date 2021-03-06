# BIDA

Code for causal effect discovery and linear causal effect estimation using Bayesian IDA (BIDA). The code package contains an C++ implementation (APS) for exact computation of parent set and ancestor relation posterior probabilities using dynamic programming (see Algorithm 1 in Pensar et al. [1]). The computationally less demanding part of the effect estimation method has been implemented in R. 

The repo also contains the supplementary material for the accepted paper. 

## APS

The APS tool can be used for two related but separate problems:

* For each variable and parent set, compute the total weight of all DAGs where a node has a particular parent set.

* For each pair `(i, j)` of variables, compute the total weight of DAGs where `j` is an ancestor of `i`.

To compile the program, go to the aps subdirectory and run `make`. 

After compilation, use `aps help` for more details about the tool.

## Calculating BIDA posteriors

After compiling the APS program, BIDA can be run from R. As an example case, we will use the included data, which contains 200 observations from a Gaussian DAG model over 10 variables. The data has been zero-centered and standardized. 

Assuming the BIDA folder is set as the current working directory in R, we read in the data and true causal effect matrix. 

``` r
data <- read.csv(file = 'example_data/data_d10.txt', header = FALSE)
true_effects <- as.matrix(read.csv(file = 'example_data/tce_d10.txt', header = FALSE))
```

In the true effect matrix, a value at row i and column j will correspond to the causal effect of i on j.

We then source the R-files.

``` r
file_paths <- list.files(pattern = "[.]R$", path = "R/", full.names = TRUE)
invisible(sapply(file_paths, source))
```

We compute the BIDA posteriors, setting the maximum parent set size to 5.

``` r
bida_post <- bida(data, max_parent_size = 5)
```

Finally, we evaluate our estimates by calculating the mean squared error (MSE) between the posterior means and true causal effects. 

``` r
# Calculate mean of BIDA posteriors
bida_mean <- calc_bida_post_mean(bida_post)

# Calculate mean squared error for the mean posterior point estimates 
mse <- mean((bida_mean-true_effects)[diag(ncol(data)) == 0]^2)
```

This should give an MSE around 0.0655. 

Alternatively, if one is interested in ranking the causal effect by their magnitude, one can calculate (numerically) the posterior absolute mean.

``` r
bida_mean_abs <- calc_bida_post_mean_abs(bida_post)
```

## Calculating Ancestor Relation Probabilities (ARPs)

Continuing with our example above, we compute the ARPs which is a more direct Bayesian approach for estimating if there exists a causal relation or not. Again, the maximum parent set size is set to 5. 

``` r
arp <- calc_arp(data, max_parent_size = 5)
```

In the resulting matrix, the value at row `i` and column `j` is the posterior probability of there existing an ancestral path from `i` to `j`.

## Reference

BIDA was developed as part of an academic project: 

[1] Pensar, J., Talvitie, T., Hyttinen, A., Koivisto, M. A Bayesian Approach for Estimating Causal Effects from Observational Data. Accepted to *AAAI-20*.

Please cite the above paper when using the method or parts of it (modified or as is).



