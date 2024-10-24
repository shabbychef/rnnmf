

# rnnmf



Implements regularized non-negative matrix factorization by a method similar to 
Lee & Seung, "Algorithms for Non-negative Matrix Factorization," 2001.

-- Steven E. Pav, shabbychef@gmail.com

## Installation

This package may be installed from CRAN; the latest version may be
found on [github](https://github.com/shabbychef/rnnmf "rnnmf")
via devtools, or installed via [drat](https://github.com/eddelbuettel/drat "drat"):


``` r
# CRAN
install.packages(c("rnnmf"))
# devtools
if (require(devtools)) {
    # latest greatest
    install_github("shabbychef/rnnmf")
}
# via drat:
if (require(drat)) {
    drat:::add("shabbychef")
    # not yet: install.packages('rnnmf')
}
```

# Basic Usage

We demonstrate the usage of the multiplicative and additive updates in
factoring a small matrix which we constructed to be the product of two
reduced rank non-negative matrices.


``` r
library(dplyr)
library(rnnmf)
library(ggplot2)

frobenius_norm_err <- function(Y, L, R) {
    sqrt(sum(abs(Y - L %*% R)^2))
}
runifmat <- function(nr, nc, ...) {
    matrix(pmax(0, runif(nr * nc, ...)), nrow = nr)
}
test_a_bunch <- function(Y_t, L_0, R_0, niter = 10000L) {
    history <<- rep(NA_real_, niter)
    on_iteration_end <- function(iteration, Y, L, R,
        ...) {
        history[iteration] <<- frobenius_norm_err(Y,
            L, R)
    }
    wuz <- aurnmf(Y_t, L_0, R_0, max_iterations = length(history),
        on_iteration_end = on_iteration_end)
    df1 <- tibble(x = seq_along(history), y = history) %>%
        mutate(method = "additive, optimal step")

    history <<- rep(NA_real_, niter)
    wuz <- murnmf(Y_t, L_0, R_0, max_iterations = length(history),
        on_iteration_end = on_iteration_end)
    df2 <- tibble(x = seq_along(history), y = history) %>%
        mutate(method = "multiplicative")

    retv <- bind_rows(df1, df2) %>%
        mutate(nr = nrow(Y_t), nc = ncol(Y_t), nd = ncol(L_0),
            max_iter = niter)
    return(retv)
}

nr <- 30
nc <- 8
nd <- 3
set.seed(1234)
L_t <- runifmat(nr, nd)
R_t <- runifmat(nd, nc)
Y_t <- L_t %*% R_t

L_0 <- runifmat(nrow(Y_t), nd + 1)
R_0 <- runifmat(ncol(L_0), ncol(Y_t))

test_a_bunch(Y_t, L_0, R_0, niter = 10000L) %>%
    ggplot(aes(x, y, color = method)) + geom_line() +
    scale_x_log10(labels = scales::comma) + scale_y_log10() +
    labs(x = "Step", y = expression(L[2] ~ ~Error),
        title = "Frobenius Norm of Error vs Step",
        color = "Method", caption = paste0("Factoring ",
            nr, " x ", nc, " matrix down to ", nd,
            " dimensions."))
```

<div class="figure">
<img src="tools/figure/basic_simulations-1.png" alt="plot of chunk basic_simulations" width="600px" height="500px" />
<p class="caption">plot of chunk basic_simulations</p>
</div>


## See also

* The original paper, by Lee, Daniel D. and Seung, H. Sebastian.
	[Algorithms for Non-negative Matrix Factorization](http://papers.nips.cc/paper/1861-algorithms-for-non-negative-matrix-factorization.pdf), 2001.
* Pav, Steven E. [System and method for unmixing spectroscopic observations with nonnegative matrix factorization](https://patentscope.wipo.int/search/en/detail.jsf?docId=US42758160), 2012.

