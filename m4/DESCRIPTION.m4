dnl divert here just means the output from basedefs does not appear.
divert(-1)
include(basedefs.m4)
divert(0)dnl
Package: PKG_NAME()
Maintainer: Steven E. Pav <shabbychef@gmail.com>
Authors@R: c(person(c("Steven", "E."), "Pav", 
    role=c("aut","cre"),
    email="shabbychef@gmail.com",
    comment = c(ORCID = "0000-0002-4197-6195")))
Version: VERSION()
Date: DATE()
License: LGPL-3
Title: Regularized Non-negative Matrix Factorization
BugReports: https://github.com/shabbychef/PKG_NAME()/issues
Description: A proof of concept implementation of regularized non-negative matrix factorization optimization.
Depends: 
    R (>= 3.0.2)
dnl Imports:
dnl LinkingTo: Rcpp
Suggests: 
    testthat, 
    dplyr,
    ggplot2,
    scales,
    viridis,
dnl cocktailApp,
    knitr
URL: https://github.com/shabbychef/PKG_NAME()
VignetteBuilder: knitr
Collate:
m4_R_FILES()
dnl vim:ts=2:sw=2:tw=79:syn=m4:ft=m4:et
