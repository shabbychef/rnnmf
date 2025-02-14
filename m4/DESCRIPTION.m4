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
Title: Regularized Non-Negative Matrix Factorization
BugReports: https://github.com/shabbychef/PKG_NAME()/issues
Description: A proof of concept implementation of regularized non-negative matrix factorization optimization.
    A non-negative matrix factorization factors non-negative matrix Y approximately as L R, for non-negative
    matrices L and R of reduced rank. This package supports such factorizations with weighted objective and
    regularization penalties. Allowable regularization penalties include L1 and L2 penalties on L and R,
    as well as non-orthogonality penalties. This package provides multiplicative update algorithms, which are
    a modification of the algorithm of Lee and Seung (2001)
    <http://papers.nips.cc/paper/1861-algorithms-for-non-negative-matrix-factorization.pdf>, as well
    as an additive update derived from that multiplicative update.  See also Pav (2024) <doi:10.48550/arXiv.2410.22698>.
Depends: 
    R (>= 3.0.2)
Imports:
    Matrix
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
