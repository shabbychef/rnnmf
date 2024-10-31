# /usr/bin/r
#
# Copyright 2024-2024 Steven E. Pav. All Rights Reserved.
# Author: Steven E. Pav 
#
# This file is part of rnnmf.
#
# rnnmf is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# rnnmf is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with rnnmf.  If not, see <http://www.gnu.org/licenses/>.
#
# Created: 2024.08.17
# Copyright: Steven E. Pav, 2024
# Author: Steven E. Pav <steven@gilgamath.com>
# Comments: Steven E. Pav

#' @title gaurnmf .
#'
#' @description 
#'
#' Additive update Non-negative matrix factorization with regularization, general form.
#'
#' @details
#'
#' Attempts to factor given non-negative matrix \eqn{Y} as the product \eqn{LR}
#' of two non-negative matrices. The objective function is Frobenius norm
#' with \eqn{\ell_1} and \eqn{\ell_2} regularization terms.
#' We seek to minimize the objective
#' \deqn{\frac{1}{2}tr((Y-LR)' W_{0R} (Y-LR) W_{0C}) + tr(W_{1L}'L) + tr(W_{1R}'R) + \frac{1}{2} \sum_j tr(L'W_{2RLj}LW_{2CLj}) + tr(R'W_{2RRj}RW_{2CRj}),}
#' subject to \eqn{L \ge 0} and \eqn{R \ge 0} elementwise, 
#' where \eqn{tr(A)} is the trace of \eqn{A}.
#'
#' The code starts from initial estimates and iteratively 
#' improves them, maintaining non-negativity.
#' This implementation uses the Lee and Seung step direction,
#' with a correction to avoid divide-by-zero.
#' The iterative step is optionally re-scaled to take the steepest 
#' descent in the step direction.
#'
#'
#' @param Y  an \eqn{r \times c} matrix to be decomposed.
#' Should have non-negative elements; an error is thrown otherwise.
#' @param L  an \eqn{r \times d} matrix of the initial estimate of L.
#' Should have non-negative elements; an error is thrown otherwise.
#' @param R  an \eqn{d \times c} matrix of the initial estimate of R.
#' Should have non-negative elements; an error is thrown otherwise.
#' @param W_0R  the row space weighting matrix.
#' This should be a positive definite non-negative symmetric \eqn{r \times r} matrix.
#' If omitted, it defaults to the properly sized identity matrix.
#' @param W_0C  the column space weighting matrix.
#' This should be a positive definite non-negative symmetric \eqn{c \times c} matrix.
#' If omitted, it defaults to the properly sized identity matrix.
#' @param W_1L  the \eqn{\ell_1} penalty matrix for the matrix \eqn{R}.
#' If a scalar, corresponds to that scalar times the all-ones matrix.
#' Defaults to all-zeroes matrix, which is no penalty term.
#' @param W_1R  the \eqn{\ell_1} penalty matrix for the matrix \eqn{L}.
#' If a scalar, corresponds to that scalar times the all-ones matrix.
#' Defaults to all-zeroes matrix, which is no penalty term.
#' @param W_2RL  the \eqn{\ell_2} row penalty matrix for the matrix \eqn{L}.
#' If a scalar, corresponds to that scalar times the identity matrix.
#' Can also be a list, in which case \code{W_2CL} must be a list of the same
#' length. The list should consist of \eqn{\ell_2} row penalty matrices.
#' Defaults to all-zeroes matrix, which is no penalty term.
#' @param W_2CL  the \eqn{\ell_2} column penalty matrix for the matrix \eqn{L}.
#' If a scalar, corresponds to that scalar times the identity matrix.
#' Can also be a list, in which case \code{W_2RL} must be a list of the same
#' length. The list should consist of \eqn{\ell_2} column penalty matrices.
#' Defaults to all-zeroes matrix, which is no penalty term.
#' @param W_2RR  the \eqn{\ell_2} row penalty matrix for the matrix \eqn{R}.
#' If a scalar, corresponds to that scalar times the identity matrix.
#' Can also be a list, in which case \code{W_2CR} must be a list of the same
#' length. The list should consist of \eqn{\ell_2} row penalty matrices.
#' Defaults to all-zeroes matrix, which is no penalty term.
#' @param W_2CR  the \eqn{\ell_2} column penalty matrix for the matrix \eqn{R}.
#' If a scalar, corresponds to that scalar times the identity matrix.
#' Can also be a list, in which case \code{W_2RR} must be a list of the same
#' length. The list should consist of \eqn{\ell_2} column penalty matrices.
#' Defaults to all-zeroes matrix, which is no penalty term.
#' @inheritParams aurnmf
#' @return a list with the elements
#' \describe{
#' \item{L}{The final estimate of L.}
#' \item{R}{The final estimate of R.}
#' \item{Lstep}{The infinity norm of the final step in L}.
#' \item{Rstep}{The infinity norm of the final step in R}.
#' \item{iterations}{The number of iterations taken.}
#' \item{converged}{Whether convergence was detected.}
#' }
#' @keywords optimization
#' @template etc
#' @template poc
#' @template ref-merritt
#' @template ref-pav
#' @template ref-leeseung
#' @seealso \code{\link{aurnmf}}
#'
#' @examples 
#'
#'  nr <- 20
#'  nc <- 5
#'  dm <- 2
#'  
#'  randmat <- function(nr,nc,...) { matrix(pmax(0,runif(nr*nc,...)),nrow=nr) }
#'  set.seed(1234)
#'  real_L <- randmat(nr,dm+2)
#'  real_R <- randmat(ncol(real_L),nc)
#'  Y <- real_L %*% real_R
#'  gram_it <- function(G) { t(G) %*% G }
#'  W_0R <- gram_it(randmat(nr+5,nr))
#'  W_0C <- gram_it(randmat(nc+5,nc))
#'  
#'  wt_objective <- function(Y, L, R, W_0R, W_0C) { 
#'    err <- Y - L %*% R
#'    0.5 * sum((err %*% W_0C) * (t(W_0R) %*% err))
#'  }
#'  matrix_trace <- function(G) {
#'    sum(diag(G))
#'  }
#'  wt_objective(Y,real_L,real_R,W_0R,W_0C)
#'  
#'  L_0 <- randmat(nr,dm)
#'  R_0 <- randmat(dm,nc)
#'  wt_objective(Y,L_0,R_0,W_0R,W_0C)
#'  out1 <- gaurnmf(Y, L_0, R_0, W_0R=W_0R, W_0C=W_0C, 
#'          max_iterations=1e4L,check_optimal_step=FALSE)
#'  wt_objective(Y,out1$L,out1$R,W_0R,W_0C)
#'  
#'  W_1L <- randmat(nr,dm)
#'  out2 <- gaurnmf(Y, out1$L, out1$R, W_0R=W_0R, W_0C=W_0C, W_1L=W_1L, 
#'          max_iterations=1e4L,check_optimal_step=FALSE)
#'  wt_objective(Y,out2$L,out2$R,W_0R,W_0C)
#'  
#'  W_1R <- randmat(dm,nc)
#'  out3 <- gaurnmf(Y, out2$L, out2$R, W_0R=W_0R, W_0C=W_0C, W_1R=W_1R, 
#'          max_iterations=1e4L,check_optimal_step=FALSE)
#'  wt_objective(Y,out3$L,out3$R,W_0R,W_0C)
#'
#' \donttest{
#' # example showing how to use the on_iteration_end callback to save iterates.
#'  max_iterations <- 1e3L
#'  it_history <<- rep(NA_real_, max_iterations)
#'  on_iteration_end <- function(iteration, Y, L, R, ...) {
#'    it_history[iteration] <<- wt_objective(Y,L,R,W_0R,W_0C)
#'  }
#'  out1b <- gaurnmf(Y, L_0, R_0, W_0R=W_0R, W_0C=W_0C, 
#'    max_iterations=max_iterations, on_iteration_end=on_iteration_end, check_optimal_step=FALSE)
#' }
#'
#' # should work on sparse matrices too.
#' if (require(Matrix)) { 
#'  real_L <- randmat(nr,dm,min=-1)
#'  real_R <- randmat(dm,nc,min=-1)
#'  Y <- as(real_L %*% real_R, "sparseMatrix")
#'  L_0 <- as(randmat(nr,dm,min=-0.5), "sparseMatrix")
#'  R_0 <- as(randmat(dm,nc,min=-0.5), "sparseMatrix")
#'  out1 <- gaurnmf(Y, L_0, R_0, max_iterations=1e2L,check_optimal_step=TRUE)
#' }
#'
#' @author Steven E. Pav \email{shabbychef@@gmail.com}
#' @export
gaurnmf <- function(Y, L, R, 
										W_0R=NULL, W_0C=NULL, 
										W_1L=0, W_1R=0,
										W_2RL=0, W_2CL=0,
										W_2RR=0, W_2CR=0,
										tau=0.1, annealing_rate=0.01, 
										check_optimal_step=TRUE, 
										zero_tolerance=1e-12,
										max_iterations=1e3L, 
										min_xstep=1e-9,
										on_iteration_end=NULL,
										verbosity=0) {
	stopifnot(all(Y >= 0))
	stopifnot(all(L >= 0))
	stopifnot(all(R >= 0))
	stopifnot((ncol(L)==nrow(R)) && (nrow(L)==nrow(Y)) && (ncol(R)==ncol(Y)))
	stopifnot(missing(W_0R) || is.null(W_0R) || all(W_0R >= 0))
	stopifnot(missing(W_0C) || is.null(W_0C) || all(W_0C >= 0))
	stopifnot(missing(W_1R) || is.null(W_1R) || all(W_1R >= 0))
	stopifnot(missing(W_1L) || is.null(W_1L) || all(W_1L >= 0))
	stopifnot((0 < tau) && (tau < 1))
	stopifnot((0 <= annealing_rate) && (annealing_rate < 1))

	fixd <- fix_LR_names(L, R)
	L <- fixd$L
	R <- fixd$R

	if (is.list(W_2RL) || is.list(W_2CL)) {
		stopifnot(is.list(W_2RL) && is.list(W_2CL) && (length(W_2RL) == length(W_2CL)))
	} else {
		W_2RL <- list(W_2RL)
		W_2CL <- list(W_2CL)
	}
	W_2L_J <- length(W_2RL)

	if (is.list(W_2RR) || is.list(W_2CR)) {
		stopifnot(is.list(W_2RR) && is.list(W_2CR) && (length(W_2RR) == length(W_2CR)))
	} else {
		W_2RR <- list(W_2RR)
		W_2CR <- list(W_2CR)
	}
	W_2R_J <- length(W_2RR)

	# precompute
	W_0R_Y <- (W_0R %**% Y)
	Y_W_0C <- (Y %**% W_0C)

	tau_k <- tau
	finished <- FALSE
	k <- 0
	while (!finished) {
		tau_k <- (1-annealing_rate) * tau_k + annealing_rate

		# update L
		WRt <- (W_0C %**% t(R))
		RWR <- R %*% WRt
		D <- W_1L - W_0R_Y %*% WRt
		F <- (W_0R %**% L) %*% RWR 
		for (jidx in 1:W_2L_J) {
			F <- F + W_2RL[[jidx]] %**% L %**% W_2CL[[jidx]]
		}
		gradfL_k <- D + F
		H_kp1 <- pick_direction(L, gradfL_k, F)
		if (check_optimal_step) {
			K_kp1 <- (W_0R %**% H_kp1) %*% RWR 
			for (jidx in 1:W_2L_J) {
				K_kp1 <- K_kp1 + W_2RL[[jidx]] %**% H_kp1 %**% W_2CL[[jidx]]
			}
		} else {
			K_kp1 <- NULL
		}
		LList <- giqpm_step(L, gradfL_k, H_kp1, K_kp1, tau_k, k, verbosity)
		L <- LList[[1]]
		Lstep <- LList[[2]]

		# update R
		LtW <- t(L) %**% W_0R
		LWL <- LtW %*% L
		D <- W_1R - LtW %*% Y_W_0C
		F <- LWL %*% (R %**% W_0C) 
		for (jidx in 1:W_2R_J) {
			F <- F + W_2RR[[jidx]] %**% R %**% W_2CR[[jidx]]
		}
		gradfR_k <- D + F
		H_kp1 <- pick_direction(R, gradfR_k, F)
		if (check_optimal_step) {
			K_kp1 <- LWL %*% (H_kp1 %**% W_0C) 
			for (jidx in 1:W_2R_J) {
				K_kp1 <- K_kp1 + W_2RR[[jidx]] %**% H_kp1 %**% W_2CR[[jidx]]
			}
		} else {
			K_kp1 <- NULL
		}
		RList <- giqpm_step(R, gradfR_k, H_kp1, K_kp1, tau_k, k, verbosity)
		R <- RList[[1]]
		Rstep <- RList[[2]]
		R[R <= zero_tolerance] <- 0
		L[L <= zero_tolerance] <- 0
		k <- k + 1
		converged <- (max(c(Lstep,Rstep)) < min_xstep)
		finished <- converged || (k >= max_iterations) 
		if (! is.null(on_iteration_end)) {
			on_iteration_end(iteration=k, Y=Y, L=L, R=R, Lstep=Lstep, Rstep=Rstep, converged=converged, finished=finished)
		}
	}
	if (verbosity > 1) { 
		print(paste0("terminated after ",k," iterations. converged: ",converged))
	}
	rownames(L) <- rownames(Y)
	colnames(R) <- colnames(Y)
	return(list(L=L,R=R,iterations=k,converged=converged,Lstep=Lstep,Rstep=Rstep))
}

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
