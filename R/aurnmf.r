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

giqpm_step <- function(X_k, gradfX_k, H_kp1, K_kp1, tau_k, k, verbosity) {
		# check for positive or zero values in hstep?
		checkus <- H_kp1 < 0
		if (any(checkus)) {
			alpha_hat_k <- min(-X_k[as.matrix(checkus)] / H_kp1[as.matrix(checkus)])
		} else if (any(H_kp1 != 0)) {
			alpha_hat_k <- Inf
		} else {
			# all hstep are zero. converged.
			converged <- TRUE
			return(list(X_k,0))
		}

		direx <- sum(gradfX_k * H_kp1)
		if (direx >= 0) {
			warning("step direction misaligned to gradient?")
		}
		if (!is.null(K_kp1)) {
			denom <- sum(H_kp1 * K_kp1)
			if (denom <= 0) {
				warning("violation of PD in G?")
			}
			alpha_star_k <- -direx / denom
		} else {
			alpha_star_k <- 1
		}
		alpha_k <- min(tau_k * alpha_hat_k, alpha_star_k)
		if (verbosity > 0) { print(paste0("iteration ", k, " taking step of size ",alpha_k,".")) }
		fullstep <- alpha_k * H_kp1
		return(list(X_k + fullstep,max(abs(fullstep))))
}
pick_direction <- function(X_k, gradfX_k, F) {
	# tweak the step:
	# if you would have division by zero, do not do that, rather just step
	# in negative gradient direction.
	# moreover, if X_k is zero and the gradient is negative, allow a step in the
	# positive direction.
	Fzer <- F <= 0
	Xzer <- X_k <= 0
	grad_neg <- gradfX_k < 0
	F[Fzer] <- 1
	H_kp1 <- - gradfX_k * X_k / F
	# H_kp1[Fzer & Xzer] <- pmax(-gradfX_k[Fzer & Xzer], 0)
	H_kp1[Xzer] <- pmax(-gradfX_k[Xzer], 0)
	return(H_kp1)
}

# like matrix multiplication, but if X or Y is not given or is null, treat as the identity.
# if one of them is a scalar, treat them as the identity times that scalar.
# 2FIX: do we want to allow a vector to stand for a diagonal matrix?
`%**%` <- function(X, Y) {
	if (missing(X) || is.null(X)) {
		return(Y)
	} else if (missing(Y) || is.null(Y)) {
		return(X)
	} else if (is.null(dim(X)) || is.null(dim(Y))) {
		return(X * Y)
	}
	return(X %*% Y)
}
tish <- function(X) {
	if (missing(X) || is.null(X)) {
		return(NULL)
	} else {
		return(t(X))
	}
}

# multiplies matrix M times (11'-I)
times_orth <- function(M) {
	# alternatively: M %*% (1 - diag(ncol(M)))
	matrix(rep(rowSums(M),ncol(M)),ncol=ncol(M),byrow=FALSE) - M
}

fix_LR_names <- function(L, R) {
	if (is.null(colnames(L)) & is.null(rownames(R))) {
		new_names <- paste0('factor',seq_len(ncol(L)))
		colnames(L) <- new_names
		rownames(R) <- new_names
	}
	if (is.null(colnames(L))) {
		colnames(L) <- rownames(R)
	}
	if (is.null(rownames(L))) {
		rownames(R) <- colnames(L)
	}
	return(list(L=L,R=R))
}


#' @title nmf .
#'
#' @description 
#'
#' Additive update Non-negative matrix factorization with regularization.
#'
#' @details
#'
#' Attempts to factor given non-negative matrix \eqn{Y} as the product \eqn{LR}
#' of two non-negative matrices. The objective function is Frobenius norm
#' with \eqn{\ell_1} and \eqn{\ell_2} regularization terms.
#' We seek to minimize the objective
#' \deqn{\frac{1}{2}tr((Y-LR)' W_{0R} (Y-LR) W_{0C}) + \lambda_{1L} |L| + \lambda_{1R} |R| + \frac{\lambda_{2L}}{2} tr(L'L) + \frac{\lambda_{2R}}{2} tr(R'R) + \frac{\gamma_{2L}}{2} tr((L'L) (11' - I)) + \frac{\gamma_{2R}}{2} tr((R'R) (11' - I)),}
#' subject to \eqn{L \ge 0} and \eqn{R \ge 0} elementwise, 
#' where \eqn{|A|} is the sum of the elements of \eqn{A} and 
#' \eqn{tr(A)} is the trace of \eqn{A}.
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
#' @param lambda_1L  the scalar \eqn{\ell_1} penalty for the matrix \eqn{L}.
#' Defaults to zero.
#' @param lambda_1R  the scalar \eqn{\ell_1} penalty for the matrix \eqn{R}.
#' Defaults to zero.
#' @param lambda_2L  the scalar \eqn{\ell_2} penalty for the matrix \eqn{L}.
#' Defaults to zero.
#' @param lambda_2R  the scalar \eqn{\ell_2} penalty for the matrix \eqn{R}.
#' Defaults to zero.
#' @param gamma_2L  the scalar \eqn{\ell_2} penalty for non-orthogonality of the matrix \eqn{L}.
#' Defaults to zero.
#' @param gamma_2R  the scalar \eqn{\ell_2} penalty for non-orthogonality of the matrix \eqn{R}.
#' Defaults to zero.
#' @param on_iteration_end  an optional function that is called at the end of
#' each iteration. The function is called as 
#' \code{on_iteration_end(iteration=iteration, Y=Y, L=L, R=R, Lstep=Lstep, Rstep=Rstep, ...)}
#'
#' @inheritParams giqpm
#' @return a list with the elements
#' \describe{
#' \item{L}{The final estimate of L.}
#' \item{R}{The final estimate of R.}
#' \item{Lstep}{The infinity norm of the final step in L.}
#' \item{Rstep}{The infinity norm of the final step in R.}
#' \item{iterations}{The number of iterations taken.}
#' \item{converged}{Whether convergence was detected.}
#' }
#' @keywords optimization
#' @template etc
#' @template poc
#' @template ref-merritt
#' @template ref-pav
#' @template ref-leeseung
#' @seealso \code{\link{gaurnmf}}, \code{\link{murnmf}}.
#'
#' @examples 
#'
#'  nr <- 100
#'  nc <- 20
#'  dm <- 4
#' 
#'  randmat <- function(nr,nc,...) { matrix(pmax(0,runif(nr*nc,...)),nrow=nr) }
#'  set.seed(1234)
#'  real_L <- randmat(nr,dm)
#'  real_R <- randmat(dm,nc)
#'  Y <- real_L %*% real_R
#' # without regularization
#'  objective <- function(Y, L, R) { sum((Y - L %*% R)^2) }
#'  objective(Y,real_L,real_R)
#' 
#'  L_0 <- randmat(nr,dm)
#'  R_0 <- randmat(dm,nc)
#'  objective(Y,L_0,R_0)
#'  out1 <- aurnmf(Y, L_0, R_0, max_iterations=5e3L,check_optimal_step=FALSE)
#'  objective(Y,out1$L,out1$R)
#' # with L1 regularization on one side
#'  out2 <- aurnmf(Y, L_0, R_0, lambda_1L=0.1, max_iterations=5e3L,check_optimal_step=FALSE)
#' # objective does not suffer because all mass is shifted to R
#'  objective(Y,out2$L,out2$R)
#' list(L1=sum(out1$L),R1=sum(out1$R),L2=sum(out2$L),R2=sum(out2$R))
#' sum(out2$L)
#' # with L1 regularization on both sides
#'  out3 <- aurnmf(Y, L_0, R_0, lambda_1L=0.1,lambda_1R=0.1,
#'      max_iterations=5e3L,check_optimal_step=FALSE)
#' # with L1 regularization on both sides, raw objective suffers
#'  objective(Y,out3$L,out3$R)
#' list(L1=sum(out1$L),R1=sum(out1$R),L3=sum(out3$L),R3=sum(out3$R))
#'
#' \donttest{
#' # example showing how to use the on_iteration_end callback to save iterates.
#' max_iterations <- 5e3L
#' it_history <<- rep(NA_real_, max_iterations)
#' quadratic_objective <- function(Y, L, R) { sum((Y - L %*% R)^2) }
#' on_iteration_end <- function(iteration, Y, L, R, ...) {
#'   it_history[iteration] <<- quadratic_objective(Y,L,R)
#' }
#' out1b <- aurnmf(Y, L_0, R_0, max_iterations=max_iterations, on_iteration_end=on_iteration_end)
#' }
#'
#' # should work on sparse matrices too.
#' if (require(Matrix)) { 
#'  real_L <- randmat(nr,dm,min=-1)
#'  real_R <- randmat(dm,nc,min=-1)
#'  Y <- as(real_L %*% real_R, "sparseMatrix")
#'  L_0 <- as(randmat(nr,dm,min=-0.5), "sparseMatrix")
#'  R_0 <- as(randmat(dm,nc,min=-0.5), "sparseMatrix")
#'  out1 <- aurnmf(Y, L_0, R_0, max_iterations=1e2L,check_optimal_step=TRUE)
#' }
#'
#' @importFrom Matrix t rowSums
#' @author Steven E. Pav \email{shabbychef@@gmail.com}
#' @export
aurnmf <- function(Y, L, R, 
									 W_0R=NULL, W_0C=NULL, 
									 lambda_1L=0, lambda_1R=0, 
									 lambda_2L=0, lambda_2R=0, 
									 gamma_2L=0,  gamma_2R=0,
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
	stopifnot(lambda_1L >= 0)
	stopifnot(lambda_1R >= 0)
	stopifnot(lambda_2L >= 0)
	stopifnot(lambda_2R >= 0)
	stopifnot(gamma_2L >= 0)
	stopifnot(gamma_2R >= 0)
	stopifnot((0 < tau) && (tau < 1))
	stopifnot((0 <= annealing_rate) && (annealing_rate < 1))

	fixd <- fix_LR_names(L, R)
	L <- fixd$L
	R <- fixd$R

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
		D <- lambda_1L - W_0R_Y %*% WRt
		F <- (W_0R %**% L) %*% RWR + lambda_2L * L
		if (gamma_2L > 0) {
			F <- F + gamma_2L * times_orth(L)
		}
		gradfL_k <- D + F
		H_kp1 <- pick_direction(L, gradfL_k, F)
		if (check_optimal_step) {
			K_kp1 <- (W_0R %**% H_kp1) %*% RWR + lambda_2L * H_kp1
		} else {
			K_kp1 <- NULL
		}
		LList <- giqpm_step(L, gradfL_k, H_kp1, K_kp1, tau_k, k, verbosity)
		L <- LList[[1]]
		Lstep <- LList[[2]]
		# update R
		LtW <- t(L) %**% W_0R
		LWL <- LtW %*% L
		D <- lambda_1R - LtW %*% Y_W_0C
		F <- LWL %*% (R %**% W_0C) + lambda_2R * R
		if (gamma_2R > 0) {
			F <- F + gamma_2R * times_orth(R)
		}
		gradfR_k <- D + F
		H_kp1 <- pick_direction(R, gradfR_k, F)
		if (check_optimal_step) {
			K_kp1 <- LWL %*% (H_kp1 %**% W_0C) + lambda_2R * H_kp1
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
