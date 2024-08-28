# /usr/bin/r
#
# Copyright 2024-2024 Steven E. Pav. All Rights Reserved.
# Author: Steven E. Pav 
#
# This file is part of rnmf.
#
# rnmf is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# rnmf is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with rnmf.  If not, see <http://www.gnu.org/licenses/>.
#
# Created: 2024.08.17
# Copyright: Steven E. Pav, 2024
# Author: Steven E. Pav <steven@gilgamath.com>
# Comments: Steven E. Pav

gipm_step <- function(X_k, gradfX_k, H_kp1, K_kp1, tau_k, k, verbosity) {
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
`%**%` <- function(X, Y) {
	if (missing(X) || is.null(X)) {
		return(Y)
	} else if (missing(Y) || is.null(Y)) {
		return(X)
	} else if (!is.matrix(X) || !is.matrix(Y)) {
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
#' @title nmf .
#'
#' @description 
#'
#' Non-negative matrix factorization with regularization.
#'
#' @details
#'
#' Attempts to factor given non-negative matrix \eqn{Y} as the product \eqn{LR}
#' of two non-negative matrices. The objective function is Frobenius norm
#' with \eqn{\ell_1} and \eqn{\ell_2} regularization terms.
#' We seek to minimize the objective
#' \deqn{\frac{1}{2}tr((Y-LR)' W_{0R} (Y-LR) W_{0C}) + \lambda_{1L} |L| + \lambda_{1R} |R| + \frac{\lambda_{2L}}{2} tr(L'L) + \frac{\lambda_{2R}}{2} tr(R'R),}
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
#' @inheritParams gipm
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
#' @seealso \code{\link{gnmf}}
#'
#' @examples 
#'
#'  nr <- 100
#'  nc <- 20
#'  dm <- 4
#' 
#'  randmat <- function(nr,nc) { matrix(runif(nr*nc),nrow=nr) }
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
#'  out1 <- nmf(Y, L_0, R_0, max_iterations=5e4L,check_optimal_step=FALSE)
#'  objective(Y,out1$L,out1$R)
#' # with L1 regularization on one side
#'  out2 <- nmf(Y, L_0, R_0, max_iterations=5e4L,lambda_1L=0.1,check_optimal_step=FALSE)
#' # objective does not suffer because all mass is shifted to R
#'  objective(Y,out2$L,out2$R)
#' list(L1=sum(out1$L),R1=sum(out1$R),L2=sum(out2$L),R2=sum(out2$R))
#' sum(out2$L)
#' # with L1 regularization on both sides
#'  out3 <- nmf(Y, L_0, R_0, max_iterations=5e4L,lambda_1L=0.1,lambda_1R=0.1,check_optimal_step=FALSE)
#' # with L1 regularization on both sides, raw objective suffers
#'  objective(Y,out3$L,out3$R)
#' list(L1=sum(out1$L),R1=sum(out1$R),L3=sum(out3$L),R3=sum(out3$R))
#'
#' @author Steven E. Pav \email{shabbychef@@gmail.com}
#' @export
nmf <- function(Y, L, R, 
								W_0R=NULL, W_0C=NULL, 
								lambda_1L=0, lambda_1R=0, 
								lambda_2L=0, lambda_2R=0, 
								tau=0.5, annealing_rate=0.10, 
								check_optimal_step=TRUE, 
								zero_tolerance=1e-12,
								max_iterations=1e3L, 
								min_xstep=1e-9,
								verbosity=0) {
	stopifnot(all(Y >= 0))
	stopifnot(all(L >= 0))
	stopifnot(all(R >= 0))
	stopifnot(missing(W_0R) || is.null(W_0R) || all(W_0R >= 0))
	stopifnot(missing(W_0C) || is.null(W_0C) || all(W_0C >= 0))
	stopifnot((0 < tau) && (tau < 1))
	stopifnot((0 <= annealing_rate) && (annealing_rate < 1))

	tau_k <- tau
	finished <- FALSE
	k <- 0
	while (!finished) {
		tau_k <- (1-annealing_rate) * tau_k + annealing_rate
		# update L
		WRt <- (W_0C %**% t(R))
		RWR <- R %*% WRt
		D <- lambda_1L - (W_0R %**% Y) %*% WRt
		F <- (W_0R %**% L) %*% RWR + lambda_2L * L
		gradfL_k <- D + F
		H_kp1 <- pick_direction(L, gradfL_k, F)
		if (check_optimal_step) {
			K_kp1 <- (W_0R %**% H_kp1) %*% RWR + lambda_2L * H_kp1
		} else {
			K_kp1 <- NULL
		}
		LList <- gipm_step(L, gradfL_k, H_kp1, K_kp1, tau_k, k, verbosity)
		L <- LList[[1]]
		Lstep <- LList[[2]]
		# update R
		LtW <- t(L) %**% W_0R
		LWL <- LtW %*% L
		D <- lambda_1R - LtW %*% (Y %**% W_0C)
		F <- LWL %*% (R %**% W_0C) + lambda_2R * R
		gradfR_k <- D + F
		H_kp1 <- pick_direction(R, gradfR_k, F)
		if (check_optimal_step) {
			K_kp1 <- LWL %*% (H_kp1 %**% W_0C) + lambda_2R * H_kp1
		} else {
			K_kp1 <- NULL
		}
		RList <- gipm_step(R, gradfR_k, H_kp1, K_kp1, tau_k, k, verbosity)
		R <- RList[[1]]
		Rstep <- RList[[2]]
		R[R <= zero_tolerance] <- 0
		L[L <= zero_tolerance] <- 0
		k <- k + 1
		converged <- (max(c(Lstep,Rstep)) < min_xstep)
		finished <- converged || (k >= max_iterations) 
	}
	if (verbosity > 1) { 
		print(paste0("terminated after ",k," iterations. converged: ",converged))
	}
	rownames(L) <- rownames(Y)
	colnames(R) <- colnames(Y)
	return(list(L=L,R=R,iterations=k,converged=converged,Lstep=Lstep,Rstep=Rstep))
}

#' @title gnmf .
#'
#' @description 
#'
#' Non-negative matrix factorization with regularization, general form.
#'
#' @details
#'
#' Attempts to factor given non-negative matrix \eqn{Y} as the product \eqn{LR}
#' of two non-negative matrices. The objective function is Frobenius norm
#' with \eqn{\ell_1} and \eqn{\ell_2} regularization terms.
#' We seek to minimize the objective
#' \deqn{\frac{1}{2}tr((Y-LR)' W_{0R} (Y-LR) W_{0C}) + tr(W_{1L}'L} + \tr{W_{1R}'R} + \frac{1}{2} \sum_j \tr(L'W_{2RLj}LW_{2CLj}) + \tr{R'W_{2RRj}RW_{2CRj}),}
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
#' Can also be a list, in which case \coode{W_2CL} must be a list of the same
#' length. The list should consist of \eqn{\ell_2} row penalty matrices.
#' Defaults to all-zeroes matrix, which is no penalty term.
#' @param W_2CL  the \eqn{\ell_2} column penalty matrix for the matrix \eqn{L}.
#' If a scalar, corresponds to that scalar times the identity matrix.
#' Can also be a list, in which case \coode{W_2RL} must be a list of the same
#' length. The list should consist of \eqn{\ell_2} column penalty matrices.
#' Defaults to all-zeroes matrix, which is no penalty term.
#' @param W_2RR  the \eqn{\ell_2} row penalty matrix for the matrix \eqn{R}.
#' If a scalar, corresponds to that scalar times the identity matrix.
#' Can also be a list, in which case \coode{W_2CR} must be a list of the same
#' length. The list should consist of \eqn{\ell_2} row penalty matrices.
#' Defaults to all-zeroes matrix, which is no penalty term.
#' @param W_2CR  the \eqn{\ell_2} column penalty matrix for the matrix \eqn{R}.
#' If a scalar, corresponds to that scalar times the identity matrix.
#' Can also be a list, in which case \coode{W_2RR} must be a list of the same
#' length. The list should consist of \eqn{\ell_2} column penalty matrices.
#' Defaults to all-zeroes matrix, which is no penalty term.
#' @inheritParams nmf
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
#' @seealso \code{\link{nmf}}
#'
#'
#' @author Steven E. Pav \email{shabbychef@@gmail.com}
#' @export
gnmf <- function(Y, L, R, 
								W_0R=NULL, W_0C=NULL, 
								W_1L=0, W_1R=0,
								W_2RL=0, W_2CL=0,
								W_2RR=0, W_2CR=0,
								tau=0.5, annealing_rate=0.10, 
								check_optimal_step=TRUE, 
								zero_tolerance=1e-12,
								max_iterations=1e3L, 
								min_xstep=1e-9,
								verbosity=0) {
	stopifnot(all(Y >= 0))
	stopifnot(all(L >= 0))
	stopifnot(all(R >= 0))
	stopifnot(missing(W_0R) || is.null(W_0R) || all(W_0R >= 0))
	stopifnot(missing(W_0C) || is.null(W_0C) || all(W_0C >= 0))
	stopifnot(missing(W_1R) || is.null(W_1R) || all(W_1R >= 0))
	stopifnot(missing(W_1L) || is.null(W_1L) || all(W_1L >= 0))
	stopifnot((0 < tau) && (tau < 1))
	stopifnot((0 <= annealing_rate) && (annealing_rate < 1))

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


	tau_k <- tau
	finished <- FALSE
	k <- 0
	while (!finished) {
		tau_k <- (1-annealing_rate) * tau_k + annealing_rate

		# update L
		WRt <- (W_0C %**% t(R))
		RWR <- R %*% WRt
		D <- W_1L - (W_0R %**% Y) %*% WRt
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
		LList <- gipm_step(L, gradfL_k, H_kp1, K_kp1, tau_k, k, verbosity)
		L <- LList[[1]]
		Lstep <- LList[[2]]

		# update R
		LtW <- t(L) %**% W_0R
		LWL <- LtW %*% L
		D <- W_1R - LtW %*% (Y %**% W_0C)
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
		RList <- gipm_step(R, gradfR_k, H_kp1, K_kp1, tau_k, k, verbosity)
		R <- RList[[1]]
		Rstep <- RList[[2]]
		R[R <= zero_tolerance] <- 0
		L[L <= zero_tolerance] <- 0
		k <- k + 1
		converged <- (max(c(Lstep,Rstep)) < min_xstep)
		finished <- converged || (k >= max_iterations) 
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
