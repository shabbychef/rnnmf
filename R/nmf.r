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
	Fzer <- F <= 0
	Xzer <- X_k <= 0
	F[Fzer] <- 1
	H_kp1 <- - gradfX_k * X_k / F
	H_kp1[Fzer & Xzer] <- pmax(-gradfX_k[Fzer & Xzer], 0)
	return(H_kp1)
}

# like matrix multiplication, but if X or Y is not given or is null, treat as the identity.
`%**%` <- function(X, Y) {
	if (missing(X) || is.null(X)) {
		return(Y)
	} else if (missing(Y) || is.null(Y)) {
		return(X)
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
#' One sentence or so that tells you some more.
#'
#' @details
#'
#' Really detailed. \eqn{\zeta}{zeta}.
#'
#' A list:
#' \itemize{
#' \item I use \eqn{n}{n} to stand for blah.
#' \item and so forth....
#' }
#' \describe{
#' \item{a}{value.}
#' \item{b}{factor.}
#' }
#'
#'
#' @param Y  an \eqn{r \times c} matrix to be decomposed.
#' Should have non-negative elements, though we do not check.
#' @param L  an \eqn{r \times d} matrix of the initial estimate of L.
#' Should have non-negative elements, though we do not check.
#' @param R  an \eqn{d \times c} matrix of the initial estimate of R.
#' Should have non-negative elements, though we do not check.
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
#' @template ref-merritt
#'
#' @examples 
#'
#' @author Steven E. Pav \email{shabbychef@@gmail.com}
#' @export
nmf <- function(Y, L, R, 
								W_0R=NULL, W_0C=NULL, 
								lambda_1L=0, lambda_1R=0, 
								lambda_2L=0, lambda_2R=0, 
								tau=0.5, annealing_rate=0.25, 
								check_optimal_step=TRUE, 
								zero_tolerance=1e-9,
								max_iterations=1e3L, 
								min_xstep=1e-9,
								verbosity=0) {
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


#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
