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
# Created: 2024.10.14
# Copyright: Steven E. Pav, 2024
# Author: Steven E. Pav <steven@gilgamath.com>
# Comments: Steven E. Pav

#' @title murnmf .
#'
#' @description 
#'
#' Multiplicative update Non-negative matrix factorization with regularization.
#'
#' @details
#'
#' This function uses multiplicative updates only, and may not optimize the
#' nominal objective. It is also unlikely to achieve optimality.
#' This code is for reference purposes and is not suited for usage other
#' than research and experimentation.
#'
#' @param epsilon  the numerator clipping value.
#' @inherit aurnmf return params details references 
#' @keywords optimization
#' @template etc
#' @template poc
#' @template ref-merritt
#' @template ref-pav
#' @template ref-leeseung
#' @seealso \code{\link{aurnmf}}, \code{\link{gaurnmf}}
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
#'  out1 <- murnmf(Y, L_0, R_0, max_iterations=5e3L)
#'  objective(Y,out1$L,out1$R)
#' # with L1 regularization on one side
#'  out2 <- murnmf(Y, L_0, R_0, max_iterations=5e3L,lambda_1L=0.1)
#' # objective does not suffer because all mass is shifted to R
#'  objective(Y,out2$L,out2$R)
#' list(L1=sum(out1$L),R1=sum(out1$R),L2=sum(out2$L),R2=sum(out2$R))
#' sum(out2$L)
#' # with L1 regularization on both sides
#'  out3 <- murnmf(Y, L_0, R_0, max_iterations=5e3L,lambda_1L=0.1,lambda_1R=0.1)
#' # with L1 regularization on both sides, raw objective suffers
#'  objective(Y,out3$L,out3$R)
#' list(L1=sum(out1$L),R1=sum(out1$R),L3=sum(out3$L),R3=sum(out3$R))
#'
#' \donttest{
#' # example showing how to use the on_iteration_end callback to save iterates.
#' max_iterations <- 1e3L
#' it_history <<- rep(NA_real_, max_iterations)
#' quadratic_objective <- function(Y, L, R) { sum((Y - L %*% R)^2) }
#' on_iteration_end <- function(iteration, Y, L, R, ...) {
#'   it_history[iteration] <<- quadratic_objective(Y,L,R)
#' }
#' out1b <- murnmf(Y, L_0, R_0, max_iterations=max_iterations, on_iteration_end=on_iteration_end)
#' }
#'
#' # should work on sparse matrices too, but beware zeros in the initial estimates
#' if (require(Matrix)) { 
#'  real_L <- randmat(nr,dm,min=-1)
#'  real_R <- randmat(dm,nc,min=-1)
#'  Y <- as(real_L %*% real_R, "sparseMatrix")
#'  L_0 <- randmat(nr,dm)
#'  R_0 <- randmat(dm,nc)
#'  out1 <- murnmf(Y, L_0, R_0, max_iterations=1e2L)
#' }
#'
#' @author Steven E. Pav \email{shabbychef@@gmail.com}
#' @export
murnmf <- function(Y, L, R, 
								W_0R=NULL, W_0C=NULL, 
								lambda_1L=0, lambda_1R=0, 
								lambda_2L=0, lambda_2R=0, 
								gamma_2L=0,  gamma_2R=0,
								epsilon=1e-7,
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

	fixd <- fix_LR_names(L, R)
	L <- fixd$L
	R <- fixd$R

	# precompute
	W_0R_Y <- (W_0R %**% Y)
	Y_W_0C <- (Y %**% W_0C)

	finished <- FALSE
	k <- 0
	while (!finished) {
		# update L
		WRt <- (W_0C %**% t(R))
		RWR <- R %*% WRt
		H <- lambda_1L - W_0R_Y %*% WRt
		H <- pmin(H, -epsilon)
		F <- (W_0R %**% L) %*% RWR + lambda_2L * L
		if (gamma_2L > 0) {
			F <- F + gamma_2L * times_orth(L)
		}
		Lprev <- L
		Fok <- (F > 0) & ! is.nan(F)
		L[Fok] <- -L[Fok] * (H[Fok] / F[Fok])
		Lstep <- max(abs(L - Lprev))
		# update R
		LtW <- t(L) %**% W_0R
		LWL <- LtW %*% L
		H <- lambda_1R - LtW %*% Y_W_0C
		H <- pmin(H, -epsilon)
		F <- LWL %*% (R %**% W_0C) + lambda_2R * R
		if (gamma_2R > 0) {
			F <- F + gamma_2R * times_orth(R)
		}
		Rprev <- R
		Fok <- (F > 0) & ! is.nan(F)
		R[Fok] <- -R[Fok] * (H[Fok] / F[Fok])
		Rstep <- max(abs(R - Rprev))
		
		k <- k + 1
		converged <- (max(c(Lstep,Rstep)) < min_xstep)
		finished <- converged || (k >= max_iterations) || all(!Fok)
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
