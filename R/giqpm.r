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
# Author: Steven E. Pav <shabbychef@gmail.com>
# Comments: Steven E. Pav

# you could imagine "multiplication" instead as doing some joins and
# summarizes or whatever
default_multfunc <- function(Gmat, x0) {
	return(Gmat %*% x0)
}

default_gradfunc <- function(Gmat, dvec, x0, Gx) {
	return(Gx + dvec)
}
# steepest descent
# note that the problem with steepest_descent is that when you have a point on
# the boundary, it may point towards infeasible solutions and then the chosen
# steplength will shrink to zero, and you will get no progress, which stinks.
steepest_stepfunc <- function(gradf, ...) {
	return(-gradf)
}
# lee and seung step function
lee_seung_stepfunc <- function(gradf, Gx, Gmat, dvec, x0, ...) {
	bad_deno <- (Gx <= 0) | is.na(Gx)
	Gx[bad_deno] <- 1
	# check Gx for zeroes?
	return(-gradf * x0 / Gx)
}
# my weird idea
# this also gives crap results.
weird_stepfunc <- function(gradf, Gx, Gmat, dvec, x0, multfunc) {
	b <- pmax(0, x0 + dvec/diag(Gmat))
	Gb <- multfunc(Gmat, b)
	# check for Gb = 0?
	raty <- b
	raty[Gb == 0] <- 0
	raty[Gb != 0] <- b[Gb != 0] / Gb[Gb != 0]
	return(-gradf * raty)
}
complementarity_stepfunc <- function(gradf, Gx, Gmat, dvec, x0, ...) {
	return(-gradf * x0)
}
default_stepfunc <- lee_seung_stepfunc


#' @title giqpm .
#'
#' @description 
#'
#' Generalized Iterative Quadratic Programming Method for non-negative quadratic optimization.
#'
#' @details
#'
#' Iteratively solves the problem
#' \deqn{\min_x \frac{1}{2}x^{\top}G x + d^{\top}x}{min x  0.5 x'Gx + d'x}
#' subject to the elementwise constraint \eqn{x \ge 0}{x >= 0}.
#'
#' This implementation allows the user to specify methods to perform matrix by
#' vector multiplication, computation of the gradient (which should be
#' \eqn{G x + d}{Gx + d}), and computation of the step direction.
#' By default we compute the optimal step in the given step direction.
#'
#' @param Gmat  a representation of the matrix \eqn{G}.
#' @param dvec  a representation of the vector \eqn{d}.
#' @param x0    the initial iterate. If none given, we spawn one of the same
#' size as \code{dvec}.
#' @param tau   the starting shrinkage factor applied to the step length.
#' Should be a value in \eqn{(0,1)}.
#' @param annealing_rate  the rate at which we scale the shrinkage factor towards 1.
#' Should be a value in \eqn{[0,1)}.
#' @param check_optimal_step  if TRUE, we attempt to take the optimal step
#' length in the given direction. If not, we merely take the longest feasible
#' step in the step direction.
#' @param mult_func  a function which takes matrix and vector and performs
#' matrix multiplication. 
#' The default does this on matrix and vector input,
#' but the user can implement this for some implicit versions of the problem.
#' @param grad_func  a function which takes matrix \eqn{G}, vector \eqn{d}, 
#' the current iterate \eqn{x} and the product \eqn{Gx} and is supposed to
#' compute \eqn{Gx + d}.
#' The default does this on matrix and vector input,
#' but the user can implement this for some implicit versions of the problem.
#' @param step_func  a function which takes the vector gradient, the product 
#' \eqn{Gx}, the matrix \eqn{G}, vector \eqn{d}, vector \eqn{x} and the
#' \code{mult_func} and produces a step vector.
#' By default this step vector is the Lee-Seung step vector, namely
#' \eqn{-(Gx + d) * x / d}, with Hadamard product and division.
#' @param zero_tolerance  values of \eqn{x} less than this will be \sQuote{snapped} to zero.
#' This happens at the end of the iteration and does not affect the measurement
#' of convergence.
#' @param max_iterations  the maximum number of iterations to perform.
#' @param min_xstep   the minimum L-infinity norm of the step taken.
#' Once the step falls under this value, we terminate.
#' @param verbosity   controls whether we print information to the console.
#' @return a list with the elements
#' \describe{
#' \item{x}{The final iterate.}
#' \item{iterations}{The number of iterations taken.}
#' \item{converged}{Whether convergence was detected.}
#' }
#' @keywords optimization
#' @template etc
#' @template poc
#' @template ref-pav
#' @template ref-merritt
#'
#' @examples 
#'
#' set.seed(1234)
#' ssiz <- 100
#' preG <- matrix(runif(ssiz*(ssiz+20)),nrow=ssiz)
#' G <- preG %*% t(preG)
#' d <- - runif(ssiz)
#' y1 <- giqpm(G, d)
#' objective <- function(G, d, x) { as.numeric(0.5 * t(x) %*% (G %*% x) + t(x) %*% d) }
#'
#' # this does not converge to an actual solution!
#' steepest_step_func <- function(gradf, ...) { return(-gradf) }
#' y2 <- giqpm(G, d, step_func = steepest_step_func)
#'
#' scaled_step_func <- function(gradf, Gx, Gmat, dvec, x0, ...) { return(-gradf * abs(x0)) }
#' y3 <- giqpm(G, d, step_func = scaled_step_func)
#'
#' sqrt_step_func <- function(gradf, Gx, Gmat, dvec, x0, ...) { return(-gradf * abs(sqrt(x0))) }
#' y4 <- giqpm(G, d, step_func = sqrt_step_func)
#' 
#' complementarity_stepfunc <- function(gradf, Gx, Gmat, dvec, x0, ...) { return(-gradf * x0) }
#' y5 <- giqpm(G, d, step_func = complementarity_stepfunc)
#'
#' objective(G, d, y1$x)
#' objective(G, d, y2$x)
#' objective(G, d, y3$x)
#' objective(G, d, y4$x)
#' objective(G, d, y5$x)
#'
#' @importFrom stats runif
#' @author Steven E. Pav \email{shabbychef@@gmail.com}
#' @export
giqpm <- function(Gmat, dvec, x0=NULL, 
								 tau=0.5, annealing_rate=0.25,
								 check_optimal_step=TRUE,
								 mult_func=NULL,
								 grad_func=NULL,
								 step_func=NULL,
								 zero_tolerance=1e-9,
								 max_iterations=1e3L, 
								 min_xstep=1e-9,
								 verbosity=0) {
	if (is.null(x0)) {
		x0 <- runif(length(dvec))
	}
	if (missing(mult_func) || is.null(mult_func)) { mult_func <- default_multfunc }
	if (missing(grad_func) || is.null(grad_func)) { grad_func <- default_gradfunc }
	if (missing(step_func) || is.null(step_func)) { step_func <- default_stepfunc }
	stopifnot((0 < tau) && (tau < 1))
	stopifnot((0 <= annealing_rate) && (annealing_rate < 1))

	tau_k <- tau
	finished <- FALSE
	k <- 0
	while (!finished) {
		Gx <- mult_func(Gmat, x0)
		gradf <- grad_func(Gmat, dvec, x0, Gx)
		hstep <- step_func(gradf, Gx, Gmat, dvec, x0, mult_func)
		# check for positive or zero values in hstep?
		checkus <- hstep < 0
		if (any(checkus)) {
			alpha_hat_k <- min(-x0[checkus] / hstep[checkus])
		} else if (any(hstep != 0)) {
			alpha_hat_k <- Inf
		} else {
			# all hstep are zero. converged.
			converged <- TRUE
			break
		}
		direx <- sum(gradf * hstep)
		if (direx >= 0) {
			warning("step direction misaligned to gradient?")
		}
		if (check_optimal_step) {
			denom <- sum(hstep * mult_func(Gmat, hstep))
			if (denom <= 0) {
				warning("violation of PD in G?")
			}
			alpha_star_k <- -direx / denom
		} else {
			alpha_star_k <- 1
		}

		tau_k <- (1-annealing_rate) * tau_k + annealing_rate
		alpha_k <- min(tau_k * alpha_hat_k, alpha_star_k)
		if (verbosity > 0) { print(paste0("iteration ", k, " taking step of size ",alpha_k,".")) }
		fullstep <- alpha_k * hstep
		x0 <- x0 + fullstep
		x0[x0 <= zero_tolerance] <- 0
		k <- k + 1
		converged <- (max(abs(fullstep)) < min_xstep)
		finished <- converged || (k >= max_iterations) 
	}
	if (verbosity > 1) { 
		print(paste0("terminated after ",k," iterations. converged: ",converged))
	}
	return(list(x=x0,iterations=k,converged=converged))
}

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
