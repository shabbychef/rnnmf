# Copyright 2024-2024 Steven E. Pav. All Rights Reserved.
# Author: Steven E. Pav

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

# env var:
# nb: 
# see also:
# todo:
# changelog: 
#
# Created: 2024.08.20
# Copyright: Steven E. Pav, 2024-2024
# Author: Steven E. Pav
# Comments: Steven E. Pav


# helpers
# random non-negative matrix; equals zero with some probability.
randmat <- function(nr,nc,min=-0.25,max=1) { matrix(pmax(0,runif(nr*nc,min=min,max=max)),nrow=nr) }

# just test if everything runs...

quadratic_objective <- function(Y, L, R) { sum((Y - L %*% R)^2) }

context("test giqpm")#FOLDUP
test_that("giqpm runs",{#FOLDUP
	set.seed(1234)
	nr <- 100
	nc <- 20
	LL <- randmat(nr,nc)
	Gmat <- t(LL) %*% LL
	dvec <- -runif(nc)
	expect_error(out0 <- giqpm(Gmat, dvec), NA)

	preG <- randmat(nr,nr+nc)
	G <- preG %*% t(preG)
	d <- - runif(nr)
	expect_error(y1 <- giqpm(G, d),NA)
	objective <- function(G, d, x) { as.numeric(0.5 * t(x) %*% (G %*% x) + t(x) %*% d) }

	# this does not converge to an actual solution!
	steepest_step_func <- function(gradf, ...) { return(-gradf) }
	expect_error(y2 <- giqpm(G, d, step_func = steepest_step_func),NA)

	scaled_step_func <- function(gradf, Gx, Gmat, dvec, x0, ...) { return(-gradf * abs(x0)) }
	expect_error(y3 <- giqpm(G, d, step_func = scaled_step_func),NA)

	sqrt_step_func <- function(gradf, Gx, Gmat, dvec, x0, ...) { return(-gradf * abs(sqrt(x0))) }
	expect_error(y4 <- giqpm(G, d, step_func = sqrt_step_func),NA)

	complementarity_stepfunc <- function(gradf, Gx, Gmat, dvec, x0, ...) { return(-gradf * x0) }
	expect_error(y5 <- giqpm(G, d, step_func = complementarity_stepfunc),NA)

	expect_lt(objective(G, d, y4$x), objective(G, d, y2$x))
	expect_lt(objective(G, d, y1$x), objective(G, d, y4$x))
	expect_lt(objective(G, d, y3$x), objective(G, d, y4$x))
	expect_lt(objective(G, d, y5$x), objective(G, d, y4$x))
})#UNFOLD
#UNFOLD
context("test aurnmf")#FOLDUP
test_that("aurnmf runs",{#FOLDUP
	nr <- 20
	nc <- 6
	dm <- 2

	set.seed(1234)
	real_L <- randmat(nr,dm)
	real_R <- randmat(dm,nc)
	Y <- real_L %*% real_R
	L_0 <- randmat(nr,dm)
	R_0 <- randmat(dm,nc)
	W_0R <- diag(runif(nr))
	W_0C <- diag(runif(nc))
	max_iter <- 20
	for (check_optimal_step in c(TRUE, FALSE)) {
		# without regularization
		expect_error(out1 <- aurnmf(Y, L_0, R_0, max_iterations=max_iter,check_optimal_step=check_optimal_step),NA)

		# with L1 regularization on one side
		expect_error(out2 <- aurnmf(Y, L_0, R_0, max_iterations=max_iter,lambda_1L=0.1,check_optimal_step=check_optimal_step),NA)

		# with L1 regularization on both sides
		expect_error(out3 <- aurnmf(Y, L_0, R_0, max_iterations=max_iter,lambda_1L=0.1,lambda_1R=0.1,check_optimal_step=check_optimal_step),NA)

		# with non-orthogonality penalty
		expect_error(out4 <- aurnmf(Y, L_0, R_0, max_iterations=max_iter,lambda_1L=0.1,lambda_1R=0.1,gamma_2L=0.1,gamma_2R=0.1,check_optimal_step=check_optimal_step),NA)

		# with weighting
		expect_error(out5 <- aurnmf(Y, L_0, R_0, W_0R=W_0R, W_0C=W_0C, max_iterations=max_iter, lambda_1L=0.1,lambda_1R=0.1,gamma_2L=0.1,gamma_2R=0.1,check_optimal_step=check_optimal_step),NA)
	}

})#UNFOLD
test_that("aurnmf callbacks",{#FOLDUP
	nr <- 100
	nc <- 20
	dm <- 4

	set.seed(1234)
	real_L <- randmat(nr,dm)
	real_R <- randmat(dm,nc)
	Y <- real_L %*% real_R

	max_iterations <- 1e2L
	it_history <<- rep(NA_real_, max_iterations)
	on_iteration_end <- function(iteration, Y, L, R, ...) {
		it_history[iteration] <<- quadratic_objective(Y,L,R)
	}
	L_0 <- randmat(nr,dm)
	R_0 <- randmat(dm,nc)
	expect_error(out1 <- aurnmf(Y, L_0, R_0, max_iterations=max_iterations, on_iteration_end=on_iteration_end),NA)

})#UNFOLD
test_that("aurnmf is good",{#FOLDUP
	nr <- 100
	nc <- 10
	dm <- 2

	set.seed(5678)
	real_L <- randmat(nr,dm)
	real_R <- randmat(dm,nc)
	Y <- real_L %*% real_R
	L_0 <- randmat(nr,dm)
	R_0 <- randmat(dm,nc)
	for (check_optimal_step in c(TRUE, FALSE)) {
		# without regularization
		expect_error(out1 <- aurnmf(Y, L_0, R_0, max_iterations=1e2L,check_optimal_step=check_optimal_step),NA)
		expect_lt(quadratic_objective(Y, out1$L, out1$R), quadratic_objective(Y, L_0, R_0))
	}
})#UNFOLD
test_that("aurnmf runs sparse input",{#FOLDUP
	nr <- 20
	nc <- 6
	dm <- 2

	set.seed(1234)
	real_L <- randmat(nr,dm,min=-1)
	real_R <- randmat(dm,nc,min=-1)
	Y <- as(real_L %*% real_R,"sparseMatrix")
	L_0 <- as(randmat(nr,dm,min=-1), "sparseMatrix")
	R_0 <- as(randmat(dm,nc,min=-1), "sparseMatrix")
	W_0R <- as(diag(runif(nr)), "sparseMatrix")
	W_0C <- as(diag(runif(nc)), "sparseMatrix")
	max_iter <- 20
	for (check_optimal_step in c(TRUE, FALSE)) {
		# without regularization
		expect_error(out1 <- aurnmf(Y, L_0, R_0, max_iterations=max_iter,check_optimal_step=check_optimal_step),NA)

		# with L1 regularization on one side
		expect_error(out2 <- aurnmf(Y, L_0, R_0, max_iterations=max_iter,lambda_1L=0.1,check_optimal_step=check_optimal_step),NA)

		# with L1 regularization on both sides
		expect_error(out3 <- aurnmf(Y, L_0, R_0, max_iterations=max_iter,lambda_1L=0.1,lambda_1R=0.1,check_optimal_step=check_optimal_step),NA)

		# with non-orthogonality penalty
		expect_error(out4 <- aurnmf(Y, L_0, R_0, max_iterations=max_iter,lambda_1L=0.1,lambda_1R=0.1,gamma_2L=0.1,gamma_2R=0.1,check_optimal_step=check_optimal_step),NA)

		# with weighting
		expect_error(out5 <- aurnmf(Y, L_0, R_0, W_0R=W_0R, W_0C=W_0C, max_iterations=max_iter, lambda_1L=0.1,lambda_1R=0.1,gamma_2L=0.1,gamma_2R=0.1,check_optimal_step=check_optimal_step),NA)
	}

})#UNFOLD
#UNFOLD
context("test murnmf")#FOLDUP
test_that("murnmf runs",{#FOLDUP
	nr <- 100
	nc <- 20
	dm <- 4

	set.seed(1234)
	real_L <- randmat(nr,dm)
	real_R <- randmat(dm,nc)
	Y <- real_L %*% real_R
	# without regularization
	objective <- function(Y, L, R) { sum((Y - L %*% R)^2) }

	max_iter <- 100
	L_0 <- randmat(nr,dm)
	R_0 <- randmat(dm,nc)
	expect_error(out1 <- murnmf(Y, L_0, R_0, max_iterations=max_iter),NA)

	# with L1 regularization on one side
	expect_error(out2 <- murnmf(Y, L_0, R_0, max_iterations=max_iter,lambda_1L=0.1),NA)

	# with L1 regularization on both sides
	expect_error(out3 <- murnmf(Y, L_0, R_0, max_iterations=max_iter,lambda_1L=0.1,lambda_1R=0.1),NA)

	# with non-orthogonality penalty
	expect_error(out4 <- murnmf(Y, L_0, R_0, max_iterations=max_iter,lambda_1L=0.1,lambda_1R=0.1,gamma_2L=0.1,gamma_2R=0.1),NA)

})#UNFOLD
test_that("murnmf is good",{#FOLDUP
	nr <- 100
	nc <- 10
	dm <- 2

	set.seed(5678)
	real_L <- randmat(nr,dm)
	real_R <- randmat(dm,nc)
	Y <- real_L %*% real_R
	L_0 <- randmat(nr,dm)
	R_0 <- randmat(dm,nc)
	# without regularization
	out1 <- murnmf(Y, L_0, R_0, max_iterations=1e2L)
	expect_error(out1 <- murnmf(Y, L_0, R_0, max_iterations=1e2L),NA)
	expect_lt(quadratic_objective(Y, out1$L, out1$R), quadratic_objective(Y, L_0, R_0))
})#UNFOLD
#UNFOLD
context("test gaurnmf")#FOLDUP
test_that("gaurnmf runs",{#FOLDUP
	nr <- 20
	nc <- 6
	dm <- 2

	set.seed(1234)
	real_L <- randmat(nr,dm)
	real_R <- randmat(dm,nc)
	Y <- real_L %*% real_R
	L_0 <- randmat(nr,dm)
	R_0 <- randmat(dm,nc)
	# without regularization
	expect_error(out1 <- gaurnmf(Y, L_0, R_0, max_iterations=1e2L,check_optimal_step=FALSE),NA)

	# with L1 regularizations
	W_1L <- randmat(nrow(L_0), ncol(L_0))
	W_1R <- randmat(nrow(R_0), ncol(R_0))
	expect_error(out2 <- gaurnmf(Y, L_0, R_0, W_1L=W_1L, W_1R=W_1R, max_iterations=5e1L,check_optimal_step=FALSE),NA)
	expect_error(out2 <- gaurnmf(Y, L_0, R_0, W_1L=W_1L, W_1R=0, max_iterations=5e1L,check_optimal_step=FALSE),NA)
	expect_error(out2 <- gaurnmf(Y, L_0, R_0, W_1L=0, W_1R=W_1R, max_iterations=5e1L,check_optimal_step=FALSE),NA)

	# with L2 regularizations
	W_2RL <- randmat(nrow(L_0), nrow(L_0))
	W_2CL <- randmat(ncol(L_0), ncol(L_0))
	W_2RR <- randmat(nrow(R_0), nrow(R_0))
	W_2CR <- randmat(ncol(R_0), ncol(R_0))
	expect_error(out3 <- gaurnmf(Y, L_0, R_0, W_1L=0, W_1R=W_1R, 
														W_2RL=W_2RL,W_2CL=W_2CL,W_2RR=W_2RR,W_2CR=W_2CR,
														max_iterations=5e1L,check_optimal_step=FALSE),NA)
	expect_error(out3 <- gaurnmf(Y, L_0, R_0, W_1L=0, W_1R=W_1R, 
														W_2RL=list(W_2RL),W_2CL=list(W_2CL),W_2RR=list(W_2RR),W_2CR=list(W_2CR),
														max_iterations=5e1L,check_optimal_step=FALSE),NA)

	# with a list of L2 regularizations
	W_2RL1 <- randmat(nrow(L_0), nrow(L_0))
	W_2CL1 <- randmat(ncol(L_0), ncol(L_0))
	W_2RL2 <- randmat(nrow(L_0), nrow(L_0))
	W_2CL2 <- randmat(ncol(L_0), ncol(L_0))
	expect_error(out4 <- gaurnmf(Y, L_0, R_0, W_1L=0, W_1R=W_1R, 
														W_2RL=list(W_2RL1,W_2RL2),W_2CL=list(W_2CL1,W_2CL2),W_2RR=list(W_2RR),W_2CR=list(W_2CR),
														max_iterations=5e1L,check_optimal_step=FALSE),NA)
	expect_error(out4 <- gaurnmf(Y, L_0, R_0, W_1L=0, W_1R=W_1R, 
														W_2RL=list(W_2RL1,W_2RL2),W_2CL=list(W_2CL1,W_2CL2),W_2RR=list(W_2RR,0.2),W_2CR=list(W_2CR,0.2),
														max_iterations=5e1L,check_optimal_step=FALSE),NA)

})#UNFOLD
test_that("gaurnmf is good",{#FOLDUP
	nr <- 100
	nc <- 10
	dm <- 2

	set.seed(5678)
	real_L <- randmat(nr,dm)
	real_R <- randmat(dm,nc)
	Y <- real_L %*% real_R
	L_0 <- randmat(nr,dm)
	R_0 <- randmat(dm,nc)
	for (check_optimal_step in c(TRUE, FALSE)) {
		# without regularization
		expect_error(out1 <- gaurnmf(Y, L_0, R_0, max_iterations=5e1L,check_optimal_step=check_optimal_step),NA)
		expect_lt(quadratic_objective(Y, out1$L, out1$R), quadratic_objective(Y, L_0, R_0))
	}
})#UNFOLD
test_that("gaurnmf runs sparse input",{#FOLDUP
	nr <- 20
	nc <- 6
	dm <- 2

	set.seed(1234)
	real_L <- randmat(nr,dm,min=-1)
	real_R <- randmat(dm,nc,min=-1)
	Y <- as(real_L %*% real_R,"sparseMatrix")
	L_0 <- as(randmat(nr,dm,min=-0.5),"sparseMatrix")
	R_0 <- as(randmat(dm,nc,min=-0.5),"sparseMatrix")
	# without regularization
	expect_error(out1 <- gaurnmf(Y, L_0, R_0, max_iterations=5e1L,check_optimal_step=FALSE),NA)

	# with L1 regularizations
	W_1L <- randmat(nrow(L_0), ncol(L_0))
	W_1R <- randmat(nrow(R_0), ncol(R_0))
	expect_error(out2 <- gaurnmf(Y, L_0, R_0, W_1L=W_1L, W_1R=W_1R, max_iterations=5e1L,check_optimal_step=FALSE),NA)
	expect_error(out2 <- gaurnmf(Y, L_0, R_0, W_1L=W_1L, W_1R=0, max_iterations=5e1L,check_optimal_step=FALSE),NA)
	expect_error(out2 <- gaurnmf(Y, L_0, R_0, W_1L=0, W_1R=W_1R, max_iterations=5e1L,check_optimal_step=FALSE),NA)

	# with L2 regularizations
	W_2RL <- randmat(nrow(L_0), nrow(L_0))
	W_2CL <- randmat(ncol(L_0), ncol(L_0))
	W_2RR <- randmat(nrow(R_0), nrow(R_0))
	W_2CR <- randmat(ncol(R_0), ncol(R_0))
	expect_error(out3 <- gaurnmf(Y, L_0, R_0, W_1L=0, W_1R=W_1R, 
														W_2RL=W_2RL,W_2CL=W_2CL,W_2RR=W_2RR,W_2CR=W_2CR,
														max_iterations=5e1L,check_optimal_step=FALSE),NA)
	expect_error(out3 <- gaurnmf(Y, L_0, R_0, W_1L=0, W_1R=W_1R, 
														W_2RL=list(W_2RL),W_2CL=list(W_2CL),W_2RR=list(W_2RR),W_2CR=list(W_2CR),
														max_iterations=5e1L,check_optimal_step=FALSE),NA)

	# with a list of L2 regularizations
	W_2RL1 <- randmat(nrow(L_0), nrow(L_0))
	W_2CL1 <- randmat(ncol(L_0), ncol(L_0))
	W_2RL2 <- randmat(nrow(L_0), nrow(L_0))
	W_2CL2 <- randmat(ncol(L_0), ncol(L_0))
	expect_error(out4 <- gaurnmf(Y, L_0, R_0, W_1L=0, W_1R=W_1R, 
														W_2RL=list(W_2RL1,W_2RL2),W_2CL=list(W_2CL1,W_2CL2),W_2RR=list(W_2RR),W_2CR=list(W_2CR),
														max_iterations=5e1L,check_optimal_step=FALSE),NA)
	expect_error(out4 <- gaurnmf(Y, L_0, R_0, W_1L=0, W_1R=W_1R, 
														W_2RL=list(W_2RL1,W_2RL2),W_2CL=list(W_2CL1,W_2CL2),W_2RR=list(W_2RR,0.2),W_2CR=list(W_2CR,0.2),
														max_iterations=5e1L,check_optimal_step=FALSE),NA)

})#UNFOLD
#UNFOLD

#for vim modeline: (do not edit)
# vim:ts=2:sw=2:tw=79:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
