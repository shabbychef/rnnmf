# Copyright 2024-2024 Steven E. Pav. All Rights Reserved.
# Author: Steven E. Pav

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
# random non-negative matrix
randmat <- function(nr,nc,zero_p=0.2) { matrix(pmax(0,runif(nr*nc)-zero_p),nrow=nr) }

# just test if everything runs...

context("test nmf")#FOLDUP
test_that("nmf runs",{#FOLDUP
	nr <- 100
	nc <- 20
	dm <- 4

	set.seed(1234)
	real_L <- randmat(nr,dm)
	real_R <- randmat(dm,nc)
	Y <- real_L %*% real_R
	# without regularization
	objective <- function(Y, L, R) { sum((Y - L %*% R)^2) }

	L_0 <- randmat(nr,dm)
	R_0 <- randmat(dm,nc)
	expect_error(out1 <- nmf(Y, L_0, R_0, max_iterations=5e3L,check_optimal_step=FALSE),NA)

	# with L1 regularization on one side
	expect_error(out2 <- nmf(Y, L_0, R_0, max_iterations=5e3L,lambda_1L=0.1,check_optimal_step=FALSE),NA)

	# with L1 regularization on both sides
	expect_error(out3 <- nmf(Y, L_0, R_0, max_iterations=5e3L,lambda_1L=0.1,lambda_1R=0.1,check_optimal_step=FALSE),NA)
})#UNFOLD
#UNFOLD

#for vim modeline: (do not edit)
# vim:ts=2:sw=2:tw=79:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
