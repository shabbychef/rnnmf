######################
# 
# Created: 2024-08-17
# Copyright: Steven E. Pav, 2024
# Author: Steven E. Pav
######################

############### FLAGS ###############

VMAJOR 						 = 0
VMINOR 						 = 1
VPATCH  					 = 0
VDEV 							 = .0001
PKG_NAME 					:= rnmf

RPKG_USES_RCPP 		:= 0

include ./rpkg_make/Makefile

#for vim modeline: (do not edit)
# vim:ts=2:sw=2:tw=129:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:tags=.tags;:syn=make:ft=make:ai:si:cin:nu:fo=croqt:cino=p0t0c5(0:
