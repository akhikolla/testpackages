# The ISOpureR package is copyright (c) 2014 Ontario Institute for Cancer Research (OICR)
# This package and its accompanying libraries is free software; you can redistribute it and/or modify it under the terms of the GPL
# (either version 1, or at your option, any later version) or the Artistic License 2.0.  Refer to LICENSE for the full license text.
# OICR makes no representations whatsoever as to the SOFTWARE contained herein.  It is experimental in nature and is provided WITHOUT
# WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE OR ANY OTHER WARRANTY, EXPRESS OR IMPLIED. OICR MAKES NO REPRESENTATION
# OR WARRANTY THAT THE USE OF THIS SOFTWARE WILL NOT INFRINGE ANY PATENT OR OTHER PROPRIETARY RIGHT.
# By downloading this SOFTWARE, your Institution hereby indemnifies OICR against any loss, claim, damage or liability, of whatsoever kind or
# nature, which may arise from your Institution's respective use, handling or storage of the SOFTWARE.
# If publications result from research using this SOFTWARE, we ask that the Ontario Institute for Cancer Research be acknowledged and/or
# credit be given to OICR scientists, as scientifically appropriate.

### FUNCTION: ISOpure.calculate.tac.R #############################################################
#
# This function performs the mathematical calculations taking bulk tumor data and deconvolved profiles
# and returning deconvolved tumour adjacent cell profiles.
# 
# Function call: ISOpure.calculate.tac <- function(tumor.profiles, deconvolved.profiles, purity.estimates)

### INPUT #########################################################################################
#  tumor.profiles: a GxD matrix representing gene expression profiles of heterogeneous (mixed) tumor 
#  samples, where G is the number of genes, D is the number of tumor samples.
#
#  deconvolved.profiles: a GxD matrix representing gene expression profiles of purified 
#  (ISOpure output) tumor samples, where G is the number of genes, D is the number of tumor samples.
#
#  purity.estimates: a vector D representing the purity estimates (output from ISOpure)
#  samples, where G is the number of genes, D is the number of tumor samples.

### OUTPUT ########################################################################################
#
# a GxD matrix representing gene expression profiles of purified (ISOpure output) tumor adjacent 
# cell signal, where G is the number of genes, D is the number of tumor samples.

ISOpure.calculate.tac <- function(tumor.profiles, deconvolved.profiles, purity.estimates) {

	# check that the dimensions match
	if (nrow(tumor.profiles) != nrow(deconvolved.profiles)) {
		stop('There are different numbers of rows in tumour.profiles and deconvolved.profiles -- make sure the rows represent the same genes');
	}
	if (ncol(tumor.profiles) != ncol(deconvolved.profiles)) {
		stop('There are different numbers of columns in tumour.profiles and deconvolved.profiles -- make sure the rows represent the same samples');
	}
	if (ncol(deconvolved.profiles) != length(purity.estimates)) {
		stop('There are different numbers of samples deconvolved.profiles and purity.estimates -- make sure the length of purity estimates is the same as the number of columns in the matrices');
	}

    # make sure data is not log transformed
	if (max(max(tumor.profiles)) < 30) {
		warning('Maximum element in matrix tumor.profiles is less than 30 -- make sure data is in normal (not log) space, otherwise output is wrong.');
	}
	if (max(max(deconvolved.profiles)) < 30) {
		warning('Maximum element in matrix deconvolved.profiles is less than 30 -- make sure data is in normal (not log) space, otherwise output is wrong.');
	}

	# make sure purity estimates are proportions of the overall bulk tumour (ie between 0 and 1)
	if (max(purity.estimates) > 1 | min(purity.estimates) < 0) {
		warning('purity.estimates are not between 0 and 1 and should be');
	}

	# set up matrix for output
	tac.profiles <- matrix(NA,ncol=ncol(deconvolved.profiles),nrow=nrow(deconvolved.profiles));
	rownames(tac.profiles) <- rownames(deconvolved.profiles);
	colnames(tac.profiles) <- colnames(deconvolved.profiles);
	
	# calculate tumour adjacent cell expression profiles
	for (i in 1:ncol(tac.profiles)) {
		tac.profiles[,i] <- (tumor.profiles[,i] - purity.estimates[i]*deconvolved.profiles[,i])/(1-purity.estimates[i]);
	}

	return(tac.profiles);
}
