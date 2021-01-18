.check_gp_file_name <- function(inputFile) {
  if (!is.character(inputFile)) stop("'inputFile' must be a character string.")
  if ( ! tools::file_ext(inputFile) %in% c("","txt") ) warning("Dangerous input file name extension:\n   the genepop code expects either extension '.txt' or no extension.")
  if( ! file.exists(inputFile) ) stop("'inputFile' not found ")
  if( file.access(inputFile,mode=4) ) stop("'inputFile' not readable (as tested by file.access(<.>,mode=4)).")
}

.remove_temp_files <- function() { ## currently this code provides documentation, but the function has no use
  if (FALSE) {
    ## remove temp HW files
    ## should not be necessary as GenClean_HW() has done it.
    abyss <- file.remove(dir(pattern = "^P[0-9]+_L[0-9]+$", full.names = TRUE ))
    abyss <- file.remove(dir(pattern = "^popc[0-9]+$", full.names = TRUE ))
    abyss <- file.remove(dir(pattern = "^locc[0-9]+$", full.names = TRUE ))
    abyss <- file.remove(dir(pattern = "^poploc$", full.names = TRUE ))
    ## remove temp Fst files
    ## should not be necessary as isolde_etc() has done it.
    abyss <- file.remove(dir(pattern = "^LOCUS[0-9]+$", full.names = TRUE ))
  }
  ## voir outloc outpop outtot
}

.returnInfo <- function(info, verbose, prestring = "Results are stored in file") {
  if (file.exists(info)) {
    resu <- paste(prestring, info)
    if (verbose) {
      cat(paste(prestring, info), "\n")
    }
    info <- structure(info,class=c(class(info),"GP_outfile"))
    invisible(info)
  } else {
    warning(paste("file",info,"*not* created. Check for an earlier problem."))
    return(NULL)
  }
}

#' @name Hardy-Weinberg
#' @title Tests of Hardy-Weinberg genotypic proportions
#' @description Compute variants of the exact conditional test for Hardy-Weinberg genotypic proportions.  The tests differ by their test statistics. \code{HWtable_analysis} handles a single table of genotype counts, and \code{test_HW} requires a standard genepop input file. See \href{../doc/all-menu-options.html#option-1-hardy-weinberg-hw-exact-tests}{this section} of the Genepop-executable documentation for more information on the statistical methods.
#' @param inputFile character: The path of the input file. For \code{test_HW}, this file should be in Genepop format. For \code{HWtable_analysis}, it should be in ad hoc format illustrated by sample file \code{Rhesus.txt} used in the Examples section, and further detailed in \href{../doc/all-menu-options.html#analyzing-a-single-genotypic-matrix}{this section} of the Genepop-executable documentation.
#' @param which character: \code{'Proba'}, \code{'excess'}, and \code{'deficit'} to perform the probability test, score test for excess, and score tests for deficit, respectively, in each population and for each locus. \code{test_HW} additionally handles \code{'global excess'} and  \code{'global deficit'} for global tests for all loci and/or all populations, and \code{HWtable_analysis} additionally handles \code{'Fis'} to report basic information (allele frequencies and Fis).
#' @param outputFile character: The path of the output file
#' @param settingsFile character: The path of the settings file
#' @param enumeration logical: whether to compute the complete enumeration test for samples with less than 5 alleles
#' @param dememorization integer: length of dememorization step of Markov chain algorithm
#' @param batches integer: Number of batches
#' @param iterations integer: Iterations per batch
#' @param verbose logical: whether to print some information
#' @return The path of the output file is returned invisibly.
#' @examples locinfile <- genepopExample('sample.txt')
#' test_HW(locinfile, which='deficit', 'sample.txt.D')
#' if ( ! interactive()) clean_workdir(otherfiles='sample.txt')
test_HW <- function(inputFile, which = "Proba", outputFile = "", settingsFile = "", enumeration = FALSE, dememorization = 10000, 
    batches = 20, iterations = 5000, verbose = interactive()) {
    .check_gp_file_name(inputFile)
    if (which == "Proba") {
        if (settingsFile == "") {
            resu <- RHWEachLocusEachPopulationProbability( inputFile, outputFile, enumeration, 
                dememorization, batches, iterations)
        } else {
            resu <- RHWEachLocusEachPopulationProbabilityWithSettingsFile( inputFile, 
                outputFile, settingsFile)
        }
    } else if (which == "excess") {
        if (settingsFile == "") {
            resu <- RHWEachLocusEachPopulationHE( inputFile, outputFile, enumeration, 
                dememorization, batches, iterations)
        } else {
            resu <- RHWEachLocusEachPopulationHEWithSettingsFile( inputFile, outputFile, 
                settingsFile)
        }
    } else if (which == "deficit") {
        if (settingsFile == "") {
            resu <- RHWEachLocusEachPopulationHD (inputFile, outputFile, enumeration, 
                dememorization, batches, iterations)
        } else {
            resu <- RHWEachLocusEachPopulationHDWithSettingsFile( inputFile, outputFile, 
                settingsFile)
        }
    } else if (which == "global excess") {
        if (settingsFile == "") {
            resu <- RHWGlobalHE(inputFile, outputFile, dememorization, batches, iterations)
        } else {
            resu <- RHWGlobalHEWithSettingsFile( inputFile, outputFile, settingsFile)
        }
    } else if (which == "global deficit") {
        if (settingsFile == "") {
            resu <- RHWGlobalHD( inputFile, outputFile, dememorization, batches, iterations)
        } else {
            resu <- RHWGlobalHDWithSettingsFile( inputFile, outputFile, settingsFile)
        }
    } else stop("invalid 'which' value")
    .returnInfo(resu, verbose = verbose)
}

#' @rdname Hardy-Weinberg
#' @examples # Example in Guo & Thompson 1992 Table 5
#' locinfile <- genepopExample('Rhesus.txt')
#' outfile <- HWtable_analysis(locinfile,which='Proba',batches = 1000,iterations = 1000)
#' readLines(outfile)[21]
#' #clean_workdir(otherfiles='Rhesus.txt')
HWtable_analysis <- function(inputFile, which = "Proba", settingsFile = "", enumeration = FALSE, dememorization = 10000, batches = 20, 
    iterations = 5000, verbose = interactive()) {
    if (which == "Fis") {
        resu <- RHWtableAlleleFrequenciesExpectedGenotypesFis( inputFile)
    } else if (which == "Proba") {
        if (settingsFile == "") {
            resu <- RHWtableProbability(inputFile, enumeration, dememorization, batches, 
                iterations)
        } else {
            resu <- RHWtableProbabilityWithSettingsFile(inputFile, settingsFile)
        }
    } else if (which == "excess") {
        if (settingsFile == "") {
            resu <- RHWtableHE(inputFile, enumeration, dememorization, batches, iterations)
        } else {
            resu <- RHWtableHEWithSettingsFile(inputFile, settingsFile)
        }
    } else {
        if (settingsFile == "") {
            resu <- RHWtableHD(inputFile, enumeration, dememorization, batches, iterations)
        } else {
            resu <- RHWtableHDWithSettingsFile(inputFile, settingsFile)
        }
    }
    .returnInfo(resu, prestring = "Results are appended to file", verbose = verbose)
}


#' @name Linkage
#' @title Tables and exact test for genotypic linkage disequilibrium
#' @description Exact test for each pair of loci in each population. See \href{../doc/all-menu-options.html#option-2-tests-and-tables-for-linkage-disequilibrium}{this section} of the Genepop executable documentation for more information on the statistical methods.
#' @param inputFile The path of the input file, in Genepop format
#' @param outputFile character: The path of the output file
#' @param settingsFile character: The path of the settings file
#' @param dememorization integer: length of dememorization step of Markov chain algorithm
#' @param batches integer: Number of batches
#' @param iterations integer: Iterations per batch
#' @param verbose logical: whether to print some information
#' @return The path of the output file is returned invisibly.
#' @examples locinfile <- genepopExample('sample.txt')
#' test_LD(locinfile,'sample.txt.DIS')
#' if ( ! interactive()) clean_workdir(otherfiles='sample.txt')
test_LD <- function(inputFile, outputFile = "", settingsFile = "", dememorization = 10000, batches = 100, iterations = 5000, 
    verbose = interactive()) {
  .check_gp_file_name(inputFile)
  if ( ! tools::file_ext(inputFile) %in% c("","txt") ) stop("Wrong input file name extension: it should either have extension '.txt' or no extension.")
    if (settingsFile == "") {
        resu <- RGDEachPairLociEachPopulation(inputFile, outputFile, dememorization, 
            batches, iterations)
    } else {
        resu <- RGDEachPairLociEachPopulationWithSettingsFile(inputFile, outputFile, 
            settingsFile)
    }
    .returnInfo(resu, verbose = verbose)
}

#' @rdname Linkage
#' @examples locinfile <- genepopExample('sample.txt')
#' write_LD_tables(locinfile,'sample.txt.TAB')
#' if ( ! interactive()) clean_workdir(otherfiles='sample.txt')
write_LD_tables <- function(inputFile, outputFile = "", verbose = interactive()) {
    resu <- RGDGenotypicContingency(inputFile, outputFile)
    .returnInfo(resu, verbose = verbose)
}

#' @name Differentiation
#' @title Tests of genic and genotypic differentiation
#' @description Exact conditional contingency-table tests for genic or genotypic differentiation. A single test for all populations, or distinct tests for all pairs of populations, may be computed. See \href{../doc/all-menu-options.html#option-3-population-differentiation}{this section} of the Genepop executable documentation for more information on the statistical methods.
#' @param inputFile The path of the input file, in Genepop format
#' @param genic logical: whether to perform genic or genotypic tests
#' @param pairs logical: whether to test differentiation between all pairs of ppulation, or to perform a single global test
#' @param outputFile character: The path of the output file
#' @param settingsFile character: The path of the settings file
#' @param dememorization integer: length of dememorization step of Markov chain algorithm
#' @param batches integer: Number of batches
#' @param iterations integer: Iterations per batch
#' @param verbose logical: whether to print some information
#' @return The path of the output file is returned invisibly.
#' @examples locinfile <- genepopExample('sample.txt')
#' test_diff(locinfile,outputFile='sample.txt.GE')
#' if ( ! interactive()) clean_workdir(otherfiles='sample.txt')
test_diff <- function(inputFile, genic = TRUE, pairs = FALSE, outputFile = "", settingsFile = "", dememorization = 10000, 
    batches = 100, iterations = 5000, verbose = interactive()) {
  .check_gp_file_name(inputFile)
  if (genic) {
        if (pairs) {
            if (settingsFile == "") {
                resu <- RPDGenicAllPairPopulationDifferentiation(inputFile, outputFile, 
                  dememorization, batches, iterations)
            } else {
                resu <- RPDGenicAllPairPopulationDifferentiationWithSettingsFile(inputFile, 
                  outputFile, settingsFile)
            }
        } else {
            if (settingsFile == "") {
                resu <- RPDGenicAllPopulationDifferentiation(inputFile, outputFile, 
                  dememorization, batches, iterations)
            } else {
                resu <- RPDGenicAllPopulationDifferentiationWithSettingsFile(inputFile, 
                  outputFile, settingsFile)
            }
        }
    } else {
        if (pairs) {
            if (settingsFile == "") {
                resu <- RPDGenotypicAllPairPopulationDifferentiation(inputFile, outputFile, 
                  dememorization, batches, iterations)
            } else {
                resu <- RPDGenotypicAllPairPopulationDifferentiationWithSettingsFile(
                  inputFile, outputFile, settingsFile)
            }
        } else {
            if (settingsFile == "") {
                resu <- RPDGenotypicAllPopulationDifferentiation(inputFile, outputFile, 
                  dememorization, batches, iterations)
            } else {
                resu <- RPDGenotypicAllPopulationDifferentiationWithSettingsFile(inputFile, 
                  outputFile, settingsFile)
            }
        }
    }
    .returnInfo(resu, verbose = verbose)
}

#' @name Contingency-test
#' @aliases struc
#' @title Exact test on a single contingency table
#' @description Performs an exact conditional contingency-table test. There are many other ways of doing this in R but this function replicates the functionality of earlier genepop code analysing a contingency table provided in a file with ad hoc format. See \href{../doc/all-menu-options.html#analyzing-a-single-contingency-table}{this section} of the Genepop executable documentation for more information on the statistical methods.
#' @param inputFile character: The path of the input file. This file should be in an ad hoc format
#' @param settingsFile character: The path of the settings file
#' @param dememorization integer: length of dememorization step of Markov chain algorithm
#' @param batches integer: Number of batches
#' @param iterations integer: Iterations per batch
#' @param verbose logical: whether to print some information
#' @return The path of the output file is returned invisibly.
#' @examples locinfile <- genepopExample('structest.txt')
#' struc(locinfile)
#' if ( ! interactive()) clean_workdir(otherfiles='structest.txt')
struc <- function(inputFile, settingsFile = "", dememorization = 10000, batches = 100, iterations = 5000, verbose = interactive()) {
    if (settingsFile == "") {
        resu <- RAnalyzingSingleContingencyTable(inputFile, dememorization, batches, 
            iterations)
    } else {
        resu <- RAnalyzingSingleContingencyTableWithSettingsFile(inputFile, settingsFile)
    }
    .returnInfo(resu, prestring = "Results are appended to file", verbose = verbose)
}

#' @name Nm_private
#' @title Private allele method
#' @description Estimation of Nm by private allele method of Slatkin and Barton. See \href{../doc/all-menu-options.html#option-4-private-alleles}{this section} of the Genepop executable documentation for more information on the statistical methods.
#' @param inputFile The path of the input file, in Genepop format
#' @param outputFile character: The path of the output file
#' @param dataType character: The haploid and diploid data
#' @param verbose logical: whether to print some information
#' @return The path of the output file is returned invisibly.
#' @examples locinfile <- genepopExample('sample.txt')
#' Nm_private(locinfile,'sample.txt.PRI')
#' if ( ! interactive()) clean_workdir(otherfiles='sample.txt')
Nm_private <- function(inputFile, outputFile = "", dataType = "Diploid", verbose = interactive()) {
  .check_gp_file_name(inputFile)
  resu <- RNmEstimates(inputFile, outputFile, dataType)
    .returnInfo(resu, verbose = verbose)
}

#' @name basic_info
#' @title Allele and genotype frequencies
#' @description Allele and genotype frequencies per locus and per sample. See \href{../doc/all-menu-options.html#sub-option-1-allele-and-genotype-frequencies}{this section} of the Genepop executable documentation for more information on the statistical methods.
#' @param inputFile The path of the input file, in Genepop format
#' @param outputFile character: The path of the output file
#' @param verbose logical: whether to print some information
#' @return The path of the output file is returned invisibly.
#' @examples locinfile <- genepopExample('sample.txt')
#' basic_info(locinfile,'sample.txt.INF')
#' if ( ! interactive()) clean_workdir(otherfiles='sample.txt')
basic_info <- function(inputFile, outputFile = "", verbose = interactive()) {
  .check_gp_file_name(inputFile)
  resu <- RDescriptifAlleleAndGenotypeFrequenciesPerLocusPerSample(inputFile, outputFile)
    .returnInfo(resu, verbose = verbose)
}

#' @name genedivFis
#' @title Gene diversities and Fis (or rho_IS)
#' @description Evaluates Fis and gene diversities, or related measures based on allele sizes. See \href{../doc/all-menu-options.html#sub-option-2-identity-based-gene-diversities-and-f_mathrmis}{this section} of the Genepop executable documentation for more information on the identity-based statistical methods, and \href{../doc/all-menu-options.html#sub-option-3-allele-size-based-gene-diversities-and-rho_mathrmis}{this one} for allele-size based ones.
#' @param inputFile The path of the input file, in Genepop format
#' @param sizes logical: whether to compute statistics based on allele size, or not.
#' @param outputFile character: The path of the output file
#' @param dataType character: The haploid and diploid data
#' @param verbose logical: whether to print some information
#' @return The path of the output file is returned invisibly.
#' @examples locinfile <- genepopExample('sample.txt')
#' genedivFis(locinfile,outputFile = 'sample.txt.DIV')
#' if ( ! interactive()) clean_workdir(otherfiles='sample.txt')
genedivFis <- function(inputFile, sizes = FALSE, outputFile = "", dataType = "Diploid", verbose = interactive()) {
  .check_gp_file_name(inputFile)
  if (sizes) {
        resu <- RDescriptifGeneDiversitiesAndFisUsingAlleleSize(inputFile, outputFile, 
            dataType)
    } else {
        resu <- RDescriptifGeneDiversitiesAndFisUsingAlleleIdentity(inputFile, outputFile, 
            dataType)
    }
    .returnInfo(resu, verbose = verbose)
}

#' @name Fst
#' @title Fst (or rho_ST) estimation
#' @description Evaluates Fst or related measures based on allele sizes, for all populations of for all pairs of populations. See \href{../doc/all-menu-options.html#sub-options-14-f-statistics-and-rho-statistics}{this section} of the Genepop executable documentation for more information on the statistical methods.
#' @param inputFile The path of the input file, in Genepop format
#' @param sizes logical: whether to estimate allele-size based statistics, or identity-based Fst
#' @param pairs whether to estimate differentiation between all pairs of populations, or to compute a global estimate for all populations
#' @param outputFile character: The path of the output file
#' @param dataType character: The haploid and diploid data
#' @param verbose logical: whether to print some information
#' @return The path of the output file is returned invisibly.
#' @examples locinfile <- genepopExample('sample.txt')
#' Fst(locinfile, outputFile= 'sample.txt.DIV')
#' if ( ! interactive()) clean_workdir(otherfiles='sample.txt')
Fst <- function(inputFile, sizes = FALSE, pairs = FALSE, outputFile = "", dataType = "Diploid", verbose = interactive()) {
  .check_gp_file_name(inputFile)
  if (sizes) {
        if (pairs) {
            resu <- REstimatingSpatialStructureAlleleSizeAllPopulationsPairs(inputFile, 
                outputFile, dataType)
        } else {
            resu <- REstimatingSpatialStructureAlleleSizeAllPopulations(inputFile, outputFile, 
                dataType)
        }
    } else {
        if (pairs) {
            resu <- REstimatingSpatialStructureAlleleIdentyAllPopulationsPairs(inputFile, 
                outputFile, dataType)
        } else {
            resu <- REstimatingSpatialStructureAlleleIdentyAllPopulations(inputFile, 
                outputFile, dataType)
        }
    }
    .returnInfo(resu, verbose = verbose)
}

#' @name IBD
#' @title Isolation by distance
#' @description Estimates isolation by distance by regression of genetic distance to geographical distance. See \href{../doc/all-menu-options.html#sub-option-5-isolation-by-distance-between-individuals}{this section} of the Genepop executable documentation for more information on individual-based analyses and \href{../doc/all-menu-options.html#sub-option-6-isolation-by-distance-between-groups}{this one} for group-based analyses.
#' @param inputFile The path of the input file, in Genepop format
#' @param outputFile character: The path of the output file
#' @param settingsFile character: The path of the settings file
#' @param dataType character: \code{'haploid'} or \code{'diploid'}
#' @param statistic character: The pairwise genetic distance, either \code{'a'} or \code{'e'} for diploid individual data, \code{'a-like'} for haploid individual data, and \code{'F/(1-F)'} or \code{'SingleGeneDiv'} for group data (haploid or diploid)
#' @param geographicScale character: gives either the scale transformation \code{'Log'} or \code{'Linear'}  for geographic distances, or the shape of the habitat \code{'2D'} or \code{'1D'}
#' @param CIcoverage numeric: The coverage probability of confidence intervals
#' @param testPoint numeric: Given value of the slope to be tested
#' @param minimalDistance numeric: The minimal geographic distance
#' @param maximalDistance numeric: The maximal geographic distance
#' @param mantelPermutations numeric: The number of permutations may be specified
#' @param mantelRankTest logical: whether to use ranks in the Mantel test
#' @param verbose logical: whether to print some information
#' @return The path of the output file is returned invisibly.
#' @examples
#' \dontrun{
#' locinfile <- genepopExample('w2.txt')
#' outfile <- ibd(locinfile,'w2.txt.ISO', geographicScale = 'Log', statistic='e')
#' if ( ! interactive()) clean_workdir(otherfiles='w2.txt')
#'
#' locinfile <- genepopExample('PEL1600withCoord.txt')
#' outfile <- ibd(locinfile,'PEL1600withCoord.ISO', statistic = 'SingleGeneDiv',
#'                geographicScale = '1D')
#' if ( ! interactive()) clean_workdir(otherfiles='PEL1600withCoord.txt')
#' }
ibd <- function(inputFile, outputFile = "", settingsFile = "", dataType = "Diploid", statistic = "F/(1-F)", geographicScale = "2D", 
    CIcoverage = 0.95, testPoint = 0, minimalDistance = 1e-04, maximalDistance = 1e+09, mantelPermutations = 1000, mantelRankTest = FALSE, 
    verbose = interactive()) {
    mc <- match.call()
    .check_gp_file_name(inputFile)
    if (statistic %in% c("a", "e", "a-like")) {
        mc[[1L]] <- quote(.GIsolationByDistanceBetweenIndividuals)
        if (statistic == "a-like") 
            mc$statistic <- "default"
    } else if (statistic %in% c("SingleGeneDiv", "F/(1-F)")) {
        mc[[1L]] <- quote(.GIsolationByDistanceBetweenGroups)
        if (statistic == "F/(1-F)") 
            mc$statistic <- "default"
    } else stop("unknown 'statistic'")
    eval(mc, parent.frame())
}

.GIsolationByDistanceBetweenGroups <- function(inputFile, outputFile = "", settingsFile = "",
                                              dataType = "Diploid", statistic = "SingleGeneDiv",
                                              geographicScale = "2D", CIcoverage = 0.95,
                                              testPoint=0, minimalDistance=1e-4,
                                              maximalDistance=1e+09, mantelPermutations=1000,
                                              mantelRankTest=FALSE, verbose = interactive()) {
  if(settingsFile == "") {
    resu <- RIsolationByDistanceBetweenGroups(inputFile, outputFile, dataType, statistic, geographicScale, CIcoverage, testPoint, minimalDistance, maximalDistance, mantelPermutations, mantelRankTest)
  } else {
    resu <- RIsolationByDistanceBetweenGroupsWithSettingsFile(inputFile, outputFile, settingsFile)
  }
  .returnInfo(resu,verbose=verbose)
}


.GIsolationByDistanceBetweenIndividuals <- function(inputFile, outputFile = "", settingsFile = "", dataType = "Diploid", statistic = "e", 
    geographicScale = "2D", CIcoverage = 0.95, testPoint = 0, minimalDistance = 1e-04, maximalDistance = 1e+09, mantelPermutations = 1000, 
    mantelRankTest = FALSE, verbose = interactive()) {
    if (settingsFile == "") {
        resu <- RIsolationByDistanceBetweenIndividuals(inputFile, outputFile, dataType, 
            statistic, geographicScale, CIcoverage, testPoint, minimalDistance, maximalDistance, mantelPermutations, mantelRankTest)
    } else {
        resu <- RIsolationByDistanceBetweenIndividualsWithSettingsFile(inputFile, outputFile, 
            settingsFile)
    }
    .returnInfo(resu, verbose = verbose)
}


#' @name conversion
#' @aliases conversion
#' @title File conversions
#' @description Converts input files from genepop format to some other formats (some maybe only of historical interest): Fstat, two Biosys formats. and linkdos. See \href{../doc/all-menu-options.html#option-7-file-conversions}{this section} of the Genepop executable documentation for more information on the statistical methods.
#' @param inputFile The path of the input file, in Genepop format
#' @param format Character string: must be one of \code{'Fstat'}, \code{'BiosysL'}, \code{'BiosysN'}, or \code{'Linkdos'}
#' @param outputFile character: The path of the output file
#' @param verbose logical: whether to print some information
#' @return The path of the output file is returned invisibly.
#' @examples locinfile <- genepopExample('sample.txt')
#' conversion(locinfile, format='Fstat', 'sample.txt.DAT')
#' if ( ! interactive()) clean_workdir(otherfiles='sample.txt')
conversion <- function(inputFile, format, outputFile = "", verbose = interactive()) {
  .check_gp_file_name(inputFile)
  resu <- switch(format, Fstat = REcumenicismFstat(inputFile, outputFile), 
        BiosysL = REcumenicismBiosysLetter(inputFile, outputFile), 
        BiosysN = REcumenicismBiosysNumber(inputFile, outputFile), 
        Linkdos = REcumenicismLinkdos(inputFile, outputFile), 
        stop("Unknown format: must be one of Fstat, BiosysL, BiosysN, or Linkdos."))
    .returnInfo(resu, verbose = verbose)
}

#' @name nulls
#' @title Estimation of allele frequencies under genotyping failure.
#' @description Estimates allele frequencies (and failure rate if relevant) under dfferent assumptions:  maximum likelihood assuming that there is null allele (default method), maximum likelihood assuming that apparent nulls are technical failures independent of genotype (\code{'ApparentNulls'}), and Brookfield's (1996) estimator (\code{'B96'}). See \href{../doc/all-menu-options.html#sub-option-1-null-alleles}{this section} of the Genepop executable documentation for more information on the statistical methods. Genepop takes the allele with the highest number for a given locus across all populations as the null allele. For example, if you have 4 alleles plus a null allele, a null homozygote individual should be indicated as e.g. \code{0505} or \code{9999} in the input file.
#' @param inputFile The path of the input file, in Genepop format
#' @param outputFile character: The path of the output file
#' @param settingsFile character: The path of the settings file
#' @param nullAlleleMethod character: \code{'ApparentNulls'}, \code{'B96'} or anything else (default method).
#' @param CIcoverage numeric: The coverage probability of confidence interval
#' @param verbose logical: whether to print some information
#' @return The path of the output file is returned invisibly.
nulls <- function(inputFile, outputFile = "", settingsFile = "", nullAlleleMethod = "", CIcoverage = 0.95, verbose = interactive()) {
  .check_gp_file_name(inputFile)
  if (settingsFile == "") {
        resu <- RNullAlleleEstimateAlleleFrequencies(inputFile, outputFile, nullAlleleMethod, 
            CIcoverage)
    } else {
        resu <- RNullAlleleEstimateAlleleFrequenciesWithSettingsFile(inputFile, outputFile, 
            settingsFile)
    }
    .returnInfo(resu, verbose = verbose)
}

#' @name manipulation
#' @title Various data manipulation utilities
#' @description Various procedures described in the linked sections of the Genepop executable documentation: \href{../doc/all-menu-options.html#sub-option-2-diploidisation-of-haploid-data}{diploidize} haploid data, \href{../doc/all-menu-options.html#sub-option-3-relabeling-alleles-names}{relabel_alleles}, \href{../doc/all-menu-options.html#sub-option-6-random-sampling-of-haploid-genotypes-from-diploid-ones}{sample_haploid}, and \href{../doc/all-menu-options.html#sub-options-4-and-5-conversion-of-population-data-to-individual-data}{pop_to_indiv}. The latter procedure converts population samples (several individuals in each population) to individual data. The names given to the individuals in the new file created (names which are to be interpreted as coordinates in a spatial analysis) may be the population coordinates (given as the name of the last individual in the original data file), or each individual's coordinates (given as the name of each individual in the original data file).
#' @param inputFile The path of the input file, in Genepop format
#' @param outputFile character: The path of the output file
#' @param coordinates character: either \code{'population'} (to use population coordinates) or any other charater string (to use individual coordinates).
#' @param verbose logical: whether to print some information
NULL

#' @rdname manipulation
#' @examples locinfile <- genepopExample('sample.txt')
#' outfile <- diploidize(inputFile = locinfile,outputFile="Dsample.txt")
#' if ( ! interactive()) clean_workdir(c("sample.txt", "Dsample.txt"))
diploidize <- function(inputFile, outputFile = "", verbose = interactive()) {
  .check_gp_file_name(inputFile)
  resu <- RDiploidisationHaploidData(inputFile, outputFile)
    .returnInfo(resu, verbose = verbose)
}

#' @rdname manipulation
relabel_alleles <- function(inputFile, outputFile = "", verbose = interactive()) {
  .check_gp_file_name(inputFile)
  resu <- RRelabelingAlleles(inputFile, outputFile)
    .returnInfo(resu, verbose = verbose)
}

#' @rdname manipulation
pop_to_indiv <- function(inputFile, coordinates, outputFile = "", verbose = interactive()) {
  .check_gp_file_name(inputFile)
  if (coordinates == "population") {
        resu <- RConversionToIndividualDataWithPopulationNames(inputFile, outputFile)
    } else resu <- RConversionToIndividualDataWithIndividualNames(inputFile, outputFile)
    .returnInfo(resu, verbose = verbose)
}

#' @rdname manipulation
sample_haploid <- function(inputFile, outputFile = "", verbose = interactive()) {
  .check_gp_file_name(inputFile)
  resu <- RRandomSamplingOfHaploidGenotypesFromDiploidOnes(inputFile, outputFile)
    .returnInfo(resu, verbose = verbose)
}

#'@rdname genepop-utils
set_restriction <- function(set=FALSE) {
  invisible(Rset_restriction(set))
}

#'@rdname clean_workdir
#' @title Removing files created by Genepop
#' @description This removes \dQuote{temporary files} created by Genepop, but also output files, so it should be used only when one no longer needs the latter files. This function asumes that the input file name contains only alphanumeric, dot, or underscore characters.  
#' @param otherfiles Character vector(s): one or more names of files to be removed and not matched by the other arguments (such as the input file, or some output files not identified by their suffix, as shown in the Example).
#' @param path character vector: path from where files should be removed.
#' @param suffixes Character vector(s): suffixes of files to be removed (useful for output files with readily identifiable suffixes).
#' @param in. boolean: whether to remove the \code{fichier.in} file created by Genepop.
#' @param cmdline boolean: whether to remove the \code{cmdline.txt} file created by Genepop.
#' @examples # Removing files possibly written by other examples in the documentation:
#' clean_workdir(otherfiles=c("sample.txt", "Dsample.txt", "w2.txt", 
#' "PEL1600withCoord.txt", "Rhesus.txt", "structest.txt"))
clean_workdir <- function(otherfiles=NULL, path=".", suffixes=c("GRA", "ISO", "MIG", "PRI", "DAT", "DG", "DIV", "D", "DIS", "FST", 
                                                "NUL", "RHO", "2G2", "G", "GE", "GE2", "INF", "MSD", "TAB", "ST2"), 
                           in.=TRUE, cmdline=TRUE) {
  oldpath <- setwd(dir = path)
  filenames <- dir()
  matchfun <- function(suffix,filenames) {
    pattern <- paste0("^[a-zA-Z0-9_\\.]*\\.",suffix,"$")
    grep(pattern,filenames)
  }
  if (length(suffixes)) {
    matches <- sapply(suffixes, matchfun, filenames=filenames)
    matches <- unlist(matches)
    if (length(matches)) file.remove(filenames[matches])
  }
  if ( ! is.null(otherfiles)) {
    otherfiles <- otherfiles[which(file.exists(otherfiles))]
    if (length(otherfiles)) file.remove(otherfiles)
  }
  if ( in. && file.exists("fichier.in")) file.remove("fichier.in")
  if ( cmdline && file.exists("cmdline.txt")) file.remove("cmdline.txt")
  setwd(oldpath)
}

#'@rdname genepopExample
#' @title Copy an example file from the Genepop package distribution
#' @description This function is used to copy an example file to the user's directory. It should not be used when analysing one's own data!
#' @param filename The name of an example file from the Genepop distribution.
#' @return Returns the \code{filename}
genepopExample <- function(filename) {
  infile <- system.file('extdata', filename, package='genepop')
  if (infile=="") stop("There is no such file in the package 'extdata' directory.")
  locinfile <- filename ## file in user's directory not in R's extdata directory
  check <- file.copy(infile,locinfile,overwrite=TRUE)
  if ( ! check) stop("The file could not be copied into the user's directory.")
  return(filename)
}