# Deprecated and defunct functions

#' @include Classes.R
#' @include Diversity.R
NULL

#### Deprecated ####

#' Generate a clonal diversity index curve
#'
#' \code{rarefyDiversity} divides a set of clones by a group annotation,
#' resamples the sequences from each group, and calculates diversity
#' scores (\eqn{D}) over an interval of diversity orders (\eqn{q}).
#' 
#' @param    data      data.frame with Change-O style columns containing clonal assignments.
#' @param    group     name of the \code{data} column containing group identifiers.
#' @param    clone     name of the \code{data} column containing clone identifiers.
#' @param    copy      name of the \code{data} column containing copy numbers for each 
#'                     sequence. If \code{copy=NULL} (the default), then clone abundance
#'                     is determined by the number of sequences. If a \code{copy} column
#'                     is specified, then clone abundances is determined by the sum of 
#'                     copy numbers within each clonal group.
#' @param    min_q     minimum value of \eqn{q}.
#' @param    max_q     maximum value of \eqn{q}.
#' @param    step_q    value by which to increment \eqn{q}.
#' @param    min_n     minimum number of observations to sample.
#'                     A group with less observations than the minimum is excluded.
#' @param    max_n     maximum number of observations to sample. If \code{NULL} then no 
#'                     maximum is set.
#' @param    ci        confidence interval to calculate; the value must be between 0 and 1.
#' @param    nboot     number of bootstrap realizations to generate.
#' @param    uniform   if \code{TRUE} then uniformly resample each group to the same 
#'                     number of observations. If \code{FALSE} then allow each group to
#'                     be resampled to its original size or, if specified, \code{max_size}.
#' @param    progress  if \code{TRUE} show a progress bar.
#' 
#' @return   A \link{DiversityCurve} object summarizing the diversity scores.
#' 
#' @details
#' Clonal diversity is calculated using the generalized diversity index (Hill numbers) 
#' proposed by Hill (Hill, 1973). See \link{calcDiversity} for further details.
#'
#' Diversity is calculated on the estimated complete clonal abundance distribution.
#' This distribution is inferred by using the Chao1 estimator to estimate the number
#' of seen clones, and applying the relative abundance correction and unseen clone
#' frequency described in Chao et al, 2015.
#'
#' To generate a smooth curve, \eqn{D} is calculated for each value of \eqn{q} from
#' \code{min_q} to \code{max_q} incremented by \code{step_q}.  When \code{uniform=TRUE}
#' variability in total sequence counts across unique values in the \code{group} column 
#' is corrected by repeated resampling from the estimated complete clonal distribution to a 
#' common number of sequences.
#' 
#' The diversity index (\eqn{D}) for each group is the mean value of over all resampling 
#' realizations. Confidence intervals are derived using the standard deviation of the 
#' resampling realizations, as described in Chao et al, 2015.
#' 
#' @references
#' \enumerate{
#'   \item  Hill M. Diversity and evenness: a unifying notation and its consequences. 
#'            Ecology. 1973 54(2):427-32.
#'   \item  Chao A. Nonparametric Estimation of the Number of Classes in a Population. 
#'            Scand J Stat. 1984 11, 265270.
#'   \item  Chao A, et al. Rarefaction and extrapolation with Hill numbers: 
#'            A framework for sampling and estimation in species diversity studies. 
#'            Ecol Monogr. 2014 84:45-67.
#'   \item  Chao A, et al. Unveiling the species-rank abundance distribution by 
#'            generalizing the Good-Turing sample coverage theory. 
#'            Ecology. 2015 96, 11891201.
#' }
#'  
#' @seealso  \link{alphaDiversity}
#' 
#' @examples
#' \dontrun{
#' # Group by sample identifier
#' div <- rarefyDiversity(ExampleDb, "sample_id", step_q=1, max_q=10, nboot=100)
#' plotDiversityCurve(div, legend_title="Sample")
#'                    
#' # Grouping by isotype rather than sample identifier
#' div <- rarefyDiversity(ExampleDb, "c_call", min_n=40, step_q=1, max_q=10, 
#'                        nboot=100)
#' plotDiversityCurve(div, legend_title="Isotype")
#' }
#' @export
rarefyDiversity <- function(data, group, clone="CLONE", copy=NULL,
                            min_q=0, max_q=4, step_q=0.05, min_n=30, max_n=NULL,
                            ci=0.95, nboot=2000, uniform=TRUE, progress=FALSE) {
    .Deprecated("alphaDiversity")
    bootstrap_obj <- estimateAbundance(data, group=group, clone=clone, copy=copy, nboot=nboot, min_n=min_n, max_n=max_n, uniform=uniform, ci=ci)
    diversity_obj <- alphaDiversity(bootstrap_obj, ci=ci, min_q=min_q, max_q=max_q, step_q)

    return(diversity_obj)
}


#' Pairwise test of the diversity index
#' 
#' \code{testDiversity} performs pairwise significance tests of the diversity index 
#' (\eqn{D}) at a given diversity order (\eqn{q}) for a set of annotation groups via
#' rarefaction and bootstrapping.
#'
#' @param    data      data.frame with Change-O style columns containing clonal assignments.
#' @param    q         diversity order to test.
#' @param    group     name of the \code{data} column containing group identifiers.
#' @param    clone     name of the \code{data} column containing clone identifiers.
#' @param    copy      name of the \code{data} column containing copy numbers for each 
#'                     sequence. If \code{copy=NULL} (the default), then clone abundance
#'                     is determined by the number of sequences. If a \code{copy} column
#'                     is specified, then clone abundances is determined by the sum of 
#'                     copy numbers within each clonal group.
#' @param    min_n     minimum number of observations to sample.
#'                     A group with less observations than the minimum is excluded.
#' @param    max_n     maximum number of observations to sample. If \code{NULL} the maximum
#'                     if automatically determined from the size of the largest group.
#' @param    nboot     number of bootstrap realizations to perform.
#' @param    ci        confidence interval to calculate; the value must be between 0 and 1.
#' @param    progress  if \code{TRUE} show a progress bar.
#' 
#' @return   A \link{DiversityCurve} object containing slot test with p-values and summary 
#'             statistics.
#' 
#' @details
#' Clonal diversity is calculated using the generalized diversity index proposed by 
#' Hill (Hill, 1973). See \link{calcDiversity} for further details.
#' 
#' Diversity is calculated on the estimated complete clonal abundance distribution.
#' This distribution is inferred by using the Chao1 estimator to estimate the number
#' of seen clones, and applying the relative abundance correction and unseen clone
#' frequency described in Chao et al, 2014.
#'
#' Variability in total sequence counts across unique values in the \code{group} column is 
#' corrected by repeated resampling from the estimated complete clonal distribution to 
#' a common number of sequences. The diversity index estimate (\eqn{D}) for each group is 
#' the mean value of over all bootstrap realizations. 
#' 
#' Significance of the difference in diversity index (\eqn{D}) between groups is tested by 
#' constructing a bootstrap delta distribution for each pair of unique values in the 
#' \code{group} column. The bootstrap delta distribution is built by subtracting the diversity 
#' index \eqn{Da} in \eqn{group-a} from the corresponding value \eqn{Db} in \eqn{group-b}, 
#' for all bootstrap realizations, yeilding a distribution of \code{nboot} total deltas; where 
#' \eqn{group-a} is the group with the greater mean \eqn{D}. The p-value for hypothesis 
#' \eqn{Da  !=  Db} is the value of \eqn{P(0)} from the empirical cumulative distribution 
#' function of the bootstrap delta distribution, multiplied by 2 for the two-tailed correction.
#' 
#' @note
#' This method may inflate statistical significance when clone sizes are uniformly small,
#' such as when most clones sizes are 1, sample size is small, and \code{max_n} is near
#' the total count of the smallest data group. Use caution when interpreting the results 
#' in such cases. We are currently investigating this potential problem.
#' 
#' @references
#' \enumerate{
#'   \item  Hill M. Diversity and evenness: a unifying notation and its consequences. 
#'            Ecology. 1973 54(2):427-32.
#'   \item  Chao A. Nonparametric Estimation of the Number of Classes in a Population. 
#'            Scand J Stat. 1984 11, 265270.
#'   \item  Wu Y-CB, et al. Influence of seasonal exposure to grass pollen on local and 
#'            peripheral blood IgE repertoires in patients with allergic rhinitis. 
#'            J Allergy Clin Immunol. 2014 134(3):604-12.
#'   \item  Chao A, et al. Rarefaction and extrapolation with Hill numbers: 
#'            A framework for sampling and estimation in species diversity studies. 
#'            Ecol Monogr. 2014 84:45-67.
#'   \item  Chao A, et al. Unveiling the species-rank abundance distribution by 
#'            generalizing the Good-Turing sample coverage theory. 
#'            Ecology. 2015 96, 11891201.
#' }
#' 
#' @seealso  \link{alphaDiversity}
#' 
#' @examples  
#' \dontrun{
#' # Groups under the size threshold are excluded and a warning message is issued.
#' testDiversity(ExampleDb, "sample_id", q=0, min_n=30, nboot=100)
#' }
#' 
#' @export
testDiversity <- function(data, q, group, clone="CLONE", copy=NULL,
                          min_n=30, max_n=NULL, nboot=2000, progress=FALSE, ci=0.95) {
    .Deprecated("alphaDiversity")

    abundance_obj <- estimateAbundance(data, group=group, clone=clone, copy=copy, nboot=nboot, min_n=min_n, max_n=max_n, ci=ci)
    diversity_obj <- alphaDiversity(abundance_obj, min_q=q, max_q=q, step_q=1, ci=ci)

    return(diversity_obj)
}

#### Defunct ####

