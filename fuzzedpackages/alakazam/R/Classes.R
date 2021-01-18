# Classes, Generics and Methods

#### Generics ####

# @exportMethod print
#setGeneric("print")

# @exportMethod plot
#setGeneric("plot")

setClassUnion("DfNULL", members=c("data.frame", "NULL"))
setClassUnion("CharNULL", members=c("character", "NULL"))

#### Diversity classes ####

#' S4 class defining a clonal abundance curve
#'
#' \code{AbundanceCurve} defines clonal abundance values.
#' 
#' @slot  abundance  data.frame with relative clonal abundance data and confidence intervals,
#'                   containing the following columns:
#'                   \itemize{
#'                     \item  \code{group}:  group identifier.
#'                     \item  \code{clone_id} or \code{CLONE}:  clone identifier. 
#'                     \item  \code{p}:      relative abundance of the clone.
#'                     \item  \code{lower}:  lower confidence inverval bound.
#'                     \item  \code{upper}:  upper confidence interval bound.
#'                     \item  \code{rank}:   the rank of the clone abundance.
#'                   }
#' @slot  bootstrap  data.frame of bootstrapped clonal distributions. 
#' @slot  clone_by   string specifying the name of the clone column. 
#' @slot  group_by   string specifying the name of the grouping column. 
#' @slot  groups     vector specifying the names of unique groups in group column. 
#' @slot  n          numeric vector indication the number of sequences sampled in each group.
#' @slot  nboot      numeric specifying the number of bootstrap iterations to use.  
#' @slot  ci         confidence interval defining the upper and lower bounds 
#'                   (a value between 0 and 1).
#' 
#' @name         AbundanceCurve-class
#' @rdname       AbundanceCurve-class
#' @aliases      AbundanceCurve
#' @exportClass  AbundanceCurve
setClass("AbundanceCurve", 
         slots=c(abundance="data.frame",
                 bootstrap="data.frame",
                 clone_by="character",
                 group_by="character",
                 groups="character",
                 n="numeric", 
                 ci="numeric",
                 nboot="numeric"))

#' S4 class defining a diversity curve 
#'
#' \code{DiversityCurve} defines diversity (\eqn{D}) scores over multiple diversity 
#' orders (\eqn{Q}).
#' 
#' @slot  diversity  data.frame defining the diversity curve with the following columns:
#'                   \itemize{
#'                     \item  \code{group}:    group label.
#'                     \item  \code{q}:        diversity order.
#'                     \item  \code{d}:        mean diversity index over all bootstrap 
#'                                             realizations.
#'                     \item  \code{d_sd}:     standard deviation of the diversity index 
#'                                             over all bootstrap realizations.
#'                     \item  \code{d_lower}:  diversity lower confidence inverval bound.
#'                     \item  \code{d_upper}:  diversity upper confidence interval bound.
#'                     \item  \code{e}:        evenness index calculated as \code{D} 
#'                                             divided by \code{D} at \code{Q=0}.
#'                     \item  \code{e_lower}:  evenness lower confidence inverval bound.
#'                     \item  \code{e_upper}:  eveness upper confidence interval bound.
#'                   }
#' @slot  tests    data.frame describing the significance test results with columns:
#'                 \itemize{
#'                   \item  \code{test}:        string listing the two groups tested.
#'                   \item  \code{delta_mean}:  mean of the \eqn{D} bootstrap delta 
#'                                              distribution for the test.
#'                   \item  \code{delta_sd}:    standard deviation of the \eqn{D} 
#'                                              bootstrap delta distribution for the test.
#'                   \item  \code{pvalue}:      p-value for the test.
#'                 }
#' @slot  group_by   string specifying the name of the grouping column in diversity calculation.
#' @slot  groups     vector specifying the names of unique groups in group column in diversity calculation.
#' @slot  method     string specifying the type of diversity calculated. 
#' @slot  q          vector of diversity hill diversity indices used for computing diversity.  
#' @slot  n          numeric vector indication the number of sequences sampled in each group.
#' @slot  ci         confidence interval defining the upper and lower bounds 
#'                   (a value between 0 and 1).
#' 
#' @name         DiversityCurve-class
#' @rdname       DiversityCurve-class
#' @aliases      DiversityCurve
#' @exportClass  DiversityCurve
setClass("DiversityCurve", 
         slots=c(diversity="data.frame", 
                 tests="DfNULL",
                 method="character",
                 group_by="character",
                 groups="character",
                 q="numeric",
                 n="numeric", 
                 ci="numeric"))
              
#### Diversity methods ####

#' @param    x    AbundanceCurve object
#' 
#' @rdname   AbundanceCurve-class
#' @aliases  AbundanceCurve-method
#' @export
setMethod("print", c(x="AbundanceCurve"), function(x) { print(x@abundance) })

#' @param    y    ignored.
#' @param    ...  arguments to pass to \link{plotDiversityCurve}.
#' 
#' @rdname   AbundanceCurve-class
#' @aliases  AbundanceCurve-method
#' @export
setMethod("plot", c(x="AbundanceCurve", y="missing"),
          function(x, y, ...) { plotAbundanceCurve(x, ...) })

#' @param    x    DiversityCurve object
#' 
#' @rdname   DiversityCurve-class
#' @aliases  DiversityCurve-method
#' @export
setMethod("print", c(x="DiversityCurve"), function(x) { print(x@diversity) })

#' @param    y    diversity order to plot (q).
#' @param    ...  arguments to pass to \link{plotDiversityCurve} or \link{plotDiversityTest}.
#' 
#' @rdname   DiversityCurve-class
#' @aliases  DiversityCurve-method
#' @export
setMethod("plot", c(x="DiversityCurve", y="missing"),
          function(x, y, ...) { plotDiversityCurve(x, ...) })

#' @rdname   DiversityCurve-class
#' @aliases  DiversityCurve-method
#' @export
setMethod("plot", c(x="DiversityCurve", y="numeric"),
          function(x, y, ...) { plotDiversityTest(x, y, ...) })

#### Lineage classes ####

#' S4 class defining a clone
#' 
#' \code{ChangeoClone} defines a common data structure for perform lineage recontruction
#' from Change-O data.
#' 
#' @slot     data      data.frame containing sequences and annotations. Contains the
#'                     columns \code{SEQUENCE_ID} and \code{SEQUENCE}, as well as any additional 
#'                     sequence-specific annotation columns.
#' @slot     clone     string defining the clone identifier.
#' @slot     germline  string containing the germline sequence for the clone.
#' @slot     v_gene    string defining the V segment gene call.
#' @slot     j_gene    string defining the J segment gene call.
#' @slot     junc_len  numeric junction length (nucleotide count).
#' 
#' @seealso  See \link{makeChangeoClone} and \link{buildPhylipLineage} for use.
#'           
#' @name         ChangeoClone-class
#' @rdname       ChangeoClone-class
#' @aliases      ChangeoClone
#' @exportClass  ChangeoClone
setClass("ChangeoClone", 
         slots=c(data="data.frame",
                 clone="character",
                 germline="character", 
                 v_gene="character", 
                 j_gene="character", 
                 junc_len="numeric"))


#### Topology classes ####

#' S4 class defining edge significance
#'
#' \code{MRCATest} defines the significance of enrichment for annotations appearing at
#' the MRCA of the tree.
#' 
#' @slot  tests         data.frame describing the significance test results with columns:
#'                      \itemize{
#'                        \item  \code{annotation}:  annotation value.
#'                        \item  \code{count}:       observed count of MRCA positions 
#'                                                   with the given annotation.
#'                        \item  \code{expected}:    expected mean count of MRCA occurance
#'                                                   for the annotation.
#'                        \item  \code{pvalue}:      one-sided p-value for the hypothesis that 
#'                                                   the observed annotation abundance is greater 
#'                                                   than expected.
#'                      }
#' @slot  permutations  data.frame containing the raw permutation test data with columns:
#'                      \itemize{
#'                        \item  \code{annotation}:  annotation value.
#'                        \item  \code{count}:       count of MRCA positions with the 
#'                                                   given annotation.
#'                        \item  \code{iter}:        numerical index define which 
#'                                                   permutation realization each 
#'                                                   observation corresponds to.
#'                      }
#' @slot  nperm         number of permutation realizations.
#' 
#' @name         MRCATest-class
#' @rdname       MRCATest-class
#' @aliases      MRCATest
#' @exportClass  MRCATest
setClass("MRCATest", 
         slots=c(tests="data.frame",
                 permutations="data.frame",
                 nperm="numeric"))


#' S4 class defining edge significance
#'
#' \code{EdgeTest} defines the significance of parent-child annotation enrichment.
#' 
#' @slot  tests         data.frame describing the significance test results with columns:
#'                      \itemize{
#'                        \item  \code{parent}:    parent node annotation.
#'                        \item  \code{child}:     child node annotation
#'                        \item  \code{count}:     count of observed edges with the given 
#'                                                 parent-child annotation set.
#'                        \item  \code{expected}:  mean count of expected edges for the 
#'                                                 given parent-child relationship.
#'                        \item  \code{pvalue}:    one-sided p-value for the hypothesis that 
#'                                                  the observed edge abundance is greater 
#'                                                  than expected.
#'                      }
#' @slot  permutations  data.frame containing the raw permutation test data with columns:
#'                      \itemize{
#'                        \item  \code{parent}:  parent node annotation.
#'                        \item  \code{child}:   child node annotation
#'                        \item  \code{count}:   count of edges with the given parent-child 
#'                                               annotation set.
#'                        \item  \code{iter}:    numerical index define which permutation
#'                                               realization each observation corresponds 
#'                                               to.
#'                      }
#' @slot  nperm         number of permutation realizations.
#' 
#' @name         EdgeTest-class
#' @rdname       EdgeTest-class
#' @aliases      EdgeTest
#' @exportClass  EdgeTest
setClass("EdgeTest", 
         slots=c(tests="data.frame",
                 permutations="data.frame",
                 nperm="numeric"))


#### Topology methods ####

#' @param    x    MRCATest object.
#' 
#' @rdname   MRCATest-class
#' @aliases  MRCATest-method
#' @export
setMethod("print", c(x="MRCATest"), function(x) { print(x@tests) })

#' @param    y    ignored.
#' @param    ...  arguments to pass to \link{plotMRCATest}.
#' 
#' @rdname   MRCATest-class
#' @aliases  MRCATest-method
#' @export
setMethod("plot", c(x="MRCATest", y="missing"),
          function(x, y, ...) { plotMRCATest(x, ...) })

#' @param    x    EdgeTest object.
#' 
#' @rdname   EdgeTest-class
#' @aliases  EdgeTest-method
#' @export
setMethod("print", c(x="EdgeTest"), function(x) { print(x@tests) })

#' @param    y    ignored.
#' @param    ...  arguments to pass to \link{plotEdgeTest}.
#' 
#' @rdname   EdgeTest-class
#' @aliases  EdgeTest-method
#' @export
setMethod("plot", c(x="EdgeTest", y="missing"),
          function(x, y, ...) { plotEdgeTest(x, ...) })