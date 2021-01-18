#' Generalized Odds Ratios
#'
#' @aliases print.Genodds
#'
#' @description Performs Agresti's Generalized Odds Ratios (GenOR) for two-group ordinal data.
#'
#' @usage
#' genodds(response, group, strata=NULL,
#'         alpha=0.05,ties="split",
#'         nnt=FALSE,verbose=FALSE,upper=TRUE)
#'
#' @param response A (non-empty) vector. Gives the outcome measure.
#'                 If a factor, level order is used to determine ranking of outcomes.
#' @param group A factor vector of equal length to \code{response}.
#'              Gives the groups for the test. Must take on exactly 2 values.
#' @param strata An optional factor vector of equal length to \code{response}.
#'               Gives treatment blocks to separate comparisons.
#' @param alpha The acceptable type 1 error used in the test.
#' @param ties A string specifying how ties should be treated. See Details.
#' @param nnt A boolean.
#'            If \code{TRUE}, then print number needed to treat in addition to generalized odds ratios.
#' @param verbose A boolean.
#'                If \code{TRUE}, then print both pooled odds and relative risk ratio matrices
#'                regardless of result of statistical test.
#' @param upper A boolean specifying if the upper triangle
#'              of relative risk ratios should be printed.
#'              If \code{FALSE}, lower triangle is used instead.
#'
#' @return A list with class "\code{Genodds}" containing the following:
#' \describe{
#'     \item{pooled_lnodds}{The pooled log(odds).}
#'     \item{pooled_lnconf.int}{(1-\code{alpha})\% Confidence intervals for pooled log(odds).}
#'     \item{pooled_SElnodds}{Standard error of pooled log(odds).}
#'     \item{pooled_SElnnull}{Standard error of pooled log(odds) under the null hypothesis.}
#'     \item{pooled_p}{The p-value of the test of pooled log(odds) = 1.}
#'     \item{pooled_rel_statistic}{Statistic of test that strata odds are equal.}
#'     \item{pooled_rel_p}{p-value for test that strata odds are equal.}
#'     \item{relative_lnodds}{A matrix giving the log of the ratio of odds between strata (generalized relative risk ratio).}
#'     \item{relative_selnodds}{A matrix containing the standard error of the log(relative risk ratio).}
#'     \item{results}{A list containing a summary of each strata measure.}
#'     \item{param.record}{A list containing parameters used in the test.}
#'}
#'
#' @details
#' Agresti's generalized odds ratios (GenOR) calculates the odds that,
#' if a pair of observations are randomly selected from
#' two groups, the outcome in one group is higher than the other.
#' This implementation determines the direction of this comparison
#' using factor levels. Odds are given with reference to
#' observations corresponding to the higher \code{group} level
#' resulting in a higher value in \code{response}.
#' The opposite direction can be calculated by either calculating 1/genodds,
#' or by specifying \code{response=1-response} in function input.
#'
#' If \code{nnt=TRUE}, the Number Needed to Treat (NNT) is printed.
#' NNT is a health economics measure and is related to generalized
#' odds ratios through the formula NNT=1+2/(GenOR-1).
#' It measures the expected number of patients required for a
#' treatment to have impacted a patient's outcome.
#' In this implementation, a positive NNT occurs when GenOR>1
#' and corresponds to the number needed to treat in the higher
#' \code{group} level to observe a higher \code{response} value,
#' while a negative NNT occurs when GenOR<1 and corresponds
#' to the number needed to treat in the higher \code{group}
#' level to observe a lower \code{response} value.
#' If the confidence interval for GenOR straddles 1,
#' the confidence interval for NNT is given as the union of disjoint
#' intervals.
#'
#' \code{ties} changes how ties are treated. If \code{"split"} is provided,
#' then ties are equally split between favoring both groups
#' (following the approach set out by O'Brien et. al. (2006)).
#' If \code{"drop"} is provided, then ties are ignored
#' (following the approach set out by Agresti (1980)).
#' By default, \code{"split"} is used.
#'
#' If \code{strata} is specified, generalized odds ratios are calculated
#' separately for each individual strata. If in-stratum odds ratios are not
#' significantly different from each other (with significance level \code{alpha}),
#' these odds are pooled to obtain a global odds ratio which is adjusted
#' for strata. If in-stratum odds ratios are significantly different, a matrix containing
#' the relative risk ratios between stratum is printed, along with Z-scores
#' corresponding to the significance of these differences.
#' If \code{verbose=TRUE} is supplied, both pooled odds and relative risk ratios
#' are printed regardless of if the between-stratum odds ratios are
#' significantly different.
#'
#' Options \code{verbose}, \code{nnt} and \code{upper}
#' may be re-specified when using print method.
#'
#' @examples
#' # Use the alteplase dataset provided by package and calculate genodds
#' df <- alteplase
#' x <- genodds(df$mRS,df$treat,df$time)
#' x
#' print(x,nnt=TRUE)
#'
#' @references
#' Agresti, A. (1980). Generalized odds ratios for ordinal data.
#' \emph{Biometrics}, 59-67.
#'
#' O'Brien, R. G., & Castelloe, J. (2006, March).
#' Exploiting the link between the Wilcoxon-Mann-Whitney test and a simple odds statistic.
#' In \emph{Thirty-first Annual SAS Users Group International Conference}.
#'
#' Churilov, L., Arnup, S., Johns, H., Leung, T., Roberts,
#' S., Campbell, B. C., Davis, S. M. & Donnan, G. A. (2014).
#' An improved method for simple, assumption-free ordinal analysis of the
#' modified Rankin Scale using generalized odds ratios.
#' \emph{International Journal of Stroke}, 9(8), 999-1005.
#'
#'
#' @export
genodds <- function(response, group, strata=NULL,alpha=0.05,ties="split",
                    nnt=FALSE, verbose=FALSE,upper=TRUE
                    )
{

  # Check inputs are non-empty
  if(length(response)==0)
  {
    stop("Response cannot be empty")
  }
  if(length(group)==0)
  {
    stop("Group cannot be empty")
  }
  if (length(response)!=length(group))
  {
    stop("Response and Group are different lengths")
  }


  # Remove NA values and warn
  if(is.null(strata))
  {
    nMissing <- sum(is.na(response) | is.na(group))
    if(sum(is.na(response) | is.na(group)) > 0)
    {
      warning(sprintf("Dropped %d observations with missing values",nMissing))
      keepList <- !(is.na(response) | is.na(group))
      response <- response[keepList]
      group <- group[keepList]
    }
  }
  else
  {
    nMissing <- sum(is.na(response) | is.na(group) | is.na(strata))
    if(nMissing > 0)
    {
      warning(sprintf("Dropped %d observations with missing values",nMissing))
      keepList <- !(is.na(response) | is.na(group) | is.na(strata) )
      response <- response[keepList]
      group <- group[keepList]
      strata <- strata[keepList]
    }
  }


  if (length(unique(group))!=2)
  {
    stop("Group must take on exactly 2 values")
  }

  # Coerce group and strata to factors
  if(!(class(group) %in% "factor"))
  {
    group=as.factor(group)
  }

  # If no strata is specified,
  # create a dummy stratum containing all data
  if (is.null(strata)){
    strata=rep("All data",length(response))
  } else {

    # Test if strata length equals response and group
    if (length(strata)!=length(group))
    {
      stop("Strata and Response/Group are different lengths")
    }

    if(!(class(strata) %in% "factor"))
    {
      strata=as.factor(strata)
    }
  }


  # Get ties treatment
  # Set up like this to allow future
  # expansion such as assume all ties favour control
  # or all ties favour treatment
  if(ties=="split")
  {
    contr_fav=0.5
  }
  else if (ties=="drop")
  {
    contr_fav=NA
  }
  else
  {
    stop("Invalid ties option specified")
  }

  # The test for an individual block/strata/layer is done in this function.
  # This will be called later on to generate odds for each stratum
  genodds_strata=function(response,group)
  {
    crosstab=as.matrix(table(response,group))

    N=sum(crosstab)
    p=crosstab/N

    Rt=p[,2:1]
    Rs=get_Rs(p)
    Rd=get_Rd(p)
    # Redistribute ties

    if(!is.na(contr_fav))
    {
      Rs=Rs+(1-contr_fav)*Rt
      Rd=Rd+contr_fav*Rt
    }

    Pc=sum(p*Rs)
    Pd=sum(p*Rd)


    odds=Pc/Pd

    SEodds=2/Pd*(sum(p*(odds*Rd-Rs)^2)/N)^0.5
    SElnodds=SEodds/odds

    conf.int=exp( qnorm(c(alpha/2,1-alpha/2),mean=log(odds),sd=SElnodds) )

    # Smooth p across groups and do this again.
    # This code is WET as hell, but whatever

    # This only seems to agree with Stata's
    # genodds routine up to 5 decimal places.
    # It's unclear if this is just a language
    # issue or an actual bug, needs more investigation

    p=apply(p,1,mean)
    p=cbind(p,p)

    Rt=p[,2:1]
    Rs=get_Rs(p)
    Rd=get_Rd(p)

    # Redistribute ties

    Rs=Rs+contr_fav*Rt
    Rd=Rd+contr_fav*Rt

    Pc=sum(p*Rs)
    Pd=sum(p*Rd)

    SEnull=2/Pd*(sum(p*(1*Rd-Rs)^2)/N)^0.5
    SElnnull=SEnull/1

    p=pnorm(abs(log(odds)),sd=SElnodds,lower.tail=FALSE)*2

    out=list(odds=odds,conf.int=conf.int,p=p,SEodds=SEodds,SEnull=SEnull,
             xtab=crosstab)

    return(out)
  }

  # Get results for each strata
  results=by(data.frame(response,group),
             strata,
             function(df) genodds_strata(df$response,df$group)
             )


  ##############
  # Pool results
  pooled_lnodds=do.call("sum",lapply(results,function(x) x$odds^2/x$SEodds^2 * log(x$odds) ))/
                do.call("sum",lapply(results,function(x) x$odds^2/x$SEodds^2 ))

  pooled_SElnodds=sqrt(1/do.call("sum",lapply(results,function(x) x$odds^2/x$SEodds^2)))

  pooled_lnconf.int=qnorm(c(alpha/2,1-alpha/2),mean=pooled_lnodds,sd=pooled_SElnodds)

#   Do we weight by lnodds under null SE?
#   pooled_lnoddsnull=do.call("sum",lapply(results,function(x) x$odds^2 * log(x$odds)/x$SEnull^2))/
#                     do.call("sum",lapply(results,function(x) x$odds^2/x$SEodds^2))

  pooled_SElnnull=sqrt(1/do.call("sum",lapply(results,function(x) 1/x$SEnull^2)))

  pooled_p=pnorm(abs(pooled_lnodds),sd=pooled_SElnodds,lower.tail=FALSE)*2



  ###########################################
  # Get measures of comparisons across layers

  pooled_rel_statistic=do.call("sum",lapply(results,function(x) ( (log(x$odds) - pooled_lnodds)/(x$SEodds/x$odds) )^2 ) )
  pooled_rel_p = pchisq(pooled_rel_statistic,length(results)-1,lower.tail = FALSE)

  crosslayer=lapply(results, function(x) {
             lapply(results, function(y) {
                return(data.frame(lnrel=log(x$odds/y$odds),selnrel= sqrt( (x$SEodds/x$odds)^2 + (y$SEodds/y$odds)^2 )) )
              })
              })

  lnrel=matrix(nrow=length(crosslayer),ncol = length(crosslayer),dimnames = list(names(crosslayer),names(crosslayer)))
  SElnrel=matrix(nrow=length(crosslayer),ncol = length(crosslayer),dimnames = list(names(crosslayer),names(crosslayer)))
  for (i in names(crosslayer) )
  {
    for (j in names(crosslayer[[i]]))
    {
      lnrel[i,j]=crosslayer[[i]][[j]]$lnrel
      SElnrel[i,j]=crosslayer[[i]][[j]]$selnrel
    }
  }
  rm(crosslayer)

  out=list(pooled_lnodds=pooled_lnodds,
           pooled_lnconf.int=pooled_lnconf.int,
           pooled_SElnodds=pooled_SElnodds,
           pooled_SElnnull=pooled_SElnnull,
           pooled_p=pooled_p,
           relative_lnodds=lnrel,
           relative_selnodds=SElnrel,
           pooled_rel_statistic=pooled_rel_statistic,
           pooled_rel_p=pooled_rel_p,
           results=results,
           param.record=list(response=response, group=group, strata=strata,
                             alpha=alpha,ties=ties,
                             nnt=nnt, verbose=verbose,upper=upper)
           )
  class(out)=c("Genodds",class(out))

  return(out)
}

print.Genodds<-function(x,...){

  args <- list(...)

  # Recover print options from x unless overridden in print command
  nnt <- x$param.record$nnt
  verbose <- x$param.record$verbose
  upper <- x$param.record$upper

  if("nnt" %in% names(args)){
    nnt <- args$nnt
  }
  if("verbose" %in% names(args)){
    verbose <- args$verbose
  }
  if("upper" %in% names(args)){
   upper <- args$upper
  }

  # It only makes sense to return either the upper or
  # lower triangle for relative risk ratios among strata.
  # This helper function will handle this
  print_triangle <- function(mat,upper=TRUE,diag=TRUE){
    blank <- sprintf("%s",paste(rep(" ",10),collapse=""))
    cat(blank)
    sapply(colnames(mat),function(x){cat(sprintf("%10s",x))})
    cat("\n")

    sapply(1:nrow(mat),function(i){
      cat(sprintf("%10s",rownames(mat)[i]))
      sapply(1:nrow(mat),function(j,i){
        if( (i<=j & upper & diag) | (i>=j & !upper & diag) |
            (i<j & upper & !diag) | (i>j & !upper & !diag)
        ){
          cat(sprintf("%10f",mat[i,j]))
        } else {
          cat(sprintf("%s",blank))
        }
      },i=i)
      cat("\n")
    })
  }

  # Header

  cat("\t Agresti's Generalized odds ratios\n\n")

  cat(sprintf("Odds that a random observation in group %s \nwill have a higher response score than a random\nobservation in group %s:\n\n", x$param.record$groupnames[2], x$param.record$groupnames[1]))

  #By layer results
  for (i in names(x$results))
  {
    if (length(x$results)>1)
    {
      cat(paste("  ",substr(paste(i,"        ",sep=""),1,10),sep=""))
    }
    else{
      cat("  ")
    }

    cat(sprintf("   Odds: %2.3f (%2.3f, %2.3f)      p=%1.4f",
                x$results[[i]]$odds,
                x$results[[i]]$conf.int[1],
                x$results[[i]]$conf.int[2],
                x$results[[i]]$p))

    if(nnt)
    {
      cat("\n")
      if (length(x$results)>1)
      {
        cat("            ")
      }
      else{
        cat("  ")
      }

      nntVals <- 1+2/(c(x$results[[i]]$odds,x$results[[i]]$conf.int)-1)

      if(sum(nntVals[2:3]>0)==1){ # Confidence intervals straddle 1

        if(nntVals[1]>0){nntVals[2:3] <- nntVals[3:2]}
        cat(sprintf("    NNT: %2.3f (-Inf, %2.3f)U(%2.3f,Inf)",
                    # ifelse(x$results[[i]]$odds>1,"B","H"),
                    nntVals[1],
                    min(nntVals[2:3]),
                    max(nntVals[2:3])
        ))


      } else {

        if(nntVals[1]>0){nntVals[2:3] <- nntVals[3:2]}
        cat(sprintf("    NNT: %2.3f (%2.3f, %2.3f)",
                    # ifelse(x$results[[i]]$odds>1,"B","H"),
                    nntVals[1],
                    nntVals[2],
                    nntVals[3]
        ))
      }

      cat("\n")




    }

    if (length(x$results)>1)
    {
      cat("\n")
    }
  }

  # Pooled results and comparison if relevant
  if (length(x$results)>1)
  {

    cat("\n")
    cat("Test of H0: odds ratios are equal among strata:\n")
    cat(sprintf("  X-squared = %2.2f, df= %d \t p=%1.4f", x$pooled_rel_statistic, length(x$results)-1, x$pooled_rel_p))
    cat("\n\n")

    if (x$pooled_rel_p>x$param.record$alpha | verbose==TRUE)
    {
      cat("Test of H0: pooled odds = 1:\n")
      cat(sprintf("  Pooled odds: %2.3f (%2.3f,%2.3f)", exp(x$pooled_lnodds), exp(x$pooled_lnconf.int[1]), exp(x$pooled_lnconf.int[2])))
      cat(sprintf("  p=%0.4f", x$pooled_p))

      if(nnt)
      {

        cat("\n")

        nntVals <- 1+2/(c(exp(x$pooled_lnodds),exp(x$pooled_lnconf.int))-1)

        if(sum(nntVals[2:3]>0)==1){ # Confidence intervals straddle 1

          if(nntVals[1]>0){nntVals[2:3] <- nntVals[3:2]}
          cat(sprintf("          NNT: %2.3f (-Inf, %2.3f)U(%2.3f,Inf)",
                      # ifelse(x$results[[i]]$odds>1,"B","H"),
                      nntVals[1],
                      min(nntVals[2:3]),
                      max(nntVals[2:3])
          ))


        } else {

          if(nntVals[1]>0){nntVals[2:3] <- nntVals[3:2]}
          cat(sprintf("          NNT: %2.3f (%2.3f, %2.3f)",
                      # ifelse(x$results[[i]]$odds>1,"B","H"),
                      nntVals[1],
                      nntVals[2],
                      nntVals[3]
          ))
        }
      }


    }

    if (x$pooled_rel_p<=x$param.record$alpha | verbose==TRUE)
    {
      cat("\n\nGeneralized relative risk ratios among strata:\n\n")

      print_triangle(exp(x$relative_lnodds),
                     upper=upper,diag=TRUE)

      cat("\n\nZ scores for relative risk ratios among strata:\n\n")

      print_triangle(x$relative_lnodds/x$relative_selnodds,
                     upper=upper,diag=TRUE)

    }

  }
  cat("\n")

  return(invisible(x))
}

