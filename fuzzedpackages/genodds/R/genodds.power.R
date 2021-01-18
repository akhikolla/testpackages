#' Power Calculations for Generalized Odds Ratios
#' @description Provides power analysis for Agresti's Generalized Odds Ratios.
#'
#' @param p0 A numeric vector containing the probabilities in control group.
#' @param p1 A numeric vector containing the probabilities in treatment group.
#' @param N  A numeric vector containing total sample sizes.
#' @param power A numeric vector containing required  total sample size.
#' @param alpha Type 1 error.
#' @param ties A string specifying how ties should be treated.
#' @param w A numeric vector of length 2 specifying the relative weighting
#'          of sample size between treatment groups.
#' @param direction Direction for hypothesis test.
#'                  Must be one of \code{"two.sided"},\code{"upper.tail"} or \code{"lower.tail"}.
#'
#' @return
#' \describe{
#'     \item{If \code{power} is supplied}{A numeric vector containing required sample sizes to achieve specified powers.}
#'     \item{If \code{N} is supplied}{A numeric vector containing power at specified sample sizes.}
#' }
#'
#' @details
#'
#' See \code{\link{genodds}} for explanation of generalized odds ratios.
#'
#' \code{N} provides the total sample size.
#' Sample size per group can be calculated by \code{N*w/sum(w)}.
#'
#' When \code{power} is supplied, if no sufficient sample size is
#' found then this function will return \code{Inf}.
#'
#' @examples
#'
#' # Provide theoretical distributions of outcomes for each group
#' # Distributions taken from Lees et. al. (2010). See ?alteplase for a citation.
#' p0 <- c(0.224,0.191,0.082,0.133,0.136,0.043,0.191)
#' p1 <- c(0.109,0.199,0.109,0.120,0.194,0.070,0.200)
#'
#' # Calculate sample size required to achieve 80% and 90%
#' # power for these distributions
#' genodds.power(p0,p1,power=c(0.8,0.9))
#'
#' # genodds.power suggests a total sample size of 619 for 80% power.
#' # Round up to 620 for even sample size per group
#'
#' # Confirm these sample sizes lead to 80% and 90% power
#' genodds.power(p0,p1,N=c(620,830))
#'
#' @references
#' O'Brien, R. G., & Castelloe, J. (2006, March).
#' Exploiting the link between the Wilcoxon-Mann-Whitney test and a simple odds statistic.
#' In \emph{Thirty-first Annual SAS Users Group International Conference}.
#'
#' @export
genodds.power <- function(p0,p1,N=NULL,power=NULL,
                          alpha=0.05,ties="split",
                          w=c(0.5,0.5),
                          direction="two.sided")
{

  if(is.null(p0) | is.null(p1))
  {
    stop("p0 and p1 must not be empty")
  }


  if(length(p0) != length(p1))
  {
    stop("p0 and p1 must be the same length")
  }

  if(is.null(N) & !is.null(power))
  {
    if(prod((power <=1) & (power >=0))==0)
    {
      stop("Power must be in the interval [0,1]")
    }

    if(length(power)>1)
    {
      # If multiple power values were requested, calculate power
      # for each

      out <- sapply(power,function(p){
        genodds.power(p0=p0,
                      p1=p1,
                      N=NULL,
                      power=p,
                      alpha=alpha,
                      ties=ties,
                      w=w,
                      direction=direction)
        })

      return(out)
    }
    else
    {
      # Recursively call this function and use
      # Bisection method to find required sample size to achieve
      # specified power

      # Find initial boundaries for search by doubling the initial window
      # until the power at upper boundary is above specified power
      nLow <- 2
      nHigh <- 4
      while(genodds.power(p0,p1,nHigh,NULL,alpha,ties,w,direction)<power)
      {
        nLow <- nHigh
        nHigh <- nHigh*2

        if(is.infinite(nHigh))
        {
          # If we never find an appropriate boundary, return Inf
          return(Inf)
        }
      }

      root <- uniroot(function(x){
        genodds.power(p0,p1,x,NULL,alpha,ties,w,direction)-power
      },
      c(nLow,nHigh)
      )

      # Round required sample size up to nearest integer
      return(ceiling(root$root))
    }
  }
  else if(!is.null(N) & is.null(power))
  {
    if(prod(N>=2)==0)
    {
      stop("N must be at least 2")
    }

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

    # Plug values into formula to get power

    # Get ~pij = P(Y_i=j | group==i)
    p <- cbind(p0/sum(p0),p1/sum(p1))

    # This is p, the joing probability P(group==i,Y_i==j)
    wp <- p*cbind(rep(w[1],nrow(p)),rep(w[2],nrow(p)))

    Rt=wp[,2:1]
    Rs=get_Rs(wp)
    Rd=get_Rd(wp)

    # Redistribute ties
    if(!is.na(contr_fav))
    {
      Rs=Rs+(1-contr_fav)*Rt
      Rd=Rd+contr_fav*Rt
    }

    Pc=sum(wp*Rs)
    Pd=sum(wp*Rd)


    odds=Pc/Pd

    SEodds=2/Pd*(sum(wp*(odds*Rd-Rs)^2)/N)^0.5
    SElnodds=SEodds/odds


    # Primary effect size measure
    d <- log(odds)/(sqrt(N)*SElnodds)



    # Smooth p across groups and do this again.
    # This code is WET as hell, but whatever

    p=apply(p,1,mean)
    p=cbind(p,p)

    # standardise w to sum to 1
    w <- w/sum(w)

    wp <- p*cbind(rep(w[1],nrow(p)),rep(w[2],nrow(p)))

    Rt=wp[,2:1]
    Rs=get_Rs(wp)
    Rd=get_Rd(wp)

    # Redistribute ties
    if(!is.na(contr_fav))
    {
      Rs=Rs+(1-contr_fav)*Rt
      Rd=Rd+contr_fav*Rt
    }

    Pc=sum(wp*Rs)
    Pd=sum(wp*Rd)

    SEnull=2/Pd*(sum(wp*(1*Rd-Rs)^2)/N)^0.5
    SElnnull=SEnull/1

    if(direction=="two.sided")
    {
      nc <-N*d^2
      out <- pchisq((SElnnull/SElnodds)^2 * qchisq(1-alpha,1,0),
                    1,nc,
                    lower.tail = FALSE)
    }
    else if(direction=="lower.tail")
    {
      out <- pnorm(SElnnull/SElnodds * qnorm(alpha) - d*sqrt(N),lower.tail =T)
    }
    else if(direction=="upper.tail")
    {
      out <- pnorm(SElnnull/SElnodds * qnorm(1-alpha) - d*sqrt(N),lower.tail =F)
    }
    else
    {
      stop("incorrect direction")
    }

    return(out)
  }
  else
  {
    stop("must provide exactly one of N or power")
  }
}


