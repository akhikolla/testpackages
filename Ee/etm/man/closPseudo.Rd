\name{closPseudo}
\alias{closPseudo}

\title{
  Pseudo Value Regression for the Extra Length-of-Stay
}
\description{
  Pseudo Value Regression for the Extra Length-of-Stay
}
\usage{
closPseudo(data, state.names, tra, cens.name, s = 0,
           formula, na.action,
           aw = FALSE, ratio = FALSE,
           ncores = 1,
           trick_ties = FALSE)
}
\arguments{
  \item{data}{
    data.frame of the form data.frame(id,from,to,time)
    or (id,from,to,entry,exit)
    \describe{
      \item{id:}{patient id}
      \item{from:}{the state from where the transition occurs}
      \item{to:}{the state to which a transition occurs}
      \item{time:}{time when a transition occurs}
      \item{entry:}{entry time in a state}
      \item{exit:}{exit time from a state}
    }
  }
  \item{state.names}{A vector of characters giving the states names.}
  \item{tra}{A quadratic matrix of logical values describing the possible
    transitions within the multistate model.}
  \item{cens.name}{ A character giving the code for censored
    observations in the column 'to' of \code{data}. If there is no
    censored observations in your data, put 'NULL'.}
  \item{s}{Starting value for computing the transition probabilities.}
  \item{formula}{A formula with the covariates at the right of a
    \code{~} operator. Leave the left part empty.}
  \item{na.action}{A function which indicates what should happen when the data
    contain 'NA's.  The default is set by the 'na.action' setting
    of 'options', and is 'na.fail' if that is unset.  The
    'factory-fresh' default is 'na.omit'.}
    \item{aw}{Logical. Whether to compute the expected change of LOS using
    alternative weighting. Default is \code{FALSE}.}
  \item{ratio}{Logical. Compute the ratio of the expected length-of-stay
    given instermediate event status instead of a difference. Default
    value is \code{FALSE}}
  \item{ncores}{Number of cores used if doing parallel computation using
    the \pkg{parallel} package}
  \item{trick_ties}{If \code{TRUE}, pseudo values are computed only one
  per subject sharing the same entry, exit times / transition types.}
}
\details{
  The function calculates the pseudo-observations for the extra
  length-of-stay for each individual. These pseudo-observations can then
  be used to fit a direct regression model using generalized estimating
  equation (e.g., package \pkg{geepack}).

  Computation of the pseudo-observations can be parallelised using the
  \code{mclapply} function of the \pkg{parallel} package. See argument
  \code{ncores}.

  Recent versions of R have changed the \code{data.frame} function,
  where the default for the \code{stringsAsFactors}
  argument from \code{TRUE} to \code{FALSE}. \code{etm} currently
  depends on the states being factors, so that the user should use
  \code{data.frame(..., stringsAsFactors=TRUE)}.
}
\value{
  An object of class \code{closPseudo} with the following components:
  \item{pseudoData}{a data.frame containing \code{id}, computed pseudo
    values (see details) and the covariates as specified in the formula}
  \item{theta}{Estimates of excess LoS in the whole sample}
  \item{aw}{like in the function call}
  \item{call}{Function call}
}

\references{
  Andersen, P.K, Klein, J.P, Rosthoj, S. (2003). Generalised
  linear  models for correlated pseudo-observations,
  with applications to multi-state models. \emph{Biometrika},
  90(1):15--27.
}

\author{
  Arthur Allignol \email{arthur.allignol@gmail.com}
}

\seealso{
  \code{\link{mclapply}}, \code{\link[etm]{clos}}
}
\examples{
if(require("kmi", quietly = TRUE)) {

    ## data in kmi package
    data(icu.pneu)
    my.icu.pneu <- icu.pneu
    
    my.icu.pneu <- my.icu.pneu[order(my.icu.pneu$id, my.icu.pneu$start), ]
    masque <- diff(my.icu.pneu$id)
    
    my.icu.pneu$from <- 0
    my.icu.pneu$from[c(1, masque) == 0] <- 1
    
    my.icu.pneu$to2 <- my.icu.pneu$event
    my.icu.pneu$to2[my.icu.pneu$status == 0] <- "cens"
    my.icu.pneu$to2[c(masque, 1) == 0] <- 1
    
    
    my.icu.pneu$to <- ifelse(my.icu.pneu$to2 \%in\% c(2, 3), 2,
                             my.icu.pneu$to2)

    my.icu.pneu <- my.icu.pneu[, c("id", "start", "stop", "from", "to",
                                   "to2", "age", "sex")]
    names(my.icu.pneu)[c(2, 3)] <- c("entry", "exit")

    ## computation of the pseudo-observations
    \dontrun{
        ps.icu.pneu <- closPseudo(my.icu.pneu, c("0", "1", "2"), tra_ill(), "cens",
                                  formula = ~ sex + age)
        
        ## regression model using geepack
        require(geepack)
        fit <- geeglm(ps.e.phi ~ sex + age, id = id, data = ps.icu.pneu$pseudoData,
                      family = gaussian)
        
        summary(fit)
    }
} else {
    print("This example requires the kmi package")
}
}

\keyword{survival}
