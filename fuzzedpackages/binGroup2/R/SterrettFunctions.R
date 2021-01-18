#####################################################################
# NAME:  Chris Bilder                                               #
# DATE:  3-13-09                                                    #
# UPDATE: 10-19-09 - made this a read-in file only and added MP and #
#         halving                                                   #
# Purpose: Summary measures for 1SIS, 2SIS, FIS, and Dorfman        #
#####################################################################
# Brianna Hitt - 03.06.2020
# Removed comparisons with halving and array testing 
# Removed associated num.perm, row.size, col.size, and seed.num arguments
# Removed associated information in the documentation details / examples

#' @title Summary measures for Sterrett algorithms
#'
#' @description Summary measures for Sterrett algorithms.
#'
#' @param p a vector of individual risk probabilities.
#' @param Sp the specificity of the diagnostic test.
#' @param Se the sensitivity of the diagnostic test.
#' @param plot logical; if TRUE, a plot of the informative Sterrett
#' CDFs will be displayed. Further details are given under 'Details'.
#' @param plot.cut.dorf logical; if TRUE, the cut-tree for Dorfman
#' testing will be displayed. Further details are given under 'Details'.
#' @param cond.prob.plot logical; if TRUE, a second axis for the
#' conditional probability plot will be displayed on the right side
#' of the plot.
#' @param font.name the name of the font to be used in plots.
#'
#' @details This function calculates summary measures for informative
#' Sterrett algorithms. Informative algorithms include one-stage informative
#' Sterrett (1SIS), two-stage informative Sterrett (2SIS), full informative
#' Sterrett (FIS), and Dorfman (two-stage hierarchical testing).
#'
#' The mean and standard deviation of the number of tests,
#' probability mass function (PMF), and cumulative distribution function (CDF)
#' are calculated for all informative Sterrett algorithms and
#' Dorfman testing. Conditional PMFs and conditional moments are calculated
#' for all informative Sterrett algorithms. Subtracting the mean number of
#' tests for two procedures gives the area difference between their CDFs.
#' This area difference is calculated for each pairwise comparison of
#' 1SIS, 2SIS, FIS, and Dorfman testing. CDF plots provide a visualization
#' of how probabilities are distributed over the number of tests. CDFs that
#' increase more rapidly to 1 correspond to more efficient retesting
#' procedures.
#'
#' Non-informative Sterrett (NIS) decodes positive groups by retesting
#' individuals at random, so there are \eqn{I!} different possible
#' NIS implementations. CDFs are found by permuting the elements in the
#' vector of individual risk probabilities and using the FIS CDF expression
#' without reordering the individual probabilities. That is, the FIS
#' procedure uses the most efficient NIS implementation, which is to
#' retest individuals in order of descending probabilities.
#' When implementing the informative Sterrett algorithms with a large
#' number of individuals, an algorithm is used to compute the PMF
#' for the number of tests under FIS. This is done automatically
#' by \kbd{Sterrett} for \eqn{I}>12. The algorithm is described in
#' detail in the Appendix of Bilder et al. (2010).
#'
#' @return A list containing:
#' \item{mean.sd}{a data frame containing the mean and standard deviation
#' of the expected number of tests for one-stage informative Sterrett (1SIS),
#' two-stage informative Sterrett (2SIS), full informative Sterrett (FIS),
#' and Dorfman testing.}
#' \item{PMF}{a data frame containing the probability mass function
#' for the number of tests possible for one-stage informative Sterrett (1SIS),
#' two-stage informative Sterrett (2SIS), full informative Sterrett (FIS),
#' and Dorfman testing.}
#' \item{CDF}{a data frame containing the cumulative distribution function
#' for the number of tests possible for one-stage informative Sterrett (1SIS),
#' two-stage informative Sterrett (2SIS), full informative Sterrett (FIS),
#' and Dorfman testing.}
#' \item{cond.PMF}{a data frame containing the conditional probability
#' mass function for the number of tests possible for one-stage informative
#' Sterrett (1SIS), two-stage informative Sterrett (2SIS), and full
#' informative Sterrett (FIS) testing.}
#' \item{cond.moments}{a data frame containing the mean and standard deviation
#' of the conditional moments for one-stage informative
#' Sterrett (1SIS), two-stage informative Sterrett (2SIS), and full
#' informative Sterrett (FIS) testing.}
#' \item{save.diff.CDF}{a data frame containing the sum of the differences
#' in the cumulative distribution function for each pairwise comparison of
#' one-stage informative Sterrett (1SIS), two-stage informative Sterrett (2SIS),
#' full informative Sterrett (FIS), and Dorfman testing.}
#'
#' @author This function was originally written as \kbd{info.gt}
#' by Christopher Bilder for Bilder et al. (2010). The function was obtained 
#' from \url{http://chrisbilder.com/grouptesting}. Minor modifications were 
#' made for inclusion of the function in the \code{binGroup2} package.
#'
#' @references
#' \insertRef{Bilder2010a}{binGroup2}
#'
#' @seealso
#' \code{\link{expectOrderBeta}} for generating a vector of individual risk 
#' probabilities for informative group testing and \code{\link{opChar1}} for 
#' calculating operating characteristics with hierarchical and array-based 
#' group testing algorithms.
#'
#' @family operating characteristic functions
#'
#' @examples
#' # Example 1: FIS provides the smallest mean
#' #   number of tests and the smallest standard
#' #   deviation. 2SIS has slightly larger mean
#' #   and standard deviation than FIS, but
#' #   its performance is comparable, indicating
#' #   2SIS may be preferred because it is
#' #   easier to implement.
#' set.seed(1231)
#' p.vec1 <- rbeta(n=8, shape1=1, shape2=10)
#' save.it1 <- Sterrett(p=p.vec1, Sp=0.90, Se=0.95)
#' save.it1$mean.sd
#'
#' # Example 2: One individual is "high risk" and
#' #   the others are "low risk". Since there is
#' #   only one high-risk individual, the three
#' #   informative Sterrett procedures perform
#' #   similarly. All three informative Sterrett
#' #   procedures offer large improvements over
#' #   Dorfman testing.
#' p.vec2 <- c(rep(x=0.01, times=9), 0.5)
#' save.it2 <- Sterrett(p=p.vec2, Sp=0.99, Se=0.99,
#'                      cond.prob.plot=TRUE)
#' save.it2$mean.sd
#'
#' # Example 3: Two individuals are at higher
#' #   risk than the others. All three informative
#' #   Sterrett procedures provide large
#' #   improvements over Dorfman testing.
#' # Due to the large initial group size, an
#' #   algorithm (described in the Appendix of
#' #   Bilder et al. (2010)) is used for FIS.
#' #   The Sterrett() function does this
#' #   automatically for I>12.
#' p.vec3 <- c(rep(x=0.01, times=98), 0.1, 0.1)
#' save.it3 <- Sterrett(p=p.vec3, Sp=0.99, Se=0.99)
#' save.it3$mean.sd

# Brianna Hitt - 04.02.2020
# Changed cat() to warning()

Sterrett <- function(p, Sp, Se, plot=TRUE, plot.cut.dorf = TRUE, 
                     cond.prob.plot = FALSE, font.name = "sans") {

  #####################################################################
  # Initial calculations

  I<-length(p)
  p<-sort(p)   #In case p's are not in order


  #####################################################################
  # 1SIS

  tests<-c(1, 3:I, I+1, I+2)
  prob<-numeric(length(tests))

  #t=1
  prob[1]<-Sp*prod(1-p) + (1-Se)*(1 - prod(1-p))

  #t = 3, 4, ..., I
  for (t in 3:I) {
    if(t==3) {
      prob[t-1]<-(1-Sp)*Sp^(t-2)*prod(1-p)*(1-Sp-Se) + Se *
        ((1-Sp)*(1-p[I+3-t]) + Se*p[I+3-t]) * (Sp*prod(1-p[1:(I+2-t)]) + (1-Se)*(1-prod(1-p[1:(I+2-t)])))
    }
    else {
      prob[t-1]<-(1-Sp)*Sp^(t-2)*prod(1-p)*(1-Sp-Se) + Se * prod(Sp*(1-p[(I+4-t):I]) + (1-Se)*p[(I+4-t):I]) *
        ((1-Sp)*(1-p[I+3-t]) + Se*p[I+3-t]) * (Sp*prod(1-p[1:(I+2-t)]) + (1-Se)*(1-prod(1-p[1:(I+2-t)])))
    }
  }

  #t = I + 1
  prob[length(tests)-1]<-(1-Sp) * Sp^(I-2) * prod(1-p) + Se*( prod( Sp*(1-p[3:I]) + (1-Se)*p[3:I]) - Sp^(I-2)*prod(1-p))

  #t = I + 2
  last.set<-numeric(length(4:(I+1)))
  for (m in 4:(I+1)) {
    if(m==I+1) {
      last.set[m-3]<-(1-Sp-Se)*(1-Sp)^2*Sp^(I+1-m)*prod(1-p) + Se *
        ((1-Sp)*(1-p[m-1]) + Se*p[m-1]) * ((1-Sp)*prod(1-p[1:(m-2)]) + Se*(1-prod(1-p[1:(m-2)]))) #First term for t=I+2
    }
    else {
      last.set[m-3]<-(1-Sp-Se)*(1-Sp)^2*Sp^(I+1-m)*prod(1-p) + Se * prod(Sp*(1-p[m:I]) + (1-Se)*p[m:I]) *
        ((1-Sp)*(1-p[m-1]) + Se*p[m-1]) * ((1-Sp)*prod(1-p[1:(m-2)]) + Se*(1-prod(1-p[1:(m-2)]))) #Other terms for t=I+2 (allows for G_ = 0 in it)
    }
  }
  prob[length(tests)]<-sum(last.set)

  mean.1stage.IS<-sum(tests*prob)
  sd.1stage.IS<-sqrt(sum(tests^2*prob) - mean.1stage.IS^2)



  #####################################################################
  # 2SIS - Additional calculations needed for 2SIS and combine with 1SIS calculations at the end

  tests5<-c(1, 3:I, I+1, I+2, I+3)

  #Can not do 2SIS when I < 4
  if (I>3) {

    prob5<-numeric(length(tests5))

    #t=5, ..., I+1 - only the new 2-stage part
    part2.5<-matrix(data = 0, nrow = I-3, ncol = I-3)

    for (t in 5:(I+1)) {
      for (m in 1:(t-4)) {

        if ((I+2-m)>I && (I+5-t)>(I-m)) {
          part2.5[t-4,m]<-(1-Sp-Se)*(1-Sp)^3*Sp^(t-4)*prod(1-p) +
            Se *                                                  ((1-Sp)*(1-p[I-m+1]) + Se*p[I-m+1])   *
            ( (1-Sp-Se)*(1-Sp)*Sp^(t-m-3)*prod(1-p[1:(I-m)]) + Se *
                ((1-Sp)*(1-p[I+4-t]) + Se*p[I+4-t]) * (Sp*prod(1-p[1:(I+3-t)]) + (1-Se)*(1-prod(1-p[1:(I+3-t)]))))
        }

        if ((I+2-m)<=I && (I+5-t)>(I-m)) {
          part2.5[t-4,m]<-(1-Sp-Se)*(1-Sp)^3*Sp^(t-4)*prod(1-p) +
            Se * prod(Sp*(1-p[(I+2-m):I]) + (1-Se)*p[(I+2-m):I]) * ((1-Sp)*(1-p[I-m+1]) + Se*p[I-m+1])   *
            ( (1-Sp-Se)*(1-Sp)*Sp^(t-m-3)*prod(1-p[1:(I-m)]) + Se *
                ((1-Sp)*(1-p[I+4-t]) + Se*p[I+4-t]) * (Sp*prod(1-p[1:(I+3-t)]) + (1-Se)*(1-prod(1-p[1:(I+3-t)]))))
        }

        if ((I+2-m)>I && (I+5-t)<=(I-m)) {
          part2.5[t-4,m]<-(1-Sp-Se)*(1-Sp)^3*Sp^(t-4)*prod(1-p) +
            Se *                                                  ((1-Sp)*(1-p[I-m+1]) + Se*p[I-m+1])   *
            ( (1-Sp-Se)*(1-Sp)*Sp^(t-m-3)*prod(1-p[1:(I-m)]) + Se * prod(Sp*(1-p[(I+5-t):(I-m)]) + (1-Se)*p[(I+5-t):(I-m)]) *
                ((1-Sp)*(1-p[I+4-t]) + Se*p[I+4-t]) * (Sp*prod(1-p[1:(I+3-t)]) + (1-Se)*(1-prod(1-p[1:(I+3-t)]))))
        }

        if ((I+2-m)<=I && (I+5-t)<=(I-m)) {
          part2.5[t-4,m]<-(1-Sp-Se)*(1-Sp)^3*Sp^(t-4)*prod(1-p) +
            Se * prod(Sp*(1-p[(I+2-m):I]) + (1-Se)*p[(I+2-m):I]) * ((1-Sp)*(1-p[I-m+1]) + Se*p[I-m+1])   *
            ( (1-Sp-Se)*(1-Sp)*Sp^(t-m-3)*prod(1-p[1:(I-m)]) + Se * prod(Sp*(1-p[(I+5-t):(I-m)]) + (1-Se)*p[(I+5-t):(I-m)]) *
                ((1-Sp)*(1-p[I+4-t]) + Se*p[I+4-t]) * (Sp*prod(1-p[1:(I+3-t)]) + (1-Se)*(1-prod(1-p[1:(I+3-t)]))))
        }
      }
    }
    prob5[4:(4+(I-3)-1)]<-rowSums(part2.5) #Each row of the matrix is for a particular t


    #t=I+2 - only the new 2-stage part
    part2.I.2<-numeric(I-3)
    for (m in 5:(I+1)) {

      if (m>I) {
        part2.I.2[m-4]<-(1-Sp-Se)*(1-Sp)^2*Sp^(I-3)*prod(1-p) +
          Se *                                       ((1-Sp)*(1-p[m-1]) + Se*p[m-1]) *
          ((1-Sp-Se)*Sp^(m-4)*prod(1-p[1:(m-2)]) + Se*(prod( Sp*(1-p[3:(m-2)]) + (1-Se)*p[3:(m-2)])))
      }

      if (m<=I) {
        part2.I.2[m-4]<-(1-Sp-Se)*(1-Sp)^2*Sp^(I-3)*prod(1-p) +
          Se * prod( Sp*(1-p[m:I]) + (1-Se)*p[m:I]) * ((1-Sp)*(1-p[m-1]) + Se*p[m-1]) *
          ((1-Sp-Se)*Sp^(m-4)*prod(1-p[1:(m-2)]) + Se*(prod( Sp*(1-p[3:(m-2)]) + (1-Se)*p[3:(m-2)])))
      }

    }
    prob5[length(tests5)-1]<-sum(part2.I.2)


    #t=I+3
    save.t3<-matrix(data = 0, nrow = I-3, ncol = I-3)
    for (m in 1:(I-3)) {
      for (r in 1:(I-m-2)) {

        #Have both sets of G_=0
        if((I+2-m)<=I && (r+3)<=(I-m)) {
          save.t3[m,r]<-(1-Sp-Se)*(1-Sp)^4*Sp^(I-r-3)*prod(1-p) +
            Se *  prod(Sp*(1-p[(I+2-m):I]) + (1-Se)*p[(I+2-m):I]) * ((1-Sp)*(1-p[I-m+1]) + Se*p[I-m+1]) *
            ((1-Sp-Se)*(1-Sp)^2*Sp^(I-m-r-2)*prod(1-p[1:(I-m)]) +
               Se * prod(Sp*(1-p[(r+3):(I-m)]) + (1-Se)*p[(r+3):(I-m)])*((1-Sp)*(1-p[r+2]) + Se*p[r+2]) *
               ((1-Sp)*prod(1-p[1:(r+1)]) + Se*(1-prod(1-p[1:(r+1)]))))
        }

        #Have NO sets of G_=0
        if((I+2-m)>I && (r+3)>(I-m)) {
          save.t3[m,r]<-(1-Sp-Se)*(1-Sp)^4*Sp^(I-r-3)*prod(1-p) +
            Se *                                                     ((1-Sp)*(1-p[I-m+1]) + Se*p[I-m+1]) *
            ((1-Sp-Se)*(1-Sp)^2*Sp^(I-m-r-2)*prod(1-p[1:(I-m)]) +
               Se *                                                           ((1-Sp)*(1-p[r+2]) + Se*p[r+2]) *
               ((1-Sp)*prod(1-p[1:(r+1)]) + Se*(1-prod(1-p[1:(r+1)]))))
        }

        #Have 1st set of G_=0 only
        if((I+2-m)<=I && (r+3)>(I-m)) {
          save.t3[m,r]<-(1-Sp-Se)*(1-Sp)^4*Sp^(I-r-3)*prod(1-p) +
            Se *   prod(Sp*(1-p[(I+2-m):I]) + (1-Se)*p[(I+2-m):I]) * ((1-Sp)*(1-p[I-m+1]) + Se*p[I-m+1]) *
            ((1-Sp-Se)*(1-Sp)^2*Sp^(I-m-r-2)*prod(1-p[1:(I-m)]) +
               Se *                                                      ((1-Sp)*(1-p[r+2]) + Se*p[r+2]) *
               ((1-Sp)*prod(1-p[1:(r+1)]) + Se*(1-prod(1-p[1:(r+1)]))))
        }

        #Have 2nd set of G_=0 only
        if((I+2-m)>I && (r+3)<=(I-m)) {
          save.t3[m,r]<-(1-Sp-Se)*(1-Sp)^4*Sp^(I-r-3)*prod(1-p) +
            Se *                                                        ((1-Sp)*(1-p[I-m+1]) + Se*p[I-m+1]) *
            ((1-Sp-Se)*(1-Sp)^2*Sp^(I-m-r-2)*prod(1-p[1:(I-m)]) +
               Se*prod(Sp*(1-p[(r+3):(I-m)]) + (1-Se)*p[(r+3):(I-m)])*((1-Sp)*(1-p[r+2]) + Se*p[r+2]) *
               ((1-Sp)*prod(1-p[1:(r+1)]) + Se*(1-prod(1-p[1:(r+1)]))))
        }
      }
    }
    prob5[length(tests5)]<-sum(save.t3)


    #To combine prob5 with 1SIS probabilities, we need to subtract off the 1-stage probabilities that are broken down further
    #  into parts (these correspond to the 1SIS probabilities with G_1..k = 1 where k>2
    IS1.sub<-prob
    IS1.sub[I+2-1]<-last.set[1] #Keep first probability for t = I+2 only for IS because this contains the G_12=1
    IS2<-c(IS1.sub,0) + prob5   #Added a 0 to IS1.sub to acccount for the t=I+3 term now in IS2
    #IS is the 2-stage PMF

    mean.2stage.IS<-sum(tests5*IS2)
    sd.2stage.IS<-sqrt(sum(tests5^2*IS2) - mean.2stage.IS^2)
  } else

  {
    warning("2SIS can not be done for I = 3")
    mean.2stage.IS<-NA
    sd.2stage.IS<-NA
    IS2<-rep(x = NA, times = length(tests)+1)
  }



  #####################################################################
  # Dorfman

  tests3<-c(1,I+1)
  prob3<-numeric(2)
  prob3[1]<-Sp*prod(1-p) + (1-Se)*(1 - prod(1-p))
  prob3[2]<-(1-Sp)*prod(1-p) + Se*(1 - prod(1-p))

  mean.dorf<-sum(tests3*prob3)
  sd.dorf<-sqrt(sum(tests3^2*prob3) - mean.dorf^2)


  #####################################################################
  # FIS - Memory problems occur for I > 12 with a full matrix approach.
  #       Peng Chen's program is run for I > 12 to get around this problem.

  # Brianna Hitt - 04.02.2020
  # changed $T from err_is_pmf() to $Tests
  
  if (I>12) {

    Se.original<-Se
    Sp.original<-Sp
    if (Se == 1) {Se<-0.99999999999999}
    if (Sp == 1) {Sp<-0.99999999999999}
    if (Se == 1 | Sp == 1) { 
      warning("If Se or Sp = 1, they are replaced with 0.99999999999999 for ",
              "FIS due to the algorithm used to find the PMF.") 
      }

    full.IS<-err_is_pmf(p,Se,Sp) #From Peng's program
    t.full.IS<-full.IS$Tests
    pmf.full.IS<-full.IS$prob

    mean.full.IS<-sum(t.full.IS*full.IS$prob)
    sd.full.IS<-sqrt(sum(t.full.IS^2*full.IS$prob) - mean.full.IS^2)

    if (Se.original == 1) {Se<-Se.original}
    if (Sp.original == 1) {Sp<-Sp.original}

  }

  else {

    #This is a way to find the PMF using matrices without paying attention to how large these matrices can be.  In the end, memory
    # problems will occur for a group size not small.  The matrix set up below is discussed somewhat at the end of the appendix (a little different version of A).
    # One could just remove the else part here and use the above I>12 for I<=12 instead.
    Se.original<-Se
    Sp.original<-Sp
    if (Se == 1) {Se<-0.99999999999999}
    if (Sp == 1) {Sp<-0.99999999999999}

    if (Se == 1 | Sp == 1) { 
      warning("If Se or Sp = 1, they are replaced with 0.99999999999999 for ",
              "FIS due to the algorithm used to find the PMF.") 
      }

    j<-function(length) { matrix(data = 1, nrow = length, ncol = 1) }
    q.vec<-kronecker(matrix(data = c(1-p[2], p[2]), nrow = 2, ncol = 1), c(1-p[1], p[1])) #q^(j) in the appendix
    A.mat<-matrix(data = c(  Sp,  (1-Sp)*Sp^2,      (1-Sp)^2*Sp,      (1-Sp)^2*Sp,    (1-Sp)^3,
                             1-Se, Se*Sp*(1-Se),          Se^2*Sp, Se*(1-Se)*(1-Sp), Se^2*(1-Sp),
                             1-Se, Se*(1-Se)*Sp, Se*(1-Se)*(1-Sp),          Se^2*Sp, Se^2*(1-Sp),
                             1-Se,  Se*(1-Se)^2,      Se^2*(1-Se),      Se^2*(1-Se),       Se^3), nrow = 4, ncol = 5, byrow = TRUE) #A matrix in the appendix
    t.vec<-c(1, rep(3, times = 4))

    for (I.iter in (3:I)) {
      A11<-kronecker(j(2^(I.iter-1)), matrix(data = c(1, rep(Sp, times = 5*2^(I.iter-3)-1)), nrow = 1, ncol = 5*2^(I.iter-3)))
      A12<-(1-Sp)*kronecker(t(j(5*2^(I.iter-3))), matrix(data = c(1-Sp, rep(Se, times = 2^(I.iter-1)-1)), nrow = 2^(I.iter-1), ncol = 1))
      A21<-cbind(c((1-Se)/Sp,j(2^(I.iter-1)-1)), kronecker(t(j(5*2^(I.iter-3)-1)), matrix(data = c(Se*(1-Se)/(1-Sp), rep(1-Se, times = 2^(I.iter-1)-1)), nrow = 2^(I.iter-1), ncol = 1)))
      A.mat<-rbind(cbind(A11*A.mat, A12*A.mat), cbind(A21*A.mat, Se^2*A.mat))
      q.vec<-kronecker(c(1-p[I.iter], p[I.iter]), q.vec)
      t.vec<-c(1, 1 + t.vec[-1], 2 + t.vec)
    }

    mean.full.IS<-sum(t(q.vec)%*%A.mat%*%t.vec)
    sd.full.IS<-sqrt(sum(t(q.vec)%*%A.mat%*%t.vec^2) - mean.full.IS^2)

    #PMF
    # Brianna Hitt - 04.02.2020
    # changed "T" in pmf1 to "Tests"
    pmf1<-data.frame(Tests = t.vec, prob = t(t(q.vec)%*%A.mat))  #S vector in the Appendix
    pmf.full.IS<-tapply(X = pmf1$prob, INDEX = pmf1$Tests, FUN = sum) #condensed over the duplicate values of t
    t.full.IS<-as.numeric(names(pmf.full.IS))

    if (Se.original == 1) {Se<-Se.original}
    if (Sp.original == 1) {Sp<-Sp.original}

  }



  ######################################################################
  # CDF plots

  if (plot == TRUE) {

    orig_par <- par(no.readonly = TRUE)
    on.exit(par(orig_par))
    
    #New graphics window
    # code below is obsolete - change to dev.new() for inclusion in binGroup2
    # win.graph(width = 10, height = 7, pointsize = 12)
    dev.new()
    par(family = font.name)
    par(lend = "square") #Line end style - default is rounded

    marg_par <- par(mar = c(5.1, 4.1, 4.1, 5.1))  #Added to help with right side y-axis

    #maximum value for x-axis
    max.x<-max(t.full.IS, I+4, na.rm = TRUE)

    #1SIS
    Ik<-I
    plot(x = tests, y = cumsum(prob), type = "s", xlab = expression(italic(t[k])), ylab = expression(italic(P)(italic(T[k]) <=
                                                                                                                 italic(t[k]))), lwd = 5, col = "darkgreen", panel.first=grid(col = "gray"),
         main = substitute(paste("CDFs with ", italic(I[k])==Ik, ", ", italic(S[e])==Se, ", and ", italic(S[p])==Sp), list(Ik=I, Se=Se, Sp=Sp)),
         xlim = c(0, max.x), ylim = c(0,1))
    abline(h = 0)
    segments(x0 = 0, y0 = 0, x1 = 1, y1 = 0, lwd = 5, col = "darkgreen")
    segments(x0 = 1, y0 = 0, x1 = 1, y1 = prob[1], lwd = 5, col = "darkgreen")
    segments(x0 = I+2, y0 = 1, x1 = max.x+1, y1 = 1, lwd = 5, col = "darkgreen")

    #2SIS
    if(I>=4) {
      lines(x = tests5, y = cumsum(IS2), type = "s", col = "purple", lty = 5, lwd = 3)
      segments(x0 = 0, y0 = 0, x1 = 1, y1 = 0, lwd = 3, col = "purple", lty = 5)
      segments(x0 = 1, y0 = 0, x1 = 1, y1 = IS2[1], lwd = 3, col = "purple", lty = 5)
      segments(x0 = I+3, y0 = 1, x1 = max.x+1, y1 = 1, lwd = 3, col = "purple", lty = 5)
    }
    else text(x=1, y=0.15, labels = "Warning: 2-stage IS can not be done with I = 3", cex = 0.75, pos = 4)


    #FIS
    lines(x = t.full.IS, y = cumsum(pmf.full.IS), type = "s", col = "blue", lty = 6, lwd = 2)
    segments(x0 = 0, y0 = 0, x1 = 1, y1 = 0, lwd = 2, col = "blue", lty = 6)
    segments(x0 = 1, y0 = 0, x1 = 1, y1 = pmf.full.IS[1], lwd = 2, col = "blue", lty = 6)
    segments(x0 = max(t.full.IS), y0 = 1, x1 = max.x+1, y1 = 1, lwd = 2, col = "blue", lty = 6)


    #Since plot can get crowded, can specify whether or not to plot cut-tree and Dorfman
    if (plot.cut.dorf == TRUE) {

      #Dorfman
      lines(x = tests3, y = cumsum(prob3), type = "s", col = "red", lwd = 2, lty = 4)
      segments(x0 = 0, y0 = 0, x1 = 1, y1 = 0, lwd = 2, col = "red", lty = 4)
      segments(x0 = 1, y0 = 0, x1 = 1, y1 = prob3[1], lwd = 2, col = "red", lty = 4)
      segments(x0 = I+1, y0 = 1, x1 = max.x+1, y1 = 1, lwd = 2, col = "red", lty = 4)
      legend(x=max.x*3/4, y=0.4, legend=c("1SIS", "2SIS", "FIS", "Dorfman"), lty=c(1, 5, 6,  4),
             col=c("darkgreen", "purple", "blue", "red"), bty="y", lwd=c(5,3,2,2), cex=1, bg = "white")

    } else
    {
      legend(x=max.x*3/4, y=0.4, legend=c("1SIS", "2SIS", "FIS"), lty=c(1, 5, 6),
             col=c("darkgreen", "purple", "cyan"), bty="n", lwd=c(5,3,2), cex=1)
    }

    p.vec<-as.character(round(p[1],3))
    for (i in 2:I) {
      p.vec<-paste(p.vec, ", ", as.character(round(p[i],3)), sep="")
    }

    if (I<=10) {
      text(x=1, y=0.1, labels = substitute(paste(italic(bolditalic(p)[k]), " = [", p.vec, "]"), list(p.vec=p.vec)), cex = 1, pos = 4)
    } else
    {
      #text(x=1, y=0.1, labels = expression(paste(italic(bolditalic(p)[k]), " = [0.01, ..., 0.01, 0.1, 0.1]")),cex = 1, pos = 4)
      text(x=1, y=0.1, labels = substitute(paste(italic(bolditalic(p)[k]), " = [", p1, ", ", p2, ", ... ,", p20, "]"),
                                           list(p1=round(p[1],4), p2=round(p[2],4), p20=round(p[length(p)],4))), cex = 1, pos = 4)
    }

    par(marg_par)

  }

  if (cond.prob.plot == TRUE) {

    #2nd y-axis for the P(T_k = t_k | Z_k = 1)
    prob.at.least1<-prob[1]
    space<-(1-prob.at.least1)/5 #space equally out from the beginning to the end.
    axis(side = 4, at = seq(from = prob.at.least1, to = 1, by = space), labels = c(0.0,0.2,0.4,0.6,0.8,1.0), line=-2)
    mtext(text = expression(paste(italic(P)(italic(T[k]) == italic(t[k]) ~ "|" ~ italic(Z[k]) == 1))), side = 4, at=0.75, line = 0.5)
    abline(h = prob.at.least1)
  }



  #####################################################################
  #Comparisons

  #Moments
  mean.sd<-data.frame(Method = c("1SIS", "2SIS", "FIS", "Dorfman"),
                      mean = c(mean.1stage.IS, mean.2stage.IS, mean.full.IS, mean.dorf),
                      sd = c(sd.1stage.IS, sd.2stage.IS, sd.full.IS, sd.dorf))

  #PMFs
  set1<-data.frame(tests, one.stage.IS.PMF = prob)
  set1.1<-data.frame(tests=tests5, two.stage.IS.PMF = IS2)
  set1.2<-data.frame(tests=t.full.IS, full.IS.PMF = pmf.full.IS)
  set3<-data.frame(tests=tests3, Dorfman.PMF = prob3)
  set1.set1.1<-merge(set1, set1.1, all = TRUE)
  set1.set1.12<-merge(set1.set1.1, set1.2, all = TRUE)
  PMF<-merge(set1.set1.12, set3, all = TRUE)

  #CDFs
  set1<-data.frame(tests, one.stage.IS.PMF = cumsum(prob))
  set1.1<-data.frame(tests=tests5, two.stage.IS.PMF = cumsum(IS2))
  set1.2<-data.frame(tests=t.full.IS, full.IS.PMF = cumsum(pmf.full.IS))
  set3<-data.frame(tests=tests3, Dorfman.PMF = cumsum(prob3))
  set1.set1.1<-merge(set1, set1.1, all = TRUE)
  set1.set1.12<-merge(set1.set1.1, set1.2, all = TRUE)
  CDF<-merge(set1.set1.12, set3, all = TRUE)


  #Fill in the NA's in the CDFs
  fill.in.CDF<-function(CDF.name) {
    temp<-is.na(CDF.name)
    CDF.new<-numeric(length(CDF.name))
    for (k in 1:(length(CDF.name)-1)) {
      if (temp[k] == FALSE)
      {CDF.new[k]<-CDF.name[k]}
      else {
        CDF.new[k]<-max(CDF.name[1:k], na.rm = TRUE)
      }
    }
    CDF.new[length(CDF.name)]<-1
    CDF.new
  }

  #Find the differences in the CDFs
  save.diff<-matrix(data = NA, nrow = 3, ncol = 3)
  for (i in 2:(ncol(CDF)-1)) {
    for (j in (i+1):ncol(CDF)) {
      CDF1<-fill.in.CDF(CDF[,i])
      CDF2<-fill.in.CDF(CDF[,j])
      diff.CDF<-CDF1-CDF2
      save.diff[i-1,j-2]<-sum(diff.CDF, na.rm = TRUE)
    }
  }
  row.names(save.diff)<-c("1SIS", "2SIS", "FIS")
  save.diff<-as.data.frame(save.diff)
  names(save.diff)<-c("2SIS", "FIS", "Dorfman")


  #Conditional comparisons
  cond.PMF<-data.frame(tests = PMF$tests[-1],
                       one.stage.IS = PMF$one.stage.IS.PMF[-1]/(1-PMF$one.stage.IS.PMF[1]),
                       two.stage.IS = PMF$two.stage.IS.PMF[-1]/(1-PMF$two.stage.IS.PMF[1]),
                       full.IS = PMF$full.IS.PMF[-1]/(1-PMF$full.IS.PMF[1]))
  cond.mean<-c(sum(cond.PMF$one.stage.IS*cond.PMF$tests, na.rm = TRUE),
               sum(cond.PMF$two.stage.IS*cond.PMF$tests, na.rm = TRUE),
               sum(cond.PMF$full.IS*cond.PMF$tests, na.rm = TRUE))
  cond.var<-c(sum((cond.PMF$tests-cond.mean[1])^2*cond.PMF$one.stage.IS, na.rm = TRUE),
              sum((cond.PMF$tests-cond.mean[2])^2*cond.PMF$two.stage.IS, na.rm = TRUE),
              sum((cond.PMF$tests-cond.mean[3])^2*cond.PMF$full.IS, na.rm = TRUE))
  cond.moments<-data.frame(Method = c("1SIS", "2SIS", "FIS"), mean = cond.mean, sd = sqrt(cond.var))


  # Brianna Hitt - 03.06.2020
  # Comparisons to non-informative procedures of matrix pooling / 
  #   array testing and halving removed
  # #####################################################################
  # #Non-informative procedures added on 10-19-09
  # 
  # set.seed(seed.num)
  # 
  # ###############
  # #Matrix pooling
  # exp.mp<-NULL
  # if (!is.null(row.size) & !is.null(col.size)) {
  #   if (I == row.size*col.size) {
  #     keep.MP<-numeric(num.perm)
  #     for(i in 1:num.perm) {
  #       keep.MP[i]<-mat.pool(p = sample(p), Se = Se, Sp = Sp, row.size = row.size, col.size = col.size)
  #     }
  #     exp.mp<-mean(keep.MP)
  #   }
  #   else warning("(row size)*(column size) does not equal group size; no array testing calculations are performed \n")
  # }
  # 
  # if (is.null(row.size) | is.null(col.size)) {
  #   warning("No row and/or column sizes given for array testing \n") #I had difficulty putting this in the if else statements above
  # }
  # 
  # ###############
  # #Halving
  # 
  # if (I <= 8) {
  #   warning("Note that 3H requires I > 2, 4H requires I > 4, and 5H requires I > 8 \n")
  # }
  # 
  # keep.halving<-matrix(data = NA, nrow = num.perm, ncol = 3)
  # for(i in 1:num.perm) {
  #   keep.halving1<-halving.comp(p = sample(p), Se = Se, Sp = Sp)
  #   keep.halving[i,1:length(keep.halving1)]<-keep.halving1 #Put results in this way in case there is a NA due to larger number of stages than possible for a halving procedure
  # }
  # 
  # halving.mean<-colMeans(keep.halving)
  # halving.mean2<-data.frame(halving.3H = halving.mean[1], halving.4H = halving.mean[2], halving.5H = halving.mean[3])
  # 

  #####################################################################
  #Return values

  # Brianna Hitt - 03.06.2020
  # Removed results for matrix pooling (MP = exp.mp) and 
  #   halving (halving = halving.mean2)
  # list(mean.sd = mean.sd, PMF = PMF, CDF = CDF, cond.PMF = cond.PMF, cond.moments = cond.moments, save.diff.CDF = save.diff,
  #      MP = exp.mp, halving = halving.mean2)
  
  list(mean.sd = mean.sd, PMF = PMF, CDF = CDF, cond.PMF = cond.PMF, 
       cond.moments = cond.moments, save.diff.CDF = save.diff)
  
}




# Start support functions for Sterrett()
##########################################################################

############################################################
#This function recursively produces the
#PMF of the number of retests T (excluding the initial test
#for the master pool) for the informative retesting procedure
#with testing errors
#Last modifed date:11-12-2008
#Author:Peng Chen
###########################################################
#3-13-09: Program was modified by Chris Bilder to match the notation
#         of the paper.

# Brianna Hitt - 04.02.2020
# changed the "T" object to "tests"
# changed "T" in returned values to

err_is_pmf<-function(p, se, sp){
  n <- length(p)

  if(se==1) se=0.99999999
  if(sp==1) sp=0.99999999

  #special case with n=2
  if(n==2){
    #T can only be 1 or 3
    tests<-c(1, 3)

    #px
    q.vec<-c((1-p[1])*(1-p[2]), 1-(1-p[1])*(1-p[2]))

    #Conditonal probability table
    A.mat<-matrix(c(sp,1-sp,1-se,se), ncol=2,byrow=2)

    #pmf
    prob<-t(q.vec)%*%A.mat

    return(list("Tests"=tests,"prob"=prob))
  }
  #special case with n=3
  #serves as a starting point for the recursive relation
  else if(n==3){
    #T can only be 1 3 4 5
    tests<-c(1,3,4,5)

    #px
    q.vec<-c(c((1-p[1])*(1-p[2]), 1-(1-p[1])*(1-p[2]))*(1-p[3]),
             c((1-p[1])*(1-p[2]), 1-(1-p[1])*(1-p[2]))*p[3])

    #Conditional prbability table
    A.mat<-matrix(0,nrow=4,ncol=4)
    A.mat[1,]<-c(sp,  (1-sp)^2*sp,      sp*(1-sp), (1-sp)^3)
    A.mat[2,]<-c(1-se,se*(1-sp)*(1-se), sp*se,     se^2*(1-sp))
    A.mat[3,]<-c(1-se,se^2*sp,          se*(1-se), se^2*(1-sp))
    A.mat[4,]<-c(1-se,se^2*(1-se),      (1-se)*se, se^3)

    #"W" vector is the first row of tb
    w<-A.mat[1,]

    #pmf
    prob<-as.vector(t(q.vec)%*%A.mat)

    return(list("Tests"=tests,"prob"=prob,"w"=w))
  }
  #general recursive step
  else{
    #obtain results from last iteration
    ret<-err_is_pmf(p[1:(n-1)], se, sp)

    #dimension of the lower table
    low_dim<-length(ret$w)

    #current "w" vector
    #the first part
    w1<-c(sp,sp*ret$w[2:(low_dim)])
    #the second part
    w2<-(1-sp)^2*ret$w

    #"merge" the two parts for "w" vector
    w<-c(w1[1],w2[1],w1[2],w1[3:low_dim]+w2[2:(low_dim-1)], w2[low_dim])

    #current pmf

    #current P{T=0}, recursive relationship is not necessary
    p0<-sp*prod(1-p) + (1-se)*(1-prod(1-p))

    #subvector 1
    prob1<- c(p0,sp*(1-p[n])*ret$prob[2:low_dim]+
                (1-se)*p[n]*ret$prob[2:low_dim]+
                (se/(1-sp)-1)*(1-se)*ret$w[2:low_dim]*p[n]*prod(1-p[1:(n-1)]))

    #subvector 2
    prob2<- se*(1-sp)*(1-p[n])*ret$prob+
      ((1-sp)^2-se*(1-sp))*ret$w*prod(1-p)+
      se^2*p[n]*ret$prob


    #combine the two pieces of prob
    prob <- c(prob1[1],prob2[1],prob1[2],prob1[3:low_dim]+prob2[2:(low_dim-1)],prob2[low_dim])
    
    #support for T
    tests <- c(1,3:(2*n-1))

    return(list("Tests"=tests,"prob"=prob,"w"=w))

  }
}





#Functions to perform non-informative matrix pooling and halving

#Matrix pooling for one matrix
mat.pool<-function(p, Se, Sp, row.size, col.size) {

  #Put p's into matrix my row - could generalize this later by allowing for a matrix of vector input and only
  #  do this if it was a vector
  p.mat<-matrix(data = p, nrow = row.size, ncol = col.size, byrow = TRUE)

  I<-nrow(p.mat)
  J<-ncol(p.mat)

  probR0<-numeric(I)
  for(i in 1:I) {
    probR0[i]<-Sp*prod(1-p.mat[i,]) + (1-Se)*(1-prod(1-p.mat[i,]))
  }

  probC0<-numeric(J)
  for(j in 1:J) {
    probC0[j]<-Sp*prod(1-p.mat[,j]) + (1-Se)*(1-prod(1-p.mat[,j]))
  }

  probR0andC0<-matrix(data = NA, nrow = I, ncol = J)
  probR0orC0<-matrix(data = NA, nrow = I, ncol = J)
  for(i in 1:I) {
    for(j in 1:J) {
      Rtilde0.Ctildle0<-Sp^2 * prod(1-p.mat[i,])*prod(1-p.mat[,j])/(1-p.mat[i,j])
      Rtilde0.Ctildle1<-Sp*(1-Se) * prod(1-p.mat[i,]) * (1 - prod(1-p.mat[-i,j]))
      Rtilde1.Ctildle0<-Sp*(1-Se) * prod(1-p.mat[,j]) * (1 - prod(1-p.mat[i,-j]))
      Rtilde1.Ctildle1<-(1-Se)^2 * (1 - prod(1-p.mat[i,]) - prod(1-p.mat[,j]) + prod(1-p.mat[i,])*prod(1-p.mat[,j])/(1-p.mat[i,j]))
      probR0andC0[i,j]<-Rtilde0.Ctildle0 + Rtilde0.Ctildle1 + Rtilde1.Ctildle0 + Rtilde1.Ctildle1
      probR0orC0[i,j]<-probR0[i] + probC0[j] - probR0andC0[i,j]
    }
  }

  exp.retests<-sum(1 - probR0orC0)
  exp.T<-nrow(p.mat) + ncol(p.mat) + exp.retests
  exp.T
}




################################################################################################################
# In order to compare matrix pooling to the other procedures in the same way as in Section 5,
#   we could use rows and columns of size I to form an IxI matrix.

#Loop for random permutations of the p vector
MP.loop<-function(p, Se, Sp, row.size = row.size, col.size = col.size, num.perm = 100, seed.num=floor(10000*runif(1))) {

  set.seed(seed.num)
  keep.MP<-numeric(num.perm)
  for(i in 1:num.perm) {
    keep.MP[i]<-mat.pool(p = sample(p), Se = Se, Sp = Sp, row.size = row.size, col.size = col.size)
    #print(i)
  }

  mean(keep.MP)

}




#################################
#Utility functions for halving

#P(G_11 = 1)
prob.g1<-function(p, Se, Sp) {
  group.size<-length(p)
  if (group.size == 1) {
    prob<-0 #Only individual testing
  }
  else {
    prob<-(1-Sp)*prod(1-p) + Se*(1 - prod(1-p))
  }
  prob
}


#P(G_11 = 1 & G_2j = 1) for j = 1,2
prob.g1.g2<-function(p, Se, Sp, range) {
  group.size<-length(range)
  if (group.size == 1) {
    prob<-0 #This sets probability to 0 like how I_sj would be set to 0 in the tech report when a group test would not be done due to only one individual remaining (easier to do it this way)
  }
  else {
    prob<-(1-Sp)^2*prod(1-p) + Se*(1-Sp)*prod(1-p[range])*(1-prod(1-p[-range])) +  Se^2*(1 - prod(1-p[range]))
  }
  prob
}


#P(G_11 = 1 & G_2k & G_3j = 1) for k = 1,2 and j = 1,2,3,4
prob.g1.g2.g3<-function(p, Se, Sp, range1, range2) {
  group.size<-length(range2)
  if (group.size == 1) {
    prob<-0 #This sets probability to 0 like how I_sj would be set to 0 in the tech report when a group test would not be done due to only one individual remaining (easier to do it this way)
  }
  else {
    prob<-(1-Sp)^3*prod(1-p) +
      Se*(1-Sp)^2 * ( prod(1-p[range1]) - prod(1-p) ) +
      Se^2*(1-Sp) * ( prod(1-p[range2]) - prod(1-p[range1]) ) +
      Se^3*(1 - prod(1-p[range2]))
  }
  prob
}


#P(G_11 = 1 & G_2k & G_3m = 1 & G_4j = 1) for k = 1,2, m = 1,2,3,4, and j = 1,...,8
prob.g1.g2.g3.g4<-function(p, Se, Sp, range1, range2, range3) {
  group.size<-length(range3)
  if (group.size == 1) {
    prob<-0 #This sets probability to 0 like how I_sj would be set to 0 in the tech report when a group test would not be done due to only one individual remaining (easier to do it this way)
  }
  else {
    prob<-(1-Sp)^4*prod(1-p) +
      Se*(1-Sp)^3   * ( prod(1-p[range1]) - prod(1-p) ) +
      Se^2*(1-Sp)^2 * ( prod(1-p[range2]) - prod(1-p[range1]) ) +
      Se^3*(1-Sp)   * ( prod(1-p[range3]) - prod(1-p[range2]) ) +
      Se^4*(1 - prod(1-p[range3]))
  }
  prob
}


#######################################
#Main halving function

# Note: Brianna Hitt 02-23-19
# Renamed from halving() to halving.comp() to avoid
#   conflicts with Michael Black's halving function
# This function is needed to calculate measures
#   for non-informative halving, for comparisons
#   with Sterrett algorithms in the Sterrett() function
# No documentation for this function is provided

halving.comp <- function(p, Se, Sp) {

  I<-length(p)

  #3-stage halving
  if (I > 2) {
    I21<-floor(I/2) #If group size is odd, the lower value goes into I21; e.g., I = 5 leads to I21 = 2 and I22=3 (one could switch the definition here if needed)
    I22<-I-I21
    exp.3H<-1 + 2*prob.g1(p=p, Se=Se, Sp=Sp) + I21*prob.g1.g2(p=p, Se=Se, Sp=Sp, range=1:I21) +
      I22*prob.g1.g2(p=p, Se=Se, Sp=Sp, range=(I-I22+1):I)

  } else
  {
    exp.3H<-NULL
  }


  #4-stage halving
  if (I > 4) {
    I31<-floor(I21/2)
    I32<-I21-I31
    I33<-floor(I22/2)
    I34<-I22-I33
    exp.4H<-1 + 2*prob.g1(p=p, Se=Se, Sp=Sp) + 2*prob.g1.g2(p=p, Se=Se, Sp=Sp, range=1:I21) +
      2*prob.g1.g2(p=p, Se=Se, Sp=Sp, range=(I-I22+1):I) +
      I31*prob.g1.g2.g3(p=p, Se=Se, Sp=Sp, range1=1:I21, range2=1:I31) +
      I32*prob.g1.g2.g3(p=p, Se=Se, Sp=Sp, range1=1:I21, range2=(I21-I32+1):I21) +
      I33*prob.g1.g2.g3(p=p, Se=Se, Sp=Sp, range1=(I-I22+1):I, range2=(I21+1):(I21+I33)) +
      I34*prob.g1.g2.g3(p=p, Se=Se, Sp=Sp, range1=(I-I22+1):I, range2=(I-I34+1):I)
  } else
  {
    exp.4H<-NULL
  }


  #5-stage halving
  if (I > 8) {
    I41<-floor(I31/2)
    I42<-I31-I41
    I43<-floor(I32/2)
    I44<-I32-I43
    I45<-floor(I33/2)
    I46<-I33-I45
    I47<-floor(I34/2)
    I48<-I34-I47
    exp.5H<-1 + 2*prob.g1(p=p, Se=Se, Sp=Sp) + 2*prob.g1.g2(p=p, Se=Se, Sp=Sp, range=1:I21) +
      2*prob.g1.g2(p=p, Se=Se, Sp=Sp, range=(I-I22+1):I) +
      2*prob.g1.g2.g3(p=p, Se=Se, Sp=Sp, range1=1:I21, range2=1:I31) +
      2*prob.g1.g2.g3(p=p, Se=Se, Sp=Sp, range1=1:I21, range2=(I21-I32+1):I21) +
      2*prob.g1.g2.g3(p=p, Se=Se, Sp=Sp, range1=(I-I22+1):I, range2=(I21+1):(I21+I33)) +
      2*prob.g1.g2.g3(p=p, Se=Se, Sp=Sp, range1=(I-I22+1):I, range2=(I-I34+1):I) +

      I41*prob.g1.g2.g3.g4(p=p, Se=Se, Sp=Sp, range1=1:I21, range2=1:I31, range3=1:I41) +
      I42*prob.g1.g2.g3.g4(p=p, Se=Se, Sp=Sp, range1=1:I21, range2=1:I31, range3=(I31-I42+1):I31) +
      I43*prob.g1.g2.g3.g4(p=p, Se=Se, Sp=Sp, range1=1:I21, range2=(I21-I32+1):I21, range3=(I31+1):(I31+I43)) +
      I44*prob.g1.g2.g3.g4(p=p, Se=Se, Sp=Sp, range1=1:I21, range2=(I21-I32+1):I21, range3=(I21-I44+1):I21) +

      I45*prob.g1.g2.g3.g4(p=p, Se=Se, Sp=Sp, range1=(I-I22+1):I, range2=(I21+1):(I21+I33), range3=(I21+1):(I21+I45)) +
      I46*prob.g1.g2.g3.g4(p=p, Se=Se, Sp=Sp, range1=(I-I22+1):I, range2=(I21+1):(I21+I33), range3=(I21+I45+1):(I21+I33)) +
      I47*prob.g1.g2.g3.g4(p=p, Se=Se, Sp=Sp, range1=(I-I22+1):I, range2=(I-I34+1):I, range3=(I21+I33+1):(I21+I33+I47)) +
      I48*prob.g1.g2.g3.g4(p=p, Se=Se, Sp=Sp, range1=(I-I22+1):I, range2=(I-I34+1):I, range3=(I-I48+1):I)
  } else
  {
    exp.5H<-NULL
  }

  c(exp.3H, exp.4H, exp.5H)
}

#
