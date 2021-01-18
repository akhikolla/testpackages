###############################################################################
#     Michael Black 04 24 2013
#     Function: contains all the functions needed to run the examples in
#               Section 4.4 of the dissertation
#
#     List of functions:
#            1. for beta.dist()
#                a. f.p.ord()        gives the pdf for ordered p_i
#                b. p.f.p.ord()      p times pdf for p_i
#                c. beta.dist()      finds vector of ordered p_i
#
#            2. for halving()
#                a. get.tests()      iterates through number of tests
#                b. sub.grp.size()   splits groups into subgroups
#                c. halving()        finds pmf, ET and Var for halving
#
#            3.  for CRC 3-stage used by hierarchical and get.CRC
#                a. Exp.Tk.3step     is for 3 steps gets the portion of E(T) 
#                     for the given grouping
#                b. grps.iter.3step   iterates through all the possible 
#                     groupings and sizes to find the optimal
#                c. Opt.grps.3step   calls grps.iter.3step for a specified 
#                     group p and number of groupings g
#                d. Opt.grps.size_number.3step   gets opt groupings and opt 
#                     sizes It calls Opt.grps.3step so looks at all possible 
#                     groupings
#                e. Opt.grps.size_number_speedg.3step  is a simple modification 
#                     to the above that lets us stop sooner
#
#            4. for ORC 3-stage used by hierarchical and get.CRC
#                a. Exp.T.3step        uses the Exp.Tk
#                b. get.fin.opt.3step  iterates until the optimal grouping is
#                     found, it assumes that
#                c. Opt.grps.speed1.3step  selects an initial estimate for 
#                     groupings Then it calls the optimizer above for a 
#                     specified number of groupings g
#                d. Opt.grps.size.speed.3step  iterates through the number of
#                     subgroups at stage 2 called g to find optimal
#                e. Opt.grps.size_number.speed.3step  is a simple modification 
#                     to the above that lets us stop sooner
#
#            5. functions for Accuracy used by hierarchical and get.CRC
#                a. ProbY         Probability of actually testing positive 
#                                 using 3-step regrouping
#                b. ProbY.0       1-Pooling specificity, obtained by changing
#                                 the pi to zero one at a time and finding 
#                                 the pooled probability of testing positive.
#                c. ProbY.1       Pooling sensitivity
#            Next same as above for 4-stages
#                d. ProbY.4s
#                e. ProbY.4s.0
#                f. ProbY.4s.1
#
#            6. 4-stage CRC  used by hierarchical and get.CRC
#                a. Exp.T.4step      first we get the ET for 4step.
#                b. steep.move       gives the descent part of steepest descent
#                c. grps.iter.both   iterates through m2 an m3 simultaneously
#                d. both.call.vec.g3  iterating the vec.g3 possibilities for 
#                     each modified vec.g2
#                e. grps.iter.size.3s  iterates through the possible g2 given g3
#                f. Opt.grps.size_number.4step iterates through the possible g3
#
#            7. 4 stage ORC
#                a. Exp.T.4step.slow   first we get the ET for 4step.
#                b. grps.iter.4step.slow   Iterates through all the possible 
#                     combinations for vec.g3 given a g3 and vec.g2 this will 
#                     be replaced with faster
#                c. grps.iter.4s.slow  iterates through all the possible vec.g2 
#                     combinations given a g3
#                d. grps.iter.size.4s.slow   iterates through the possible g2 
#                     given a g3
#                e. Opt.grps.size_number.4step.slow   iterates through the 
#                     possible g3
#
#            8. converting notation
#                a. get.mc.3s
#                b. get.mc.4s
#
#            9. function called for hierarchical.desc()
#                a. get.g2_g3    to covert I2 and I3 to the form used in the
#                     programs
#                b. hierarchical.desc    gets descriptive values for given 
#                     vector of p's and configuration
#            10. get.CRC         gets ORC or CRC for give vector of p's
#
#
#
###############################################################################



#Start beta.dist() functions
###############################################################################
#    Michael Black 03 12 2013
#     Function: Get E(p(i)) from a Beta distn. given a p and alpha or 
#               beta and alpha
#       calls: f.p.ord and p.f.p.ord which provide functions for beta 
#              distribution
#       inputs: p average probability = alpha/(alpha + beta),
#               or a b = beta parameter for the beta distribution,
#               a = alpha the alpha parameter for the beta distribution
#               grp.sz the number of individuals in the group
#               plot = FALSE produces plot of the distribution with p(i) marked
#       if alpha = 0 uses an extreme value distribution based on a bernoulli(p) 
#         for the p_i
#       if alpha = "inf" uses p_i = p for all i
#       note: including a value for b means function ignores p unless a = "inf"
#
#       demo(plotmath) will produce formatting notes for R labels creation
#
#
###############################################################################

#PDF of ordered p_(i)
f.p.ord <- function(p, i, grp.sz, a, b) {
  factorial(grp.sz)/(factorial(i - 1)*factorial(grp.sz - i))*
    dbeta(x = p, shape1 = a, shape2 = b)*
    pbeta(q = p, shape1 = a, shape2 = b)^(i - 1)*
    (1 - pbeta(q = p, shape1 = a, shape2 = b))^(grp.sz - i)
}

#p * PDF of ordered p_(i) - used to find E(p_(i))
p.f.p.ord <- function(p, i, grp.sz, a, b) {
  p*factorial(grp.sz)/(factorial(i - 1)*factorial(grp.sz - i))*
    dbeta(x = p, shape1 = a, shape2 = b)*
    pbeta(q = p, shape1 = a, shape2 = b)^(i - 1)*
    (1 - pbeta(q = p, shape1 = a, shape2 = b))^(grp.sz - i)
}



# Start beta.dist()
###############################################################################

# Brianna Hitt - 04.02.2020
# changed print() to message() for easy suppression 

# function to get vector of p(i)
# beta.dist <- function(p = 0.05,alpha = 1,beta = NULL, grp.sz = 10, 
#                       simul = FALSE, plot = FALSE, 
#                       rel.tol = ifelse(a >= 1, .Machine$double.eps^0.25, .Machine$double.eps^0.1)) {
#   a <- alpha
#   b <- beta
#   if (a == "inf"|| a == "hom") p.vec = rep(p, grp.sz)
#   else if (a == 0) {
#     p.vec <- numeric(grp.sz)
#     save.err <- numeric(grp.sz)
#     save.int <- 0
#     for(i in 1:grp.sz) {
#       save.int <- save.int + choose(grp.sz, (i - 1))*(p^(grp.sz - i + 1))*((1 - p)^(i - 1))
#       p.vec[i] <- save.int
#       save.err[i] <- 0
#     }
#   }
# 
#   else {
#     if (is.null(b)) b = a*(1/p - 1)
#     if (simul == FALSE) {
#       p.vec <- numeric(grp.sz)
#       save.err <- numeric(grp.sz)
# 
#       for(i in 1:grp.sz) {
#         save.int <- integrate(f = p.f.p.ord, lower = 0, upper = 1, i = i, 
#                               grp.sz = grp.sz, a = a, b = b, rel.tol = rel.tol)
#         p.vec[i] <- save.int$value
#         save.err[i] <- save.int$abs.error
#       }
#     }
#     else  {
#       message("Using simulation")
#       presort <- matrix(rbeta(10000*grp.sz, a, b), ncol = grp.sz, byrow = 10000)
#       sorted <- t(apply(presort,1,sort))
#       p.vec <- colMeans(sorted)
#     }
# 
#   }
#   if (plot == TRUE && is.numeric(a) && a != 0) {
#     get.beta  <- c(sort(rbeta(1000, a, b)), 1)
#     y.beta <- dbeta(get.beta, a, b)
#     max.y <- min(max(y.beta), 100)
# 
# 
#     plot(x = get.beta, y = y.beta, type = "l",
#          lwd = 1, col = "darkgreen", panel.first=grid(col = "gray"),
#          xlim = c(0, 1), ylim = c(0, max.y),
#          xlab = substitute(paste("PDF and ", italic(E)(italic(p)[(italic(i))]) , 
#                                  " for beta distribution with I = ", grp.sz, ", ", 
#                                  alpha, " = ", a, " and ", beta, " = ", b), 
#                            list(grp.sz = grp.sz, a = a, b = b)),
#          ylab = substitute(paste(italic(f)(italic(p)[italic(i)])), list()))
#     points(x = p.vec, y = rep(0.5, length(p.vec)), type = "h", lwd = 1)
# 
#   }
#   if (plot == TRUE && !is.numeric(a) || a == 0){
#     warning("Plot not available for homogeneous or extreme")
#   } 
#   p.vec
# }



# Start beta.dist2()
###############################################################################
# function to get vector of p(i)

# Brianna Hitt - 04.02.2020
# changed print() to message() and warning(), as appropriate, for easy suppression 

beta.dist2 <- function(p = 0.05, alpha = 1, beta = NULL, grp.sz = 10, 
                       num.sim = 10000, simul = FALSE, plot = FALSE, 
                       rel.tol = ifelse(a >= 1, .Machine$double.eps^0.25,  .Machine$double.eps^0.1)) {
  a <- alpha
  b <- beta
  if (a == "inf"|| a == "hom") p.vec = rep(p, grp.sz)
  else if (a == 0) {
    p.vec <- numeric(grp.sz)
    save.err <- numeric(grp.sz)
    save.int <- 0
    for(i in 1:grp.sz) {
      save.int <- save.int + choose(grp.sz, (i - 1))*(p^(grp.sz - i + 1))*((1 - p)^(i - 1))
      p.vec[i] <- save.int
      save.err[i] <- 0
    }
  }

  else {
    if (is.null(b)) b = a*(1/p - 1)
    if (simul == FALSE) {
      p.vec <- numeric(grp.sz)
      save.err <- numeric(grp.sz)

      for(i in 1:grp.sz) {
        save.int <- integrate(f = p.f.p.ord, lower = 0, upper = 1, i = i, 
                              grp.sz = grp.sz, a = a, b = b, rel.tol = rel.tol)
        p.vec[i] <- save.int$value
        save.err[i] <- save.int$abs.error
      }
    }
    else  {
      message("Using simulation")
      presort <- matrix(rbeta(num.sim*grp.sz, a, b), ncol = grp.sz, 
                        byrow = num.sim)
      sorted <- t(apply(presort,1,sort))
      p.vec <- colMeans(sorted)
    }

  }
  if (plot == TRUE && is.numeric(a) && a != 0) {
    get.beta  <- c(sort(rbeta(1000, a, b)), 1)
    y.beta <- dbeta(get.beta, a, b)
    max.y <- min(max(y.beta), 100)


    plot(x = get.beta, y = y.beta, type = "l",
         lwd = 1, col = "darkgreen", panel.first=grid(col = "gray"),
         xlim = c(0, 1), ylim = c(0, max.y),
         xlab = substitute(paste("PDF and ", italic(E)(italic(p)[(italic(i))]) , 
                                 " for beta distribution with I = ",  grp.sz, 
                                 ", ", alpha, " = ", a, " and ", beta,  " = ", b), 
                           list(grp.sz = grp.sz, a = a, b = b)),
         ylab = substitute(paste(italic(f)(italic(p)[italic(i)])), list()))
    points(x = p.vec, y = rep(0.5, length(p.vec)), type = "h", lwd = 1)

  }
  if (plot == TRUE && !is.numeric(a) || a == 0){
    warning("Plot not available for homogeneous or extreme")
  } 
  p.vec
}
# End beta.dist() functions



# Start halving() functions
###############################################################################
#    Michael Black 03 12 2013
#     Function: Halving PMF
#       calls:  get.tests combines possible number of tests
#               sub.grp.size  gets the subgroup sizes 
#                 (could be added to main function)
#       inputs: se and sp sensitivity and specificity,
#               p a vector of individual probabilities
#               stages the number of stages for halving
#               order.p = TRUE if FALSE p.vec is not sorted for halving
###############################################################################

#get.tests is called by halving
get.tests <- function(t1, t2) {
  t1 <- data.matrix(t1)
  t2 <- data.matrix(t2)
  one <- data.matrix(c(rep(x = 1, times = length(t1))))
  Tests <- c(1, t(one%x%t1 + t2%x%one + one%x%one))
  Tests
}



#sub.grp.size called by halving
sub.grp.size <- function(I.0, stages) {

  if (stages == 2) {
    return(I.0)
  }
  else if (stages == 3) {
    return(c(floor(I.0/2), I.0-floor(I.0/2)))
    #return(c(ceiling(I.0/2), I.0 - ceiling(I.0/2)))
  }
  else {
    return(c(sub.grp.size(I.0 = floor(I.0/2), stages = stages - 1), 
             sub.grp.size(I.0 = (I.0 - floor(I.0/2)), stages=stages - 1)))
    #return(c(sub.grp.size(I.0 = ceiling(I.0/2), stages = stages - 1), 
    #         sub.grp.size(I.0 = (I.0 - ceiling(I.0/2)), stages=stages - 1)))
  }

}



# Start halving
###############################################################################

#' @title Probability mass function for halving
#'
#' @description Calculate the probability mass function for the number of tests 
#' from using the halving algorithm.
#' @param p a vector of individual risk probabilities.
#' @param sp the specificity of the diagnostic test.
#' @param se the sensitivity of the diagnostic test.
#' @param stages the number of stages for the halving algorithm.
#' @param order.p logical; if TRUE, the vector of individual risk probabilities 
#' will be sorted.
#'
#' @details Halving algorithms involve successively splitting a positive
#' testing group into two equal-sized halves (or as close to equal as possible) 
#' until all individuals have been identified as positive or negative. 
#' \eqn{S}-stage halving begins by testing the whole group of \eqn{I}
#' individuals. Positive groups are split in half until the final stage of the 
#' algorithm, which consists of individual testing. For example, consider an 
#' initial group of size \eqn{I}=16 individuals. Three-stage halving (3H) 
#' begins by testing the whole group of 16 individuals. If this group tests 
#' positive, the second stage involves splitting into two groups of size 8. 
#' If either of these groups test positive, a third stage involves testing each 
#' individual rather than halving again. Four-stage halving (4H) would continue 
#' with halving into groups of size 4 before individual testing. Five-stage 
#' halving (5H) would continue with halving into groups of size 2 before 
#' individual testing. 3H requires more than 2 individuals, 4H requires more 
#' than 4 individuals, and 5H requires more than 8 individuals.
#'
#' This function calculates the probability mass function, expected testing 
#' expenditure, and variance of the testing expenditure for halving algorithms 
#' with 3 to 5 stages.
#'
#' @return A list containing:
#' \item{pmf}{the probability mass function for the halving algorithm.}
#' \item{et}{the expected testing expenditure for the halving algorithm.}
#' \item{vt}{the variance of the testing expenditure for the halving 
#' algorithm.}
#'
#' @author This function was originally written by Michael Black for Black 
#' et al. (2012). The function was obtained from 
#' \url{http://chrisbilder.com/grouptesting}. Minor modifications have been 
#' made for inclusion of the function in the \code{binGroup2} package.
#'
#' @references
#' \insertRef{Black2012}{binGroup2}
#'
#' @seealso
#' \code{\link{expectOrderBeta}} for generating a vector of individual risk 
#' probabilities for informative group testing.
#'
#' @family operating characteristic functions
#'
#' @examples
#' # Equivalent to Dorfman testing (two-stage hierarchical)
#' halving(p=rep(0.01, 10), sp=1, se=1, stages=2, 
#'         order.p=TRUE)
#'
#' # Halving over three stages; each individual has a 
#' #   different probability of being positive
#' set.seed(12895)
#' p.vec <- expectOrderBeta(p=0.05, alpha=2, grp.sz=20)
#' halving(p=p.vec, sp=0.95, se=0.95, stages=3, 
#'         order.p=TRUE)

halving <- function(p, sp = 1, se = 1, stages = 2, order.p = TRUE) {

  if (order.p == TRUE) p <- sort(p)
  N <- length(p)

  if(stages < 2) {
    stages <- 2
    warning("stages out of range, using Dorfman")
  }

  #these make sure that stages is in bounds for the program
  if(stages > 5) {
    stages = 5
    warning("Too many stages, go back to 5")
  }
  max.stages <- floor(log(N, 2)) + 1
  if (stages > max.stages) {
    stages <- max.stages
    warning("Too many stages using max.stages")
  }

  # This splits it into final sub-groups
  fin.cnt <- sub.grp.size(I.0 = N, stages = stages)
  f <- length(fin.cnt)

  tb <- matrix(c(sp, 1 - sp, 1 - se, se), ncol = 2, byrow = 2)

  if (stages > 2){

    #the for loop below gets the matrix for the testing error part
    for (i in 2:(stages - 1)) {
      c1 <- c(sp, rep(1 - se, (2^(2^(i - 1)) -1)))
      ident <- diag(1 - c1)
      tb2 <- (tb%x%tb)
      tb3 <- t(t(tb2)%*%ident)
      tb <- cbind(c1, tb3)

    }
    # below gets the tests matrix for given final sub-groups
    tests.start <- NULL
    for (i in 1:f) {
      tests.start <- rbind(tests.start, c(1, 1 + fin.cnt[i]))
    }

    for (i in 1:(stages - 2)){
      tests.fin <- NULL
      for (k in 1:(f/(2^i))) {
        tests.fin <- rbind(tests.fin, get.tests(t1 = tests.start[(2*k - 1), ], 
                                                t2 = tests.start[(2*k), ]))
      }
      tests.start <- tests.fin
    }
    Tests <- tests.start

  }
  if (stages == 2) {
    Tests <- t(c(1, 1 + N))
  }
  # gets the probability vector
  prob.p <- 1
  p.start <- 0
  for (j in 1:f) {
    p.end <- p.start + fin.cnt[j]
    p.part <- p[(p.start + 1):p.end]
    p.part.prod <- prod(1 - p.part)
    add.prob <- t(c(p.part.prod, 1 - p.part.prod))
    prob.p <- add.prob%x%prob.p
    p.start <- p.end
  }
  prob.ts <- prob.p%*%tb
  test.prop <- rbind(Tests, prob.ts)
  #finalizes the pmf
  pmf <- rbind((by(test.prop[1, ], test.prop[1, ], sum)/by(test.prop[1, ],
                                                           test.prop[1, ], length)),
               by(test.prop[2, ], test.prop[1, ], sum))
  #get the E(T)
  et <- sum(pmf[1, ]*pmf[2, ])
  #get variance for T
  et2 <- sum(((pmf[1, ])^2)*pmf[2, ])
  vt <- et2 - et^2

  pmf <- data.frame(num.tests = round(pmf[1, ]), 
                    prob.tests = round(pmf[2, ], 4), row.names = NULL)

  list(pmf = pmf, et = et, vt = vt)
}
# End halving() functions



# Start  hierarchical.desc functions.
# First auxilary functions that are called by hierarchical
    #Start three-stage optimal functions()
###############################################################################
# Function for grouping N individuals into g groups of optimized sizes
# 2/1/2011 Michael Black
#
# The idea here is to set up a recursive call to the function
# I1 goes from 1 to N-g+1
# I2 goes from ceiling1 to N-I1-g+2 , etc
# Ig is N-I1-I2-...-I(g-1)
#     this will need a conditional for final status
#
###############################################################################

# Exp.Tk.3step is for 3 steps gets the portion of E(T) for the given grouping
## Brianna Hitt 11-13-19 - allow Se, Sp values to vary across stages of testing

Exp.Tk.3step <- function(p.group, p.rest, sp = 1, se = 1) {
  group.size <- length(p.group)
  if (group.size > 1) {
    # below is Michael Black's original code
    # ETk <- group.size*(((1 - sp)^2)*prod(1 - p.group)*prod(1 - p.rest) + 
    #                     se*(1 - sp)*(prod(1 - p.group)*(1 - prod(1 - p.rest))) + 
    #                     (se^2)*(1 - prod(1 - p.group)))
    
    # below allows Se, Sp to vary across stages of testing
    ETk <- group.size*(prod(1-sp[1:2])*prod(1 - p.group)*prod(1 - p.rest) + 
                         se[1]*(1 - sp[2])*prod(1 - p.group)*(1 - prod(1 - p.rest)) + 
                         prod(se[1:2])*(1 - prod(1 - p.group)))
  }
  else if (group.size == 1){
    ETk <- 0
  }
  ETk
}



###############################################################################
#grps.iter.3step iterates through all the possible groupings and sizes to 
#   find the optimal
# grps.iter.3step <- function(p, N, g, grps.size, opt.grps, fin.opt, sp, se) {
#   N.all <- length(p)
#   g.all <- length(opt.grps) - 1
#   grps.start <- 1
#   grps.fin <- N - g + 1
#   if (g == 1) {
# 
#     opt.grps[g.all - g + 1] <- N.all - sum(opt.grps[1:g.all - 1])
#     count <- 0
#     opt.grps[g.all + 1] <- Exp.T.3step(p = p, opt.grps[1:g.all], 
#                                        sp = sp, se = se)
# 
#     if (opt.grps[g.all + 1] < fin.opt[g.all + 1]) {
#       fin.opt <- opt.grps
#     }
#   }
#   else if (g > 1) {
#     for (i in grps.start:grps.fin) {
#       grps.size <- i
#       opt.grps[g.all - g + 1] <- i
#       fin.opt <- grps.iter.3step(p = p, N = N - i, g = g - 1, grps.size = i, 
#                                  opt.grps = opt.grps, fin.opt = fin.opt, 
#                                  sp = sp, se = se)
#     }
#   }
#   fin.opt
# }



#################################
#Opt.grps.3step calls grps.iter.3step for a specified group p and 
#  number of groupings g
# Opt.grps.3step <- function(p, g = 2, sp = 1, se = 1){
#   N <- length(p)
#   opt.grps <- rep(0, g + 1)
#   grps.size <- N - g + 1
#   fin.opt <- c(rep(0, g), N + 2)
#   out.this <- grps.iter.3step(p = p, N = N, g = g, grps.size = grps.size, 
#                               opt.grps = opt.grps, fin.opt = fin.opt, 
#                               sp = sp, se = se)
#   out.this
# }



###
#This next program gets opt groupings and opt sizes
#   It calls Opt.grps.3step so looks at all possible groupings
###

# Opt.grps.size_number.3step <- function(p, sp = 1, se = 1) {
#   start.time <- proc.time()
#   N <- length(p)
#   expt.t <- 2*N
#   Opt.all <- NULL
#   for (g in 1:N) {
#     Opt.sizes <- Opt.grps.3step(p = p, g = g, sp = sp, se = se)
#     if (Opt.sizes[g + 1] < expt.t) {
#       expt.t <- Opt.sizes[g + 1]
#       Opt.all <- Opt.sizes
#     }
#   }
# 
#   end.time <- proc.time()
#   save.time <- end.time-start.time
#   print("\n Number of minutes running for I = ", N, ":", save.time[3]/60, "\n \n")
# 
#   Opt.all
# }



# Opt.grps.size_number_speedg.3step is a simple modification to the above that
#   lets us stop sooner
# Opt.grps.size_number_speedg.3step<-function(p, sp = 1, se = 1) {
#   start.time <- proc.time()
#   N <- length(p)
#   expt.t <- 2*N
#   Opt.all <- NULL
#   for (g in 1:N) {
#     Opt.sizes <- Opt.grps.3step(p = p, g = g, sp = sp, se = se)
#     if (Opt.sizes[g + 1] < expt.t) {
#       expt.t <- Opt.sizes[g + 1]
#       Opt.all <- Opt.sizes
#     }
#     if (Opt.sizes[g + 1] > expt.t) break
#   }
# 
#   end.time <- proc.time()
#   save.time <- end.time - start.time
# 
#   Opt.all
# }



###############################################################################
#
#   The following improve the speed of the program
#       1 Set an initial grouping
#       2 add and subtract one from different groupings
#           to find the best groupings
#
###############################################################################
# First get a function for Exp.T given a grouping
# this uses the Exp.Tk
## Brianna Hitt 11-13-19 - allow Se, Sp values to vary across stages of testing

Exp.T.3step <- function(p, grping, sp = 1, se = 1){
  g.all <- length(grping) # number of groups in second stage of testing
  N <- length(p)
  g.ord <- c(0, cumsum(grping))
  if (g.all > 1) {
    # below is Michael Black's original code
    # Exp.T<- 1 + g.all*((1 - sp)*prod(1 - p) + se*(1 - prod(1 - p)))

    # below allows Se, Sp to vary across stages of testing    
    Exp.T <- 1 + g.all*((1 - sp[1])*prod(1 - p) + se[1]*(1 - prod(1 - p)))
    for (i in 1:g.all) {
      Exp.T <- Exp.T + Exp.Tk.3step(p.group = p[(g.ord[i] + 1):g.ord[i + 1]], 
                                    p.rest = p[-((g.ord[i] + 1):g.ord[i + 1])], 
                                    sp = sp, se = se)
    }
  }
  # Brianna 10.04.19 - initial group is retested
  if (g.all == 1) {
    # below is Michael Black's original code
    # Exp.T <- 1 + grping*((1 - sp)*prod(1 - p) + se*(1 - prod(1 - p)))
    
    # below allows Se, Sp to vary across stages of testing
    Exp.T <- 1 + grping*((1 - sp[1])*prod(1 - p) + se[1]*(1 - prod(1 - p)))
  }
  Exp.T
}



#This Function iterates until the optimal grouping is found, it assumes that
#     the Exp.t is a convex function of the groupings
# get.fin.opt.3step <- function(p, fin.opt, sp = 1, se = 1) {
#   g <- length(fin.opt) - 1
#   N <- length(p)
#   init.grps <- fin.opt[1:g]
#   chk.chg <- 0
#   for (k in 1:g) {
#     for (j in 1:2) {
#       for (i in 1:g) {
#         if (j == 1) {
#           if (init.grps[i] > 1) {
#             init.grps.mod <- init.grps
#             init.grps.mod[i] <- init.grps.mod[i] - 1
#             init.grps.mod[k] <- init.grps.mod[k] + 1
#             chk.opt <- c(init.grps.mod, Exp.T.3step(p, grping = init.grps.mod, 
#                                                     sp = sp, se = se))
#             if ((chk.opt[g + 1]) < (fin.opt[g + 1])) {
#               fin.opt <- chk.opt
#               chk.chg <- 1
#             }
#           }
#         }
#         if (j == 2) {
#           if (init.grps[k] > 1) {
#             init.grps.mod <- init.grps
#             init.grps.mod[i] <- init.grps.mod[i] + 1
#             init.grps.mod[k] <- init.grps.mod[k] - 1
#             chk.opt <- c(init.grps.mod, Exp.T.3step(p, grping = init.grps.mod, 
#                                                     sp = sp, se = se))
#             if (chk.opt[g + 1] < fin.opt[g + 1]) {
#               fin.opt <- chk.opt
#               chk.chg <- 1
#             }
#           }
#         }
# 
#       }
#     }
#   }
# 
#   if (chk.chg == 1) return(get.fin.opt.3step(p, fin.opt, sp = sp, se = se))
#   if (chk.chg == 0) return(fin.opt)
# }



# This function selects an initial estimate for groupings
#   Then it calls the optimizer above for a specified number of groupings g
# 6/15/11 Updated the starting point to be more efficient with ordering
# Opt.grps.speed1.3step <- function(p, g, sp = 1, se = 1) {
#   N <- length(p)
#   init.grps <- c((N - (g - 1)), rep(1, g - 1))
#   fin.opt <- c(init.grps, Exp.T.3step(p, grping = init.grps, sp = sp, se = se))
#   fin.opt.out <- get.fin.opt.3step(p, fin.opt, sp = sp, se = se)
#   fin.opt.out
# 
# }



# This iterates through the number of subgroups at stage 2 called g to find 
#   optimal
# Opt.grps.size.speed.3step <- function(p, sp = 1, se = 1) {
#   start.time <- proc.time()
#   N <- length(p)
#   expt.t <- 2*N
#   Opt.all <- NULL
#   for (g in 1:N) {
#     Opt.sizes <- Opt.grps.speed1.3step(p = p, g = g, sp = sp, se = se)
#     if (Opt.sizes[g + 1] < expt.t) {
#       expt.t <- Opt.sizes[g + 1]
#       Opt.all <- Opt.sizes
#     }
#   }
# 
#   end.time <- proc.time()
#   save.time <- end.time - start.time
#   print("\n Number of minutes running for I = ", N, ":", save.time[3]/60, "\n \n")
# 
#   Opt.all
# }



#Below is a simple modification to the above that lets us stop sooner
# Opt.grps.size_number.speed.3step <- function(p, sp = 1, se = 1) {
#   start.time <- proc.time()
#   N <- length(p)
#   expt.t <- 2*N
#   Opt.all <- NULL
#   stop.count<-0
#   for (g in 1:N) {
#     Opt.sizes <- Opt.grps.speed1.3step(p = p, g = g, sp = sp, se = se)
#     if (Opt.sizes[g + 1] < expt.t) {
#       expt.t <- Opt.sizes[g + 1]
#       Opt.all <- Opt.sizes
#     }
#     if (Opt.sizes[g + 1] > expt.t) stop.count <- stop.count + 1
#     if (stop.count > 1) break
#   }
# 
#   end.time <- proc.time()
#   save.time <- end.time - start.time
# 
#   Opt.all
# }
#End three-stage functions



  #start accuracy measure functions
###############################################################################
#
#   Author: Michael Black
#   Last Updated: 5-25-2012
#   Purpose: These functions were added to help get the
#             Expected individual testing probabilities and
#             diagnostic statistics
#         1. ProbY  gets the probability an individual tests positive
#         2. ProbY.0 gets individual 1 - specificity
#         3. ProbY.1 gets individual sensitivity
#         4-6 same for 4s
#
###############################################################################

# Probability of actually testing using 3-step regrouping
ProbY <- function(p.vec, g3, sp=1, se=1) {
  g <- length(g3)
  N <- length(p.vec)
  p.ord <- p.vec
  g.ord <- c(0, cumsum(g3))
  Prob.Y <- NULL
  for (i in 1:g) {
    p.g2 <- p.ord[(g.ord[i] + 1):g.ord[i + 1]]
    I.g2 <- length(p.g2)
    if (I.g2 > 1) {
      for (j in 1:I.g2) {
        p.j <- p.g2[j]
        # below is Michael Black's original code
        probYis <- ((1 - sp)^3)*prod(1 - p.ord) + 
          se*((1 - sp)^2)*(prod(1 - p.g2) - prod(1 - p.ord)) + 
          (se^2)*(1 - sp)*((1 - p.j) - prod(1 - p.g2)) + (se^3)*p.j
        Prob.Y <- rbind(Prob.Y, probYis)
      }
    }
    else if (I.g2 == 1) {
      p.j <- p.g2
      # below is Michael Black's original code
      probYis <- ((1 - sp)^2)*prod(1 - p.ord) + 
        se*(1 - sp)*((1 - p.j) - prod(1 - p.ord)) + (se^2)*p.j
      Prob.Y<-rbind(Prob.Y, probYis)
    }
  }
  Prob.Y
}



# 1-Pooling specificity: Is specific for each individual obtained by changing
#                        the pi to zero one at a time and finding the pooled 
#                        probability of testing positive.
# Could be done on all the positives and negatives to see how they fall out.
## Brianna Hitt 11-13-19 - allow Se, Sp values to vary across stages of testing

ProbY.0 <- function(p.vec, g3, sp=1, se=1) {
  g <- length(g3) #Number of subgroups at stage 2
  N <- length(p.vec)
  p.ord <- p.vec
  g.ord <- c(0, cumsum(g3))
  Prob.Y <- NULL
  for (i in 1:g) {
    p.g2 <- p.ord[(g.ord[i] + 1):g.ord[i + 1]]
    I.g2 <- length(p.g2)
    # Brianna 09.30.19 - individuals who are decoded in 3 stages
    if (I.g2 > 1) {
      for (j in 1:I.g2) {
        p.g2.mod <- p.g2
        p.ord.mod <- p.ord
        p.j <- 0
        p.g2.mod[j] <- 0
        p.ord.mod[g.ord[i] + j] <- 0
        # below is Michael Black's original code
        # probYis <- ((1 - sp)^3)*prod(1 - p.ord.mod) + 
        #   se*((1 - sp)^2)*(prod(1 - p.g2.mod) - prod(1 - p.ord.mod)) + 
        #   (se^2)*(1 - sp)*((1 - p.j) - prod(1 - p.g2.mod)) + (se^3)*p.j
        
        # below allows Se, Sp to vary across stages of testing
        probYis <- prod(1 - sp[1:3])*prod(1 - p.ord.mod) + 
          se[1]*prod(1-sp[2:3])*(prod(1 - p.g2.mod) - prod(1 - p.ord.mod)) + 
          prod(se[1:2])*(1 - sp[3])*((1 - p.j) - prod(1 - p.g2.mod)) + 
          prod(se[1:3])*p.j
        Prob.Y <- rbind(Prob.Y, probYis)
      }
    }
    # Brianna 09.30.19 - individuals who are decoded in 2 stages
    else if (I.g2 == 1) {
      # Michael Black's original code - p.j is the individual of interest
      p.j <- 0
      p.ord.mod <- p.ord
      p.ord.mod[g.ord[i] + 1] <- 0
      
      # below is Michael Black's original code
      # probYis <- ((1 - sp)^2)*prod(1 - p.ord.mod) + 
      #   se*(1 - sp)*((1 - p.j) - prod(1 - p.ord.mod)) + (se^2)*p.j

      # Brianna 09.30.19 
      # below allows Se, Sp to vary across stages of testing
      probYis <- prod(1 - sp[1:2])*prod(1 - p.ord.mod) + 
        se[1]*(1 - sp[2])*((1 - p.j) - prod(1 - p.ord.mod)) + prod(se[1:2])*p.j
      
      Prob.Y <- rbind(Prob.Y, probYis)
    }
  }
  # Brianna 09.30.19 - added this to remove "probYis" from row names of the 
  #   resulting matrix
  row.names(Prob.Y) <- NULL
  
  Prob.Y
}



###############################################################################
# GSE
## Brianna Hitt 11-13-19 - allow Se, Sp values to vary across stages of testing

ProbY.1 <- function(p.vec, g3, sp=1, se=1) {
  g <- length(g3)
  N <- length(p.vec)
  p.ord <- p.vec
  g.ord <- c(0, cumsum(g3))
  Prob.Y <- NULL
  for (i in 1:g) {
    p.g2 <- p.ord[(g.ord[i] + 1):g.ord[i + 1]]
    I.g2 <- length(p.g2)
    for (j in 1:I.g2) {
      # Brianna 09.30.19 - 3-stage testing
      if (I.g2 > 1) {
        # below is Michael Black's original code
        # Prob.Y <- rbind(Prob.Y, se^3)
        
        # Brianna 09.30.19 
        # below allows Se, Sp to vary across stages of testing
        Prob.Y <- rbind(Prob.Y, prod(se[1:3]))
      }
      # Brianna 09.30.19 - 2-stage testing
      else if (I.g2 == 1) {
        # below is Michael Black's original code
        # Prob.Y <- rbind(Prob.Y, se^2)
        
        # Brianna 09.30.19
        # below allows Se, Sp to vary across stages of testing
        Prob.Y <- rbind(Prob.Y, prod(se[1:2]))
      }
    }
  }
  Prob.Y
}

###############################################################################
###############################################################################

#  4 step method for these tests

# 1-Pooling specificity: Is specific for each individual obtained by changing 
#                        the pi to zero one at a time and finding the pooled 
#                       probability of testing positive.
# Could be done on all the positives and negatives to see how they fall out.

# Brianna Hitt - 09.30.19
# Did not edit this function because it is not called anywhere

ProbY.4s <- function(p.vec, vec.g2, vec.g3, sp = 1, se = 1) {
  g2 <- length(vec.g2)
  g3 <- length(vec.g3)
  N <- length(p.vec)
  p.ord <- p.vec
  g2.ord <- c(0,cumsum(vec.g2))
  g3.ord <- c(0,cumsum(vec.g3))
  Prob.Y <- NULL
  p.ord.prob <- prod(1-p.ord)
  for (i in 1:g2) {
    p.g2 <- p.ord[(g3.ord[(g2.ord[i] + 1)] + 1):g3.ord[g2.ord[i + 1] + 1]]
    p.g2.prob <- prod(1 - p.g2)
    for (j in (g2.ord[i] + 1):g2.ord[i + 1]) {
      p.g3 <- p.ord[(g3.ord[j] + 1):g3.ord[j + 1]]
      p.g3.prob <- prod(1 - p.g3)
      use.lgth <- length(p.g3)
      if (vec.g2[i] > 1) {
      if (use.lgth > 1) {
      for (k in 1:use.lgth) {
        p.k <- p.g3[k]
        probYis<-((1 - sp)^4)*p.ord.prob + 
          se*((1 - sp)^3)*(p.g2.prob - p.ord.prob) +
          (se^2)*((1 - sp)^2)*(p.g3.prob - p.g2.prob) +
          (se^3)*(1 - sp)*((1 - p.k) - p.g3.prob) + (se^4)*p.k
        Prob.Y<-rbind(Prob.Y, probYis)
      }
      }
      else if (use.lgth == 1){
        probYis <- ((1 - sp)^3)*p.ord.prob + 
          se*((1 - sp)^2)*(p.g2.prob - p.ord.prob) +
          (se^2)*((1 - sp))*(p.g3.prob - p.g2.prob) +
          (se^3)*(1 - p.g3.prob)
        Prob.Y <- rbind(Prob.Y, probYis)

      }
      }

      else if (vec.g2[i] == 1) {
        if (use.lgth > 1) {
          for (k in 1:use.lgth) {
            p.k <- p.g3[k]
            probYis <- ((1 - sp)^3)*p.ord.prob + 
              se*((1 - sp)^2)*(p.g2.prob - p.ord.prob) +
              (se^2)*(1 - sp)*((1 - p.k) - p.g3.prob) + (se^3)*p.k
            Prob.Y<-rbind(Prob.Y, probYis)
            
          }
        }
        else if (use.lgth == 1) {
          probYis <- ((1 - sp)^2)*p.ord.prob + 
            se*((1 - sp))*(p.g2.prob - p.ord.prob) +
            (se^2)*(1 - p.g3.prob)
          Prob.Y <- rbind(Prob.Y, probYis)
          
        }
      }
    }

  }
  #The following two if statements deal with the case of 4-stage optimal
  #     being a 3-stage configuration
  if (g2 == 1) {
   Prob.Y <- ProbY(p.vec = p.vec, g3 = vec.g3, sp = sp, se = se)
  }
  if (g3 == N) {
   Prob.Y <- ProbY(p.vec = p.vec, g3 = vec.g2, sp = sp, se = se)
  }

  Prob.Y
}


## Brianna Hitt 11-13-19 - allow Se, Sp values to vary across stages of testing

ProbY.4s.0 <- function(p.vec, vec.g2, vec.g3, sp = 1, se = 1) {
  g2 <- length(vec.g2)
  g3 <- length(vec.g3)
  N <- length(p.vec)
  p.ord <- p.vec
  g2.ord <- c(0, cumsum(vec.g2))
  g3.ord <- c(0, cumsum(vec.g3))
  Prob.Y <- NULL
  for (i in 1:g2) {
    p.g2 <- p.ord[(g3.ord[(g2.ord[i] + 1)] + 1):g3.ord[g2.ord[i + 1] + 1]]
    for (j in (g2.ord[i] + 1):g2.ord[i + 1]) {
      p.g3 <- p.ord[(g3.ord[j] + 1):g3.ord[j + 1]]
      use.lgth <- length(p.g3)
      for (k in 1:use.lgth) {
        p.k <- 0
        p.g3.0 <- p.g3
        p.g3.0[k] <- 0
        p.ord.0 <- p.ord
        p.ord.0[(g3.ord[j] + k)] <- 0
        p.g2.0 <- p.ord.0[(g3.ord[(g2.ord[i] + 1)] + 1):g3.ord[g2.ord[i + 1] + 1]]
        p.ord.prob <- prod(1 - p.ord.0)
        p.g2.prob <- prod(1 - p.g2.0)
        p.g3.prob <- prod(1 - p.g3.0)
        if (vec.g2[i] > 1) {
          if (use.lgth > 1) {
            # below is Michael Black's original code
            # probYis <- ((1 - sp)^4)*p.ord.prob + 
            #   se*((1 - sp)^3)*(p.g2.prob - p.ord.prob) +
            #   (se^2)*((1 - sp)^2)*(p.g3.prob - p.g2.prob) +
            #   (se^3)*(1 - sp)*((1 - p.k) - p.g3.prob) + (se^4)*p.k
            
            # below allows Se, Sp to vary across stages of testing
            probYis <- prod(1 - sp[1:4])*p.ord.prob + 
              se[1]*prod(1 - sp[2:4])*(p.g2.prob - p.ord.prob) +
              prod(se[1:2])*prod(1 - sp[3:4])*(p.g3.prob - p.g2.prob) +
              prod(se[1:3])*(1 - sp[4])*((1 - p.k) - p.g3.prob) + 
              prod(se[1:4])*p.k
            Prob.Y<-rbind(Prob.Y, probYis)
          }
          else if (use.lgth == 1){
            # below is Michael Black's original code
            # probYis<-((1 - sp)^3)*p.ord.prob + 
            #   se*((1 - sp)^2)*(p.g2.prob - p.ord.prob) +
            #   (se^2)*((1 - sp))*(p.g3.prob - p.g2.prob) +
            #   (se^3)*(1 - p.g3.prob)
            
            # below allows Se, Sp to vary across stages of testing
            probYis<-prod(1 - sp[1:3])*p.ord.prob + 
              se[1]*prod(1 - sp[2:3])*(p.g2.prob - p.ord.prob) +
              prod(se[1:2])*(1 - sp[3])*(p.g3.prob - p.g2.prob) +
              prod(se[1:3])*(1 - p.g3.prob)
            Prob.Y<-rbind(Prob.Y, probYis)
          }
        }
        
        else if (vec.g2[i] == 1) {
          if (use.lgth > 1) {
            # below is Michael Black's original code
            # probYis <- ((1 - sp)^3)*p.ord.prob + 
            #   se*((1 - sp)^2)*(p.g2.prob - p.ord.prob) +
            #   (se^2)*(1 - sp)*((1 - p.k) - p.g3.prob) + (se^3)*p.k
            
            # below allows Se, Sp to vary across stages of testing
            probYis <- prod(1 - sp[1:3])*p.ord.prob + 
              se[1]*prod(1 - sp[2:3])*(p.g2.prob - p.ord.prob) +
              prod(se[1:2])*(1 - sp[3])*((1 - p.k) - p.g3.prob) + 
              prod(se[1:3])*p.k
            Prob.Y<-rbind(Prob.Y, probYis)
          }
          else if (use.lgth == 1){
            # below is Michael Black's original code
            # probYis<-((1 - sp)^2)*p.ord.prob + 
            #   se*((1 - sp))*(p.g2.prob - p.ord.prob) +
            #   (se^2)*(1 - p.g3.prob)
            
            # below allows Se, Sp to vary across stages of testing
            probYis<-prod(1 - sp[1:2])*p.ord.prob + 
              se[1]*(1 - sp[2])*(p.g2.prob - p.ord.prob) +
              prod(se[1:2])*(1 - p.g3.prob)
            Prob.Y<-rbind(Prob.Y, probYis)
          }
        }
      }
    }
  }
  #The following two if statements deal with the case of 4-stage optimal
  #     being a 3-stage configuration
  if (g2 == 1) {
    # Brianna 10.06.19 - typo in Michael Black's original program
    # Prob.Y.0 <- ProbY.0(p.vec = p.vec, g3 = vec.g3, sp = sp,se = se)
    Prob.Y <- ProbY.0(p.vec = p.vec, g3 = vec.g3, sp = sp,se = se)
  }
  if (g3 == N) {
    Prob.Y <- ProbY.0(p.vec = p.vec, g3 = vec.g2, sp = sp, se = se)
  }
  
  Prob.Y
}

###############################################################################
## Brianna Hitt 11-13-19 - allow Se, Sp values to vary across stages of testing
ProbY.4s.1 <- function(p.vec, vec.g2, vec.g3, sp = 1, se = 1) {
  g2 <- length(vec.g2)
  g3 <- length(vec.g3)
  N <- length(p.vec)
  p.ord <- p.vec
  g2.ord <- c(0, cumsum(vec.g2))
  g3.ord <- c(0, cumsum(vec.g3))
  Prob.Y <- NULL
  for (i in 1:g2) {
    p.g2 <- p.ord[(g3.ord[(g2.ord[i] + 1)] + 1):g3.ord[g2.ord[i + 1] + 1]]
    for (j in (g2.ord[i] + 1):g2.ord[i + 1]) {
      p.g3 <- p.ord[(g3.ord[j] + 1):g3.ord[j + 1]]
      use.lgth <- length(p.g3)
      for (k in 1:use.lgth) {
        p.k <- 1
        p.g3.1 <- p.g3
        p.g3.1[k] <- 1
        p.ord.1 <- p.ord
        p.ord.1[(g3.ord[j] + k)] <- 1
        p.g2.1 <- p.ord.1[(g3.ord[(g2.ord[i] + 1)] + 1):g3.ord[g2.ord[i + 1] + 1]]
        p.ord.prob <- prod(1 - p.ord.1)
        p.g2.prob <- prod(1 - p.g2.1)
        p.g3.prob <- prod(1 - p.g3.1)
        if (vec.g2[i] > 1) {
          if (use.lgth > 1) {
            # below is Michael Black's original code
            # probYis <- ((1 - sp)^4)*p.ord.prob + 
            #   se*((1 - sp)^3)*(p.g2.prob - p.ord.prob) +
            #   (se^2)*((1 - sp)^2)*(p.g3.prob - p.g2.prob) +
            #   (se^3)*(1 - sp)*((1 - p.k) - p.g3.prob) + (se^4)*p.k
            
            # below allows Se, Sp to vary across stages of testing
            probYis <- prod(1 - sp[1:4])*p.ord.prob + 
              se[1]*prod(1 - sp[2:4])*(p.g2.prob - p.ord.prob) +
              prod(se[1:2])*prod(1 - sp[3:4])*(p.g3.prob - p.g2.prob) +
              prod(se[1:3])*(1 - sp[4])*((1 - p.k) - p.g3.prob) + 
              prod(se[1:4])*p.k
            Prob.Y <- rbind(Prob.Y, probYis)
            
          }
          else if (use.lgth == 1) {
            # below is Michael Black's original code
            # probYis <- ((1 - sp)^3)*p.ord.prob + 
            #   se*((1 - sp)^2)*(p.g2.prob - p.ord.prob) +
            #   (se^2)*((1 - sp))*(p.g3.prob - p.g2.prob) +
            #   (se^3)*(1 - p.g3.prob)
            
            # below allows Se, Sp to vary across stages of testing
            probYis <- prod(1 - sp[1:3])*p.ord.prob + 
              se[1]*prod(1 - sp[2:3])*(p.g2.prob - p.ord.prob) +
              prod(se[1:2])*(1 - sp[3])*(p.g3.prob - p.g2.prob) +
              prod(se[1:3])*(1 - p.g3.prob)
            Prob.Y <- rbind(Prob.Y, probYis)
            
          }
        }
        
        else if (vec.g2[i] == 1) {
          if (use.lgth > 1) {
            # below is Michael Black's original code
            # probYis <- ((1 -sp)^3)*p.ord.prob + 
            #   se*((1 - sp)^2)*(p.g2.prob - p.ord.prob) +
            #   (se^2)*(1 - sp)*((1 - p.k) - p.g3.prob) + (se^3)*p.k
            
            # below allows Se, Sp to vary across stages of testing
            probYis <- prod(1 - sp[1:3])*p.ord.prob + 
              se[1]*prod(1 - sp[2:3])*(p.g2.prob - p.ord.prob) +
              prod(se[1:2])*(1 - sp[3])*((1 - p.k) - p.g3.prob) + 
              prod(se[1:3])*p.k
            Prob.Y <- rbind(Prob.Y, probYis)
          }
          else if (use.lgth == 1){
            # below is Michael Black's original code
            # probYis <- ((1 - sp)^2)*p.ord.prob + 
            #   se*((1 - sp))*(p.g2.prob - p.ord.prob) +
            #   (se^2)*(1 - p.g3.prob)
            
            # below allows Se, Sp to vary across stages of testing
            probYis <- prod(1 - sp[1:2])*p.ord.prob + 
              se[1]*(1 - sp[2])*(p.g2.prob - p.ord.prob) +
              prod(se[1:2])*(1 - p.g3.prob)
            Prob.Y <- rbind(Prob.Y, probYis)
          }
        }
      }
    }
  }
  #The following two if statements deal with the case of 4-stage optimal
  #     being a 3-stage configuration
  if (g2 == 1) {
    Prob.Y <- ProbY.1(p.vec = p.vec, g3 = vec.g3, sp = sp, se = se)
  }
  if (g3 == N) {
    Prob.Y <- ProbY.1(p.vec = p.vec, g3 = vec.g2, sp = sp, se = se)
  }
  
  Prob.Y
}
#end accuracy functions



  #start four-stage CRC programs
###############################################################################
#   Author: Michael Black
#   4-step extension of optimization
#     Attempt 1 to program 4 step will try a slow program first then a fast
#     1. need to get a E(T|pvec) for 4-step set up
#         This will need final group sizes and size of previous groups
#     2. g3 will be the final number of groups at the 3rd step
#     3. g2 will be the number of groups at the 2nd step
#     4. Step4 is individual testing and step 1 is the initial group.
#     5. For the final group sizes we can have a vector of length g3 with the
#         number of individuals in each group
#     6. For the g2 groupings the vector of size g2 will have the number of
#         groupings from g3.
#         Ex for an I = 12 we could have g3 = 5, g2 = 3,
#           with vec.g3 = c(3,3,2,2,2), and vec.g2 = c(2,2,1)
#     7. If optimal vec.g2 or vec.g3 is a vector of 1s then 3 step is optimal.
#
###############################################################################

#first we get the ET for 4step.
## Brianna Hitt 11-13-19 - allow Se, Sp values to vary across stages of testing
Exp.T.4step <- function(p, vec.g2, vec.g3, se = 1, sp = 1){
  
  g2 <- length(vec.g2)
  g3 <- length(vec.g3)
  N <- length(p)
  if (g2 == 1) {
    get.Exp.T <- N - 1 #Exp.T.3step(p,grping=vec.g3,sp=sp,se=se)
  }
  if (g2 > 1) {
    if (g2 < g3) {
      if (g3 < N) {
        # below is Michael Black's original code
        # get.Exp.T <- 1 + g2*((1 - sp)*prod(1 - p) + se*(1 - prod(1 - p)))
        
        # below allows Se, Sp to vary across stages of testing
        get.Exp.T <- 1 + g2*((1 - sp[1])*prod(1 - p) + se[1]*(1 - prod(1 - p)))
        for (i in 1:g2) {
          start.g2 <- sum(vec.g3[0:(sum(vec.g2[0:(i - 1)]))]) + 1
          end.g2 <- sum(vec.g3[0:(sum(vec.g2[0:i]))])
          p.group <- p[start.g2:end.g2]
          p.rest <- p[-(start.g2:end.g2)]
          if(vec.g2[i] > 1) {
            # below is Michael Black's original code
            # get.Exp.T <- get.Exp.T + vec.g2[i]*(((1 - sp)^2)*prod(1 - p) +
            #                                       se*(1 - sp)*(prod(1 - p.group)*(1 - prod(1 - p.rest))) +
            #                                       se^2*(1 - prod(1 - p.group)))

            # below allows Se, Sp to vary across stages of testing            
            get.Exp.T <- get.Exp.T + vec.g2[i]*(prod(1 - sp[1:2])*prod(1 - p) +
                                                  se[1]*(1 - sp[2])*(prod(1 - p.group)*(1 - prod(1 - p.rest))) + 
                                                  prod(se[1:2])*(1 - prod(1 - p.group)))
            for (j in (sum(vec.g2[0:(i - 1)]) + 1):(sum(vec.g2[1:i]))) {
              start.g3 <- sum(vec.g3[0:(j - 1)]) + 1
              end.g3 <- sum(vec.g3[1:j])
              p.group3 <- p[start.g3:end.g3];
              p.rest3 <- p.group[-((start.g3 - start.g2 + 1):(end.g3 - start.g2 + 1))]
              if(vec.g3[j] > 1) {
                # below is Michael Black's original code
                # get.Exp.T <- get.Exp.T + vec.g3[j]*(((1 - sp)^3)*prod(1 - p) +
                #                                       se*((1 - sp)^2)*(prod(1 - p.group)*(1 - prod(1 - p.rest))) +
                #                                       (se^2)*(1 - sp)*(prod(1 - p.group3)*(1 - prod(1 - p.rest3))) + 
                #                                       (se^3)*(1 - prod(1 - p.group3)))
                
                # below allows Se, Sp to vary across stages of testing
                get.Exp.T <- get.Exp.T + vec.g3[j]*(prod(1-sp[1:3])*prod(1 - p) +
                                                      se[1]*prod(1 - sp[2:3])*prod(1 - p.group)*(1 - prod(1 - p.rest)) +
                                                      prod(se[1:2])*(1 - sp[3])*prod(1 - p.group3)*(1 - prod(1 - p.rest3)) + 
                                                      prod(se[1:3])*(1 - prod(1 - p.group3)))
              }
              if(vec.g3[j] == 1) {
                get.Exp.T <- get.Exp.T + 0
              }
            }
          }
          if(vec.g2[i] == 1) {
            if(vec.g3[sum(vec.g2[1:i])] > 1) {
              # below is Michael Black's original code
              # get.Exp.T <- get.Exp.T + vec.g3[sum(vec.g2[1:i])]*(((1 - sp)^2)*prod(1 - p) +
              #                                                      se*(1 - sp)*(prod(1 - p.group)*(1 - prod(1 - p.rest))) + 
              #                                                      (se^2)*(1 - prod(1 - p.group)))
              
              # below allows Se, Sp to vary across stages of testing
              get.Exp.T <- get.Exp.T + vec.g3[sum(vec.g2[1:i])]*(prod(1 - sp[1:2])*prod(1 - p) +
                                                                   se[1]*(1 - sp[2])*prod(1 - p.group)*(1 - prod(1 - p.rest)) + 
                                                                   prod(se[1:2])*(1 - prod(1 - p.group)))
            }
            if(vec.g3[sum(vec.g2[1:i])] == 1) {
              get.Exp.T = get.Exp.T + 0
            }
            
          }
        }
      }
    }
    
    get.mc <- get.mc.4s(vec.g2, vec.g3)
    # get.mc restructures the coding variables to variables used in the paper 
    #   and gives the testing number of stages, that is if a 4 stage 
    #   construction is really a 3 stage it returns 3 stage.
    # an example would be if vec.g2 = c(1,3,1) and vec.g3 = c(3,1,1,1,3) this 
    #   is a 3 stage testing pattern.
    # this allows us to get a strictly 4 stage construction for the optimal, 
    #   the code below just gives a bad value if actual stages less than 4
    if (get.mc$S < 4) {
      get.Exp.T <- N - 1
    }
    if (g3 == N) {
      get.Exp.T <- N - 1 #Exp.T.3step(p,grping=vec.g2,sp=sp,se=se)
    }
    # when g2 equals g3 then individual testing is taking place immediately 
    #   so is not a 4 stage situation, returning N-1 insures it is not optimal
    if (g2 == g3) {   
      get.Exp.T <- N - 1 #Exp.T.3step(p,grping=vec.g3,sp=sp,se=se)
    }
  }
  
  list(vec.g2 = vec.g2, vec.g3 = vec.g3, Exp.T = get.Exp.T)
}



####################################################
# This function gives the descent part of steepest descent
#  the movement continues in the steepest direction
#   until improvement stops or we hit a wall

# steep.move <- function(chg.g2, chg.g3, p, N, g3, g2, opt.grps3,
#                   opt.grps2, iter.opt.3s, iter.opt.2s, iter.opt.ET, sp, se) {
# 
#   new.vec.g2 <- iter.opt.2s + chg.g2
#   new.vec.g3 <- iter.opt.3s + chg.g3
#   if (min(new.vec.g2) == 0){
#     return(grps.iter.both(p, N, g3, g2, opt.grps3 = iter.opt.3s,
#                           opt.grps2 = iter.opt.2s, iter.opt.3s, 
#                           iter.opt.2s, iter.opt.ET, sp, se))
#   }   
#   if (min(new.vec.g3) == 0){
#     return(grps.iter.both(p, N, g3, g2, opt.grps3 = iter.opt.3s,
#                           opt.grps2 = iter.opt.2s, iter.opt.3s, 
#                           iter.opt.2s, iter.opt.ET, sp, se))
#   } 
#   else {
#       Opt.chk <- Exp.T.4step(p, new.vec.g2, new.vec.g3, se, sp)
#           if(Opt.chk$Exp.T < iter.opt.ET) {
#             iter.opt.ET <- Opt.chk$Exp.T
#             iter.opt.3s <- Opt.chk$vec.g3
#             iter.opt.2s <- Opt.chk$vec.g2
#             Iter.all <- Opt.chk
#             return(steep.move(chg.g2, chg.g3, p, N, g3, g2, 
#                               opt.grps3 = iter.opt.3s, opt.grps2 = iter.opt.2s, 
#                               iter.opt.3s, iter.opt.2s, iter.opt.ET, sp, se))
#           }
#          else return(grps.iter.both(p, N, g3, g2, opt.grps3 = iter.opt.3s,
#                             opt.grps2 = iter.opt.2s, iter.opt.3s, iter.opt.2s, 
#                             iter.opt.ET, sp, se))
# 
#   }
# }



###############################################################################
# Trying to iterate both vec.g2 and vec.g3 simultaneously.

# grps.iter.both<-function(p, N, g3, g2, opt.grps3, opt.grps2, iter.opt.3s, 
#                          iter.opt.2s, iter.opt.ET, sp, se) {
#   N.all <- length(p)
#   g3.all <- sum(opt.grps2)
#   g2.all <- length(opt.grps2)
#   Iter.all <- list(vec.g2=iter.opt.2s,vec.g3=iter.opt.3s,Exp.T=iter.opt.ET)
#   init.vec.g3 <- opt.grps3
#   init.vec.g2 <- opt.grps2
#  # init.chk<-Exp.T.4step(p,vec.g2=opt.grps2,vec.g3=init.vec.g3,sp,se)
#   init.chk <- both.call.vec.g3(p, N, g3, g2, opt.grps3 = init.vec.g3, 
#                                opt.grps2 = opt.grps2, iter.opt.3s, iter.opt.2s, 
#                                iter.opt.ET, sp, se)
# 
#   chk.chg <- 0
# 
#           if(init.chk$Exp.T < iter.opt.ET) {
#             iter.opt.ET <- init.chk$Exp.T
#             iter.opt.3s <- init.chk$vec.g3
#             iter.opt.2s <- init.chk$vec.g2
#             chk.chg <- 1
#             Iter.all <- init.chk
#           }
# 
#   for (k in 1:g2.all) {
#     for(j in 1:2) {
#       for (i in 1:g2.all) {
#         if(j == 1) {
#           if (init.vec.g2[i] > 1) {
#             init.vec.g2.mod <- init.vec.g2
#             init.vec.g2.mod[i] <- init.vec.g2.mod[i]-1
#             init.vec.g2.mod[k] <- init.vec.g2.mod[k]+1
#             chk.opt <- both.call.vec.g3(p, N, g3, g2, opt.grps3 = init.vec.g3, 
#                                         opt.grps2 = init.vec.g2.mod,
#                                         iter.opt.3s, iter.opt.2s, iter.opt.ET, 
#                                         sp, se)
# 
#           if(chk.opt$Exp.T < iter.opt.ET) {
#             iter.opt.ET <- chk.opt$Exp.T
#             iter.opt.3s <- chk.opt$vec.g3
#             iter.opt.2s <- chk.opt$vec.g2
#             chk.chg <- 1
#             Iter.all <- chk.opt
#           }
# 
#           }
# 
# 
# 
#         }
#         if(j == 2) {
#           if (init.vec.g2[k] > 1) {
#             init.vec.g2.mod <- init.vec.g2
#             init.vec.g2.mod[k] <- init.vec.g2.mod[k]-1
#             init.vec.g2.mod[i] <- init.vec.g2.mod[i]+1
#             chk.opt <- both.call.vec.g3(p, N, g3, g2, opt.grps3 = init.vec.g3, 
#                                         opt.grps2 = init.vec.g2.mod,
#                               iter.opt.3s, iter.opt.2s, iter.opt.ET, sp, se)
# 
#           if(chk.opt$Exp.T < iter.opt.ET) {
#             iter.opt.ET <- chk.opt$Exp.T
#             iter.opt.3s <- chk.opt$vec.g3
#             iter.opt.2s <- chk.opt$vec.g2
#             chk.chg <- 1
#             Iter.all <- chk.opt
#           }
# 
#           }
#         }
#       }
#     }
#   }
#   chg.g2 <- iter.opt.2s - init.vec.g2
#   chg.g3 <- iter.opt.3s - init.vec.g3
#   if (chk.chg == 1) return(steep.move(chg.g2, chg.g3, p, N, g3, g2, 
#                                       opt.grps3 = iter.opt.3s,
#                             opt.grps2 = iter.opt.2s, iter.opt.3s, iter.opt.2s, 
#                             iter.opt.ET, sp, se))
#   if (chk.chg == 0) return(Iter.all)
# }



#################
# iterating the vec.g3 possibilities for each modified vec.g2

#function(p,N,g3,g2,opt.grps3,opt.grps2,fin.opt.3s,fin.opt.2s,fin.opt.ET,sp,se)
# original function call changing fin. to iter. since need to check for this 
#   g2, g3 values
# both.call.vec.g3 <- function(p, N, g3, g2, opt.grps3, opt.grps2, iter.opt.3s, 
#                              iter.opt.2s, iter.opt.ET, sp, se) {
#   N.all <- length(p)
#   g3.all <- sum(opt.grps2)
#   g2.all <- length(opt.grps2)
#   Iter.all <- list(vec.g2 = iter.opt.2s, vec.g3 = iter.opt.3s, 
#                    Exp.T = iter.opt.ET)
#   init.vec.g3 <- opt.grps3
#   init.chk <- Exp.T.4step(p, vec.g2 = opt.grps2, vec.g3 = init.vec.g3, 
#                           sp, se)
#   chk.chg <- 0
#           if(init.chk$Exp.T < iter.opt.ET) {
#           iter.opt.ET <- init.chk$Exp.T
#           iter.opt.3s <- init.chk$vec.g3
#           iter.opt.2s <- init.chk$vec.g2
#           chk.chg <- 1
#           Iter.all <- init.chk
#           }
# 
#   for (k in 1:g3.all) {
#     for (j in 1:2)    {
#       for (i in 1:g3.all) {
#       if (j == 1) {
#         if (init.vec.g3[i] > 1) {
#         init.vec.g3.mod <- init.vec.g3
#         init.vec.g3.mod[i] <- init.vec.g3.mod[i]-1
#         init.vec.g3.mod[k] <- init.vec.g3.mod[k]+1
#         chk.opt <- Exp.T.4step(p, vec.g2 = opt.grps2, 
#                                vec.g3 = init.vec.g3.mod, sp, se)
#           if(chk.opt$Exp.T < iter.opt.ET) {
#           iter.opt.ET <- chk.opt$Exp.T
#           iter.opt.3s <- chk.opt$vec.g3
#           iter.opt.2s <- chk.opt$vec.g2
#           chk.chg <- 1
#           Iter.all <- chk.opt
#           }
#         }
#       }
#       if (j == 2) {
#         if (init.vec.g3[k] > 1) {
#         init.vec.g3.mod <- init.vec.g3
#         init.vec.g3.mod[k] <- init.vec.g3.mod[k] - 1
#         init.vec.g3.mod[i] <- init.vec.g3.mod[i] + 1
#         chk.opt <- Exp.T.4step(p, vec.g2 = opt.grps2, 
#                                vec.g3 = init.vec.g3.mod, sp, se)
#           if(chk.opt$Exp.T < iter.opt.ET) {
#           iter.opt.ET <- chk.opt$Exp.T
#           iter.opt.3s <- chk.opt$vec.g3
#           iter.opt.2s <- chk.opt$vec.g2
#           chk.chg <- 1
#           Iter.all <- chk.opt
#           }
#         }
#       }
# 
#       }
#     }
#   }
#   return(Iter.all)
# }
###############################################################################



# This iterates through the possible g2 given a g3
# grps.iter.size.3s <- function(p, N, g3, fin.opt.3s, fin.opt.2s, fin.opt.ET, 
#                               sp = 1, se = 1, init.config = "hom") {
#   Opt.all <- list(vec.g2 = fin.opt.2s, vec.g3 = fin.opt.3s, Exp.T = fin.opt.ET)
# 
#    for (g2 in 1:(g3 - 1)) {
#   # gets front loaded initial groups
#  if (init.config == "ord" || init.config == "both") {
#   opt.grps3 <- c(N - g3 + 1, rep(1, g3 - 1))
#   opt.grps2 <- c(g3 - g2 + 1, rep(1, g2 - 1))
#   iter.opt.3s <- opt.grps3
#   iter.opt.2s <- opt.grps2
#      Iter.all <- Exp.T.4step(p, vec.g2 = opt.grps2, vec.g3 = opt.grps3, se, sp)
#      iter.opt.ET <- Iter.all$Exp.T
# 
#      Opt.sizes <- grps.iter.both(p, N, g3, g2 = g2, opt.grps3, opt.grps2, 
#                                  iter.opt.3s, iter.opt.2s, iter.opt.ET, sp, se)
#      Opt.ET <- Opt.sizes$Exp.T
#      if (Opt.all$Exp.T > Opt.sizes$Exp.T) {
#         Opt.all <- Opt.sizes
#      }
#  }
#   #### Gets even sized groups. is the default method
#  if (init.config == "hom" || init.config == "both") {
#   if (N%%g3 == 0) {
#     opt.grps3 <- c(rep(N%/%g3,g3))
#     }
#   else  {
#     opt.grps3 <- c(rep(N%/%g3 + 1, N%%g3), rep(N%/%g3, g3 - N%%g3))
#     }
#   if  (g3%%g2 == 0) {
#     opt.grps2 <- c(rep(g3%/%g2, g2))
#     }
#   else {
#     opt.grps2<- c(rep(g3%/%g2 + 1, g3%%g2), rep(g3%/%g2, g2 - g3%%g2))
#     }
# 
#     iter.opt.3s <- opt.grps3
#     iter.opt.2s <- opt.grps2
#      Iter.all <- Exp.T.4step(p, vec.g2 = opt.grps2, vec.g3 = opt.grps3, se, sp)
#      iter.opt.ET <- Iter.all$Exp.T
# 
#      Opt.sizes <- grps.iter.both(p, N, g3, g2 = g2, opt.grps3, opt.grps2, 
#                                  iter.opt.3s, iter.opt.2s, iter.opt.ET, sp, se)
#      Opt.ET <- Opt.sizes$Exp.T
#      if (Opt.all$Exp.T > Opt.sizes$Exp.T) {
#         Opt.all <- Opt.sizes
#      }
# 
#         fin.opt.ET <- Opt.all$Exp.T
#         fin.opt.3s <- Opt.all$vec.g3
#         fin.opt.2s <- Opt.all$vec.g2
# 
#    }
#   #could possibly add an early break here don't know if it is justified
#   }
#         fin.opt.ET <- Opt.all$Exp.T
#         fin.opt.3s <- Opt.all$vec.g3
#         fin.opt.2s <- Opt.all$vec.g2
#   Opt.all
# }



#this program iterates through the possible g3.
# Opt.grps.size_number.4step <- function(p, sp = 1, se = 1, init.config = "hom"){
#   start.time <- proc.time()
#   N <- length(p)
#   use.N <- N - floor((N - 4)/2)
#   fin.opt.ET <- 2*N
#   fin.opt.3s <- c(N)
#   fin.opt.2s <- c(1)
#   Opt.all <- NULL
#   for (g3 in 2:use.N) {
#      Opt.sizes2 <- grps.iter.size.3s(p, N, g3 = g3, fin.opt.3s, fin.opt.2s, 
#                                      fin.opt.ET, sp, se, 
#                                      init.config = init.config)
#      Opt.ET2 <- Opt.sizes2$Exp.T
#      if (g3 > 5) {
#       if (Opt.sizes2$Exp.T >= fin.opt.ET) break
#      }
#      if (fin.opt.ET > Opt.sizes2$Exp.T) {
#         fin.opt.ET <- Opt.sizes2$Exp.T
#         fin.opt.3s <- Opt.sizes2$vec.g3
#         fin.opt.2s <- Opt.sizes2$vec.g2
#         Opt.all <- Opt.sizes2
#      }
#         fin.opt.ET <- Opt.all$Exp.T
#         fin.opt.3s <- Opt.all$vec.g3
#         fin.opt.2s <- Opt.all$vec.g2
# 
#   }
# 
#   end.time <- proc.time()
#   save.time <- end.time-start.time
#   print("\n Number of minutes running for I = ",N,":", save.time[3]/60, "\n \n")
# 
# 
#   Opt.all
# }
  #End four-stage CRC functions



  #start four stage ORC functions
###############################################################################
#
#   4-step extension of optimization
#     Attempt 1 to program 4 step will try a slow program first then a fast
#     1. need to get a E(T|pvec) for 4-step set up
#         This will need final group sizes and size of previous groups
#     2. g3 will be the final number of groups at the 3rd step
#     3. g2 will be the number of groups at the 2nd step
#     4. Step4 is individual testing and step 1 is the initial group.
#     5. For the final group sizes we can have a vector of length g3 with the
#         number of individuals in each group
#     6. For the g2 groupings the vector of size g2 will have the number of
#         groupings from g3.
#         Ex for an I = 12 we could have g3 = 5, g2 = 3,
#           with vec.g3 = c(3,3,2,2,2), and vec.g2 = c(2,2,1)
#     7. If optimal vec.g2 is a vector of 1s then 3 step is optimal.
#
###############################################################################

#first we get the ET for 4step.
# Exp.T.4step.slow<-function(p,vec.g2,vec.g3,se=1,sp=1){
#   g2<-length(vec.g2)
#   g3<-length(vec.g3)
#   N<-length(p)
#   if (g2==1){
#     get.Exp.T<-N-1#Exp.T.3step(p,grping=vec.g3,sp=sp,se=se)
#   }
#   if (g2>1) {
#   if (g2<g3) {
#   if (g3 < N) {
#   get.Exp.T<-1+g2*((1-sp)*prod(1-p)+se*(1-prod(1-p)))
#   for (i in 1:g2) {
#      start.g2<-sum(vec.g3[0:(sum(vec.g2[0:(i-1)]))])+1
#      end.g2<- sum(vec.g3[0:(sum(vec.g2[0:i]))])
#      p.group<-p[start.g2:end.g2]
#      p.rest<-p[-(start.g2:end.g2)]
#      if(vec.g2[i]>1) {
#        get.Exp.T<-get.Exp.T + vec.g2[i]*(((1-sp)^2)*prod(1-p) + 
#                                            se*(1-sp)*(prod(1-p.group)*(1-prod(1-p.rest))) + 
#                                            (se^2)*(1 - prod(1-p.group)))
#        for (j in (sum(vec.g2[0:(i-1)])+1):(sum(vec.g2[1:i]))) {
#           start.g3<-sum(vec.g3[0:(j-1)])+1
#           end.g3<- sum(vec.g3[1:j])
#           p.group3<-p[start.g3:end.g3];
#           p.rest3<-p.group[-((start.g3-start.g2+1):(end.g3-start.g2+1))]
#           if(vec.g3[j]>1) {
#             get.Exp.T<-get.Exp.T + vec.g3[j]*(((1-sp)^3)*prod(1-p) +
#                                                 se*((1-sp)^2)*(prod(1-p.group)*(1-prod(1-p.rest))) +
#                                                 (se^2)*(1-sp)*(prod(1-p.group3)*(1-prod(1-p.rest3))) + 
#                                                 (se^3)*(1 - prod(1-p.group3)))
#           }
#           if(vec.g3[j]==1) {
#             get.Exp.T<-get.Exp.T+0
#           }
#         }
#      }
#      if(vec.g2[i]==1) {
#           if(vec.g3[sum(vec.g2[1:i])] > 1) {
#             get.Exp.T=get.Exp.T+vec.g3[sum(vec.g2[1:i])]*(((1-sp)^2)*prod(1-p)+
#                 se*(1-sp)*(prod(1-p.group)*(1-prod(1-p.rest)))+(se^2)*(1 - prod(1-p.group)))
#           }
#           if(vec.g3[sum(vec.g2[1:i])]==1) {
#             get.Exp.T=get.Exp.T + 0
#           }
# 
#      }
#   }
#   }
#   }
# 
#   get.mc <- get.mc.4s(vec.g2,vec.g3)
#   # get.mc restructures the coding variables to variables used in the paper 
#   #   and gives the testing number of stages, that is if a 4 stage construction 
#   #   is really a 3 stage it returns 3 stage.
#   # an example would be if vec.g2 = c(1,3,1) and vec.g3 = c(3,1,1,1,3) 
#   # this is a 3 stage testing pattern.
#   if (get.mc$S < 4) {
#     get.Exp.T<- N-1
#   }
#   if (g3 == N) {
#     get.Exp.T<-N-1#Exp.T.3step(p,grping=vec.g2,sp=sp,se=se)
#   }
#   # when g2 equals g3 then individual testing is taking place immediately so is 
#   #   not a 4 stage situation, returning N-1 insures it is not optimal
#   if (g2 == g3) {   
#     get.Exp.T<-N-1#Exp.T.3step(p,grping=vec.g3,sp=sp,se=se)
#   }
#   }
#   list(vec.g2=vec.g2,vec.g3=vec.g3,Exp.T=get.Exp.T)
# }



# Iterates through all the possible combinations for vec.g3 given a g3 and 
#   vec.g2 this will be replaced with faster
# grps.iter.4step.slow<-function(p,N,g3,g2,opt.grps3,opt.grps2,fin.opt.3s,
#                                fin.opt.2s,fin.opt.ET,sp,se) {
# 
#   N.all<-length(p)
#   g3.all<-sum(opt.grps2)
#   g2.all<-length(opt.grps2)
#   grps.start<-1
#   grps.fin<-N-g3+1
#   Opt.all<-list(vec.g2=fin.opt.2s,vec.g3=fin.opt.3s,Exp.T=fin.opt.ET)
#   if (g3==1) {
# 
#       opt.grps3[g3.all-g3+1]<-N.all-sum(opt.grps3[1:(g3.all-1)])
#       opt.grps.ET<-Exp.T.4step.slow(p,vec.g2=opt.grps2,vec.g3=opt.grps3,sp,se)
#     if (fin.opt.ET >= opt.grps.ET$Exp.T) {
#       fin.opt.3s<-opt.grps.ET$vec.g3
#       fin.opt.2s<-opt.grps.ET$vec.g2
#       fin.opt.ET<-opt.grps.ET$Exp.T
#       Opt.all<-opt.grps.ET
#     }
#         fin.opt.ET<-Opt.all$Exp.T
#         fin.opt.3s<-Opt.all$vec.g3
#         fin.opt.2s<-Opt.all$vec.g2
# 
#   }
#   if (g3>1) {
#     for (i in grps.start:grps.fin) {
#       opt.grps3[g3.all-g3+1]<-i
#       Opt.all<-grps.iter.4step.slow(p,N=N-i,g3=g3-1,g2,opt.grps3,opt.grps2,
#                                     fin.opt.3s,fin.opt.2s,fin.opt.ET,sp,se)
#         fin.opt.ET<-Opt.all$Exp.T
#         fin.opt.3s<-Opt.all$vec.g3
#         fin.opt.2s<-Opt.all$vec.g2
#     }
#   }
#         fin.opt.ET<-Opt.all$Exp.T
#         fin.opt.3s<-Opt.all$vec.g3
#         fin.opt.2s<-Opt.all$vec.g2
# 
#   Opt.all
# }



#This iterates through all the possible vec.g2 combinations given a g3

# grps.iter.4s.slow<-function(p,N,g3,g2,opt.grps3,opt.grps2,fin.opt.3s,
#                             fin.opt.2s,fin.opt.ET,sp=1,se=1){
#   g3.all<-length(opt.grps3)
#   g2.all<-length(opt.grps2)
#   Opt.all<-list(vec.g2=fin.opt.2s,vec.g3=fin.opt.3s,Exp.T=fin.opt.ET)
#   grps.start<-1
#   grps.fin<-g3-g2+1
#     if (g2==1) {
#       opt.grps2[g2.all-g2+1]<-g3.all-sum(opt.grps2[1:(g2.all-1)])
#       opt.grps<-grps.iter.4step.slow(p,N,g3=g3.all,g2=g2.all,opt.grps3,
#                                      opt.grps2,fin.opt.3s,fin.opt.2s,fin.opt.ET,
#                                      sp,se)
#     if (fin.opt.ET >= opt.grps$Exp.T) {
#       fin.opt.3s<-opt.grps$vec.g3
#       fin.opt.2s<-opt.grps$vec.g2
#       fin.opt.ET<-opt.grps$Exp.T
#       Opt.all<-opt.grps
#     }
#         fin.opt.ET<-Opt.all$Exp.T
#         fin.opt.3s<-Opt.all$vec.g3
#         fin.opt.2s<-Opt.all$vec.g2
#   }
#   else if (g2>1) {
#     for (i in grps.start:grps.fin) {
#       opt.grps2[g2.all-g2+1]<-i
#       Opt.all<-grps.iter.4s.slow(p,N,g3=g3-i,g2=g2-1,opt.grps3,opt.grps2,
#                                  fin.opt.3s,fin.opt.2s,fin.opt.ET,sp,se)
#         fin.opt.ET<-Opt.all$Exp.T
#         fin.opt.3s<-Opt.all$vec.g3
#         fin.opt.2s<-Opt.all$vec.g2
#     }
#   }
#         fin.opt.ET<-Opt.all$Exp.T
#         fin.opt.3s<-Opt.all$vec.g3
#         fin.opt.2s<-Opt.all$vec.g2
#   Opt.all
# }



# This iterates through the possible g2 given a g3
# grps.iter.size.4s.slow<-function(p,N,g3,fin.opt.3s,fin.opt.2s,
#                                  fin.opt.ET,sp=1,se=1) {
#   Opt.all<-list(vec.g2=fin.opt.2s,vec.g3=fin.opt.3s,Exp.T=fin.opt.ET)
#   opt.grps3<-rep(0,g3)
#   use.g3<-min(g3,7)
#    for (g2 in 1:g3) {
#      opt.grps2<-rep(0,g2)
#      Opt.sizes<-grps.iter.4s.slow(p,N,g3,g2=g2,opt.grps3,opt.grps2,fin.opt.3s,
#                                   fin.opt.2s,fin.opt.ET,sp,se)
#      Opt.ET<-Opt.sizes$Exp.T
#      if (fin.opt.ET>=Opt.sizes$Exp.T) {
#         fin.opt.2s<-Opt.sizes$vec.g2
#         fin.opt.3s<-Opt.sizes$vec.g3
#         fin.opt.ET<-Opt.sizes$Exp.T
#         Opt.all<-Opt.sizes
#      }
#         fin.opt.ET<-Opt.all$Exp.T
#         fin.opt.3s<-Opt.all$vec.g3
#         fin.opt.2s<-Opt.all$vec.g2
# 
#   }
#         fin.opt.ET<-Opt.all$Exp.T
#         fin.opt.3s<-Opt.all$vec.g3
#         fin.opt.2s<-Opt.all$vec.g2
#   Opt.all
# }



#this program iterates through the possible g3.
# Opt.grps.size_number.4step.slow<-function(p,sp=1,se=1) {
#   start.time<-proc.time()
#   N<-length(p)
#   use.N<-min(N,10)
#   fin.opt.ET<-2*N
#   fin.opt.3s<-c(N)
#   fin.opt.2s<-c(1)
#   Opt.all<-NULL
#   for (g3 in 1:use.N) {
#      Opt.sizes2<-grps.iter.size.4s.slow(p,N,g3=g3,fin.opt.3s,fin.opt.2s,
#                                         fin.opt.ET,sp,se)
#      Opt.ET2<-Opt.sizes2$Exp.T
#      if (fin.opt.ET>=Opt.sizes2$Exp.T) {
#         fin.opt.ET<-Opt.sizes2$Exp.T
#         fin.opt.3s<-Opt.sizes2$vec.g3
#         fin.opt.2s<-Opt.sizes2$vec.g2
#         Opt.all<-Opt.sizes2
#      }
#         fin.opt.ET<-Opt.all$Exp.T
#         fin.opt.3s<-Opt.all$vec.g3
#         fin.opt.2s<-Opt.all$vec.g2
# 
#   }
# 
#   end.time<-proc.time()
#   save.time<-end.time-start.time
#   print("\n Number of minutes running for I = ",N,":", save.time[3]/60, "\n \n")
# 
# 
#   Opt.all
# }
  #end four-stage ORC functions



  #start notation adjustment functions
###############################################################################
#   Author: Michael Black
###############################################################################
get.mc.3s <- function(grping){
  c2 <- length(grping)
  N <- sum(grping)
  g.ord <- c(0, cumsum(grping))
  vec.g2.start <- g.ord[1:c2]
  vec.g2.end <- g.ord[2:(c2 + 1)]
  vec.I2 <- vec.g2.end - vec.g2.start
  vec.m2 <- grping
  S <- 3

    for (i in 1:c2) {
      if (vec.m2[i] == 1) vec.m2[i] = 0
    }

  if (c2 == 1) {        #this is the case for 2-stage testing
    S <- 2
    vec.m2 <- NULL
    vec.I2 <- NULL
    c2 <- N
  }
  if (c2 == N) {        #this is the case for 2-stage testing
    S <- 2
    vec.m2 <- NULL
    vec.I2 <- NULL
  }


list(S = S, c2 = c2, vec.m2 = vec.m2, vec.I2 = vec.I2)

}



#first we get the ET for 4step.
get.mc.4s <- function(vec.g2, vec.g3){
  g2 <- length(vec.g2)
  g3 <- length(vec.g3)
  N <- sum(vec.g3)
  S  <- 4
  c2 <- g2
  c3 <- g3
  vec.m2 <- vec.g2
  vec.m3 <- vec.g3
    g2.ord <- c(0, cumsum(vec.g2))
    g3.ord <- c(0, cumsum(vec.g3))
    vec.I2 <- g3.ord[(g2.ord[2:(c2 + 1)] + 1)] - g3.ord[(g2.ord[1:c2] + 1)]
    vec.I3 <- g3.ord[2:(c3 + 1)] - g3.ord[1:c3]
    condition.1 <- c2
    condition.2 <- c3
    if (condition.1 == 1){           # signifies 3 stages or less
      get.mc <- get.mc.3s(grping = vec.g3)
      S <- get.mc$S
      c2 <- get.mc$c2
      c3 <- NULL
      vec.m2 <- get.mc$vec.m2
      vec.I2 <- get.mc$vec.I2
      vec.m3 <- NULL
      vec.I3 <- NULL

    }
    if (condition.2 == N) {          #signifies 3 stages or less
      get.mc <- get.mc.3s(grping = vec.g2)
      S <- get.mc$S
      c2 <- get.mc$c2
      c3 <-NULL
      vec.m2 <- get.mc$vec.m2
      vec.I2 <- get.mc$vec.I2
      vec.m3 <- NULL
      vec.I3 <- NULL
    }

  if (condition.1 > 1) {

  if (g3 < N) {

  for (i in condition.1:1) {
    # moves testing up a stage if it was a place holder
     if (vec.g2[i] == 1) {           

        condition.3 <- vec.I3[g2.ord[i + 1]]
        if (condition.3 > 1) {
          vec.I2[i] <- vec.I3[g2.ord[i + 1]]
          vec.m2[i] <- vec.m3[g2.ord[i + 1]]

          if (g2.ord[i + 1] < c3) {
            vec.m3 <- c(vec.m3[0:(g2.ord[i + 1] - 1)], 
                        rep(1, vec.I3[g2.ord[i + 1]]), 
                        vec.m3[(g2.ord[i + 1] + 1):c3])
            vec.I3 <- c(vec.I3[0:(g2.ord[i + 1] - 1)], 
                        rep(1, vec.I3[g2.ord[i + 1]]), 
                        vec.I3[(g2.ord[i + 1] + 1):c3])
          }
          if (g2.ord[i + 1] == c3) {
            vec.m3 <- c(vec.m3[1:(g2.ord[i + 1] - 1)], 
                        rep(1, vec.I3[g2.ord[i + 1]]))
            vec.I3 <- c(vec.I3[1:(g2.ord[i + 1] - 1)], 
                        rep(1, vec.I3[g2.ord[i + 1]]))
          }
        }
        #eliminates space for already tested individual
        if (condition.3 == 1) {   

           vec.I3 <- vec.I3[-g2.ord[i + 1]]
           vec.m3 <- vec.m3[-g2.ord[i + 1]]
           vec.m2[i] = 0
        }
        if (vec.m2[i] == 1) vec.m2[i] = 0
     }
   c3 <- length(vec.I3)
  }

  c3 <- length(vec.I3)

  for (j in 1:c3)   {
     if (vec.m3[j] == 1) vec.m3[j] = 0
  }

  }
  }
  if(!is.null(vec.m3)) {
  if(sum(vec.m3) == 0) {
    vec.m3 <- NULL
    vec.I3 <- NULL
    S <- 3
  }
  }
  list(S = S, c2 = c2, m11 = c2, c3 = c3, vec.g2 = vec.g2, vec.g3 = vec.g3,
        vec.m2 = vec.m2, vec.m3 = vec.m3, I11 = N, vec.I2 = vec.I2,
        vec.I3 = vec.I3)
}
  #end notation adjustment functions



  #start hierarchical.desc() functions
###############################################################################
#    Michael Black 03 12 2013
#     Function: hierarchical.desc gets descriptive/diagnostic info given a p 
#               vector and vectors for subgroup sizes for 2, 3 or 4 stages
#       calls:
#       inputs: p vector of probabilities,
#               se and sp sensitivity and specificity,,
#               g2 vector of number of subgroups stage 2 subgroups split into
#               g3 vector of number of subgroups stage 3 subgroups split into
#               m?  if we add more stages we will need more m vectors
#               order.p = TRUE if False p vector not sorted
#
#        note: if g2 is null, then stages = 2
#              if g3 is null but g2 has values then stages = 3
#              if g2 and g3 both have values then stages =4
#              g vectors should be entered using notation that keeps track of 
#                 all individuals through all stages (eg if a group of 10 
#                 splits into 5, 4 and 1 individual at stage 2 then into 
#                 3, 2, 2, 1 and 1 individuals at stage 3
#                 before individual testing at stage 4, then g2 = c(2,3,1) and 
#                 g3 = c(3,2,2,1,1,1) the one that tested individually at 
#                 stage 2 is still numbered at stage 3,
#                 even though it won't be tested again
#
###################################################################
# function called by hierarchical.desc() to convert I2 and I3 to the form 
#   used in the programs
get.g2_g3 <- function(I2, I3) {
  extra.g3  <- ifelse(I2 == 1, 1, 0)
  length.g3 <- length(I3) + sum(extra.g3)
  g2 <- rep(0, length(I2))

  for (i in 1:length(I2)) {
    if (I2[i] == 1)   I3 <- c(I3[0:(j - 1)], 1, I3[j:length(I3)])
    for (j in 1:length(I3)) {
      if (cumsum(I3)[j] == cumsum(I2)[i]) {
        g2[i] = j
      }
    }
  }
  g2 <- g2 - c(0, g2[1:(length(g2) - 1)])
  g3 <- I3
  list(g2 = g2, g3 = g3)
}



##Note Mike Black 3/22/13: I've added the call to the function with I2 and I3 
#  as In section 3.2,

# not used due to creation of hierarchical.desc2() with suppresses printing of stages 

# hierarchical.desc <- function(p, I2 = NULL, I3 = NULL, se = 1, sp = 1, 
#                               order.p = TRUE) {
#   if (order.p == TRUE) p = sort(p)
#   N <- length(p)
# 
#   g2 <- NULL
#   g3 <- NULL
# 
#   #converts I2 to what is used in program
#   if (!is.null(I2) && is.null(I3)) g2 <- I2 
# 
#   if (!is.null(I2) && !is.null(I3)) {
#     #converts I2, I3 to what is used in the programs
#     get.g.vec <- get.g2_g3(I2 = I2, I3 = I3) 
#     g2 <- get.g.vec$g2
#     g3 <- get.g.vec$g3
#   }
# 
#   if (is.null(g2) || (is.null(g3) && max(g2) == 1 || max(g2) == N)) {
#     print("Two stage procedure")
#     stages = 2
# 
#     g2 <- "individual testing"
#     g3 <- NULL
# 
#     et <- 1 + N*(se*(1 - prod(1 - p)) + (1 - sp)*prod(1 - p))
# 
#     pse.vec  <- ProbY.1(p.vec = p, g3 = rep(1, N), sp = sp, se = se)
#     psp.vec  <- 1 - ProbY.0(p.vec = p, g3 = rep(1,N), sp = sp, se = se)
#     pppv.vec <- p*pse.vec/(p*pse.vec + (1 - p)*(1 - psp.vec))
#     pnpv.vec <- (1 - p)*psp.vec/((1 - p)*psp.vec + p*(1 - pse.vec))
# 
#     pse  <- sum(p*pse.vec)/sum(p)
#     psp  <- sum((1 - p)*psp.vec)/sum(1 - p)
#     pppv <- sum(p*pse.vec)/sum((p*pse.vec + (1 - p)*(1 - psp.vec)))
#     pnpv <- sum((1 - p)*psp.vec)/sum(((1 - p)*psp.vec + p*(1 - pse.vec)))
# 
#     I2 <- "individual testing"
#     I3 <- NULL
# 
#     m1 <- N
#     m2 <- "individual testing"
#     m3 <- NULL
# 
#   }
#   else if (!is.null(g2) && 
#            (is.null(g3) || length(g3)==length(g2) || max(g3)==1)){
#     print("Three stage procedure")
#     stages = 3
# 
#     g3 <- "individual testing"
# 
#     et <- Exp.T.3step(p = p, grping = g2, se = se, sp = sp)
# 
#     pse.vec  <- ProbY.1(p.vec = p,g3 = g2,sp = sp,se = se)
#     psp.vec  <- 1 - ProbY.0(p.vec = p,g3 = g2,sp = sp,se = se)
#     pppv.vec <- p*pse.vec/(p*pse.vec + (1 - p)*(1 - psp.vec))
#     pnpv.vec <- (1 - p)*psp.vec/((1 - p)*psp.vec + p*(1 - pse.vec))
# 
#     pse  <- sum(p*pse.vec)/sum(p)
#     psp  <- sum((1 - p)*psp.vec)/sum(1 - p)
#     pppv <- sum(p*pse.vec)/sum((p*pse.vec + (1 - p)*(1 - psp.vec)))
#     pnpv <- sum((1 - p)*psp.vec)/sum(((1 - p)*psp.vec + p*(1 - pse.vec)))
# 
#     test.pattern <- get.mc.3s(grping = g2)
# 
#     I2 <- test.pattern$vec.I2
#     I3 <- "individual testing"
# 
#     m2 = test.pattern$vec.m2
#     m1 = length(m2)
#     m3 = "individual testing"
# 
# 
#   }
#   else if (!is.null(g2) && !is.null(g3)) {
#     print("Four stage procedure")
#     stages = 4
#     et <- Exp.T.4step(p = p, vec.g2 = g2, vec.g3 = g3, se = se, sp = sp)$Exp.T
# 
#     pse.vec  <- ProbY.4s.1(p.vec = p, vec.g2 = g2, vec.g3 = g3, 
#                            sp = sp, se = se)
#     psp.vec  <- 1 - ProbY.4s.0(p.vec = p, vec.g2 = g2, vec.g3 = g3, 
#                                sp = sp, se = se)
#     pppv.vec <- p*pse.vec/(p*pse.vec + (1 - p)*(1 - psp.vec))
#     pnpv.vec <- (1 - p)*psp.vec/((1 - p)*psp.vec + p*(1 - pse.vec))
# 
#     pse  <- sum(p*pse.vec)/sum(p)
#     psp  <- sum((1 - p)*psp.vec)/sum(1 - p)
#     pppv <- sum(p*pse.vec)/sum((p*pse.vec + (1 - p)*(1 - psp.vec)))
#     pnpv <- sum((1 - p)*psp.vec)/sum(((1 - p)*psp.vec + p*(1 - pse.vec)))
# 
#     test.pattern <- get.mc.4s(vec.g2 = g2, vec.g3 = g3)
# 
#     I2 <- test.pattern$vec.I2
#     I3 <- test.pattern$vec.I3
# 
#     m2 = test.pattern$vec.m2
#     m1 = length(m2)
#     m3 = test.pattern$vec.m3
# 
#   }
# 
# 
#   c("PSe.individual", "PSp.individual", "PPPV.individual", "PNPV.individual")
# 
# 
#   list(ET = et, stages = stages, group.size = N, I2 = I2, I3 = I3, m1 = m1, 
#        m2 = m2, m3 = m3, 
#        individual.testerror = data.frame(p, pse.vec, psp.vec, pppv.vec, 
#                                          pnpv.vec, row.names=NULL),
#        group.testerror = c(PSe = pse, PSp = psp, PPPV = pppv, PNPV = pnpv),
#        individual.probabilities = p)
# 
# }


# start hierarchical.desc2() function
###############################################################################
## Chris 5-16-14 - turn off printing of number of stages
## Brianna Hitt 4-11-16 - turn off printing for two, three, and four stages
## Brianna Hitt 11-13-19 - allow Se, Sp values to vary across stages of testing
##                         Users may specify a single Se, Sp value if 
##                         the value is assumed to be the same for all stages

hierarchical.desc2 <- function(p, I2 = NULL, I3 = NULL, se = 1, sp = 1, 
                               order.p = TRUE) {
  if (order.p == TRUE) p = sort(p)
  N <- length(p)
  
  g2 <- NULL
  g3 <- NULL
  
  # Brianna 09.24.19 - stage 2 subgroup sizes
  # converts I2 to what is used in program
  if (!is.null(I2) && is.null(I3)) g2 <- I2 
  
  # Brianna 09.24.19 - stage 2 and 3 subgroup sizes 
  # g2 is the number of subgroups that will result if a second stage group 
  #   tests positive
  if (!is.null(I2) && !is.null(I3)) {
    #converts I2, I3 to what is used in the programs
    get.g.vec <- get.g2_g3(I2 = I2, I3 = I3) 
    g2 <- get.g.vec$g2
    g3 <- get.g.vec$g3
  }
  
  # Brianna 09.24.19 - two-stage testing procedure
  if (is.null(g2) || (is.null(g3) && max(g2) == 1 || max(g2) == N)) {
    #print("Two stage procedure")
    stages = 2
    
    g2 <- "individual testing"
    g3 <- NULL
    
    # Brianna 09.30.19 - need to specify the first element of se/sp
    et <- 1 + N*(se[1]*(1 - prod(1 - p)) + (1 - sp[1])*prod(1 - p))
    
    # Brianna 09.30.19 - edited ProbY.1() and ProbY.0() for se/sp
    pse.vec  <- ProbY.1(p.vec = p, g3 = rep(1, N), sp = sp, se = se)
    psp.vec  <- 1 - ProbY.0(p.vec = p, g3 = rep(1,N), sp = sp, se = se)
    
    # Brianna 09.30.19 - the predictive values do not change
    pppv.vec <- p*pse.vec/(p*pse.vec + (1 - p)*(1 - psp.vec))
    pnpv.vec <- (1 - p)*psp.vec/((1 - p)*psp.vec + p*(1 - pse.vec))
    
    # Brianna 09.30.19 - the overall accuracy measures do not change
    pse  <- sum(p*pse.vec)/sum(p)
    psp  <- sum((1 - p)*psp.vec)/sum(1 - p)
    pppv <- sum(p*pse.vec)/sum((p*pse.vec + (1 - p)*(1 - psp.vec)))
    pnpv <- sum((1 - p)*psp.vec)/sum(((1 - p)*psp.vec + p*(1 - pse.vec)))
    
    I2 <- "individual testing"
    I3 <- NULL
    
    m1 <- N
    m2 <- "individual testing"
    m3 <- NULL
    
  }
  
  # Brianna 09.30.19 - three-stage testing procedure
  else if (!is.null(g2) && 
           (is.null(g3) || length(g3) == length(g2) || max(g3) == 1)) {
    # print("Three stage procedure")
    stages = 3
    
    g3 <- "individual testing"
    
    et <- Exp.T.3step(p = p, grping = g2, se = se, sp = sp)
    
    pse.vec  <- ProbY.1(p.vec = p,g3 = g2,sp = sp,se = se)
    psp.vec  <- 1 - ProbY.0(p.vec = p,g3 = g2,sp = sp,se = se)
    pppv.vec <- p*pse.vec/(p*pse.vec + (1 - p)*(1 - psp.vec))
    pnpv.vec <- (1 - p)*psp.vec/((1 - p)*psp.vec + p*(1 - pse.vec))
    
    pse  <- sum(p*pse.vec)/sum(p)
    psp  <- sum((1 - p)*psp.vec)/sum(1 - p)
    pppv <- sum(p*pse.vec)/sum((p*pse.vec + (1 - p)*(1 - psp.vec)))
    pnpv <- sum((1 - p)*psp.vec)/sum(((1 - p)*psp.vec + p*(1 - pse.vec)))
    
    test.pattern <- get.mc.3s(grping = g2)
    
    I2 <- test.pattern$vec.I2
    I3 <- "individual testing"
    
    m2 = test.pattern$vec.m2
    m1 = length(m2)
    m3 = "individual testing"
    
    
  }
  else if (!is.null(g2) && !is.null(g3)) {
    # print("Four stage procedure")
    stages = 4
    et <- Exp.T.4step(p = p, vec.g2 = g2, vec.g3 = g3, se = se, sp = sp)$Exp.T
    
    pse.vec  <- ProbY.4s.1(p.vec = p, vec.g2 = g2, vec.g3 = g3, 
                           sp = sp, se = se)
    psp.vec  <- 1 - ProbY.4s.0(p.vec = p, vec.g2 = g2, vec.g3 = g3, 
                               sp = sp, se = se)
    pppv.vec <- p*pse.vec/(p*pse.vec + (1 - p)*(1 - psp.vec))
    pnpv.vec <- (1 - p)*psp.vec/((1 - p)*psp.vec + p*(1 - pse.vec))
    
    pse  <- sum(p*pse.vec)/sum(p)
    psp  <- sum((1 - p)*psp.vec)/sum(1 - p)
    pppv <- sum(p*pse.vec)/sum((p*pse.vec + (1 - p)*(1 - psp.vec)))
    pnpv <- sum((1 - p)*psp.vec)/sum(((1 - p)*psp.vec + p*(1 - pse.vec)))
    
    test.pattern <- get.mc.4s(vec.g2 = g2, vec.g3 = g3)
    
    I2 <- test.pattern$vec.I2
    I3 <- test.pattern$vec.I3
    
    m2 = test.pattern$vec.m2
    m1 = length(m2)
    m3 = test.pattern$vec.m3
    
  }
  
  
  c("PSe.individual", "PSp.individual", "PPPV.individual", "PNPV.individual")
  
  
  list(ET = et, stages = stages, group.size = N, 
       I2 = I2, I3 = I3, m1 = m1, m2 = m2, m3 = m3,
       individual.testerror = data.frame(p, pse.vec, psp.vec, 
                                         pppv.vec, pnpv.vec, row.names = NULL),
       group.testerror = c(PSe = pse, PSp = psp, PPPV = pppv, PNPV = pnpv),
       individual.probabilities = p)
  
}
#end hierarchical.desc() functions

#end hierarchical and additional functions



#Start get.CRC, uses same additional functions as above for hierarchical.desc()
###############################################################################
#    Michael Black 03 12 2013
#     Function: get.CRC gets CRC or ORC the optimal retesting configuration 
#               for 2, 3 or 4 stages
#       calls:
#       inputs: p vector of probabilities,
#               se and sp sensitivity and specificity,,
#               stages is number of stages to optimize for can be 2, 3 or 4
#               order.p = TRUE if False p vector not sorted
#
#        note:
#
###############################################################################

# get.CRC <- function(p, se = 1, sp = 1, stages = 2, order.p = TRUE,
#                     everycase = FALSE, init.config = "hom")   {
#   if (order.p == TRUE) p = sort(p)
# 
#   if (stages == 2) {
#     print("Two stage immediately retests individually if initial group is positive")
#     CRC.desc <- hierarchical.desc2(p = p, se = se, sp = sp, I2 = NULL, I3 = NULL, 
#                                   order.p = order.p)
#   }
# 
#   if (everycase != TRUE) {
#     if (stages == 3) {
#       #in "integer_programming.R"
#       CRC <- Opt.grps.size_number.speed.3step(p = p, se = se, sp = sp) 
#       CRC.desc <- hierarchical.desc2(p=p, se = se, sp = sp, 
#                                     I2 = CRC[1:length(CRC) - 1], I3 = NULL, 
#                                     order.p = order.p)
#     }
#     if (stages == 4) {
#       #in "ET4stepOptimizing_steeper_4s_only.R"
#       CRC <- Opt.grps.size_number.4step(p = p, se = se, sp =sp) 
#       change.notation <- get.mc.4s(vec.g2 = CRC$vec.g2, vec.g3 = CRC$vec.g3)
#       CRC.desc <- hierarchical.desc2(p = p, se = se, sp = sp,
#                                     I2 = change.notation$vec.I2, 
#                                     I3 = change.notation$vec.I3,
#                                     order.p = order.p)
#     }
#   }
# 
#   if (everycase == TRUE) {
#     if (stages == 3) {
# 
#       warning("If group size is large (>18) program may take excessive time")
# 
#       #in "integer_programming.R"
#       CRC <- Opt.grps.size_number_speedg.3step(p = p, sp = sp, se = se)
#       CRC.desc <- hierarchical.desc2(p = p, se = se, sp = sp, 
#                                     I2 = CRC[1:length(CRC) - 1], I3 = NULL, 
#                                     order.p = order.p)
#     }
#     if (stages == 4) {
# 
#       warning("If group size is large (>13) program may take excessive time")
# 
#       #in "ET4stepOptimizing_steeper_4s_only.R"
#       CRC <- Opt.grps.size_number.4step.slow(p = p, sp = sp, se = se) 
#       change.notation <- get.mc.4s(vec.g2 = CRC$vec.g2, vec.g3 = CRC$vec.g3)
#       CRC.desc <- hierarchical.desc2(p = p, se = se, sp = sp,
#                                     I2 = change.notation$vec.I2, 
#                                     I3 = change.notation$vec.I3, order.p = order.p)
#     }
#     print("ORC is the optimal looking at every possible configuration")
#   }
#   CRC.desc
# }


