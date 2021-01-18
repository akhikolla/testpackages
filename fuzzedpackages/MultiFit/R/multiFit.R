#' @title Multiscale Multivariate Framework for Independence Tests
#' @description Perform multiscale test of independence for multivariate vectors.
#' See vignettes for further examples.
#' @param xy A list, whose first element corresponds to the matrix x as below, and
#' its second element corresponds to the matrix y as below. If \code{xy} is not
#' specified, \code{x} and \code{y} need to be assigned.
#' @param x A matrix, number of columns = dimension of random vector,
#' number of rows = number of observations.
#' @param y A matrix, number of columns = dimension of random vector,
#' number of rows = number of observations.
#' @param p_star Numeric, cuboids associated with tests whose p-value is below \code{p_star}
#' will be halved and further tested. 
#' @param R_star A positive integer, if set to an integer
#' between 0 and \code{R_max}, all tests up to and including resolution \code{R_star}
#' will be performed (algorithm will stop at a lower resolution than requested if
#' all tables in it do not meet the criteria specified at \code{min.tbl.tot},
#' \code{min.row.tot} and \code{min.col.tot}). For higher resolutions only the children of 
#' tests with p-value lower than \code{p_star} will be considered.
#' @param R_max A positive integer (or Inf), the maximal number of
#' resolutions to scan (algorithm will stop at a lower resolution if
#' all tables in it do not meet the criteria specified at \code{min.tbl.tot},
#' \code{min.row.tot} and \code{min.col.tot})
#' @param rank.transform Logical, if \code{TRUE}, marginal rank transform is
#' performed on all margins of \code{x} and \code{y}. If \code{FALSE}, all
#' margins are scaled to 0-1 scale. When \code{FALSE}, the average and top
#' statistics of the negative logarithm of the \code{p}-values are only computed
#' for the univariate case.
#' @param test.method String, choose "Fisher" for Fisher's exact test (slowest), "chi.sq" for
#' Chi-squared test, "LR" for likelihood-ratio test and "norm.approx" for approximating
#' the hypergeometric distribution with a normal distribution (fastest).
#' @param correct Logical, if \code{TRUE} compute mid-p corrected \code{p}-values for
#' Fisher's exact test, or Yates corrected \code{p}-values for the Chi-squared test,
#' or Williams corrected \code{p}-values for the likelihood-ratio test.
#' @param min.tbl.tot Non-negative integer, the minimal number of observations
#' per table below which a \code{p}-value for a given table will not be computed.
#' @param min.row.tot Non-negative integer, the minimal number of observations
#' for row totals in the 2x2 contingency tables below which a contingency table will not be tested.
#' @param min.col.tot Non-negative integer, the minimal number of observations
#' for column totals in the 2x2 contingency tables below which a contingency table will not be tested.
#' @param p.adjust.methods String, choose between "H" for Holm, "Hcorrected" for Holm with
#' the correction as specified in \code{correct}, or "MH" for Modified-Holm (for Fisher's
#' exact test only). See documentation for further details.
#' @param compute.all.holm Logical, if \code{FALSE}, only global \code{p}-value is
#' computed (may be faster, especially when Modified-Holm correction is used). If \code{TRUE}
#' adjusted \code{p}-values are computed for all tests.
#' @param cutoff Numerical between 0 and 1, an upper limit for the \code{p}-values
#' that are to be adjusted (the lower the cutoff - the fewer computations are required
#' for the Modified Holm method).
#' @param top.max.ps Positive integer, report the mean of the top \code{top.max.ps}
#' order statistics of the negative logarithm of all \code{p}-values.
#' @param return.all.pvs Logical, if TRUE, a data frame with all \code{p}-values
#' is returned.
#' @param save.all.pvs Logical, if \code{TRUE} a data frame with all \code{p}-values
#' is saved to a RData file named according to \code{all.pvs.fname}
#' @param all.pvs.fname String, file name to which all \code{p}-values are saved if
#' \code{save.all.pvs} is \code{TRUE}.
#' @param uv.approx.null Logical, in a univariate case, if \code{TRUE} and the
#' testing method is either Fisher's exact test or the normal approximation of
#' the hypergeometric distribution, an approximate null distribution for the global
#' test statistics is simulated. See documentation for further details.
#' @param uv.exact.null Logical, in a univariate case, if \code{TRUE} and the
#' testing method is either Fisher's exact test or the normal approximation of
#' the hypergeometric distribution, an exact null distribution for the global
#' test statistics is simulated. See documentation for further details.
#' @param uv.null.sim Positive integer, the number of simulated values to be computed
#' in a univariate case when an exact or approximate null distribution is simulated.
#' @param plot.marginals Logical, if \code{TRUE} plots the marginal scatter plots between
#' all pairs of margins of x and y, before and after rank transforming or scaling them.
#' @param rk Logical, if \code{FALSE}, select only tests with \code{p}-values 
#' more extreme than \code{p_star} to halve and further test. FWER control guaranteed. 
#' If \code{TRUE}, choose at each resolution the \code{M} tests with the most extreme 
#' \code{p}-values to further halve and test.
#' @param M A positive integer (or Inf), the number of top ranking
#' tests to continue to split at each resolution. FWER control not guaranteed 
#' for this method.
#' @param verbose Logical.
#' @return \code{test.stats}, a named numerical vector containing the test
#' statistics for the global null hypothesis (i.e. x independent of y)
#' @return \code{p.values}, a named numerical vector containing the \code{p}-values of
#' for the global null hypothesis (i.e. x independent of y). These are not computed
#' if \code{p.adjust.methods} is \code{NULL}.
#' @return \code{pvs}, a data frame that contains all \code{p}-values and adjusted
#' \code{p}-values that are computed. Returned if \code{return.all.pvs} is \code{TRUE}.
#' @return \code{all}, a nested list. Each entry is named and contains data about a resolution
#' that was tested. Each resolution is a list in itself, with \code{cuboids}, a summary of
#' all tested cuboids in a resolution, \code{tables}, a summary of all 2x2
#' contingency tables in a resolution, \code{pv}, a numerical vector containing the
#' \code{p}-values from the tests of independence on 2x2 contingency table in \code{tables}
#' that meet the criteria defined by \code{min.tbl.tot}, \code{min.row.tot} and \code{min.col.tot}.
#' The length of \code{pv} is equal to the number of rows of \code{tables}. \code{pv.correct},
#' similar to the above \code{pv}, corrected \code{p}-values are computed and returned when
#' \code{correct} is \code{TRUE}. \code{rank.tests}, logical vector that indicates
#' whether or not a test was ranked among the top \code{M} tests in a resolution. The
#' length of \code{rank.tests} is equal to the number of rows of \code{tables}. \code{parent.cuboids},
#' an integer vector, indicating which cuboids in a resolution are associated with
#' the ranked tests, and will be further halved in the next higher resolution.
#' \code{parent.tests}, a logical vector of the same length as the
#' number of rows of \code{tables}, indicating whether or not a test was chosen as a parent
#' test (same tests may have multiple children).
#' @return \code{approx.nulls}, in a univariate case, a list of numerical vectors
#' whose values are the simulated approximate null values.
#' @return \code{exact.nulls}, in a univariate case, a list of numerical vectors
#' whose values are the simulated theoretical null values.
#' @examples
#' set.seed(1)
#' n = 300
#' Dx = Dy = 2
#' x = matrix(0, nrow=n, ncol=Dx)
#' y = matrix(0, nrow=n, ncol=Dy)
#' x[,1] = rnorm(n)
#' x[,2] = runif(n)
#' y[,1] = rnorm(n)
#' y[,2] = sin(5*pi*x[,2]) + 1/5*rnorm(n)
#' fit = multiFit(x=x, y=y, verbose=TRUE)
#' w = multiSummary(x=x, y=y, fit=fit, alpha=0.0001)
#' @export
multiFit = function(xy, x = NULL, y = NULL, 
                    p_star = NULL, R_max = NULL,
                    R_star = 1,
                    rank.transform = TRUE,
                    test.method = "Fisher", correct=TRUE,
                    min.tbl.tot = 25L, min.row.tot = 10L, min.col.tot = 10L,
                    p.adjust.methods = c("H", "Hcorrected", "MH"),
                    compute.all.holm  = TRUE,
                    cutoff = 0.05, top.max.ps = 4L,
                    return.all.pvs = TRUE, save.all.pvs = FALSE, all.pvs.fname = NULL,
                    uv.approx.null = FALSE, uv.exact.null = FALSE,
                    uv.null.sim = 10000L, plot.marginals = FALSE,
                    rk = FALSE, M = 10, verbose = FALSE) {
  
# :::::::::::::::::::::::::::::::::
  
  s0=Sys.time()
  
  # Pre-processing:
  if (!missing(xy)) {
    x=xy[[1]]
    y=xy[[2]]
  }
  
  if (is.vector(x)) x=matrix(x,ncol=1)
  if (is.vector(y)) y=matrix(y,ncol=1)
  
  Dx = ncol(x); Dy = ncol(y); n = nrow(x)
  if (is.null(R_max)) R_max = as.integer(floor(log(n/10,2)))
  if (!is.null(R_star)) {
    if (R_star>=R_max) {
      warning(paste0("Full coverage was asked for resolution ",
                     R_star,
                     " whereas the maximal resolution is ",
                     R_max,
                     ". All cuboids are scanned."))
      M=Inf
    }
  } else { R_star = -1 }
  
  if (is.null(p_star)) p_star = 1/(Dx*Dy*log(n,2))
  if (is.null(R_star)) R_star = 1
  if (!rk) M=Inf
  
  if (n<50 & min.tbl.tot==25L & min.row.tot==10L & min.col.tot==10L) {
    min.tbl.tot = as.integer(floor(n/4))
    min.row.tot = as.integer(floor(0.4*min.tbl.tot))
    min.col.tot = as.integer(floor(0.4*min.tbl.tot))
  }
  
  if (!is.null(p.adjust.methods)) {
    if (sum(!is.na(p.adjust.methods))==0) {
      p.adjust.methods=NULL
    } else {
      switch (test.method,
              Fisher = {
                test.name = "Fisher's exact tests"
                choices = c("H", "Hcorrected", "MH")
                if ("Hcorrected"%in%p.adjust.methods | "Hcorrected"%in%p.adjust.methods) correct=TRUE
              },
              chi.sq = {
                test.name ="Chi-squared tests"
                p.adjust.methods = setdiff(p.adjust.methods, c("MH"))
                choices = c("H", "Hcorrected")
                if ("Hcorrected"%in%p.adjust.methods | "Hcorrected"%in%p.adjust.methods) correct=TRUE
              },
              LR = {
                test.name ="likelihood-ratio tests"
                p.adjust.methods = setdiff(p.adjust.methods, c("MH"))
                choices = c("H", "Hcorrected")
                if ("Hcorrected"%in%p.adjust.methods | "Hcorrected"%in%p.adjust.methods) correct=TRUE
              },
              norm.approx = {
                test.name = "normal approximation based tests"
                correct=FALSE
                p.adjust.methods = setdiff(p.adjust.methods, c("Hcorrected", "MH"))
                choices = c("H")
              })
      
      if (any(!(p.adjust.methods%in%choices))) {
        stop(paste0("p.adjust.methods not specified correctly. For ",test.name,
                    " p.adjust.methods must be NULL/NA or include ",
                    paste(paste0(choices[1:(length(choices)-1)],", "),collapse=""),
                    "or ",choices[length(choices)],".",sep=""))
      }}
  }
  
  names.test.stats=NULL
  # if (rank.transform | (Dx==1 & Dy==1)) names.test.stats = c("avg.log.pv", paste0("avg.log.top",top.max.ps,".pv"))
  # if ((rank.transform | (Dx==1 & Dy==1)) & correct) {
  if (!correct) names.test.stats = c("avg.log.pv", paste0("avg.log.top",top.max.ps,".pv"))
  if (correct) {
    names.test.stats=c(names.test.stats, "avg.log.pv.corrected",
                       paste0("avg.log.top",top.max.ps,".pv.corrected"))
  }
  if (length(names.test.stats)>0) {
    test.stats = rep(NA, length(names.test.stats))
    names(test.stats) = names.test.stats
  }
  
  names.p.values = NULL
  if (uv.exact.null & Dx==1 & Dy==1 & (test.method!="Fisher" | correct)) names.p.values = c("pv.avg.exact", paste0("pv.top",top.max.ps,".exact"))
  if (uv.approx.null & Dx==1 & Dy==1 & (test.method!="Fisher" | correct)) names.p.values = c(names.p.values, "pv.avg.approx", paste0("pv.top",top.max.ps,".approx"))
  if (!is.null(p.adjust.methods)) names.p.values = c(names.p.values, p.adjust.methods)
  if (length(names.p.values)>0) {
    p.values = rep(NA, length(names.p.values))
    names(p.values) = names.p.values }
  
  #if (!("MH"%in%p.adjust.methods)) compute.all.holm = FALSE
  
  if (plot.marginals) {
    for (i in 1:Dx) for (j in 1:Dy) plot(x[,i],y[,j],xlab=paste0("x",i),ylab=paste0("y",j))
  }
  
  if (rank.transform) {
    if (verbose) cat("Applying rank transformation\n")
    x.tran=apply(x,2,function(x) ecdf(x)(x)-1/n)
    y.tran=apply(y,2,function(y) ecdf(y)(y)-1/n)
  }

  if (!rank.transform) {
    scale01 = function(z){(1-10^(-7))*(z-min(z))/(max(z)-min(z))}
    if (verbose) cat("Scaling variables\n")
    x.tran = apply(x,2,scale01)
    y.tran = apply(y,2,scale01)
  }
  
  if (plot.marginals) {
    for (i in 1:Dx) for (j in 1:Dy) plot(x.tran[,i],y.tran[,j],xlab=paste0("x.tran",i),ylab=paste0("y.tran",j))
  }
  
  #-------------------------------------------------------------------------------------------
  
  # Get to work:
  kx.header = paste0("kx.",1:Dx)
  ky.header = paste0("ky.",1:Dy)
  k.header = c(kx.header, ky.header)
  lx.header = paste0("lx.",1:Dx)
  ly.header = paste0("ly.",1:Dy)
  l.header = c(lx.header, ly.header)
  kl.header = c(k.header, l.header)
  cuboid.header = c("n00", "n01", "n10", "n11")
  ij = expand.grid(1L:Dy,1L:Dx)
  ij = as.matrix(ij[,c(2,1)])
  
  cuboids.colnames = c("res", kl.header, "cuboid.idx", "child.of.cuboid", "child.of.test")
  W=as.data.frame(matrix(0L, nrow=1, ncol=length(cuboids.colnames)))
  
  colnames(W) = cuboids.colnames
  
  # Initial values for resolution 0:
  W[,"cuboid.idx"] = 1L
  W[,l.header] = 1L
  
  rng=1L:(Dx*Dy) # rng: range of tests(=rows in pvs) defined for each resolution
  
  if (verbose) cat("Testing")
  if (verbose & correct & test.method=="Fisher") cat(" and computing mid-p corrected p-values")
  if (verbose & correct & test.method=="chi.sq") cat(" and computing Yates corrected p-values")
  if (verbose & correct & test.method=="LR") cat(" and computing Williams corrected p-values")
  if (verbose) cat(paste0(":\nResolution ",0,"/",R_max, ": performing ",length(rng)," tests\n"))
  
  discCpp = discretizeCpp(a=x.tran, b=y.tran,
                          w=as.numeric(W[,c(kl.header,"cuboid.idx")]),
                          mask=rep(TRUE, n),
                          ij=ij, Dx=Dx, Dy=Dy);
  colnames(discCpp$tables) =  c("i", "j",
                                "n00", "n01", "n10", "n11",
                                "cuboid.idx", "test.idx");
  disc = list();
  disc$tables = discCpp$tables;
  
  disc$`2`$mask = as.logical(discCpp$mask);
  disc$`2`$x0.sub.mask = matrix(as.logical(discCpp$x0.sub.mask), nrow = nrow(discCpp$x0.sub.mask));
  disc$`2`$y0.sub.mask = matrix(as.logical(discCpp$y0.sub.mask), nrow = nrow(discCpp$y0.sub.mask));
  names(disc)[2] = "";
  
  
  # **********************************
  disc$tables[,"test.idx"] = rng
  if (return.all.pvs) tx = disc$tables[,"test.idx"]
  
  t = matrix(integer(), nrow=nrow(disc$tables), ncol=4+5+2*(test.method=="Fisher"))
  if (test.method=="Fisher") {
    colnames(t) = c(cuboid.header, "row0.tot", "row1.tot", "col0.tot", "col1.tot", "grand.tot", "hi", "lo")
  } else {
    colnames(t) = c(cuboid.header, "row0.tot", "row1.tot", "col0.tot", "col1.tot", "grand.tot")
  }
  
  if (is.matrix(disc$tables[,cuboid.header])) {
    t[,cuboid.header] = disc$tables[,cuboid.header]
    t[,"col0.tot"] = rowSums(t[,c("n00","n10")])
    t[,"col1.tot"] = rowSums(t[,c("n01","n11")])
    t[,"grand.tot"] = rowSums(t[,c("col0.tot","col1.tot")])
    t[,"row1.tot"] = rowSums(t[,c("n10","n11")])
    t[,"row0.tot"] = rowSums(t[,c("n00","n01")])
  } else {
    t[,cuboid.header] = matrix(disc$tables[,cuboid.header], nrow=1)
    t[,"col0.tot"] = sum(t[,c("n00","n10")])
    t[,"col1.tot"] = sum(t[,c("n01","n11")])
    t[,"grand.tot"] = sum(t[,c("col0.tot","col1.tot")])
    t[,"row1.tot"] = sum(t[,c("n10","n11")])
    t[,"row0.tot"] = sum(t[,c("n00","n01")])
  }
  
  do.tests = ((t[,"grand.tot"]>min.tbl.tot) & (t[,"row0.tot"]>min.row.tot)
              & (t[,"row1.tot"]>min.row.tot) & (t[,"col0.tot"]>min.col.tot)
              & (t[,"col1.tot"]>min.col.tot))
  
  if (sum(do.tests)==0) {
    if (verbose) cat("No pairs of margins in resolution 0 had enough observations to be tested,
                     no tests performed.\n")
    return(NULL)
  }
  
  if (test.method=="Fisher") {
    t[,"hi"] = t[,"n00"]+pmin(t[,"n10"],t[,"n01"])
    t[,"lo"] = pmax(0, t[,"col0.tot"]-t[,"row1.tot"])
  }
  
  if (test.method=="chi.sq" | test.method=="LR") {
    if (is.matrix(disc$tables[,cuboid.header]) & sum(do.tests)>1) {
      E = exp(c(rowSums(log(t[do.tests,c("row0.tot","col0.tot")])),
                rowSums(log(t[do.tests,c("row0.tot","col1.tot")])),
                rowSums(log(t[do.tests,c("row1.tot","col0.tot")])),
                rowSums(log(t[do.tests,c("row1.tot","col1.tot")])))-log(t[do.tests,"grand.tot"]))
    } else {
      E = exp(c(sum(log(t[do.tests,c("row0.tot","col0.tot")])),
                sum(log(t[do.tests,c("row0.tot","col1.tot")])),
                sum(log(t[do.tests,c("row1.tot","col0.tot")])),
                sum(log(t[do.tests,c("row1.tot","col1.tot")])))-log(t[do.tests,"grand.tot"]))
    }
  }
  
  all = list(list(cuboids=W[,cuboids.colnames], tables=disc$tables))
  names(all) = '0'
  
  switch (test.method,
          Fisher = {
            ret.all.probs = ("MH"%in%p.adjust.methods)
            if (sum(do.tests)==1) {
              ST = single_Fisher_test(t[do.tests,],  correct, ret.all.probs)
              if (ret.all.probs) {
                ALL.PROBS = list(ST$all.probs)
              }
              pv = rep(NA, nrow(t))
              pv[do.tests] = ST$pv
              if (correct) {
                pv.correct = rep(NA, nrow(t))
                pv.correct[do.tests] = ST$pv.correct
              }
            } else {
              ST = apply(t[do.tests,], 1, single_Fisher_test, correct, ret.all.probs)
              if (ret.all.probs) {
                ALL.PROBS = lapply(ST, `[[`, "all.probs")
              }
              pv = rep(NA, nrow(t))
              pv[do.tests] = vapply(ST, `[[`, double(1), "pv")
              if (correct) {
                pv.correct = rep(NA, nrow(t))
                pv.correct[do.tests] = vapply(ST, `[[`, double(1), "pv.correct")
              }
            }
          },
          chi.sq = {
            U = abs(c(t[do.tests, cuboid.header]) - E)
            stat = rowSums(matrix(U^2/E,ncol=4))
            if (correct) {stat.yates = rowSums(matrix((U - pmin(0.5, U))^2/E,ncol=4))} else {stat.yates = NULL}
            p = pchisq(c(stat, stat.yates), 1, lower.tail = FALSE)
            pv = rep(NA, nrow(t))
            pv[do.tests] = p[1:((length(p)/2)+as.integer(!correct)*(length(p)/2))]
            pv[pv==0] = .Machine$double.xmin
            if (correct) {
              pv.correct = rep(NA, nrow(t))
              pv.correct[do.tests] = p[(length(p)/2+1):length(p)]
              pv.correct[pv.correct==0] = .Machine$double.xmin
            }
          },
          LR = {
            t[do.tests,cuboid.header][t[do.tests,cuboid.header]==0]=NA
            U = c(t[do.tests,cuboid.header])
            stat = 2*rowSums(matrix(U*log(U/E),ncol=4), na.rm=TRUE)
            if (correct) {
              row.tot = 1/t[do.tests,"row0.tot"]+1/t[do.tests,"row1.tot"]
              col.tot = 1/t[do.tests,"col0.tot"]+1/t[do.tests,"col1.tot"]
              gt = t[do.tests,"grand.tot"]
              q = 1 + exp(log(gt*row.tot-1)+log(gt*col.tot-1)-log(6*gt))
              stat.williams = stat/q
            } else {
              stat.williams = NULL
            }
            p = pchisq(c(stat,stat.williams), 1, lower.tail = FALSE)
            pv = rep(NA, nrow(t))
            pv[do.tests] = p[1:((length(p)/2)+as.integer(!correct)*(length(p)/2))]
            pv[pv==0] = .Machine$double.xmin
            if (correct) {
              pv.correct = rep(NA, nrow(t))
              pv.correct[do.tests] = p[(length(p)/2+1):length(p)]
              pv.correct[pv.correct==0] = .Machine$double.xmin
            }
          },
          norm.approx = {
            c0 = t[do.tests,"col0.tot"]
            r0 = t[do.tests,"row0.tot"]
            gt = t[do.tests,"grand.tot"]
            EA00 = exp(log(c0)+log(r0)-log(gt))
            VarA00 = exp(log(EA00) + log(gt-r0)-log(gt) +log(gt-c0)-log(gt-1))
            pv = rep(NA, nrow(t))
            pv[do.tests] = 2*pnorm(-abs((t[do.tests,"n00"]-EA00)/sqrt(VarA00)))
            pv[pv==0] = .Machine$double.xmin
          }
  )
  
  if (correct) {
    all$`0` = c(all$`0`, list(pv=pv, pv.correct=pv.correct))
  } else {
    all$`0` = c(all$`0`, list(pv=pv))
  }
  
  disc$tables=NULL
  selector = disc
  names(selector) = "1"
  p.rng = !is.na(all$`0`$pv)
  
  if (sum(is.na(all$`0`$pv))==length(all$`0`$pv)) {
    warning("No tests were performed, too few observations per table. Consider changing min.tbl.tot, min.row.tot or min.col.tot.")
    if (return.all.pvs) { return(invisible(all)) } else { return(NULL) }
  }
  
  # Rank the tests that have p-values:
  pvs.per.res = all$`0`$pv
  n.pvs.per.res = sum(p.rng)
  if (rk) {
    mM = min(n.pvs.per.res, M)
    cond = (mM>=n.pvs.per.res | 0<R_star)
  } else {
    mM = n.pvs.per.res
    cond = (0<R_star)
  }
  
  if (cond) {
    rank.tests = rep(TRUE, n.pvs.per.res)
  } else {
    if (!rk) {
      rank.tests =
        as.logical(rep(FALSE, n.pvs.per.res) +
                     (pvs.per.res<=p_star))
    } else {
      rank.tests =
        as.logical(rep(FALSE, n.pvs.per.res) +
                     (pvs.per.res<=sort(pvs.per.res,partial=1:mM)[mM]))
    }
  }

  parent.tests = rank.tests
  parent.tests[is.na(parent.tests)] = 0
  parent.tests = as.logical(parent.tests)
  num.parent.tests = sum(parent.tests)
  
  # On the annoying occasion where several different tests are ranked last with the
  # same p-value (not negligible chance, because of the rank transform),
  # omit until down to M ranked tests:
  if (0>=R_star) {
    while (num.parent.tests>M) {
      rank.tests[rank.tests][tail(which(pvs.per.res[rank.tests]==max(pvs.per.res[rank.tests],na.rm=TRUE)),1)]=FALSE
      parent.tests = rank.tests
      num.parent.tests = sum(parent.tests,na.rm=TRUE)
    }
  }
  parent.cuboids = unname(all$`0`$tables[parent.tests,"cuboid.idx"])
  num.parent.cuboids = length(unique(parent.cuboids))
  
  all$`0` = c(all$`0`, list(rank.tests=rank.tests, parent.cuboids=parent.cuboids, parent.tests=parent.tests))
  unq = !duplicated(W$cuboid.idx)
  
  if (sum(rank.tests)==0) {
    if (verbose) cat("No potential parent tests in resolution 0 have p-values below threshold.\n")
    R_max=0
  }
  
  if (R_max>0) {
    # Resolutions 1 and above:
    res=1
    old.kx.header = paste0("parent.",kx.header)
    old.ky.header = paste0("parent.",ky.header)
    old.k.header = c(old.kx.header, old.ky.header)
    lx.add.header = paste0("add.",lx.header)
    ly.add.header = paste0("add.",ly.header)
    cuboids.colnames.extended = c("res", kl.header, "cuboid.idx", "child.of.cuboid", "child.of.test",
                               old.k.header, lx.add.header, ly.add.header)
    while (res <= R_max) {
      # Find unique tests:
      r.prev=as.character(res-1)
      curr.res = as.character(res)
      min.cuboid.idx = max(all[[r.prev]]$cuboid[,"cuboid.idx"])+1L
      min.test.idx = max(all[[r.prev]]$tables[,"test.idx"])+1L
      mM = ifelse(res<=R_star,
                  num.parent.tests,
                  min(num.parent.tests,M))
      all.max.cuboid.idx = min.cuboid.idx-1L+4L*mM
      all.cuboid.rng = min.cuboid.idx:all.max.cuboid.idx
      if (verbose) cat(paste0("Resolution ",res,"/",R_max, ": scanning ",length(all.cuboid.rng)," cuboids"))
      
      W=as.data.frame(matrix(0L, nrow=length(all.cuboid.rng), ncol=length(cuboids.colnames.extended)+1))
      colnames(W) = c(cuboids.colnames.extended, "id")
      
      W[,"res"] = res
      
      parent.cuboids.row.num = match(parent.cuboids, all[[r.prev]]$cuboids$cuboid.idx)
      old.rep.idx = rep(parent.cuboids.row.num, each=4)
      
      W[,"child.of.cuboid"] = all[[r.prev]]$cuboids[old.rep.idx,"cuboid.idx"]
      
      # This is gross, but all options differ dependeing on whether this is
      # a 1x1 case, there is only one parent selected, or there is only one
      # p-value in the previous resolution:
      if (Dx==1 & Dy==1 & nrow(all[[r.prev]]$cuboids)==1) {
        all.parent.tests.tables =
          matrix(all[[r.prev]]$tables[p.rng,][all[[r.prev]]$parent.tests],nrow=1,byrow=TRUE)
        colnames(all.parent.tests.tables) = c("i", "j", "n00", "n01", "n10", "n11", "cuboid.idx", "test.idx")
      } else {
        if (sum(p.rng)>1 & num.parent.tests>1) {
          all.parent.tests.tables = all[[r.prev]]$tables[p.rng,][all[[r.prev]]$parent.tests,]
        } else {
          if (sum(p.rng)==1) {
            all.parent.tests.tables = matrix(all[[r.prev]]$tables[p.rng,],nrow=1)
          } else {
            all.parent.tests.tables = matrix(all[[r.prev]]$tables[p.rng,][all[[r.prev]]$parent.tests,],nrow=1)
          }
          colnames(all.parent.tests.tables) = c("i", "j", "n00", "n01", "n10", "n11", "cuboid.idx", "test.idx")
        }
      }
      
      W[,"child.of.test"] = rep(as.integer(all.parent.tests.tables[,"test.idx"]), each=4)
      
      kx.add = do.call("rbind", lapply(all.parent.tests.tables[,"i"],
                                       function(i) {
                                         a=matrix(0L,nrow=2,Dx); a[1,i]=1L;
                                         if (is.vector(a[rep(1L:2L, each=2),])) {
                                           return(matrix(a[rep(1L:2L, each=2),],ncol=1))
                                         } else {
                                           return(a[rep(1L:2L, each=2),])
                                         }}))
      ky.add = do.call("rbind", lapply(all.parent.tests.tables[,"j"],
                                       function(j) {
                                         a=matrix(0L,nrow=2,ncol=Dy); a[2,j]=1L;
                                         if (is.vector(a[rep(1L:2L, each=2),])) {
                                           return(matrix(a[rep(1L:2L, each=2),],ncol=1))
                                         } else {
                                           return(a[rep(1L:2L, each=2),])
                                         }}))
      W[,old.k.header] = all[[r.prev]]$cuboids[old.rep.idx,k.header]
      W[,kx.header]=W[,old.kx.header]+kx.add
      W[,ky.header]=W[,old.ky.header]+ky.add
      
      W[,lx.add.header] = do.call("rbind", lapply(all.parent.tests.tables[,"i"],
                                                  function(i) {
                                                    a=matrix(1L,nrow=4,Dx); a[,i]=c(1L,2L,1L,1L);
                                                    return(a)
                                                  }))
      W[,ly.add.header] = do.call("rbind", lapply(all.parent.tests.tables[,"j"],
                                                  function(j) {
                                                    a=matrix(1L,nrow=4,ncol=Dy); a[,j]=c(1L,1L,1L,2L);
                                                    return(a)
                                                  }))
      old.l = all[[r.prev]]$cuboids[old.rep.idx,l.header]
      W[,lx.header] = kx.add*(2L*(old.l[,lx.header]-1L)+W[,lx.add.header]) + (1L-kx.add)*old.l[,lx.header]
      W[,ly.header] = ky.add*(2L*(old.l[,ly.header]-1L)+W[,ly.add.header]) + (1L-ky.add)*old.l[,ly.header]
      
      W[,"cuboid.idx"] = all.cuboid.rng
      
      id = NULL;
      invisible(data.table::setDT(W)[, id := .GRP, by = kl.header]$id)
      W$cuboid.idx = min.cuboid.idx - 1L + W$id
      W=data.frame(W)
      W$id=NULL
      
      unq = !duplicated(W$cuboid.idx)
      rng = all.cuboid.rng[unq]
      if (verbose) cat(paste0(", performing ",length(rng)*Dx*Dy," tests\n"))
      
      disc = apply(W[unq,],1, function(G) {
        which.old.cuboid = as.character(G["child.of.cuboid"])
        pm = selector[[which.old.cuboid]]$mask
        if (which(as.logical(G[k.header]-G[old.k.header]))<=Dx) {
          icol = all.parent.tests.tables[all.parent.tests.tables[,"test.idx"]==G["child.of.test"],]["i"]
          if (G[lx.add.header[icol]]==1) {
            pm[pm==TRUE] = selector[[which.old.cuboid]]$x0.sub.mask[,icol]
          } else {
            pm[pm==TRUE] = !selector[[which.old.cuboid]]$x0.sub.mask[,icol]
          }
        } else {
          jcol = all.parent.tests.tables[all.parent.tests.tables[,"test.idx"]==G["child.of.test"],]["j"]
          if (G[ly.add.header[jcol]]==1) {
            pm[pm==TRUE] = selector[[which.old.cuboid]]$y0.sub.mask[,jcol]
          } else {
            pm[pm==TRUE] = !selector[[which.old.cuboid]]$y0.sub.mask[,jcol]
          }
        }
        {
          discCpp = discretizeCpp(a=x.tran, b=y.tran,
                                  w=as.numeric(G[c(kl.header,"cuboid.idx")]),
                                  mask=pm,
                                  ij=ij, Dx=Dx, Dy=Dy);
          colnames(discCpp$tables) =  c("i", "j",
                                        "n00", "n01", "n10", "n11",
                                        "cuboid.idx", "test.idx");
          disc = list();
          disc$tables = discCpp$tables;
          
          disc$`2`$mask = as.logical(discCpp$mask);
          disc$`2`$x0.sub.mask = matrix(as.logical(discCpp$x0.sub.mask), nrow = nrow(discCpp$x0.sub.mask));
          disc$`2`$y0.sub.mask = matrix(as.logical(discCpp$y0.sub.mask), nrow = nrow(discCpp$y0.sub.mask));
          names(disc)[2] = "";
          disc  
        }
      })
      
      dt = do.call("rbind", lapply(disc, `[[`,"tables"))
      dt[,"test.idx"] = min.test.idx:(min.test.idx+nrow(dt)-1)
      if (return.all.pvs) tx = c(tx, dt[,"test.idx"])
      
      all.add = list(list(cuboids=W[,cuboids.colnames], tables=dt))
      names(all.add) = res
      all=c(all, all.add)
      
      t = matrix(integer(), nrow=nrow(dt), ncol=4+5+2*(test.method=="Fisher"))
      if (test.method=="Fisher") {
        colnames(t) = c(cuboid.header, "row0.tot", "row1.tot", "col0.tot", "col1.tot", "grand.tot", "hi", "lo")
      } else {
        colnames(t) = c(cuboid.header, "row0.tot", "row1.tot", "col0.tot", "col1.tot", "grand.tot")
      }
      
      if (is.matrix(dt[,cuboid.header])) {
        t[,cuboid.header] = dt[,cuboid.header]
        t[,"col0.tot"] = rowSums(t[,c("n00","n10")])
        t[,"col1.tot"] = rowSums(t[,c("n01","n11")])
        t[,"grand.tot"] = rowSums(t[,c("col0.tot","col1.tot")])
        t[,"row1.tot"] = rowSums(t[,c("n10","n11")])
        t[,"row0.tot"] = rowSums(t[,c("n00","n01")])
      } else {
        t[,cuboid.header] = matrix(dt[,cuboid.header], nrow=1)
        t[,"col0.tot"] = sum(t[,c("n00","n10")])
        t[,"col1.tot"] = sum(t[,c("n01","n11")])
        t[,"grand.tot"] = sum(t[,c("col0.tot","col1.tot")])
        t[,"row1.tot"] = sum(t[,c("n10","n11")])
        t[,"row0.tot"] = sum(t[,c("n00","n01")])
      }
      
      do.tests = ((t[,"grand.tot"]>min.tbl.tot) & (t[,"row0.tot"]>min.row.tot)
                  & (t[,"row1.tot"]>min.row.tot) & (t[,"col0.tot"]>min.col.tot)
                  & (t[,"col1.tot"]>min.col.tot))
      
      if (sum(do.tests)==0) {
        if (correct) {
          all[[curr.res]] = c(all[[curr.res]], list(pv=rep(NA,nrow(t)), pv.correct=rep(NA,nrow(t))))
        } else {
          all[[curr.res]] = c(all[[curr.res]], list(pv=rep(NA,nrow(t))))
        }
        if (verbose) cat("No pairs of margins in resolution", res, "had enough observations to be tested.\n")
        break
      }
      
      if (test.method=="Fisher") {
        t[,"hi"] = t[,"n00"]+pmin(t[,"n10"],t[,"n01"])
        t[,"lo"] = pmax(0, t[,"col0.tot"]-t[,"row1.tot"])
      }
      
      if (test.method=="chi.sq" | test.method=="LR") {
        if (is.matrix(dt[,cuboid.header]) & sum(do.tests)>1) {
          E = exp(c(rowSums(log(t[do.tests,c("row0.tot","col0.tot")])),
                    rowSums(log(t[do.tests,c("row0.tot","col1.tot")])),
                    rowSums(log(t[do.tests,c("row1.tot","col0.tot")])),
                    rowSums(log(t[do.tests,c("row1.tot","col1.tot")])))-log(t[do.tests,"grand.tot"]))
        } else {
          E = exp(c(sum(log(t[do.tests,c("row0.tot","col0.tot")])),
                    sum(log(t[do.tests,c("row0.tot","col1.tot")])),
                    sum(log(t[do.tests,c("row1.tot","col0.tot")])),
                    sum(log(t[do.tests,c("row1.tot","col1.tot")])))-log(t[do.tests,"grand.tot"]))
        }
      }
      
      switch (test.method,
              Fisher = {
                ret.all.probs = ("MH"%in%p.adjust.methods)
                if (sum(do.tests)==1) {
                  ST = single_Fisher_test(t[do.tests,],  correct, ret.all.probs)
                  if (ret.all.probs) {
                    ALL.PROBS = c(ALL.PROBS, ST$all.probs)
                  }
                  pv = rep(NA, nrow(t))
                  pv[do.tests] = ST$pv
                  if (correct) {
                    pv.correct = rep(NA, nrow(t))
                    pv.correct[do.tests] = ST$pv.correct
                  }
                } else {
                  ST = apply(t[do.tests,], 1, single_Fisher_test, correct, ret.all.probs)
                  if (ret.all.probs) {
                    ALL.PROBS.add = lapply(ST, `[[`, "all.probs")
                    ALL.PROBS = c(ALL.PROBS, ALL.PROBS.add)
                  }
                  pv = rep(NA, nrow(t))
                  pv[do.tests] = vapply(ST, `[[`, double(1), "pv")
                  if (correct) {
                    pv.correct = rep(NA, nrow(t))
                    pv.correct[do.tests] = vapply(ST, `[[`, double(1), "pv.correct")
                  }
                }
              },
              chi.sq = {
                U = abs(c(t[do.tests, cuboid.header]) - E)
                stat = rowSums(matrix(U^2/E,ncol=4))
                if (correct) {stat.yates = rowSums(matrix((U - pmin(0.5, U))^2/E,ncol=4))} else {stat.yates = NULL}
                p = pchisq(c(stat, stat.yates), 1, lower.tail = FALSE)
                pv = rep(NA, nrow(t))
                pv[do.tests] = p[1:((length(p)/2)+as.integer(!correct)*(length(p)/2))]
                pv[pv==0] = .Machine$double.xmin
                if (correct) {
                  pv.correct = rep(NA, nrow(t))
                  pv.correct[do.tests] = p[(length(p)/2+1):length(p)]
                  pv.correct[pv.correct==0] = .Machine$double.xmin
                }
              },
              LR = {
                t[do.tests,cuboid.header][t[do.tests,cuboid.header]==0]=NA
                U=c(t[do.tests,cuboid.header])
                stat = 2*rowSums(matrix(U*log(U/E),ncol=4), na.rm=TRUE)
                if (correct) {
                  row.tot = 1/t[do.tests,"row0.tot"]+1/t[do.tests,"row1.tot"]
                  col.tot = 1/t[do.tests,"col0.tot"]+1/t[do.tests,"col1.tot"]
                  gt = t[do.tests,"grand.tot"]
                  q = 1 + exp(log(gt*row.tot-1)+log(gt*col.tot-1)-log(6*gt))
                  stat.williams = stat/q
                } else {
                  stat.williams = NULL
                }
                p = pchisq(c(stat,stat.williams), 1, lower.tail = FALSE)
                pv = rep(NA, nrow(t))
                pv[do.tests] = p[1:((length(p)/2)+as.integer(!correct)*(length(p)/2))]
                pv[pv==0] = .Machine$double.xmin
                if (correct) {
                  pv.correct = rep(NA, nrow(t))
                  pv.correct[do.tests] = p[(length(p)/2+1):length(p)]
                  pv.correct[pv.correct==0] = .Machine$double.xmin
                }
              },
              norm.approx = {
                c0 = t[do.tests,"col0.tot"]
                r0 = t[do.tests,"row0.tot"]
                gt = t[do.tests,"grand.tot"]
                EA00 = exp(log(c0)+log(r0)-log(gt))
                VarA00 = exp(log(EA00) + log(gt-r0)-log(gt) +log(gt-c0)-log(gt-1))
                pv = rep(NA, nrow(t))
                pv[do.tests] = 2*pnorm(-abs((t[do.tests,"n00"]-EA00)/sqrt(VarA00)))
                pv[pv==0] = .Machine$double.xmin
              }
      )
      
      
      if (correct) {
        all[[curr.res]] = c(all[[curr.res]], list(pv=pv, pv.correct=pv.correct))
      } else {
        all[[curr.res]] = c(all[[curr.res]], list(pv=pv))
      }
      
      selector = lapply(disc, `[[`, 2)
      names(selector) = W$cuboid.idx[unq]
      
      # Rank the tests that have p-values:
      p.rng = !is.na(all[[curr.res]]$pv)
      pvs.per.res = all[[curr.res]]$pv[p.rng]
      n.pvs.per.res = length(pvs.per.res)
      if (rk) {
        mM = min(n.pvs.per.res, M)
        cond = (mM>=n.pvs.per.res | res<R_star)
      } else {
        cond = (res<R_star)
      }
      
      if (cond) {
      # if (mM>=n.pvs.per.res | res<R_star) {
        rank.tests = rep(TRUE, n.pvs.per.res)
      } else {
        if (!rk) {
        rank.tests =
          as.logical(rep(FALSE, n.pvs.per.res) +
                       (pvs.per.res<=p_star))
        } else {
          rank.tests =
            as.logical(rep(FALSE, n.pvs.per.res) +
                         (pvs.per.res<=sort(pvs.per.res,partial=1:mM)[mM]))
        }
      }
      
      all[[curr.res]] = c(all[[curr.res]], list(rank.tests=rank.tests))
      
      parent.tests = rank.tests
      parent.tests[is.na(parent.tests)] = 0
      parent.tests = as.logical(parent.tests)
      num.parent.tests = sum(parent.tests)
      
      if (res>=R_star) {
        while (num.parent.tests>M) {
          rank.tests[rank.tests][tail(which(pvs.per.res[rank.tests]==max(pvs.per.res[rank.tests],na.rm=TRUE)),1)]=FALSE
          parent.tests = rank.tests
          num.parent.tests = sum(parent.tests,na.rm=TRUE)
        }
      }
      
      parent.cuboids = all[[curr.res]]$tables[p.rng,"cuboid.idx"][parent.tests]
      num.parent.cuboids = length(unique(parent.cuboids))
      all[[curr.res]] = c(all[[curr.res]], list(parent.cuboids=parent.cuboids, parent.tests=parent.tests))
      
      if (sum(rank.tests)==0) {
        if (verbose) cat("No potential parents in resolution", res, "have p-values below threshold.\n")
        break
      }
      
      res=res+1
    } # end res loop
  } # end res>0 if
  
  if (verbose) print(Sys.time()-s0)
  
  #---------------------------------------------------------------------
  
  if (verbose) cat("Individual tests completed, post-processing...\n")
  
  pvs = unlist(lapply(all, `[[`, "pv"), use.names=FALSE)
  p.rng = !is.na(pvs)
  p = pvs[p.rng]
  lp = length(p)
  min.p = min(p)
  if (correct) {
    pvs.correct = unlist(lapply(all, `[[`, "pv.correct"), use.names=FALSE)
    pv.correct = pvs.correct[p.rng]
    min.pv.correct = min(pv.correct)
  }
  if (return.all.pvs) {
    if (!compute.all.holm) holm.names=NULL
    holm.methods = c("H"%in%p.adjust.methods, "Hcorrected"%in%p.adjust.methods, "MH"%in%p.adjust.methods)
    if (compute.all.holm) holm.names=c("H","Hcorrected","MH")[holm.methods]
    all.pvs=as.data.frame(matrix(NA, nrow=length(pvs),
                                 ncol=2L+correct+length(holm.names)+
                                   length(setdiff(p.adjust.methods,c("H","Hcorrected","MH")))))
    if (correct) correct.p.name="pv.correct" else correct.p.name =NULL
    colnames(all.pvs) = c("test.idx", "pv", correct.p.name, holm.names, setdiff(p.adjust.methods,c("H","Hcorrected","MH")))
    all.pvs$test.idx = tx
    all.pvs$pv = pvs
    if (correct) all.pvs$pv.correct = pvs.correct
  }
  
  if (!is.null(p.adjust.methods)) {
    s1=Sys.time()
    # (The following code was adapted from the R
    #  packages 'discreteMTP' and 'MHTdiscrete')
    if ("H"%in%p.adjust.methods) {
      if (compute.all.holm) {
        try({
          if (verbose) cat("H: Computing all adjusted p-values...\n")
          pvs.H = rep(NA, length(pvs))
          s = seq_len(lp)
          o = order(p)
          ro = order(o)
          pvs.H[p.rng] = pmin(1, cummax((lp - s + 1L) * p[o]))[ro]
          if (return.all.pvs) all.pvs$H = pvs.H
          p.values["H"] = min(pvs.H[p.rng])
        })
      } else {
        if (verbose) cat("H: Computing global p-value...\n")
        p.values["H"] = min(1,lp*min.p)
      }
    }
    
    if ("Hcorrected"%in%p.adjust.methods) {
      if (compute.all.holm) {
        try({
          if (verbose) cat("Hcorrected: Computing all adjusted p-values...\n")
          pvs.Hcorrected = rep(NA, length(pvs))
          s = seq_len(lp)
          o = order(pv.correct)
          ro = order(o)
          pvs.Hcorrected[p.rng] = pmin(1, cummax((lp - s + 1L) * pv.correct[o]))[ro]
          if (return.all.pvs) all.pvs$Hcorrected = pvs.Hcorrected
          p.values["Hcorrected"] = min(pvs.Hcorrected[p.rng])
        })
      } else {
        if (verbose) cat("Hcorrected: Computing global p-value...\n")
        p.values["Hcorrected"] = min(1,lp*min.pv.correct)
      }
    }
    
    
    if ("MH"%in%p.adjust.methods |
        (Dx==1 & Dy==1 & uv.exact.null)) {
      ttables = do.call(rbind, lapply(all, `[[`, "tables"))
      ttables = ttables[p.rng,-c(1,2,7,8)]
      col0.tot = ttables[,"n00"]+ttables[,"n10"]
      col1.tot = ttables[,"n01"]+ttables[,"n11"]
      grand.tot = col0.tot + col1.tot
      row1.tot = ttables[,"n10"]+ttables[,"n11"]
      row0.tot = ttables[,"n00"]+ttables[,"n01"]
      min.margin = pmin(col0.tot, row0.tot)
    }
    
    # Compute MH:
    # Step 1: compute CDF for each unique test.
    if ("MH"%in%p.adjust.methods) {
      if (verbose) cat("MH: Computing CDF...\n")
      ppp = make_CDF( lp, col0.tot, row1.tot,
                      min.margin, ALL.PROBS, 
                      compute.all.holm, min.p )
      pCDFlist = lapply(ppp, `[[`, "pvals")
      if (!compute.all.holm) {
        if (verbose) cat("MH: Computing global p-value...\n")
        mp = vapply(ppp, `[[`, numeric(1), "mp")
      }
      
      lc = sum(p <= cutoff)
      o = order(p)
      ro = order(o)
      
      if (compute.all.holm) {
        if (verbose) cat("MH: Computing all adjusted p-values...\n")
        pvs.MH = rep(NA, length(pvs))
        sort.p.set = pCDFlist[o]
        sort.p = p[o][1L:lc]
        p.sums.MH = genStepFunCpp(sort.p.set, sort.p, lc, lp)
        p.sums.MH = c(pmin(1, p.sums.MH),rep(1,lp-lc))
        pvs.MH[p.rng] = cummax(p.sums.MH)[ro]
        if (return.all.pvs) all.pvs$MH[p.rng] = pvs.MH[p.rng]
        p.values["MH"] = min(pvs.MH[p.rng])
      }
      
      if ("MH"%in%p.adjust.methods & !compute.all.holm) p.values["MH"] = min(1,sum(mp))
    }
    
    if (verbose) print(Sys.time()-s1)
  }
  
  attributes(all)$parameters =
    list(Dx=Dx, Dy=Dy, n=n, M=M, R_max=R_max,
         R_star=R_star,
         rank.transform=rank.transform,
         test.method=test.method,
         correct=correct,
         min.tbl.tot=min.tbl.tot,
         min.row.tot=min.row.tot,
         min.col.tot=min.col.tot,
         p.adjust.methods=p.adjust.methods,
         compute.all.holm=compute.all.holm,
         cutoff=cutoff, top.max.ps=top.max.ps)
  
  if (return.all.pvs) attributes(all.pvs)$parameters = attributes(all)$parameters
  
  try(if (save.all.pvs) {cat("Saving all p-values\n");
    save(list(all=all, pvs=all.pvs), file=all.pvs.fname)})
  
  test.stats["avg.log.pv"] = sum(-log(p))/lp
  test.stats[paste0("avg.log.top",top.max.ps,".pv")] =
    sum(-log(sort(p, partial=1:min(lp,top.max.ps))[1:min(lp,top.max.ps)]))/min(lp,top.max.ps)
  if (correct) {
    test.stats["avg.log.pv.corrected"] = sum(-log(pv.correct))/lp
    test.stats[paste0("avg.log.top",top.max.ps,".pv.corrected")] =
      sum(-log(sort(pv.correct, partial=1:min(lp,top.max.ps))[1:min(lp,top.max.ps)]))/min(lp,top.max.ps)
  }
  
  if (Dx==1 & Dy==1) {
    if (uv.approx.null) {
      if (test.method=="Fisher" | test.method=="norm.approx") {
        uv.approx.nulls = uvApproxNull(uv.null.sim, sum(p.rng), top.max.ps, verbose)
        if (!correct) {
          p.values["pv.avg.approx"] =
            mean(uv.approx.nulls$log.avg>=test.stats["avg.log.pv"])
          p.values[paste0("pv.top",top.max.ps,".approx")] =
            mean(uv.approx.nulls$log.top>=test.stats[paste0("avg.log.top",top.max.ps,".pv")])
        } else {
          p.values["pv.avg.approx"] =
            mean(uv.approx.nulls$log.avg>=test.stats["avg.log.pv.corrected"])
          p.values[paste0("pv.top",top.max.ps,".approx")] =
            mean(uv.approx.nulls$log.top>=test.stats[paste0("avg.log.top",top.max.ps,".pv.corrected")])
        }
      }
    } # end uv.approx.null if
    
    if (uv.exact.null) {
      if (test.method=="Fisher" | test.method=="norm.approx") {
        uv.exact.nulls = uvExactNull(uv.null.sim, test.method, correct,
                                     col0.tot, row0.tot, grand.tot, top.max.ps,
                                     verbose)
        if (correct) {
          p.values["pv.avg.exact"] = mean(uv.exact.nulls$log.avg>test.stats["avg.log.pv.corrected"])
          p.values[paste0("pv.top",top.max.ps,".exact")] =
            mean(uv.exact.nulls$log.top>test.stats[paste0("avg.log.top",top.max.ps,".pv.corrected")])
        } else {
          p.values["pv.avg.exact"] = mean(uv.exact.nulls$log.avg>test.stats["avg.log.pv"])
          p.values[paste0("pv.top",top.max.ps,".exact")] =
            mean(uv.exact.nulls$log.top>test.stats[paste0("avg.log.top",top.max.ps,".pv")])
        }
      } else {
        warning("Function only generates approximate or theoretical null distributions in the univariate case, where normal approximation for the hypergeometric is used (Fisher's exact test or
                normal approximation).")
      } # end (correct) if-else
      } # end uv.exact.null if
    }
  
  if (verbose) {
    cat("\n")
    switch(test.method,
           Fisher = cat("Fisher's Exact Tests:\n"),
           norm.approx = cat("Two-tailed p-values from Normal Approximations to HG Distribution:\n"),
           chi.sq = cat("Chi Squared Tests:\n"),
           LR = cat("Likelihood Ratio Tests:\n")
    )
    cat(paste0("Mean of -log(p-values): ",signif(test.stats["avg.log.pv"],6),"\n"))
    cat(paste0("Mean of top ",top.max.ps, " -log(p-values): ",signif(test.stats[paste0("avg.log.top",top.max.ps,".pv")],6),"\n"))
    if (correct) {
      switch(test.method,
             Fisher = {correction.name = "mid-p"},
             chi.sq = {correction.name = "Yates"},
             LR = {correction.name = "Williams"})
      cat(paste0("Mean of -log(p-values with ",correction.name," correction): ",signif(test.stats["avg.log.pv.corrected"],6),"\n"))
      cat(paste0("Mean of top ",top.max.ps, " -log(p-values with ",correction.name," correction): ",signif(test.stats[paste0("avg.log.top",top.max.ps,".pv.corrected")],6),"\n"))
    }
    if (uv.exact.null & Dx==1 & Dy==1) {
      if (test.method=="Fisher" & correct) {
        cat(paste0("\nGlobal p-value, theoretical null, mean -log(p-values with mid-p correction): ",signif(p.values["pv.avg.exact"],6),"\n"))
        cat(paste0("Global p-value, theoretical null, mean top ",top.max.ps," -log(p-values with mid-p correction): ",signif(p.values[paste0("pv.top",top.max.ps,".exact")],6),"\n"))
      }
    }
    if (uv.approx.null & Dx==1 & Dy==1) {
      if (test.method=="Fisher" & correct) {
        cat(paste0("\nGlobal p-value, approximate null, mean -log(p-values with mid-p correction): ",signif(p.values["pv.avg.approx"],6),"\n"))
        cat(paste0("Global p-value, approximate null, mean top ",top.max.ps," -log(p-values with mid-p correction): ",signif(p.values[paste0("pv.top",top.max.ps,".approx")],6),"\n"))
      }
    }
    if (Dx==1 & Dy==1) {
      if (test.method!="Fisher" & uv.exact.null) {
        cat(paste0("\nGlobal p-value, exact null, mean -log(p-values): ",signif(p.values["pv.avg.exact"],6),"\n"))
        cat(paste0("Global p-value, exact null, mean top ",top.max.ps," -log(p-values): ",signif(p.values[paste0("pv.top",top.max.ps,".exact")],6),"\n"))
      }
      if (test.method!="Fisher" & uv.approx.null) {
        cat(paste0("\nGlobal p-value, approximate null, mean -log(p-values): ",signif(p.values["pv.avg.approx"],6),"\n"))
        cat(paste0("Global p-value, approximate null, mean top ",top.max.ps," -log(p-values): ",signif(p.values[paste0("pv.top",top.max.ps,".approx")],6),"\n"))
      }
    }
    cat("\n")
    if ("H"%in%p.adjust.methods) {
      cat(paste0("Global p-value, Holm on p-values: ",signif(p.values["H"],6),"\n"))
    }
    if ("Hcorrected"%in%p.adjust.methods) {
      cat(paste0("Global p-value, Holm on p-values with ",correction.name," correction: ",signif(p.values["Hcorrected"],6),"\n"))
    }
    if ("MH"%in%p.adjust.methods & test.method=="Fisher") {
      cat(paste0("Global p-value, Modified Holm step-down: ",signif(p.values["MH"],6),"\n"))
    }
  }
  
  val = list(test.stats=test.stats, all=all)
  if (!is.null(p.adjust.methods)) val = c(val, list(p.values=p.values))
  if (return.all.pvs) val = c(val, list(pvs=all.pvs))
  if (uv.approx.null & Dx==1 & Dy==1) val = c(val, list(approx.nulls=uv.approx.nulls))
  if (uv.exact.null  & Dx==1 & Dy==1) val = c(val, list(exact.nulls=uv.exact.nulls))
  return(invisible(val))
} # end multiFit

#==========================================================================
#' @title Approximate univariate null distribution
#' @description In a univariate case, simulate an approximate null distribution for the global
#' test statistics.
#'
#' @param uv.null.sim Positive integer, the number of simulated values to be computed
#' in a univariate case when an exact or approximate null distribution is simulated.
#' @param num.wins Positive integer, the number of windows that are tested in each
#' simulation.
#' @param top.max.ps Positive integer, report the mean of the top \code{top.max.ps}
#' order statistics of the negative logarithm of all \code{p}-values.
#' @param verbose Logical.
#' @return List of two numerical vectors for each of the global test statistics. Each
#' such vector is of length \code{uv.null.sim}.
#' @export
uvApproxNull = function(uv.null.sim, num.wins, top.max.ps=4L, verbose=FALSE) {
  t0=Sys.time()
  if (verbose) cat("Simulating an approximate null distribution...\n")
  approx.null = matrix(rexp(num.wins*uv.null.sim),nrow=num.wins)
  approx.log.avg = apply(approx.null,2,mean)
  approx.log.top = apply(-apply(-approx.null,2,sort,partial=1:top.max.ps)[1:top.max.ps,],2,mean)
  if (verbose) print(Sys.time()-t0)
  return(list(log.avg=approx.log.avg,
              log.top=approx.log.top))
} # end uvApproxNull

#==========================================================================
#' @title Exact univariate null distribution
#' @description In a univariate case, simulate an exact null distribution for the global
#' test statistics.
#'
#' @param uv.null.sim Positive integer, the number of simulated values to be computed
#' in a univariate case when an exact or approximate null distribution is simulated.
#' @param test.method String, "Fisher" and "norm.approx" are applicable here.
#' @param correct Logical, relates to Fisher's exact test, if \code{TRUE} the exact
#' null is simulated for the mid-\code{p} corrected \code{p}-values.
#' @param col0.tot Numerical vector, containing the column-0 totals of all 2x2 contingency
#' tables for which simulated values are asked for. Has to be the same length as \code{row0.tot}
#' and \code{grand.tot}.
#' @param row0.tot Numerical vector, containing the row-0 totals of all 2x2 contingency
#' tables for which simulated values are asked for. Has to be the same length as \code{col0.tot}
#' and \code{grand.tot}.
#' @param grand.tot Numerical vector, containing the totals of all 2x2 contingency
#' tables for which simulated values are asked for. Has to be the same length as \code{col0.tot}
#' and \code{row0.tot}.
#' @param top.max.ps Positive integer, report the mean of the top \code{top.max.ps}
#' order statistics of the negative logarithm of all \code{p}-values.
#' @param verbose Logical.
#' @return List of two numerical vectors for each of the global test statistics. Each
#' such vector is of length \code{uv.null.sim}.
#' @export
uvExactNull = function(uv.null.sim, test.method, correct=TRUE,
                         col0.tot, row0.tot, grand.tot, top.max.ps=4L,
                         verbose=FALSE) {
  t1=Sys.time()
  row1.tot=grand.tot-row0.tot
  col1.tot=grand.tot-col0.tot
  if (verbose) cat("Simulating from the theoretical null distribution...\n")
  if (test.method=="norm.approx") {
    EA00 = col0.tot * row0.tot / grand.tot
    VarA00 =  EA00 * ((grand.tot-row0.tot)/grand.tot) *
      ((grand.tot-col0.tot)/(grand.tot-1))
  }
  min.margin = pmin(col0.tot, row0.tot)
  sample.win = function(s) {
    A00 = rhyper(uv.null.sim, row0.tot[s], row1.tot[s], col0.tot[s])
    lo = max(0L, col0.tot[s]-row1.tot[s])
    hi = min.margin[s]
    support = lo:hi
    switch (test.method,
            Fisher = {
              all.probs = { d = dhyper(support, row0.tot[s], row1.tot[s], col0.tot[s], log=TRUE)
              d = exp(d - max(d))
              d/sum(d) }
              if (!correct) {
              pvals = sapply(A00,
                             function(a00) { relErr=1 + 10^(-7);
                             pv=sum(all.probs[all.probs<=all.probs[a00-lo+1]*relErr])
                             ifelse(pv==0, .Machine$double.xmin, pv)
                             })
              }
              if (correct) {
                pvals = sapply(A00,
                               function(a00) { relErr=1 + 10^(-7);
                               pv.correct=(sum(all.probs[all.probs<=all.probs[a00-lo+1]*relErr])
                                           +sum(all.probs[all.probs<all.probs[a00-lo+1]*relErr]))/2
                               ifelse(pv.correct==0, .Machine$double.xmin, pv.correct)
                               })
              }
            },
            norm.approx = {
              all.probs = 2*pnorm(-abs((support-EA00[s])/sqrt(VarA00[s])))
              pvals = sapply(A00,
                             function(a00) { 2*pnorm(-abs((a00-EA00[s])/sqrt(VarA00[s]))) })
            })
    return(-log(pvals))
  }
  sample.wins = sapply(1:length(grand.tot), sample.win)
  exact.log.avg = apply(sample.wins,1,mean)
  exact.log.top = apply(t(-apply(-sample.wins,1,sort,partial=1:top.max.ps))[,1:top.max.ps],1,mean)
  if (verbose) print(Sys.time()-t1)
  return(list(log.avg=exact.log.avg,
              log.top=exact.log.top))
}

#==========================================================================
#' @title Permutation null distribution
#' @description Simulate a permutation null distribution (per input data and testing
#' parameters) for the global test statistics. See vignettes for examples of usage.
#'
#' @param perm.null.sim Positive integer, the number of simulated values to be computed
#' when a permutation null distribution is simulated.
#' @param xy A list, whose first element corresponds to the matrix x as below, and
#' its second element corresponds to the matrix y as below. If \code{xy} is not
#' specified, \code{x} and \code{y} need to be assigned.
#' @param x A matrix, number of columns = dimension of random vector,
#' number of rows = number of observations.
#' @param y A matrix, number of columns = dimension of random vector,
#' number of rows = number of observations.
#' @param fit An object generated by \code{multiFit}.
#' @param parts Positive integer, divide computation to four parts. Useful for getting
#' a sense of progress (when \code{verbose=TRUE}) or for distributing memory requirement.
#' @param save.perm.null Logical, if \code{TRUE}, save the permutation null into an RData
#' file named \code{perm.null.fname}
#' @param perm.null.fname String, file name to which to save an RData file containing
#' the permutation null if \code{save.perm.null=TRUE}.
#' @param verbose Logical.
#' @export
permNullTest = function(perm.null.sim=10000L, xy, x=x, y=y, fit,
                        parts=4L, save.perm.null=FALSE,
                        perm.null.fname=NULL, verbose=FALSE) {
  s0=Sys.time()
  if (!missing(xy)) {
    x=xy[[1]]
    y=xy[[2]]
  }
  a = attributes(fit$all)
  Dx=a$parameters$Dx
  Dy=a$parameters$Dy
  n=a$parameters$n
  M=a$parameters$M
  top.max.ps=a$parameters$top.max.ps
  R_max=a$parameters$R_max
  R_star=a$parameters$R_star
  rank.transform=a$parameters$rank.transform
  test.method=a$parameters$test.method
  correct=a$parameters$correct
  min.tbl.tot=a$parameters$min.tbl.tot
  min.row.tot=a$parameters$min.row.tot
  min.col.tot=a$parameters$min.col.tot

  p.adjust.methods=NULL

  if (perm.null.sim<100) { parts=1 }

  names.test.stats=NULL
  if (rank.transform | (Dx==1 & Dy==1)) names.test.stats = c("avg.log.pv", paste0("avg.log.top",top.max.ps,".pv"))
  if ((rank.transform | (Dx==1 & Dy==1)) & correct) {
    names.test.stats=c(names.test.stats, "avg.log.pv.corrected",
                       paste0("avg.log.top",top.max.ps,".pv.corrected"))
  }
  if (length(names.test.stats)>0) {
    test.stats = rep(NA, length(names.test.stats))
    names(test.stats) = names.test.stats
  }

  null.dist=matrix(as.numeric(NA), nrow=perm.null.sim, ncol=length(test.stats))
  colnames(null.dist) = names.test.stats

  parts.vec = c(floor(perm.null.sim/parts)+perm.null.sim-floor(perm.null.sim/parts)*parts,rep(floor(perm.null.sim/parts),parts-1))
  if (verbose) cat("\nSimulating ",perm.null.sim," observations from a permutation null distribution in ",parts," part",ifelse(parts>1,"s",""),"...\n",sep="")
  for (prt in 1:parts) {
    if (verbose) {
      if (prt==1) cat(paste0("Running simulations ",1," to ",parts.vec[1],"\n"))
      if (prt>1) cat(paste0("Running simulations ",(1+sum(parts.vec[1:(prt-1)])),
                            " to ",(sum(parts.vec[1:prt])),"\n"))
    }
    s0=Sys.time()
    l = t(sapply( 1:parts.vec[prt], function(sim) 
                        multiFit(x=x[sample(n),], y=y, M=M,
                                 R_star=R_star, R_max=R_max,
                                 test.method=test.method, correct=correct,
                                 rank.transform=rank.transform, p.adjust.methods=NULL,
                                 return.all.pvs=FALSE,
                                 min.tbl.tot=min.tbl.tot, min.row.tot=min.row.tot,
                                 min.col.tot=min.col.tot, top.max.ps=top.max.ps,
                                 verbose=FALSE)$test.stats ))
    s1=Sys.time()
    if (verbose) print(s1-s0)
    if (prt==1) {
      null.dist[1:parts.vec[prt],]=l
    } else {
      null.dist[(1+sum(parts.vec[1:(prt-1)])):sum(parts.vec[1:prt]),]=l
    }
  }

  try(if (save.perm.null) save(null.dist, file=perm.null.fname) )

  p = sapply(names(fit$test.stats), function(nm) mean(null.dist[,nm]>fit$test.stats[nm]))
  if (!correct) names(p) = c("pv.avg", paste0("pv.top",top.max.ps))
  if (correct) names(p) = c("pv.avg", paste0("pv.top",top.max.ps),"pv.avg.correct", paste0("pv.top",top.max.ps,".correct"))
  if (verbose) {
    cat(paste0("\nGlobal p-value, permutation null, mean -log(p-values): ",p["pv.avg"],"\n"))
    cat(paste0("Global p-value, permutation null, mean top ",top.max.ps," -log(p-values): ",p[paste0("pv.top",top.max.ps)],"\n"))
    if (correct) {
      cat(paste0("\nGlobal p-value, permutation null, mean -log(corrected p-values): ",p["pv.avg.correct"],"\n"))
      cat(paste0("Global p-value, permutation null, mean top ",top.max.ps," -log(p-values with mid-p correction): ",p[paste0("pv.top",top.max.ps,".correct")],"\n"))
    }
  }
  names(p) = paste0(names(p),".","permutation")

  return(list(null.dist=null.dist, p.values=p))
}
