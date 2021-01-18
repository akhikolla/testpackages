# CPFunChisq.R -- Comparative functional chisquare test
#
# MS
# Created: May 3, 2016. Extracted from FunChisq.R

cp.fun.chisq.test <- function(x,
                              method=c("fchisq", "nfchisq", "default", "normalized"),
                              log.p=FALSE)
{
  if(mode(x)!="list" || length(x)<2 ) {
    stop("only accept list of 2 or more matrices as input!")
  }

  method <- match.arg(method)
  if(method == "default") {
    warning(paste0("method=\"", method, "\" is deprecated. Use \"fchisq\" instead."))
    method <- "fchisq"
  } else if(method == "normalized") {
    warning(paste0("method=\"", method, "\" is deprecated. Use \"nfchisq\" instead."))
    method <- "nfchisq"
  }

  finalStat <- 0
  finalDf <- 0

  for(i in 1:nrow(x[[1]]))
  {
    oneT <- c() # one table for each row
    for(j in 1:length(x))
    {
      oneT <- rbind(oneT, x[[j]][i,])
    }
    oneresult <- fun.chisq.test(oneT)
    finalStat <- finalStat + oneresult$statistic
    finalDf <- finalDf + oneresult$parameter
  }

  DNAME <- deparse(substitute(x))

  if(method=="fchisq") {

    names(finalStat) <- "statistic"
    names(finalDf) <- "parameter"
    p.value <- pchisq(finalStat, df = finalDf, lower.tail=FALSE, log.p=log.p)
    return(structure(list( statistic=finalStat, parameter=finalDf, p.value=p.value,
                           method = "Comparative functional chi-square test for heterogeneity",
                           data.name= DNAME),
                     class = "htest"))

  } else if(method=="nfchisq") {

    finalStat <- as.numeric((finalStat-finalDf)/sqrt(2*finalDf))
    names(finalStat) <- "statistic"
    names(finalDf) <- "parameter"
    return(structure(list(statistic = finalStat, parameter = finalDf,
                          p.value = pnorm( finalStat, lower.tail=FALSE, log.p=log.p),
                          method = "Nomalized comparative functional chi-square test for heterogeneity",
                          data.name= DNAME),
                     class = "htest"))

  }
}
