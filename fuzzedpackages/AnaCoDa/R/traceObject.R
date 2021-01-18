#### TODO, lets move it into parameterObject.R and use a parameter instead of trace. thats how it is done for the acf function

# see mcmc Object.R convergence.test function for documentation
convergence.test.Rcpp_Trace <- function(object, samples = 10, frac1 = 0.1,
                                        frac2 = 0.5, thin = 1, plot = FALSE, what = "Mutation", mixture = 1)
{
  current.trace <- 0
  if(what[1] == "Mutation" || what[1] == "Selection")
  {
    names.aa <- aminoAcids()
    numCodons <- 0
    for(aa in names.aa)
    {
      if (aa == "M" || aa == "W" || aa == "X") next
      codons <- AAToCodon(aa, T)
      numCodons <- numCodons + length(codons)
    }
    index <- 1
    cur.trace <- vector("list", numCodons)
    for(aa in names.aa)
    {
      if (aa == "M" || aa == "W" || aa == "X") next
      codons <- AAToCodon(aa, T)
      for(i in 1:length(codons))
      {
        if(what[1] == "Mutation"){
          cur.trace[[index]]<- object$getCodonSpecificParameterTraceByMixtureElementForCodon(mixture, codons[i], 0, T)
        }else{
          cur.trace[[index]] <- object$getCodonSpecificParameterTraceByMixtureElementForCodon(mixture, codons[i], 1, T)
        }
        index <- index + 1
      }
    }
    current.trace <- do.call("rbind", cur.trace)
    ## Transpose matrix to get in correct format for coda::mcmc. Transposing results in same output from coda::geweke.test as performing the test separately on each codon specific parameter
    current.trace <- t(current.trace)
  }
  
  if(what[1] == "Alpha" || what[1] == "Lambda" || what[1] == "NSERate" || what[1] == "LambdaPrime")
  {
    codon.list <- codons()
    codon.list <- codon.list[1:(length(codon.list)-3)]
    cur.trace <- vector("list",length(codon.list))
    for (i in 1:length(codon.list))
    {
      if (what[1]=="Alpha")
      {
        cur.trace[[i]]<- object$getCodonSpecificParameterTraceByMixtureElementForCodon(mixture, codon.list[i], 0, F)
      } else if (what[1]=="Lambda" || what[1]=="LambdaPrime"){
        cur.trace[[i]]<- object$getCodonSpecificParameterTraceByMixtureElementForCodon(mixture, codon.list[i], 1, F)
      } else if (what[1]=="NSERate"){
        cur.trace[[i]]<- object$getCodonSpecificParameterTraceByMixtureElementForCodon(mixture, codon.list[i], 2, F)
      }
    }
    current.trace <- do.call("rbind", cur.trace)
    ## Transpose matrix to get in correct format for coda::mcmc. Transposing results in same output from coda::geweke.test as performing the test separately on each codon specific parameter
    current.trace <- t(current.trace)
  }
  
  if(what[1] == "MixtureProbability")
  {
    numMixtures <- object$getNumberOfMixtures()
    cur.trace <- vector("list", numMixtures)
    for(i in 1:numMixtures)
    {
      cur.trace[[i]] <- object$getMixtureProbabilitiesTraceForMixture(i)
    }
    current.trace <- do.call("rbind", cur.trace)
    current.trace <- t(current.trace)
  }
  if(what[1] == "Sphi")
  {
    sphi <- object$getStdDevSynthesisRateTraces()
    current.trace <- do.call("rbind", sphi)
    current.trace <- t(current.trace)
  }
  if(what[1] == "Mphi")
  {
    sphi <- object$getStdDevSynthesisRateTraces()
    sphi <- do.call("rbind", sphi)
    mphi <- -(sphi * sphi) / 2;
    current.trace <- t(mphi)
  }
  if(what[1] == "Aphi")
  {
    # TODO need way to determine number of Aphi traces
  }
  if(what[1] == "Sepsilon")
  {
    # TODO need way to determine number of Sepsilon traces
  }
  if(what[1] == "ExpectedPhi")
  {
    current.trace <- object$getExpectedSynthesisRateTrace()
  }
  if(what[1] == "InitiationCost")
  {
    current.trace <- object$getInitiationCostTrace()
  }
  if(what[1] == "Expression")
  {
    # TODO need way to determine number of expression traces
  }
  if(what[1] == "AcceptanceCSP")
  {
    names.aa <- aminoAcids()
    index <- 1
    cur.trace <- vector("list", length(names.aa) - length(c("M","W","X")))
    for(aa in names.aa)
    {
      if (aa == "M" || aa == "W" || aa == "X") next
      cur.trace[[index]] <- object$getCodonSpecificAcceptanceRateTraceForAA(aa)
      index <- index + 1
    }
    current.trace <- do.call("rbind", cur.trace)
    current.trace <- t(current.trace)
  }
  trace.length <- length(current.trace)
  start <- max(0, trace.length - samples)
  mcmcobj <- coda::mcmc(data=current.trace, start=start, thin=thin)
  if(plot){
    coda::geweke.plot(mcmcobj, frac1=frac1, frac2=frac2)
  } else{
    diag <- coda::geweke.diag(mcmcobj, frac1=frac1, frac2=frac2)
    return(diag)
  }
}
