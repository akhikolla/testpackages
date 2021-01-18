tolatt2Ns2 <- function(parvec) { ##util. to convert parvec : c(2Nmu, 2Nm, g) to c(2Nmu, 2Ds2[lattice])
  if(inherits(parvec,c("data.frame","list"))) parvec <- unlist(parvec)
  ##tmp <- condaxialS2fromg(parvec["g"])*parvec["twoNm"]
  ##return(as.numeric(c(parvec["twoNmu"], tmp)))
  return(c(twoNmu=parvec[["twoNmu"]],
           latt2Ns2 = condaxialS2fromg(parvec[["g"]], D2bool=("2D" %in% blackbox.getOption("DemographicModel"))) * parvec[["twoNm"]] ))
}


toDgmuFromCanonical <- function(parvec) { ##utility to convert parvec c(twoNmu, D,T,..., twoNancmu) to c(TwoNmu,Dgmu,T, ..., twoNancmu)
  if(inherits(parvec,c("data.frame","list"))) parvec <- unlist(parvec) ## keep names
  return(c(Dgmu=parvec[["D"]]*parvec[["twoNmu"]]))
}

tom1overmuFromCanonical <- function(parvec) { ##utility to convert parvec c(twoNmu,Q1,M1,M1) to c(TwoNmu,Q1,m1overmun,M2)
  if(inherits(parvec,c("data.frame","list"))) parvec <- unlist(parvec) ## keep names
  return(c(m1overmu=parvec[["M1"]]/parvec[["twoNmu"]]/parvec[["Q1"]]))
}

tom2overmuFromCanonical <- function(parvec) { ##utility to convert parvec c(twoNmu,Q1,M1,M2,...) to c(TwoNmu,Q1,M1,m2overmu,...)
  if(inherits(parvec,c("data.frame","list"))) parvec <- unlist(parvec) ## keep names
  return(c(m2overmu=parvec[["M2"]]/parvec[["twoNmu"]]/(1.0-parvec[["Q1"]])))
}

tomratioFromCanonical <- function(parvec) { ##utility to convert parvec c(twoNmu,Q1,M1,M2,...) to c(twoNmu,Q1,mratio,M2,...)
  if(inherits(parvec,c("data.frame","list"))) parvec <- unlist(parvec) ## keep names
  return(c(mratio=parvec[["M1"]]/parvec[["M2"]]*(1.0-parvec[["Q1"]])/parvec[["Q1"]]))
}

toNactNfounderratioFromCanonical <- function(parvec) { ##utility to convert parvec c(twoNmu, ..., twoNfoundermu, twoNancmu) to c(NactNfounderratio, ..., twoNfoundermu, twoNancmu)
  if(inherits(parvec,c("data.frame","list"))) parvec <- unlist(parvec) ## keep names
  return(c(NactNfounderratio=parvec[["twoNmu"]]/parvec[["twoNfoundermu"]]))
}

toNfounderNancratioFromCanonical <- function(parvec) { ##utility to convert parvec c(twoNmu, ..., twoNfoundermu, twoNancmu) to c(NfounderNancratio, ..., twoNfoundermu, twoNancmu)
  if(inherits(parvec,c("data.frame","list"))) parvec <- unlist(parvec) ## keep names
  return(c(NfounderNancratio=parvec[["twoNfoundermu"]]/parvec[["twoNancmu"]]))
}

toNMratioFromCanonical <- function(parvec) { ##utility to convert parvec c(twoNmu,Q1,M1,M2,...) to c(twoNmu,Q1,NMratio,M2,...)
  if(inherits(parvec,c("data.frame","list"))) parvec <- unlist(parvec) ## keep names
  return(c(NMratio=parvec[["M1"]]/parvec[["M2"]]))
}

toNratioFromCanonical <- function(parvec) { ##utility to convert parvec c(twoNmu, ..., twoNancmu) to c(Nratio, ..., twoNancmu)
  if(inherits(parvec,c("data.frame","list"))) parvec <- unlist(parvec) ## keep names
  return(c(Nratio=parvec[["twoNmu"]]/parvec[["twoNancmu"]]))
}

toTgmuFromCanonical <- function(parvec) { ##utility to convert parvec c(twoNmu,D,T, ..., twoNancmu) to c(TwoNmu,D,Tgmu ..., twoNancmu)
  if(inherits(parvec,c("data.frame","list"))) parvec <- unlist(parvec) ## keep names
  return(c(Tgmu=parvec[["T"]]*parvec[["twoNmu"]]))
}

