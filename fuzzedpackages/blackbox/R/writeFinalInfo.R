writeFinalInfo <- function(cleanResu="") {
  rosglobal <- blackbox.getOption("rosglobal")
  plotOptions <- blackbox.getOption("plotOptions")
  oneDimCIvars <- blackbox.getOption("oneDimCIvars")
  message.redef("...done.")
  message.redef("\n*** Final likelihood estimates, and predicted logL: ***")
  message.redef(prettyNamedUserValues(c(rosglobal$canonVP, "ln(L)"=rosglobal$value), extradigits=2))
  returncode <- rosglobal$convergence
  tmp <- rosglobal$edgelevel
  if (tmp>0) returncode <- returncode+tmp/(10^ceiling(log(tmp, 10))) ## second summand goes in decimal part of returcode
  write("\n*** Point estimates *** \n", file=cleanResu)
  writeCleanNAMED(prettyNamedUserValues(rosglobal$canonVP), file=cleanResu)
  DemographicModel <- blackbox.getOption("DemographicModel")
  if ("IBD" %in% DemographicModel) {
    D1IBDbool <- "1D" %in% DemographicModel
    Nbfactor <- blackbox.getOption("Nbfactor")
    write(paste("\n      Neighborhood: ", prettynum(rosglobal$latt2Ns2*Nbfactor), " ",
                if (D1IBDbool) {blackbox.getOption("GeoUnit")} else {""}, sep=""), file=cleanResu)
    if (D1IBDbool) write(paste("## Conversion factor from Nb in lattice units to Nb in geographic distance units as deduced from input:", Nbfactor),
                         file=blackbox.getOption("estimOutf"))
  } else if (length(intersect(DemographicModel, c("OnePopVarSize", "IM")))>0) {
    write(paste("\n      N ratio: ", prettynum(rosglobal$Nratio), " ", sep=""), file=cleanResu)
  } else if ("OnePopFounderFlush" %in% DemographicModel) {
    write(paste("\n      Nanc ratio: ", prettynum(rosglobal$Nratio), " ", sep=""), file=cleanResu)
    write(paste("\n      NactNfounder ratio: ", prettynum(rosglobal$NactNfounderratio), " ", sep=""), file=cleanResu)
    write(paste("\n      NfounderNanc ratio: ", prettynum(rosglobal$NfounderNancratio), " ", sep=""), file=cleanResu)
  }
  if (length(intersect(DemographicModel, c("OnePopVarSize", "OnePopFounderFlush", "IM")))>0) {
    if ( !("IM" %in% DemographicModel) && ( ("DgmuProf" %innc% plotOptions) || ("Dgmu" %innc% oneDimCIvars) ) ) write(paste("\n      Dg*mu: ", prettynum(rosglobal$Dgmu), " ", sep=""), file=cleanResu)
    if ( ("TgmuProf" %innc% plotOptions) || ("Tgmu" %innc% oneDimCIvars) ) write(paste("\n      Tg*mu: ", prettynum(rosglobal$Tgmu), " ", sep=""), file=cleanResu)
  }
  if (length(intersect(DemographicModel, c("Npop", "IM")))>0) {
    write(paste("\n      NMratio: ", prettynum(rosglobal$NMratio), " ", sep=""), file=cleanResu)
    write(paste("\n      mratio: ", prettynum(rosglobal$mratio), " ", sep=""), file=cleanResu)
    if (  ("movermuProf" %innc% plotOptions) || ("m1overmu" %innc% oneDimCIvars) ) write(paste("\n      m1/mu: ", prettynum(rosglobal$m1overmu), " ", sep=""), file=cleanResu)
    if ( ("movermuProf" %innc% plotOptions) || ("m2overmu" %innc% oneDimCIvars) ) write(paste("\n      m2/mu: ", prettynum(rosglobal$m2overmu), " ", sep=""), file=cleanResu)
  }
  ## note that the C codes seeks estimates in the VERY LAST line of the output.txt file: do not write comments after the following output:
  upperPred_crits <- blackbox.getOption("upperPred_crits")
  if (is.null(upperPred_crits)) upperPred_crits <- list(RMSpred=NA,GOP=NA) ## if sampleByResp has not been run (which should not occur)
  writeoutput(paste(blackbox.getOption("dataFile"), "(final)", sep=""),
              returncode=returncode,
              NA,
              upperPred_crits$RMSpred, ## RMSpred for upper points only
              upperPred_crits$GOP) ## GOP root mean relative squared error of prediction
  if ( ! blackbox.getOption("interactiveGraphics")) { ##
    plotFiles <- blackbox.getOption("plotFiles")
    if(!is.null(plotFiles) & length(plotFiles)>0)
      message.redef(paste("See file(s) ", paste(names(plotFiles), sep="", collapse=", "), " for figures", sep=""))
    graphics.off() ## close all graphic files, but note that plotFiles is not cleaned.
    close(blackbox.getOption("estimOutf")) ## ! leaves connection open if interactive run !
  }
  write("\nNormal ending.", file=cleanResu)
  close(cleanResu) ##
  blackbox.options(cleanResu="") ## so that all 'cleanResu' output now goes to the standard output connection; see ?write
  alarm() ## Sounds when execution of the R source file reaches this point. Except in Rstudio...
  invisible(NULL)
}
