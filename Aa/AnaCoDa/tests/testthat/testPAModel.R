# library(testthat)
# library(AnaCoDa)
# 
# context("PA Model")
# 
# test_that("PA Model testing simulated versus actual accuracy", {
#   # Skip unless manually run or changed
# #  if (F)
#     skip("PA Model testing is optional.")
#   
#   #####################
#   ### Initial Setup ###
#   #####################
# 
#   # Test with only 1500 genes
#   fileName = file.path("UnitTestingData", "testPAModelFiles", "testPAModel.csv")
#   fileTable = file.path("UnitTestingData", "testPAModelFiles", "codonTranslationRates.csv")
#   
#   # Ensure the input files exist.
#   test_that("file exists: testPAModel.csv", {
#     expect_equal(file.exists(fileName), T)
#   })
#   
#   test_that("file exists: codonTranslationRates.csv", {
#     expect_equal(file.exists(fileTable), T)
#   })
#   
#   genome <- initializeGenomeObject(file = fileName, fasta = FALSE, append = FALSE)
#   
#   sphi_init <- c(2)
#   numMixtures <- 1
#   mixDef <- "allUnique"
#   geneAssignment <- c(rep(1, length(genome))) 
#   parameter <- initializeParameterObject(genome, sphi_init, numMixtures, geneAssignment, model= "PA", split.serine = TRUE, mixture.definition = mixDef)
#   #parameter <- initializeParameterObject(model="PA", restart.file="30restartFile.rst")
#   
#   samples <- 20000
#   thinning <- 10
#   adaptiveWidth <- 10
#   mcmc <- initializeMCMCObject(samples = samples, thinning = thinning, adaptive.width = adaptiveWidth, 
#                                est.expression=TRUE, est.csp=TRUE, est.hyper=TRUE)
#   
#   model <- initializeModelObject(parameter, "PA")
#   setRestartSettings(mcmc, "restartFile.rst", adaptiveWidth, TRUE)
#   
#   outFile = file.path("UnitTestingOut", "testPAModelLog20000.txt")
#   
#   sink(outFile)
#   system.time(
#     runMCMC(mcmc, genome, model, 9)
#   )
#   sink()
#   
#   #########################################################################################
#   ### Output File 1: MCMC, Log Posterior, Mixture Probability, Mphi, Sphi, Expected Phi ###
#   #########################################################################################
#   
#   # plots different aspects of trace
#   trace <- parameter$getTraceObject()
#   writeParameterObject(parameter, file = file.path("UnitTestingOut", "PAObject.Rdat"))
#   writeMCMCObject(mcmc, file = file.path("UnitTestingOut", "MCMCObject.Rdat"))
#   
#   
#   pdf(file.path("UnitTestingOut", "RFP_Genome_allUnique_startCSP_startPhi_adaptSphi_true.pdf"))
#   plot(mcmc, main = "MCMC Trace") #plots the whole log posterior trace
#   
#   
#   # take a subset of the trace values for the log posterior trace and plot them.
#   # The primary reason for doing this is the "jump" that throws the scale of the graph
#   # at the beginning is removed by taking out the beginning values.
#   logPost.trace <- mcmc$getLogPosteriorTrace()
#   start <- length(logPost.trace) * 0.7 
#   # the multiplier (currently 0.7) determines how much of the beginning trace is eliminated.
#   
#   logL <- logL <- mean(logPost.trace[start:length(logPost.trace)]) #get the mean for the subset
#   plot(logPost.trace[start:length(logPost.trace)], type="l", main=paste("logL:", logL), xlab="Sample", ylab="log(Posterior)")
#   grid (NULL,NULL, lty = 6, col = "cornsilk2")
#   
#   
#   plot(trace, what = "MixtureProbability", main = "Mixture Probability")
#   plot(trace, what = "Mphi", main = "Mphi")
#   plot(trace, what = "Sphi", main = "Sphi")
#   plot(trace, what = "ExpectedPhi", main = "Expected Phi")
#   
#   #logPost.trace <- mcmc$getLogPosteriorTrace()
#   #acf(logPost.trace)
#   #acf(logPost.trace[start:length(logPost.trace)])
#   dev.off()
#   
#   #################################
#   ### Output File 2: CSP Traces ###
#   #################################
#   
#   pdf(file.path("UnitTestingOut", "RFP_CSP_Values_Mixture1.pdf"), width = 11, height = 20)
#   #plot(trace, what = "Expression", geneIndex = 905) #used to make sure gene stabalized, not really needed now
#   plot(trace, what = "Alpha", mixture = 1)
#   plot(trace, what = "LambdaPrime", mixture = 1)
#   dev.off()
#   
#   #####################
#   ### Output File 3 ###
#   #####################
#   
#   pdf(file.path("UnitTestingOut", "RFPConfidenceIntervalsForAlphaAndLambdaPrime.pdf"))
#   
#   #eventually this will need loop over all categories if there are multiple mixtures
#   cat <- 1
#   proposal <- FALSE
#   alphaList <- numeric (61)
#   lambdaPrimeList <- numeric (61)
#   waitingTimes <- numeric(61)
#   alpha.ci <- matrix(0, ncol=2, nrow=61)
#   lambdaPrime.ci <- matrix(0, ncol=2, nrow=61)
#   psiList <- numeric(length(genome))
#   ids <- numeric(length(genome))
#   codonList <- codons()
#   
#   waitRates <- numeric(61)
#   for (i in 1:61)
#   {
#     codon <- codonList[i]
#     alphaList[i] <- parameter$getCodonSpecificPosteriorMean(cat, samples * 0.5, codon, 0, FALSE)
#     alphaTrace <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(1, codon, 0, FALSE)
#     alpha.ci[i,] <- quantile(alphaTrace[(samples * 0.5):samples], probs = c(0.025,0.975))
#     
#     
#     lambdaPrimeList[i] <- parameter$getCodonSpecificPosteriorMean(cat, samples * 0.5, codon, 1, FALSE)
#     lambdaPrimeTrace <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(1, codon, 1, FALSE)
#     lambdaPrime.ci[i,] <- quantile(lambdaPrimeTrace[(samples * 0.5):samples], probs = c(0.025,0.975))
#     waitingTimes[i] <- alphaList[i] * lambdaPrimeList[i]
#     waitRates[i] <- (1.0/waitingTimes[i])
#   }
#   
#   for (geneIndex in 1:length(genome))
#   {
#     psiList[geneIndex] <- parameter$getSynthesisRatePosteriorMeanByMixtureElementForGene(samples * 0.5, geneIndex, 1)
#     
#     g <- genome$getGeneByIndex(geneIndex, F)
#     ids[geneIndex] <- g$id
#   }
#   
#   #Plot confidence intervals for alpha and lambda prime
#   plot(NULL, NULL, xlim=range(1:61, na.rm = T), ylim=range(alpha.ci), 
#        main = "Confidence Intervals for Alpha Parameter", xlab = "Codons", 
#        ylab = "Estimated values", axes=F) 
#   confidenceInterval.plot(x = 1:61, y = alphaList, sd.y = alpha.ci)
#   axis(2)
#   axis(1, tck = 0.02, labels = codonList[1:61], at=1:61, las=2, cex.axis=.6)
#   
#   plot(NULL, NULL, xlim=range(1:61, na.rm = T), ylim=range(lambdaPrime.ci), 
#        main = "Confidence Intervals for LambdaPrime Parameter", xlab = "Codons", 
#        ylab = "Estimated values", axes=F) 
#   confidenceInterval.plot(x = 1:61, y = lambdaPrimeList, sd.y = lambdaPrime.ci)
#   axis(2)
#   axis(1, tck = 0.02, labels = codonList[1:61], at=1:61, las=2, cex.axis=.6)
#   
#   
#   # correlation between PAModel and Pop's wait rates
#   # load Pop's data
#   X <- read.csv(fileTable)
#   X <- X[order(X[,1]) , ]
#   XM <- matrix(c(X[,1], X[,2]), ncol = 2, byrow = FALSE)
#   
#   
#   Y <- data.frame(codonList[-c(62,63,64)], waitRates)
#   colnames(Y) <- c("Codon", "PausingTimeRates")
#   Y <- Y[order(Y[,1]) , ]
#   
#   
#   plot(NULL, NULL, xlim=range(XM[,2], na.rm = T), ylim=range(Y[,2]), 
#        main = "Correlation Between Pop and PA Model Pausing Time Rates", xlab = "Pop's Rates", ylab = "PA's Rates")
#   upper.panel.plot(XM[,2], Y[,2])
#   
#   # correlation between PAModel WAIT RATES (inverse) and Pop's wait rates
#   Y <- data.frame(codonList[-c(62,63,64)], waitingTimes)
#   colnames(Y) <- c("Codon", "WaitingTimeRates")
#   Y <- Y[order(Y[,1]) , ]
#   
#   
#   plot(NULL, NULL, xlim=range(XM[,2], na.rm = T), ylim=range(Y[,2]), 
#        main = "Correlation Between Pop and PA Model Waiting Times", xlab = "Pop's Rates", ylab = "PA's Waiting Times")
#   upper.panel.plot(XM[,2], Y[,2])
# 
#   dev.off()
#   
#   #################
#   ### CSV Files ###
#   #################
#   # Will write csv files based off posterior for alpha, lambda prime, and psi
#   m <- matrix(c(codonList[-c(62,63,64)], alphaList), ncol = 2, byrow = FALSE)
#   colnames(m) <- c("Codon", "Alpha")
#   write.csv(m, file.path("UnitTestingOut", "RFPAlphaValues.csv"), quote = F, row.names = F)
#   
#   
#   m <- matrix(c(codonList[-c(62,63,64)], lambdaPrimeList), ncol = 2, byrow = FALSE)
#   colnames(m) <- c("Codon", "LambdaPrime")
#   write.table(m, file.path("UnitTestingOut", "RFPLambdaPrimeValues.csv"), quote = F, row.names = F)
#   
#   
#   m <- matrix(c(ids, psiList, psiList), ncol = 3, byrow = FALSE)
#   colnames(m) <- c("Gene", "PsiValue", "PsiValue")
#   write.table(m, file.path("UnitTestingOut", "RFPPsiValues.csv"), quote = F, row.names = F)
# })
