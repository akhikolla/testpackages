#' @title Accessor for 'nb' slot
#' @description This function outputs number of blocks of the variable block structure.
#' @name getNb
#' @param object Object of class "VB" or "HMMVB".
#' @export getNb
setGeneric("getNb", function(object) {standardGeneric("getNb")})

#' @title Accessor for 'dim' slot
#' @description This function outputs dimensionality.
#' @name getDim
#' @param object Object of class "VB", "HMM" or "HMMVB".
#' @export getDim
setGeneric("getDim", function(object) {standardGeneric("getDim")})


#' @title Accessor for 'bdim' slot
#' @description This function outputs dimensionality of blocks of variable block structure.
#' @name getBdim
#' @param object Object of class "VB" or "HMMVB".
#' @export getBdim
setGeneric("getBdim", function(object) {standardGeneric("getBdim")})

#' @title Accessor for 'numst' slot
#' @description This function outputs the number of states for each variable block in the variable block structure, 
#' the number of states of the HMM, or the number of states for each variable block of the HMM-VB.
#' @name getNumst
#' @param object Object of class "VB", "HMM" or "HMMVB".
#' @export getNumst
setGeneric("getNumst", function(object) {standardGeneric("getNumst")})

#' @title Accessor for 'varorder' slot
#' @description This function outputs the ordering of the variable blocks.
#' @name getVarorder
#' @param object Object of class "VB" or "HMMVB".
#' @export getVarorder
setGeneric("getVarorder", function(object) {standardGeneric("getVarorder")})

#' @title Accessor for 'prenumst' slot
#' @description This function outputs the number of states in the HMM for the preceding block of HMM-VB.
#' @name getPrenumst
#' @param object Object of class "HMM".
#' @export getPrenumst
setGeneric("getPrenumst", function(object) {standardGeneric("getPrenumst")})

#' @title Accessor for parameters of HMM
#' @description This function outputs a list with means, covariance matrices, inverse covarince matrices and
#' logarithms of the determinants of the covariance matrices for all states of the HMM.
#' @name getHmmParam
#' @param object Object of class "HMM".
#' @export getHmmParam
setGeneric("getHmmParam", function(object) {standardGeneric("getHmmParam")})

#' @title Accessor for 'VbStructure' slot.
#' @description This function outputs the variable block structure in the HMM-VB.
#' @name getVb
#' @param object Object of class "HMMVB".
#' @export getVb
setGeneric("getVb", function(object) {standardGeneric("getVb")})

#' @title Accessor for 'HmmChain' slot.
#' @description This function outputs a list with trained HMMs.
#' @name getHmmChain
#' @param object Object of class "HMMVB".
#' @export getHmmChain
setGeneric("getHmmChain", function(object) {standardGeneric("getHmmChain")})

#' @title Accessor for 'BIC' slot. 
#' @description This function outputs BIC for a trained HMM-VB model or a vector with BIC values calculated in model selection.
#' @name getBIC
#' @param object Object of class "HMMVB" or "HMMVBBIC".
#' @export getBIC
setGeneric("getBIC", function(object) {standardGeneric("getBIC")})

#' @title Accessor for 'Loglikehd' slot.
#' @description This function outputs Loglikelihood for each data point in a trained HMM-VB model or Loglikelihood for a new dataset in a HMM-VB model.
#' @name getLoglikehd
#' @param object Object of class "HMMVB", "HMMVBBIC" "HMMVBclust".
#' @export getLoglikehd
setGeneric("getLoglikehd", function(object) {standardGeneric("getLoglikehd")})

#' @title Accessor for 'diagCov' slot.
#' @description This function outputs diagCov logical indicator of diagonal covariance 
#' matrices for HMM-VB model.
#' @name getDiagCov
#' @param object Object of class "HMMVB".
#' @export getDiagCov
setGeneric("getDiagCov", function(object) {standardGeneric("getDiagCov")})

#' @title Accessor for 'clustParam' slot.
#' @description This function outputs clusterPar for the object of class HMMVBclust.
#' @name getClustParam
#' @param object Object of class "HMMVBclust".
#' @export getClustParam
setGeneric("getClustParam", function(object) {standardGeneric("getClustParam")})

#' @title Accessor for 'clsid' slot.
#' @description This function outputs the cluster labels for the object of class HMMVBclust.
#' @name getClsid
#' @param object Object of class "HMMVBclust".
#' @export getClsid
setGeneric("getClsid", function(object) {standardGeneric("getClsid")})

#' @title Accessor for 'size' slot.
#' @description This function outputs the number of points in each cluster for the object of class HMMVBclust.
#' @name getSize
#' @param object Object of class "HMMVBclust".
#' @export getSize
setGeneric("getSize", function(object) {standardGeneric("getSize")})

#' @title Accessor for 'optHMMVB' slot.
#' @description This function outputs the optimal HMM-VB found via BIC model selection.
#' @name getOptHMMVB
#' @param object Object of class "HMMVBBIC".
#' @export getOptHMMVB
setGeneric("getOptHMMVB", function(object) {standardGeneric("getOptHMMVB")})
