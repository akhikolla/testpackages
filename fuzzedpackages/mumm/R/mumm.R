#' Fit Multiplicative Mixed Models with TMB
#'
#' Fit a multiplicative mixed-effects model to data with use of the Template Model Builder.
#'
#' @param formula a two-sided formula object describing the linear fixed-effects and random-effects part
#' together with the multiplicative part. The response is on the left of a ~ operator and the terms which are
#' separated by + operators are on the right. The random-effect terms are recognized by vertical bars "|",
#' separating an expression for a model matrix and a grouping factor. The syntax for the multiplicative term
#' is 'mp("random effect","fixed effect")'.
#'
#' @param data a data frame containing the variables in the formula.
#'
#' @param cor logical. If FALSE the random effect in the multiplicative term is assumed to be independent of
#' the corresponding random main effect.
#'
#' @param start a numeric vector of starting values for the parameters in the model.
#'
#' @param control a list of control parameters passed on to the \code{nlminb} function used for the optimization.
#'
#' @details Fit a multiplicative mixed model via maximum likelihood with use of the Template Model Builder.
#' A multiplicative mixed model is here considered as a model with a linear mixed model part and one
#' multiplicative term. A multiplicative term is here defined as a product of a random effect and a fixed effect,
#' i.e. a term that models a part of the interaction as a random coefficient model based on linear regression
#' on a fixed main effect.
#'
#' @return An object of class \code{mumm}.
#'
#' @examples
#' set.seed(100)
#' sigma_e <- 1.5
#' sigma_a <- 0.8
#' sigma_b <- 0.5
#' sigma_d <- 0.7
#' nu <- c(8.2, 6.2, 2.3, 10.4, 7.5, 1.9)

#' nA <- 15
#' nP <- 6
#' nR <- 5

#' a <- rnorm(nA, mean = 0, sd = sigma_a)
#' b <- rnorm(nA, mean = 0, sd = sigma_b)
#' d <- rnorm(nA*nP, mean = 0, sd = sigma_d)
#' e <- rnorm(nA*nP*nR, mean = 0, sd = sigma_e)

#' Assessor <- factor(rep(seq(1,nA),each = (nP*nR)))
#' Product <- factor(rep(rep(seq(1,nP),each = nR), nA))
#' AssessorProduct <- (Assessor:Product)

#' y <- nu[Product] + a[Assessor] + b[Assessor]*(nu[Product]-mean(nu)) + d[AssessorProduct] + e

#' sim_data <- data.frame(y, Assessor, Product)

#' fit <- mumm(y ~ 1 + Product + (1|Assessor) + (1|Assessor:Product) +
#'              mp(Assessor,Product) ,data = sim_data)

#'
#' @useDynLib mumm
#' @importFrom lme4 subbars findbars mkReTrms nobars
#' @importFrom TMB MakeADFun sdreport tmbprofile
#' @importFrom Rcpp sourceCpp
#' @importFrom Matrix t
#' @importFrom methods as
#' @export
mumm <- function(formula, data, cor = TRUE, start = c(), control = list()) {

  #Checiking input:


  #-------------Building fixed effect design matrix, X-----------------------

  fixedform <- formula   #formula_obj
  fixedform[[3]]<-nobars(fixedform[[3]])
  terms_fix_mult = attr(stats::terms(fixedform),"term.labels")

  #Seperating the fixed effect terms and the multiplicative terms
  TFind = rep(FALSE,length(terms_fix_mult))

  for(i in 1:length(terms_fix_mult)) {
    TFind[i] = substr(terms_fix_mult[i],1,3) == "mp(";
  }
  NUMind = 1:length(terms_fix_mult);
  NUMind1 = NUMind[TFind]  #index for mp terms
  NUMind2 = NUMind[!TFind] #index for fixed terms

  terms_fix = stats::drop.terms(stats::terms(fixedform),NUMind1,keep.response = TRUE)
  terms_mult = stats::drop.terms(stats::terms(fixedform),NUMind2,keep.response = TRUE)
  formula_mult = formula(terms_mult)

  #The number of multiplicative terms in the model
  n_mult = length(attr(terms_mult,"term.labels"))

  #Finding the random and fixed effects that are part of the mulitiplicative terms
  randomef = vector(length = 0)
  fixedef = vector(length = 0)
  for(i in 1:n_mult) {
    mult_term = attr(terms_mult,"term.labels")[i];
    pos = stringr::str_locate_all(pattern = ",",mult_term)
    randomef[i] = substr(mult_term,4,pos[[1]][1]-1);
    fixedef[i] = substr(mult_term,pos[[1]][1]+1,nchar(mult_term)-1);

  }
  #Removing potential blank space from the strings
  randomef = stringr::str_trim(randomef)
  fixedef = stringr::str_trim(fixedef)


  #Seperating the fixed effects that is and isn't part of the multiplicative term
  remove = rep(FALSE,length(attr(terms_fix,"term.labels")))
  for(i in 1:length(fixedef)){

    if(sum((attr(terms_fix,"term.labels")==fixedef[i]))==0) {
      stop(sprintf("%s in multipliactive term is not a fixed effect",fixedef[i]))
      }

    remove = remove + (attr(terms_fix,"term.labels")==fixedef[i])
  }
  NUMremove = 1:length(attr(terms_fix,"term.labels"));
  NUMremove1 = NUMremove[as.logical(remove)];  #part of scaling
  NUMremove2 = NUMremove[!as.logical(remove)]; #not part of scaling


  #------------Building fixed effect design matrices, X and Xnu -------------------

  Xbig = stats::model.matrix(terms_fix,data = data)

  Xindex = as.logical(attr(Xbig,"assign")%in%NUMremove2 + (attr(Xbig,"assign")==0))
  X = Xbig[,Xindex, drop= F]
  Xnu = Xbig[,!Xindex, drop = F]

  t = as.data.frame(table(attr(Xbig,"assign")))

  #The number of parameters to be estimated for each fixed effect in the multiplicative term
  sizenu = t[t[["Var1"]]%in%NUMremove1,]$Freq



  #--------------- Building random effect design matrix, Z ----------------------------
  randform <- formula

  if (is.null(findbars(randform[[3]]))) {

    Z = matrix(0, nrow = nrow(data), ncol = 0)
    Z = as(Z,"dgTMatrix")
    rterms = NULL
    npar= 0

  } else {

    rterms = mkReTrms(findbars(randform[[3]]),data)
    Z = t(rterms$Zt)
    npar = sapply(rterms$flist,nlevels)
  }

  #number of random effects in the model (excluding random regression coef.)
  n_rand = length(rterms$cnms)

  #finding the index of the random effects that is related to the regression coef
  #(due to the correlation)
  TFindex = attr(npar,"names") == randomef
  NUMindex = which(TFindex)
  npar_cumsum = c(0,cumsum(npar))
  indexIna = c(npar_cumsum[NUMindex]+1,npar_cumsum[NUMindex+1])

  if (length(NUMindex) == 0){
    NUMindex = 0
    indexIna = c(1,0)
    cor = FALSE
  }

  data_named = data

  colnames(data_named)[colnames(data)==formula[[2]]] <- "y"

  #If the fixed effect in the mp term is an interaction, add it to data_named
  for(i in 1:length(fixedef)){
    if(regexpr(":",fixedef[i]) != -1){
      data_named[fixedef[i]] = with(data = data_named,eval(parse(text=fixedef[i])))
    }
  }

  nlevelsf = sapply(data_named[fixedef],nlevels);
  nlevelsr = sapply(data_named[randomef],nlevels);

  ffac = data.matrix(data_named[fixedef]);
  ffac = matrix(ffac,dimnames = NULL, ncol = length(fixedef))
  rfac = data.matrix(data_named[randomef]);
  rfac = matrix(rfac,dimnames = NULL, ncol = length(randomef))

  dataTMB = list(y = data_named[['y']],ffac = ffac, rfac = rfac,
                 X = X, Z = Z, Xnu = Xnu,  npar = npar, nlevelsf = nlevelsf, nlevelsr = nlevelsr, sizenu = sizenu,
                 indexIna = indexIna, indexInSiga = NUMindex)


  par_ix = c(ncol(dataTMB$X),ncol(Xnu),n_rand,n_mult,1,1)
  par_cix = cumsum(par_ix)


  if(length(start)==0){
    start = c(rep(0,ncol(dataTMB$X)+ncol(Xnu)),rep(1,n_rand+n_mult+1),0)
  }else if(length(start)<sum(par_ix)){
    start = c(start,0)
  }

  TF_beta = as.logical(par_cix[1])
  TF_sigma_a = (par_cix[2]+1)<=par_cix[3]

  parameters = list(
    beta  =  start[(0+as.logical(par_cix[1])):par_cix[1]],
    a = rep(0,ncol(dataTMB$Z)),
    b  = rep(0,sum(nlevelsr)),
    nu  = start[(par_cix[1]+1):par_cix[2]],
    log_sigma_a     = log(start[((par_cix[2]+1)*TF_sigma_a):(par_cix[3]*TF_sigma_a)]),
    log_sigma_b     = log(start[(par_cix[3]+1):par_cix[4]]),
    log_sigma       = log(start[length(start)-1]),
    transf_rho      = sign(start[length(start)])*sqrt(start[length(start)]^2/(1-start[length(start)]^2))
  )


  if(cor == TRUE){

    obj <- MakeADFun(
      data=dataTMB,
      parameters= parameters,
      random = c("a","b"),
      DLL    = "mumm",
      silent = TRUE
    )
  }else{

    obj <- MakeADFun(
      data=dataTMB,
      parameters= parameters,
      random = c("a","b"),
      DLL    = "mumm",
      silent = TRUE,
      map=list(transf_rho=factor(NA)) #fixerer korrelationen til 0
    )
  }



  opt = stats::nlminb(obj$par,obj$fn,obj$gr, control = control);

  sdr = sdreport(obj)


  #---------------------------- Finalizing output ---------------------------------------

  ## The estimated fixed effect coefficients
  par_fix = opt$par[1:(length(opt$par)-(n_rand+n_mult+1+cor))]


  names_fixed_ef = c(colnames(X),colnames(Xnu))

  names(par_fix) = names_fixed_ef

  est_cor = opt$par['transf_rho']/ sqrt(1. + opt$par['transf_rho']^2);
  names(est_cor) = "Correlation"

  ##The estimated variance components
  sigmas = exp(opt$par[(length(opt$par)-(n_rand+n_mult+(cor==TRUE))):(length(opt$par)-(cor==TRUE))])

  ##Creating a vector with names for the variance components
  names_random_ef = c(names(rterms$flist),paste0("mp ",randomef,":",fixedef),"Residual")
  names(sigmas) = names_random_ef

  ##Random effects
  par_rand = sdr$par.random
  ## Creating a vector with names for the random parameters
  mp_rdata = data_named[randomef]

  if(is.null(findbars(randform[[3]]))){
    names_random_par = sapply(mp_rdata,levels)
  } else {
    names_random_par = c(rterms$Zt@Dimnames[[1]],sapply(mp_rdata,levels))
  }

  names(par_rand) = names_random_par

  nlevels_par_rand = c(sapply(rterms$flist,nlevels),sapply(mp_rdata,nlevels))


  res = list(par = opt$par, objective = opt$objective, convergence = opt$convergence,
             iterations = opt$iterations, evaluations = opt$evaluations, convmessage = opt$message,
             par_fix = par_fix, sigmas = sigmas, est_cor = est_cor , par_rand = par_rand, nlevels_par_rand =  nlevels_par_rand,
             call = match.call(), nobs = nrow(data), df = length(opt$par), sdreport = sdr, obj = obj,
             index_num = NUMind1, data = data)

  class(res) <- "mumm"
  return(res)
}


