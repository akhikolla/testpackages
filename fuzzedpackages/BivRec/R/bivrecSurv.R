#################### CREATE A BIVREC OBJECT ######################

#####
#' Create a Bivariate Alternating Recurrent Event Object
#'
#' @description
#' This function creates a bivariate recurrent survival object to be used as a response variable in a model formula.
#'
#' @importFrom stats na.omit
#'
#' @param id Numeric or character vector of subject's unique identifier (i).
#' @param episode Vector indicating the observation or episode (j) for a subject (i). This will determine order of events for each subject.
#' @param xij Vector with the lengths of time spent in event of Type I for individual i in episode j.
#' @param yij Vector with the lengths of time spent in event of Type II for individual i in episode j.
#' @param d1 Vector of censoring indicator corresponding to Type I gap times (xij): = 1 for uncensored, and = 0 for censored gap times.
#' @param d2 Vector of censoring indicator corresponding to Type II gap times (yij): = 1 for uncensored, and = 0 for censored gap times. Note that in the last episode (or pair) of gap times, yij is always censored, i.e., d2 = 0.
#'
#' @details
#' All the following data requirements must be met to avoid errors or warnings:
#' \enumerate{
#'  \item All variables should have the same length.
#'  \item The Type I (xij) and Type II (yij) gap times must be numeric and meet the requirements xij\eqn{>}0, yij\eqn{\ge}0 respectively for all observations.
#'  \item Only unique, positive, consecutive integers starting from 1, are allowed within each subject for episode (i.e., 1, 2, ...), but they should be allowed to be in any order in the raw data as long as they are unique and have no gap within each subject.
#'  \item Censoring indicators d1 and d2 must be numeric. Only the last episode is allowed to have zero value for either, and for the last episode, the only combinations allowed are (d1, d2) = (0, 0) or (1, 0).
#'  \item No missing values. If a subject has missing values in any variables, the whole subject's data will not be used.
#' }
#'
#' @return A bivrecSurv object ready to be used as the response for analysis using \verb{bivrecReg} or \verb{bivrecNP}.
#'
#' @rdname BivRec
#' @export
#' @examples
#' library(BivRec)
#' set.seed(28)
#' sim_data <- simBivRec(nsize=100, beta1=c(0.5,0.5), beta2=c(0,-0.5))
#' bivrecsurv_data <- with(sim_data, bivrecSurv(id, epi, xij, yij, d1, d2))
#' class(bivrecsurv_data)

bivrecSurv <- function(id, episode, xij, yij, d1, d2) {

  #Check if anything is missing
  if (missing(xij)) stop("Missing gap times for Type I event (xij).")
  if (missing(yij)) stop("Missing gap times for Type II event (yij).")
  if (missing(id)) stop("Missing subject identifiers (id).")
  if (missing(episode)) stop("Missing episodes for each subject (episode).")
  if (missing(d2)) stop("Missing censoring indicator for Type I event (d1).")
  if (missing(d1)) stop("Missing censoring indicator for Type II event (d2).")

  Xcind <- d1
  Ycind <- d2

  #Check all vectors have same length
  all_lengths <- c(length(id),length(episode),length(xij),length(yij),length(Ycind),length(Xcind))
  if (length(unique(all_lengths)) != 1) stop("Data not cleaned. One or more input vectors (id, episode, xij, yij, d1, d2) differs in length from the rest.")

  #Check xij > 0 and yij >=0 both numeric vectors
  if (!is.numeric(xij)) stop("Data not cleaned. Time arguments (xij and yij) must be numeric.")
  if (!is.numeric(yij)) stop("Data not cleaned. Time arguments (xij and yij) must be numeric.")
  if (any(xij <= 0)) stop("Data not cleaned. Time arguments for event Type I (xij) must be positive.")
  if (any(yij < 0)) stop("Data not cleaned. Time arguments for event Type II (yij) must be non-negative")

  #Check censoring indicators are made of only 0 or 1 values
  if (any(Xcind!=0 & Xcind!=1)) stop("Data not cleaned. Indicator vector for Type I gap times (d1) must be made of 0 or 1 values only.")
  if (any(Ycind!=0 & Ycind!=1)) stop("Data not cleaned. Indicator vector for Type II gap times (d2) must be made of 0 or 1 values only.")

  #ensure id's are numeric
  if (!is.numeric(id)) {
    if (is.character(id)) {id = as.numeric(as.factor(id))} else {
      if (is.factor(id)) {id = as.numeric((id))} else {
        stop("Data not cleaned. The id vector must be numeric, character or factor.")}
    }
  }

  id_ref = id
  inputdf <- data.frame(id=id, epi=episode, xij=xij, yij=yij, d1=Xcind, d2=Ycind)

  #Checks for each subject
  err_xind = err_yind = err_epi = NULL
  unique_id <- unique(inputdf$id)

  for (i in 1:length(unique_id)) {
    sub_id <- unique_id[i]
    temp_by_subject <- subset(inputdf, inputdf$id==sub_id)
    temp_by_subject <- temp_by_subject[order(temp_by_subject$epi),]
    sub_n <- nrow(temp_by_subject)

    #Check last (d1,d2) = (1,0) or (0,0)
    last_cx <- temp_by_subject$d1[sub_n]
    last_cy <- temp_by_subject$d2[sub_n]
    if (last_cy!=0) {err_yind <- c(err_yind, sub_id)}
    if (last_cx==0) {if (last_cy==1) {err_yind <- c(err_yind, sub_id)}

    if (sub_n > 1) {
      other_cx <- temp_by_subject$d1[-sub_n]
      other_cy <- temp_by_subject$d2[-sub_n]
      #check all d1, d2 besides last are 0 or 1
      if (sum(other_cx)!=(sub_n-1)) {err_xind <- c(err_xind, sub_id)}
      if (sum(other_cy)!=(sub_n-1)) {err_yind <- c(err_yind, sub_id)}
      ind_check1 <- unique(other_cx >= other_cy)
      if (length(ind_check1) > 1) {
        err_xind <- c(err_xind, sub_id)
        err_yind <- c(err_yind, sub_id)
      } else {
        if (ind_check1==FALSE) {
          err_xind <- c(err_xind, sub_id)
          err_yind <- c(err_yind, sub_id)
        }
      }
    }

    #Check episodes don't have gaps
    if (!is.numeric(episode) & !is.integer(episode)) {stop("Episode variable must be numeric.")}
    for (j in 1:sub_n){
      if (temp_by_subject$epi[j]!=j) {
        err_epi <- c(err_epi, sub_id)}
    }
  }

  error_subjects <- unique(c(err_xind, err_yind, err_epi))
  if (length(error_subjects>0)){
    errmsg <- paste(error_subjects, collapse = ", ")
    msg <- paste("Warning: Data not cleaned. Subjects with id", errmsg,
                 "removed because of gaps in episodes or incorrect values for d1, d2.",
                 sep=" ")
    print(msg)
    df4mdat <- inputdf[-which(inputdf$id %in% error_subjects), ]
  } else {df4mdat <- inputdf}
  }

  #calculate censoring time
  ci=id2=NULL
  j=1
  df4mdat$zij <- df4mdat$xij + df4mdat$yij
  for (i in unique(df4mdat$id)){
    tempi=df4mdat[df4mdat$id == i,]
    if (nrow(tempi) == 1){
      ci=c(ci,tempi$zij)
      id2=c(id2,j)
    } else {
      ci=c(ci,rep(sum(tempi$zij),nrow(tempi)))
      id2=c(id2,rep(j,nrow(tempi)))
    }
    j=j+1
  }

  df4mdat <- cbind(id=id2, df4mdat[-1], ci)

  result <- list()
  result$id_ref = id_ref
  result$error_ids <- error_subjects
  d2check <- unique(Ycind)

  if (length(d2check)==1) {
    if (d2check==0) {
      stop("Data not cleaned. All episodes provided are censored (all d2=0).")
    }
  }

  result$data4Lreg <- mdat(dat=df4mdat) #data for Lee regression
  result$data4Creg <- df4mdat #data for Chang regression (this is also the df that is used in bivrecPlot)
  #####ADD data for cdf and marginal of NP model
  df4np <- df4mdat
  colnames(df4np)=c("id", "epi", "vij", "wij", "d1", "d2", "x0ij", "ci")
  df4np=df4np[,c("id","vij","wij","d2","d1","epi","x0ij","ci")] #change order of columns
  forcdf1 <- np.dat(df4np, ai=1)
  forcdf2 <- np.dat(df4np, ai=2)
  marg1 <- formarginal(dat = df4np) #this is from the reformat code
  marg2 <- formarginal(dat = df4np)
  formarg1 <- np.dat(dat=marg1, ai=1)
  formarg2 <- np.dat(dat=marg2, ai=2)
  #two np objects that have data for cdf and marg depending on ai
  result$dat4np1 <- list(forcdf=forcdf1, formarg=formarg1,refdata = df4np) #for ai=1
  result$dat4np2 <- list(forcdf=forcdf2, formarg=formarg2,refdata = df4np) #for ai=2

  result$call <- match.call()

  class(result) <- "bivrecSurv"
  return(result)

}

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

