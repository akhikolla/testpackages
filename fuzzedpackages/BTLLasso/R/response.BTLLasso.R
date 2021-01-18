#' Create response object for BTLLasso
#' 
#' Create a response object for \code{BTLLasso} and \code{cv.BTLLasso}
#' 
#' 
#' @param response Vector containing results (binary or ordinal) of single paired
#' comparisons. Alternatively, also a  \code{\link[psychotools]{paircomp}} object as defined 
#' in the package \code{psychotools} could be used. In this case, none of the further 
#' arguments are needed. 
#' @param first.object Vector (character or factor, same length as response) indicating the first
#' object of the respective paired comparison from response.
#' @param second.object Vector (character or factor, same length as response) indicating the second
#' object of the respective paired comparison from response.
#' @param subject Vector (character, same length as response) indicating the subject that
#' generated the respective paired comparison from response.
#' @param with.order Boolean vector containing indicators for each paired comparison if an order effect was 
#' present. By default, an order effect is assumed for each comparison. This option is relevant whenever 
#' only some of the paired comparisons had an order effect and others did not, for example if some matches are
#' played on neutral ground. This option is only effective if either \code{order.effect = TRUE} or \code{object.order.effect = TRUE}.
#' @return Object of class \code{response.BTLLasso}
#' @author Gunther Schauberger\cr \email{gunther.schauberger@@tum.de}
#' @seealso \code{\link{BTLLasso}}, \code{\link{cv.BTLLasso}}
#' @references Schauberger, Gunther and Tutz, Gerhard (2019): BTLLasso - A Common Framework and Software 
#' Package for the Inclusion  and Selection of Covariates in Bradley-Terry Models, \emph{Journal of 
#' Statistical Software}, 88(9), 1-29, \url{https://doi.org/10.18637/jss.v088.i09}
#' 
#' Schauberger, Gunther and Tutz, Gerhard (2017): Subject-specific modelling 
#' of paired comparison data: A lasso-type penalty approach, \emph{Statistical Modelling},
#' 17(3), 223 - 243
#' 
#' Schauberger, Gunther, Groll Andreas and Tutz, Gerhard (2018): 
#' Analysis of the importance of on-field covariates in the German Bundesliga, 
#' \emph{Journal of Applied Statistics}, 45(9), 1561 - 1578
#' @examples
#'
#' \dontrun{
#' ##############################
#' ##### Example how response object for Bundesliga data Buli1516 was created
#' ##############################
#' 
#' data(BuliResponse)
#' 
#' Y.Buli <- response.BTLLasso(response = BuliResponse$Result, 
#'                             first.object = BuliResponse$TeamHome,
#'                             second.object = BuliResponse$TeamAway,
#'                             subject = BuliResponse$Matchday)
#' 
#' 
#' ##############################
#' ##### Example to create response object from paircomp object
#' ##############################
#' data("Topmodel2007", package = "psychotree")
#' 
#' Y.models <- response.BTLLasso(Topmodel2007$preference)
#' X.models <- scale(model.matrix(preference~., data = Topmodel2007)[,-1])
#' rownames(X.models) <- paste0("Subject",1:nrow(X.models))
#' colnames(X.models) <- c("Gender","Age","KnowShow","WatchShow","WatchFinal")
#' 
#' set.seed(5)
#' m.models <- cv.BTLLasso(Y = Y.models, X = X.models)
#' }
response.BTLLasso <- function(response, first.object = NULL, second.object = NULL, 
  subject = NULL, with.order = rep(TRUE, length(response))) {
  
  if(inherits(response, "paircomp")){
    response <- as.matrix(response)

    model_names <- str_split(colnames(response),pattern=":")
    
    model_names <- matrix(unlist(model_names),nrow=2)

    subject <- paste0("Subject",rep(1:nrow(response),ncol(response)))
    first.object <- rep(model_names[1,],each=nrow(response))
    second.object <- rep(model_names[2,],each=nrow(response))
  }

  withS <- FALSE
  if (!is.null(subject)) {
    withS <- TRUE
    if (!is.character(subject)) 
      stop("Argument subject has to be a character vector")
  }
  
  if (!withS) {
    subject <- 1:length(response)
  }
  
  ly <- length(response)
  lo1 <- length(first.object)
  lo2 <- length(second.object)
  ls <- length(subject)
  lorder <- length(with.order)
  
  if (!all(sapply(list(lo1, lo2, ls,lorder), identical, ly))) 
    stop("The arguments response, first.object, second.object and (if specified) subject and with.order
     have to be of the same length")

  all.objects <- as.factor(as.character(unlist(list(first.object, 
    second.object))))
  object.names <- levels(all.objects)
  
  first.object <- as.numeric(all.objects[1:ly])
  second.object <- as.numeric(all.objects[(ly + 1):(2 * ly)])
  
  m <- length(object.names)
  
  ## make response ordered
  response <- as.ordered(response)
  
  # number of response categories
  q <- length(levels(response)) - 1
  k <- q + 1
  
  
  ## everything about the subjects
  subject.names <- levels(as.factor(subject))
  n <- length(subject.names)

  
  RET <- list(response = response, first.object = first.object, 
    second.object = second.object, subject = subject, withS = withS, 
    subject.names = subject.names, object.names = object.names, 
    n = n, m = m, k = k, q = q, with.order = with.order)
  
  class(RET) <- "responseBTLLasso"
  
  RET
}

