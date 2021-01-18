#' @title Make predictions from a probabilistic network
#' 
#' @description Makes predictions (either as the most likely state or
#'     as the conditional distributions) of variables conditional on
#'     finding (evidence) on other variables in an independence
#'     network.
#'
#' @name grain_predict
#' 
#' @param object A grain object
#' @param response A vector of response variables to make predictions
#'     on
#' @param predictors A vector of predictor variables to make
#'     predictions from.  Defaults to all variables that are note
#'     responses.
#' @param newdata A data frame
#' @param type If "class", the most probable class is returned; if
#'     "distribution" the conditional distrubtion is returned.
#' @param ... Not used
#' @return A list with components \item{pred}{A list with the
#'     predictions} \item{pFinding}{A vector with the probability of
#'     the finding (evidence) on which the prediction is based}
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @seealso \code{\link{grain}}
#' @references Søren Højsgaard (2012). Graphical Independence Networks
#'     with the gRain Package for R. Journal of Statistical Software,
#'     46(10), 1-26.  \url{http://www.jstatsoft.org/v46/i10/}.
#' @keywords models
#'
#' @examples
#' data(chest_cpt)
#' data(chestSim500)
#'
#' chest.bn <- grain(compileCPT(chest_cpt))
#' nd <- chestSim500[1:4]
#'
#' predict(chest.bn, response="bronc", newdata=nd)
#' predict(chest.bn, response="bronc", newdata=nd, type="distribution")
#'

#' @export
predict.grain <- function(object, response, predictors=setdiff(names(newdata), response),
                          newdata, type="class", ...){

    ## cat("+++ predict.grain\n")
    if (!isCompiled(object)){
        object <- compile(object, propagate=TRUE)
    }

    ## if (!inherits(object, "compgrain"))
    ##     object <- compile(object, propagate=TRUE)

    type <- match.arg(type, c("class", "distribution"))
    nstate <- nodeStates(object, response)
    if (missing(predictors))
        predictors  <- setdiff(names(newdata),response)

    p.evec <- rep(NA, nrow(newdata))
    ans <- lapply(nstate,
                  function(a){
                      v <- matrix(NA, ncol=length(a),nrow=nrow(newdata))
                      colnames(v) <- a
                      v
                  })

    nd <- do.call(cbind, lapply(newdata, as.character))
    nd <- nd[, predictors, drop=FALSE]

    ## print(head(nd))
    vn <- colnames(nd)
    ## print(vn)

    #for (i in 1:nrow(newdata)){
        ## case      <- newdata[i, predictors, drop=FALSE]
        ## objecttmp1    <- setFinding(object, nodes=names(case), states=case)

    for (i in 1:nrow(nd)){
        objecttmp1    <- setFinding(object, nodes=vn, states=nd[i,,drop=FALSE])
        p.e       <- pEvidence(objecttmp1)
        ##cat(sprintf("pEvidence=%20.18f\n", p.e))
        if (p.e < .Machine$double.xmin){
            cat(sprintf("The evidence for row %i has probability smaller than %f in then model.\
                         Consider using the 'smooth' argument when building the network. \
                         Exiting...\n",
                        i, .Machine$double.xmin))
            return(NULL)
        }

        p.evec[i] <- p.e
        for (j in 1:length(response)){
            pj   <- .nodeMarginal(objecttmp1, response[j])[[1]] ## BRIS
                                        #print(pj)
            ans[[j]][i,] <- pj
        }
    }

    if (type=="class"){
        ns <- nodeStates(object, response)
        for (i in 1:length(ans)){
            a<-ans[[i]]
                                        #print(a)
            mlc <- apply(a,1,which.max)
                                        #print(mlc)

            ans[[i]] <- ns[[i]][mlc]
        }
    }

    value <- list(pred=ans, pEvidence=p.evec)
    return(value)
}
