## #' Stepwise model selection in (graphical) interaction models
## #' 
## #' Stepwise model selection in (graphical) interaction models
## #' 
## #' 
## #' @aliases stepwise.iModel backward forward
## #' @param object An \code{iModel} model object
## #' @param criterion Either \code{"aic"} or \code{"test"} (for
## #'     significance test)
## #' @param alpha Critical value for deeming an edge to be significant/
## #'     insignificant. When \code{criterion="aic"}, \code{alpha}
## #'     defaults to 0; when \code{criterion="test"}, \code{alpha}
## #'     defaults to 0.05.
## #' @param type Type of models to search. Either \code{"decomposable"}
## #'     or \code{"unrestricted"}. If \code{type="decomposable"} and the
## #'     initial model is decompsable, then the search is among
## #'     decomposable models only.
## #' @param search Either \code{'all'} (greedy) or \code{'headlong'}
## #'     (search edges randomly; stop when an improvement has been
## #'     found).
## #' @param steps Maximum number of steps.
## #' @param k Penalty term when \code{criterion="aic"}. Only k=2 gives
## #'     genuine AIC.
## #' @param fixin Matrix (p x 2) of edges. If those edges are in the
## #'     model, they are not considered for removal.
## #' @param fixout Matrix (p x 2) of edges. If those edges are not in
## #'     the model, they are not considered for addition.
## #' @param direction Direction for model search. Either
## #'     \code{"backward"} or \code{"forward"}.
## #' @param details Controls the level of printing on the screen.
## #' @param trace For debugging only.
## #' @param \dots Further arguments to be passed on to \code{testdelete}
## #'     (for \code{testInEdges}) and \code{testadd} (for
## #'     \code{testOutEdges}).
## #' @return An \code{iModel} model object.
## #' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
## #' @seealso \code{\link{cmod}} \code{\link{dmod}} \code{\link{mmod}}
## #'     \code{\link{testInEdges}} \code{\link{testOutEdges}}
## #' @keywords models
## #' @examples
## #' 
## #' data(reinis)
## #' ## The saturated model
## #' m1 <- dmod(~.^., data=reinis)
## #' m2 <- stepwise(m1)
## #' m2
## #' 
## #' 














