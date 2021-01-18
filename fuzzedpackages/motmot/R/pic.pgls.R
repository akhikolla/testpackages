#' @title Fast PLGS estimation based on contrasts
#' @description Estimates regression parameters for a phylogenetic generalised least-squares analysis using the fast constrasts method (Felsenstein 1973; 1985; Freckleton 2012). This implementation is applicable for continuous traits only and not factors
#' @param formula A model formula with continuous variables only
#' @param phy An object of class \code{phylo} (see \pkg{ape}).
#' @param y A matrix of trait values with rownames corresponding to the phy tip labels, and column names corresponding to the formula variable names
#' @param lambda Default is "ML" meaning the phylogenetic signal of the response variable will be estimated using restricted Maximum likelihood. If a numeric value between 0-1 is provided this will be used in the calculation of regression coefficients
#' @param return.intercept.stat Logical. If \code{TRUE} the standard errors, t value, and p value of the estimated Intercept is provided for comparison with output from \code{\link[caper]{pgls}} from \pkg{caper} etc.,. Default is \code{FALSE} as this slows the function as it involves constructing and calculating the inverse of the phy variance-covariance matrix, and based on contrasts the design matrix column of ones would have zero contrasts.
#' @param meserr A vector (or matrix) of measurement error for each tip. This is only applicable to univariate analyses.
#' @return A list containing the model, model summary, intercept, estimate of Lambda, model log-Likelihood, model AIC
#' @useDynLib motmot
#' @importFrom Rcpp sourceCpp
#' @seealso \code{\link[caper]{pgls}}
#' @author Mark N Puttick and Rob Freckleton
#' @examples
#' # Data and phylogeny
#' data(anolis.tree)
#' anolis.tree$node.label <- NULL
#' lm.data <- transformPhylo.sim(phy=anolis.tree, n=2, model="bm")
#' dat <- data.frame(x = lm.data[,1], y = lm.data[,2],
#' names = anolis.tree$tip, row.names = anolis.tree$tip)
#' picModel <- pic.pgls(formula=y ~  x,
#' phy=anolis.tree, y = dat, lambda=1,
#' return.intercept.stat=FALSE)
#' @export

pic.pgls <- function(formula, phy, y, lambda="ML", return.intercept.stat=FALSE, meserr=NULL) {
	
	formula.break <- strsplit(as.character(formula), "~")
	
	response <- attr(terms.formula(formula, data=y), "term.labels")
	model.data <- model.frame(formula, y, na.action = na.pass)
	class.data <- !all(sapply(model.data, function(x) is.numeric(x)[1]))
	if(class.data) stop("please provide numerical variables only")
	
	
	yy <- strsplit(as.character(formula), "~")[[2]]
	response.xx <- paste(attr(terms.formula(formula, data=y), "term.labels"), collapse=" + ")
	if(response.xx == "") response.xx <- 1
	formula <- as.formula(paste(yy, " ~ 0 + ", response.xx))	
	
	response.var <- as.matrix(model.data[,1])
	explanatory.var <- as.matrix(model.data[,-1])
	rownames(explanatory.var) <- rownames(response.var) <- rownames(model.data)
	model.var <- cbind(response.var, explanatory.var)
	
	if(lambda == "ML") {	
		lambda.model <- transformPhylo.ML(y=response.var, phy=phy, model='lambda', returnPhy=TRUE, modelCIs=FALSE, meserr=meserr)
		phy.lambda <- lambda.model$lambdaPhy
		lambda <- lambda.model$Lambda[1]
	} else {
		lambda.model <- transformPhylo.ll(y=response.var, phy=phy, model='lambda', lambda=lambda, meserr=meserr)		
		phy.lambda <- transformPhylo(phy=phy, model='lambda', lambda=lambda)		
		}
	
	pic.data <- data.frame(apply(model.var, 2, function(x) pic(x, phy.lambda)))
	colnames(pic.data) <- colnames(model.data)
	model <- lm(formula, pic.data, y = TRUE, x=TRUE)
	model.y <- transformPhylo.ML(y=as.matrix(model.var[,1]), phy.lambda, model="bm", meserr=meserr)
	rst.y <- model.y$root.state
	rst.x <- 0
	if(dim(model.var)[2] > 1) for(u in 2:dim(model.var)[2]) rst.x[u-1] <- transformPhylo.ML(y=as.matrix(model.var[,u]), phy.lambda, model="bm", meserr=meserr)$root.state
	intercept <- rst.y - sum(rst.x *  model[[1]])
	vars <- pic(model.var[,1], phy.lambda, var.contrasts = TRUE)[,2]
	sigma2 <- model.y$brownianVariance
	V <- pic.motmot(x=model.var[,1], phy=phy.lambda)$V
	u <- resid(model)
	n <- length(u)
	k <- length(coef(model)) 
	logLik <- -0.5 * ((n + 1) * log( 2 * pi * sigma2) + sum(log(vars)) + sum(u ^ 2) / sigma2)
	ll.model <- logLik - 0.5 * log(V)
	
	
	output <- list()
	output$model <- model
	output$model.summary <- summary(model)
	output$intercept <- intercept
	output$lambda <- lambda
	output$logLikelihood <- ll.model
	output$AIC <- (-2 * ll.model) + k
	
	if(return.intercept.stat) {	
		n <- dim(y)[1]
		x.values <- matrix(c(rep(1, n), model.var[,-1]), ncol=dim(y)[2])
		errMat <- t(x.values) %*% solve(vcv(phy.lambda)) %*% x.values
		errMat <- solve(errMat) * summary(model)[[6]] ^ 2
		sterr <- diag(errMat)
		sterr <- sqrt(sterr)
		intercept.error <- sterr[1]
		t.value <- abs((0-intercept)/intercept.error)
		p.value <- pt(t.value, dim(pic.data)[1] - dim(pic.data)[2], lower.tail=FALSE) * 2
		output$intercept <- c("intercept"=intercept, "std.error"=intercept.error, "t.value"=t.value, "p.value"=p.value)
		
	}
	return(output)
	print(output$model.summary)
}