############ assess ############################

assess<-function(d1, d2, ...){
  UseMethod("assess")}

assess.ace<-function(d1, d2, B = NULL, n.assess = 20, relative = TRUE,...){

RU<-d1$utility	

if(is.null(B)){
B<-d1$B}

if(inherits(d2,"ace")){
    d2a<-d2$phase2.d}
if(inherits(d2,"pace")){
    d2a<-d2$d}
if(inherits(d2,"matrix")){
    d2a<-d2}

if(d1$deterministic){
U1<-RU(d = d1$phase2.d, B = B)
U2<-RU(d = d2a, B = B)
} else{
U1<-rep(0,n.assess)
U2<-rep(0,n.assess)
for(j in 1:n.assess){
U1[j]<-mean(RU(d = d1$phase2.d, B = B[1]))
U2[j]<-mean(RU(d = d2a, B = B[1]))
}
}

if(!relative){
oldU1<-U1	
oldU2<-U2	
U1<-oldU2
U2<-oldU1}

eff<-NULL
if(d1$deterministic){
if(d1$glm){
if(d1$criterion=="D"){
p <- ncol(model.matrix(object = d1$formula, data = data.frame(d1$phase2.d)))
eff <- 100 * exp((U1 - U2) / p)}
if(d1$criterion=="E"){
eff <- 100 * U1 / U2}
if(d1$criterion=="A"){
eff <- 100 * U2/U1}
}
if(d1$nlm){
if(d1$criterion=="D"){
p <- length(setdiff(all.vars(d1$formula), dimnames(d1$phase2.d)[[2]]))		
eff <- 100 * exp((U1 - U2) / p)}
if(d1$criterion=="E"){
eff <- 100 * U1 / U2}
if(d1$criterion=="A"){
eff <- 100 * U2/U1}
}
}

if(!relative){
U1<-oldU1
U2<-oldU2}

out<-list(U1 = U1, U2 = U2, eff = eff, d1 = d1, d2 = d2)
class(out)<-"assess"

out}

assess.pace<-function(d1, d2, B = NULL, n.assess = 20, relative = TRUE,...){

RU<-d1$utility	

if(is.null(B)){
B<-d1$B}

if(inherits(d2,"ace")){
    d2a<-d2$phase2.d}
if(inherits(d2,"pace")){
    d2a<-d2$d}
if(inherits(d2,"matrix")){
    d2a<-d2}

if(d1$deterministic){
U1<-RU(d = d1$d, B = B)
U2<-RU(d = d2a, B = B)
} else{
U1<-rep(0,n.assess)
U2<-rep(0,n.assess)
for(j in 1:n.assess){
U1[j]<-mean(RU(d = d1$d, B = B[1]))
U2[j]<-mean(RU(d = d2a, B = B[1]))
}
}

if(!relative){
oldU1<-U1	
oldU2<-U2	
U1<-oldU2
U2<-oldU1}

eff<-NULL
if(d1$deterministic){
if(d1$glm){
if(d1$criterion=="D"){
p <- ncol(model.matrix(object = d1$formula, data = data.frame(d1$d)))
eff <- 100 * exp((U1 - U2) / p)}
if(d1$criterion=="E"){
eff <- 100 * U1 / U2}
if(d1$criterion=="A"){
eff <- 100 * U2/U1}
}
if(d1$nlm){
if(d1$criterion=="D"){
p <- length(setdiff(all.vars(d1$formula), dimnames(d1$d)[[2]]))		
eff <- 100 * exp((U1 - U2) / p)}
if(d1$criterion=="E"){
eff <- 100 * U1 / U2}
if(d1$criterion=="A"){
eff <- 100 * U2/U1}
}
}

if(!relative){
U1<-oldU1
U2<-oldU2}

out<-list(U1 = U1, U2 = U2, eff = eff, d1 = d1, d2 = d2)
class(out)<-"assess"

out}

print.assess<-function(x, ...){

if(x$d1$deterministic){
cat("Approximate expected utility of d1 =", x$U1, "\n")
cat("Approximate expected utility of d2 =", x$U2, "\n")
if(x$d1$glm|x$d1$nlm){
if(x$d1$criterion=="D"|x$d1$criterion=="E"|x$d1$criterion=="A"){
cat("Approximate relative ",x$d1$criterion,"-efficiency = ", x$eff,"% \n",sep="")
}
}
} else{
cat("Mean (sd) approximate expected utility of d1 = ", mean(x$U1)," (",sd(x$U1) ,") \n",sep="")
cat("Mean (sd) approximate expected utility of d2 = ", mean(x$U2)," (",sd(x$U2) ,") \n",sep="")	
}	
	
}

summary.assess<-function(object, ...){

print.assess(x=object)
}

plot.assess<-function(x, ...){

if(x$d1$deterministic){
warning("No meaningful plot produced when deterministic = TRUE")
} else{
boxplot(x$U1,x$U2, names = c("d1","d2"), xlab= "Object", ylab="Approximate expected utility")}	
	
}	

############ efficiency ############################

# efficiency<-function(d1, d2, ...){
#   UseMethod("efficiency")}
# 
# efficiency.matrix<-function(d1, d2, relative=TRUE, model = "nlm", B, formula, family = NULL, prior, criterion, ...){
#   
#   if(model=="glm" & is.null(family)){
#     stop("The family argument must be specified when model argument is glm")}
#   
#   if(!any(c(criterion=="A"|criterion=="D"|criterion=="E"))){
#     stop("The criterion argument must be one of A, D or E")}
#     
#   if(model=="glm"){
#   RU <- utilityglm(formula = formula, family=family,prior = prior, criterion = criterion, method = "quadrature", nrq = B)} 
#   if(model=="nlm"){
#   RU <- utilitynlm(formula = formula, prior = prior, desvars = dimnames(d1)[[2]], criterion = criterion, method = "quadrature", nrq = B)} 
#   
#   U<-function(d) RU$utility(d = d, B = B)
# 
#   U1 <- mean(U(d1))
#   if(model=="glm"){
#   p <- ncol(model.matrix(object = formula, data = data.frame(d1)))}
#   if(model=="nlm"){
#   allvars <- all.vars(formula)
#   paravars <- setdiff(allvars, dimnames(d1)[[2]])
#   p <- length(paravars)}
# 
#   U2 <- mean(U(d2))
#   
#   if(relative){
#     out<-switch(EXPR=criterion,
#                 D = 100 * exp((U1 - U2) / p),
#                 E = 100 * U1 / U2,
#                 100 * U2/U1)} else{
#                   out<-switch(EXPR=criterion,
#                               D = 100 * exp((U2 - U1) / p),
#                               E = 100 * U2 / U1,
#                               100 * U1/U2)}
#   out}
# 
# efficiency.ace<-function(d1, d2, relative=TRUE, B, formula = NULL, family = NULL, prior = NULL, criterion = NULL, ...){
#   
# if(!(d1$nlm|d1$glm)){
# stop("efficiency function should only be used when either a) d1 is a result of a call to (p)aceglm or (p)acenlm with the argument criterion being A, D or E; or b) d1 is a matrix object and the call to efficiency has specified arguments formula, prior and criterion")}
# 
# if(d1$glm){
#  
#   if(is.null(formula)){
#     formula<-d1$formula}
#   if(is.null(family)){
#     family<-d1$family}
#   if(is.null(prior)){
#     prior<-d1$prior}
#   if(is.null(criterion)){
#     criterion<-d1$criterion}
#     
#   if(!any(c(criterion=="A"|criterion=="D"|criterion=="E"))){
#     stop("The criterion argument must be one of A, D or E")}
#   
#   RU <- utilityglm(formula = formula, family=family,prior = prior, criterion = criterion, method = "quadrature", nrq = B)
#   U<-function(d) RU$utility(d = d, B = B)
#   
#   U1 <- mean(U(d1$phase2.d))
#   p <- ncol(model.matrix(object = formula, data = data.frame(d1$phase2.d)))
#   if(inherits(d2,"ace")){
#     d2a<-d2$phase2.d}
#   if(inherits(d2,"pace")){
#     d2a<-d2$d}
#   if(inherits(d2,"matrix")){
#     d2a<-d2}
#   U2 <- mean(U(d2a))
#   
#   if(relative){
#     out<-switch(EXPR=criterion,
#                 D = 100 * exp((U1 - U2) / p),
#                 E = 100 * U1 / U2,
#                 100 * U2/U1)} else{
#                   out<-switch(EXPR=criterion,
#                               D = 100 * exp((U2 - U1) / p),
#                               E = 100 * U2 / U1,
#                               100 * U1/U2)}
#   out}
#   
#   
# if(d1$nlm){
# 
#   if(is.null(formula)){
#     formula<-d1$formula}
#   if(is.null(prior)){
#     prior<-d1$prior}
#   if(is.null(criterion)){
#     criterion<-d1$criterion}
#     
#   if(!any(c(criterion=="A"|criterion=="D"|criterion=="E"))){
#     stop("The criterion argument must be one of A, D or E")}
#   
#   RU <- utilitynlm(formula = formula, prior = prior, desvars = dimnames(d1$phase2.d)[[2]], 
#                          criterion = criterion, method = "quadrature", nrq = B)
#     U<-function(d) RU$utility(d = d, B = B)
# 
#     U1 <- mean(U(d1$phase2.d))
#     p <- length(setdiff(all.vars(formula), dimnames(d1$phase2.d)[[2]]))
#     if(inherits(d2,"ace")){
#     d2a<-d2$phase2.d}
#     if(inherits(d2,"pace")){
#     d2a<-d2$d}
#     if(inherits(d2,"matrix")){
#     d2a<-d2}
#     U2 <- mean(U(d2a))
#   
#   if(relative){
#     out<-switch(EXPR=criterion,
#          D = 100 * exp((U1 - U2) / p),
#          E = 100 * U1 / U2,
#          100 * U2/U1)} else{
#            out<-switch(EXPR=criterion,
#                        D = 100 * exp((U2 - U1) / p),
#                        E = 100 * U2 / U1,
#                        100 * U1/U2)}
#     out}
# 
# out}
# 
# efficiency.pace<-function(d1, d2, relative=TRUE, B, formula = NULL, family = NULL, prior = NULL, criterion = NULL, ...){
#   
# if(!(d1$nlm|d1$glm)){
# stop("efficiency function should only be used when either a) d1 is a result of a call to (p)aceglm or (p)acenlm with the argument criterion being A, D or E; or b) d1 is a matrix object and the call to efficiency has specified arguments formula, prior and criterion")}
# 
# if(d1$glm){
#  
#   if(is.null(formula)){
#     formula<-d1$formula}
#   if(is.null(family)){
#     family<-d1$family}
#   if(is.null(prior)){
#     prior<-d1$prior}
#   if(is.null(criterion)){
#     criterion<-d1$criterion}
#     
#   if(!any(c(criterion=="A"|criterion=="D"|criterion=="E"))){
#     stop("The criterion argument must be one of A, D or E")}
#   
#   RU <- utilityglm(formula = formula, family=family,prior = prior, criterion = criterion, method = "quadrature", nrq = B)
#   U<-function(d) RU$utility(d = d, B = B)
#   
#   U1 <- mean(U(d1$d))
#   p <- ncol(model.matrix(object = formula, data = data.frame(d1$d)))
#   if(inherits(d2,"ace")){
#     d2a<-d2$phase2.d}
#   if(inherits(d2,"pace")){
#     d2a<-d2$d}
#   if(inherits(d2,"matrix")){
#     d2a<-d2}
#   U2 <- mean(U(d2a))
#   
#   if(relative){
#     out<-switch(EXPR=criterion,
#                 D = 100 * exp((U1 - U2) / p),
#                 E = 100 * U1 / U2,
#                 100 * U2/U1)} else{
#                   out<-switch(EXPR=criterion,
#                               D = 100 * exp((U2 - U1) / p),
#                               E = 100 * U2 / U1,
#                               100 * U1/U2)}
#   out}
#   
#   
# if(d1$nlm){
# 
#   if(is.null(formula)){
#     formula<-d1$formula}
#   if(is.null(prior)){
#     prior<-d1$prior}
#   if(is.null(criterion)){
#     criterion<-d1$criterion}
#     
#   if(!any(c(criterion=="A"|criterion=="D"|criterion=="E"))){
#     stop("The criterion argument must be one of A, D or E")}
#   
#   RU <- utilitynlm(formula = formula, prior = prior, desvars = dimnames(d1$d)[[2]], 
#                          criterion = criterion, method = "quadrature", nrq = B)
#     U<-function(d) RU$utility(d = d, B = B)
# 
#     U1 <- mean(U(d1$d))
#     p <- length(setdiff(all.vars(formula), dimnames(d1$d)[[2]]))
#     if(inherits(d2,"ace")){
#     d2a<-d2$phase2.d}
#     if(inherits(d2,"pace")){
#     d2a<-d2$d}
#     if(inherits(d2,"matrix")){
#     d2a<-d2}
#     U2 <- mean(U(d2a))
#   
#   if(relative){
#     out<-switch(EXPR=criterion,
#          D = 100 * exp((U1 - U2) / p),
#          E = 100 * U1 / U2,
#          100 * U2/U1)} else{
#            out<-switch(EXPR=criterion,
#                        D = 100 * exp((U2 - U1) / p),
#                        E = 100 * U2 / U1,
#                        100 * U1/U2)}
#     out}
# 
# out}

############ utilityglm ############################

############ utilityglm ############################

utilityglm <- function (formula, family, prior, criterion = c("D", "A", "E", "SIG", "NSEL", "SIG-Norm", "NSEL-Norm"), 
                        method = c("quadrature", "MC"), nrq) 
{
  
  criterion <- match.arg(criterion)
  if(length(method) > 1) {
    method = switch(EXPR=criterion,
                    D = "quadrature",
                    A = "quadrature",
                    E = "quadrature",
                    "MC")
  } 
  method <- match.arg(method)
  if(identical(method, "MC") && !is.function(prior)) stop("For method = \"MC\", argument prior must specify a function.") 
  
  if(criterion %in% c("SIG", "NSEL", "SIG-Norm", "NSEL-Norm") && pmatch(method, "quadrature", nomatch = 0)) {
    warning("method = \"quadrature\" is not available for SIG and NSEL utilities. Using method = \"MC\"", immediate. = TRUE)
    method <- "MC"
  }
  
  if(missing(nrq)) {
    nrq <- switch(method,
                  quadrature = c(2, 8),
                  NULL)
  }  
  nr <- nrq[[1]]
  nq <- nrq[[2]]
  
  
  if (criterion == "A" | criterion == "D" | criterion == "E") {
    if(identical(method, "quadrature")) {
      #no.terms <- 1 + length(attr(terms(formula), "term.labels"))    #### Problem here if you don't want an interecpt
      no.terms <- attr(terms(formula), "intercept") + length(attr(terms(formula), "term.labels")) 
      if(identical(names(prior)[1:2], c("mu", "sigma2"))) {
        if(identical(length(prior$mu), as.integer(1))) {
          qmu <- rep(prior$mu, no.terms)
        } else {
          qmu <- prior$mu 
        }
        if(identical(dim(prior$sigma2), as.integer(2))) {
          qsigma <- prior$sigma
        } else {
          qsigma2 <- diag(prior$sigma2, nrow = no.terms)
        }
        quad <- RSquadrature(no.terms, qmu, qsigma2, nr, nq)
        abscissas <- as.matrix(quad$a)
        weights <- quad$w
      } else if(identical(names(prior)[1], "support")) {
        if(!identical(dim(t(prior$support)), as.integer(c(no.terms, 2)))) {
          stop("Uniform prior: the prior support must be specified as a 
               matrix with dimension c(2, number of terms); see help file ")
        }
        quad <- RSquadrature.uniform(no.terms, t(prior$support), nr, nq)
        abscissas <- as.matrix(quad$a)
        weights <- quad$w
        } else {
          stop("Argument \"prior\" must correctly specify a normal or uniform prior for the model parameters (see help file).")
        }
      ## dcw 28-11-2018 check the weights and abscissas are all OK
      if(anyNA(weights) || any(is.infinite(weights)) || any(is.nan(weights))) stop("The quadrature scheme is not valid for large values of nr and/or nq. Try making these values smaller.")
      if(anyNA(abscissas) || any(is.infinite(abscissas)) || any(is.nan(abscissas))) stop("The quadrature scheme is not valid for large values of nr and/or nq. Try making these values smaller.")
      inte <- function(d, B) {
        x <- model.matrix(object = formula, data = data.frame(d))
        v <- utilglm(x = x, beta = abscissas, family = family, criterion = criterion)
        weighted.mean(v, weights)
      }
    } else {
      inte <- function(d, B) {
        x <- model.matrix(object = formula, data = data.frame(d))
        beta <- prior(B)
        utilglm(x = x, beta = beta, family = family, criterion = criterion)
      }
    }
  }
  else {
    if (criterion == "SIG") {
      if (!is.function(family)) {
        if (is.list(family)) {
          stuff <- family
        }
        else {
          family2 <- get(family, mode = "function", envir = parent.frame())
          stuff <- family2()
        }
      }
      else {
        stuff <- family()
      }
      if (stuff$family != "binomial" & stuff$family != 
          "poisson") {
        stop("Family or link not implemented for Shannon information gain (SIG) utility yet")
      }
      else {
        if (stuff$family == "binomial") {
          if (stuff$link == "logit") {
            inte <- function(d, B) {
              sam <- prior(2 * B)
              x <- model.matrix(object = formula, data = data.frame(d))
              n1 <- dim(x)[1]
              rho <- 1/(1 + exp(-sam %*% t(x)))
              Z <- log(1 - rho[-(1:B), ])
              frho <- as.vector(.Call("rowSumscpp", Z, 
                                      PACKAGE = "acebayes"))
              y <- matrix(rbinom(n = n1 * B, size = 1, 
                                 prob = as.vector(rho[1:B, ])), ncol = n1)
              Z <- dbinom(x = y, size = 1, prob = rho[1:B, 
                                                      ], log = TRUE)
              rsll <- as.vector(.Call("rowSumscpp", Z, 
                                      PACKAGE = "acebayes"))
              sam <- sam[-(1:B), ]
              rsll4 <- as.vector(.Call("siglrcpp", y, 
                                       x, sam, frho, PACKAGE = "acebayes"))
              MY3 <- log(rsll4/B)
              eval <- rsll - MY3
              eval
            }
          }
          else {
            stop("Family or link not implemented for Shannon information gain (SIG) utility yet")
          }
        }
        if (stuff$family == "poisson") {
          if (stuff$link == "log") {
            inte <- function(d, B) {
              sam <- prior(2 * B)
              x <- model.matrix(object = formula, data = data.frame(d))
              n1 <- dim(x)[1]
              rho <- exp(sam %*% t(x))
              frho <- as.vector(.Call("rowSumscpp", rho[-(1:B), 
                                                        ], PACKAGE = "acebayes"))
              y <- matrix(rpois(n = n1 * B, lambda = as.vector(rho[1:B, 
                                                                   ])), ncol = n1)
              yfrac <- as.vector(.Call("rowSumscpp", 
                                       lfactorial(y), PACKAGE = "acebayes"))
              Z <- dpois(x = y, lambda = rho[1:B, ], 
                         log = TRUE)
              rsll <- as.vector(.Call("rowSumscpp", Z, 
                                      PACKAGE = "acebayes"))
              sam <- sam[-(1:B), ]
              rsll4 <- as.vector(.Call("sigprcpp", y, 
                                       x, sam, frho, yfrac, PACKAGE = "acebayes"))
              rsll - rsll4
            }
          }
          else {
            stop("Family or link not implemented for Shannon information gain (SIG) utility yet")
          }
        }
      }
    }
    if (criterion == "SIG-Norm") {
      if (!is.function(family)) {
        if (is.list(family)) {
          stuff <- family
        }
        else {
          family2 <- get(family, mode = "function", envir = parent.frame())
          stuff <- family2()
        }
      }
      else {
        stuff <- family()
      }
      if (stuff$family != "binomial" & stuff$family != 
          "poisson") {
        stop("Family or link not implemented for Shannon information gain (SIG) utility yet")
      }
      else {
        if (stuff$family == "binomial") {
          if (stuff$link == "logit") {
            inte <- function(d, B) {
              sam <- prior(2 * B)
              pm <- apply(sam[1:B, ], 2, mean)
              pv <- apply(sam[1:B, ], 2, var)
              sam <- sam[-(1:B), ]
              X <- model.matrix(object = formula, data = data.frame(d))
              n1 <- dim(X)[1]
              rho <- 1/(1 + exp(-sam %*% t(X)))
              y <- matrix(rbinom(n = n1 * B, size = 1, 
                                 prob = as.vector(rho)), ncol = n1)
              out <- as.vector(.Call("LRLAP2cpp", y, 
                                     sam, X, pm, pv, PACKAGE = "acebayes"))
              out
            }
          }
          else {
            stop("Family or link not implemented for Shannon information gain (SIG) utility yet")
          }
        }
        if (stuff$family == "poisson") {
          if (stuff$link == "log") {
            inte <- function(d, B) {
              sam <- prior(2 * B)
              pm <- apply(sam[1:B, ], 2, mean)
              pv <- apply(sam[1:B, ], 2, var)
              sam <- sam[-(1:B), ]
              X <- model.matrix(object = formula, data = data.frame(d))
              n1 <- dim(X)[1]
              rho <- exp(sam %*% t(X))
              y <- matrix(rpois(n = n1 * B, lambda = as.vector(rho)), 
                          ncol = n1)
              out <- as.vector(.Call("PRLAP2cpp", y, 
                                     sam, X, pm, pv, PACKAGE = "acebayes"))
              out
            }
          }
          else {
            stop("Family or link not implemented for Shannon information gain (SIG) utility yet")
          }
        }
      }
    }
    if (criterion == "NSEL") {
      if (!is.function(family)) {
        if (is.list(family)) {
          stuff <- family
        }
        else {
          family2 <- get(family, mode = "function", envir = parent.frame())
          stuff <- family2()
        }
      }
      else {
        stuff <- family()
      }
      if (stuff$family != "binomial" & stuff$family != 
          "poisson") {
        stop("Family or link not implemented for negative squared error loss (NSEL) utility yet")
      }
      else {
        if (stuff$family == "binomial") {
          if (stuff$link == "logit") {
            inte <- function(d, B) {
              sam <- prior(2 * B)
              x <- model.matrix(object = formula, data = data.frame(d))
              n1 <- dim(x)[1]
              rho <- 1/(1 + exp(-sam %*% t(x)))
              Z <- log(1 - rho[-(1:B), ])
              frho <- as.vector(.Call("rowSumscpp", Z, 
                                      PACKAGE = "acebayes"))
              y <- matrix(rbinom(n = n1 * B, size = 1, 
                                 prob = as.vector(rho[1:B, ])), ncol = n1)
              fsam <- sam[1:B, ]
              sam <- sam[-(1:B), ]
              rsll4 <- .Call("nsellrcpp", y, x, sam, 
                             frho, PACKAGE = "acebayes")
              Z <- (rsll4 - fsam)^2
              eval <- as.vector(.Call("rowSumscpp", Z, 
                                      PACKAGE = "acebayes"))
              -eval
            }
          }
          else {
            stop("Family or link not implemented for negative squared error loss (NSEL) utility yet")
          }
        }
        if (stuff$family == "poisson") {
          if (stuff$link == "log") {
            inte <- function(d, B) {
              sam <- prior(2 * B)
              x <- model.matrix(object = formula, data = data.frame(d))
              n1 <- dim(x)[1]
              rho <- exp(sam %*% t(x))
              frho <- as.vector(.Call("rowSumscpp", rho, 
                                      PACKAGE = "acebayes"))
              y <- matrix(rpois(n = n1 * B, lambda = as.vector(rho[1:B, 
                                                                   ])), ncol = n1)
              yfrac <- as.vector(.Call("rowSumscpp", 
                                       lfactorial(y), PACKAGE = "acebayes"))
              fsam <- sam[1:B, ]
              sam <- sam[-(1:B), ]
              rsll4 <- .Call("nselprcpp", y, x, sam, 
                             frho, yfrac, PACKAGE = "acebayes")
              Z <- (rsll4 - fsam)^2
              eval <- as.vector(.Call("rowSumscpp", Z, 
                                      PACKAGE = "acebayes"))
              -eval
            }
          }
          else {
            stop("Family or link not implemented for negative squared error loss (NSEL) utility yet")
          }
        }
      }
    }
    if (criterion == "NSEL-Norm") {
      if (!is.function(family)) {
        if (is.list(family)) {
          stuff <- family
        }
        else {
          family2 <- get(family, mode = "function", envir = parent.frame())
          stuff <- family2()
        }
      }
      else {
        stuff <- family()
      }
      if (stuff$family != "binomial" & stuff$family != 
          "poisson") {
        stop("Family or link not implemented for negative squared error loss (NSEL) utility yet")
      }
      else {
        if (stuff$family == "binomial") {
          if (stuff$link == "logit") {
            inte <- function(d, B) {
              sam <- prior(2 * B)
              pm <- apply(sam[1:B, ], 2, mean)
              pv <- apply(sam[1:B, ], 2, var)
              sam <- sam[-(1:B), ]
              X <- model.matrix(object = formula, data = data.frame(d))
              n1 <- dim(X)[1]
              rho <- 1/(1 + exp(-sam %*% t(X)))
              y <- matrix(rbinom(n = n1 * B, size = 1, 
                                 prob = as.vector(rho)), ncol = n1)
              out <- as.vector(.Call("LRNSELLAP2cpp", 
                                     y, sam, X, pm, pv, PACKAGE = "acebayes"))
              out
            }
          }
          else {
            stop("Family or link not implemented for negative squared error loss (NSEL) utility yet")
          }
        }
        if (stuff$family == "poisson") {
          if (stuff$link == "log") {
            inte <- function(d, B) {
              sam <- prior(2 * B)
              pm <- apply(sam[1:B, ], 2, mean)
              pv <- apply(sam[1:B, ], 2, var)
              sam <- sam[-(1:B), ]
              X <- model.matrix(object = formula, data = data.frame(d))
              n1 <- dim(X)[1]
              rho <- exp(sam %*% t(X))
              y <- matrix(rpois(n = n1 * B, lambda = as.vector(rho)), 
                          ncol = n1)
              out <- as.vector(.Call("PRNSELLAP2cpp", 
                                     y, sam, X, pm, pv, PACKAGE = "acebayes"))
              out
            }
          }
          else {
            stop("Family or link not implemented for negative squared error loss (NSEL) utility yet")
          }
        }
      }
    }
  }
  output <- list(utility = inte)
  output
}

############ utilitynlm ############################

############ utilitynlm ############################

utilitynlm <- function (formula, prior, desvars, criterion = c("D", "A", "E", "SIG", "NSEL"), 
                        method = c("quadrature", "MC"), nrq) 
{
  
  criterion <- match.arg(criterion)
  if(length(method) > 1) {
    method = switch(EXPR=criterion,
                    D = "quadrature",
                    A = "quadrature",
                    E = "quadrature",
                    "MC")
  } 
  method <- match.arg(method)
  if(identical(method, "MC") && !is.function(prior)) stop("For method = \"MC\", argument prior must specify a function.") 
  
  if(criterion %in% c("SIG", "NSEL") && pmatch(method, "quadrature", nomatch = 0)) {
    warning("method = \"quadrature\" is not available for SIG and NSEL utilities. Using method = \"MC\"", immediate. = TRUE)
    method <- "MC"
  }
  
  if(missing(nrq)) {
    nrq <- switch(method,
                  quadrature = c(2, 8),
                  NULL)
  }  
  nr <- nrq[[1]]
  nq <- nrq[[2]]
  
  
  allvars <- all.vars(formula)
  paravars <- setdiff(allvars, desvars)
  p <- length(paravars)
  k <- length(desvars)
  DD <- deriv(expr = formula, namevec = paravars)
  aDD <- as.character(DD)
  grad <- function() {
  }
  gradtext <- "grad<-function(d,paras){ \n"
  gradtext <- paste(gradtext, desvars[1], "<-d[,1]", " \n ", 
                    sep = "")
  if (k > 1) {
    for (j in 2:k) {
      gradtext <- paste(gradtext, desvars[j], "<-d[,", 
                        j, "]", " \n ", sep = "")
    }
  }
  gradtext <- paste(gradtext, paravars[1], "<-paras[,1]", " \n ", 
                    sep = "")
  if (p > 1) {
    for (j in 2:p) {
      gradtext <- paste(gradtext, paravars[j], "<-paras[,", 
                        j, "]", " \n ", sep = "")
    }
  }
  gradtext <- paste(gradtext, substr(x = aDD, start = 2, stop = nchar(aDD)), 
                    sep = "")
  eval(parse(text = gradtext))
  
  if (criterion == "A" | criterion == "D" | criterion == "E") {
    if(identical(method, "quadrature")) {
      no.terms <- p
      
      if(identical(names(prior)[1:2], c("mu", "sigma2"))) {
        if(is.null(names(prior$mu))) stop("The names attribute for mu must correspond to the named parameters in the model formula.")
        if(!identical(length(prior$mu), p)) {
          stop("Normal prior: mean must have the same length as the number of model parameters")
        } else {
          qmu <- prior$mu 
        }
        if(identical(dim(prior$sigma2), 2)) {
          qsigma <- prior$sigma
        } else {
          qsigma2 <- diag(prior$sigma2, nrow = no.terms)
        }
        if(!identical(TRUE, compare(paravars, names(prior$mu), ignoreOrder = TRUE)$result)) {
          stop("Normal prior: the names attribute for mu must correspond to the named parameters in the model formula.")
        }
        quad <- RSquadrature(no.terms, qmu, qsigma2, nr, nq)
        abscissas <- as.matrix(quad$a)
        colnames(abscissas) <- names(prior$mu)
        weights <- quad$w
      } else if(identical(names(prior)[1], "support")) {
        if(!identical(dim(t(prior$support)), as.integer(c(no.terms, 2)))) {
          stop("Uniform prior: the prior support must be specified as a 
               matrix with dimension c(2, number of terms); see help file ")
        }
        if(!identical(TRUE, compare(paravars, colnames(prior$support), ignoreOrder = TRUE)$result)) {
          stop("Uniform prior: the column names attribute for support must correspond to the named parameters in the model formula.")
        }
        quad <- RSquadrature.uniform(no.terms, t(prior$support), nr, nq)
        abscissas <- as.matrix(quad$a)
        colnames(abscissas) <- colnames(prior$support)
        weights <- quad$w
        } else {
          stop("Argument \"prior\" must correctly specify a normal or uniform prior for the model parameters (see help file).")
        }
      ## dcw 28-11-2018 check the weights and abscissas are all OK
      if(anyNA(weights) || any(is.infinite(weights)) || any(is.nan(weights))) stop("The quadrature scheme is not valid for large values of nr and/or nq. Try making these values smaller.")
      if(anyNA(abscissas) || any(is.infinite(abscissas)) || any(is.nan(abscissas))) stop("The quadrature scheme is not valid for large values of nr and/or nq. Try making these values smaller.")
  }
  }
  if (criterion == "D") {
    if(identical(method, "MC")) {
      inte <- function(d, B) {
        n1 <- dim(d)[1]
        sam <- prior(B)
        #        d2 <- 0.5 * (d + 1) * diff + lower
        d3 <- matrix(0, ncol = k, nrow = B * n1)
        for (i in 1:k) {
          d3[, i] <- rep(d[, i], B)
        }
        sam3 <- matrix(0, ncol = p, nrow = B * n1)
        for (i in 1:p) {
          sam3[, paravars == paravars[i]] <- rep(sam[, colnames(sam) == paravars[i]], each = n1)
        }
        jac <- attributes(grad(d = d3, paras = sam3))$gradient
        as.vector(.Call("Dnlmcpp", jac, c(n1, B), PACKAGE = "acebayes"))
      }
    } else {
      inte <- function(d, B) {
        n1 <- dim(d)[1]
        sam <- abscissas
        B <- nrow(sam)
        #        d2 <- 0.5 * (d + 1) * diff + lower
        d3 <- matrix(0, ncol = k, nrow = B * n1)
        for (i in 1:k) {
          d3[, i] <- rep(d[, i], B)
        }
        sam3 <- matrix(0, ncol = p, nrow = B * n1)
        for (i in 1:p) {
          # cat("p = ", p, ", rep ... = ", rep(sam[, colnames(sam) == paravars[i]], each = n1), "\n")
          sam3[, paravars == paravars[i]] <- rep(sam[, colnames(sam) == paravars[i]], each = n1)
        }
        jac <- attributes(grad(d = d3, paras = sam3))$gradient
        v <- as.vector(.Call("Dnlmcpp", jac, c(n1, B), PACKAGE = "acebayes"))
        weighted.mean(v, weights)
      }
    }
  }
  if (criterion == "A") {
    if(identical(method, "MC")) {
      inte <- function(d, B) {
        n1 <- dim(d)[1]
        sam <- prior(B)
        #       d2 <- 0.5 * (d + 1) * diff + lower
        d3 <- matrix(0, ncol = k, nrow = B * n1)
        for (i in 1:k) {
          d3[, i] <- rep(d[, i], B)
        }
        sam3 <- matrix(0, ncol = p, nrow = B * n1)
        for (i in 1:p) {
          sam3[, paravars == paravars[i]] <- rep(sam[, 
                                                     colnames(sam) == paravars[i]], each = n1)
        }
        jac <- attributes(grad(d = d3, paras = sam3))$gradient
        as.vector(.Call("Anlmcpp", jac, c(n1, B), PACKAGE = "acebayes"))
      }
    } else {
      inte <- function(d, B) {
        n1 <- dim(d)[1]
        sam <- abscissas
        B <- nrow(sam)
        #        d2 <- 0.5 * (d + 1) * diff + lower
        d3 <- matrix(0, ncol = k, nrow = B * n1)
        for (i in 1:k) {
          d3[, i] <- rep(d[, i], B)
        }
        sam3 <- matrix(0, ncol = p, nrow = B * n1)
        for (i in 1:p) {
          sam3[, paravars == paravars[i]] <- rep(sam[, 
                                                     colnames(sam) == paravars[i]], each = n1)
        }
        jac <- attributes(grad(d = d3, paras = sam3))$gradient
        v <- as.vector(.Call("Anlmcpp", jac, c(n1, B), PACKAGE = "acebayes"))
        weighted.mean(v, weights)
      }
    }
  }
  if (criterion == "E") {
    if(identical(method, "MC")) {
      inte <- function(d, B) {
        n1 <- dim(d)[1]
        sam <- prior(B)
        #        d2 <- 0.5 * (d + 1) * diff + lower
        d3 <- matrix(0, ncol = k, nrow = B * n1)
        for (i in 1:k) {
          d3[, i] <- rep(d[, i], B)
        }
        sam3 <- matrix(0, ncol = p, nrow = B * n1)
        for (i in 1:p) {
          sam3[, paravars == paravars[i]] <- rep(sam[, 
                                                     colnames(sam) == paravars[i]], each = n1)
        }
        jac <- attributes(grad(d = d3, paras = sam3))$gradient
        as.vector(.Call("Enlmcpp", jac, c(n1, B), PACKAGE = "acebayes"))
      } 
    } else {
      inte <- function(d, B) {
        n1 <- dim(d)[1]
        sam <- abscissas
        B <- nrow(sam)
        #        d2 <- 0.5 * (d + 1) * diff + lower
        d3 <- matrix(0, ncol = k, nrow = B * n1)
        for (i in 1:k) {
          d3[, i] <- rep(d[, i], B)
        }
        sam3 <- matrix(0, ncol = p, nrow = B * n1)
        for (i in 1:p) {
          sam3[, paravars == paravars[i]] <- rep(sam[, 
                                                     colnames(sam) == paravars[i]], each = n1)
        }
        jac <- attributes(grad(d = d3, paras = sam3))$gradient
        v <- as.vector(.Call("Enlmcpp", jac, c(n1, B), PACKAGE = "acebayes"))
        weighted.mean(v, weights)
      }
    }
  }
  if (criterion == "SIG") {
    inte <- function(d, B) {
      n1 <- dim(d)[1]
      B2 <- B * 2
      sam <- prior(B2)
      d3 <- matrix(0, ncol = k, nrow = B2 * n1)
      for (i in 1:k) {
        d3[, i] <- rep(d[, i], B2)
      }
      sam3 <- matrix(0, ncol = p, nrow = B2 * n1)
      for (i in 1:p) {
        sam3[, paravars == paravars[i]] <- rep(sam[, 
                                                   colnames(sam) == paravars[i]], each = n1)
      }
      mu1 <- matrix(grad(d = d3, paras = sam3)[1:(B2 * 
                                                    n1)], ncol = n1, byrow = TRUE)
      mu2 <- matrix(mu1[-(1:B), ],nrow=B)
      mu1 <- matrix(mu1[1:B, ],nrow=B)
      y <- matrix(rnorm(n = n1 * B, mean = as.vector(mu1), 
                        sd = rep(sqrt(sam[1:B, colnames(sam) == "sig2"]), 
                                 n1)), ncol = n1)
      as.vector(.Call("SIGnlmcpp", y, mu1, mu2, sam[1:B, 
                                                    colnames(sam) == "sig2"], sam[-(1:B), colnames(sam) == 
                                                                                    "sig2"], PACKAGE = "acebayes"))
    }
  }
  if (criterion == "NSEL") {
    inte <- function(d, B) {
      n1 <- dim(d)[1]
      B2 <- B * 2
      sam <- prior(B2)
      d3 <- matrix(0, ncol = k, nrow = B2 * n1)
      for (i in 1:k) {
        d3[, i] <- rep(d[, i], B2)
      }
      sam3 <- matrix(0, ncol = p, nrow = B2 * n1)
      for (i in 1:p) {
        sam3[, paravars == paravars[i]] <- rep(sam[, 
                                                   colnames(sam) == paravars[i]], each = n1)
      }
      mu1 <- matrix(grad(d = d3, paras = sam3)[1:(B2 * 
                                                    n1)], ncol = n1, byrow = TRUE)
      mu2 <- matrix(mu1[-(1:B), ],nrow=B)
      mu1 <- matrix(mu1[1:B, ],nrow=B)
      y <- matrix(rnorm(n = n1 * B, mean = as.vector(mu1), 
                        sd = rep(sqrt(sam[1:B, colnames(sam) == "sig2"]), 
                                 n1)), ncol = n1)
      as.vector(.Call("NSELnlmcpp", y, mu2, sam[-(1:B), 
                                                colnames(sam) == "sig2"], sam[1:B, colnames(sam) != 
                                                                                "sig2"], sam[-(1:B), colnames(sam) != "sig2"], 
                      PACKAGE = "acebayes"))
    }
  }
  output <- list(utility = inte)
  output
  }

############ acenlm ############################

## if method="quadrature", prior is a list with (a) for normal dist, mean, var, (optionally) nr and ns; 
## (b) for uniform dist, support, (optionally) nr and ns
## mu must have named elements, support must have named columns

acenlm <- function (formula, start.d, prior, B, criterion = c("D", "A", "E", "SIG", "NSEL"), 
                     method = c("quadrature", "MC"), Q = 20, N1 = 20, N2 = 100, lower = -1, upper = 1, 
                     progress = FALSE, limits = NULL) 
{
  
  criterion <- match.arg(criterion)
  if(length(method) > 1) {
    method = switch(EXPR=criterion,
                    D = "quadrature",
                    A = "quadrature",
                    E = "quadrature",
                    "MC")
  } 
  method <- match.arg(method)
  if(identical(method, "MC") && !is.function(prior)) stop("For method = \"MC\", argument prior must specify a function.") 
  
  if(missing(B)) {
    B <- switch(method,
                MC = c(20000, 1000),
                quadrature = c(2, 8))
  }
  
  utilobj <- utilitynlm(formula = formula, prior = prior, desvars = dimnames(start.d)[[2]], 
                         criterion = criterion, method = method, nrq = B) 
  
 inte <- function(d, B) {
    d2 <- 0.5 * (d + 1) * diff + lower
    utilobj$utility(d2, B)
  }
  
  deterministic <- FALSE
  if(identical(method, "quadrature")) deterministic <- TRUE
  mainupper <- upper
  mainlower <- lower
  mainlimits <- limits
  limits2 <- mainlimits
  if (length(mainlimits) != 0) {
    limits2 <- function(d, i, j) {
      d1 <- 0.5 * (d + 1) * (mainupper - mainlower) + mainlower
      out <- mainlimits(d = d1, i = i, j = j)
      if (is.matrix(mainupper)) {
        out1 <- 2 * (out - mainlower[i, j])/(mainupper[i, 
                                                       j] - mainlower[i, j]) - 1
      }
      else {
        out1 <- 2 * (out - mainlower)/(mainupper - mainlower) - 
          1
      }
      out1
    }
  }
  
  diff <- upper - lower
  start.d2 <- 2 * (start.d - lower)/diff - 1
  output <- ace(utility = inte, start.d = start.d2, B = B, 
                 Q = Q, N1 = N1, N2 = N2, lower = -1, upper = 1, progress = progress, 
                 limits = limits2, deterministic = deterministic)
  inte2 <- function(d, B) {
    d1 <- 2 * (d - mainlower)/(mainupper - mainlower) - 1
    inte(d = d1, B = B)
  }
  output$utility <- inte2
  output$glm <- FALSE
  output$nlm <- TRUE
  output$criterion <- criterion
  output$method <- method
  output$prior <- prior
  output$formula <- formula
  output$phase1.d <- 0.5 * (output$phase1.d + 1) * (mainupper - 
                                                      mainlower) + mainlower
  output$phase2.d <- 0.5 * (output$phase2.d + 1) * (mainupper -                 
                         mainlower) + mainlower
output
}

############ utilglm ############################

utilglm<-function(x , beta, family, criterion = "D"){

eta<-beta%*%t(x)

if(!is.function(family)){
	if(is.list(family)){
	stuff<-family} else{
		family2 <- get(family, mode = "function", envir = parent.frame())
		stuff<-family2()}
		} else{
	stuff<-family()}
mu<-stuff$linkinv(eta)
w<-(stuff$mu.eta(eta)^2)/stuff$variance(mu)
if(stuff$family=="gaussian" & stuff$link=="identity"){    ## To fix error with Gaussian("identity")
w<-matrix(1,nrow=nrow(mu),ncol=ncol(mu))}	

if(criterion=="D"){
eval<-.Call("Dcpp", x, w, PACKAGE = "acebayes")}
if(criterion=="A"){
eval<-.Call("Acpp", x, w, PACKAGE = "acebayes")}
if(criterion=="E"){
eval<-.Call("Ecpp", x, w, PACKAGE = "acebayes")}

as.vector(eval)}

############ pval ############################

pval <- function(oldeval, neweval, binary){
if(binary){
old_n<-length(oldeval)
new_n<-length(neweval)
old_sum<-sum(oldeval)
new_sum<-sum(neweval)
new_beta_sam<-rbeta(n = 10000, shape1 = 1 + new_sum, shape2 = 1 + new_n - new_sum)
out<-mean(pbeta(q = new_beta_sam, shape1 = 1 + old_sum, shape2 = 1 + old_n - old_sum))} else{
details<-as.vector(.Call( "pvalcpp", oldeval, neweval, PACKAGE = "acebayes" ))
out<-1-pt(details[1],df=details[2])}
out}

############ utilglm ############################

noisyexpectutil<-function(utility, B, d, i, j, Dij){
temp<-d
z<-c()
for(L in 1:length(Dij)){
temp[i,j]<-Dij[L]
z[L]<-mean(utility(d=temp, B=B))}
list(z=z,Dij=Dij)}

############ distmat ############################

distmat <- function(Dij){
	.Call( "distcpp",Dij,PACKAGE = "acebayes" )
}

############ GPpred ############################

GPpred <- function(paras, dist, z, newDij, Dij){
as.vector(.Call( "GPpredcpp",paras, dist, z, newDij, Dij, PACKAGE = "acebayes" ))
}

############ FisherScoring2D ############################

FisherScoring2D<-function(par, Dij, z, dist, tol = 1e-10, max.iter = 25){

singular<-0
QQ<-length(Dij)

theta<-par

e1<-exp(theta[1])
e2<-exp(theta[2])

G<-exp(-e2*dist)
R<-dist*G
B<-e1*diag(QQ)+G
iB<-c()
try(iB<-solve(B),silent=TRUE)
if(!is.null(iB)){
M1<-iB%*%iB
M2<-iB%*%R
M3<-M2%*%iB
M1z<-as.vector(M1%*%matrix(z,ncol=1))
M3z<-as.vector(M3%*%matrix(z,ncol=1))

D1<-0.5*e1*(sum(M1z*z) - sum(diag(iB)))
D2<-0.5*e2*(sum(diag(M2)) - sum(M3z*z))
A11<--0.5*e1*e1*sum(diag(M1))
A12<-0.5*e1*e2*sum(diag(M3))
A22<--0.5*e2*e2*sum(M2*t(M2))
F<-A11*A22-A12*A12
counter<-1

	while(max(abs(c(D1,D2)))>tol & counter<max.iter & abs(F)<Inf & abs(F)>0 & !any(is.na(c(D1,D2,F)))){
	last<-theta
	theta<-theta-c(A22*D1 - A12*D2,A11*D2 - A12*D1)/F
	e1<-exp(theta[1])
	e2<-exp(theta[2])

	G<-exp(-e2*dist)
	R<-dist*G
	B<-e1*diag(QQ)+G
	iB<-c()
	try(iB<-solve(B),silent=TRUE)
	if(!is.null(iB)){
		M1<-iB%*%iB
		M2<-iB%*%R
		M3<-M2%*%iB
		M1z<-as.vector(M1%*%matrix(z,ncol=1))
		M3z<-as.vector(M3%*%matrix(z,ncol=1))

		D1<-0.5*e1*(sum(M1z*z) - sum(diag(iB)))
		D2<-0.5*e2*(sum(diag(M2)) - sum(M3z*z))
		A11<--0.5*e1*e1*sum(diag(M1))
		A12<-0.5*e1*e2*sum(diag(M3))
		A22<--0.5*e2*e2*sum(M2*t(M2))
		F<-A11*A22-A12*A12

		counter<-counter+1}   else{
		
		counter<-max.iter+2}

	}

if(counter==(max.iter+2)){
theta<-last}

} else{

singular<-1}

list(par=theta, singular=singular)}

############ acephase1 ############################

acephase1 <- function (utility, start.d, B, Q = 20, N1 = 20, 
          lower, upper, limits = NULL, progress = FALSE, binary = FALSE, deterministic = FALSE) 
{
  if(missing(B) && identical(deterministic, FALSE)) B <- c(20000, 1000)
  if(missing(B) && identical(deterministic, TRUE)) B <- NULL

  ptm <- proc.time()[3]
  if (!is.matrix(upper)) {
    UPPER <- upper + 0 * start.d
  }
  else {
    UPPER <- upper
  }
  if (!is.matrix(lower)) {
    LOWER <- lower + 0 * start.d
  }
  else {
    LOWER <- lower
  }
  DIFF <- UPPER - LOWER
  if (length(limits) == 0) {
    limits2 <- function(i, j, d) {
      c(LOWER[i, j], sort(runif(9998)) * DIFF[i, j] + LOWER[i, 
                                                            j], UPPER[i, j])
    }
  }
  else {
    limits2 <- limits
  }
  n <- dim(start.d)[1]
  k <- dim(start.d)[2]
  DESIGN <- start.d
  eval <- utility(d = start.d, B = B[[1]])
  curr <- mean(eval)
  curr2 <- curr
  counter <- 1
  best <- DESIGN
  inner <- DESIGN
  inner_eval <- curr
  best_eval <- inner_eval + 1
  while (counter <= N1) { 
    for (i in 1:n) {
      for (j in 1:k) {
        xxx <- as.vector(optimumLHS(n = Q, k = 1)) * 
          2 - 1
        xxx2 <- 0.5 * (xxx + 1) * DIFF[i, j] + LOWER[i, 
                                                     j]
        yyy <- noisyexpectutil(utility = utility, B = B[[2]], 
                               d = DESIGN, i = i, j = j, Dij = xxx2)$z 
        A.array <- distmat(xxx)
        meany <- mean(yyy)
        sdy <- sd(yyy)
        zzz <- (yyy - meany)/sdy
        optgp <- FisherScoring2D(par = c(0, 0), Dij = xxx, 
                                 z = zzz, dist = A.array)
        opt <- NULL
        if (optgp$singular != 1) {
          xxxz <- limits2(i = i, j = j, d = DESIGN)
          xxxz2 <- 2 * (xxxz - LOWER[i, j])/DIFF[i, j] - 
            1
          yyyz <- meany + sdy * GPpred(paras = optgp$par, 
                                       dist = A.array, z = matrix(zzz, ncol = 1), 
                                       newDij = xxxz2, Dij = xxx)
          opt <- NULL
          start <- xxxz[max(yyyz) == yyyz]
          if (length(start) == 1) {
            opt <- list(best.x = start)
          }
        }
        if (!is.null(opt)) {
          old_DESIGN <- DESIGN
          new_DESIGN <- DESIGN
          new_DESIGN[i, j] <- opt$best.x
          old_eval <- utility(d = old_DESIGN, B = B[[1]])
          new_eval <- utility(d = new_DESIGN, B = B[[1]])
          if(deterministic) { 
            the.p.val <- ifelse(new_eval > old_eval, 1, 0)
          } else {
            the.p.val <- pval(old_eval, new_eval, binary)
            the.p.val <- ifelse(is.na(the.p.val), 0, the.p.val) 
          }
          if (the.p.val >= runif(1)) { 
            DESIGN <- new_DESIGN
            curr <- c(curr, mean(new_eval))
            if (curr[length(curr)] > inner_eval) {
              inner_eval <- curr[length(curr)]
              inner <- new_DESIGN
            }
          }
          else {
            DESIGN <- old_DESIGN
            curr <- c(curr, curr[length(curr)])
          }
        }
        else {
          curr <- c(curr, curr[length(curr)])
        }
      }
    }
    old_DESIGN <- best
    new_DESIGN <- inner
    old_eval <- utility(d = old_DESIGN, B = B[[1]])
    new_eval <- utility(d = new_DESIGN, B = B[[1]])
    if(deterministic) {
      the.p.val <- ifelse(new_eval > old_eval, 1, 0)
#      if(identical(the.p.val, 0)) break
    } else {
      the.p.val <- pval(old_eval, new_eval, binary)
      the.p.val <- ifelse(is.na(the.p.val), 0, the.p.val)
    }
    if (the.p.val >= runif(1)) {
      best_eval <- mean(new_eval)
      best <- inner
    }
    else {
      inner <- best
      best_eval <- mean(old_eval)
      inner_eval <- best_eval
    }
    curr2 <- c(curr2, best_eval)
    if (progress) {
      cat("Phase I iteration ", counter, " out of ", N1, 
          " (Current value = ", best_eval, ") \n", sep = "")
#      print(best)
    }
#    if(deterministic && identical(the.p.val, 0)) break
    counter <- counter + 1
  
  }
  ptm <- proc.time()[3] - ptm
  output <- list(utility = utility, start.d = start.d, phase1.d = best, 
                 phase2.d = best, phase1.trace = curr2, phase2.trace = NULL, B=B, 
                 Q = Q, N1 = N1, N2 = 0, glm = FALSE, nlm = FALSE, criterion = "NA", 
                 family = "NA", prior = "NA", time = ptm, binary = binary, deterministic = deterministic)
  class(output) <- "ace"
  output
}

############ acephase2 ############################

acephase2 <- function (utility, start.d, B, N2 = 100, progress = FALSE, binary = FALSE, 
         deterministic = FALSE) 
{
 if (missing(B) && identical(deterministic, FALSE)) 
   B <- c(20000, 1000)
 if (missing(B) && identical(deterministic, TRUE)) 
   B <- NULL
 ptm <- proc.time()[3]
 n <- dim(start.d)[1]
 k <- dim(start.d)[2]
 DESIGN <- start.d
 CAND <- DESIGN
 eval <- utility(d = DESIGN, B = B[[1]])
 best <- DESIGN
 best_ob <- eval
 curr2 <- mean(eval)
 counter2 <- 0 ## changed from 1 to ensure we go in the while loop
 # if (progress) {
 #   cat("Phase II iteration ", counter2, " out of ", N2, 
 #       " (Current value = ", curr2, ") \n", sep = "")
 # }
 while (counter2 < N2) {
   crt <- c()
   for (j in 1:n) {
     INTER <- rbind(DESIGN, CAND[j, ])
     crt[j] <- mean(utility(d = INTER, B = B[[2]]))
   }
   potpts <- (1:n)[max(crt) == crt]
   if (length(potpts) > 1) {
     potpts <- sample(x = potpts, size = 1)
   }
   INTER <- rbind(DESIGN, CAND[potpts, ])
   crt <- c()
   for (j in 1:(n + 1)) {
     INTER2 <- matrix(INTER[-j, ], nrow = n, dimnames = dimnames(CAND)) ## adjusted to fix missing colnames
     crt[j] <- mean(utility(d = INTER2, B = B[[2]]))
   }
   potpts <- (1:(n + 1))[max(crt) == crt]
   if (length(potpts) > 1) {
     potpts <- sample(x = potpts, size = 1)
   }
   new_DESIGN <- matrix(INTER[-potpts, ], nrow = n, dimnames = dimnames(CAND))
   old_DESIGN <- DESIGN
   old_eval <- utility(d = old_DESIGN, B = B[[1]])
   new_eval <- utility(d = new_DESIGN, B = B[[1]])
   if (deterministic) {
     the.p.val <- ifelse(new_eval > old_eval, 1, 0)
   }
   else {
     the.p.val <- pval(old_eval, new_eval, binary)
     the.p.val <- ifelse(is.na(the.p.val), 0, the.p.val)
   }
   if (the.p.val > runif(1)) {
     curr2 <- c(curr2, mean(new_eval))
     DESIGN <- new_DESIGN
     the.p.val <- pval(best_ob, new_eval, binary)
     the.p.val <- ifelse(is.na(the.p.val), 0, the.p.val)
     if (the.p.val >= runif(1)) {
       best <- DESIGN
       best_ob <- new_eval
     }
   }
   else {
     curr2 <- c(curr2, mean(old_eval))
   }
   counter2 <- counter2 + 1
   if (progress) {
     cat("Phase II iteration ", counter2, " out of ", 
         N2, " (Current value = ", curr2[length(curr2)], 
         ") \n", sep = "")
   }
 }

curr2<-curr2[-1]  # Delete 1st element of curr2
 ptm <- proc.time()[3] - ptm
 output <- list(utility = utility, start.d = start.d, phase1.d = start.d, 
                phase2.d = best, phase1.trace = NULL, phase2.trace = curr2, 
                B = B, Q = NULL, N1 = 0, N2 = N2, glm = FALSE, nlm = FALSE, 
                criterion = "NA", family = "NA", prior = "NA", time = ptm, 
                binary = binary, deterministic = deterministic)
 class(output) <- "ace"
 output
}


# acephase2 <- function (utility, start.d, B, N2 = 100, progress = FALSE, 
#           binary = FALSE, deterministic = FALSE) 
# {
#   if(missing(B) && identical(deterministic, FALSE)) B <- c(20000, 1000)
#   if(missing(B) && identical(deterministic, TRUE)) B <- NULL
#   
#   ptm <- proc.time()[3]
#   n <- dim(start.d)[1]
#   k <- dim(start.d)[2]
#   DESIGN <- start.d
#   CAND <- DESIGN
#   eval <- utility(d = DESIGN, B = B[[1]])
#   best <- DESIGN
#   best_ob <- eval
#   curr2 <- mean(eval)
#   counter2 <- 1
#   if (progress) {
#     cat("Phase II iteration ", counter2, " out of ", N2, 
#         " (Current value = ", curr2, ") \n", sep = "")
#   }
#   while (counter2 < N2) {
#     crt <- c()
#     for (j in 1:n) {
#       INTER <- rbind(DESIGN, CAND[j, ])
#       crt[j] <- mean(utility(d = INTER, B = B[[2]]))
#     }
#     potpts <- (1:n)[max(crt) == crt]
#     if (length(potpts) > 1) {
#       potpts <- sample(x = potpts, size = 1)
#     }
#     INTER <- rbind(DESIGN, CAND[potpts, ])
#     crt <- c()
#     for (j in 1:(n + 1)) {
#       INTER2 <- as.matrix(INTER[-j, ], nrow = n)
#       crt[j] <- mean(utility(d = INTER2, B = B[[2]]))
#     }
#     potpts <- (1:(n + 1))[max(crt) == crt]
#     if (length(potpts) > 1) {
#       potpts <- sample(x = potpts, size = 1)
#     }
#     new_DESIGN <- matrix(INTER[-potpts, ], nrow = n, dimnames = dimnames(CAND))
#     old_DESIGN <- DESIGN
#     old_eval <- utility(d = old_DESIGN, B = B[[1]])
#     new_eval <- utility(d = new_DESIGN, B = B[[1]])
#     if(deterministic) {
#       the.p.val <- ifelse(new_eval > old_eval, 1, 0) 
#     } else {
#       the.p.val <- pval(old_eval, new_eval, binary)
#       the.p.val <- ifelse(is.na(the.p.val), 0, the.p.val)
#     }
#     if (the.p.val > runif(1)) {
#       curr2 <- c(curr2, mean(new_eval))
#       DESIGN <- new_DESIGN
#       the.p.val <- pval(best_ob, new_eval, binary)
#       the.p.val <- ifelse(is.na(the.p.val), 0, the.p.val)
#       if (the.p.val >= runif(1)) {
#         best <- DESIGN
#         best_ob <- new_eval
#       }
#     }
#     else {
#       curr2 <- c(curr2, mean(old_eval))
#     }
#     counter2 <- counter2 + 1
#     if (progress) {
#       cat("Phase II iteration ", counter2, " out of ", 
#           N2, " (Current value = ", curr2[length(curr2)], 
#           ") \n", sep = "")
#     }
# #    if(deterministic && counter2 > 2) {
# #      if(!(curr2[counter2 - 1] > curr2[counter2 - 2] + tolerence)) break
# #    }
#   }
#   ptm <- proc.time()[3] - ptm
#   output <- list(utility = utility, start.d = start.d, phase1.d = start.d, 
#                  phase2.d = best, phase1.trace = NULL, phase2.trace = curr2, B=B, 
#                  Q = NULL, N1 = 0, N2 = N2, glm = FALSE, nlm = FALSE, 
#                  criterion = "NA", family = "NA", prior = "NA", time = ptm, 
#                  binary = binary, deterministic = deterministic)
#   class(output) <- "ace"
#   output
# }

############ ace ############################

ace <- function (utility, start.d, B, Q = 20, N1 = 20, 
          N2 = 100, lower = -1, upper = 1, limits = NULL, progress = FALSE, 
          binary = FALSE, deterministic = FALSE) 
{
  if(missing(B) && identical(deterministic, FALSE)) B <- c(20000, 1000)
  if(missing(B) && identical(deterministic, TRUE)) B <- NULL
  
  ptm <- proc.time()[3]
  if (N1 > 0) {
    interim <- acephase1(utility = utility, start.d = start.d, 
                         B = B, Q = Q, N1 = N1, lower = lower, upper = upper, 
                         limits = limits, progress = progress, binary = binary, 
                         deterministic = deterministic)
    interim.d <- interim$phase1.d
    interim.trace <- interim$phase1.trace
  }
  else {
    interim.d <- start.d
    interim.trace <- NULL
  }
  if (N2 > 0) {
    last <- acephase2(utility = utility, start.d = interim.d, 
                      B = B, N2 = N2, progress = progress, binary = binary, 
                      deterministic = deterministic)
    last.d <- last$phase2.d
    last.trace <- last$phase2.trace
  }
  else {
    last.d <- interim.d
    last.trace <- NULL
  }
  ptm <- proc.time()[3] - ptm
  output <- list(utility = utility, start.d = start.d, phase1.d = interim.d, 
                 phase2.d = last.d, phase1.trace = interim.trace, phase2.trace = last.trace, 
                 B = B, Q = Q, N1 = N1, N2 = N2, glm = FALSE, nlm = FALSE, criterion = "NA", 
                 prior = "NA", time = ptm, binary = binary, deterministic = deterministic)
  class(output) <- "ace"
  output
}

############ aceglm ############################

## if method="quadrature", prior is a list with (a) for normal dist, mean, var, (optionally) nr and ns; 
## (b) for uniform dist, support, (optionally) nr and ns

aceglm <- function (formula, start.d, family, prior, B, criterion = c("D", "A", "E", "SIG", "NSEL", "SIG-Norm", "NSEL-Norm"), 
                     method = c("quadrature", "MC"), Q = 20, N1 = 20, N2 = 100, lower = -1, upper = 1, 
                     progress = FALSE, limits = NULL) 
{
	
# if (is.character(family)) 
#         family <- get(family, mode = "function", envir = parent.frame())
#     if (is.function(family)) 
#         family <- family()
#     if (is.null(family$family)) {
#         print(family)
#         stop("'family' not recognized")
#     }	
	
  criterion <- match.arg(criterion)
  if(length(method) > 1) {
    method = switch(EXPR=criterion,
                    D = "quadrature",
                    A = "quadrature",
                    E = "quadrature",
                    "MC")
  } 
  method <- match.arg(method)
  if(identical(method, "MC") && !is.function(prior)) stop("For method = \"MC\", argument prior must specify a function.") 
  
  if(missing(B)) {
   B <- switch(method,
          MC = c(20000, 1000),
          quadrature = c(2, 8))
  }
  
  utilobj <- utilityglm(formula = formula, family = family, prior = prior, 
                         criterion = criterion, method = method, nrq = B) 
  
  inte <- function(d, B) {
    utilobj$utility(d, B)
  }
  
  deterministic <- FALSE
  if(identical(method, "quadrature")) deterministic <- TRUE
  
  output <- ace(utility = inte, start.d = start.d, B = B, Q = Q, 
                N1 = N1, N2 = N2, lower = lower, upper = upper, progress = progress, 
                limits = limits, deterministic = deterministic)
  output$glm <- TRUE
  output$criterion <- criterion
  output$method <- method
  output$prior <- prior
  output$phase1.d <- output$phase1.d
  output$phase2.d <- output$phase2.d
  output$family <- family
  output$formula <- formula
output

}

############ plot.ace ############################

plot.ace<-function(x,...){

if(length(x$phase1.trace)>0 & length(x$phase2.trace)>0){
ulim<-max(c(x$phase1.trace,x$phase2.trace))
llim<-min(c(x$phase1.trace,x$phase2.trace))
plot(0:(length(x$phase1.trace)-1),x$phase1.trace,xlab="Phase I iteration",ylab="Observation of expected utility",ylim=c(llim,ulim),type="l",xlim=c(0,length(x$phase1.trace)))
new_z<-axTicks(side=1)
new_x<-(1:length(x$phase2.trace))*((length(x$phase1.trace)-1)/length(x$phase2.trace))
lines(new_x,x$phase2.trace,col=8)
legend(x="bottomright",legend=c("Phase I","Phase II"),col=c(1,8),lty=c(1,1),bty="n")
axis(side=3,labels=new_z/((length(x$phase1.trace)-1)/length(x$phase2.trace)),at=new_z)
mtext("Phase II iteration", side=3, line = par("mgp")[1]) }

if(length(x$phase1.trace)>0 & length(x$phase2.trace)==0){
plot(0:(length(x$phase1.trace)-1),x$phase1.trace,type="l",xlab="Phase I iteration",ylab="Observation of expected utility")
legend(x="bottomright",legend=c("Phase I"),col=1,lty=1,bty="n") }

if(length(x$phase1.trace)==0 & length(x$phase2.trace)>0){
plot(1:length(x$phase2.trace),x$phase2.trace,type="l",xlab="Phase II iteration",ylab="Observation of expected utility",col=8)
legend(x="bottomright",legend=c("Phase II"),col=8,lty=1,bty="n") }

}

############ print.ace ############################

print.ace<-function(x,...){

hrs<-round(x$time%/%3600,0)
mins<-round((x$time%%3600)%/%60,0)
secs<-round((x$time%%3600)%%60,0)
hrs<-ifelse(hrs<10,paste("0",hrs,sep=""),hrs)
mins<-ifelse(mins<10,paste("0",mins,sep=""),mins)
secs<-ifelse(secs<10,paste("0",secs,sep=""),secs)

if(x$glm==TRUE){
  cat("Generalised Linear Model \n")
  cat("Criterion = Bayesian ",x$criterion,"-optimality \n",sep="")
  cat("Formula: "); print(x$formula)
  if(is.function(x$family)){
  print(x$family())} else{
  print(x$family)}	
  cat("Method: ", x$method, "\n\n")
  if(identical(x$method, "MC")) cat("B: ", x$B, "\n\n")
  else {
    cat("nr = ", x$B[[1]], ", nq = ", x$B[[2]],"\n")
    if(identical(names(x$prior)[1:2], c("mu", "sigma2"))) cat("Prior: normal\n\n")
    if(identical(names(x$prior)[1], "support")) cat("Prior: uniform\n\n")       
  }
} 
if(x$nlm==TRUE){
  cat("Non Linear Model \n")
  cat("Criterion = Bayesian ",x$criterion,"-optimality \n",sep="")
  cat("Formula: "); print(x$formula)
  cat("Method: ", x$method, "\n\n")
  if(identical(x$method, "MC")) cat("B: ", x$B, "\n\n")
  else {
    cat("nr = ", x$B[[1]], ", nq = ", x$B[[2]],"\n")
    if(identical(names(x$prior)[1:2], c("mu", "sigma2"))) cat("Prior: normal\n\n")
    if(identical(names(x$prior)[1], "support")) cat("Prior: uniform\n\n")       
  }
} 
if(x$nlm==FALSE & x$glm==FALSE){
cat("User-defined model & utility \n")}
#cat("\n")
cat("Number of runs = ",dim(x$phase2.d)[1],"\n",sep="")
cat("\n")
cat("Number of factors = ",dim(x$phase2.d)[2],"\n",sep="")
cat("\n")
cat("Number of Phase I iterations = ",x$N1,"\n",sep="")
cat("\n")
cat("Number of Phase II iterations = ",x$N2,"\n",sep="")
cat("\n")
cat("Computer time = ",paste(hrs,":",mins,":",secs,sep=""),"\n",sep="")

}

############ summary.ace ############################

summary.ace<-function(object,...){

print.ace(x=object)}

################################ Utility Functions ###########################################################

################################################################################################################################
### 4.1 Linear model
################################################################################################################################

utillinmod<-function(d, B){
x<-cbind(1,d,d^2,d[,1]*d[,2])
log_deter<-as.vector(.Call("LMcpp", x, PACKAGE = "acebayes"))
log_deter+rnorm(B)}

optdeslinmod<-function(n, type = "ACE"){
if(type=="ACE"){
des<-linmodoptdesigns[linmodoptdesigns$n==n,]} else{
des<-truelinmodoptdesigns[truelinmodoptdesigns$n==n,]}
des<-as.matrix(des)
des<-des[,-1]
rownames(des)<-NULL
des}

################################################################################################################################
### 4.2 Compartmental model
################################################################################################################################

### n=18 Gotwalt et al model Bayesian D-optimality

utilcomp18bad<-function(d, B){
low<-c(0.01884,0.298)
upp<-c(0.09884,8.298)
sam<-cbind(runif(n=B,min=low[1],max=upp[1]),runif(n=B,min=low[2],max=upp[2]))
as.vector(.Call("utilcomp18badcpp", d , sam, PACKAGE = "acebayes"))}

optdescomp18bad<-function(type = "ACE"){
if(type=="ACE"){
des<-comp18badoptdesign}
if(type=="Gotwalt"){
des<-comp18badoptdesign_got}
if(type=="Atkinson"){
des<-comp18badoptdesign_atk}
des<-as.matrix(des)
rownames(des)<-NULL
des}

### n=15 Ryan et al model SIG

utilcomp15bad<-function(d, B){
d2<-(d+1)*12
theta<-exp(cbind(rnorm(n=B,mean=log(0.1),sd=sqrt(0.05)),rnorm(n=B,mean=log(1),sd=sqrt(0.05)),rnorm(n=B,mean=log(20),sd=sqrt(0.05))))
as.vector(.Call("utilcomp15badcpp", d2 , theta, PACKAGE = "acebayes"))}

optdescomp15bad<-function(){
des<-comp15badoptdesign
des<-as.matrix(des)
rownames(des)<-NULL
des}

### n=15 Ryan et al model SIG

#inidescomp15sig<-function(rep){
#des<-comp15badinidesign$time[comp15badinidesign$rep==rep]
#des<-as.matrix(des)
#rownames(des)<-NULL
#des}

utilcomp15sig<-function(d, B){
D<-400
sigadd<-0.1
sigpro<-0.01
d2<-12*(as.vector(d)+1)
n1<-length(d2)
sam<-cbind(rnorm(n=2*B,mean=log(0.1),sd=sqrt(0.05)),rnorm(n=2*B,mean=log(1),sd=sqrt(0.05)),rnorm(n=2*B,mean=log(20),sd=sqrt(0.05)))
sam<-exp(sam)
mu<-(D/matrix(rep(sam[,3],n1),ncol=n1))*(matrix(rep(sam[,2],n1),ncol=n1)/(matrix(rep(sam[,2],n1),ncol=n1)-matrix(rep(sam[,1],n1),ncol=n1)))*(exp(-matrix(rep(sam[,1],n1),ncol=n1)*matrix(rep(d2,2*B),ncol=n1,byrow=TRUE))-exp(-matrix(rep(sam[,2],n1),ncol=n1)*matrix(rep(d2,2*B),ncol=n1,byrow=TRUE)))
vv<-sigadd+sigpro*(mu^2)
y<-matrix(rnorm(n=n1*B,mean=as.vector(mu[1:B,]),sd=sqrt(vv[1:B,])),ncol=n1)

frho<-as.vector(.Call("rowSumscpp", log(vv[-(1:B),]), PACKAGE = "acebayes"))
loglik<-as.vector(.Call("rowSumscpp", matrix(dnorm(x=as.vector(y),mean=as.vector(mu[1:B,]),sd=sqrt(vv[1:B,]),log=TRUE),ncol=n1), PACKAGE = "acebayes"))

rsll4<-as.vector(.Call("utilcomp15sigcpp", y, mu[-(1:B),], vv[-(1:B),], frho, PACKAGE = "acebayes"))
MY3<-log(rsll4/B)

eval<-loglik-MY3

eval}

optdescomp15sig<-function(){
des<-comp15sigoptdesign
des<-as.matrix(des)
rownames(des)<-NULL
des}

### n=15 Ryan et al model SIG DRS

utilcomp15sigDRS<-function(d, B){
D<-400
sigadd<-0.1
sigpro<-0.01
d2<-24*qbeta(p=seq(from=0,to=1,length=17)[-c(1,17)],shape1=d[1,1],shape2=d[2,1])
n1<-length(d2)
sam<-cbind(rnorm(n=2*B,mean=log(0.1),sd=sqrt(0.05)),rnorm(n=2*B,mean=log(1),sd=sqrt(0.05)),rnorm(n=2*B,mean=log(20),sd=sqrt(0.05)))
sam<-exp(sam)
mu<-(D/matrix(rep(sam[,3],n1),ncol=n1))*(matrix(rep(sam[,2],n1),ncol=n1)/(matrix(rep(sam[,2],n1),ncol=n1)-matrix(rep(sam[,1],n1),ncol=n1)))*(exp(-matrix(rep(sam[,1],n1),ncol=n1)*matrix(rep(d2,2*B),ncol=n1,byrow=TRUE))-exp(-matrix(rep(sam[,2],n1),ncol=n1)*matrix(rep(d2,2*B),ncol=n1,byrow=TRUE)))
vv<-sigadd+sigpro*(mu^2)
y<-matrix(rnorm(n=n1*B,mean=as.vector(mu[1:B,]),sd=sqrt(vv[1:B,])),ncol=n1)

frho<-as.vector(.Call("rowSumscpp", log(vv[-(1:B),]), PACKAGE = "acebayes"))
loglik<-as.vector(.Call("rowSumscpp", matrix(dnorm(x=as.vector(y),mean=as.vector(mu[1:B,]),sd=sqrt(vv[1:B,]),log=TRUE),ncol=n1), PACKAGE = "acebayes"))

rsll4<-as.vector(.Call("utilcomp15sigcpp", y, mu[-(1:B),], vv[-(1:B),], frho, PACKAGE = "acebayes"))
MY3<-log(rsll4/B)

eval<-loglik-MY3

eval}

optdescomp15sigDRS<-function(){
des<-comp15sigDRSoptdesign
des<-as.matrix(des)
rownames(des)<-NULL
des}

################################################################################################################################
### 4.3 Logistic Regression
################################################################################################################################

### Standard Logistic Regression - Bayesian D-optimality

utillrbad<-function(d, B){
low<-c(-3,4,5,-6,-2.5)
upp<-c(3,10,11,0,3.5)
beta<-t(t(6*matrix(runif(n=5*B),ncol=5))+low)
x<-cbind(1,d)
eval<-as.vector(.Call("LRDcpp", x, beta, PACKAGE = "acebayes"))
eval}

optdeslrbad<-function(n, type = "ACE"){
if(type=="ACE"){
des<-LRBADoptdesign[LRBADoptdesign$n==n,]
des<-as.matrix(des)
des<-des[,-1]}
if(type=="Gotwalt1"){
des<-LRBADoptdesign2[LRBADoptdesign2$n==n & LRBADoptdesign2$type==0,]}
if(type=="Gotwalt2"){
des<-LRBADoptdesign2[LRBADoptdesign2$n==n & LRBADoptdesign2$type==1,]}
if(type=="Woods"){
des<-LRBADoptdesign2[LRBADoptdesign2$n==n & LRBADoptdesign2$type==2,]}
if(type!="ACE"){
des<-as.matrix(des)
des<-des[,-c(1,2)]}
rownames(des)<-NULL
des}

### Hierarchical Logistic Regression - Bayesian D-optimality

utilhlrbad<-function(d, B){
x<-cbind(1,d)
m<-6
n<-dim(x)[1]
G<-n/m
low<-c(-3,4,5,-6,-2.5)
upp<-c(3,10,11,0,3.5)
bet<-c(3,3,3,1,1)
beta<-t(t(6*matrix(runif(n=5*B),ncol=5))+low)
marker<-rep(1:G,each=m)
S<-t(bet*t(matrix(1-sqrt(runif(5*B)),ncol=5)))
gam<-matrix(0,ncol=G*5,nrow=B)
z<-matrix(0,nrow=n,ncol=G*5)
for(i in 1:G){
for(j in 1:5){
z[marker==i,5*(i-1)+j]<-x[marker==i,j]
gam[,5*(i-1)+j]<-runif(n=B,min=-S[,j],max=S[,j])}}

#eval<-as.vector(test.cpp(x=x,z=z,beta=beta,gam=gam,S=S))
eval<-as.vector(.Call("HLRDcpp", x, z, beta, gam, S, PACKAGE = "acebayes"))
eval}

optdeshlrbad<-function(n){
des<-HLRBADoptdesign[HLRBADoptdesign$n==n,] 
des<-as.matrix(des)
des<-des[,-1]
rownames(des)<-NULL
des}

### Standard Logistic Regression - SIG

inideslrsig<-function(n, rep){
des<-LRSIGinidesign[LRSIGinidesign$n==n & LRSIGinidesign$rep==rep,-c(1,2)]
des<-as.matrix(des)
rownames(des)<-NULL
des}

utillrsig<-function(d, B){
low<-c(-3,4,5,-6,-2.5)
upp<-c(3,10,11,0,3.5)

sam<-matrix(runif(n=2*B*5),ncol=5)
for(jj in 1:5){
sam[,jj]<-sam[,jj]*(upp[jj]-low[jj])+low[jj]}

x<-cbind(1,d)
n1<-dim(x)[1]

rho<-1/(1+exp(-sam%*%t(x)))
Z<-log(1-rho[-(1:B),])
frho<-as.vector(.Call("rowSumscpp", Z, PACKAGE = "acebayes"))

y<-matrix(rbinom(n=n1*B,size=1,prob=as.vector(rho[1:B,])),ncol=n1)
Z<-dbinom(x=y,size=1,prob=rho[1:B,],log=TRUE)                        ## loglik
rsll<-as.vector(.Call("rowSumscpp", Z, PACKAGE = "acebayes"))
sam<-sam[-(1:B),]

rsll4<-as.vector(.Call("siglrcpp", y, x, sam, frho, PACKAGE = "acebayes"))

MY3<-log(rsll4/B)
eval<-rsll-MY3

eval}

optdeslrsig<-function(n){
des<-LRSIGoptdesign[LRSIGoptdesign$n==n,]
des<-as.matrix(des)
des<-des[,-1]
rownames(des)<-NULL
des}

### Hierarchical Logistic Regression - SIG

inideshlrsig<-function(n, rep){
des<-LRSIGinidesign[LRSIGinidesign$n==n & LRSIGinidesign$rep==rep,-c(1,2)]
des<-as.matrix(des)
rownames(des)<-NULL
des}

utilhlrsig<-function(d, B){
B2d<-B
B4d<-B
x<-cbind(1,d)
m<-6
n1<-dim(x)[1]
G<-n1/m
low<-c(-3,4,5,-6,-2.5)
upp<-c(3,10,11,0,3.5)
bet<-c(3,3,3,1,1)
sam<-t(t(6*matrix(runif(n=5*(B2d+2*B4d)),ncol=5))+low)
marker<-rep(1:G,each=m)
zam<-matrix(0,ncol=G*5,nrow=B2d+2*B4d)
z<-matrix(0,nrow=n1,ncol=G*5)
for(i in 1:G){
for(j in 1:5){
z[marker==i,5*(i-1)+j]<-x[marker==i,j]
zam[,5*(i-1)+j]<-bet[j]*(1-sqrt(runif(B2d+2*B4d)))*(2*runif(B2d+2*B4d)-1)}}
etax<-sam%*%t(x)
etaz<-zam%*%t(z)
rho<-1/(1+exp(-etax-etaz))
y<-matrix(rbinom(n=B2d*n1,size=1,prob=as.vector(rho[1:B2d,])),ncol=n1)

frho<-as.vector(.Call("rowSumscpp", log(1-rho[((B2d+1):(B2d+B4d)),]), PACKAGE="acebayes"))

rsll4<-.Call("sighlrcpp", cbind(x,z), y, cbind(sam[-(1:B2d),],zam[-(1:B2d),]), frho, sam[1:B2d,], exp(etax[1:B2d,]), exp(etaz[-(1:(B2d+B4d)),]), PACKAGE="acebayes")

rsll4<-log(rsll4)

eval<-rsll4[,2]-rsll4[,1]

eval}

optdeshlrsig<-function(n){
des<-HLRSIGoptdesign[HLRSIGoptdesign$n==n,]
des<-as.matrix(des)
des<-des[,-1]
rownames(des)<-NULL
des}

### Standard Logistic Regression - Bayesian A-optimality

utillrbaa<-function(d, B){
low<-c(-3,4,5,-6,-2.5)
upp<-c(3,10,11,0,3.5)
beta<-t(t(6*matrix(runif(n=5*B),ncol=5))+low)
x<-cbind(1,d)
eval<-as.vector(.Call("LRAcpp", x, beta, PACKAGE = "acebayes"))
eval}

optdeslrbaa<-function(n){
des<-LRBAAoptdesign[LRBAAoptdesign$n==n,]
des<-as.matrix(des)
des<-des[,-1]
rownames(des)<-NULL
des}

### Hierarchical Logistic Regression - Bayesian A-optimality

utilhlrbaa<-function(d, B){
x<-cbind(1,d)
m<-6
n<-dim(x)[1]
G<-n/m
low<-c(-3,4,5,-6,-2.5)
upp<-c(3,10,11,0,3.5)
bet<-c(3,3,3,1,1)
beta<-t(t(6*matrix(runif(n=5*B),ncol=5))+low)
marker<-rep(1:G,each=m)
S<-t(bet*t(matrix(1-sqrt(runif(5*B)),ncol=5)))
gam<-matrix(0,ncol=G*5,nrow=B)
z<-matrix(0,nrow=n,ncol=G*5)
for(i in 1:G){
for(j in 1:5){
z[marker==i,5*(i-1)+j]<-x[marker==i,j]
gam[,5*(i-1)+j]<-runif(n=B,min=-S[,j],max=S[,j])}}

#eval<-as.vector(test.cpp(x=x,z=z,beta=beta,gam=gam,S=S))
eval<-as.vector(.Call("HLRAcpp", x, z, beta, gam, S, PACKAGE = "acebayes"))
eval}

optdeshlrbaa<-function(n){
des<-HLRBAAoptdesign[HLRBAAoptdesign$n==n,]
des<-as.matrix(des)
des<-des[,-1]
rownames(des)<-NULL
des}

### Standard Logistic Regression - NSEL

inideslrnsel<-function(n, rep){
des<-LRNSELinidesign[LRNSELinidesign$n==n & LRNSELinidesign$rep==rep,-c(1,2)]
des<-as.matrix(des)
rownames(des)<-NULL
des}

utillrnsel<-function(d, B){
low<-c(-3,4,5,-6,-2.5)
upp<-c(3,10,11,0,3.5)

sam<-matrix(runif(n=2*B*5),ncol=5)
for(jj in 1:5){
sam[,jj]<-sam[,jj]*(upp[jj]-low[jj])+low[jj]}

x<-cbind(1,d)
n1<-dim(x)[1]

rho<-1/(1+exp(-sam%*%t(x)))
Z<-log(1-rho[-(1:B),])
frho<-as.vector(.Call("rowSumscpp", Z, PACKAGE = "acebayes"))

y<-matrix(rbinom(n=n1*B,size=1,prob=as.vector(rho[1:B,])),ncol=n1)
fsam<-sam[1:B,]
sam<-sam[-(1:B),]

rsll4<-.Call("nsellrcpp", y, x, sam, frho, PACKAGE = "acebayes")

Z<-(rsll4-fsam)^2
eval<-as.vector(.Call("rowSumscpp", Z, PACKAGE = "acebayes"))

-eval}

optdeslrnsel<-function(n){
des<-LRNSELoptdesign[LRNSELoptdesign$n==n,]
des<-as.matrix(des)
des<-des[,-1]
rownames(des)<-NULL
des}

### Hierarchical Logistic Regression - NSEL

inideshlrnsel<-function(n, rep){
des<-LRSIGinidesign[LRSIGinidesign$n==n & LRSIGinidesign$rep==rep,-c(1,2)]
des<-as.matrix(des)
rownames(des)<-NULL
des}

utilhlrnsel<-function(d, B){
x<-cbind(1,d)
m<-6
n1<-dim(x)[1]
G<-n1/m
low<-c(-3,4,5,-6,-2.5)
upp<-c(3,10,11,0,3.5)
bet<-c(3,3,3,1,1)
sam<-t(t(6*matrix(runif(n=5*2*B),ncol=5))+low)
marker<-rep(1:G,each=m)
zam<-matrix(0,ncol=G*5,nrow=2*B)
z<-matrix(0,nrow=n1,ncol=G*5)
for(i in 1:G){
for(j in 1:5){
z[marker==i,5*(i-1)+j]<-x[marker==i,j]
zam[,5*(i-1)+j]<-bet[j]*(1-sqrt(runif(2*B)))*(2*runif(2*B)-1)}}
etax<-sam%*%t(x)
etaz<-zam%*%t(z)
rho<-1/(1+exp(-etax-etaz))
y<-matrix(rbinom(n=B*n1,size=1,prob=as.vector(rho[1:B,])),ncol=n1)
frho<-as.vector(.Call("rowSumscpp", log(1-rho[-(1:B),]), PACKAGE = "acebayes"))
rsll4<-.Call("nselhlrcpp", cbind(x,z), y, cbind(sam[-(1:B),],zam[-(1:B),]), frho, PACKAGE = "acebayes")

MY3<-(rsll4-sam[1:B,])^2
eval<-as.vector(.Call("rowSumscpp", MY3, PACKAGE = "acebayes"))

-eval
}

optdeshlrnsel<-function(n){
des<-HLRNSELoptdesign[HLRNSELoptdesign$n==n,]
des<-as.matrix(des)
des<-des[,-1]
rownames(des)<-NULL
des}

################################################################################################################################
### 4.4 Beetles
################################################################################################################################

utilbeetle<-function(d, B){
nd<-dim(d)[1]
L<-1.6907
U<-1.8839
R<-U-L
ck<-log(-log(0.5));lambda<-60
Xda<-cbind(1,(d+1)/2)
Xda[,2]<-Xda[,2]*R+L
Xda2<-cbind(Xda,Xda[,2]^2)
msam<-sample(x=1:6,size=B,prob=as.vector(probs),replace=TRUE)
tmsam<-c(length(msam[msam==1]),length(msam[msam==2]),length(msam[msam==3]),length(msam[msam==4]),length(msam[msam==5]),length(msam[msam==6]))
wr1<-sample(x=1:10000,size=tmsam[1],replace=TRUE)
wr2<-sample(x=1:10000,size=tmsam[2],replace=TRUE)
wr3<-sample(x=1:10000,size=tmsam[3],replace=TRUE)
wr4<-sample(x=1:10000,size=tmsam[4],replace=TRUE)
wr5<-sample(x=1:10000,size=tmsam[5],replace=TRUE)
wr6<-sample(x=1:10000,size=tmsam[6],replace=TRUE)
sam1<-as.matrix(lrlin[wr1,],ncol=2)
sam2<-as.matrix(lrquad[wr2,],ncol=3)
sam3<-as.matrix(colin[wr3,],ncol=2)
sam4<-as.matrix(coquad[wr4,],ncol=3)
sam5<-as.matrix(prlin[wr5,],ncol=2)
sam6<-as.matrix(prquad[wr6,],ncol=3)
eta1<-sam1%*%t(Xda)
eta2<-sam2%*%t(Xda2)
eta3<-sam3%*%t(Xda)
eta4<-sam4%*%t(Xda2)
eta5<-sam5%*%t(Xda)
eta6<-sam6%*%t(Xda2)
phi0<-c(-sam1[,1]/sam1[,2],
0.5*(-sam2[,2]+sqrt(sam2[,2]^2 - 4*sam2[,1]*sam2[,3]))/sam2[,3],
(ck-sam3[,1])/sam3[,2],
0.5*(-sam4[,2]+sqrt(sam4[,2]^2 - 4*(sam4[,1]-ck)*sam4[,3]))/sam4[,3],
-sam5[,1]/sam5[,2],
0.5*(-sam6[,2]+sqrt(sam6[,2]^2 - 4*sam6[,1]*sam6[,3]))/sam6[,3])
mu<-lambda*rbind(1/(1+exp(-eta1)),1/(1+exp(-eta2)),1-exp(-exp(eta3)),1-exp(-exp(eta4)),pnorm(eta5),pnorm(eta6))
y<-matrix(rpois(n=prod(dim(mu)),lambda=as.vector(mu)),ncol=nd)

wr1<-sample(x=1:10000,size=tmsam[1],replace=TRUE)
wr2<-sample(x=1:10000,size=tmsam[2],replace=TRUE)
wr3<-sample(x=1:10000,size=tmsam[3],replace=TRUE)
wr4<-sample(x=1:10000,size=tmsam[4],replace=TRUE)
wr5<-sample(x=1:10000,size=tmsam[5],replace=TRUE)
wr6<-sample(x=1:10000,size=tmsam[6],replace=TRUE)
sam1<-as.matrix(lrlin[wr1,],ncol=2)
sam2<-as.matrix(lrquad[wr2,],ncol=3)
sam3<-as.matrix(colin[wr3,],ncol=2)
sam4<-as.matrix(coquad[wr4,],ncol=3)
sam5<-as.matrix(prlin[wr5,],ncol=2)
sam6<-as.matrix(prquad[wr6,],ncol=3)
eta1<-sam1%*%t(Xda)
eta2<-sam2%*%t(Xda2)
eta3<-sam3%*%t(Xda)
eta4<-sam4%*%t(Xda2)
eta5<-sam5%*%t(Xda)
eta6<-sam6%*%t(Xda2)
mu1<-1/(1+exp(-eta1))
mu2<-1/(1+exp(-eta2))
mu3<-1-exp(-exp(eta3))
mu4<-1-exp(-exp(eta4))
mu5<-pnorm(eta5)
mu6<-pnorm(eta6)

phi<-c(-sam1[,1]/sam1[,2],
0.5*(-sam2[,2]+sqrt(sam2[,2]^2 - 4*sam2[,1]*sam2[,3]))/sam2[,3],
(ck-sam3[,1])/sam3[,2],
0.5*(-sam4[,2]+sqrt(sam4[,2]^2 - 4*(sam4[,1]-ck)*sam4[,3]))/sam4[,3],
-sam5[,1]/sam5[,2],
0.5*(-sam6[,2]+sqrt(sam6[,2]^2 - 4*sam6[,1]*sam6[,3]))/sam6[,3])

mu<-lambda*rbind(mu1,mu2,mu3,mu4,mu5,mu6)
Lmu<-log(mu)
fy<-lfactorial(y)

frho<--as.vector(.Call("rowSumscpp", mu, PACKAGE = "acebayes"))
ncr<--as.vector(.Call("rowSumscpp", fy, PACKAGE = "acebayes"))

rsll<-as.vector(.Call("beetlecpp", phi, y, Lmu, frho, ncr, PACKAGE = "acebayes"))

eval<-(phi0-rsll)^2

names(eval)<-NULL

-eval}

optdesbeetle<-function(n){
des<-beetleoptdesign[beetleoptdesign$n==n,]
des<-as.matrix(des)
des<-matrix(des[,-1],nrow=n)
rownames(des)<-NULL
des}

################# Quadrature Functions

# A simplex function to numerically approximate the inner integral. 
simplex <- function(p) {
  # cretate matrix v
  v <- matrix(0, nrow = p+1, ncol = p)
  #extended simplex method begins by creating a p+1 vertex-centered simplex in Rp. 
  for (i in 1:(p+1)) {
    for (j in 1:p) {
      if (j<i) { v[i,j] = (-1) * sqrt( (p+1) / ( p * (p-j+2) * (p-j+1) ) ) }
      if (j==i) { v[i,j] = sqrt( (p+1) * (p-i+1) / ( p * (p-i+2)  ) ) }
      if (j>i) { v[i,j] = 0 }
    }
  }
  # create midpoints of the simplex vertices. 
  midpoints <- matrix(ncol=p, nrow=p*(p+1)/2)
  k<-1
  for (i in 1:p) {
    for (j in (i+1):(p+1)) {
      midpoints[k,] = 0.5 * (v[i,] + v[j,])
      k<-k+1
    }
  }
  proj.pts <- midpoints # gets correct dimensions
  for (k in 1:(p*(p+1)/2)) 
  {
    norm <- (sum(midpoints[k,]^2))^0.5
    proj.pts[k,] <- midpoints[k,]/norm
    if(identical(proj.pts[k, ], NaN)) proj.pts[k,] <- 0
  }
  
  # Form extended simplex by adding in negative images of these points
  simplex <- rbind(v, -v, proj.pts, -proj.pts)
  
  # simplex weights
  w.s <- rep(0,(p+1)*(p+2))
  w.s[1:(2*(p+1))] <- p*(7-p) / (2*((p+1)^2) * (p+2))
  w.s[-(1:(2*(p+1)))] <- 2*((p-1)^2) / (p*((p+1)^2) * (p+2))
  return(list(simplex=simplex,w.s=w.s))
}
# sm1 <- simplex(p=3)

# A function which finds the roots of gaussian laugerre and return them as a vector. 
# Numerical values of abscissas a. 
# Code found at Press et al 1992 as C++ code.
gaulag<-function(p, Nr, its, precision) { 
  alpha <- p/2 
  x <- c(0)
  for(i in 1:Nr){ 
    if(i==1) { z <- (1 + alpha) * (3 + 0.92*alpha) / (1 + 2.4*Nr + 1.8*alpha) } 
    if(i==2) { z <- z + (15 + 6.25*alpha) / (1 + 2.4*Nr + 1.8*alpha) } 
    if(i>2) { ai <- i-2 } 
    if(i>2)  { z <- z + ((1 + 2.55*ai) / (1.9*ai) + 1.26*ai*alpha / (1+3.5*ai)) * (z - x[i-2]) / (1 + 0.3*alpha) } 
    
    for(its in 1:its){ 
      p1 <- 1
      p2 <- 0 
      for(j in 1:Nr){ 
        p3 <- p2
        p2 <- p1 
        p1 <- ( ( 2*j - 1 + alpha - z ) * p2 - ( j -1 + alpha ) * p3 ) / (j)
      } 
      pp <- ( Nr*p1 - (Nr +  alpha)*p2 ) / z 
      z1 <- z 
      z <- z1 - p1 / pp 
      if(abs(z-z1) < precision) break 
      # if(its >= its) warning("too many iterations in gaussian laguerre")
    } 
    x[i] <- z 
  } 
  x<-c(0,x) 
  return(x) 
} 
# a <-gaulag(p=3, Nr=10, its=1e+06, precision=1e-06)

# A function for the associated Laguerre polynomials.
laguerre <- function( a, n, s ) {
  laguerre.matrix <- matrix(0, nrow = length(a), ncol = (n+1))
  for (i in 0:n) {
    laguerre.matrix[,i+1] <- (-1)^i * exp( lfactorial(n) + lchoose(n+s, n-i) + log(a^i) - lfactorial(i) )
  }
  laguerre.vector <- rowSums(laguerre.matrix)
  return(laguerre.vector)
}

# A function for the weights.
# Numerical values for the coefficients H.
# Formula found at Cassity 1965.
cassity.weight <- function(a,p) {
  s <- (p/2) - 1
  n <- length(a)
  constant <- exp( lgamma(n) + lgamma(n+s) - log(n+s) )
  weight <- constant/(laguerre(a,n-1,s)^2) 
  weight[1] <- weight[1] * (1+s)
  
  return(weight)
}
# w <- laguerre(p=3, a=a)
# l <- cassity.weight(a=a, p=3)

# normal prior quadrature
library(randtoolbox)
RSquadrature <- function(p, mu, Sigma, Nr, Nq) {
  
  if(!identical(length(mu), as.integer(p))) stop("length(mu) must equal p")
  if(!identical(length(Sigma), as.integer(p * p))) stop("Sigma must have size p * p")
  
  Q <- array(dim=c(p, p, Nq))
  big.ran.mat <- matrix(halton(p * p * Nq, dim = 1, normal = T), ncol = p)
  for (i in 1:Nq) {
    ran.mat <- big.ran.mat[(1 + (i - 1) * p):(p + (i - 1) * p), ]
    qr <- qr(ran.mat)
    Q[,,i]<-qr.Q(qr)
  }
  
  # cholesky root of Sigma
  L <- t(chol(Sigma))
  # transpose because we want the lower triangular cholensky root. (entries above diagonal to be zero)
  
  # Spherical inner integral: simplex and simplex weights 
  simp <- simplex(p)
  simplex <- simp$simplex
  w.s <- simp$w.s
  
  # Radial outer integral
  # find abscissas
  r <- gaulag(p=p, Nr=Nr, its=1e+06, precision=1e-06)
  # find coefficients/weights
  w.R <- exp( log(cassity.weight(p=p, a=r)) - lgamma((p)/2) )
  r <- 2*r
  
  #create arrays to store abscissas and weights, and put values for the zero abscissas
  a <- matrix(nrow=( 1+(p+1)*(p+2)*Nr*Nq) , ncol=p)
  a[1,] <- mu
  w.a <- rep(0, 1+(p+1)*(p+2)*Nr*Nq)
  w.a[1] <- w.R[1]
  l <- 2
  
  # calculate remaining abscissas and weights
  for (i in 1:Nr) {
    for (j in 1:((p+1)*(p+2))) {
      for (k in 1:Nq) {
        # compute the abscissa in beta-log(Sigma^2) space
        a[l,] <- mu + (r[i+1]^0.5)*L%*%Q[,,k]%*%simplex[j,]
        w.a[l] <- w.R[i+1]*w.s[j]/Nq
        l <- l+1 
      }
    }
  }
  return(list(a = a, w = w.a))
}

# uniform prior quadrature
RSquadrature.uniform <- function(p, limits, Nr, Nq) {
  
  #generate orthogonal matrices
  Q <- array(dim=c(p, p, Nq))
  big.ran.mat <- matrix(halton(p * p * Nq, dim = 1, normal = T), ncol = p)
  for (i in 1:Nq) {
    #		ran.mat<-matrix(rnorm(P*P),ncol=P)
    ran.mat <- big.ran.mat[(1 + (i - 1) * p):(p + (i - 1) * p), ]
    qr <- qr(ran.mat)
    Q[,,i]<-qr.Q(qr)
  }
  
  
  # compute simplex + weights
  simp <- simplex(p)
  simplex <- simp$simplex
  w.s <- simp$w.s
  
  #print(simplex)
  #print(w.s)
   
  # compute radial abscissae, store in vector r
  
  r <- gaulag(p=p, Nr=Nr, its=1e6, precision=1e-6)
  w.R <- exp( log(cassity.weight(p=p, a=r)) - lgamma((p)/2) )
  r <- 2*r
  #print(r)
  
  d1<-limits[,2]-limits[,1]
  d2<-limits[,1]
  #create arrays to store abscissae and weights, and put values for the zero abscissa
  a <- matrix(nrow=( 1+(p+1)*(p+2)*Nr*Nq), ncol=p)
  a[1,] <-d1*pnorm(0)+d2  
  
  w.a <- rep(0, 1+(p+1)*(p+2)*Nr*Nq)
  w.a[1] <- w.R[1]
  
  l <- 2
  
  # calculate remaining abscissae and weights
  for (i in 1:Nr) {
    for (j in 1:((p+1)*(p+2))) {
      for (k in 1:Nq) {
        # compute the abscissa in beta-log(Sigma^2) space
        a[l,] <- d1*pnorm((r[i+1]^0.5)*Q[,,k]%*%simplex[j,])+d2
        w.a[l] <- w.R[i+1]*w.s[j]/Nq
        l <- l+1
        # exponentiate log(Sigma^2) to put on beta-Sigma^2 scale
        #ONLY FOR GLMM case
      }
    }
  }
  return(list(a = a, w = w.a))
}

# RSquadrature <- function(p, mu, Sigma, Nr, Nq) {
# 
#   if(!identical(length(mu), as.integer(p))) stop("length(mu) must equal p")
#   if(!identical(length(Sigma), as.integer(p * p))) stop("Sigma must have size p * p")
# 	#generate orthogonal matrices
# 	#P<- p+1
# 	P<-p
# 	Q <- array(dim=c(Nq,P,P))
# 	big.ran.mat <- matrix(halton(P * P * Nq, dim = 1, normal = T), ncol = P)
# 	for (i in 1:Nq) {
# #	  ran.mat<-matrix(rnorm(P*P),ncol=P)
# 	  ran.mat <- big.ran.mat[(1 + (i - 1) * P):(P + (i - 1) * P), ]
# 		qr <- qr(ran.mat)
# 		Q[i,,]<-qr.Q(qr)
# 	}
# 
# 
# 	# compute L, Cholesky root of Sigma
# 	# L <- sqrt(Sigma)
#   L <- t(chol(Sigma))
# 
# 	# compute simplex + weights
# 	simp <-simplex(P)
# 	simplex <- simp$simplex
# 	w.s <- simp$w.s
# 
# 	#print(simplex)
# 	#print(w.s)
# 
# 	# compute radial abscissae, store in vector r
# 
# 	r <- gaulag(p=P,Nr=Nr,its=1e6,precision=1e-6)
# 	w.R <- cassity.weight(r,P)/gamma((P)/2)
# 	r <- 2*r
# 	#print(r)
# 
# 	#create arrays to store abscissae and weights, and put values for the zero abscissa
# 	a<- matrix(mu,ncol=P)
# 	#a[1,P] <<- exp(mu[P])
# 
# 	w.a <- NULL
# 	w.a <- w.R[1]
# 
# 	# calculate remaining abscissae and weights
# 	for (i in 1:Nr) {
# 			for (j in 1:((P+1)*(P+2))) {
# 				for (k in 1:Nq) {
# 					# compute the abscissa in beta-log(sigma^2) space
# 					theta <- mu + (r[i+1]^0.5)*L%*%Q[k,,]%*%simplex[j,]
# 
# 					# exponentiate log(sigma^2) to put on beta-sigma^2 scale
# 					#ONLY FOR GLMM case
# 					#theta[P] <- exp(theta[P])
# 					a <- rbind(a,t(theta))
# 					w.a <- c(w.a,w.R[i+1]*w.s[j]/Nq)
# 				}
# 			}
# 		}
# 	return(list(a = a, w = w.a))
# }
# 
# 
# RSquadrature.uniform <- function(p, limits, Nr, Nq) {
# 
# 	#generate orthogonal matrices
# 	#P<- p+1
# 	P<-p
# 	Q <- array(dim=c(Nq,P,P))
# 	big.ran.mat <- matrix(halton(P * P * Nq, dim = 1, normal = T), ncol = P)
# 	for (i in 1:Nq) {
# #		ran.mat<-matrix(rnorm(P*P),ncol=P)
#     ran.mat <- big.ran.mat[(1 + (i - 1) * P):(P + (i - 1) * P), ]
# 	  qr <- qr(ran.mat)
# 		Q[i,,]<-qr.Q(qr)
# 	}
# 
# 
# 
# 	# compute simplex + weights
# 	simp <-simplex(P)
# 	simplex <- simp$simplex
# 	w.s <- simp$w.s
# 
# 	#print(simplex)
# 	#print(w.s)
# 
# 	# compute radial abscissae, store in vector r
# 
# 	r <- gaulag(p=P,Nr=Nr,its=1e6,precision=1e-6)
# 	w.R <- cassity.weight(r,P)/gamma((P)/2)
# 	r <- 2*r
# 	#print(r)
# 
# 	d1<-limits[,2]-limits[,1]
# 	d2<-limits[,1]
# 
# 	#create arrays to store abscissae and weights, and put values for the zero abscissa
# 	a <- matrix(d1*pnorm(0)+d2,ncol=P)
# 	#a[1,P] <<- exp(mu[P])
# 
# 	w.a <- NULL
# 	w.a <- w.R[1]
# 
# 	# calculate remaining abscissae and weights
# 	for (i in 1:Nr) {
# 			for (j in 1:((P+1)*(P+2))) {
# 				for (k in 1:Nq) {
# 					# compute the abscissa in beta-log(sigma^2) space
# 					theta <- d1*pnorm((r[i+1]^0.5)*Q[k,,]%*%simplex[j,])+d2
# 
# 					# exponentiate log(sigma^2) to put on beta-sigma^2 scale
# 					#ONLY FOR GLMM case
# 					#theta[P] <- exp(theta[P])
# 					a <- rbind(a,t(theta))
# 					w.a <- c(w.a,w.R[i+1]*w.s[j]/Nq)
# 				}
# 			}
# 		}
# 	return(list(a = a, w = w.a))
# }
# 
# 
# 
# 
# simplex <- function(p) {
# 	V <- matrix(ncol=p,nrow=p+1)
# 	for (i in 1:(p+1))
# 	{
# 		for (j in 1:p) {
# 			if (j<i) { V[i,j] = -((p+1)/(p*(p-j+2)*(p-j+1)))^0.5 }
# 			if (j==i) { V[i,j] = ((p+1)*(p-i+1)/(p*(p-i+2)))^0.5}
# 			if (j>i) { V[i,j] = 0 }
# 		}
# 	}
# 
# 	# Form midpoints and project onto the sphere
# 
# 	midpoints <- matrix(ncol=p,nrow=p*(p+1)/2)
# 	k<-1
# 	for (i in 1:p) {
# 		for (j in (i+1):(p+1)) {
# 			midpoints[k,] = 0.5*(V[i,]+V[j,])
# 			k<-k+1
# 		}
# 	}
# 	proj.pts <- midpoints # gets correct dimensions
# 	for (k in 1:(p*(p+1)/2))
# 	{
# 		norm <- (sum(midpoints[k,]^2))^0.5
# 		proj.pts[k,] <- midpoints[k,]/norm
# 		if(identical(proj.pts[k, ], NaN)) proj.pts[k,] <- 0
# 	}
# 
# 	# Form extended simplex by adding in negative images of these points
# 	simplex <- rbind(V,-V,proj.pts,-proj.pts)
# 
# 	# Compute the simplex weights
# 	w.s <- vector(length=(p+1)*(p+2))
# 	w.s[1:(2*(p+1))] <- p*(7-p)/(2*(p+1)^2*(p+2))
# 	w.s[-(1:(2*(p+1)))] <- 2*((p-1)^2)/(p*(p+1)^2*(p+2))
# 	return(list(simplex=simplex,w.s=w.s))
# }
# 
# # John's Laguerre poly root finder
# # translated from the C code in Press et al.
# # I have checked this against the latter, TW 20.1.2010.
# 
# gaulag<-function(p, Nr, its,precision)
# {
# alpha<-p/2
# a<-vector()
# w<-vector()
# for(i in 1:Nr){
# if(i==1){z<-(1+alpha)*(3+0.92*alpha)/(1+2.4*Nr+1.8*alpha)}
# if(i==2){z<-z+(15+6.25*alpha)/(1+2.4*Nr+1.8*alpha)}
# if(i>2){ai<-i-2}
# if(i>2){z<-z+((1+2.55*ai)/(1.9*ai)+1.26*ai*alpha/(1+3.5*ai))*(z-a[i-2])/(1+0.3*alpha)}
# 
# for(its in 1:its){
# p1<-1;p2<-0
# for(j in 1:Nr){
# p3<-p2;p2<-p1
# p1<-((2*j-1+alpha-z)*p2-(j-1+alpha)*p3)/j}
# 
# pp<-(Nr*p1-(Nr+alpha)*p2)/z
# z1<-z
# z<-z1-p1/pp
# if(abs(z-z1)< precision) break
# }
# a[i]<-z
# }
# a<-c(0,a)
# return(a)
# }
# #Evaluates the laguerre polynomial with parameters n and s at the value a.
# laguerre<-function(a,n,s)
# {
# laguerre.matrix<-matrix(nrow=length(a), ncol=n+1)
# laguerre.vector<-vector(length=length(a))
# #Loops up the recurrence relation
# for(j in 1:length(a)){
# for(i in 0:n){
# laguerre.matrix[j,i+1]<-factorial(n)*choose(n+s, n-i)*(-a[j])^i/factorial(i)
# 
# 
# }
# laguerre.vector[j]<-sum(laguerre.matrix[j,])
# }
# return(laguerre.vector)
# }
# 
# #Reproduces the weights formula given by Cassity(1965). But takes
# #the number of parameters as input to make it easier for the function user.
# cassity.weight<-function(a,p)
# {
# n<-length(a);s<-p/2-1
# constant<-gamma(n)*gamma(n+s)/(n+s)
# weight<-constant/(laguerre(a,n-1,s)^2)
# weight[1]<-weight[1]*(1+s)#as described by cassity et al 'incorporate factor 1+s
# 					# TW: Yeah, I find that sentence a bit confusing.
# 					# but multiplying by 1+s here gives us the right answers
# 					# (checking against the tables...)
# return(weight)
# }

##### pace ######################

pace<-function(utility, start.d, B, Q = 20, N1 = 20, N2 = 100, lower = -1, upper = 1,
	limits = NULL, binary = FALSE, deterministic = FALSE, mc.cores = 1, n.assess = 20){

ptm<-proc.time()[3]	
	
C<-length(start.d)

if (missing(B) && identical(deterministic, FALSE)){ 
        B <- c(20000, 1000)}
if (missing(B) && identical(deterministic, TRUE)){ 
        B <- NULL}

innerace<-function(d){
out<-ace(utility = utility, start.d = d, B = B, Q = Q, N1 = N1, N2 = N2, lower = lower,
	upper = upper, limits = limits, binary = binary, deterministic = deterministic, progress = FALSE)
list(phase2.d=out$phase2.d,phase1.trace=out$phase1.trace,phase2.trace=out$phase2.trace)}

if(.Platform$OS.type=="unix"){
rout<-mclapply(start.d, innerace, mc.cores = mc.cores)} else{
if(mc.cores==1){
rout<-lapply(start.d, innerace)}
if(mc.cores>1){
warning("mc.cores > 1 not currently supported under a non-Unix OS. Proceeding with mc.cores = 1 \n")
rout<-lapply(start.d, innerace)}
}

# cl <- makeCluster(getOption("cl.cores", mc.cores))
# clusterExport(cl, list("ace","utility", "start.d", "B", "Q", "N1", "N2", "lower","upper","limits", "binary","deterministic"), envir=environment())
# rout<-parLapply(cl = cl, X = start.d, fun = innerace)
# stopCluster(cl)

routd<-list()
routphase1<-list()
routphase2<-list()
for(i in 1:C){
routd[[i]]<-rout[[i]]$phase2.d	
routphase1[[i]]<-rout[[i]]$phase1.trace	
routphase2[[i]]<-rout[[i]]$phase2.trace}
	
if(!deterministic){
inner<-function(d){
evals<-rep(0,n.assess)	
for(i in 1:n.assess){
evals[i]<-mean(utility(d = d, B = B[1]))}
evals}

if(.Platform$OS.type=="unix"){
fout<-mclapply(routd, inner, mc.cores = mc.cores)} else{
if(mc.cores==1){
fout<-lapply(routd, inner)}
if(mc.cores>1){
#warning("mc.cores > 1 not currently supported under a non-Unix OS. Proceeding with mc.cores = 1 \n")
fout<-lapply(routd, inner)}
}

# cl <- makeCluster(getOption("cl.cores", mc.cores))
# clusterExport(cl, list("utility", "routd", "B"), envir=environment())
# fout<-parLapply(cl = cl, X = routd, fun = inner)
# stopCluster(cl)

mout<-lapply(fout, mean)
besti<-which.max(mout)}	
if(deterministic){
inner<-function(d){
utility(d = d, B = B)}

if(.Platform$OS.type=="unix"){
fout<-mclapply(routd, inner, mc.cores = mc.cores)} else{
if(mc.cores==1){
fout<-lapply(routd, inner)}
if(mc.cores>1){
#warning("mc.cores > 1 not currently supported under a non-Unix OS. Proceeding with mc.cores = 1 \n")
fout<-lapply(routd, inner)}
}

# cl <- makeCluster(getOption("cl.cores", mc.cores))
# clusterExport(cl, list("utility", "routd", "B"), envir=environment())
# fout<-parLapply(cl = cl, X = routd, fun = inner)
# stopCluster(cl)

besti<-which.max(fout)}	

ptm<-proc.time()[3]-ptm

phase1.trace<-NULL
phase2.trace<-NULL
if(N1>0){phase1.trace<-routphase1[[besti]]}
if(N2>0){phase2.trace<-routphase2[[besti]]}

output<-list(d = routd[[besti]], phase1.trace = phase1.trace, 
phase2.trace = phase2.trace, eval = fout[[besti]],
utility = utility, start.d = start.d, final.d = routd, besti = besti,
B = B, Q = Q, N1 = N1, N2 = N2, glm = FALSE, nlm = FALSE, criterion = "NA",
prior = "NA", time = ptm, binary = binary, deterministic = deterministic)	
class(output)<-"pace"
output
}

####### paceglm ##############

paceglm<-function(formula, start.d, family, prior, B, criterion = c("D", "A", "E", "SIG", "NSEL", "SIG-Norm", "NSEL-Norm"), 
	method = c("quadrature", "MC"), Q = 20, N1 = 20, N2 = 100, lower = -1, upper = 1,
	limits = NULL, mc.cores = 1, n.assess = 20){

ptm<-proc.time()[3]	
# 
# if (is.character(family)) 
#         family <- get(family, mode = "function", envir = parent.frame())
#     if (is.function(family)) 
#         family <- family()
#     if (is.null(family$family)) {
#         print(family)
#         stop("'family' not recognized")
#     }
	
C<-length(start.d)

criterion <- match.arg(criterion)
  if(length(method) > 1) {
    method = switch(EXPR=criterion,
                    D = "quadrature",
                    A = "quadrature",
                    E = "quadrature",
                    "MC")
  } 
  method <- match.arg(method)
  if(identical(method, "MC") && !is.function(prior)) stop("For method = \"MC\", argument prior must specify a function.") 

if(missing(B)) {
   B <- switch(method,
          MC = c(20000, 1000),
          quadrature = c(2, 8))
}
  
utilobj <- utilityglm(formula = formula, family = family, prior = prior, 
                         criterion = criterion, method = method, nrq = B) 
  
inte <- function(d, B) {
    utilobj$utility(d, B)
}

 deterministic <- FALSE
  if(identical(method, "quadrature")) deterministic <- TRUE
    
innerace<-function(d){
out<-aceglm(formula=formula, start.d = d, family=family, prior = prior, B = B, 
	criterion = criterion, method = method, Q = Q, N1 = N1, N2 = N2, lower = lower,
	upper = upper, limits = limits, progress = FALSE)
list(phase2.d=out$phase2.d,phase1.trace=out$phase1.trace,phase2.trace=out$phase2.trace)}

if(.Platform$OS.type=="unix"){
rout<-mclapply(start.d, innerace, mc.cores = mc.cores)} else{
if(mc.cores==1){
rout<-lapply(start.d, innerace)}
if(mc.cores>1){
warning("mc.cores > 1 not currently supported under a non-Unix OS. Proceeding with mc.cores = 1 \n")
rout<-lapply(start.d, innerace)}
}
# cl <- makeCluster(getOption("cl.cores", mc.cores))
# clusterExport(cl, list("aceglm","formula", "start.d", "family", "prior", "B", "criterion", "method", "Q", "N1", "N2", "lower","upper","limits"), envir=environment())
# rout<-parLapply(cl = cl, X = start.d, fun = innerace)
# stopCluster(cl)

routd<-list()
routphase1<-list()
routphase2<-list()
for(i in 1:C){
routd[[i]]<-rout[[i]]$phase2.d	
routphase1[[i]]<-rout[[i]]$phase1.trace	
routphase2[[i]]<-rout[[i]]$phase2.trace}
	
if(!deterministic){
inner<-function(d){
evals<-rep(0,n.assess)	
for(i in 1:n.assess){
evals[i]<-mean(inte(d = d, B = B[1]))}
evals}

if(.Platform$OS.type=="unix"){
fout<-mclapply(routd, inner, mc.cores = mc.cores)} else{
if(mc.cores==1){
fout<-lapply(routd, inner)}
if(mc.cores>1){
#warning("mc.cores > 1 not currently supported under a non-Unix OS. Proceeding with mc.cores = 1 \n")
fout<-lapply(routd, inner)}
}
# cl <- makeCluster(getOption("cl.cores", mc.cores))
# clusterExport(cl, list("inte", "routd", "B"), envir=environment())
# fout<-parLapply(cl = cl, X = routd, fun = inner)
# stopCluster(cl)
mout<-lapply(fout, mean)
besti<-which.max(mout)}	
if(deterministic){
inner<-function(d){
inte(d = d, B = B)}
if(.Platform$OS.type=="unix"){
fout<-mclapply(routd, inner, mc.cores = mc.cores)} else{
if(mc.cores==1){
fout<-lapply(routd, inner)}
if(mc.cores>1){
#warning("mc.cores > 1 not currently supported under a non-Unix OS. Proceeding with mc.cores = 1 \n")
fout<-lapply(routd, inner)}
}
# cl <- makeCluster(getOption("cl.cores", mc.cores))
# clusterExport(cl, list("inte", "routd", "B"), envir=environment())
# fout<-parLapply(cl = cl, X = routd, fun = inner)
# stopCluster(cl)
besti<-which.max(fout)}	

ptm<-proc.time()[3]-ptm

phase1.trace<-NULL
phase2.trace<-NULL
if(N1>0){phase1.trace<-routphase1[[besti]]}
if(N2>0){phase2.trace<-routphase2[[besti]]}

output<-list(d = routd[[besti]], phase1.trace = phase1.trace, 
phase2.trace = phase2.trace, eval = fout[[besti]],
utility = inte, start.d = start.d, final.d = routd, besti = besti,
B = B, Q = Q, N1 = N1, N2 = N2, glm = TRUE, nlm = FALSE, criterion = criterion,
prior = prior, time = ptm, binary = FALSE, deterministic = deterministic,
method = method, family = family, formula = formula)	
class(output)<-"pace"
output
}

##### pacenlm #####

pacenlm<-function(formula, start.d, prior, B, criterion = c("D", "A", "E", "SIG", "NSEL"), 
	method = c("quadrature", "MC"), Q = 20, N1 = 20, N2 = 100, lower = -1, upper = 1,
	limits = NULL, mc.cores = 1, n.assess = 20){

ptm<-proc.time()[3]	
	
C<-length(start.d)

criterion <- match.arg(criterion)
  if(length(method) > 1) {
    method = switch(EXPR=criterion,
                    D = "quadrature",
                    A = "quadrature",
                    E = "quadrature",
                    "MC")
  } 
  method <- match.arg(method)
  if(identical(method, "MC") && !is.function(prior)) stop("For method = \"MC\", argument prior must specify a function.") 

if(missing(B)) {
   B <- switch(method,
          MC = c(20000, 1000),
          quadrature = c(2, 8))
}
  
utilobj <- utilitynlm(formula = formula, prior = prior, desvars = dimnames(start.d[[1]])[[2]],
                         criterion = criterion, method = method, nrq = B) 
  
 deterministic <- FALSE
  if(identical(method, "quadrature")) deterministic <- TRUE
    
innerace<-function(d){
out<-acenlm(formula=formula, start.d = d, prior = prior, B = B, 
	criterion = criterion, method = method, Q = Q, N1 = N1, N2 = N2, lower = lower,
	upper = upper, limits = limits, progress = FALSE)
list(phase2.d=out$phase2.d,phase1.trace=out$phase1.trace,phase2.trace=out$phase2.trace, utility = out$utility)}

if(.Platform$OS.type=="unix"){
rout<-mclapply(start.d, innerace, mc.cores = mc.cores)} else{
if(mc.cores==1){
rout<-lapply(start.d, innerace)}
if(mc.cores>1){
warning("mc.cores > 1 not currently supported under a non-Unix OS. Proceeding with mc.cores = 1 \n")
rout<-lapply(start.d, innerace)}
}

# cl <- makeCluster(getOption("cl.cores", mc.cores))
# clusterExport(cl, list("acenlm","formula", "start.d", "prior", "B", "criterion", "method", "Q", "N1", "N2", "lower","upper","limits"), envir=environment())
# rout<-parLapply(cl = cl, X = start.d, fun = innerace)
# stopCluster(cl)

routd<-list()
routphase1<-list()
routphase2<-list()
for(i in 1:C){
routd[[i]]<-rout[[i]]$phase2.d	
routphase1[[i]]<-rout[[i]]$phase1.trace	
routphase2[[i]]<-rout[[i]]$phase2.trace}

inte<-function(d, B){
rout[[1]]$utility(d=d,B=B)}	
	
if(!deterministic){
inner<-function(d){
evals<-rep(0,n.assess)	
for(i in 1:n.assess){
evals[i]<-mean(inte(d = d, B = B[1]))}
evals}

if(.Platform$OS.type=="unix"){
fout<-mclapply(routd, inner, mc.cores = mc.cores)} else{
if(mc.cores==1){
fout<-lapply(routd, inner)}
if(mc.cores>1){
#warning("mc.cores > 1 not currently supported under a non-Unix OS. Proceeding with mc.cores = 1 \n")
fout<-lapply(routd, inner)}
}

# cl <- makeCluster(getOption("cl.cores", mc.cores))
# clusterExport(cl, list("inte", "routd", "B"), envir=environment())
# fout<-parLapply(cl = cl, X = routd, fun = inner)
# stopCluster(cl)

mout<-lapply(fout, mean)
besti<-which.max(mout)}	
if(deterministic){
inner<-function(d){
inte(d = d, B = B)}

if(.Platform$OS.type=="unix"){
fout<-mclapply(routd, inner, mc.cores = mc.cores)} else{
if(mc.cores==1){
fout<-lapply(routd, inner)}
if(mc.cores>1){
#warning("mc.cores > 1 not currently supported under a non-Unix OS. Proceeding with mc.cores = 1 \n")
fout<-lapply(routd, inner)}
}

# cl <- makeCluster(getOption("cl.cores", mc.cores))
# clusterExport(cl, list("inte", "routd", "B"), envir=environment())
# fout<-parLapply(cl = cl, X = routd, fun = inner)
# stopCluster(cl)
besti<-which.max(fout)}	

ptm<-proc.time()[3]-ptm

phase1.trace<-NULL
phase2.trace<-NULL
if(N1>0){phase1.trace<-routphase1[[besti]]}
if(N2>0){phase2.trace<-routphase2[[besti]]}

output<-list(d = routd[[besti]], phase1.trace = phase1.trace, 
phase2.trace = phase2.trace, eval = fout[[besti]],
utility = inte, start.d = start.d, final.d = routd, besti = besti,
B = B, Q = Q, N1 = N1, N2 = N2, glm = FALSE, nlm = TRUE, criterion = criterion,
prior = prior, time = ptm, binary = FALSE, deterministic = deterministic,
method = method, formula = formula)	
class(output)<-"pace"
output
}

### pace generics #####

print.pace<-function(x, ...){
	
hrs<-round(x$time%/%3600,0)
mins<-round((x$time%%3600)%/%60,0)
secs<-round((x$time%%3600)%%60,0)
hrs<-ifelse(hrs<10,paste("0",hrs,sep=""),hrs)
mins<-ifelse(mins<10,paste("0",mins,sep=""),mins)
secs<-ifelse(secs<10,paste("0",secs,sep=""),secs)

if(x$glm==TRUE){
  cat("Generalised Linear Model \n")
  cat("Criterion = Bayesian ",x$criterion,"-optimality \n",sep="")
  cat("Formula: "); print(x$formula)
   if(is.function(x$family)){
  print(x$family())} else{
  print(x$family)}	
  cat("Method: ", x$method, "\n\n")
  if(identical(x$method, "MC")) cat("B: ", x$B, "\n\n")
  else {
    cat("nr = ", x$B[[1]], ", nq = ", x$B[[2]],"\n")
    if(identical(names(x$prior)[1:2], c("mu", "sigma2"))) cat("Prior: normal\n\n")
    if(identical(names(x$prior)[1], "support")) cat("Prior: uniform\n\n")       
  }
} 
if(x$nlm==TRUE){
  cat("Non Linear Model \n")
  cat("Criterion = Bayesian ",x$criterion,"-optimality \n",sep="")
  cat("Formula: "); print(x$formula)
  cat("Method: ", x$method, "\n\n")
  if(identical(x$method, "MC")) cat("B: ", x$B, "\n\n")
  else {
    cat("nr = ", x$B[[1]], ", nq = ", x$B[[2]],"\n")
    if(identical(names(x$prior)[1:2], c("mu", "sigma2"))) cat("Prior: normal\n\n")
    if(identical(names(x$prior)[1], "support")) cat("Prior: uniform\n\n")       
  }
} 
if(x$nlm==FALSE & x$glm==FALSE){
cat("User-defined model & utility \n")}
#cat("\n")
cat("Number of repetitions = ",length(x$start.d),"\n",sep="")
cat("\n")
cat("Number of runs = ",dim(x$d)[1],"\n",sep="")
cat("\n")
cat("Number of factors = ",dim(x$d)[2],"\n",sep="")
cat("\n")
cat("Number of Phase I iterations = ",x$N1,"\n",sep="")
cat("\n")
cat("Number of Phase II iterations = ",x$N2,"\n",sep="")
cat("\n")
cat("Computer time = ",paste(hrs,":",mins,":",secs,sep=""),"\n",sep="")

}

summary.pace<-function(object,...){

print.pace(x=object)}

plot.pace<-function(x,...){

if(length(x$phase1.trace)>0 & length(x$phase2.trace)>0){
ulim<-max(c(x$phase1.trace,x$phase2.trace))
llim<-min(c(x$phase1.trace,x$phase2.trace))
plot(0:(length(x$phase1.trace)-1),x$phase1.trace,xlab="Phase I iteration",ylab="Observation of expected utility",ylim=c(llim,ulim),type="l",xlim=c(0,length(x$phase1.trace)))
new_z<-axTicks(side=1)
new_x<-(1:length(x$phase2.trace))*((length(x$phase1.trace)-1)/length(x$phase2.trace))
lines(new_x,x$phase2.trace,col=8)
legend(x="bottomright",legend=c("Phase I","Phase II"),col=c(1,8),lty=c(1,1),bty="n")
axis(side=3,labels=new_z/((length(x$phase1.trace)-1)/length(x$phase2.trace)),at=new_z)
mtext("Phase II iteration", side=3, line = par("mgp")[1]) }

if(length(x$phase1.trace)>0 & length(x$phase2.trace)==0){
plot(0:(length(x$phase1.trace)-1),x$phase1.trace,type="l",xlab="Phase I iteration",ylab="Observation of expected utility")
legend(x="bottomright",legend=c("Phase I"),col=1,lty=1,bty="n") }

if(length(x$phase1.trace)==0 & length(x$phase2.trace)>0){
plot(1:length(x$phase2.trace),x$phase2.trace,type="l",xlab="Phase II iteration",ylab="Observation of expected utility",col=8)
legend(x="bottomright",legend=c("Phase II"),col=8,lty=1,bty="n") }

}




