#' @title SPstats
#' @description Calculates the deviance information criterion (DIC)
#'      and Log-likelihood for fitted model outputs of pooled,
#' exchangeable, and spatial Split Population survival models for which a log-likelihood can be obtained using the formula
#' \emph{DIC = -2 * (L - P)}, where \emph{L} is the log likelihood of the data given the posterior means of the parameter and
#' \emph{P} is the  estimate of the effective number of parameters in the model.
#' @param object An object of the output of pooled, exchangeable, or spatial Split Population survival model .
#
#' @return List.
#
#' @export

SPstats = function(object){

    #Calculate L
    X  <- as.matrix(object$spstats$X)
    Z  <- as.matrix(object$spstats$Z)
    Y  <- as.matrix(object$spstats$Y)
    Y0 <- as.matrix(object$spstats$Y0)
    C  <- as.matrix(object$spstats$C)
    form <- object$spstats$form
    if(length(class(object)) == 2){

        data <- as.data.frame(cbind(Y, Y0, C, X, Z))

        theta_post = cbind(object$gammas, object$betas, object$rho)
        theta_hat = apply(theta_post, 2, mean)
        L = rllFun(theta_hat,Y, Y0, C, X, Z, data, form)$llik
        #Calculate P
        S = nrow(theta_post) #S = number of iterations
        #Add up the log likelihoods of each iteration
        llSum = 0
        sum1 = 0
        sum2 = 0
        sum3 = 0
        #l <- as.matrix(NA, nrow=S, ncol=1)
        for (s in 1:S) {
            theta_s = as.matrix(theta_post[s,])
            ll <- rllFun(theta_s, Y, Y0, C, X, Z, data, form)
            llSum <- llSum + ll$llik
            sum1 <- sum1 + ll$one
            sum2 <- sum2 + ll$two
            sum3 <- sum3 + ll$three
            #l[s,] <- llFun(theta_s,Y,Y0,C,X,Z,data)
        }
        P = 2 * (L - (1 / S * llSum))
        #Calculate DIC
        DIC <- -2 * (L - P)
        all <- sum1/S
        finite <- sum2/S
        small <- sum3/S
        #Return the results
        list(DIC = DIC, Loglik = L)

    } else {

        S <- object$spstats$S
        SP <- as.matrix(S)
        data <- as.data.frame(cbind(Y, Y0, C, X, Z, S))

        W <- matrix(NA, ncol = 1, nrow = nrow(X))
        V <- matrix(NA, ncol = 1, nrow = nrow(X))

        for(i in 1:length(SP)){
            sid <- SP[i]
            W[i] <- mean(object$W[,sid])
            V[i] <- mean(object$V[,sid])
        }

        data <- as.data.frame(cbind(Y, Y0, C, X, Z))
        theta_post = cbind(object$gammas, object$betas, object$rho)
        theta_hat = apply(theta_post, 2, mean)
        L = llFun(theta_hat,Y, Y0, C, X, Z, W, V, data, form)$llik
        #Calculate P
        S = nrow(theta_post) #S = number of iterations
        #Add up the log likelihoods of each iteration
        llSum = 0
        sum1 = 0
        sum2 = 0
        sum3 = 0
        #l <- as.matrix(NA, nrow=S, ncol=1)
        for (s in 1:S) {
            theta_s = as.matrix(theta_post[s,])
            ll <- llFun(theta_s, Y, Y0, C, X, Z, W, V, data, form)
            llSum <- llSum + ll$llik
            sum1 <- sum1 + ll$one
            sum2 <- sum2 + ll$two
            sum3 <- sum3 + ll$three
            #l[s,] <- llFun(theta_s,Y,Y0,C,X,Z,data)
        }
        P = 2 * (L - ((1 / S )* llSum))
        #Calculate DIC
        DIC <- -2 * (L - P)
        all <- sum1/S
        finite <- sum2/S
        small <- sum3/S
        #Return the results
        list(DIC = DIC, Loglik = L)}
}









