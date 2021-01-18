library(Countr)
data(football)
POISSON <- FALSE

home_poiss <- glm(formula = homeTeamGoals ~ 1, family = poisson, data = football)
home_wei <- renewalCount(formula = homeTeamGoals ~ 1, data = football,
                         dist = "weibull", weiMethod = "conv_dePril",
                         computeHessian = FALSE, 
                         control = renewal.control(trace = 0,
                                                   method = "nlminb")
                         )

print(rbind(
    c(logLik(home_poiss), logLik(home_wei)),
    c(AIC(home_poiss), AIC(home_wei))
    )
)

away_poiss <- glm(formula = awayTeamGoals ~ 1, family = poisson, data = football)
away_wei <- renewalCount(formula = awayTeamGoals ~ 1, data = football,
                         dist = "weibull", weiMethod = "conv_dePril",
                         computeHessian = FALSE, 
                         control = renewal.control(trace = 0,
                                                   method = "nlminb")
                         )

print(rbind(
    c(logLik(away_poiss), logLik(away_wei)),
    c(AIC(away_poiss), AIC(away_wei))
)
)

## we can add a likelihood ratio test using the library(lmtest)
library(lmtest)
print(lrtest(away_poiss, away_wei))

breaks_ <- 0:5
pears <- compareToGLM(poisson_model = away_poiss,
                      breaks = breaks_, weibull = away_wei)

frequency_plot(pears$Counts, pears$Actual,
               dplyr::select(pears, contains("_predicted"))
               )

## run the formal chi-sq test gof
test_wei <- chiSq_gof(away_wei, breaks = breaks_)
test_pois <- chiSq_gof(away_poiss, breaks = breaks_)


## ---- Old code: not useful anymore
## ================================================================================
## ----------------------- formal chi-square g-o-f test for the Poisson model -----
## ================================================================================
if (POISSON) {
    ## compute the same for Poisson
    library(pscl)
    ppoiss <- predprob(away_poiss)
    ppoissave <- colMeans(ppoiss)
    relfreqs <- table(football$awayTeamGoals) / nrow(football)
    pears2 <- data.frame(Counts = c(0:4, ">= 5"),
                         Actual = c(relfreqs[1:5], sum(relfreqs[6:length(relfreqs)])),
                         poisson = c(ppoissave[1:5], sum(ppoissave[6:length(relfreqs)]))
                         )
    
    ## now estimate the poisson model using the score argument
    ## .obj_fit <- function(theta) {
    ##     th <- c(theta, 0)
    ##     names(th) <- names(coef(away_wei))
    ##     -sum(away_wei$score_fct(th))
    ## }
    
    ## init <- coef(away_wei)[-length(coef(away_wei))]
    ## res <- nlminb(start = init, objective = .obj_fit)
    ## print(res$par - coef(away_poiss))
    
    .obj <- function(theta) {
        th <- c(theta, 0)
        names(th) <- names(coef(away_wei))
        away_wei$score_fct(th)
    }
    
    library(numDeriv)
    si <- jacobian(.obj, coef(away_poiss))
    colnames(si) <- paste0("si_", names(coef(away_poiss)))
    ppoiss2 <- cbind(ppoiss[, 0:5], rowSums(ppoiss[, -(0:5)]))
    
    .ff <- function(ind) {
        res <- rep(0, ncol(ppoiss2))
        res[ind] <- 1
        res
    }
    
    cl <- ifelse(football$awayTeamGoals > 4, 6, football$awayTeamGoals + 1)
    dij <- t(sapply(cl, .ff))
    dij_pij <- dij - ppoiss2
    colnames(dij_pij) <- paste0("di_", as.character(pears2$Counts))
    
    mat <- data.frame(y = 1, si, dij_pij[, 1:5])
    
    ## compute the X2 statistic
    mod <- lm(y~. -1, data = mat)
    stat <- nrow(mat) * summary(mod)$r.squared
    
    rval <- matrix(NA, nrow = 1, ncol = 3)
    colnames(rval) <- c("DF", "Chisq", "Pr(>Chisq)")
    
    rval[, 1] <- ncol(dij_pij) - 1 ## number of cells - 1
    rval[, 2] <- stat
    rval[, 3] <- pchisq(stat, rval[, 1], lower.tail = FALSE)
    
    title <- "chi-square goodness-of-fit test\n"
    topnote <- paste("Cells considered ", paste(colnames(dij_pij), collapse = " "),
                     sep = "", collapse = "\n")
    
    test_pois <- structure(as.data.frame(rval), heading = c(title, topnote),
                           class = c("anova", "data.frame"))
}
