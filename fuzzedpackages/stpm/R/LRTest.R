#'Likelihood-ratio test
#'@param LA Log-likelihood for alternative hyphotesis
#'@param L0 Log-likelihood for null hyphotesis
#'@param df Degrees of freedom for Chi-square test
#'@return p-value of LR test.
#'
LRTest <- function(LA, L0, df=1) {
    d.log.lik.chisq <- 2*(LA - L0)
    return(1-pchisq(d.log.lik.chisq, df = df))
}