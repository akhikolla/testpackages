#'Plot of a Dirichlet process mixture of skew normal distribution partition
#'
#'@param z data matrix \code{d x n} with \code{d} dimensions in rows
#'and \code{n} observations in columns.
#'
#'@param c allocation vector of length \code{n} indicating which observation belongs to which
#'clusters.
#'
#'@param alpha current value of the DP concentration parameter.
#'
#'@param U_SS a list containing \code{"xi"}, \code{"psi"}, \code{"S"}, and \code{"df"}.
#'
#'@param i current MCMC iteration number.
#'
#'@param dims2plot index vector, subset of \code{1:d} indicating which dimensions should be drawn.
#'Default is all of them.
#'
#'@param ellipses a logical flag indicating whether ellipses should be drawn around clusters. Default
#'is \code{TRUE} if only 2 dimensions are plotted, \code{FALSE} otherwise.
#'
#'@param gg.add
#'A list of instructions to add to the \code{ggplot2} instruction (see \code{\link[ggplot2]{gg-add}}).
#'Default is \code{list(theme())}, which adds nothing to the plot.
#'
#'@param nbsim_dens number of simulated points used for computing clusters density contours in 2D
#'plots. Default is \code{1000} points.
#'
#'@author Boris Hejblum
#'
#'@import ellipse
#'
#'@import reshape2
#'
#'@importFrom stats dnorm pnorm rnorm
#'
#'@export

plot_DPMsn <- function(z, c, i="", alpha="?", U_SS,
                       dims2plot=1:nrow(z),
                       ellipses=ifelse(length(dims2plot)<3,TRUE,FALSE),
                       gg.add=list(theme()), nbsim_dens=1000){

    mean_sn01 <- (stats::dnorm(0)-stats::dnorm(Inf))/(stats::pnorm(Inf)-stats::pnorm(0))

    z <- z[dims2plot,]

    n <- ncol(z)
    p <- nrow(z)
    m <- numeric(n) # number of observations in each cluster
    m[unique(c)] <- table(c)[as.character(unique(c))]

    fullCl <- which(m!=0)

    U_xi2plot=sapply(U_SS, "[[", "xi")
    U_psi2plot=sapply(U_SS, "[[", "psi")
    U_Sigma2plot=lapply(U_SS, "[[", "S")
    U_SS2plot <- U_SS
    U_mu2plot <- U_xi2plot + U_psi2plot*mean_sn01
    rownames(U_mu2plot) <- rownames(z)
    zClusters <- factor(c, levels=as.character(fullCl), ordered=TRUE)

    expK <- ifelse(is.numeric(alpha), round(alpha*(digamma(alpha+n)-digamma(alpha))), NA)
    alpha2print <- ifelse(is.numeric(alpha), formatC(alpha, digits=2), alpha)

    if(p>2){
        zDplot <- melt(cbind.data.frame("ID"=as.character(1:n),
                                        t(z),
                                        "Cluster"=zClusters
        ),
        id.vars=c("ID", "Cluster"),
        variable.name = "dimensionX",
        value.name="X"
        )
        zDplotfull <- zDplot
        zDplotfull$Y <- zDplot$X
        zDplotfull$dimensionY <- zDplot$dimensionX

        lev <- as.character(1:length(levels(zDplot$dimensionX)))
        for(l in 2:length(lev)){
            move <- which(as.numeric(zDplot$dimensionX)<l)
            zDplottemp <- rbind.data.frame(zDplot[-move,], zDplot[move,])
            zDplottemp$Y <- zDplot$X
            zDplottemp$dimensionY <- zDplot$dimensionX
            zDplotfull <- rbind.data.frame(
                zDplotfull, zDplottemp)
        }

        UDplot <- melt(cbind.data.frame(t(U_mu2plot),
                                        "Cluster"=factor(as.character(fullCl),
                                                         levels=as.character(fullCl),
                                                         ordered=TRUE)
        ),
        id.vars=c("Cluster"),
        variable.name = "dimensionX",
        value.name="X"
        )
        UDplotfull <- UDplot
        UDplotfull$Y <- UDplotfull$X
        UDplotfull$dimensionY <- UDplotfull$dimensionX

        lev <- levels(UDplotfull$dimensionX)
        for(l in 2:length(lev)){
            move <- which(as.numeric(UDplotfull$dimensionX)<l)
            UDplottemp <- rbind.data.frame(UDplotfull[-move,], UDplotfull[move,])
            UDplottemp$Y <- UDplotfull$X
            UDplottemp$dimensionY <- UDplotfull$dimensionX
            UDplotfull <- rbind.data.frame(
                UDplotfull, UDplottemp)
        }

        #         ellipse95 <- data.frame()
        #         f