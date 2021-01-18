#' Loop through basicplot for each categorical variable
#'
#' @param pred_levels pass from plot fcts
#' @param plotdat pass from plot fcts
#' @param cov_name pass from plot fcts
#' @param args pass from plot fcts
#'
#' @import graphics
#' @keywords internal
#' @noRd
#'
ploteach <- function(pred_levels, plotdat, cov_name, args) {
  dfs = args_new = list()
  for (p in 1:length(pred_levels)) {
    ##EXTRACT VECTORS FOR PLOTTING FUNCTION
    parameters <- plotdat[which(plotdat[ ,6] == pred_levels[p]), 1:5]
    new_main = paste(cov_name, " = ", pred_levels[p], sep="")
    unik_ids <- unique(parameters$id)
    nsubject2 <- length(unik_ids)
    parameters$id2 <- rep(NA, nrow(parameters))
    for (i in 1:nsubject2){
      index <- which(parameters$id == unik_ids[i])
      parameters$id2[index]=rep(i, length(index))
    }
    parameters2 <- cbind(parameters[,6], parameters[,2:5])
    colnames(parameters2) <- c("id", "episode", "xij", "yij", "ci")
    dfs[[p]] <- data.frame(parameters2)
    args2 = args
    args2[1] = new_main
    args_new[[p]] <- args2
  }
  rdim <- ceiling(length(pred_levels)/2) + 1
  layoutvect <- c(seq(1, length(pred_levels)),
                  rep(length(pred_levels)+1, 2))
  layout(matrix(layoutvect, nrow=rdim, byrow=TRUE),
         heights = (c(rep(4,rdim-1), 1)))
  par(mar=c(5,4,4,2)+0.1)
  for (p_iter in 1:length(pred_levels)) {
    #draw p_iter plot in mfrow
    basicplot(parameters = dfs[[p_iter]], ctimes = unique(dfs[[p_iter]]$ci),
              nsubject=max(unique(dfs[[p_iter]]$id)), temp=NULL,
              args = args_new[[p_iter]], c=0.75, cm=0.9, byp=TRUE)
  }
  legendtext = c(args[4], args[5])
  #xlim2 = round(max(ctimes), digits = -1) + 10
  par(mar=c(1,1,1,1)+0.1)
  plot(0:1, 0:1, xaxt='n',yaxt='n',bty='n',ylab='',xlab='',
       col="white")
  legend("center", legend=legendtext, col = c("blue", "red"),
         lty = 1, bg = "white", bty='n', horiz=TRUE)
}

#' Plot by function
#'
#' @param df passed from plot.bivrecSurv
#' @param predictors passed from plot.bivrecSurv
#' @param args passed from plot.bivrecSurv
#'
#' @keywords internal
#' @noRd

plotBy <- function(df, args) {

  #number of levels for each predictor
  num_levs <- apply(df[, 6:ncol(df)], 2, function(x) length(unique(x)))
  to_delete <- which(num_levs > 6) + 5
  message1 <- paste(colnames(df)[to_delete], " not used - either continuous or had more than 6 levels.", sep="")
  print(message1)
  df <- df[, -to_delete]

  if (ncol(df)==5) {stop("Cannot break by covariate. All covariates are continuous or have more than 6 levels.")}

  cov_names <- colnames(df)[6:ncol(df)]
  nsubject <-length(unique(df$id))

  message <- paste("Subjects for plots: ", nsubject, ".", sep="")
  print(message)

  if (length(cov_names)==1) {
    pred_levels = unique(df[,6])
    plotdat = na.omit(df[ , 1:6])
    ploteach(pred_levels, plotdat, cov_name = cov_names, args)
  } else {
    for (k in 1:length(cov_names)) {
      pred_levels = unique(df[ ,5+k])
      plotdat = na.omit(df[, c(1:5, 5+k)])
      ploteach(pred_levels, plotdat, cov_name = cov_names[k], args)
    }
  }
  par(mfrow=c(1, 1))
}


