#' Basis for all bivrec plot functions
#'
#' @param parameters pass from plot fcts
#' @param ctimes pass from plot fcts
#' @param nsubject pass from plot fcts
#' @param temp pass from plot fcts
#' @param args pass from plot fcts
#' @param c optional
#' @param cm optional
#'
#' @keywords internal
#' @noRd

basicplot <- function(parameters, ctimes, nsubject, temp, args, c, cm, byp) {

  if (missing(c)) {c = 1}
  if (missing(cm)) {cm = 2}

  #### Reformat data to start-stop times ########
  for (iter in 1:nsubject) {
    subject <- parameters[parameters$id==iter,] #pull dat for each subject
    #print(subject)
    if (nrow(subject)==1) {
      start_time <- c(0, subject$xij)
      stop_time <- c(subject$xij, subject$xij+subject$yij)
      temp <- rbind(temp, cbind(sub_id=rep(iter, 2), sub_epi=rep(1, 2),
                                start_time, stop_time, sub_ci=rep(subject$ci, 2)))
    } else {
      subject <- subject[order(subject$episode),]
      times <- c(rbind(subject$xij, subject$yij))
      start_time <- 0
      stop_time <- times[1]
      for (j in 2:length(times)) {
        start_time[j] <- stop_time[j-1]
        stop_time[j] <- start_time[j]+times[j]
      }
      sub_id <- rep(iter, length(times))
      sub_ci <- rep(subject$ci[1], length(times))
      sub_epi <- rep(subject$episode, each=2)
      temp <- rbind(temp, cbind(sub_id, sub_epi, start_time, stop_time, sub_ci))
    }
  }

  data4plot <- data.frame(temp)
  data4plot <- data4plot[-which(data4plot$start_time==data4plot$stop_time), ]
  colnames(data4plot) <- c("id", "epi", "start_time", "stop_time", "ci")
  data4plot  <- data4plot[order(data4plot$ci),]

  ###### PLOT ########
  # get the range for the x and y axis
  xrange <- c(0, round(max(ctimes), digits = -1) + 10) #Time
  yrange <- c(0, nsubject) #subjects

  # set up the plot
  plot.new()
  legendtext = c(args[4], args[5])
  plot.window(xlim=xrange, ylim=yrange)
  axis(side=1)
  title(main=args[1], xlab=args[2], ylab=args[3], cex.main = cm)
  if (byp == FALSE) {
    legend("bottomright", legend=legendtext, col = c("blue", "red"),
         lty = 1, cex=c, inset=.02, bg = "white", bty='n', seg.len = 0.7,
         x.intersp=0.9, y.intersp = 0.85, xjust = 0, yjust = 0)
  }

  # add line segments
  newid = 0
  for (c_iter in sort(ctimes)) {
    newid = newid + 1
    subject <- subset(data4plot, data4plot$ci == c_iter)
    subject$newid = newid
    if (nrow(subject)==1) {
      segments(subject$start_time[1], subject$newid,
               subject$stop_time[1], subject$newid, col="blue")
    } else {
      for (j in 1:length(subject$id)) {
        if (j %% 2 == 1) colors <- "blue" else colors <- "red"
        segments(subject$start_time[j], subject$newid,
                 subject$stop_time[j], subject$newid, col=colors)
      }
    }
  }

}
