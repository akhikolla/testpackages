plotCounts <- function(data, paired = TRUE, points = TRUE, 
                       points.method = "jitter", xlabel = "", ylabel = "Faecal egg counts [epg]", 
                       ...){
  epgsL <- reshape(data, direction="long", varying = list(names(data)))
  epgsL$time <- factor(epgsL$time, levels=1:2, labels=c("before treatment","after treatment"))
  if (paired){
  xyplot(epgsL[,2] ~ time, group=epgsL$id, data=epgsL, type=c("p","l"), col=1, xlab=xlabel, ylab=ylabel, ...)
  } else {
    boxplot(epgsL[,2] ~ time, data=epgsL ,xlab=xlabel, ylab=ylabel, ...)
    if (points){
    stripchart(epgsL[,2] ~ time, vertical = TRUE, data = epgsL, 
               method = points.method, add = TRUE, pch = 20, col = 'blue')
    }
  }
}