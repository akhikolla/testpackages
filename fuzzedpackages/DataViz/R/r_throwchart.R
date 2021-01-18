## Code donne par Christophe, edite par Timothy
## Dans le cas de la version statique une serie de courbes quadratiques est envoyee a l'utilisateur

quad <- function(point,before_point,after_point,offSet){
  polyAux <- function(point){return((point - before_point)*(point - after_point))}
  result <- (point - before_point)*(point - after_point) / polyAux((before_point+after_point)/2) * (after_point-before_point-offSet)
  return(result)
}
r_throwchart <- function(before, after, xlim, ylim, col, lwd, offSet){
  ## If the xlims are equal then xmin and xmax are the lims
  if((xlim[1] == xlim[2]))
  {
    max_before = max(before)
    max_after = max(after)
    min_before = min(before)
    min_after = min(after)
    if(min_before <= min_after)
    {
      definite_min = min_before
    }
    else
    {
      definite_min = min_after
    }
    if(max_before >= max_after)
    {
      definite_max = max_before
    }
    else
    {
      definite_max = max_after
    }
    xlim <- c(definite_min, definite_max)
  }
  ## If the ylims are the equal then they are set as maximum value
  if(ylim[1] == ylim[2]) 
  {
    max_difference = max(abs(after-before))
    ylim <- c(-max_difference, max_difference)
  }
  if(offSet > 0)
  {
    ylim <- c(offSet, max_difference)
  }
  if(offSet < 0)
  {
    ylim <- c(-max_difference, offSet)
  }
  plot(x=xlim,y=c(offSet,offSet),xlim=xlim,ylim=ylim,type="l", ann = FALSE)
  ## since values were transformed into tibbles in the following functions to iterate through the before values is equivalent to before[[1]][i]
  for(i in 1:NROW(before)){
    points <- seq(before[[1]][i],after[[1]][i],length.out = 100)
    lines(points,quad(points,before[[1]][[i]],after[[1]][[i]],offSet)+offSet,xlim=before[[1]][[i]],after[[1]][[i]],type="l",col="blue",lwd=2.5)  
  }
}
