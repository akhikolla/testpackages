print.tladle <- function(x, ...){
  cat("The tensorial Ladle based on", x$method, "gives for the data set", x$data.name, "the dimension estimates:\n")
  print(sapply(x$ResMode, function(x) x$k))
}