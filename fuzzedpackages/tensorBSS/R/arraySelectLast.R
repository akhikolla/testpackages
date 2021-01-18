arraySelectLast <- function(x, take)
{
  DIMS <- dim(x)
  nDIMS <- length(DIMS)
  COMMAS <- paste(rep(",", nDIMS - 1), collapse = "")
  eval(parse(text = paste("x[", COMMAS, "take]", collapse = "")))
}