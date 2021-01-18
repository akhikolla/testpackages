#'Returns string w/o leading whitespace
#'@param x a string to trim
trim.leading <- function (x)  sub("^\\s+", "", x)

#'Returns string w/o trailing whitespace
#'@param x a string to trim
trim.trailing <- function (x) sub("\\s+$", "", x)

#'Returns string w/o leading or trailing whitespace
#'@param x a string to trim
trim <- function (x) gsub("^\\s+|\\s+$", "", x)