################################################################################

ERROR_NLINE <- "Error when getting the number of lines of the file."
ERROR_FILE  <- "Error with file. Lines not equal? No last empty line?"
ERROR_DIM   <- "Error when guessing the dimensions of the file."

################################################################################

CODE_012    <- rep(NA_integer_, 256); CODE_012[49:51]    <- 0:2
CODE_DIGITS <- rep(NA_integer_, 256); CODE_DIGITS[49:58] <- 0:9

################################################################################

is_int <- function(x) floor(x) == x

################################################################################

stop2 <- function(msg) stop(msg, call. = FALSE) 
warning2 <- function(msg) warning(msg, call. = FALSE) 

################################################################################