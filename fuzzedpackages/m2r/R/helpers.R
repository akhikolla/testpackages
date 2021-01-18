is.mac <- function() grepl("darwin", R.version$platform)
is.win <- function() .Platform$OS.type == "windows"
is.linux <- function() (.Platform$OS.type == "unix") && (is.mac() == FALSE)
is.unix  <- function() .Platform$OS.type == "unix"
is.solaris <- function() grepl("solaris", R.version$os)


file.path2 <- function(...){
	sep <- if(is.unix()) "/" else "\\"
	paste0(list(...), collapse = sep)
}

write_to_temp <- function(contents) {
  tempname <- tempfile()
  writeLines(contents, tempname)
  tempname
}

name_and_increment <- function(prefix, option) {
  # grab # of current rings, set ring number, and increment
  num <- getOption(option)
  num <- ifelse(is.null(num), 1L, strtoi(num))
  setOption(option, num + 1)

  # make name
  sprintf("m2rint%s%08d", prefix, num)
}




listify <- function (strings) paste0("{", paste(strings, collapse = ","), "}")

listify_mat <- function (mat) listify(apply(mat, 1, listify))
# (mat <- matrix(1:9, nrow = 3))
# listify_mat(mat)




delistify <- function (string, f, collapse, ...) {

  # dewhiten
  s <- str_replace_all(string, "[\\s]+", "")

  # kill outside {}
  s <- str_sub(s, 3, -3)

  # break it up
  s <- str_split(s, fixed("},{"))[[1]]
  s <- str_split(s, fixed(","))

  # f, if available
  if ( !missing(f) ) s <- lapply(s, f, ...)

  # collapse, if available
  if ( !missing(collapse) ) s <- do.call(collapse, s)

  # return
  s

}
# string <- "{{1,2,3}, {1,34,45}, {2213,1123,6543}, {0,0,0}}"
# delistify(string)
# str(delistify(string))
# delistify(string, as.integer)
# delistify(string, as.integer, rbind)






mpolyList_to_m2_str <- function(mpolyList) {

  # allow for character vectors
  if (inherits(mpolyList, "mpoly")) mpolyList <- mpolyList(mpolyList)

  # parse it if it's a numeric or character string
  if(is.numeric(mpolyList)) mpolyList <- as.character(mpolyList)
  if(is.character(mpolyList)) mpolyList <- mp(mpolyList)

  # convert mpolylist to strings readable by m2
  vec <- vapply(mpolyList, print, character(1),
    silent = TRUE, stars = TRUE
  )
  vec <- str_replace_all(vec, "[\\s]+", "")
  vec <- str_replace_all(vec, "\\*\\*", "^")

  # return
  vec
}
# mpolyList_to_m2_str( mp(c("x^3","x + y z", "1")) )
# mpolyList_to_m2_str( c(1, 2, 3) )






`%notin%` <- function(x, y) { !(x %in% y) }
# 2 %notin% 1:5
# 2 %notin% (1:5)[-2]





`%:%` <- function(x, y) {
  # do type checking
  stopifnot(
    is.character(x), length(x) == 1,
    is.character(y), length(y) == 1
  )

  # check letters
  if(all(c(x,y) %in% letters)) { # lower case letters given
    x_ndx <- which(x == letters)
    y_ndx <- which(y == letters)
    letters[x_ndx:y_ndx]
  } else if (all(c(x,y) %in% LETTERS)) { # upper case letters given
    x_ndx <- which(x == LETTERS)
    y_ndx <- which(y == LETTERS)
    LETTERS[x_ndx:y_ndx]
  } else {
    stop("letters must be of the same case")
  }

}


# copy pryr::dots so that we don't have to import pryr
dots <- function(...) {
  eval(substitute(alist(...)))
}



