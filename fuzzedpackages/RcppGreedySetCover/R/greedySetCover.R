#' Greedy Set Cover
#'
#' Fast greedy set cover algorithm.

#' @param X Two-column data.frame in long format: Column 1 identifies the sets, column 2 the elements.
#' @param data.table If \code{TRUE} returns a \code{data.table} with keys given by sets and elements. 
#' If FALSE returns a \code{data.frame}, sorted by sets and elements. 
#' @keywords greedy set cover
#' @return If \code{data.table == TRUE} a \code{data.table}, keyed by sets and elements. 
#' Else a \code{data.frame}, sorted by sets and elements. 
#' Column names are derived from input.
#' @examples
#' # Create some data.
#' set.seed(333)
#' X <- data.table::rbindlist(
#'   lapply(
#'     seq_len(1e4L),
#'     function(x) list(element=sample.int(n=1e3L,size=sample.int(50L,1L)))
#'   ),
#'   idcol="set"
#' )
#' # Elements are integers 1,2,...,1000.
#' 
#' # Run set cover
#' res <- greedySetCover(X,FALSE)
#' head(res)
#' 
#' # Check if all elements are covered.
#' identical(sort(unique(res$element)),sort(unique(X$element)))


greedySetCover <- function(X,data.table=TRUE) {
  
  stopifnot(ncol(X)==2L)
  # Input: Two column dataframe. First column represents the sets,
  # Second column represents the elements in the sets.
  X <- setDT(copy(X)) # This might be improved in the future by looking at 
  # shallow copy.
  names_ <- copy(names(X))
  
  n1 <- names_[1L]
  n2 <- names_[2L]
  
  X[,"i0":=.GRP-1L,by=n1] 
  X[,"i1":=.GRP-1L,by=n2]
  
  # Get group sizes for pre allocation of space in greedy_set_cover.
  # This data.tables can convientently by reused as lookup table later.
  
  ex_text0 <- sprintf(".(.N,'orvar'=%s[1L] )  ",n1)
  ex_text1 <- sprintf(".(.N,'orvar'=%s[1L] )  ",n2)
  
  Group_size_i0 <- X[,eval(parse(text=ex_text0)),keyby="i0"]
  Group_size_i1 <- X[,eval(parse(text=ex_text1)),keyby="i1"]
  
  Res_set_cover <- greedy_set_cover2(
    X[["i0"]],X[["i1"]],Group_size_i0[["N"]],Group_size_i1[["N"]]
  )
  
  # Recover original cols.
  # Index using i0 / i1 + offset to translate
  
  Out <- data.table(
    var1=Group_size_i0[[3L]][Res_set_cover[,1L] + 1L],
    var2=Group_size_i1[[3L]][Res_set_cover[,2L] + 1L]
    
  )
  
  setnames(Out,names_)
  setkey(Out)
  
  if(!data.table) setDF(Out)
  Out  
  
}