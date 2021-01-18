#' Dermatology Database
#'
#' We have removed 8 observations with missing values.
#' Data contains 12 clinical attributes and 21 histopathological attributes. The age
#' attribute has been discretized. The class variable has six levels; each
#' describing a skin disease.
#' 
#' @docType data
#' 
#' @references \url{https://archive.ics.uci.edu/ml/datasets/dermatology}
"derma"


## #' A data frame with observations corresponding to handwritten digits
## #'
## #' The original data contained 32x32 bitmaps for each digit. Cell values are 1 if the grayscale color is more
## #' towards black and zero otherwise. These bitmaps were then converted into blocks of 4x4 non-overlapping blocks
## #' in which all the cell values were added (values ranging from 0 to 16).
## #' The resulting matrix for an observations was then of dimension 8x8 (with each cell having values between 0 and 16).
## #'
## #' In the final version of the data in molic, we have a data frame (of characters)
## #' where the last column is the true class of the handwritten digit. All the remaining cells in the first 64 columns was
## #' converted to values in {a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p} corresponding to the 16 possibilities. This is needed
## #' in the outlier model. The to_single_chars() function can be exploited to achieve this for other data sets.
## #' 
## #' @docType data
## #' 
## #' @references \href{https://archive.ics.uci.edu/ml/datasets/optical+recognition+of+handwritten+digits}{Optical Recognition of Handwritten Digits}
## "digits"


#' A data frame with genetic data from the 1000 genomes project
#'
#' The data consists of \code{2504} DNA profiles, each genotyped on 304 SNPs (binary variables).
#' The data frame has \code{5008} rows, since each profile has two copies.
#' 
#' @docType data
#' 
#' @references \href{https://www.ncbi.nlm.nih.gov/pubmed/26432245}{1000 Genomes Project}
"tgp_dat"


#' A named list of character vectors.
#'
#' Every element in the list is a character vector that forms a haplotype
#' from the 1000 genomes project. If the list is unlisted, it should correspond
#' to colnames of tgp_dat. In other words, tgp_haps is a "haplotype-grouping"
#' of the variables in tgp_dat.
#' 
#' @docType data
#' 
#' @references \href{https://www.ncbi.nlm.nih.gov/pubmed/26432245}{1000 Genomes Project}
"tgp_haps"
