## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(lazyarray)

## -----------------------------------------------------------------------------
# Sample data (~24 MB)
x <- rnorm(3e6); dim(x) <- c(10, 100, 100, 30)

# Save array to a path
path <- tempfile()
arr <- create_lazyarray(path, 'double', dim(x), multipart = TRUE)
arr[] <- x

## -----------------------------------------------------------------------------
# Load existing array
arr <- load_lazyarray(path)

## -----------------------------------------------------------------------------
# Make loaded array writable
arr$make_writable()
arr$can_write

## -----------------------------------------------------------------------------
arr$make_writable()
dimnames(arr) <- list(
  A = 1:10,
  B = 1:100,
  C = 1:100,
  D = 1:30
)

## -----------------------------------------------------------------------------
# Subset/read array
y1 <- arr[]              
y2 <- arr[,,,3]          

# Write to slice of data, writing to slices along the 
# last dimension is optimized
arr[,,,1] <- seq_len(1e5)

## -----------------------------------------------------------------------------
sub <- subset(arr, A ~ A <= 2, B ~ B == 10)
dim(sub)

## -----------------------------------------------------------------------------
arr$remove_data()

