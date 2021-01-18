context("coef_gen")

# set up input variables
edge_list <- vector("list", length = 5)
edge_list[[1]] <- integer(0)
edge_list[[2]] <- integer(0)
edge_list[[3]] <- 1
edge_list[[4]] <- c(1, 3)
edge_list[[5]] <- c(1, 2)
edge_list <- sparsebnUtils::as.edgeList(edge_list)
names(edge_list) <- c("V1", "V2", "V3", "V4", "V5")
nlevels <- c(3, 5, 2, 2, 3)

# test
test_that("coef_gen run as expected", {
  ### coef_gen runs with default setting
  expect_error(coef_gen(edge_list = edge_list), NA)

  ### coef_gen runs with manual setting
  expect_error(coef_gen(edge_list = edge_list, n_levels = nlevels, FUN = runif, flip = FALSE), NA)
})

test_that("Check input edge_list", {
  edge_list_list <- as.list(edge_list)
  expect_error(coef_gen(edge_list = edge_list_list))

})

test_that("Check input n_levels", {
  n_levels_wrongLevel <- nlevels
  n_levels_wrongLevel[1] <- 1
  expect_error(coef_gen(edge_list = edge_list, n_levels = n_levels_wrongLevel))

})

test_that("coef_gen output", {
  out <- coef_gen(edge_list = edge_list, n_levels = nlevels)
  expect_equal(is.list(out), TRUE)
  expect_equal(length(out), 5)
  for (i in 1:5) {
    n_parents <- length(edge_list[[i]])
    if (n_parents) {
      col_length <- sum(nlevels[edge_list[[i]]]-1)+1
      expect_equal(dim(out[[i]]), c(nlevels[i], col_length))
    }
    else
      expect_equal(out[[i]], NULL)
  }
})
