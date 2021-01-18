# test_simulate.R
# Created by: Sajal Kumar
#
# Modified by : Ruby Sharma
# Date : February 27 2017
#
# Modified by : Ruby Sharma
# Date : December 2 2018

# Modified Update : Dr. Joe Song and Ruby Sharma
# Date            : May 26 2020
#                 : Added new test cases to test marginal
#                 : distribution

library(testthat)
library(FunChisq)

context("Testing simulate_tables()")

#Attributes to test
#1) y = f(x)
#2) Not a constant function - All populated samples are not supposed to be on the same column.

Test_Functional_table = function(iter)
{

  func.flag = FALSE

  for(i in seq(iter))
  {
    Get.Stats = Construct_Table("functional")

    conti.table = Get.Stats$conti.table
    noise.table = Get.Stats$noise.table

    # check if all samples are populated

    check.all.samples = All.sample.check(conti.table,Get.Stats$sample.size)

    if(check.all.samples$flag)
    {
      func.flag = TRUE
      failure.summary = check.all.samples$failure.table
      break
    }

    # Check for non zero cells

    check.non.zero = No.Non.Zero.Check(conti.table)

    if(check.non.zero$flag)
    {
      func.flag = TRUE
      failure.summary = check.non.zero$failure.table
      break
    }

    # Check y = f(x)

    check.functional = Functional.check(conti.table)

    if(check.functional$flag)
    {
      func.flag = TRUE
      failure.summary = check.functional$failure.table
      break
    }

    # check constant functions

    check.constant = Constant.check(conti.table)

    if(check.constant$flag)
    {
      func.flag = TRUE
      failure.summary = check.constant$failure.table
      break
    }

    # check.margin = margin.check(conti.table, noise.table)
    # if(check.margin$flag)
    # {
    #  func.flag = TRUE
    #  failure.summary = check.margin$failure.table
    #  break
    # }
  }

  #if any table was flagged, the test failed.
  expect_identical(func.flag, FALSE)
}


#Attributes to test
#1) y = f(x)
#2) Not a constant function - All populated samples are not supposed to be on the same column.
Test_Functional_Discontinuous_table = function(iter)
{

  func.flag = FALSE

  for(i in seq(iter))
  {
    Get.Stats = Construct_Table("discontinuous")

    conti.table = Get.Stats$conti.table

    # check if all samples are populated

    check.all.samples = All.sample.check(conti.table,Get.Stats$sample.size)

    if(check.all.samples$flag)
    {
      func.flag = TRUE
      failure.summary = check.all.samples$failure.table
      break
    }

    # Check for non zero cells

    check.non.zero = No.Non.Zero.Check(conti.table)

    if(check.non.zero$flag)
    {
      func.flag = TRUE
      failure.summary = check.non.zero$failure.table
      break
    }

    # Check y = f(x)

    check.functional = Functional.check(conti.table)

    if(check.functional$flag)
    {
      func.flag = TRUE
      failure.summary = check.functional$failure.table
      break
    }

    # check constant functions

    check.constant = Constant.check(conti.table)

    if(check.functional$flag)
    {
      func.flag = TRUE
      failure.summary = check.constant$failure.table
      break
    }
  }

  #if any table was flagged, the test failed.
  expect_identical(func.flag, FALSE)
}

#Attributes to test
#1) y = f(x)
#2) x ! = f(y)
#2) Not a constant function - All populated samples are not supposed to be on the same column.

Test_Functional_Many_to_one_table = function(iter)
{

  non.mono.func.flag = FALSE

  for(i in seq(iter))
  {

    Get.Stats = Construct_Table("many.to.one")

    conti.table = Get.Stats$conti.table

    # check if all samples are populated

    check.all.samples = All.sample.check(conti.table,Get.Stats$sample.size)

    if(check.all.samples$flag)
    {

      non.mono.func.flag = TRUE
      failure.summary = check.all.samples$failure.table
      break
    }

    # Check for non zero cells

    check.non.zero = No.Non.Zero.Check(conti.table)

    if(check.non.zero$flag)
    {
      non.mono.func.flag = TRUE
      failure.summary = check.non.zero$failure.table
      break
    }

    # Check y = f(x)

    check.functional = Functional.check(conti.table)

    if(check.functional$flag)
    {
      non.mono.func.flag = TRUE
      failure.summary = check.functional$failure.table
      break
    }

    # check x != f(y)

    check.non.functional.invert = Non.functional.check(t(conti.table))

    if(check.non.functional.invert$flag)
    {
      non.mono.func.flag = TRUE
      failure.summary = conti.table
      break
    }

    # check constant functions

    check.constant = Constant.check(conti.table)

    if(check.constant$flag)
    {
      non.mono.func.flag = TRUE
      failure.summary = check.constant$failure.table
      break
    }
  }

  #if any table was flagged, the test failed.

  expect_identical(non.mono.func.flag, FALSE)
}


#Attributes to test
#1) y ! = f(x)

Test_Non_Functional_table = function(iter)
{

  non.func.flag = FALSE

  for(i in seq(iter))
  {

    Get.Stats = Construct_Table("dependent.non.functional")

    conti.table = Get.Stats$conti.table

    # check if all samples are populated

    check.all.samples = All.sample.check(conti.table,Get.Stats$sample.size)

    if(check.all.samples$flag)
    {

      non.func.flag = TRUE
      failure.summary = check.all.samples$failure.table
      break
    }

    # Check for non zero cells

    check.non.zero = No.Non.Zero.Check(conti.table)

    if(check.non.zero$flag)
    {
      non.func.flag = TRUE
      failure.summary = check.non.zero$failure.table
      break
    }

    # check y != f(x)

    check.non.functional = Non.functional.check(conti.table)

    if(check.non.functional$flag)
    {
      non.func.flag = TRUE
      failure.summary = conti.table
      break
    }
  }

  #if any table was flagged, the test failed.
  expect_identical(non.func.flag, FALSE)
}

#Attributes to test
#1) No non zero cell
#2) No dependency
#3) y != f(x) AND x != f(y)

Test_Independent_table = function(iter)
{

  ind.flag = FALSE

  for(i in seq(iter))
  {

    Get.Stats = Construct_Table("independent")

    conti.table = Get.Stats$conti.table

    # check if all samples are populated

    check.all.samples = All.sample.check(conti.table,Get.Stats$sample.size)

    if(check.all.samples$flag)
    {
      ind.flag = TRUE
      failure.summary = check.all.samples$failure.table
      break
    }
  }

  #if any table was flagged, the test failed.
  expect_identical(ind.flag, FALSE)
}


#Setup random parameters and construct a table using simulate.tables
Construct_Table = function(type)
{
  row_col.equality = sample(c(TRUE,FALSE),1)
  if(row_col.equality) {
    nrows = sample(c(2:15),1)
    ncols = nrows
  } else {
    nrows = sample(c(2:15),1)
    ncols = sample(c(2:15),1)
  }

  if(type=="many.to.one" && nrows == 2)
    nrows = 3

  sample.size = ((sample(c(1:100),1) * sample(c(1:100),1))) + nrows
  row.marginal = rep(0,nrows)

  if(type=="independent" && sample.size < (nrows * ncols))
    sample.size = nrows * ncols

  if(type=="dependent.non.functional" && sample.size < (nrows * ncols))
    sample.size = nrows * ncols

  row.marginal.set = sample(c(TRUE,FALSE),1)
  if(row.marginal.set){
    row.marginal = runif(nrows)
    row.marginal = row.marginal/sum(row.marginal)
    conti.table = simulate_tables(n=sample.size, nrow = nrows, ncol = ncols, type = type, n.tables = 1, row.marginal = row.marginal, col.marginal = NULL, noise = 0.2)
  } else {
    conti.table = simulate_tables(n=sample.size, nrow = nrows, ncol = ncols, type = type, n.tables = 1, noise =0.2)
  }

  list(conti.table = conti.table$sample.list[[1]], nrows = nrows, ncols = ncols, sample.size = sample.size, row.marginal = row.marginal, noise.table = conti.table$noise.list[[1]])
}

# The table should contain all the samples specified in sample size
All.sample.check = function(conti.table, sample.size)
{
  flag=FALSE
  failure.summary = NULL

  if(sum(conti.table)!=sample.size){
    flag = TRUE
    failure.summary = conti.table
  }

  list(flag=flag,failure.table = failure.summary)
}

# If there are no non-zero elements in the table then something's wrong.
No.Non.Zero.Check = function(conti.table)
{
  index = which(conti.table!=0, arr.ind=TRUE)

  flag=FALSE
  failure.summary = NULL

  if(length(index)==0){
    flag = TRUE
    failure.summary = conti.table
  }

  list(flag=flag, failure.table = failure.summary)
}

# If there are any zero element in the table then something's wrong.
All.Non.Zero.Check = function(conti.table)
{
  index = which(conti.table==0, arr.ind=TRUE)

  flag=FALSE
  failure.summary = NULL

  if(length(index)!=0){
    flag = TRUE
    failure.summary = conti.table
  }

  list(flag=flag, failure.table = failure.summary)
}

# Testing y = f(x)
#If there exists a row such that it contains samples for more than one y then the table is invalid
Functional.check = function(conti.table)
{
  index = which(conti.table!=0, arr.ind=TRUE)

  flag=FALSE
  failure.summary = NULL

  if(length(index[,1])!=length(unique(index[,1]))){
    flag = TRUE
    failure.summary = conti.table
  }

  list(flag=flag, failure.table = failure.summary)
}

# Testing y != f(X)
Non.functional.check = function(conti.table)
{

  status = Functional.check(conti.table)

  flag=FALSE
  failure.summary = NULL

  if(!status$flag)
    flag = TRUE

  list(flag=flag, failure.table = failure.summary)
}

# If there exists a column (or y) such that it contains all the samples then the table is invalid
Constant.check = function(conti.table)
{
  index = which(conti.table!=0, arr.ind=TRUE)

  flag=FALSE
  failure.summary = NULL

  if(length(unique(index[,2]))==1){
    flag=TRUE
    failure.summary = conti.table
  }

  list(flag=flag, failure.table = failure.summary)
}


#If there exists a zero cell after removing all zero rows and columns then the table is invalid
Dependency.check = function(conti.table)
{
  index.zero = which(conti.table==0, arr.ind=TRUE)
  index.non.zero = which(conti.table!=0, arr.ind=TRUE)

  flag=FALSE
  failure.summary = NULL

  if(length(index.zero) == 0)
  {
    flag = TRUE
    failure.summary = conti.table

  } else {

    index.zero.copy = index.zero
    for(j in 1:length(index.zero.copy[,1]))
    {

      if(length(which(index.zero[,1]==index.zero.copy[j,1]))==ncol(conti.table)){
        index.zero = index.zero[-which(index.zero[,1]==index.zero.copy[j,1]),,drop=FALSE]
      }

      if(length(which(index.zero[,2]==index.zero.copy[j,2]))==nrow(conti.table)){
        index.zero = index.zero[-which(index.zero[,2]==index.zero.copy[j,2]),,drop=FALSE]
      }
    }

    if(nrow(index.zero)==0){
      flag = TRUE
      failure.summary = conti.table
    }
  }

  list(flag=flag, failure.table = failure.summary)
}

margin.check = function(conti.table, noise.table)
{
  flag=FALSE
  failure.summary = NULL

  if(any(rowSums(conti.table)!= rowSums(noise.table)))
    {
      flag = TRUE
      failure.summary = noise.table
    }
  list(flag=flag, failure.table = failure.summary)

}

test_that("Testing simulate_tables()", {
  Test_Functional_table(2)
  Test_Independent_table(2)
  Test_Non_Functional_table(2)
  Test_Functional_Many_to_one_table(2)
  Test_Functional_Discontinuous_table(2)
})

test_that("Testing marginals in simulate_tables()", {

  col.mar <- c(4, 2, 15)
  col.mar <- col.mar / sum(col.mar)
  row.mar <- c(15, 1, 1)
  row.mar <- row.mar / sum(row.mar)

  # Testing for row marginals for functional tables
  tables <- simulate_tables(
    100, type="fun", n.tables = 10,
    row.marginal=row.mar
  )

  p.values <- vector("numeric", 10)
  for(i in 1:10) {
    p.values[i] <- chisq.test(rowSums(tables$sample.list[[i]]), p=row.mar)$p.value
  }
  expect_equal(sum(p.values > 0.05) > 5, TRUE)



  # Testing for column marginals, this would fail everytime
  # because of zero column, so only comparing probabities where the observed
  # colsum is non-zero.

  tables <- simulate_tables(
    100, type="fun", n.tables = 10,
    col.marginal = col.mar
  )

  for(i in 1:10) {
    sumcols = colSums(tables$sample.list[[i]])
    nonzero.sum = which(sumcols!=0)
    sumcols = sumcols[nonzero.sum]
    scaled.mar = col.mar[nonzero.sum]/sum(col.mar[nonzero.sum])
    p.values[i] <- chisq.test(sumcols, p=scaled.mar)$p.value
  }
  expect_equal(sum(p.values > 0.05) > 5, TRUE)


  # Testing for indepdendent tables, both row and column row marginal are begin tested
  col.mar <- c(10, 20, 15)
  col.mar <- col.mar / sum(col.mar)
  row.mar <- c(15, 1, 30)
  row.mar <- row.mar / sum(row.mar)

  tables <- simulate_tables(
    100, type="in", n.tables = 10,
    row.marginal=row.mar, col.marginal = col.mar
  )

  p.values <- vector("numeric", 10)
  for(i in 1:10) {
    p.values[i] <- suppressWarnings(
      chisq.test(
        rowSums(tables$sample.list[[i]]), p=row.mar
      )$p.value
    )
  }
  expect_equal(sum(p.values > 0.05) > 5, TRUE)

  for(i in 1:10) {
    p.values[i] <- chisq.test(colSums(tables$sample.list[[i]]), p=col.mar)$p.value
  }
  expect_equal(sum(p.values > 0.05) > 5, TRUE)


  # Testing for many.to.one tables, row marginal test
  row.mar <- c(15, 1, 30)
  row.mar <- row.mar / sum(row.mar)
  col.mar <- c(20, 40, 30)
  col.mar <- col.mar / sum(col.mar)

  tables <- simulate_tables(
    100, type="many.to.one", n.tables = 10,
    row.marginal=row.mar
  )
  p.values <- vector("numeric", 10)
  for(i in 1:10) {
    p.values[i] <- suppressWarnings(
      chisq.test(
        rowSums(tables$sample.list[[i]]), p=row.mar
      )$p.value
    )
  }
  expect_equal(sum(p.values > 0.05) > 5, TRUE)


  # Testing for many.to.one tables, column marginal test
  tables <- simulate_tables(
    100, type="many.to.one", n.tables = 10,
    col.marginal = col.mar
  )
  for(i in 1:10) {
    sumcols = colSums(tables$sample.list[[i]])
    nonzero.sum = which(sumcols!=0)
    sumcols = sumcols[ nonzero.sum]
    scaled.mar = col.mar[nonzero.sum]/sum(col.mar[nonzero.sum])
    p.values[i] <- chisq.test(sumcols, p=scaled.mar)$p.value
  }
  expect_equal(sum(p.values > 0.05) > 5, TRUE)

  # Testing for discontinuous tables, row marginal test
  row.mar <- c(10, 25, 30)
  row.mar <- row.mar / sum(row.mar)
  col.mar <- c(20, 40, 30)
  col.mar <- col.mar / sum(col.mar)

  tables <- simulate_tables(
    100, type="discontinuous", n.tables = 10,
    row.marginal=row.mar
  )
  p.values <- vector("numeric", 10)
  for(i in 1:10) {
    p.values[i] <- chisq.test(rowSums(tables$sample.list[[i]]), p=row.mar)$p.value
  }
  expect_equal(sum(p.values > 0.05) > 5, TRUE)

  # Testing for discontinuous tables, column marginal test
  tables <- simulate_tables(
    100, type="discontinuous", n.tables = 10,
    col.marginal = col.mar
  )
  for(i in 1:10) {
    sumcols = colSums(tables$sample.list[[i]])
    nonzero.sum = which(sumcols!=0)
    sumcols = sumcols[ nonzero.sum]
    scaled.mar = col.mar[nonzero.sum]/sum(col.mar[ nonzero.sum])
    p.values[i] <- chisq.test(sumcols, p=scaled.mar)$p.value
  }
  expect_equal(sum(p.values > 0.05) > 5, TRUE)


  # Testing for dependent.non.functional tables, row marginal test
  row.mar <- c(10, 25, 30)
  row.mar <- row.mar / sum(row.mar)
  col.mar <- c(20, 40, 30)
  col.mar <- col.mar / sum(col.mar)

  tables <- simulate_tables(
    100, type="dep", n.tables = 10,
    row.marginal=row.mar
  )
  p.values <- vector("numeric", 10)
  for(i in 1:10) {
    p.values[i] <- chisq.test(rowSums(tables$sample.list[[i]]), p=row.mar)$p.value
  }
  expect_equal(sum(p.values > 0.05) > 5, TRUE)

})
