context("test-alphapart")

test_that("Test input for AlphaPart ped", {

  ## Example - not really interesting but set up in a way that we can test behaviour of AlphaPart()
  ped <- data.frame(
        id=c( "A",   "B",   "C",   "F",   "G",   "J",   "K",   "E",   "D",   "I",   "H",   "L",   "M"),
       fid=c(  NA,    NA,    NA,   "A",   "C",   "F",   "F",   "A",   "A",   "A",   "C",   "K",    NA),
       mid=c(  NA,    NA,    NA,   "B",   "B",   "G",   "G",   "B",   "B",   "B",   "B",   "H",   "H"),
        by=c(   0,     0,     0,     1,     1,     2,     2,     1,     1,     1,     2,     3,     3) + 2000,
      dum1=13:1,
       pat=c( "A",   "A",   "B",   "A",   "A",   "A",   "A",   "A",   "B",   "C",    NA,   "B",   "B"),
      trt1=c(0.56,  0.04, -0.60,  0.47, -0.31,  0.03,     0,  1.00,  1.00,  1.00,  0.10,  0.50,  0.00),
      dum2=(1:13)+6,
      trt2=c(0.10,  0.24, -0.30,  0.17, -0.21,  0.13,     0,     0,     0,     0,  0.20,     0,  0.00))

  ## ... to test unknown argument
  ped2 <- ped
  ped2$fid <- as.character(ped2$fid)
  ped2[is.na(ped2$fid), "fid"] <- "0"
  ped2$fid <- as.factor(ped2$fid)

  ## ... to test recode argument
  ped3 <- ped[order(orderPed(ped=ped[, c("id", "fid", "mid")])), ]
  ped3$idI  <- 1:nrow(ped3)
  ped3$fidI <- match(ped3$fid, ped3$id)
  ped3$midI <- match(ped3$mid, ped3$id)

  ## ... to test recode and unknown argument
  ped4 <- ped[order(orderPed(ped=ped[, c("id", "fid", "mid")])), ]
  ped4$idI  <- 1:nrow(ped4)
  ped4$fidI <- match(ped4$fid, ped3$id, nomatch=99)
  ped4$midI <- match(ped4$mid, ped3$id, nomatch=99)

  ## .. to test that colBV columns are not numeric
  ped5 <- ped
  ped5$trt1 <- as.character(ped5$trt1)

  ## --- Run ---

  ## Error when path column contains NA
  expect_error(AlphaPart(x=ped[, c("id", "fid", "mid", "pat", "trt1", "trt2")]))

  ## Error if recode=FALSE and input is not numeric
  expect_error(AlphaPart(x=ped[, c("id", "fid", "mid", "pat", "trt1", "trt2")], pathNA=TRUE, recode=FALSE, verbose=0))

  ## Error if coBV columns are not numeric
  expect_error(AlphaPart(x=ped[, c("id", "fid", "mid", "pat", "trt1", "trt2")], pathNA=TRUE, recode=FALSE, verbose=0))

  ## Error if NAs are present
  pedX <- ped
  pedX[1, "trt1"] <- NA
  expect_error(AlphaPart(x=pedX[, c("id", "fid", "mid", "pat", "trt1", "trt2")], pathNA=TRUE, verbose=0))
})

test_that("Test the output of AlphaPart function", {
  ped <- data.frame(
        id=c( "A",   "B",   "C",   "F",   "G",   "J",   "K",   "E",   "D",   "I",   "H",   "L",   "M"),
       fid=c(  NA,    NA,    NA,   "A",   "C",   "F",   "F",   "A",   "A",   "A",   "C",   "K",    NA),
       mid=c(  NA,    NA,    NA,   "B",   "B",   "G",   "G",   "B",   "B",   "B",   "B",   "H",   "H"),
        by=c(   0,     0,     0,     1,     1,     2,     2,     1,     1,     1,     2,     3,     3) + 2000,
      dum1=13:1,
       pat=c( "A",   "A",   "B",   "A",   "A",   "A",   "A",   "A",   "B",   "C",    NA,   "B",   "B"),
      trt1=c(0.56,  0.04, -0.60,  0.47, -0.31,  0.03,     0,  1.00,  1.00,  1.00,  0.10,  0.50,  0.00),
      dum2=(1:13)+6,
      trt2=c(0.10,  0.24, -0.30,  0.17, -0.21,  0.13,     0,     0,     0,     0,  0.20,     0,  0.00))

  ## ... to test unknown argument
  ped2 <- ped
  ped2$fid <- as.character(ped2$fid)
  ped2[is.na(ped2$fid), "fid"] <- "0"
  ped2$fid <- as.factor(ped2$fid)

  ## ... to test recode argument
  ped3 <- ped[order(orderPed(ped=ped[, c("id", "fid", "mid")])), ]
  ped3$idI  <- 1:nrow(ped3)
  ped3$fidI <- match(ped3$fid, ped3$id)
  ped3$midI <- match(ped3$mid, ped3$id)
    ## ... to test recode and unknown argument

  ret   <- AlphaPart(x=ped[, c("id", "fid", "mid", "pat", "trt1", "trt2")],  pathNA=TRUE, verbose=0)
  ret2  <- AlphaPart(x=ped,                                                  pathNA=TRUE, verbose=0, colId=1,     colFid=2,      colMid=3,      colPath=6,     colBV=c(7, 9))
  ret3  <- AlphaPart(x=ped,                                                  pathNA=TRUE, verbose=0, colId="id",  colFid="fid",  colMid="mid",  colPath="pat", colBV=c("trt1", "trt2"))
  ped$idI <- ped$id
  ret3a <- AlphaPart(x=ped,                                                  pathNA=TRUE, verbose=0, colId="idI", colFid="fid",  colMid="mid",  colPath="pat", colBV=c("trt1", "trt2"))
  ped$idI <- NULL
  ret4  <- AlphaPart(x=ped[, c("id", "fid", "mid", "pat", "trt1", "trt2")],  pathNA=TRUE, verbose=0, colId="id",  colFid="fid",  colMid="mid",  colPath="pat", colBV=c("trt1", "trt2"))

  ## ... to test recode argument
  ret5  <- AlphaPart(x=ped3,                                                 pathNA=TRUE, verbose=0, colId="idI", colFid="fidI", colMid="midI", colPath="pat", colBV=c("trt1", "trt2"))
  ret6  <- AlphaPart(x=ped3,                                                 pathNA=TRUE, verbose=0, colId="idI", colFid="fidI", colMid="midI", colPath="pat", colBV=c("trt1", "trt2"), recode=FALSE)
  ## ... to test recode and unknown argument
  ret7  <- AlphaPart(x=ped3,                                                 pathNA=TRUE, verbose=0, colId="idI", colFid="fidI", colMid="midI", colPath="pat", colBV=c("trt1", "trt2"), recode=FALSE, unknown=99)


  ## --- Overall result ---

  ## List component names
  expect_equal(names(ret), c(ret$info$lT, "info"))

  ## The same input (so should be the output), but different column definitions
  expect_equal(ret2, ret3)
  ret3a$trt1$"idI" <- NULL
  expect_equal(ret3$trt1, ret3a$trt1)
  expect_equal(ret,  ret4)

  ## Use of argument recode and unknown
  expect_equal(ret5, ret6)
  expect_equal(ret6, ret7)

  ## Use of recode=FALSE
  tmp <- ret6$trt1[order(ret6$trt1$id), ]
  tmp$idI <- tmp$fidI <- tmp$midI <- NULL
  expect_equal(tmp, ret3$trt1[order(ret3$trt1$id), ])

  ## --- Check the meta info component ---

  expect_equal(as.character(ret$info$path), "pat", checkNames=FALSE)
  expect_equal(ret$info$nP, 4)
  expect_equal(ret$info$lP, c("A", "B", "C", "XXX"))
  expect_equal(ret$info$nT, 2)
  expect_equal(ret$info$lT, c("trt1", "trt2"))

  ## --- Check partition components ---

  ## Make sure we have all individuals for all traits
  expect_equal(nrow(ret$trt1), 13)
  expect_equal(nrow(ret$trt2), 13)

  ## Test column names
  expect_equal(colnames(ret$trt1),  c("id", "fid", "mid", "pat", "trt1", "trt1_pa", "trt1_w", "trt1_A", "trt1_B", "trt1_C", "trt1_XXX"))
})

test_that("Test computation", {
    ped <- data.frame(
        id=c( "A",   "B",   "C",   "F",   "G",   "J",   "K",   "E",   "D",   "I",   "H",   "L",   "M"),
       fid=c(  NA,    NA,    NA,   "A",   "C",   "F",   "F",   "A",   "A",   "A",   "C",   "K",    NA),
       mid=c(  NA,    NA,    NA,   "B",   "B",   "G",   "G",   "B",   "B",   "B",   "B",   "H",   "H"),
        by=c(   0,     0,     0,     1,     1,     2,     2,     1,     1,     1,     2,     3,     3) + 2000,
      dum1=13:1,
       pat=c( "A",   "A",   "B",   "A",   "A",   "A",   "A",   "A",   "B",   "C",    NA,   "B",   "B"),
      trt1=c(0.56,  0.04, -0.60,  0.47, -0.31,  0.03,     0,  1.00,  1.00,  1.00,  0.10,  0.50,  0.00),
      dum2=(1:13)+6,
      trt2=c(0.10,  0.24, -0.30,  0.17, -0.21,  0.13,     0,     0,     0,     0,  0.20,     0,  0.00))


  ret   <- AlphaPart(x=ped[, c("id", "fid", "mid", "pat", "trt1", "trt2")],  pathNA=TRUE, verbose=0)

  ## --- Check computations ---

  ## ret   <- AlphaPart(x=ped[, c("id", "fid", "mid", "pat", "trt1", "trt2")], pathNA=TRUE, verbose=2)
  ## Gene flow (T)
  ##  [1,] 1.000 .    .     .    .    . .   . . . .   . .
  ##  [2,] .     1.00 .     .    .    . .   . . . .   . .
  ##  [3,] .     .    1.000 .    .    . .   . . . .   . .
  ##  [4,] 0.500 0.50 .     1.00 .    . .   . . . .   . .
  ##  [5,] .     0.50 0.500 .    1.00 . .   . . . .   . .
  ##  [6,] 0.250 0.50 0.250 0.50 0.50 1 .   . . . .   . .
  ##  [7,] 0.250 0.50 0.250 0.50 0.50 . 1.0 . . . .   . .
  ##  [8,] 0.500 0.50 .     .    .    . .   1 . . .   . .
  ##  [9,] 0.500 0.50 .     .    .    . .   . 1 . .   . .
  ## [10,] 0.500 0.50 .     .    .    . .   . . 1 .   . .
  ## [11,] .     0.50 0.500 .    .    . .   . . . 1.0 . .
  ## [12,] 0.125 0.50 0.375 0.25 0.25 . 0.5 . . . 0.5 1 .
  ## [13,] .     0.25 0.250 .    .    . .   . . . 0.5 . 1
  ## Gene flow inverse (inv(T))
  ## 1   1.0  .    .    .    .   .  .   . . .  .   . .
  ## 2   .    1.0  .    .    .   .  .   . . .  .   . .
  ## 3   .    .    1.0  .    .   .  .   . . .  .   . .
  ## 4  -0.5 -0.5  .    1.0  .   .  .   . . .  .   . .
  ## 5   .   -0.5 -0.5  .    1.0 .  .   . . .  .   . .
  ## 6   .    .    .   -0.5 -0.5 1  .   . . .  .   . .
  ## 7   .    .    .   -0.5 -0.5 .  1.0 . . .  .   . .
  ## 8  -0.5 -0.5  .    .    .   .  .   1 . .  .   . .
  ## 9  -0.5 -0.5  .    .    .   .  .   . 1 .  .   . .
  ## 10 -0.5 -0.5  .    .    .   .  .   . . 1  .   . .
  ## 11  .   -0.5 -0.5  .    .   .  .   . . .  1.0 . .
  ## 12  .    .    .    .    .   . -0.5 . . . -0.5 1 .
  ## 13  .    .    .    .    .   .  .   . . . -0.5 . 1
  ## Selected gene flow (TP) for path 1
  ##  [1,] 1.000 .    . .    .    . .   . . . . . .
  ##  [2,] .     1.00 . .    .    . .   . . . . . .
  ##  [3,] .     .    0 .    .    . .   . . . . . .
  ##  [4,] 0.500 0.50 . 1.00 .    . .   . . . . . .
  ##  [5,] .     0.50 0 .    1.00 . .   . . . . . .
  ##  [6,] 0.250 0.50 0 0.50 0.50 1 .   . . . . . .
  ##  [7,] 0.250 0.50 0 0.50 0.50 . 1.0 . . . . . .
  ##  [8,] 0.500 0.50 . .    .    . .   1 . . . . .
  ##  [9,] 0.500 0.50 . .    .    . .   . 0 . . . .
  ## [10,] 0.500 0.50 . .    .    . .   . . 0 . . .
  ## [11,] .     0.50 0 .    .    . .   . . . 0 . .
  ## [12,] 0.125 0.50 0 0.25 0.25 . 0.5 . . . 0 0 .
  ## [13,] .     0.25 0 .    .    . .   . . . 0 . 0
  ## Selected gene flow (TP) for path 2
  ##  [1,] 0 . .     . . . . . . . . . .
  ##  [2,] . 0 .     . . . . . . . . . .
  ##  [3,] . . 1.000 . . . . . . . . . .
  ##  [4,] 0 0 .     0 . . . . . . . . .
  ##  [5,] . 0 0.500 . 0 . . . . . . . .
  ##  [6,] 0 0 0.250 0 0 0 . . . . . . .
  ##  [7,] 0 0 0.250 0 0 . 0 . . . . . .
  ##  [8,] 0 0 .     . . . . 0 . . . . .
  ##  [9,] 0 0 .     . . . . . 1 . . . .
  ## [10,] 0 0 .     . . . . . . 0 . . .
  ## [11,] . 0 0.500 . . . . . . . 0 . .
  ## [12,] 0 0 0.375 0 0 . 0 . . . 0 1 .
  ## [13,] . 0 0.250 . . . . . . . 0 . 1
  ## Selected gene flow (TP) for path 3
  ##  [1,] 0 . . . . . . . . . . . .
  ##  [2,] . 0 . . . . . . . . . . .
  ##  [3,] . . 0 . . . . . . . . . .
  ##  [4,] 0 0 . 0 . . . . . . . . .
  ##  [5,] . 0 0 . 0 . . . . . . . .
  ##  [6,] 0 0 0 0 0 0 . . . . . . .
  ##  [7,] 0 0 0 0 0 . 0 . . . . . .
  ##  [8,] 0 0 . . . . . 0 . . . . .
  ##  [9,] 0 0 . . . . . . 0 . . . .
  ## [10,] 0 0 . . . . . . . 1 . . .
  ## [11,] . 0 0 . . . . . . . 0 . .
  ## [12,] 0 0 0 0 0 . 0 . . . 0 0 .
  ## [13,] . 0 0 . . . . . . . 0 . 0
  ## Selected gene flow (TP) for path 4
  ##  [1,] 0 . . . . . . . . . .   . .
  ##  [2,] . 0 . . . . . . . . .   . .
  ##  [3,] . . 0 . . . . . . . .   . .
  ##  [4,] 0 0 . 0 . . . . . . .   . .
  ##  [5,] . 0 0 . 0 . . . . . .   . .
  ##  [6,] 0 0 0 0 0 0 . . . . .   . .
  ##  [7,] 0 0 0 0 0 . 0 . . . .   . .
  ##  [8,] 0 0 . . . . . 0 . . .   . .
  ##  [9,] 0 0 . . . . . . 0 . .   . .
  ## [10,] 0 0 . . . . . . . 0 .   . .
  ## [11,] . 0 0 . . . . . . . 1.0 . .
  ## [12,] 0 0 0 0 0 . 0 . . . 0.5 0 .
  ## [13,] . 0 0 . . . . . . . 0.5 . 0

  ## ret$trt1 --> trait 1
  ##    id  fid  mid   path  trt1 trt1_pa trt1_w trt1_A trt1_B trt1_C trt1_XXX
  ## 1   A <NA> <NA>      A  0.56    0.00   0.56   0.56   0.00    0.0    0.000
  ## --> base animal from path 1 and all AGV is from path A
  expect_equal(as.vector(unlist(ret$trt1[1, -(1:4)])), c(0.56, 0, 0.56, 0.56, 0, 0, 0), checkNames=FALSE)

  ##    id  fid  mid   path  trt1 trt1_pa trt1_w trt1_A trt1_B trt1_C trt1_XXX
  ## 2   B <NA> <NA>      A  0.04    0.00   0.04   0.04   0.00    0.0    0.000
  ## --> ditto

  ##    id  fid  mid   path  trt1 trt1_pa trt1_w trt1_A trt1_B trt1_C trt1_XXX
  ## 3   C <NA> <NA>      B -0.60    0.00  -0.60   0.00  -0.60    0.0    0.000
  ## --> ditto for path B
  expect_equal(as.vector(unlist(ret$trt1[3, -(1:4)])), c(-0.6, 0, -0.6, 0, -0.6, 0, 0), checkNames=FALSE)

  ##    id  fid  mid   path  trt1 trt1_pa trt1_w trt1_A trt1_B trt1_C trt1_XXX
  ## 4   F    A    B      A  0.47    0.30   0.17   0.47   0.00    0.0    0.000
  ## --> ditto for path 1 and both parents are from path A
  expect_equal(as.vector(unlist(ret$trt1[4, -(1:4)])), c(0.47, 0.3, 0.17, 0.47, 0, 0, 0), checkNames=FALSE)

  ##    id  fid  mid   path  trt1 trt1_pa trt1_w trt1_A trt1_B trt1_C trt1_XXX
  ## 5   G    C    B      A -0.31   -0.28  -0.03  -0.01  -0.30    0.0    0.000
  ##
  ## id fid mid path
  ##  A   0   0    A 1
  ##  B   0   0    A 2
  ##  C   0   0    B 3
  ##  F   A   B    A 4
  ##  G   C   B    A 5
  ##
  ## a_._5 = 1/2(a_3  + a_2)  + w_5
  ##       = 1/2(w_3  + w_2)  + w_5 --> this corresponds to T[5,] .    0.5    0.50 .   1.0 . . . . .
  ##       = 1/2(-0.6 + 0.04) + -0.03
  ##
  ## a_1_5 = T_5 P_1_5 w
  ##     T[5,] .    0.5    0.50 .   1.0 . . . . .
  ##    TP[5,] .    0.5       . .   1.0 . . . . . --> selects genes from path 1 - parent 1 genes and Mendelian sampling (parent 3 is from path 2!!!)
  ##  TInv[5,] .   -0.5   -0.50 .   1.0 . . . . . --> when computing Mendelian sampling we substract parent average -1/2(a_f(i) + a_m(i))
  ##       = 1/2w_2 + w_5
  ##       = 1/2*0.04 + -0.03
  ##       = -0.01
  ##
  ## a_2_5 = T_5 P_2_5 w
  ##     T[5,] .    0.5    0.50 .   1.0 . . . . .
  ##    TP[5,] .    0.0    0.50 .     0 . . . . . --> selects genes from path 2 (parent 2 and animals is from path 1!!!)
  ##
  ##       = 1/2w_3
  ##       = 1/2*-0.60
  ##       = -0.30
  expect_equal(as.vector(unlist(ret$trt1[5, -(1:4)])), c(-0.31, -0.28, -0.03, -0.01, -0.3, 0, 0), checkNames=FALSE)

  ##    id  fid  mid   path  trt1 trt1_pa trt1_w trt1_A trt1_B trt1_C trt1_XXX
  ## 6   J    F    G      A  0.03    0.08  -0.05   0.18  -0.15    0.0    0.000
  ##
  ## id fid mid path
  ##  A   0   0    1 1
  ##  B   0   0    1 2
  ##  C   0   0    2 3
  ##  F   A   B    1 4
  ##  G   3   B    1 5
  ##  J   F   G    1 6
  ##
  ## a_._6 = 1/2(a_4  + a_5)  + w_6
  ##       = 1/2(1/2a_1 + 1/2a_2 + w_4) + 1/2(1/2a_3  + 1/2a_2 + w_5) + w_6
  ##       = 1/4w_1 + 1/4w_2 + 1/2w_4 + 1/4w_3 + 1/4w_2 + 1/2w_5 + w_6
  ##       = 1/4w_1 + 1/2w_2 + 1/4w_3 + 1/2w_4 + 1/2w_5 + w_6 --> this corresponds to T[6,] 0.25 0.5 0.25 0.5 0.5 1 . . . .
  ##       = 0.03
  ##
  ## a_1_6 = T_6 P_1_6 w
  ##              1   2    3   4   5 6
  ##     T[6,] 0.25 0.5 0.25 0.5 0.5 1 . . . .
  ##    TP[6,] 0.25 0.5 0    0.5 0.5 1 . . . . --> selects genes from path 1 - Mendelian sampling (6) and both parents (4,5) as well as three grandparents (1, 2x2)
  ##       = 1/4w_1 + 1/2w_2 + 0w_3 + 1/2w_4 + 1/2w_5 + w_6
  ##       = 0.18
  ## ...
  expect_equal(as.vector(unlist(ret$trt1[6, -(1:4)])), c(0.03, 0.08, -0.05, 0.18, -0.15, 0, 0), checkNames=FALSE)

  ##    id  fid  mid   path  trt1 trt1_pa trt1_w trt1_A trt1_B trt1_C trt1_XXX
  ## 7   K    F    G      A     0    0.08  -0.08   0.15  -0.15     0        0
  expect_equal(as.vector(unlist(ret$trt1[7, -(1:4)])), as.numeric(c(0, 0.08, -0.08, 0.15, -0.15, 0, 0)), checkNames=FALSE)

  ##    id  fid  mid   path  trt1 trt1_pa trt1_w trt1_A trt1_B trt1_C trt1_XXX
  ## 8   E    A    B      A  1.00    0.30   0.70   1.00   0.00    0.0    0.000
  ## --> reference for test animals, both parents path 1 and animal path 1 = all OK
  expect_equal(as.vector(unlist(ret$trt1[8, -(1:4)])), c(1, 0.3, 0.7, 1, 0, 0, 0), checkNames=FALSE)

  ##    id  fid  mid   path  trt1 trt1_pa trt1_w trt1_A trt1_B trt1_C trt1_XXX
  ## 9   D    A    B      B  1.00    0.30   0.70   0.30   0.70    0.0    0.000
  ## --> the same parents as for 8, but different partitioning, due to different claimed path - path 2 was doing "selection" of Mendelian sampling!!!
  expect_equal(as.vector(unlist(ret$trt1[9, -(1:4)])), c(1, 0.3, 0.7, 0.3, 0.7, 0, 0), checkNames=FALSE)

  ##    id  fid  mid   path  trt1 trt1_pa trt1_w trt1_A trt1_B trt1_C trt1_XXX
  ## 10  I    A    B      C  1.00    0.30   0.70   0.30   0.00    0.7    0.000
  ## --> ditto
  expect_equal(as.vector(unlist(ret$trt1[10, -(1:4)])), c(1, 0.3, 0.7, 0.3, 0, 0.7, 0), checkNames=FALSE)

  ##    id  fid  mid   path  trt1 trt1_pa trt1_w trt1_A trt1_B trt1_C trt1_XXX
  ## 11  H    C    B    XXX  0.10  -0.065  0.165  0.235  -0.30    0.0    0.165
  ## --> argument pathNA=TRUE sets unknown path to dummy path called XXX
  expect_equal(as.vector(unlist(ret$trt1[11, -(1:4)])), c(0.1, -0.28, 0.38, 0.02, -0.3, 0, 0.38), checkNames=FALSE)

  ##    id  fid  mid   path  trt1 trt1_pa trt1_w trt1_A trt1_B trt1_C trt1_XXX
  ## 12  L    K    H      B  0.50    0.05   0.45  0.085  0.225    0.0     0.19
  expect_equal(as.vector(unlist(ret$trt1[12, -(1:4)])), c(0.5, 0.05, 0.45, 0.085, 0.225, 0.0, 0.19), checkNames=FALSE)

  ##    id  fid  mid   path  trt1 trt1_pa trt1_w trt1_A trt1_B trt1_C trt1_XXX
  ## 13  M   NA    H      B  0.00    0.05  -0.05   0.01  -0.20      0     0.19
  ##
  ## id fid mid path
  ##  A   0   0    A 1
  ##  B   0   0    A 2
  ##  C   0   0    B 3
  ##  H   C   B  XXX 11
  ##  M   0   H    B 13
  ##
  ## a_._13 = 1/2(0 + a_11)                   + w_13
  ##        = 1/2(0 + 1/2a_3 + 1/2a_2 + w_11) + w_13 --> this corresponds to T[13,]
  ##        = 1/2*(0 + -0.3   + 0.02   + 0.38) + -0.05
  ##        = 1/2*0.1 + -0.05
  ##        = 0
  ##
  ## a_1_5 = T_5 P_1_5 w
  ##     T[13,] .     0.25 0.250 .    .    . .   . . .  0.5 . 1
  ##    TP[13,] .     0.25 0     .    .    . .   . . .  0   . 0 --> selects genes from path 1
  ##  TInv[13,] .     .    .     .    .    . .   . . . -0.5 . 1 --> when computing Mendelian sampling we substract parent average -1/2(a_f(i) + a_m(i))
  ##       = 1/4*w_2
  ##       = 1/2*0.04
  ##       = 0.01
  ##
  ## ...
  expect_equal(as.vector(unlist(ret$trt1[13, -(1:4)])), c(0, 0.05, -0.05, 0.01, -0.20, 0, 0.19), checkNames=FALSE)

})
