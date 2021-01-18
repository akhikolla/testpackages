context("get_self() etc.")

set.seed(20171204)
n <- 5
m <- 7
d <- matrix(runif(n*m, 0, 10), nrow=n, ncol=m)
dimnames(d) <- list(LETTERS[1:n], sample(LETTERS[3:(m+2)]))


test_that("get_self(), etc. work", {

    expect_equal(get_self(d),
                 c(A=NA, B=NA,
                   C=d["C","C"], D=d["D", "D"], E=d["E", "E"],
                   F=NA, I=NA, G=NA, H=NA))


    # get_best() and # which_best()
    expect_equal(get_best(d, "row"),
                 c(A=d["A","F"], B=d["B","C"], C=d["C","F"], D=d["D","E"], E=d["E","G"],
                   F=NA, I=NA, G=NA, H=NA))

    expect_equal(which_best(d, "row"),
                 c(A="F", B="C", C="F", D="E", E="G", F=NA, I=NA, G=NA, H=NA))

    expect_equal(get_best(d, "row", FALSE),
                 c(A=d["A","H"], B=d["B","D"], C=d["C","G"], D=d["D","G"], E=d["E","D"],
                   F=NA, I=NA, G=NA, H=NA))

    expect_equal(which_best(d, "row", FALSE),
                 c(A="H", B="D", C="G", D="G", E="D", F=NA, I=NA, G=NA, H=NA))

    expect_equal(get_best(d, "col"),
                 c(A=NA, B=NA, C=d["B","C"], D=d["C","D"], E=d["D","E"], F=d["C","F"],
                   I=d["B","I"], G=d["B","G"], H=d["C","H"]))

    expect_equal(which_best(d, "col"),
                 c(A=NA, B=NA, C="B", D="C", E="D", F="C", I="B", G="B", H="C"))

    expect_equal(get_best(d, "col", FALSE),
                 c(A=NA, B=NA, C=d["A","C"], D=d["E","D"], E=d["E","E"],
                   F=d["D","F"], I=d["A","I"], G=d["A","G"], H=d["A","H"]))

    expect_equal(which_best(d, "col", FALSE),
                 c(A=NA, B=NA, C="A", D="E", E="E", F="D", I="A", G="A", H="A"))

    # get_2ndbest() and # which_2ndbest()
    expect_equal(get_2ndbest(d, "row"),
                 c(A=d["A","E"], B=d["B","F"], C=d["C","C"], D=d["D","C"], E=d["E","C"],
                   F=NA, I=NA, G=NA, H=NA))

    expect_equal(which_2ndbest(d, "row"),
                 c(A="E", B="F", C="C", D="C", E="C", F=NA, I=NA, G=NA, H=NA))

    expect_equal(get_2ndbest(d, "row", FALSE),
                 c(A=d["A","G"], B=d["B","I"], C=d["C","E"], D=d["D","D"], E=d["E","E"],
                   F=NA, I=NA, G=NA, H=NA))

    expect_equal(which_2ndbest(d, "row", FALSE),
                 c(A="G", B="I", C="E", D="D", E="E", F=NA, I=NA, G=NA, H=NA))

    expect_equal(get_2ndbest(d, "col"),
                 c(A=NA, B=NA, C=d["C","C"], D=d["A","D"], E=d["B","E"],
                   F=d["B","F"], I=d["C","I"], G=d["E","G"], H=d["B","H"]))

    expect_equal(which_2ndbest(d, "col"),
                 c(A=NA, B=NA, C="C", D="A", E="B", F="B", I="C", G="E", H="B"))

    expect_equal(get_2ndbest(d, "col", FALSE),
                 c(A=NA, B=NA, C=d["E","C"], D=d["B","D"], E=d["C","E"],
                   F=d["E","F"], I=d["E","I"], G=d["D","G"], H=d["D","H"]))

    expect_equal(which_2ndbest(d, "col", FALSE),
                 c(A=NA, B=NA, C="E", D="B", E="C", F="E", I="E", G="D", H="D"))


    expected <- d
    for(i in LETTERS[3:5]) expected[i,i] <- NA
    expect_equal(get_nonself(d), expected)

})
