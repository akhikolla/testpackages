context("Correlation between two matrices")

test_that("corr_betw_matrices works", {

    set.seed(20171201)
    n <- 100
    p <- 5
    q <- 7
    x <- matrix(rnorm(p*n), ncol=p)
    y <- matrix(rnorm(q*n), ncol=q)

    # paired (use first 5 columns of y)
    result <- corr_betw_matrices(x, y[,1:p], "paired")
    expect_equivalent(result, diag(cor(cbind(x,y[,1:p]))[1:p,p+1:p]))

    # all
    full_result <- corr_betw_matrices(x, y, "all")
    expect_equivalent(full_result, cor(cbind(x,y))[1:p,p+1:q])

    # best right
    expected <- data.frame(corr=apply(full_result, 1, max),
                           yindex=apply(full_result, 1, which.max))
    expected$ycol <- paste0("V", 1:ncol(y))[expected$yindex]

    result <- corr_betw_matrices(x, y, "bestright")
    expect_equivalent(result, expected)

    # all best
    expect_equivalent(corr_betw_matrices(x, y, "bestpairs"),
                      data.frame(corr=numeric(0),
                                 xindex=integer(0),
                                 yindex=integer(0),
                                 xcol=character(0),
                                 ycol=character(0),
                                 stringsAsFactors=FALSE))

    expected <- data.frame(corr=full_result[full_result > 0.15],
                           xindex=row(full_result)[full_result > 0.15],
                           yindex=col(full_result)[full_result > 0.15])
    expected$xcol <- paste0("V", 1:ncol(x))[expected$xindex]
    expected$ycol <- paste0("V", 1:ncol(y))[expected$yindex]
    expected <- expected[order(expected$xindex, expected$yindex),]

    expect_equivalent(corr_betw_matrices(x, y, "bestpairs", corr_threshold=0.15),
                      expected)

})


test_that("corr_betw_matrices works with some missing values", {

    set.seed(20171201)
    n <- 100
    p <- 5
    q <- 7
    x <- matrix(rnorm(p*n), ncol=p)
    y <- matrix(rnorm(q*n), ncol=q)

    x[sample(1:prod(dim(x)), 5)] <- NA
    y[sample(1:prod(dim(y)), 6)] <- NA

    # paired (use first 5 columns of y)
    result <- corr_betw_matrices(x, y[,1:p], "paired")
    expect_equivalent(result, diag(cor(cbind(x,y[,1:p]), use="pairwise.complete.obs")[1:p,p+1:p]))

    # all
    full_result <- corr_betw_matrices(x, y, "all")
    expect_equivalent(full_result, cor(cbind(x,y), use="pairwise.complete.obs")[1:p,p+1:q])

    # best right
    expected <- data.frame(corr=apply(full_result, 1, max),
                           yindex=apply(full_result, 1, which.max))
    expected$ycol <- paste0("V", 1:ncol(y))[expected$yindex]

    result <- corr_betw_matrices(x, y, "bestright")
    expect_equivalent(result, expected)

    # all best
    expect_equivalent(corr_betw_matrices(x, y, "bestpairs"),
                      data.frame(corr=numeric(0),
                                 xindex=integer(0),
                                 yindex=integer(0),
                                 xcol=character(0),
                                 ycol=character(0),
                                 stringsAsFactors=FALSE))

    expected <- data.frame(corr=full_result[full_result > 0.15],
                           xindex=row(full_result)[full_result > 0.15],
                           yindex=col(full_result)[full_result > 0.15])
    expected$xcol <- paste0("V", 1:ncol(x))[expected$xindex]
    expected$ycol <- paste0("V", 1:ncol(y))[expected$yindex]
    expected <- expected[order(expected$xindex, expected$yindex),]

    expect_equivalent(corr_betw_matrices(x, y, "bestpairs", corr_threshold=0.15),
                      expected)

})



test_that("corr_betw_matrices works, multi-core", {

    if(isnt_karl()) skip("this test only runs locally")

    set.seed(20171201)
    n <- 100
    p <- 5
    q <- 7
    x <- matrix(rnorm(p*n), ncol=p)
    y <- matrix(rnorm(q*n), ncol=q)

    # paired (use first 5 columns of y)
    result <- corr_betw_matrices(x, y[,1:p], "paired", cores=0)
    expect_equivalent(result, diag(cor(cbind(x,y[,1:p]))[1:p,p+1:p]))

    # all
    full_result <- corr_betw_matrices(x, y, "all", cores=0)
    expect_equivalent(full_result, cor(cbind(x,y))[1:p,p+1:q])

    # best right
    expected <- data.frame(corr=apply(full_result, 1, max),
                           yindex=apply(full_result, 1, which.max))
    expected$ycol <- paste0("V", 1:ncol(y))[expected$yindex]

    result <- corr_betw_matrices(x, y, "bestright", cores=0)
    expect_equivalent(result, expected)

    # all best
    expect_equivalent(corr_betw_matrices(x, y, "bestpairs", cores=0),
                      data.frame(corr=numeric(0),
                                 xindex=integer(0),
                                 yindex=integer(0),
                                 xcol=character(0),
                                 ycol=character(0),
                                 stringsAsFactors=FALSE))

    expected <- data.frame(corr=full_result[full_result > 0.15],
                           xindex=row(full_result)[full_result > 0.15],
                           yindex=col(full_result)[full_result > 0.15])
    expected$xcol <- paste0("V", 1:ncol(x))[expected$xindex]
    expected$ycol <- paste0("V", 1:ncol(y))[expected$yindex]
    expected <- expected[order(expected$xindex, expected$yindex),]

    expect_equivalent(corr_betw_matrices(x, y, "bestpairs", corr_threshold=0.15, cores=0),
                      expected)

})


test_that("corr_betw_matrices works with some missing values, multi-core", {

    if(isnt_karl()) skip("this test only runs locally")

    set.seed(20171201)
    n <- 100
    p <- 5
    q <- 7
    x <- matrix(rnorm(p*n), ncol=p)
    y <- matrix(rnorm(q*n), ncol=q)

    x[sample(1:prod(dim(x)), 5)] <- NA
    y[sample(1:prod(dim(y)), 6)] <- NA

    # paired (use first 5 columns of y)
    result <- corr_betw_matrices(x, y[,1:p], "paired", cores=0)
    expect_equivalent(result, diag(cor(cbind(x,y[,1:p]), use="pairwise.complete.obs")[1:p,p+1:p]))

    # all
    full_result <- corr_betw_matrices(x, y, "all", cores=0)
    expect_equivalent(full_result, cor(cbind(x,y), use="pairwise.complete.obs")[1:p,p+1:q])

    # best right
    expected <- data.frame(corr=apply(full_result, 1, max),
                           yindex=apply(full_result, 1, which.max))
    expected$ycol <- paste0("V", 1:ncol(y))[expected$yindex]

    result <- corr_betw_matrices(x, y, "bestright", cores=0)
    expect_equivalent(result, expected)

    # all best
    expect_equivalent(corr_betw_matrices(x, y, "bestpairs", cores=0),
                      data.frame(corr=numeric(0),
                                 xindex=integer(0),
                                 yindex=integer(0),
                                 xcol=character(0),
                                 ycol=character(0),
                                 stringsAsFactors=FALSE))

    expected <- data.frame(corr=full_result[full_result > 0.15],
                           xindex=row(full_result)[full_result > 0.15],
                           yindex=col(full_result)[full_result > 0.15])
    expected$xcol <- paste0("V", 1:ncol(x))[expected$xindex]
    expected$ycol <- paste0("V", 1:ncol(y))[expected$yindex]
    expected <- expected[order(expected$xindex, expected$yindex),]

    expect_equivalent(corr_betw_matrices(x, y, "bestpairs", corr_threshold=0.15, cores=0),
                      expected)

})
