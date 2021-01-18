
context("get_geno")

test_that("simple test of get_geno", {

    dat <- list("1"=list(mat=list(alleles=c(1, 2, 1), locations=c(9.22, 83.6, 100)),
                         pat=list(alleles=c(2, 1), locations=c(48.7, 100))),
                "2"=list(mat=list(alleles=c(2, 1, 2), locations=c(18.64, 55.21, 100)),
                         pat=list(alleles=c(2, 1, 2), locations=c(31.55, 95.42, 100))),
                "3"=list(mat=list(alleles=c(2, 1, 2), locations=c(5.5, 39.4, 100)),
                         pat=list(alleles=c(1, 2),    locations=c(89.1, 100))))

    ex_geno <- rbind(c(2, 2), c(1, 2), c(1, 1))
    geno <- get_geno(dat, 30)
    expect_equivalent(geno, ex_geno)

    ex_geno <- rbind(c(1,2), c(2,2), c(2, 1))
    geno <- get_geno(dat, 0)
    expect_equivalent(geno, ex_geno)

    ex_geno <- rbind(c(1,1), c(2,2), c(2, 2))
    geno <- get_geno(dat, 100)
    expect_equivalent(geno, ex_geno)

    expect_error(get_geno(dat, -1))
    expect_error(get_geno(dat, 101))

})

