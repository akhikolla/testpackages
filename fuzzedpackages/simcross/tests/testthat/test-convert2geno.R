context("convert2geno")

test_that("convert2geno works for 2-allele case", {

    # data as list with alleles in intervals
    #    (as produced by sim_from_pedigree)
    dat <- list("1" = list(mat = list(alleles = c(2L, 1L, 2L, 1L, 2L, 1L),
                                      locations = c(1.63593688048422, 50.6488258950412,
                                                    54.441562271677, 72.7342922240496,
                                                    80.7659703074023, 100)),
                           pat = list(alleles = c(2L, 1L, 2L, 1L, 2L, 1L),
                                      locations = c(8.94253242295235, 43.9252462703735,
                                                    82.6503853080794, 90.0170957203954,
                                                    91.0740879364312, 100))),
                "2" = list(mat = list(alleles = c(1L, 2L, 1L, 2L, 1L),
                                      locations = c(1.72073859721422, 2.41351707372814,
                                                    72.7342922240496, 80.7659703074023, 100)),
                           pat = list(alleles = c(2L, 1L, 2L, 1L, 2L, 1L, 2L, 1L),
                                      locations = c(8.94253242295235, 43.9252462703735,
                                                    54.1366793680936, 75.1830249791965,
                                                    82.6503853080794, 90.0170957203954,
                                                    91.0740879364312, 100))),
                "3" = list(mat = list(alleles = c(1L, 2L, 1L, 2L, 1L, 2L, 1L),
                                      locations = c(1.72073859721422, 2.41351707372814,
                                                    50.6488258950412, 54.441562271677,
                                                    89.9764276109636, 97.2906407434493, 100)),
                           pat = list(alleles = c(2L, 1L, 2L, 1L),
                                      locations = c(8.94253242295235, 43.9252462703735,
                                                    97.3880127770826, 100))),
                "4" = list(mat = list(alleles = c(2L, 1L, 2L, 1L),
                                      locations = c(1.63593688048422, 72.7342922240496,
                                                    80.7659703074023, 100)),
                           pat = list(alleles = c(2L, 1L, 2L, 1L, 2L, 1L, 2L, 1L),
                                      locations = c(20.5045772017911, 38.2374913664535,
                                                    48.7418097676709, 75.1830249791965,
                                                    82.6503853080794, 88.3461387362331,
                                                    97.3880127770826, 100))))

    # marker map
    map <- seq(0, 100, by=5)
    names(map) <- paste0("marker", seq(along=map))

    # expected genotype matrix (force integers)
    expected <- rbind(c(3, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 1, 1, 1, 1),
                      c(2, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 2, 3, 1, 1, 1, 1),
                      c(2, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 1),
                      c(3, 2, 2, 2, 2, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 2, 3, 1, 2, 2, 1))
    storage.mode(expected) <- "integer"
    dimnames(expected) <- list(as.character(1:4), paste0("marker", 1:21))

    expect_equal(expected, convert2geno(dat, map))

    # with founder genotypes
    founder_geno <- rbind(rep(1, length(map)),
                          rep(2, length(map)))
    expect_equal(expected, convert2geno(dat, map, founder_geno))

    # switch founder alleles
    founder_geno <- rbind(rep(2, length(map)),
                          rep(1, length(map)))
    expected2 <- 4 - expected
    expect_equal(expected2, convert2geno(dat, map, founder_geno))

    # random founder alleles
    founder_geno <- rbind(c(2, 1, 2, 2, 2, 2, 1, 1, 1, 2, 1, 2, 1, 2, 1,
                            2, 2, 2, 1, 1, 2),
                          c(2, 1, 2, 1, 2, 2, 2, 2, 2, 1, 1, 1, 2, 1, 2,
                            2, 2, 1, 2, 2, 2))
    expected3 <- expected
    for(i in 1:ncol(expected3)) {
        if(founder_geno[1,i]==1 & founder_geno[2,i]==1)
            expected3[,i] <- 1
        else if(founder_geno[1,i]==2 & founder_geno[2,i]==2)
            expected3[,i] <- 3
        else if(founder_geno[1,i]==2 & founder_geno[2,i]==1)
            expected3[,i] <- 4-expected[,i]
    }
    expect_equal(expected3, convert2geno(dat, map, founder_geno))

})


test_that("convert2geno works for 8-allele case", {

    # data as list with alleles in intervals
    #    (as produced by sim_from_pedigree)
    dat <- list("1" = list(mat = list(alleles=c(8L, 2L, 8L, 5L, 3L, 5L, 3L),
                                      locations=c(8.20657967124134, 18.4184361249208,
                                                  55.8109006844461, 86.3054546294734,
                                                  88.7538079405203, 90.0979985482991, 100)),
                           pat = list(alleles = c(8L, 2L, 1L, 8L, 5L, 3L, 5L, 3L),
                                      locations=c(21.6286054579541, 27.2114308550954,
                                                  32.0876344339922, 55.8109006844461,
                                                  86.3054546294734, 88.7538079405203,
                                                  90.0979985482991, 100))),
                "2" = list(mat = list(alleles=c(8L, 2L, 1L, 8L, 5L, 3L, 5L, 3L),
                                      locations=c(21.6286054579541, 27.2114308550954,
                                                  32.7756949001923, 55.8109006844461,
                                                  86.3054546294734, 88.7538079405203,
                                                  90.0979985482991, 100)),
                           pat = list(alleles = c(8L, 2L, 1L, 8L, 5L, 3L, 5L, 3L),
                                      locations=c(21.6286054579541, 27.2114308550954,
                                                  29.4709661742672, 55.8109006844461,
                                                  86.3054546294734, 88.7538079405203,
                                                  90.0979985482991, 100))),
                "3" = list(mat = list(alleles=c(8L, 2L, 8L, 5L, 3L, 5L, 3L),
                                      locations=c(8.20657967124134, 18.4184361249208,
                                                  55.8109006844461, 86.3054546294734,
                                                  88.7538079405203, 90.0979985482991, 100)),
                           pat = list(alleles = c(8L, 2L, 1L, 8L, 5L, 3L, 5L, 3L),
                                      locations=c(21.6286054579541, 27.2114308550954,
                                                  32.7756949001923, 55.8109006844461,
                                                  86.3054546294734, 88.7538079405203,
                                                  90.0979985482991, 100))),
                "4" = list(mat = list(alleles=c(8L, 2L, 8L, 2L, 1L, 8L, 5L, 3L, 5L, 3L),
                                      locations=c(8.20657967124134, 12.613788805902,
                                                  21.6286054579541, 27.2114308550954,
                                                  32.0876344339922, 55.8109006844461,
                                                  86.3054546294734, 88.7538079405203,
                                                  90.0979985482991, 100)),
                           pat = list(alleles = c(8L, 2L, 1L, 8L, 5L, 3L, 5L, 3L),
                                      locations=c(21.6286054579541, 27.2114308550954,
                                                  29.4709661742672, 55.8109006844461,
                                                  86.3054546294734, 88.7538079405203,
                                                  90.0979985482991, 100))))


    # marker map
    map <- seq(0, 100, by=5)
    names(map) <- paste0("marker", seq(along=map))

    # Use get_geno to construct expected array
    expected <- array(dim=c(length(dat), length(map), 2))
    dimnames(expected) <- list(names(dat), names(map), c("mat", "pat"))
    for(i in seq(along=map))
        expected[,i,] <- get_geno(dat, map[i])
    storage.mode(expected) <- "integer"

    expect_equal(expected, convert2geno(dat, map))

    # now test with use of founder genotypes
    founder_geno <- rbind(c(2, 2, 2, 2, 2, 1, 2, 2, 2, 1, 1, 2, 1, 1, 1,
                            2, 1, 1, 1, 2, 2),
                          c(2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 1, 1, 2, 2,
                            1, 1, 2, 1, 1, 2),
                          c(2, 2, 1, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 2,
                            1, 2, 2, 2, 1, 1),
                          c(1, 2, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 2,
                            2, 1, 1, 2, 1, 1),
                          c(2, 2, 1, 1, 2, 2, 1, 2, 2, 1, 2, 1, 2, 1, 2,
                            1, 2, 2, 2, 1, 2),
                          c(1, 2, 2, 1, 2, 1, 1, 2, 2, 2, 2, 2, 2, 1, 1,
                            1, 2, 2, 1, 1, 2),
                          c(2, 2, 1, 2, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2, 1,
                            1, 2, 2, 2, 1, 1),
                          c(2, 2, 2, 1, 1, 2, 1, 1, 2, 2, 2, 2, 2, 1, 2,
                            2, 2, 2, 2, 1, 1))
    storage.mode(founder_geno) <- "integer"

    expected2 <- rbind(c(3, 3, 3, 2, 1, 3, 2, 1, 3, 3, 3, 3, 3, 1, 3, 1,
                         3, 3, 3, 1, 1),
                       c(3, 3, 3, 1, 1, 3, 2, 1, 3, 3, 3, 3, 3, 1, 3, 1,
                         3, 3, 3, 1, 1),
                       c(3, 3, 3, 2, 1, 3, 2, 1, 3, 3, 3, 3, 3, 1, 3, 1,
                         3, 3, 3, 1, 1),
                       c(3, 3, 3, 1, 1, 3, 2, 1, 3, 3, 3, 3, 3, 1, 3, 1,
                         3, 3, 3, 1, 1))
    storage.mode(expected2) <- "integer"
    dimnames(expected2) <- list(as.character(1:4), paste0("marker", 1:21))

    expect_equal(expected2, convert2geno(dat, map, founder_geno))
})

test_that("convert2geno and get_geno give same answer when shifted", {

    # marker map starting at 5
    map <- seq(20, 80, by=5)
    names(map) <- paste0("marker", seq(along=map))

    # simulate large ail pedigree
    set.seed(37513076)
    tab <- sim_ail_pedigree(6, 100)
    dat <- sim_from_pedigree(tab, max(map))[tab[,"gen"]==6]

    g <- convert2geno(dat, map)
    g_shifted <- convert2geno(dat, map, shift_map=TRUE)

    expect_equal(g[,map==30], g_shifted[,map==50])

    qtl_g <- get_geno(dat, 30)
    qtl_g <- (qtl_g[,1] + qtl_g[,2] - 1)
    expect_equal(qtl_g, g[,map==30])
    expect_equal(qtl_g, g_shifted[,map==50])

})
