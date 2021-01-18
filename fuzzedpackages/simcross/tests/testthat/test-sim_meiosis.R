
context("meiosis")

test_that("create_parent works", {

    expected_p1 <- list(mat=list(alleles=1L, locations=100.0),
                        pat=list(alleles=1L, locations=100.0))

    expect_equal(create_parent(100, 1), expected_p1)

    expected_p2 <- list(mat=list(alleles=2L, locations=100.0),
                        pat=list(alleles=2L, locations=100.0))

    expect_equal(create_parent(100, 2), expected_p2)

    expected_f1 <- list(mat=list(alleles=1L, locations=100.0),
                        pat=list(alleles=2L, locations=100.0))

    expect_equal(create_parent(100, 1:2), expected_f1)

    # numeric alleles converted to integers
    expect_equal(create_parent(100, c(1.9, 2.9)), expected_f1)

    # non-numeric alleles not allowed
    expect_error(create_parent(100, "a"))

    # can't have more than 2 alleles
    expect_error(create_parent(100, 1:3))

    # must have at least 1 allele
    expect_error(create_parent(100, NULL))
    expect_error(create_parent(100, numeric(0)))

})


test_that("check_individual works", {

    expected_p1 <- list(mat=list(alleles=1L, locations=100.0),
                        pat=list(alleles=1L, locations=100.0))

    expected_p2 <- list(mat=list(alleles=2L, locations=100.0),
                        pat=list(alleles=2L, locations=100.0))

    expected_f1 <- list(mat=list(alleles=1L, locations=100.0),
                        pat=list(alleles=2L, locations=100.0))

    # should be clean
    expect_true(check_individual(expected_p1))
    expect_true(check_individual(expected_p2))
    expect_true(check_individual(expected_f1))

    # no data
    expect_error(check_individual(NULL))
    expect_error(check_individual(list(mat=NULL, pat=NULL)))

    # alleles not integers
    z <- list(mat=list(alleles=1, locations=100),
              pat=list(alleles=2, locations=100))
    expect_error(check_individual(z))
    z <- list(mat=list(alleles="a", locations=100),
              pat=list(alleles="b", locations=100))
    expect_error(check_individual(z))

    # wrong orders
    z <- list(mat=list(locations=100, alleles=1L),
              pat=list(alleles=1L, locations=100))
    expect_error(check_individual(z))
    z <- list(mat=list(alleles=1L, locations=100),
              pat=list(locations=100, alleles=1L))
    expect_error(check_individual(z))
    z <- list(pat=list(alleles=2L, locations=100),
              mat=list(alleles=1L, locations=100))
    expect_error(check_individual(z))

    # locations not in order
    z <- list(mat=list(alleles=c(2L, 1L), locations=c(5, 100)),
              pat=list(alleles=c(1L, 2L), locations=c(100, 85)))
    expect_error(check_individual(z))
    z <- list(mat=list(alleles=c(2L, 1L), locations=c(500, 100)),
              pat=list(alleles=c(1L, 2L), locations=c(100, 105)))
    expect_error(check_individual(z))
})

test_that("simulations of meiosis and crosses work", {

    # these are all rather contrived, but it's a start

    seed <- 79251223
    set.seed(seed)

    expected <- c(36.1079982714728, 51.9751604413614, 204.360923892818)
    expect_equal(sim_crossovers(300, 10, 0.3), expected)

    allele <- ifelse(runif(1) < 0.5, 1L, 2L)
    expect_equal(allele, 1L)

    another_set <- sim_crossovers(300, 3, 0.01)
    another_allele <- ifelse(runif(1) < 0.5, 1L, 2L)
    expect_equal(another_allele, 2L)

    p1 <- create_parent(300, 1)
    p2 <- create_parent(300, 2)
    f1 <- cross(p1, p2)
    expect_equal(create_parent(300, 1:2), f1)

    set.seed(seed)
    expected2 <- list(alleles=c(1L, 2L, 1L, 2L), locations=c(expected, 300))
    expect_equal(sim_meiosis(f1, m=10, p=0.3), expected2)

    expected3 <- list(alleles=c(2L, 1L, 2L), locations=c(another_set, 300))
    expect_equal(sim_meiosis(f1, m=3, p=0.01), expected3)

    set.seed(seed)
    f2 <- cross(f1, f1, m=10, p=0.3)
    expected <- list(mat=expected2,
                     pat=list(alleles=c(1L, 2L, 1L),
                              locations=c(22.7139970986173, 234.89673093427, 300)))
    expect_equal(f2, expected)

    set.seed(seed)
    junk <- sim_meiosis(f1, m=10, p=0.3)
    f2 <- cross(f1, f1, m=3, p=0.01)
    expected <- list(mat=expected3,
                     pat=list(alleles=c(2L, 1L, 2L, 1L),
                     locations=c(72.8879569564015, 194.954296317883, 285.253091622144, 300)))
    expect_equal(f2, expected)

    set.seed(seed)
    p1 <- create_parent(300, 1:2)
    p2 <- create_parent(300, 3:4)
    f1a <- cross(p1, p2)
    f2a <- cross(p1, p2)
    sib <- cross(f1a, f2a)
    expected <- list(mat=list(alleles=c(4L, 3L, 2L, 3L, 2L, 1L, 4L),
                                    locations=c(113.656340399757,
                                                119.926706166007, 151.954430225305, 188.945353589952,
                                                214.686310663819, 242.297911760397, 300)),
                     pat=list(alleles=c(4L, 3L, 4L, 3L, 4L),
                                    locations=c(130.161284795031, 178.185045230202, 224.519506352954,
                                                250.761620071717, 300)))

    expect_equal(sib, expected)

})


test_that("average no. crossovers match expectation", {

    set.seed(82271836)
    tol <- 0.01
    n <- 10000

    x <- replicate(n, sim_crossovers(L=100, m=0, p=0))
    expect_equal(mean(sapply(x, length)), 1, tolerance=tol)

    Lstar <- calc_Lstar(100, 0, 0)
    x <- replicate(n, sim_crossovers(L=100, m=0, p=0, obligate_chiasma=TRUE, Lstar=Lstar))
    expect_equal(mean(sapply(x, length)), 1, tolerance=tol)

    x <- replicate(n, sim_crossovers(L=100, m=1, p=0))
    expect_equal(mean(sapply(x, length)), 1, tolerance=tol)

    Lstar <- calc_Lstar(100, 1, 0)
    x <- replicate(n, sim_crossovers(L=100, m=1, p=0, obligate_chiasma=TRUE, Lstar=Lstar))
    expect_equal(mean(sapply(x, length)), 1, tolerance=tol)

    set.seed(82271836)
    x <- replicate(n, sim_crossovers(L=100, m=10, p=0))
    expect_equal(mean(sapply(x, length)), 1, tolerance=tol)

    Lstar <- calc_Lstar(100, 10, 0)
    x <- replicate(n, sim_crossovers(L=100, m=10, p=0, obligate_chiasma=TRUE, Lstar=Lstar))
    expect_equal(mean(sapply(x, length)), 1, tolerance=tol)

    x <- replicate(n, sim_crossovers(L=100, m=10, p=0.1))
    expect_equal(mean(sapply(x, length)), 1, tolerance=tol)

    Lstar <- calc_Lstar(100, 10, 0.1)
    x <- replicate(n, sim_crossovers(L=100, m=10, p=0.1, obligate_chiasma=TRUE, Lstar=Lstar))
    expect_equal(mean(sapply(x, length)), 1, tolerance=tol)

})
