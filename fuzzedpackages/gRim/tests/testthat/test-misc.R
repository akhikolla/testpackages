context("gRim model")

##############################

setsetequal <- function(lst1, lst2){
    r1 <- all(sapply(lst1, function(x) is_inset(x, lst2)))
    r2 <- all(sapply(lst2, function(x) is_inset(x, lst1)))
    r1 && r2
}

set.seed(1212)
Y1 <- rbinom(500, 1, 0.5)
Y2 <- Y1
Y2[3:10] <- rep(1,8)
Y3 <- rbinom(500, 1, 0.5)
df <- data.frame(Y1, Y2, Y3)

test_that("dmod()", {

    m1 <- dmod(~.^1, data=df)
    ms <- dmod(~.^., data=df)

    detail = 0
    m2 <- forward(m1, criterion="test", details=detail,type = "decomposable") #OK
    m3 <- forward(m1, criterion="test", details=detail,type = "unrestricted") #OK
    m4 <- backward(ms, k=2, details=detail,type = "decomposable") ## OK
    m5 <- backward(ms, k=2, details=detail,type = "unrestricted") ## OK

    
    expect_true(setsetequal(list("Y1", "Y2", "Y3"), terms(m1)))
    expect_true(setsetequal(list(c("Y1", "Y2", "Y3")), terms(ms)))
    
    t0 <- list(c("Y1", "Y2"), "Y3")
    
    expect_true(setsetequal(t0, terms(m2)))
    expect_true(setsetequal(t0, terms(m3)))
    expect_true(setsetequal(t0, terms(m4)))
    expect_true(setsetequal(t0, terms(m5)))
})
