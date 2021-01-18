## Increasing coverage of spray_ops.R 

test_that("test suite aab",{

    checker1 <- function(S){
        expect_error(!S)
        expect_true(+S == S)
        expect_error(S&S)
        expect_silent(jj <- 4*S)
        expect_silent(jj <- S*4)
        expect_silent(jj <- 4+S)
        expect_silent(jj <- S+4)
        expect_silent(jj <- S/4)
        expect_silent(jj <- -(S*0))
        expect_silent(jj <- -(0*S))
        expect_silent(jj <- S*(0*S))
        expect_silent(jj <- S*(S*0))
        expect_silent(jj <- (0*S)*S)
        expect_silent(jj <- (0*S)*S)
        expect_true(S != S+3)
        expect_true(S != 3+S)
        expect_error(4/S)
        expect_error(S^S)
        expect_error(S^(-9))
        expect_error(S == 6)
        expect_true(S == S + (S*0))
        expect_true(S == (S*0) + S)
        expect_error(S != 6)
        expect_false(1000+S == S)
        expect_false(13+S == S*0)

        SS1 <- rspray(6,arity=119)
        SS2 <- rspray(6,arity=119)
        expect_false(SS1==SS1+SS2)
        expect_false(SS2==SS1+SS2)
        expect_false(SS1+SS2==SS1)
        expect_false(SS1+SS2==SS2)

        return(TRUE)
    }

    checker2 <- function(S1,S2){
        expect_false(S1==S2)
        return(TRUE)
    }

    checker2a <- function(S1,S2){
        expect_false(S1==S2)
        return(TRUE)
    }

    
    for(i in 1:3){
        checker1(rspray(8))
        checker2(rspray(8,arity=3),rspray(9,arity=4))
        checker2(rspray(8,arity=99),rspray(8,arity=99))
        checker2a(rspray(4,arity=119),rspray(5,arity=119))

    }

})

