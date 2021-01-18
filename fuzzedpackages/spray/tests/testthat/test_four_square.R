## Four-square identity for sprays

test_that("four-square identity",{
    # first a cute test:
    LHS <- spray(cbind(diag(2,4),diag(0,4)))*spray(cbind(diag(0,4),diag(2,4)))

    RHS <-
  spray(kronecker(t(rep(1,2)),diag(4)))^2+
    spray(cbind(diag(4),magic::adiag(1-diag(2),1-diag(2))),c(1,-1,1,-1))^2 + 
      spray(cbind(diag(4),kronecker(1-diag(2),diag(2))), c(1,-1,-1,1))^2 +
        spray(cbind(diag(4),magic::arev(diag(4),1)), c(1,1,-1,-1))^2

    expect_true(is.zero(LHS-RHS))
    expect_true(is.zero(RHS-LHS))
    expect_true(LHS == RHS)


    # Now a very tough test:
    foo <- function(
                    a1,a2,a3,a4,
                    b1,b2,b3,b4)
    {
        LHS <- (a1^2+a2^2+a3^2+a4^2)*(b1^2+b2^2+b3^2+b4^2)
        RHS <- (
            (a1*b1-a2*b2-a3*b3-a4*b4)^2 +
            (a1*b2+a2*b1+a3*b4-a4*b3)^2 +
            (a1*b3-a2*b4+a3*b1+a4*b2)^2 +
            (a1*b4+a2*b3-a3*b2+a4*b1)^2
        )

    expect_true(is.zero(LHS-RHS))
    expect_true(is.zero(RHS-LHS))
    expect_true(LHS == RHS)
    }

    for(i in 1:4){
        foo(
            rspray(2,powers=0:9),
            rspray(2,powers=0:9),
            rspray(2,powers=0:9),
            rspray(2,powers=0:9),
            rspray(2,powers=0:9),
            rspray(2,powers=0:9),
            rspray(2,powers=0:9),
            rspray(2,powers=0:9)
        )
    }
})

