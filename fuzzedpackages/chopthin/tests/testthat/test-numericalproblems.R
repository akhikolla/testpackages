context("extreme cases")

test_that("When calibrating a, h=N does not hold due to numerical reasons",{
    wlist <- list( c(rep(0.0010152336491610143739,214),rep(0.0029576782456277052845,246),rep(0.002,40)),
                   c(rep(0.00105,214),rep(0.0029,246),rep(0.002,40)))
    for (w in wlist){
        N <- length(w)
        expect_true(length(chopthin(w=w,N=N)$weights)==length(w))
        expect_true(all(sort(chopthin(w=w,N=N,normalise=FALSE)$weights)==sort(w)))
        expect_true(all(sort(chopthin(w=w,N=N)$indices)==1:N))
        expect_true(all(sort(chopthin(w=c(0.,w),N=N)$indices)==2:(N+1)))
        expect_true(all(sort(chopthin(w=c(w,0.,0.),N=N)$indices)==1:N))
        expect_true(all(sort(chopthin(w=c(0.,w,0.),N=N)$indices)==2:(N+1)))
        expect_true(all(sort(chopthin(w=c(w[1:2],0.,w[3:N]),N=N)$indices)==c(1,2,4:(N+1))))
    }
})
