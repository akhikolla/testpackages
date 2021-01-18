test_that("benchmarks", {
  if (requireNamespace("survival", quietly = TRUE)) {
    library(fastlogranktest)
    size<-1000
    repititions<-100
    maxvalue<-50

    T1 <- runif(size, 0.0, maxvalue)
    E1 <- sample(0:1, size, replace=T)
    T2 <- runif(size, 0.0, maxvalue)
    E2 <- sample(0:1, size, replace=T)
    t1s<-list(T1)[rep(1,repititions)]
    e1s<-list(E1)[rep(1,repititions)]
    t2s<-list(T2)[rep(1,repititions)]
    e2s<-list(E2)[rep(1,repititions)]



    time1<-system.time({result1=logrank_test(T1,T2,E1,E2)})[["elapsed"]]
    # show(result1)
    time2<-system.time({result2=multi_logrank_test(t1s,t2s,e1s,e2s)})[["elapsed"]]
    #show(result2)

    data<-data.frame(time=c(T1,T2),censored=c(E1,E2),group=c(rep(0,length(T1)),rep(1,length(T2))))
    object<-survival::Surv(data$time, data$censored)
    time3<-system.time({result3=survival::survdiff(object ~ data$group)})[["elapsed"]]
    #show(result3)

    print("")
    print("normal time: ")
    print(time1)
    print("multiple time: ")
    print(time2/repititions)
    print("survdiff time: ")
    print(time3)
    expect_equal(result1[3], pchisq(result3$chisq, length(result3$n)-1, lower.tail = FALSE), tolerance=1e-4)
  }
  else{
    print("")
    print("Package 'survival' is needed for this test")
    return(TRUE)
  }
})


