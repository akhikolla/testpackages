test_that("state functions work", {
  f1 <- function(x,newkeyv){
    stopifnot(is.data.table(x))
    state <- savestate(x)
    tryCatch(expr = {
      #do some code which may alter the key of x
      setkeyv(x,newkeyv)
      ##add a temporary column
      x[, temp:=NA]
      ##ooops but the function broke and now the key of x is altered
      c(1,2)%*%c(1,2,3) #generate an error
      out <- 3

    },
    error=function(e){
      setstate(x,state)
      stop(e)
    })

    setstate(x,state)
    out
  }

  x <- data.table(id=c(1,1,2,2,2),
                  id2=c(5,4,3,2,1),
                  value=c(9,2,4,2,8),
                  row_n=c(1,2,3,4,5)
  )
  setkey(x,id)
  x_original <- copy(x)
  expect_error(f1(x,newkeyv="id2"),"non-conformable")
  expect_equal(x,x_original)
})

