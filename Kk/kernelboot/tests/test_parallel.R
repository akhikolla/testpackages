
library("kernelboot")

if ( future::availableCores() > 1L ) {

  # simply check if it fails

  stopifnot( kernelboot(mtcars, function(data) coef(lm(mpg ~ disp + hp + drat, data = data)),
                        R = 10, parallel = TRUE, workers = 2L)$param$parallel )

  # check reproducibility

  set.seed(0xBEEF)
  s1 <- kernelboot(mtcars, function(data) coef(lm(mpg ~ disp + hp + drat, data = data)),
                   R = 10, parallel = TRUE, workers = 2L)

  set.seed(0xBEEF)
  s2 <- kernelboot(mtcars, function(data) coef(lm(mpg ~ disp + hp + drat, data = data)),
                   R = 10, parallel = TRUE, workers = 2L)

  stopifnot( all.equal(s1, s2) )

}
