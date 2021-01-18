test_that("Test suite aae.R",{  # quaternions

  c1 <- clifford(list(numeric(0),c(1,2),c(1,3),c(2,3)),1:4)
  c2 <- clifford(list(numeric(0),c(1,2),c(1,3),c(2,3)),4:1)

  q1 <- clifford_to_quaternion(c1)
  q2 <- clifford_to_quaternion(c2)

  ## q1 = 1-2i-3j-4k; q2 = 4-3i-2i-k

  expect_true(quaternion_to_clifford(clifford_to_quaternion(c1)) == c1)
  expect_true(quaternion_to_clifford(clifford_to_quaternion(c2)) == c2)


  transcript_using_onion_package <- 
'
> q1 <- as.quaternion(c(1,-2,-3,-4),single=TRUE)
> q2 <- as.quaternion(c(4,-3,-2,-1),single=TRUE)
> q1*q2
   [1]
Re -12
i  -16
j   -4
k  -22
> q2*q1
   [1]
Re -12
i   -6
j  -24
k  -12
> 
'

  expect_true(
      quaternion_to_clifford(clifford_to_quaternion(c1*c2)) ==
      quaternion_to_clifford(c(-12,-16,-4,-22))
  )
  
  expect_true(
      quaternion_to_clifford(clifford_to_quaternion(c2*c1)) ==
      quaternion_to_clifford(c(-12,-6,-24,-12))
      )
    
})
