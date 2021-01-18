## This file follows the structure of aaa.R in the free group package.

## Define some checker functions, and call them at the end.  They
## should all return TRUE if the package works, and stop with error if
## a test is failed.  Function checker1() has one argument, checker2()
## two, and checker3() has three.  Equation numbers are from Hestenes.

test_that("Test suite aab.R",{

checker1 <- function(A){

  expect_true(A == +A)
  expect_true(A == -(-A))
  expect_error(!A)

  expect_true(A == A+0) # 1.6
  expect_false(A == 1+A) # 1.7
  expect_false(1+A == A) # 1.7
  expect_false(A == A+1)
    
  expect_true(A+A == 2*A)
  expect_true(A+A == A*2)

  expect_true(A-A == as.clifford(0))  # 1.8
  expect_true(is.zero(A-A))    # 1.8
  expect_true(A-A == as.clifford(0)) # 1.8
  expect_true(A+A+A == 3*A)
  expect_true(A+A+A == A*3)

  expect_true(A/2 + A/2 == A)

  expect_error(A&A)
  expect_true(A*A == A^2)

  expect_true(is.zero(A %^% as.clifford(0)))
  expect_true(is.zero(A %.% as.clifford(0)))

  expect_true(A^0 == as.clifford(1))
  expect_true(A^1 ==     A)
  expect_true(A^2 ==   A*A)
  expect_true(A^3 == A*A*A)

  expect_true(is.homog(grade(A,0)))
  expect_true(is.homog(grade(A,1)))
  expect_true(is.homog(grade(A,2)))
  expect_true(is.homog(grade(A,3)))

  expect_true(rev(rev(A)) == A)

  
  for(r in rstloop){
    expect_true(grade(grade(A,r,drop=FALSE),r) == grade(A,r)) # 1.12; grade() is idempotent
    expect_true(rev(grade(A,0)) == grade(A,0))
    for(lam in lamloop){
      expect_true(grade(lam*A,r,drop=FALSE) == lam*grade(A,r,drop=FALSE))  # 1.11
    }
  }

  total <- as.clifford(0)
  for(r in unique(grades(A))){
    total <- total + grade(A,r)
  }
  expect_true(A == total)  # 1.9
  if(signature()==0){
      expect_true(grade(grade(A,1,drop=FALSE)*grade(A,1,drop=FALSE),0)>=0) # 1.13
  }
  expect_true(grade(rev(A),0) == grade(A,0)) # 1.17c
  expect_true(rev(grade(A,1)) == grade(A,1)) # 1.17d

  for(lambda in lamloop){
    Ar <- grade(A,r,drop=FALSE)
    expect_equal(lambda*Ar, Ar %^% lambda) # 1.22b
  }

  expect_true(is.even(evenpart(A)))
  expect_equal(evenpart(evenpart(A)),evenpart(A))

  expect_true(is.odd(oddpart(A)))
  expect_equal(oddpart(oddpart(A)),oddpart(A))

  expect_equal(A,evenpart(A)+oddpart(A))

  expect_true(is.odd(A - evenpart(A)))
  expect_true(is.even(A - oddpart(A)))

  expect_visible(summary(A))
  expect_visible(as.character(A))
  expect_visible(as.character(-A))
  
}   # checker1() closes
  
checker2 <- function(A,B){
  expect_true(A+B == B+A) # 1.1
  expect_true(A+2*B == B+B+A)
  for(r in rstloop){
    Ar <- grade(A,r,drop=FALSE)
    Br <- grade(B,r,drop=FALSE)
    expect_true(grade(A+B,r,drop=FALSE) == Ar+Br)  # 1.10
    expect_equal(as.clifford(Ar %star% Br), Ar %.% Br) # 1.45b
  }

  expect_true(rev(A*B) == rev(B)*rev(A))  # 1.17a
  expect_true(rev(A + B) == rev(B) + rev(A)) # 1.17b

  for(r in rstloop){
    for(s in rstloop){
      LHS <- grade(A,r,drop=FALSE) %.% grade(B,s,drop=FALSE)
      RHS <- grade(grade(A,r,drop=FALSE)*grade(B,s,drop=FALSE),abs(r-s),drop=FALSE)
      expect_true(LHS == RHS) # 1.21a

      if((r==0) | (s==0)){
        LHS <-  grade(A,r,drop=FALSE) %.% grade(B,s,drop=FALSE)
        RHS <- as.clifford(0)
        expect_true(LHS == RHS)
      }
      
      Ar <- grade(A,r,drop=FALSE)
      Bs <- grade(B,s,drop=FALSE)
      expect_equal(Ar %^% Bs , grade(Ar*Bs,r+s,drop=FALSE)) # 1.22a
      if(r<=s){
        expect_equal(Ar %.% Bs , (-1)^(r*(s-1))*Bs %.% Ar)  # 1.23a
      }
      expect_equal(Ar %^% Bs, (-1)^(r*s)*Bs %^% Ar)         # 1.23b
    } # s loop closes
  } # r loop closes

   dotprod <- as.clifford(0)
     cprod <- as.clifford(0)
  starprod <- 0
  for(r in unique(grades(A))){
    for(s in unique(grades(B))){
      Ar <- grade(A,r,drop=FALSE)
      Bs <- grade(B,s,drop=FALSE)
      dotprod <-   dotprod + Ar %.% Bs
      cprod   <-     cprod + Ar %^% Bs
      if(r !=s){
        expect_true(Ar %star% Bs == 0) # 1.45a
      } else {
        starprod <- starprod + Ar %star% Bs
      }
    } # s loop closes
  } # r loop closes
  expect_true(dotprod == A %.% B)     # 1.21c
  expect_true(  cprod == A %^% B)     # 1.22c
  expect_true(starprod == A %star% B) # 1.46
  
  
  expect_true(A %star% B == grade(A*B,0))  # 1.44

  expect_true(A %star% B == grade(A*B,0)) # 1.47a
  expect_true(A %star% B == grade(B*A,0)) # 1.47a
  expect_true(A %star% B == B %star% A)   # 1.47a

  expect_true(A %star% B == rev(A) %star% rev(B)) # 1.48

  expect_true(A %X% B + B %X% A == as.clifford(0))

  ## Now some checks of the Dorst products:
  expect_true(Conj(A %|_% B) == Conj(B) %_|% Conj(A))
  expect_true(A %_|% B + A %|_%B == A %star% B + A %o% B)

}   # checker2() closes

checker3 <- function(A,B,C){
  expect_true(A+(B+C) == (A+B)+C)  # addition is associative; 1.2
  expect_true(A*(B*C) == (A*B)*C)  # geometric product is associative; 1.3


  expect_true(A*(B+C) == A*B + A*C) # left distributive; 1.4
  expect_true((A+B)*C == A*C + B*C) # right distributive; 1.5

  expect_equal( A %.% (B+C) , A %.% B + A %.% C)  # 1.24a

  expect_true(A %^% (B %^% C) == (A %^% B) %^% C) # 1.25a
  expect_true(A %^% (B + C) == A %^% B + A %^% C) # 1.24b

  expect_true(A %star% (2*B + 3*C) == 2*A %star% B + 3*A %star% C) # 1.47b


  expect_true(is.zero(A %X% (B %X% C) + B %X% (C %X% A) + C %X% (A %X% B))) # 1.56c
  expect_true(A %X% (B*C) == (A %X% B)*C + B*(A %X% C))  # 1.57

  for(r in rstloop){
    for(s in rstloop){
      for(t in rstloop){
        if((r+s <=t) & (r>0) & (s>0)){
          Ar <- grade(A,r,drop=FALSE)
          Bs <- grade(B,s,drop=FALSE)
          Ct <- grade(C,t,drop=FALSE)
          expect_equal(Ar %.% (Bs %.% Ct) , (Ar %^% Bs) %.% Ct) # 1.25b
        }
      }
    }
  }

  ## Following checks from Chisholm, "Geometric algebra", arXiv:1205.5935v1, 27 May 2012
  expect_true(A %_|% (B %|_% C) == (A %_|% B) %|_% C)  # eqn 83
  expect_true(A %_|% (B %_|% C) == (A %^%  B) %_|% C)
  expect_true(A %|_% (B %^%  C) == (A %|_% B) %|_% C)
  
  ## Following checks from Dorst 
  expect_true((A %^% B) %star% C == A %star% (B %_|% C))   # 2.2.4  NB the LHS takes much longer to evaluate than the RHS
  expect_true(C %star% (B %^% C) == (C %|_% B) %star% C)


}  # checker3() closes
  
for(i in 1:10){
    for(sigs in 0:2){
        signature(sigs)
        rstloop <- 0:4
        lamloop <- 0:2
        A <- rcliff(include.fewer=TRUE)
        B <- rcliff(5)
        C <- rcliff(5)
        
        checker1(A)
        checker2(A,B)
        checker3(A,B,C)
    }
}

})
