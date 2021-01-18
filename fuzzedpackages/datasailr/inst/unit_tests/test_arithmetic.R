test_arithmetic <- function(){
  data(iris)
  code = "
Sepal.Area = Sepal.Length * Sepal.Width
Petal.Area = Petal.Length * Petal.Width
 
Sepal.Petal.Ratio = Sepal.Area / Petal.Area
 
exp = 2 ^ 5
exp2 = 2 ** 5 
exp3 = 2.2 ^ 3 
"
iris_result = iris
  iris_result = datasailr::sail(iris, code)

  iris2 = iris;
  iris2[,"Sepal.Area"] = iris2[,"Sepal.Length"] * iris2[,"Sepal.Width"]
  iris2[,"Petal.Area"] = iris2[,"Petal.Length"] * iris2[,"Petal.Width"]
  iris2[,"Sepal.Petal.Ratio"] = iris2[,"Sepal.Area"] / iris2[,"Petal.Area"]
  
  n_rows = nrow(iris2)
  iris2$exp = rep(2 ^ 5, n_rows)
  iris2$exp2 = rep(2 ^ 5, n_rows)
  iris2$exp3 = rep(2.2 ^ 3, n_rows) 

  RUnit::checkEqualsNumeric( iris_result[,"Sepal.Area"] , iris2[,"Sepal.Area"])
  RUnit::checkEqualsNumeric( iris_result[,"Petal.Area"] , iris2[,"Petal.Area"])
  RUnit::checkEqualsNumeric( iris_result[,"Sepal.Petal.Ratio"] , iris2[,"Sepal.Petal.Ratio"])
  RUnit::checkEqualsNumeric( iris_result[,"exp"] , iris2[,"exp"])
  RUnit::checkEqualsNumeric( iris_result[,"exp2"] , iris2[,"exp2"])
  RUnit::checkEqualsNumeric( iris_result[,"exp3"] , iris2[,"exp3"])
}

