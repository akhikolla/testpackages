test_if <- function(){
  data(mtcars)
  mtcars[,"name"] = rownames(mtcars)

  code = '
power = \'\'
if( hp > 145 ){
  power = "powerful"
}else if( 145 >= hp && hp > 0){
  power = "low power"
}else{
  print("hp variable has missing value")
}

efficient = ""
if( mpg > 20){
  efficient = "efficient"
}else if( 20 >= mpg && mpg > 0 ){
  efficient = "inefficient"
}else{
  print("mpg variable has missing value")
}

print(cyl)
description = name + " is " + power + " " + efficient + " car"
  '

  mtcars_result = datasailr::sail(mtcars, code)
  mtcars2 = mtcars

  mtcars2$power = ifelse( mtcars2$hp > 145 , "powerful" , "low power" )
  mtcars2$efficient = ifelse( mtcars2$mpg > 20 , "efficient" , "inefficient" )
  mtcars2$description = paste(mtcars2$name, "is", mtcars2$power, mtcars2$efficient, "car" )
  
  RUnit::checkEquals( mtcars_result[,"power"] , mtcars2[,"power"] )
  RUnit::checkEquals( mtcars_result[,"efficient"] ,   mtcars2[,"efficient"] )
  RUnit::checkEquals( mtcars_result[,"description"] , mtcars2[,"description"] )
}

