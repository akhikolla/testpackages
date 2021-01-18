test_regexp  <- function(){
  data(mtcars)
  mtcars[,"name"] = rownames(mtcars)
  
  code = '
germany = re/(^Merc|^Porsche|^Volvo)/
usa = re/(^Hornet|^Cadillac|^Lincoln|^Chrysler|^Dodge|^AMC|^Camaro|^Chevrolet|^Pontiac|^Ford)/
japan = re/(^Mazda|^Datsun|^Honda|^Toyota)/
austoralia = re/(^Valiant)/
france=re/(^Duster)/
italy=re/(^Fiat|^Ferrari|^Maserati)/
uk = re/(^Lotus)/

power = ""
if( hp > 145 ){
  power = "powerful"
}else if( 145 >= hp && hp > 0){
  power = "low power"
}else{
  print("hp variable has missing value")
  power = "unknown power"
}

country = ""
if(name =~ germany){
  country = "Germany"
}else if(name =~ usa){
  country = "USA"
}else if(name =~ japan){
  country = "Japan"
}else if(name =~ austoralia){
  country = "Austoralia"
}else if(name =~ france){
  country = "France"
}else if(name =~ italy){
  country = "Italy"
}else if(name =~ uk){
  country = "UK"
}else{
  country = "other country"
}

description = name + " is " +  power + " " + country + " made car."
'

  mtcars_result = datasailr::sail(mtcars, code)
  mtcars2 = mtcars

  mtcars2$power = ifelse( mtcars2$hp > 145 , "powerful" , "low power" )
  country_list = c( rep("Japan", 3), rep("USA", 2), "Austoralia", "France", rep("Germany", 7), rep("USA", 3),  "Italy", rep("Japan",3), rep("USA",4), "Italy", "Germany", "UK", "USA", rep("Italy",2), "Germany")
  mtcars2$country = country_list
  mtcars2$description = paste0(mtcars2$name , " is ", mtcars2$power, " ", mtcars2$country, " made car.")

  RUnit::checkEquals( mtcars_result[,"power"] , mtcars2[,"power"])
  RUnit::checkEquals( mtcars_result[,"country"] , mtcars2[,"country"])
  RUnit::checkEquals( mtcars_result[,"description"] , mtcars2[,"description"])

}
