# global variables added so that there are no notes when running R CMD check

if(getRversion() >= "2.15.1"){
  # global variables necessary because of processSlots()
  fromProcessSlots <- c('Trace','Year','Variable','Value','Month')
  utils::globalVariables(c(fromProcessSlots, "."))
}
