context('PhenoFlex')

test_that('PhenoFlex_1', {
  data(KA_weather)
  hourtemps <- chillR::stack_hourly_temps(KA_weather, latitude=50.4)
  iSeason <- chillR::genSeason(hourtemps, years=c(2009))
  zc <- 190
  yc <- 40
  x <- chillR::PhenoFlex(temp=hourtemps$hourtemps$Temp[iSeason[[1]]],
                         times=c(1: length(hourtemps$hourtemps$Temp[iSeason[[1]]])),
                         zc=zc, stopatzc=TRUE, yc=yc, basic_output=FALSE)
  DBreakDay <- x$bloomindex
  expect_equal(DBreakDay, 6170)
})
