context('chifull')

test_that('chifull_1', {
  hourtemps <- chillR::stack_hourly_temps(KA_weather, latitude=50.4)
  par <- c(40, 190, 0.5, 25, 3372.8, 9900.3, 6319.5, 5.939917e13, 4, 36, 4, 1.6)
  SeasonList <- chillR::genSeasonList(hourtemps$hourtemps, years=c(2007,2008))
  
  x <- chillR:::chifull(par=par,
                        modelfn=chillR::PhenoFlex_GDHwrapper,
                        bloomJDays=KA_bloom$pheno[c(24,25)],
                        SeasonList=SeasonList
                        )
  expect_equal(object=x, expected=2.3767361111111, tolerance = 0.0000001)
})
