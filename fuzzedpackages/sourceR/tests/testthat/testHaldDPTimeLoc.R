context("Test time/location data structure")

testdata = function()
{
  # Test priors
  priors <-
    list(
      a_theta = 0.01,
      b_theta = 0.00001,
      a_alpha = 1,
      a_r = 0.1
    )
  # Data
  set.seed(1)
  types = unique(sim_SA$cases$Type)[1:10]

  X = filter(sim_SA$sources, Type %in% types)
  y = filter(sim_SA$cases, Type %in% types)
  list(
    y = y,
    X = X,
    prior = priors,
    prev = sim_SA$prev
  )
}

test_that("Detect y/X data mismatch", {
  dat = testdata()
  y = Y(
    data = dat$y,
    y = 'Human',
    type = 'Type',
    time = 'Time',
    location = 'Location'
  )
  x = X(
    data = dat$X,
    x = 'Count',
    type = 'Type',
    time = 'Time',
    source = 'Source'
  )
  k = Prev(dat$prev,
           prev = 'Value',
           time = 'Time',
           source = 'Source')

  expect_silent(HaldDP(
    y = y,
    x = x,
    k = k,
    priors = dat$prior,
    a_q = 0.1
  ))

  dat$X = dat$X[dat$X$Time == '1', ]
  x = X(
    data = dat$X,
    x = 'Count',
    type = 'Type',
    time = 'Time',
    source = 'Source'
  )
  expect_error(HaldDP(
    y = y,
    x = x,
    k = k,
    priors = dat$prior,
    a_q = 0.1
  ),
  "Times in x and y do not match")
})
