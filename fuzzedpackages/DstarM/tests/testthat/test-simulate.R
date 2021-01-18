context("simulate, estDstarM, estND, estObserved and methods.")

set.seed(42)
tt <- seq(0, 5, .2)
pdfND <- dbeta(tt, 10, 30)
n <- 100
pars <- c(1, 2, .5, .5, .5)
dat <- simData(n, pars, tt, pdfND)
fixed <- matrix(c('z1', .5,  'sz1', .5, 'sv1', .5), 2, 3)

test_that("simulation works", {

	sum <- summary(dat)
	expect_equal2(
		object = summary(dat),
		expected = structure(c("Min.   :0.400  ", "1st Qu.:0.400  ", "Median :0.400  ",
							   "Mean   :0.474  ", "3rd Qu.:0.600  ", "Max.   :1.200  ", "lower:20  ",
							   "upper:80  ", NA, NA, NA, NA, "Min.   :1  ", "1st Qu.:1  ", "Median :1  ",
							   "Mean   :1  ", "3rd Qu.:1  ", "Max.   :1  "), .Dim = c(6L, 3L
							   ), .Dimnames = list(c("", "", "", "", "", ""), c("      rt",
							   												 " response", "  condition")), class = "table")


	)
})

test_that("estDstarM useRcpp = FALSE works", {
  skip_on_cran()
  fitD <- estDstarM(data = dat, tt = tt, fixed = fixed, verbose = 0)
	expect_equal2(
		object = fitD$Bestvals,
		expected = structure(c(0.949695690096916, 1.57237466724544, 0.5, 0.5, 0.5
							   ), .Names = c("a1", "v1", "z1", "sz1", "sv1")),
		tol = 1e-4
	)
})

set.seed(42)
fitD <- estDstarM(data = dat, tt = tt, fixed = fixed, verbose = 0, useRcpp = TRUE)

test_that("estDstarM useRcpp = TRUE works", {

	expect_equal2(
		object = fitD$Bestvals,
		expected = structure(c(0.949695690096916, 1.57237466724544, 0.5, 0.5, 0.5
		), .Names = c("a1", "v1", "z1", "sz1", "sv1")),
		tol = 1e-4
	)
})
set.seed(42)
fitND <- estND(fitD, verbose = 0)

test_that("estND works", {

  expect_equal2(
    object = fitND$r.hat,
    expected = structure(
      c(3.25486872725596e-08, 4.94098337050759, 0.0543616639789932,
        1.49766131273096e-07, 1.55420026887307e-08, 0.00465478393093598,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
      .Dim = c(26L, 1L)),
    tol = 1e-3
  )
})

set.seed(42)
fitObs <- estObserved(fitD, fitND)

test_that("estObserved works", {

  expect_equal2(
    object = fitObs$obs,
    expected = structure(
      c(0, 4.7382332893854e-09, 0.719277307661969, 0.201508092222449,
        0.0556723249170272, 0.015627024291952, 0.00512304445828108, 0.00146200655907774,
        0.000422680654569136, 0.00012347483645932, 3.64029951054436e-05,
        1.0821044663176e-05, 3.24052403193299e-06, 9.76932256131796e-07,
        2.96309977694946e-07, 9.03699329878131e-08, 2.77006028471892e-08,
        8.53017568066561e-09, 2.63795161695005e-09, 8.18974200683555e-10,
        2.55174920375501e-10, 7.97724370220762e-11, 2.50155796132101e-11,
        7.8671162041282e-12, 2.480741374849e-12, 7.84241748339073e-13,
        0, 1.92923593266303e-08, 2.92863509257847, 0.787951959082506,
        0.205422214951634, 0.0545938185464467, 0.0175272397822968, 0.00477018417294568,
        0.0013165712433021, 0.000368488648507474, 0.000104426157198195,
        2.99237654415752e-05, 8.66056960077416e-06, 2.52913094609005e-06,
        7.44594172428541e-07, 2.20834998634084e-07, 6.59375130215366e-08,
        1.98090558146495e-08, 5.98467245259715e-09, 1.8174640757515e-09,
        5.54580763528876e-10, 1.69972748294523e-10, 5.23081757647044e-11,
        1.61591328553273e-11, 5.00931757947714e-12, 1.55797896370483e-12
      ), .Dim = c(26L, 2L)),
    # higher tolerance since the floating point error
    # of estDstarM and estND accumulates here.
    tol = 1e-3,
    label = "estObserved same time grid"
  )

})

test_that("rtDescriptives works", {

	obj <- rtDescriptives(data = dat, plot = FALSE, verbose = FALSE)

	expect_equal2(
		object = obj$table,
		expected = structure(list(counts = structure(list(conditionResponse = structure(c(20L, 80L), .Dim = 1:2, .Dimnames = list("1", c("lower", "upper"))),
														  condition = 100, response = c(20, 80)), .Names = c("conditionResponse",
														  												   "condition", "response")), props = structure(list(conditionResponse = structure(c(0.2,
														  												   																				  0.8), .Dim = 1:2, .Dimnames = list("1", c("lower", "upper"))),
														  												   												  condition = 1, response = c(0.2, 0.8)), .Names = c("conditionResponse",
														  												   												  												   "condition", "response")), responses = c("lower", "upper")), .Names = c("counts",
														  												   												  												   																		"props", "responses"), class = "DstarM")
	)

	f <- tempfile()
	grDevices::pdf(f)
	on.exit({
	  if (file.exists(f))
	    file.remove(f)
	  dev.off()
	})
	expect_error(capture.output(print(obj)), NA, label = "print(rtDescriptives)")

})

test_that("getPdf works", {

  pdf1 <- getPdfs(fitD)

  expect_equal2(
    object = pdf1,
    expected = structure(
      c(0, 0.727865247447104, 0.1959059183997, 0.0541816283045266,
        0.0152174888025591, 0.00433108329889229, 0.00124725384076482,
        0.000362962123803565, 0.000106619908216455, 3.1584465676862e-05,
        9.42777367605277e-06, 2.83356391486915e-06, 8.56981811418903e-07,
        2.60665880685095e-07, 7.9699946339339e-08, 2.44852424993832e-08,
        7.55536277036991e-09, 2.34078046456242e-09, 7.2792393724709e-10,
        2.2714876405139e-10, 7.11090434038004e-11, 2.23270255960626e-11,
        7.02973917409224e-12, 2.2190585660137e-12, 7.02175520399019e-13,
        2.22690957517935e-13, 0, 2.96361903380849, 0.764758044187432,
        0.199461950577259, 0.0530514141430213, 0.0143609673281179, 0.0039487031638906,
        0.00110094583345708, 0.000310799361306402, 8.87250963722666e-05,
        2.558512050515e-05, 7.44537730216164e-06, 2.18463815970493e-06,
        6.45869910520832e-07, 1.92264868928976e-07, 5.75960759637361e-08,
        1.73540295666976e-08, 5.25681387124611e-09, 1.60022430952626e-09,
        4.89344204250317e-10, 1.50272178616129e-10, 4.63278840057676e-11,
        1.43346671355407e-11, 4.45047576720503e-12, 1.38612062853476e-12,
        4.32994168631387e-13),
      .Dim = c(26L, 2L)),
    label = "getPdf without time grid."
  )

  pdf2 <- getPdfs(fitD, tt = tt*2)
  expect_equal2(
    object = pdf2,
    expected = structure(
      c(0, 0.473392700802263, 0.0367719780112711, 0.00301390008575954,
        0.000257639415501978, 2.2781543709942e-05, 2.07083551921233e-06,
        1.92589244672013e-07, 1.82569960959983e-08, 1.75897635684979e-09,
        1.71829939510607e-10, 1.69868641067257e-11, 1.6967571411535e-12,
        1.7102193963543e-13, 1.73754499873704e-14, 1.77776054796384e-15,
        1.83030977536413e-16, 1.89496182249514e-17, 1.9717497879796e-18,
        2.06092976174104e-19, 2.16295410816017e-20, 2.27845494943933e-21,
        2.40823518593757e-22, 2.55326528650017e-23, 2.71468467434776e-24,
        2.89380693458856e-25, 0, 1.84798335321093, 0.128194964336426,
        0.00954175999730018, 0.000751024523701583, 6.18246217124132e-05,
        5.27902253870849e-06, 4.64594364045617e-07, 4.19347766187439e-08,
        3.86682808750321e-09, 3.63122018322027e-10, 3.46387023211313e-11,
        3.34946178931164e-12, 3.2774935792496e-13, 3.24066967539408e-14,
        3.23389233223661e-15, 3.2536151433629e-16, 3.29741823717026e-17,
        3.36372411048145e-18, 3.45160483936721e-19, 3.56065010678367e-20,
        3.6908766632802e-21, 3.84266668313421e-22, 4.01672676799234e-23,
        4.21406209559239e-24, 4.43596200398678e-25), .Dim = c(26L, 2L
        )),
    label = "getPdf with new time grid."
  )

  pars <- fitD$Bestvals[fitD$restr.mat]
  dim(pars) <- dim(fitD$restr.mat)
  pdf3 <- getPdfs(tt = fitD$tt, pars = pars)
  expect_equal2(
    object = pdf3,
    expected = structure(
      c(0, 0.727865247447104, 0.1959059183997, 0.0541816283045266,
        0.0152174888025591, 0.00433108329889229, 0.00124725384076482,
        0.000362962123803565, 0.000106619908216455, 3.1584465676862e-05,
        9.42777367605277e-06, 2.83356391486915e-06, 8.56981811418903e-07,
        2.60665880685095e-07, 7.9699946339339e-08, 2.44852424993832e-08,
        7.55536277036991e-09, 2.34078046456242e-09, 7.2792393724709e-10,
        2.2714876405139e-10, 7.11090434038004e-11, 2.23270255960626e-11,
        7.02973917409224e-12, 2.2190585660137e-12, 7.02175520399019e-13,
        2.22690957517935e-13, 0, 2.96361903380849, 0.764758044187432,
        0.199461950577259, 0.0530514141430213, 0.0143609673281179, 0.0039487031638906,
        0.00110094583345708, 0.000310799361306402, 8.87250963722666e-05,
        2.558512050515e-05, 7.44537730216164e-06, 2.18463815970493e-06,
        6.45869910520832e-07, 1.92264868928976e-07, 5.75960759637361e-08,
        1.73540295666976e-08, 5.25681387124611e-09, 1.60022430952626e-09,
        4.89344204250317e-10, 1.50272178616129e-10, 4.63278840057676e-11,
        1.43346671355407e-11, 4.45047576720503e-12, 1.38612062853476e-12,
        4.32994168631387e-13),
      .Dim = c(26L, 2L)),
    label = "getPdf with time grid and parameters."
  )
})

test_that("methods work", {

  expect_error(capture.output(print(fitD)), NA, label = "print(fitD)")
  expect_error(capture.output(print(fitND)), NA, label = "print(fitND)")
  expect_error(capture.output(print(fitObs)), NA, label = "print(fitObs)")
  expect_error(capture.output(coef(fitD)), NA, label = "coef(fitD)")
  expect_error(capture.output(summary(fitD)), NA, label = "summary(fitD)")

  f <- tempfile()
  grDevices::pdf(f)
  on.exit({
    if (file.exists(f))
      file.remove(f)
    dev.off()
  })
  expect_error(capture.output(plot(fitD)), NA, label = "plot(fitD)")
  expect_error(capture.output(plot(fitND)), NA, label = "plot(fitND)")
  expect_error(capture.output(plot(fitObs)), NA, label = "plot(fitObs)")

})

