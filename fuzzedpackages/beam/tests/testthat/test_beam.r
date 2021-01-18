data("TCPAprad")

test_that("NgreaterthanP", {

  ##################################################
  ## beam estimation
  fit <- beam(X = TCPAprad, type="both", verbose=TRUE, D=diag(ncol(TCPAprad)))
  fit <- beam(X = TCPAprad, type="both", verbose=TRUE, D=NULL)

  ## Test beam-class object
  # Test table
  expect_equal(dim(fit@table), c(17766, 6))
  expect_equal(colnames(fit@table), c("m_cor", "m_logBF", "m_tail_prob", "p_cor", "p_logBF", "p_tail_prob"))
  colnames(fit@table) <- NULL
  expect_equal(fit@table[1,1], -0.3499239, tolerance=1e-5)
  expect_equal(fit@table[1,2], 10.04931, tolerance=1e-5)
  expect_equal(fit@table[1,3], 3.897943e-08, tolerance=1e-5)
  expect_equal(fit@table[1,4], 0.00415224, tolerance=1e-5)
  expect_equal(fit@table[1,5], -0.8375797, tolerance=1e-5)
  expect_equal(fit@table[1,6], 0.8700933, tolerance=1e-5)
  colnames(fit@table) <- c("m_cor", "m_logBF", "m_tail_prob", "p_cor", "p_logBF", "p_tail_prob")
  
  # Test deltaOpt
  expect_equal(length(fit@deltaOpt), 1)
  
  # Test alphaOpt
  expect_equal(length(fit@alphaOpt), 1)
  
  # Test dimX
  expect_equal(fit@dimX, c(164, 189))
  
  # Test type
  expect_equal(fit@type, "both")
  
  # Test varlabs
  expect_equal(length(fit@varlabs), 189)
  
  # Test gridAlpha
  expect_equal(dim(fit@gridAlpha), c(100, 3))
  
  # Test valOpt
  expect_equal(length(fit@valOpt), 1)
  
  # Test return.only
  expect_equal(fit@return.only, c("cor","BF","prob"))
  
  
  
  ## Test accessors beam-class object
  # mcor
  expect_equal(dim(mcor(fit)), c(189, 189))
  
  # mcor
  expect_equal(dim(pcor(fit)), c(189, 189))
  
  # postExpOmega
  expect_equal(dim(postExpOmega(fit)), c(189, 189))
  
  # bgraph
  expect_equal(class(bgraph(fit)), "igraph")
  
  # ugraph
  expect_equal(class(ugraph(fit)), "igraph")
  
  # summary
  expect_output(summary(fit))
  
  # print
  expect_output(print(fit))
  
  # show
  expect_output(show(fit))
  
  # Test marg
  expect_equal(is.data.frame(marg(fit)), TRUE)
  
  # Test cond
  expect_equal(is.data.frame(cond(fit)), TRUE)
  
  # Test plotML
  expect_identical(plotML(fit), box())
  
  # Test plotCor
  expect_identical(plotCor(fit, type="both", order = "original", by = 'marginal'), box())
  expect_identical(plotCor(fit, type="both", order = "clust", by = 'marginal'), box())
  expect_identical(plotCor(fit, type="both", order = "clust", by = 'conditional'), box())
  expect_identical(plotCor(fit, type="marginal", order = "original", by = 'marginal'), box())
  expect_identical(plotCor(fit, type="marginal", order = "clust", by = 'marginal'), box())
  expect_identical(plotCor(fit, type="conditional", order = "original", by = 'conditional'), box())
  expect_identical(plotCor(fit, type="conditional", order = "clust", by = 'conditional'), box())
  
  ##################################################
  ## beam selection
  sel <- beam.select(fit, thres=0.01, method = "BH") 
  
  ## Test beam.select-class object
  # Test marginal
  expect_equal(dim(sel@marginal), c(6101, 4))
  expect_equal(colnames(sel@marginal), c("m_cor", "m_logBF", "m_tail_prob", "m_tail_prob_BH"))
  
  # Test conditional
  expect_equal(dim(sel@conditional), c(154, 4))
  expect_equal(colnames(sel@conditional), c("p_cor", "p_logBF", "p_tail_prob", "p_tail_prob_BH"))
  
  # Test method
  expect_equal(sel@method, "BH")
  
  # Test thres
  expect_equal(sel@thres, 0.01)
  
  ## Test accessors beam.select-class object
  # Test marg
  expect_equal(is.data.frame(marg(sel)), TRUE)
  
  # Test cond
  expect_equal(is.data.frame(cond(sel)), TRUE)
  
  # Test mcor
  expect_equal(dim(mcor(sel)), c(189, 189))
  
  # Test pcor
  expect_equal(dim(pcor(sel)), c(189, 189))
  
  # summary
  expect_output(summary(sel))
  
  # print
  expect_output(print(sel))
  
  # show
  expect_output(show(sel))
  
  # Test plotML
  expect_identical(plotML(sel), box())
  
  # Test plotAdj
  expect_identical(plotAdj(sel, type="both", order = "original"), box())
  expect_identical(plotAdj(sel, type="both", order = "clust"), box())
  expect_identical(plotAdj(sel, type="marginal", order = "original"), box())
  expect_identical(plotAdj(sel, type="marginal", order = "clust"), box())
  expect_identical(plotAdj(sel, type="conditional", order = "original"), box())
  expect_identical(plotAdj(sel, type="conditional", order = "clust"), box())
  
  # bgraph
  expect_equal(class(bgraph(sel)), "igraph")
  
  # ugraph
  expect_equal(class(ugraph(sel)), "igraph")
  
  
  ## beam selection
  sel <- beam.select(fit, thres=0.01, method = "HC") 
  
  # Test method
  expect_equal(sel@method, "HC")
  
  # Test marg
  expect_equal(is.data.frame(marg(sel)), TRUE)
  
  # Test cond
  expect_equal(is.data.frame(cond(sel)), TRUE)
  
  # Test mcor
  expect_equal(dim(mcor(sel)), c(189, 189))
  
  # Test pcor
  expect_equal(dim(pcor(sel)), c(189, 189))
  
})



test_that("NlowerthanP", {

  ##################################################
  ## beam estimation
  fit <- beam(X = t(TCPAprad), type="both", verbose=TRUE, D=NULL)
  fit <- beam(X = t(TCPAprad), type="both", verbose=TRUE, D=diag(nrow(TCPAprad)))

  ## Test beam-class object
  # Test table
  expect_equal(dim(fit@table), c(13366, 6))
  expect_equal(colnames(fit@table), c("m_cor", "m_logBF", "m_tail_prob", "p_cor", "p_logBF", "p_tail_prob"))
  colnames(fit@table) <- NULL
  expect_equal(fit@table[1,1], 0.4994693, tolerance=1e-5)
  expect_equal(fit@table[1,2], 26.13615, tolerance=1e-5)
  expect_equal(fit@table[1,3], 7.699685e-16, tolerance=1e-5)
  expect_equal(fit@table[1,4], 0.1178105, tolerance=1e-5)
  expect_equal(fit@table[1,5], 1.506616, tolerance=1e-5)
  expect_equal(fit@table[1,6], 0.0001768815, tolerance=1e-5)
  colnames(fit@table) <- c("m_cor", "m_logBF", "m_tail_prob", "p_cor", "p_logBF", "p_tail_prob")

  # Test deltaOpt
  expect_equal(length(fit@deltaOpt), 1)

  # Test alphaOpt
  expect_equal(length(fit@alphaOpt), 1)

  # Test dimX
  expect_equal(fit@dimX, c(189, 164))

  # Test type
  expect_equal(fit@type, "both")

  # Test varlabs
  expect_equal(length(fit@varlabs), 164)

  # Test gridAlpha
  expect_equal(dim(fit@gridAlpha), c(100, 3))

  # Test valOpt
  expect_equal(length(fit@valOpt), 1)

  # Test return.only
  expect_equal(fit@return.only, c("cor","BF","prob"))



  ## Test accessors beam-class object
  # mcor
  expect_equal(dim(mcor(fit)), c(164, 164))

  # mcor
  expect_equal(dim(pcor(fit)), c(164, 164))

  # bgraph
  expect_equal(class(bgraph(fit)), "igraph")

  # ugraph
  expect_equal(class(ugraph(fit)), "igraph")

  # summary
  expect_output(summary(fit))

  # print
  expect_output(print(fit))

  # show
  expect_output(show(fit))

  # Test marg
  expect_equal(is.data.frame(marg(fit)), TRUE)

  # Test cond
  expect_equal(is.data.frame(cond(fit)), TRUE)

  # Test plotML
  expect_identical(plotML(fit), box())

  # Test plotCor
  expect_identical(plotCor(fit, type="both", order = "original", by = 'marginal'), box())
  expect_identical(plotCor(fit, type="both", order = "clust", by = 'marginal'), box())
  expect_identical(plotCor(fit, type="both", order = "clust", by = 'conditional'), box())
  expect_identical(plotCor(fit, type="marginal", order = "original", by = 'marginal'), box())
  expect_identical(plotCor(fit, type="marginal", order = "clust", by = 'marginal'), box())
  expect_identical(plotCor(fit, type="conditional", order = "original", by = 'conditional'), box())
  expect_identical(plotCor(fit, type="conditional", order = "clust", by = 'conditional'), box())

  ##################################################
  ## beam selection
  sel <- beam.select(fit, thres=0.01, method = "BY")

  ## Test beam.select-class object
  # Test marginal
  expect_equal(dim(sel@marginal), c(4922, 4))
  expect_equal(colnames(sel@marginal), c("m_cor", "m_logBF", "m_tail_prob", "m_tail_prob_BY"))

  # Test conditional
  expect_equal(dim(sel@conditional), c(54, 4))
  expect_equal(colnames(sel@conditional), c("p_cor", "p_logBF", "p_tail_prob", "p_tail_prob_BY"))

  # Test method
  expect_equal(sel@method, "BY")

  # Test thres
  expect_equal(sel@thres, 0.01)

  ## Test accessors beam.select-class object
  # Test marg
  expect_equal(is.data.frame(marg(sel)), TRUE)

  # Test cond
  expect_equal(is.data.frame(cond(sel)), TRUE)

  # Test mcor
  expect_equal(dim(mcor(sel)), c(164, 164))

  # Test pcor
  expect_equal(dim(pcor(sel)), c(164, 164))

  # summary
  expect_output(summary(sel))

  # print
  expect_output(print(sel))

  # show
  expect_output(show(sel))

  # Test plotML
  expect_identical(plotML(sel), box())

  # Test plotAdj
  expect_identical(plotAdj(sel, type="both", order = "original"), box())
  expect_identical(plotAdj(sel, type="both", order = "clust"), box())
  expect_identical(plotAdj(sel, type="marginal", order = "original"), box())
  expect_identical(plotAdj(sel, type="marginal", order = "clust"), box())
  expect_identical(plotAdj(sel, type="conditional", order = "original"), box())
  expect_identical(plotAdj(sel, type="conditional", order = "clust"), box())

  # bgraph
  expect_equal(class(bgraph(sel)), "igraph")

  # ugraph
  expect_equal(class(ugraph(sel)), "igraph")


  ## beam selection
  sel <- beam.select(fit, thres=0.01, method = "holm")

  # Test method
  expect_equal(sel@method, "holm")

  # Test marg
  expect_equal(is.data.frame(marg(sel)), TRUE)

  # Test cond
  expect_equal(is.data.frame(cond(sel)), TRUE)

  # Test mcor
  expect_equal(dim(mcor(sel)), c(164, 164))

  # Test pcor
  expect_equal(dim(pcor(sel)), c(164, 164))

})

