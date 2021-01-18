library(benchr)
library(ggplot2)

options(benchr.use_ggplot = TRUE)

b <- benchr::benchmark(1 + 1, 2 + 2)

bp <- benchr:::boxplot_ggplot(b, xlab = "expr", ylab = "time", log = FALSE)
expect_equal(class(bp), c("gg", "ggplot"))
expect_identical(b, bp$data)
expect_equal(bp$labels, list(x = "expr", y = "time"))
expect_true(is(bp$layers[[1]]$geom, "GeomBoxplot"))
expect_true(is(bp$layers[[1]]$stat, "StatBoxplot"))

bpp <- ggplot2::ggplot_build(bp)
expect_equal(bpp$layout$panel_params[[1]]$x$limits, c("1 + 1", "2 + 2"))
expect_equal(bpp$plot$scales$scales[[2]]$trans$name, "identity")


bph <- benchr:::boxplot_ggplot(b, xlab = "expr", ylab = "time", horizontal = TRUE)
expect_equal(class(bph), c("gg", "ggplot"))
expect_identical(b, bph$data)
expect_equal(bph$labels, list(x = "expr", y = "time"))
expect_true(is(bph$layers[[1]]$geom, "GeomBoxplot"))
expect_true(is(bph$layers[[1]]$stat, "StatBoxplot"))

bphp <- ggplot2::ggplot_build(bph)
expect_equal(bphp$layout$panel_params[[1]]$y$limits, c("1 + 1", "2 + 2"))
expect_equal(bphp$plot$scales$scales[[1]]$trans$name, "log-10")


bpv <- benchr:::boxplot_ggplot(b, xlab = "expr", ylab = "time", violin = TRUE)
expect_equal(class(bpv), c("gg", "ggplot"))
expect_identical(b, bpv$data)
expect_equal(bpv$labels, list(x = "expr", y = "time"))
expect_true(is(bpv$layers[[1]]$geom, "GeomViolin"))
expect_true(is(bpv$layers[[1]]$stat, "StatYdensity"))

bpvp <- ggplot2::ggplot_build(bpv)
expect_equal(bpvp$layout$panel_params[[1]]$x$limits, c("1 + 1", "2 + 2"))
expect_equal(bpvp$plot$scales$scales[[1]]$trans$name, "log-10")


pp <- benchr:::plot_ggplot(b, xlab = "replications", ylab = "time")
expect_equal(class(pp), c("gg", "ggplot"))
expect_identical(b, pp$data)
expect_equal(names(pp$labels), c("x", "y", "colour"))
expect_equal(pp$labels, list(x = "replications", y = "time", colour = NULL))
expect_true(is(pp$layers[[1]]$geom, "GeomPoint"))
expect_true(is(pp$layers[[1]]$stat, "StatIdentity"))
expect_equal(pp$layers[[1]]$aes_params, list(shape = 19))

ppp <- ggplot2::ggplot_build(pp)
expect_equal(ppp$plot$scales$scales[[1]]$trans$name, "log-10")
expect_equal(ppp$plot$scales$scales[[2]]$trans$name, "identity")


pp2 <- benchr:::plot_ggplot(b, xlab = "replications", ylab = "time", log = FALSE)
expect_equal(class(pp2), c("gg", "ggplot"))
expect_identical(b, pp2$data)
expect_equal(names(pp2$labels), c("x", "y", "colour"))
expect_equal(pp2$labels, list(x = "replications", y = "time", colour = NULL))
expect_true(is(pp2$layers[[1]]$geom, "GeomPoint"))
expect_true(is(pp2$layers[[1]]$stat, "StatIdentity"))
expect_equal(pp2$layers[[1]]$aes_params, list(shape = 19))

pp2p <- ggplot2::ggplot_build(pp2)
expect_equal(pp2p$plot$scales$scales[[1]]$trans$name, "identity")
expect_equal(pp2p$plot$scales$scales[[2]]$trans$name, "identity")
