library(benchr)

b <- benchr::benchmark(1 + 1, times = 1L)
us <- benchr:::convert_units(b, units = "us")
ms <- benchr:::convert_units(b, units = "ms")

expect_equal(attr(us, "units"), "us")
expect_equal(attr(ms, "units"), "ms")
expect_equal(us$time / ms$time, 1000)
