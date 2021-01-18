# Test that type = fwd returns a decomposable graph
d1 <- derma[, sample(1:ncol(derma), 10)]
g1 <- fit_graph(d1, type = "fwd", trace = FALSE, q = 0)$G_adj
expect_true(is_decomposable(g1))

# Test that type = tree returns a decomposable graph
d2 <- derma[, sample(1:ncol(derma), 10)]
g2 <- fit_graph(d2, type = "tree", trace = FALSE, q = 0)$G_adj
expect_true(is_decomposable(g2))


# Test that type = bwd returns a decomposable graph
d3<- derma[, sample(1:ncol(derma), 10)]
g3 <- fit_graph(d3, type = "bwd", trace = FALSE, q = 0)$G_adj
expect_true(is_decomposable(g3))
