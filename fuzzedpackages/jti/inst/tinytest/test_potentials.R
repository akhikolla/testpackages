el <- matrix(c(
  "A", "T",
  "T", "E",
  "S", "L",
  "S", "B",
  "L", "E",
  "E", "X",
  "E", "D",
  "B", "D"),
  nc = 2,
  byrow = TRUE
)

g  <- igraph::graph_from_edgelist(el)
cp <- compile(cpt_list(asia, g))
j  <- jt(cp)

## ---------------------------------------------------------
## Test that potentials sum to one when propagation has finished        
## ---------------------------------------------------------
lapply(j$charge$C, function(x) expect_equal(sum(x), 1))
lapply(Filter(jti:::neq_null, j$charge$S), function(x) expect_equal(sum(x), 1))
## ---------------------------------------------------------
