# Relative D-efficiency with respect to the minimax criterion
meff(formula = ~1/(1 + exp(-b * (x-a))), predvars = "x",
     parvars = c("a", "b"), family = "binomial",
     lp = c(-3, .5), up = c(3, 2),
     x2 = c(-3, -1.608782, 0, 1.608782, 3),
     w2 = c(0.22291601, 0.26438449, 0.02539899, 0.26438449, 0.22291601),
     x1 = c(-1, 1), w1 = c(.5, .5))



# A function to calculate the locally D-optimal design for the 2PL model
Dopt_2pl <- function(a, b){
  x <- c(a + (1/b) * 1.5434046, a - (1/b) * 1.5434046)
  return(list(x = x, w = c(.5, .5)))
}
# Relative D-efficiency with respect to the standardized maximin criterion
meff (formula = ~1/(1 + exp(-b * (x-a))), predvars = "x",
      parvars = c("a", "b"), family = "binomial",
      lp = c(-3, .5), up = c(3, 2),
      x2 = c(-3, -1.611255, 0, 1.611255, 3),
      w2 = c(0.22167034, 0.26592974, 0.02479984, 0.26592974, 0.22167034),
      x1 = c(0, -1), w1 = c(.5, .5),
      standardized = TRUE,
      localdes = Dopt_2pl)


