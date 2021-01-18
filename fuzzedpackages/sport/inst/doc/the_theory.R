## ----message=FALSE, warning=FALSE---------------------------------------------
library(sport)

# example taken from Glickman (1999)
data <- data.frame(id = 1, 
                   name = c("A", "B", "C", "D"), 
                   rank = c(3, 4, 1, 2))
r <- setNames(c(1500, 1400, 1550, 1700), c("A","B","C","D"))
rd <- setNames(c(200, 30, 100, 300), c("A","B","C","D"))

model <- glicko_run(rank | id ~ player(name), data = data, r = r, rd = rd)
print(model$final_r)

## ----message=FALSE, warning=FALSE---------------------------------------------
# example taken from Glickman (2013)
data <- data.frame(id = 1, 
                   name = c("A", "B", "C", "D"), 
                   rank = c(3, 4, 1, 2))
r <- setNames(c(1500, 1400, 1550, 1700), c("A","B","C","D"))
rd <- setNames(c(200, 30, 100, 300), c("A","B","C","D"))

model <- glicko2_run(rank | id ~ player(name), data = data, r = r, rd = rd)
print(model$final_r)


## ----message=FALSE, warning=FALSE---------------------------------------------
data <- data.frame(
  id = c(1, 1, 1, 1),
  team = c("A", "A", "B", "B"),
  player = c("a", "b", "c", "d"),
  rank_player = c(3, 4, 1, 2)
)

model <- bbt_run(
  data = data, 
  formula = rank_player | id ~ player(player),
  r = setNames(c(25, 23.3, 25.83, 28.33), 
               c("a", "b", "c", "d")),
  rd = setNames(c(4.76, 0.71, 2.38, 7.14), 
                c("a", "b", "c", "d"))
)
print(model$final_r)

## ----message=FALSE, warning=FALSE---------------------------------------------
data <- data.frame(
  id = c(1, 1, 1, 1),
  name = c("A", "B", "C", "D"),
  rank = c(3, 4, 1, 2),
  gate = c(1, 2, 3, 4),
  factor1 = c("a", "a", "b", "b"),
  factor2 = c("a", "b", "a", "b")
)

dbl <- dbl_run(
  data = data, 
  formula = rank | id ~ player(name) + gate * factor1)

print(dbl$final_r)

