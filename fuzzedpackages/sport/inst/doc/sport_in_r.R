## ----eval=FALSE, warning=FALSE , message=FALSE--------------------------------
#  install.packages("sport")
#  devtools::install_github("gogonzo/sport")

## ----echo=TRUE----------------------------------------------------------------
library(sport)
data <- gpheats[1:1002, ]
str(data)

## ----echo=FALSE---------------------------------------------------------------
data[1:8, c("id","rider","rank")]

## ----message=FALSE------------------------------------------------------------
glicko <- glicko_run(formula = rank | id ~ player(rider), data = data)
glicko2 <- glicko2_run(formula = rank | id ~ player(rider), data = data)
bbt <- bbt_run(formula = rank | id ~ player(rider), data = data)
dbl <- dbl_run(formula = rank | id ~ player(rider), data = data)

print(glicko)

## -----------------------------------------------------------------------------
summary(dbl)

## ----message=FALSE, fig.show='hold', out.width = "50%"------------------------
plot(glicko, n = 15)
plot(glicko, players = c("Greg HANCOCK","Tomasz GOLLOB","Tony RICKARDSSON"))

## -----------------------------------------------------------------------------
names(glicko)

## -----------------------------------------------------------------------------
tail(glicko$r)
tail(glicko$pairs)

## -----------------------------------------------------------------------------
glicko2 <- glicko2_run(
  data = data.frame(
    id = c(1, 1, 1, 1),
    player = c("a", "b", "c", "d"),
    rank_player = c(3, 4, 1, 2)
  ), 
  formula = rank_player | id ~ player(player)
)

## -----------------------------------------------------------------------------
glicko2 <- glicko2_run(
  data = data.frame(
    id = c(1, 1, 1, 1),
    team = c("A", "A", "B", "B"),
    player = c("a", "b", "c", "d"),
    rank_team = c(1, 1, 2, 2)
  ), 
  formula = rank_team | id ~ player(player | team)
 )

## -----------------------------------------------------------------------------
dbl <- dbl_run(
  data = data.frame(
    id = c(1, 1, 1, 1),
    name = c("A", "B", "C", "D"),
    rank = c(3, 4, 1, 2),
    gate = c(1, 2, 3, 4),
    factor1 = c("a", "a", "b", "b"),
    factor2 = c("a", "b", "a", "b")
  ), 
  formula = rank | id ~ player(name) + gate * factor1)


## -----------------------------------------------------------------------------
model <- glicko_run(data = gpheats[1:16, ], 
                    rank | id ~ player(rider))

## -----------------------------------------------------------------------------
glicko_run(
  formula = rank | id ~ player(rider), 
  data = gpheats[17:20, ], 
  r    = model$final_r, 
  rd   = model$final_rd
  )$final_r

## -----------------------------------------------------------------------------
library(dplyr)
data <- mutate(data, 
               weight = ifelse(heat >= (max(heat) - 3), 2, 1))

glicko <- glicko_run(formula = rank | id ~ player(rider), 
                     data = data, 
                     weight = "weight")

## -----------------------------------------------------------------------------
bbt1 <- bbt_run(formula = rank | id ~ player(rider),
                data = data, 
                kappa = 0.99) # RD decreases at most 1%

bbt2 <- bbt_run(formula = rank | id ~ player(rider), 
                data = data, 
                kappa = 0.8)  # RD decreases at most 20%

all(bbt1$final_rd >= bbt2$final_rd)

## -----------------------------------------------------------------------------
# bbt example
data <- data %>%
 group_by(rider) %>%
 mutate(idle_30d = if_else(as.integer(date - lag(date)) > 30, 1.0, 2.0)) %>%
 filter(!is.na(idle_30d))

bbt <- bbt_run(rank | id ~ player(rider),
               data = data, 
               lambda = "idle_30d")

## -----------------------------------------------------------------------------
glicko2 <- glicko2_run(
    data = data.frame(
    id = c(1, 1, 1, 1),
    team = c("A", "A", "B", "B"),
    player = c("a", "b", "c", "d"),
    rank_team = c(1, 1, 2, 2),
    share = c(0.4, 0.6, 0.5, 0.5)
  ), 
  formula = rank_team | id ~ player(player | team),
  share = "share"
 )

glicko2$final_r

## -----------------------------------------------------------------------------
glicko2$pairs
glicko2$r

