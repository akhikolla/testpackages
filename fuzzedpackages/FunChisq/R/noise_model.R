## Hua created at Mar 3 2014
## Adding noise to contingency table by using House noise model
## House noise model was designed by Yang Zhang and Joe Song
## Note: The adding noise code only supports on parent and one child, does not support combinatorial parents.
##
## Hua modified at Feb 14, 2017
## Updates:
## 1. Parameter "all.tables" must be A list of tables or one table (matrix/data frame).
## 2. u is the noise level between [0, 1].
## 3. Allow adding noise to one margin only (X/Y).
##    Added parameter margin with values between [0,1,2]. Default is 0.
##       0: add noise to both row and column margins of a table;
##       1: add noise to column margin. Column sums are fixed;
##       2: add noise to row margin. Row sums are fixed.
##
## Hua modified at Apr 17, 2017
## Updates:
## 1. Speed-up for house noise model:
##    (1) r*s (r>1 & s>1) tables, Only generate the probability table once rather than for each sample.
##    (2) r*1 and 1*s tables, Directly generate random table based on multinomial distribution with probabilities.
## 2. Added candle noise model.
## 3. Added new function "add.noise", which can specify with noise model to use with parameter "noise.model".
## 4. "add.house.noise" & "add.candle.noise" two function work as well with the same parameters as before.
## 5. Renamed the file name to be "noise_model.R"

add.house.noise <- function(tables, u, margin=0){
  add.noise(tables, u, "house", margin)
}

add.candle.noise <- function(tables, u, margin=0){
  add.noise(tables, u, "candle", margin)
}

add.noise <- function(tables, u, noise.model, margin=0){
  if(is.matrix(tables) || is.data.frame(tables)){
    tables.noised <- add.noise.one.table(tables, u, noise.model, margin)
  }else if(is.list(tables)){
    tables.noised <- lapply(1:length(tables),
                            function(k){
                              t <- tables[[k]]
                              t <- add.noise.one.table(t, u, noise.model, margin)
                              return (t)
                            }
    )
  }else{
    stop("Wrong input format!")
  }
  return (tables.noised)
}

add.noise.one.table <- function(one.table, u, noise.model, margin=0){
  # margin: 0, both X and Y; 1, change along X (in a column); 2, change along Y (in a row).
  tab <- one.table
  table.result <- NULL

  t.row <- nrow(tab)
  t.col <- ncol(tab)

  if(t.row > 1) {
    if(noise.model == "house"){
      col.prob.mat <- house.noise.model.prob.matrix(t.row, u)
    }else if(noise.model == "candle"){
      col.prob.mat <- candle.noise.model.prob.matrix(t.row, u)
    }else{
      stop("Wrong noise.model!")
    }
  }
  if(t.col > 1) {
    if(noise.model == "house"){
      row.prob.mat <- house.noise.model.prob.matrix(t.col, u)
    }else if(noise.model == "candle"){
      row.prob.mat <- candle.noise.model.prob.matrix(t.col, u)
    }else{
      stop("Wrong noise.model!")
    }
  }

  if(margin == 0){
    tables <- lapply(1:(t.row * t.col), function(x) {
      row <- ceiling(x/t.col)
      col <- (x-1)%%t.col+1

      t.one.sample <- matrix(0, nrow=t.row, ncol=t.col)
      if(tab[row, col]==0)return (t.one.sample)

      if(t.row==1 && t.col>1){
        prvector.row <- row.prob.mat$prvector[col,]
        row.index <- sample.int(n=t.col, size=tab[row, col], prob=prvector.row, replace = TRUE)
        t.one.sample[row,] <- as.numeric(table(factor(row.index, levels=c(1:t.col))))
        return(t.one.sample)
      }else if(t.col==1 && t.row>1){
        prvector.col <- col.prob.mat$prvector[row,]
        col.index <- sample.int(n=t.row, size=tab[row, col], prob=prvector.col, replace = TRUE)
        t.one.sample[,col] <- as.numeric(table(factor(col.index, levels=c(1:t.row))))
        return(t.one.sample)
      }else if(t.row>1 && t.col>1){
        valuevector.row <- row.prob.mat$valuevector[col,]
        valuevector.col <- col.prob.mat$valuevector[row,]

        randomNumber <- runif(tab[row, col], 0, 1)
        col.index <- sapply(randomNumber, function(y){
          which(valuevector.row >= y)[1]
        })

        randomNumber <- runif(tab[row, col], 0, 1)
        row.index <- sapply(randomNumber, function(y){
          which(valuevector.col >= y)[1]
        })

        t.one.sample <- Reduce('+', lapply(c(1:tab[row, col]), function(y){
          t.one.sample.tmp <- matrix(0, nrow=t.row, ncol=t.col)
          t.one.sample.tmp[row.index[y], col.index[y]] <- 1
          return(t.one.sample.tmp)
        }))

        return(t.one.sample)
      }else{
        stop("Wrong table size 1*1!")
      }
    })
    table.result <- Reduce('+', tables)
  }else if(margin == 1){
    table.result <- t(apply(tab, 1,function(x){
      t.tmp <- as.matrix(x)
      add.noise.one.table(t.tmp, u, noise.model, 0)
    }))
  }else if(margin == 2){
    table.result <- apply(tab, 2, function(x){
      t.tmp <- as.matrix(x)
      add.noise.one.table(t.tmp, u, noise.model, 0)
    })
  }else{
    stop("Wrong margin values! ([0,1,2])")
  }
  return (table.result)
}


house.noise.model.prob.matrix <- function(baseNum, u){#u: mnoise level
  #Refer to the C++ code in GLN, the code is translated to here in R
  prvector <- matrix (0, nrow=baseNum, ncol=baseNum)
  tempsum <- sapply(1:baseNum, simplify="array", USE.NAMES=FALSE, function(x){
    a <- c(1:baseNum)
    a <- abs(a-x)
    return (sum(a))
  })
  prvector <- t(apply(matrix(1:baseNum, ncol=1), 1, function(x){
    a <- c(1:baseNum)
    res <- (1 - abs(x - a)/tempsum[x]) * u/(baseNum-1)
    res[x] <- res[x]+1-u
    return (res)
  }))
  prvector <- prvector * (1-u) + u/baseNum

  valuevector <- matrix (0, nrow=baseNum, ncol=baseNum)
  valuevector <- apply(matrix(1:baseNum, nrow=1), 2, function(x){
    if(x==1){
      return (prvector[,1:x])
    }else{
      return (rowSums(prvector[,1:x]))
    }

  })
  valuevector[,baseNum] <- 1

  return (list(prvector=prvector, valuevector=valuevector))
}


candle.noise.model.prob.matrix <- function(baseNum, u){#u: mnoise level
  #Refer to the C++ code in GLN, the code is translated to here in R
  prvector <- matrix (u/(baseNum-1), nrow=baseNum, ncol=baseNum)

  for(i in c(1:baseNum)){
    prvector[i,i] <- 1-u
  }

  valuevector <- matrix (0, nrow=baseNum, ncol=baseNum)
  valuevector <- apply(matrix(1:baseNum, nrow=1), 2, function(x){
    if(x==1){
      return (prvector[,1:x])
    }else{
      return (rowSums(prvector[,1:x]))
    }

  })
  valuevector[,baseNum] <- 1

  return (list(prvector=prvector, valuevector=valuevector))
}

####Test: random k 5*5 test cases.
test.noise.model.case <- function(k=3){ # Number of randomed tables
  # A list of k tables
  test <- function(noise.model, k){
    row.num <- sample(x = c(2:5), size = k, replace = TRUE)
    col.num <- sample(x = c(2:5), size = k, replace = TRUE)

    t <- lapply(1:k, function(x){matrix(sample.int(5, row.num[x]*col.num[x], TRUE), row.num[x], col.num[x])})
    t.XY <- add.noise(t, 0.5, noise.model, 0)
    t.X <- add.noise(t, 0.5, noise.model, 1)
    t.Y <- add.noise(t, 0.5, noise.model, 2)

    t.res <- lapply(1:length(t), function(x){
      return(
        sum(t[[x]]) == sum(t.XY[[x]]) &
          all(colSums(t[[x]]) == colSums(t.Y[[x]])) &
          all(rowSums(t[[x]]) == rowSums(t.X[[x]])) &
          length(unique(nrow(t[[x]]), nrow(t.XY[[x]]), nrow(t.X[[x]]), nrow(t.Y[[x]]))) == 1 &
          length(unique(ncol(t[[x]]), ncol(t.XY[[x]]), ncol(t.X[[x]]), ncol(t.Y[[x]]))) == 1
      )
    })
    t.res <- unlist(t.res)


    # One table
    t <- t[[1]]
    t.XY <- add.noise(t, 0.5, noise.model, 0)
    t.X <- add.noise(t, 0.5, noise.model, 1)
    t.Y <- add.noise(t, 0.5, noise.model, 2)

    t.res <- c(t.res,
               sum(t) == sum(t.XY) &
                 all(colSums(t) == colSums(t.Y)) &
                 all(rowSums(t) == rowSums(t.X)) &
                 length(unique(nrow(t), nrow(t.XY), nrow(t.X), nrow(t.Y))) == 1 &
                 length(unique(ncol(t), ncol(t.XY), ncol(t.X), ncol(t.Y))) == 1
               )

    t.res <- unlist(t.res)
    if(all(t.res)){
      message(paste("All test cases passed for noise.model: \"", noise.model, "\"!", sep=""))
    }else{
      message(paste(sum(t.res), '/', length(t.res), " test cases passed for noise.model: \"", noise.model, "\"!", sep=''))
    }
  }
  test("house", k)
  test("candle", k)
}
