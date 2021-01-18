######################################################################################################*
######################################################################################################*
#' @importFrom stats runif
GenerateNewCountry <- function(NumOfCountries,
                               lower,
                               upper,
                               sym,
                               sym_point,
                               x_id,
                               w_id,
                               npred,
                               equal_weight,
                               is.only.w){
  # sym: TRUE if a symetric design should be found
  # sym_point: a point that the design is symmetric around
  # npred: number of independnt variables. We must know the dimesnion of the design space.
  # when 'sym = TRUE', you need 'x_id' and 'w_id'. 'x_id' is the locations of the design points and weights
  # for example let k = 5 and  point 3 is the sym_point (known). So we have 4 decision variables for points and 5 decision variables for the weights
  # the indices in this case are x_id <- 1:4 and w_id <- 5:9

  # in the current code x_id and w_id is calculated in updtae_ICA as follows
  ##   if(AlgorithmParams$Symetric)
  ##    x_id <- 1:floor(ProblemParams$k/2) else
  ##       x_id <- 1:(ProblemParams$k * AlgorithmParams$npred)
  ##    w_id <- (x_id[length(x_id)] + 1):length(ProblemParams$VarMinOuter)

  # output is a matrix of the initial solutions. each row is points and weights of the design.
  # each column is one country
  VarMinMatrix <-  matrix(lower, NumOfCountries, length(lower), byrow = TRUE)
  VarMaxMatrix <-  matrix(upper, NumOfCountries, length(upper), byrow = TRUE)
  NewCountry <- (VarMaxMatrix - VarMinMatrix) * matrix(runif (length(VarMaxMatrix)), dim(VarMaxMatrix)[1],  dim(VarMaxMatrix)[2]) + VarMinMatrix
  if (!is.only.w){
    ## we are sorting even a non symetric design!
    npoint <- length(x_id)
    #  sort only if the number of independent variables is 1, otherwise it's a bug!
    if (npred == 1)
      NewCountry[,1:npoint] <- t(apply(NewCountry[,1:npoint, drop = FALSE], 1, sort))
    # we dont sort the weights based on the points becuase it is an inital random country!
    if(!equal_weight){
      w_mat <- NewCountry[, w_id, drop = FALSE]
      # ony for symmetric case is useful!
      if (sym)
        even_odd <- ifelse(length(lower)%% 2 == 0, "even", "odd") else
          even_odd <- NA
        ## the sum of weight will be one!!
        w_sum_to_one <- SumToOne(w_mat = w_mat, sym = sym, even_odd = even_odd)
        NewCountry[, w_id] <- w_sum_to_one
    }
  }else{
    ## here the sum of weights will be equal to 1 for each row (country)
    NewCountry <- SumToOne(w_mat =  NewCountry, sym = sym, even_odd = NA)
  }
  return(NewCountry)
}
######################################################################################################*
######################################################################################################*
CreateInitialEmpires <- function(sorted_Countries,
                                 sorted_Cost,
                                 sorted_InnerParam,
                                 Zeta,
                                 NumOfInitialImperialists,
                                 NumOfAllColonies
){
  # create initial empires, as a list

  # sorted_Countries: matrix of countries, must be sorted while the 'nimp' fisrt of them will be choosen as imperialists
  #  sorted_Cost: corresponding cost of sorted_count
  # sorted_InnerParam: corresponding inner parameters for each countries only for minimax
  # Zeta
  # NumOfInitialImperialists: number of imperialists
  # ncolon: number of colonies
  # output: a list contains positions and costs for the corresponding colonies and imperialists.
  AllImperialistsPosition <- sorted_Countries[1:NumOfInitialImperialists, , drop = FALSE]
  AllImperialistsCost = sorted_Cost[1:NumOfInitialImperialists]
  AllImperialistsInnerParam <- sorted_InnerParam[1:NumOfInitialImperialists, , drop = FALSE]


  ## other colonies, the imperialists are excluded
  AllColoniesPosition = sorted_Countries[-(1:NumOfInitialImperialists), , drop = FALSE]
  AllColoniesCost = sorted_Cost[-c(1:NumOfInitialImperialists)]
  AllColoniesInnerParam= sorted_InnerParam[-(1:NumOfInitialImperialists), , drop = FALSE]


  # calculate imperialists power
  if (max(AllImperialistsCost)>0)
    AllImperialistsPower <- 1.3 * max(AllImperialistsCost) - AllImperialistsCost else
      AllImperialistsPower <- 0.7 * max(AllImperialistsCost) - AllImperialistsCost
  AllImperialistsPower <- AllImperialistsPower/sum(AllImperialistsPower)

  # finding the number of colonies for each empire
  ## number of colonies of each impire. each one should at leats has one colony
  AllImperialistNumOfColonies <- rep(1, length(AllImperialistsPower))
  NumOfRemainColonies <- NumOfAllColonies - sum(AllImperialistNumOfColonies)
  AllImperialistNumOfColonies <- floor(AllImperialistsPower * NumOfRemainColonies)  + AllImperialistNumOfColonies
  ##if still some colonies remain add it to the strongest one
  diff_col <- NumOfAllColonies - sum(AllImperialistNumOfColonies)
  # we add or to the strongest imperialist
  AllImperialistNumOfColonies[1] <-  AllImperialistNumOfColonies[1] + diff_col
  ## modfy the number of colonies for last Impire to be sure that the total number of colonies not less or more than ncolon
  ## warning: bug prone, becasue it may delet all the colonies when it only has one
  last_empire_id <- length(AllImperialistNumOfColonies)
  AllImperialistNumOfColonies[last_empire_id] <- NumOfAllColonies - sum(AllImperialistNumOfColonies[-last_empire_id])



  RandomIndex <- sample(x = NumOfAllColonies, size = NumOfAllColonies, replace = FALSE)
  RandonIndexList <- vector("list", NumOfInitialImperialists)
  InitialInd <- 1
  LastInd <- AllImperialistNumOfColonies[1]

  for(jj in 1:NumOfInitialImperialists){
    if( InitialInd > length(RandomIndex))
      RandonIndexList [[jj]] <- 0 else
        RandonIndexList [[jj]] <- RandomIndex[InitialInd:LastInd]
      InitialInd  <- 1 + LastInd
      LastInd <- LastInd +   AllImperialistNumOfColonies[jj + 1]
  }

  ##now we create the empires. we sent the Imperialist  position and cost and to a list
  Empires <- vector("list", NumOfInitialImperialists)
  for(i in 1:length(AllImperialistNumOfColonies)){
    Empires[[i]]$ImperialistPosition <- AllImperialistsPosition[i, , drop = FALSE]
    Empires[[i]]$ImperialistCost <- AllImperialistsCost[i]
    Empires[[i]]$ImperialistInnerParam <- AllImperialistsInnerParam[i, , drop = FALSE]
    if(all(RandonIndexList[[i]]== 0)){
      stop("One initial empire has no colony.\nIncrease the number of colonies or decrease the number of imperialists.\nUsually the number of imperialists set to be 10 percent of the colonies.")
    }else{
      Empires[[i]]$ColoniesPosition <- AllColoniesPosition[RandonIndexList[[i]], , drop = FALSE]
      Empires[[i]]$ColoniesCost <- AllColoniesCost[RandonIndexList[[i]]]
      Empires[[i]]$ColoniesInnerParam <- AllColoniesInnerParam[RandonIndexList[[i]], , drop = FALSE]
      Empires[[i]]$TotalCost <- Empires[[i]]$ImperialistCost + Zeta * mean(Empires[[i]]$ColoniesCost)
    }

  }

  return(Empires)

}

# CreateInitialEmpires_bayes <- function(sorted_Countries,
#                                        sorted_Cost,
#                                        Zeta,
#                                        NumOfInitialImperialists,
#                                        NumOfAllColonies
# ){
#   # create initial empires, as a list
#
#   # sorted_Countries: matrix of countries, must be sorted while the 'nimp' fisrt of them will be choosen as imperialists
#   # sorted_Cost: corresponding cost of sorted_count
#   # Zeta
#   # NumOfInitialImperialists: number of imperialists
#   # ncolon: number of colonies
#   # output: a list contains positions and costs for the corresponding colonies and imperialists.
#   AllImperialistsPosition <- sorted_Countries[1:NumOfInitialImperialists, , drop = FALSE]
#   AllImperialistsCost = sorted_Cost[1:NumOfInitialImperialists]
#   ## other colonies, the imperialists are excluded
#   AllColoniesPosition = sorted_Countries[-(1:NumOfInitialImperialists), , drop = FALSE]
#   AllColoniesCost = sorted_Cost[-c(1:NumOfInitialImperialists)]
#   # calculate imperialists power
#   if (max(AllImperialistsCost)>0)
#     AllImperialistsPower <- 1.3 * max(AllImperialistsCost) - AllImperialistsCost else
#       AllImperialistsPower <- 0.7 * max(AllImperialistsCost) - AllImperialistsCost
#   AllImperialistsPower <- AllImperialistsPower/sum(AllImperialistsPower)
#   # finding the number of colonies for each empire
#   # number of colonies of each impire. each one should at leats has one colony
#   AllImperialistNumOfColonies <- rep(1, length(AllImperialistsPower))
#   NumOfRemainColonies <- NumOfAllColonies - sum(AllImperialistNumOfColonies)
#   AllImperialistNumOfColonies <- floor(AllImperialistsPower * NumOfRemainColonies)  + AllImperialistNumOfColonies
#   ##if still some colonies remain add them to the strongest one
#   diff_col <- NumOfAllColonies - sum(AllImperialistNumOfColonies)
#   # we add them to the strongest imperialist
#   AllImperialistNumOfColonies[1] <-  AllImperialistNumOfColonies[1] + diff_col
#   ## modfy the number of colonies for the last Impire to be sure that the total number of colonies not less or more than number of colonies
#   ## warning: bug prone, becasue it may delet all the colonies when it only has one
#   last_empire_id <- length(AllImperialistNumOfColonies)
#   AllImperialistNumOfColonies[last_empire_id] <- NumOfAllColonies - sum(AllImperialistNumOfColonies[-last_empire_id])
#   RandomIndex <- sample(x = NumOfAllColonies, size = NumOfAllColonies, replace = FALSE)
#   RandonIndexList <- vector("list", NumOfInitialImperialists)
#   InitialInd <- 1
#   LastInd <- AllImperialistNumOfColonies[1]
#   for(jj in 1:NumOfInitialImperialists){
#     if( InitialInd > length(RandomIndex))
#       RandonIndexList [[jj]] <- 0 else
#         RandonIndexList [[jj]] <- RandomIndex[InitialInd:LastInd]
#       InitialInd  <- 1 + LastInd
#       LastInd <- LastInd +   AllImperialistNumOfColonies[jj + 1]
#   }
#   # now we create the empires. we send the Imperialist position and cost and to a list
#   Empires <- vector("list", NumOfInitialImperialists)
#   for(i in 1:length(AllImperialistNumOfColonies)){
#     Empires[[i]]$ImperialistPosition <- AllImperialistsPosition[i, , drop = FALSE]
#     Empires[[i]]$ImperialistCost <- AllImperialistsCost[i]
#     if(all(RandonIndexList[[i]]== 0)){
#       stop("One initial empire has no colony.\nIncrease the number of colonies or decrease the number of imperialists.\nUsually the number of imperialists set to be 10 percent of the colonies.")
#     }else{
#       Empires[[i]]$ColoniesPosition <- AllColoniesPosition[RandonIndexList[[i]], , drop = FALSE]
#       Empires[[i]]$ColoniesCost <- AllColoniesCost[RandonIndexList[[i]]]
#       Empires[[i]]$TotalCost <- Empires[[i]]$ImperialistCost + Zeta * mean(Empires[[i]]$ColoniesCost)
#     }
#   }
#   return(Empires)
#}
######################################################################################################*
######################################################################################################*
LocalSearch <- function(TheEmpire, lower, upper, l, fixed_arg = fixed_arg){
  # the local search is done for the half of the positions for 'l' times
  # when the design points are given, then the positions are only weights
  nfeval <- 0
  n_success <- 0
  imperialist <- as.vector(TheEmpire$ImperialistPosition)
  for(k in 1:(length(imperialist)/2)){
    counter <- 1
    while(counter <= l){
      diff_vec <- -sweep(x = TheEmpire$ColoniesPosition, MARGIN = 2, STATS = TheEmpire$ImperialistPosition, FUN = "-")
      ## 'd' is the maximal range is set to the distance between the imperialist and its closest colony in
      # the same emipire divided by  the square root of the number of veriables
      d <- apply(diff_vec, 1, function(x)sqrt(sum(x^2)))
      d <- min(d)/sqrt(length(imperialist)/2) #because the local search is only done for the half of the dimension
      if (round(d, 8) == 0)
        d <- .05 ## we set the d be equal to .05 if d is eqaul to zero!
      ## because if d is equal to zero then the local search is not useful anymore!!!
      NewPos <-  TheEmpire$ImperialistPosition
      lambda <- runif(1, -1, 1)
      NewPos[k] <-  NewPos[k] + lambda * d * .05
      NewPos[k] <- (NewPos[k] <= upper[k] & NewPos[k] >= lower[k]) * NewPos[k] +
        (NewPos[k] > upper[k]) * (upper[k] - .25 * (upper[k] - lower[k]) * runif(1)) +
        (NewPos[k] < lower[k]) * (lower[k] + .25 * (upper[k] - lower[k]) * runif(1))
      if (fixed_arg$is.only.w)
        NewPos <- NewPos/sum(NewPos)
      output <- fixed_arg$Calculate_Cost(mat = matrix(NewPos, nrow = 1), fixed_arg = fixed_arg)
      NewPos_cost <- output$cost
      fneval_candidate <- output$nfeval


      ### here make it with tolerance
      if(NewPos_cost < TheEmpire$ImperialistCost){
        n_success <- n_success + 1
        nfeval <- nfeval + fneval_candidate
        TheEmpire$ImperialistPosition <- matrix(NewPos, 1, length(lower))
        TheEmpire$ImperialistCost <- NewPos_cost
        TheEmpire$ImperialistInnerParam <- output$inner_optima
        counter <- l

      }
      counter <- counter + 1
    }
  }
  return(list(TheEmpire = TheEmpire, nfeval = nfeval, n_success = n_success))
}
######################################################################################################*
######################################################################################################*
AssimilateColonies2 <- function(TheEmpire,
                                AssimilationCoefficient,
                                VarMin,
                                VarMax,
                                ExceedStrategy,
                                AsssimilationStrategy,
                                sym,
                                MoveOnlyWhenImprove,
                                fixed_arg,
                                w_id,
                                equal_weight){
  #   TheEmpire: a list  contains the position of the colonies and imperialists and the cost values
  #   AssimilationCoefficient
  #   VarMin: lower bound of the positions
  #   VarMax: upper bounds of the positions
  #   ExceedStrategy: passed to 'CheckBoxConstraints' function. can be 'random', 'perturbed'
  #   AsssimilationStrategy: 'PICA' or original 'ICA' or 'OICA'
  #   sym: if the positions must be symetric
  #   MoveOnlyWhenImprove.
  #   fixed_arg. arguments that are fixed and will be passed to the optimization over the inner problem, passed to Calculate cost
  #   w_id: column index of weights in the position matrix. you can take w_id from fixed_arg, but is avoided to not make any confusion
  # if  ExceedStrategy = Random, the points that exceed lower and upper bound take the random values between the bounds.
  # if ExceedStrategy = Lower_Upper, the points that exceed lower nad upper bound take the lower and upper bounds.

  # Rerurn the assimilated positions of colonies and update  TheEmpire$ColoniesPosition, TheEmpire$ColoniesCost and   TheEmpire$ColoniesInnerParam

  NumOfColonies <- dim(TheEmpire$ColoniesPosition)[1]
  ##MARGIN is 2 since it should be subtracted columnwise!
  diff_vec <- -sweep(x = TheEmpire$ColoniesPosition, MARGIN = 2, STATS = TheEmpire$ImperialistPosition, FUN = "-")


  if(AsssimilationStrategy =="ICA")
    NewPosition = TheEmpire$ColoniesPosition +
    AssimilationCoefficient * matrix(runif(length(diff_vec)), dim(diff_vec)[1], dim(diff_vec)[2]) * diff_vec

  if(AsssimilationStrategy == "PICA")
    NewPosition = TheEmpire$ColoniesPosition +
    (AssimilationCoefficient * matrix(runif(length(diff_vec)), dim(diff_vec)[1], dim(diff_vec)[2]) - 1) * diff_vec


  NewPosition <- CheckBoxConstraints(Position = NewPosition, Lower = VarMin, Upper = VarMax, Strategy = "perturbed")


  ### for weights position
  if (!equal_weight){
    ##if you want to only move in feasible region!
    w_mat <- NewPosition[, w_id, drop = FALSE]
    even_odd <- ifelse(length(VarMin)%% 2 == 0, "even", "odd")
    NewPosition[, w_id] <- SumToOne(w_mat = w_mat, sym = sym, even_odd = even_odd)
  }
  temp <- fixed_arg$Calculate_Cost(mat=NewPosition, fixed_arg=fixed_arg)
  nfeval <- temp$nfeval
  NewCost <- temp$cost
  nimprove <- 0 ## if there is no improvement this will be returnEd!!

  if(MoveOnlyWhenImprove){
    OldCost <- TheEmpire$ColoniesCost
    ### here make it with tolerance
    replace_id <- which(NewCost < OldCost)
    if(length(replace_id) != 0){
      TheEmpire$ColoniesPosition[replace_id,] <- NewPosition[replace_id, , drop = FALSE]
      TheEmpire$ColoniesCost[replace_id] <- NewCost[replace_id]
      TheEmpire$ColoniesInnerParam[replace_id, ] <- temp$inner_optima[replace_id, , drop = FALSE]
      nimprove <- length(replace_id)
    }##else dont change anything
  }
  if(!MoveOnlyWhenImprove){
    TheEmpire$ColoniesPosition <- NewPosition
    TheEmpire$ColoniesCost <- NewCost
    TheEmpire$ColoniesInnerParam <- temp$inner_optima
  }


  return(list(TheEmpire = TheEmpire, nfeval = nfeval, nimprove = nimprove))

}
######################################################################################################*
######################################################################################################*
ImperialisticCompetition <- function(Empires, Zeta){
  # Empires: all empires
  # perform an imperialistic competition
  # Zeta: a cofficient needed to compute the total cost
  # updating 'ColoniesPosition', "ColoniesCost"
  NumOfEmpires <- length(Empires)
  # competition pressure
  if (NumOfEmpires == 1 || runif(1)>.11){
    return(Empires)
  }else{
    #we first select the the empire that possess the weakest colony from the weakest empire.
    TotalCosts <- sapply(Empires, "[[", "TotalCost", simplify =TRUE )
    # maybe you should move the following syntax in a better position!
    if (any(is.infinite(TotalCosts)))
      stop("Please check the inputs of your design problem.\n  The value of the criterion is not finite for some imperialists from the design space. It may be a sign that the FIM for the given design setting is (computationaly) singular.")
    WeakestEmpireInd <- which.max(TotalCosts)

    TotalPowers <- TotalCosts[WeakestEmpireInd] - TotalCosts
    if (all(TotalPowers == 0))
      return(Empires)
    PossessionProbability <- TotalPowers/sum(TotalPowers)

    SelectedEmpireInd <- sample(1:NumOfEmpires, 1, prob = PossessionProbability)

    ## nn is the number of colonies in the weakest Empire. we should choose a colony by random!
    nn <- length(Empires[[WeakestEmpireInd]]$ColoniesCost)
    jj <- sample(1:nn,1)
    ## jj is the index of the colony that should be moved!

    # add the colony to the selected empire
    Empires[[SelectedEmpireInd]]$ColoniesPosition <- rbind(Empires[[SelectedEmpireInd]]$ColoniesPosition,
                                                           Empires[[WeakestEmpireInd]]$ColoniesPosition[jj, ])

    Empires[[SelectedEmpireInd]]$ColoniesCost <- c(Empires[[SelectedEmpireInd]]$ColoniesCost,
                                                   Empires[[WeakestEmpireInd]]$ColoniesCost[jj])


    Empires[[SelectedEmpireInd]]$ColoniesInnerParam <- rbind(Empires[[SelectedEmpireInd]]$ColoniesInnerParam,
                                                             Empires[[WeakestEmpireInd]]$ColoniesInnerParam[jj, ])
    # now we should remove the colony from the weakest empire (the jjth colony must be removed).
    Empires[[WeakestEmpireInd]]$ColoniesPosition <- Empires[[WeakestEmpireInd]]$ColoniesPosition[-jj, , drop = FALSE]
    Empires[[WeakestEmpireInd]]$ColoniesCost <- Empires[[WeakestEmpireInd]]$ColoniesCost[-jj]

    # collapsing the weakest eimpire when it only has one colony
    # if the weakest empire has only one colony it should be collapsed and be combined with the strongest empire
    nn <- length(Empires[[WeakestEmpireInd]]$ColoniesCost)
    if (nn<=1){
      Empires[[SelectedEmpireInd]]$ColoniesPosition <- rbind(Empires[[SelectedEmpireInd]]$ColoniesPosition,
                                                             Empires[[WeakestEmpireInd]]$ImperialistPosition)

      Empires[[SelectedEmpireInd]]$ColoniesCost <- c(Empires[[SelectedEmpireInd]]$ColoniesCost,
                                                     Empires[[WeakestEmpireInd]]$ImperialistCost)

      Empires[[SelectedEmpireInd]]$ColoniesInnerParam <- rbind(Empires[[SelectedEmpireInd]]$ColoniesInnerParam,
                                                               Empires[[WeakestEmpireInd]]$ColoniesInnerParam)
      Empires[[SelectedEmpireInd]]$TotalCost <- Empires[[SelectedEmpireInd]]$ImperialistCost + Zeta * mean(Empires[[SelectedEmpireInd]]$ColoniesCost)
      # remove the empire from the list (collpase the empire because it has only one colony)
      Empires[[WeakestEmpireInd]] <- NULL
    }
    return(Empires)
  }
}

######################################################################################################*
######################################################################################################*
PossesEmpire <- function(TheEmpire){
  ## Exchanging the position of the Imperialist and  the best Colony, if the best colony is better
  ColoniesCost <- TheEmpire$ColoniesCost
  BestColonyInd <- which.min(ColoniesCost)
  if(ColoniesCost[BestColonyInd] < TheEmpire$ImperialistCost){
    OldImperialistPosition <- TheEmpire$ImperialistPosition
    OldImperialistCost <- TheEmpire$ImperialistCost
    OldImperialistInnerParam <- TheEmpire$ImperialistInnerParam

    TheEmpire$ImperialistPosition <- TheEmpire$ColoniesPosition[BestColonyInd, ,drop = FALSE]
    TheEmpire$ImperialistCost <- TheEmpire$ColoniesCost[BestColonyInd]
    TheEmpire$ImperialistInnerParam <- TheEmpire$ColoniesInnerParam[BestColonyInd, ,drop = FALSE]

    TheEmpire$ColoniesPosition[BestColonyInd,] <- OldImperialistPosition
    TheEmpire$ColoniesCost[BestColonyInd] <- OldImperialistCost
    TheEmpire$ColoniesInnerParam[BestColonyInd,] <- OldImperialistInnerParam
  }
  return(TheEmpire)
}
######################################################################################################*
######################################################################################################*
RevolveColonies <- function(TheEmpire, RevolutionRate, NumOfCountries,
                            lower,
                            upper,
                            sym,
                            sym_point,
                            fixed_arg,
                            w_id,
                            equal_weight){
  # Revolve the colonies and update 'TheEmpire'
  # the position is changed if and only if improvement happens
  nrevol <- 0 ## number of revolutions
  nfeval <- 0
  # the number of design points is odd or even
  even_odd <- ifelse(length(lower)%% 2 == 0, "even", "odd")
  sigma <- 0.1 * (upper - lower)
  n_col <- dim(TheEmpire$ImperialistPosition)[2]

  ## we choose the columns by random
  DimSize <- sample(x = 1:n_col, size = 1)
  DimIndex <- sample(x = 1:n_col, size = DimSize )
  new_imp <- TheEmpire$ImperialistPosition

  new_imp[, DimIndex]  <- new_imp[, DimIndex, drop = FALSE] + sigma[DimIndex] * runif(DimSize, -1, 1)

  new_imp[, DimIndex] <- new_imp[, DimIndex]  * (( new_imp [, DimIndex] <= upper[DimIndex]) & ( new_imp[, DimIndex] >= lower[DimIndex])) +
    (new_imp[, DimIndex] > upper[DimIndex]) * (upper[DimIndex] - .25 * (upper[DimIndex] - lower[DimIndex]) * runif(1)) +
    (new_imp[, DimIndex] < lower[DimIndex]) * (lower[DimIndex] + .25 * (upper[DimIndex] - lower[DimIndex]) * runif(1))
  if(!equal_weight)
    new_imp[, w_id] <- SumToOne(w_mat =  new_imp[, w_id, drop = FALSE], sym=sym, even_odd=even_odd)
  temp <- fixed_arg$Calculate_Cost(mat = new_imp, fixed_arg = fixed_arg)
  nfeval <- nfeval + temp$nfeval
  new_imp_cost <- temp$cost
  ## here make it with tolerance
  if(new_imp_cost < TheEmpire$ImperialistCost){
    nrevol <- nrevol + 1
    TheEmpire$ImperialistCost <- new_imp_cost
    TheEmpire$ImperialistPosition <- new_imp
    TheEmpire$ImperialistInnerParam <- temp$inner_optima ## new_imp_inner_param
  }

  # now we revolve the colonies
  NumOfRevolvingColonies <- round(RevolutionRate * length(TheEmpire$ColoniesCost))

  if(NumOfRevolvingColonies != 0){


    # choose the colonies that should be revolved.
    RandomIndex <- sample(x = 1:length(TheEmpire$ColoniesCost),
                          size = length(TheEmpire$ColoniesCost),
                          replace = FALSE)[1:NumOfRevolvingColonies]

    ######################################################*
    ### a local search for revolved colonies
    DimSize <- sample(x = 1:n_col, size = 1)
    DimIndex <- sample(x = 1:n_col, size = DimSize )

    new_colonies <- TheEmpire$ColoniesPosition[RandomIndex, , drop = FALSE]
    new_colonies[, DimIndex] <-  new_colonies[, DimIndex, drop = FALSE] +
      matrix(sigma[DimIndex], nrow = length(RandomIndex), ncol = DimSize, byrow = TRUE) * matrix(runif(length(new_colonies[, DimIndex, drop = FALSE])), length(  RandomIndex) , DimSize )
    if(!equal_weight)
      new_colonies[,w_id] <- SumToOne(w_mat = new_colonies[, w_id, drop = FALSE], sym=sym, even_odd=even_odd)
    # check the lower upper bound

    new_colonies[, DimIndex]<- CheckBoxConstraints(Position = new_colonies[, DimIndex, drop = FALSE],
                                                   Lower = lower[DimIndex] , Upper = upper[DimIndex], Strategy = "perturbed")

    temp1 <- fixed_arg$Calculate_Cost(mat = new_colonies, fixed_arg = fixed_arg)
    nfeval <- nfeval + temp$nfeval
    new_colonies_cost <- temp1$cost

    change_id <-  which(new_colonies_cost < TheEmpire$ColoniesCost[RandomIndex])
    if(length(change_id) != 0){
      nrevol <- nrevol + length(change_id)
      TheEmpire$ColoniesPosition[RandomIndex[change_id], ] <- new_colonies[change_id, , drop = FALSE]
      TheEmpire$ColoniesCost[RandomIndex[change_id]] <- new_colonies_cost[change_id]
      TheEmpire$ColoniesInnerParam[RandomIndex[change_id], ] <- temp1$inner_optima[change_id, , drop = FALSE]
    }
  }
  return(list(TheEmpire = TheEmpire, nfeval = nfeval, nrevol = nrevol))
}

######################################################################################################*
######################################################################################################*
UniteSimilarEmpires <- function(Empires, Zeta, UnitingThreshold, SearchSpaceSize){
  # Unite the similar empires
  # it can be removed

  TheresholdDistance <- UnitingThreshold * base::norm(as.matrix(SearchSpaceSize), type = '2')
  NumOfEmpires <- length(Empires)

  # because the length(Empire) = NumOfEmpires is dynamic (it will be decreased if uniting happens we should use while)
  counter1 <- 1
  counter2 <- 2

  while(counter1 <= (NumOfEmpires-1)){
    while(counter2 <= NumOfEmpires){
      DistanceVector <- Empires[[counter1]]$ImperialistPosition - Empires[[counter2]]$ImperialistPosition
      Distance <- base::norm(DistanceVector, "2")

      if(Distance <= TheresholdDistance){
        if (Empires[[counter1]]$ImperialistCost < Empires[[counter2]]$ImperialistCost){
          BetterEmpireInd <- counter1
          WorseEmpireInd <- counter2} else{
            BetterEmpireInd <- counter2
            WorseEmpireInd <- counter1
          }

        Empires[[BetterEmpireInd]]$ColoniesPosition <- rbind(Empires[[BetterEmpireInd]]$ColoniesPosition,
                                                             Empires[[WorseEmpireInd]]$ImperialistPosition,
                                                             Empires[[WorseEmpireInd]]$ColoniesPosition)

        Empires[[BetterEmpireInd]]$ColoniesCost <- c(Empires[[BetterEmpireInd]]$ColoniesCost,
                                                     Empires[[WorseEmpireInd]]$ImperialistCost,
                                                     Empires[[WorseEmpireInd]]$ColoniesCost)

        Empires[[BetterEmpireInd]]$ColoniesInnerParam <- rbind(Empires[[BetterEmpireInd]]$ColoniesInnerParam,
                                                               Empires[[WorseEmpireInd]]$ImperialistInnerParam,
                                                               Empires[[WorseEmpireInd]]$ColoniesInnerParam)
        # Update TotalCost for new United Empire
        Empires[[BetterEmpireInd]]$TotalCost = Empires[[BetterEmpireInd]]$ImperialistCost + Zeta * mean(Empires[[BetterEmpireInd]]$ColoniesCost)
        Empires[[WorseEmpireInd]] <- NULL
        NumOfEmpires <- NumOfEmpires - 1
      }
      counter2 <- counter2 + 1
    }
    counter1 <- counter1 + 1
  }
  return(Empires)
}
######################################################################################################*
######################################################################################################*
CheckBoxConstraints <- function(Position, Lower, Upper, Strategy){
  # Position: matrix of positions
  # Lower: vector of lower bound corresponding to each row
  # Upper: upper bound
  # Startegy: character can be 'perturbed' or 'random'
  # output: the corrected positions
  Position <- Position
  n_row <- dim(Position)[1]

  if(Strategy == "perturbed"){
    ExceedLower <- t(sapply(1:n_row, function(j) Position[j,, drop = FALSE] < Lower))
    if(dim(Position)[2] == 1) # because it gives us wrong answer when we have a matrix whit one column
      ExceedLower <- t(ExceedLower)

    if(any(ExceedLower)){
      ExceedLowerIndex <- which(ExceedLower, arr.ind = TRUE)
      L <-  Lower[ExceedLowerIndex[, 2]]
      U <- Upper[ExceedLowerIndex[, 2]]
      x <- Position[ExceedLowerIndex]
      Position[ExceedLowerIndex] <- 2*L - (x + floor((L-x)/(U-L))*(U-L))
      L <- U <- x <- NA
    }

    ExceedUpper <- t(sapply(1:n_row, function(j) Position[j, , drop = FALSE] > Upper))
    if(dim(Position)[2] == 1) # because it gives us wrong answer when we have a matrix whit one column
      ExceedUpper <- t(ExceedUpper)
    if(any(ExceedUpper)) {
      ExceedUpperIndex <- which(ExceedUpper, arr.ind = TRUE)
      L <-  Lower[ExceedUpperIndex[, 2]]
      U <- Upper[ExceedUpperIndex[, 2]]
      x <- Position[ExceedUpperIndex]
      Position[ExceedUpperIndex] <- 2*U - (x - floor((x-U)/(U-L))*(U-L))
    }
  }

  if(Strategy == "random"){
    ExceedLower <- t(sapply(1:n_row, function(j) Position[j, , drop = FALSE] < Lower))
    if(dim(Position)[2] == 1) ##because it gives us wrong answer when we have a matrix whit one colmun
      ExceedLower <- t(ExceedLower)
    if(any(ExceedLower)) {
      ExceedLowerIndex <- which(ExceedLower, arr.ind = TRUE)
      Position[ExceedLowerIndex] <- Lower[ExceedLowerIndex[, 2]] + runif(dim(ExceedLowerIndex)[1])
    }

    ExceedUpper <- t(sapply(1:n_row, function(j) Position[j, , drop = FALSE] > Upper))
    if(dim(Position)[2] == 1) ##because it gives us wrong answer when we have a matrix whit one colmun
      ExceedUpper <- t(ExceedUpper)
    if(any(ExceedUpper)) {
      ExceedUpperIndex <- which(ExceedUpper, arr.ind = TRUE)
      Position[ExceedUpperIndex] <- Upper[ExceedUpperIndex[, 2]] - runif(dim(ExceedUpperIndex)[1])
    }
  }
  return(Position)
}
######################################################################################################*
######################################################################################################*
find_on_points <- function(fn, ..., points){
  # ... is the further arguments to be passed
  # x and w is the arguments that will be passed
  # output calculate fn(points)
  if(!is.function(fn))
    stop("'fn' must be a function.")


  cost <- sapply(1:dim(points)[1], function(j) fn(points[j, ],...))
  #cost <- apply(points, 1, fn,...)
  points <- cbind(points, cost, deparse.level = 0)


  # warnings:this is not minima!
  return(list(minima = points, counts = nrow(points)))
}
######################################################################################################*
######################################################################################################*
make_vertices <- function(lower, upper){
  ## give the lower and upper of region of uncertainty and the output is the vertices. each row is a vertex
  par_list <- vector("list",length(lower))
  for(i in 1:length(lower))
    par_list[[i]] <- c(lower[i], upper[i])
  vertices <- matrix(unlist(expand.grid(par_list)), nrow = 2^length(upper))
  return(vertices)
}
######################################################################################################*
######################################################################################################*
SumToOne <- function(w_mat, sym, even_odd){
  #even_odd is NA when sym == FALSE, otherwise is a character string "even" or "odd"
  # even_odd: is the number of support points is odd or even. needed when sym = TRUE
  # because when sym = TRUE the sum ofweights must be one, while they are symmetric
  if(!sym)
    w_mat_out <- w_mat/apply(w_mat, 1, sum) else{
      n_col <- dim(w_mat)[2]

      if(even_odd == "even"){
        all_w <- cbind(w_mat, w_mat[, rev(1:n_col), drop = FALSE])
        all_w <- all_w/apply( all_w, 1, sum)
        w_mat_out <- all_w[, 1:n_col, drop = FALSE]
      }
      if(even_odd == "odd"){
        sym_w <- w_mat[,n_col, drop = FALSE]
        ##now we remove the symetric weight!
        w_mat <- w_mat[, -n_col, drop = FALSE]
        all_w <- cbind(w_mat, w_mat[, rev(1:dim(w_mat)[2]), drop = FALSE])
        all_w <- cbind(all_w, sym_w)
        all_w <- all_w/apply( all_w, 1, sum)
        w_mat_out <- cbind(all_w[,1:(dim(all_w)[2]/2), drop =FALSE], all_w[, dim(all_w)[2], drop = FALSE])
      }
    }
  return(w_mat_out)
}
######################################################################################################*
######################################################################################################*
ICA_extract_x_w <- function(x, w, sym, sym_point){
  # when the deign is symetrix, this function extract x and w.
  # x and w are the decision variables while the output are the real design points and weights
  include_symetric <- length(x) != length(w)
  x <- c(x, rev(2 * sym_point - x))
  if(include_symetric)
    x <- c(x, sym_point)
  if(include_symetric){
    w_sym <- w[length(w)]
    w <- w[-length(w)]
    w <- c(w, rev(w))
    w <- c(w, w_sym)
  }else
    w <- c(w, rev(w))
  return(list(x=x, w=w))
}
#############################################################################################################*
#############################################################################################################*
add_fixed_ICA.control <- function(ICA.control.list){
  ## the followings are not fixed and not controled by the user!
  ICA.control.list$zeta <- .1
  ICA.control.list$l <-  2
  ICA.control.list$assim_strategy <- "PICA"
  ICA.control.list$lsearch <- 2
  ICA.control.list$only_improve <- TRUE
  return(ICA.control.list)
}
######################################################################################################*
######################################################################################################*
#' @importFrom graphics abline axis legend lines mtext plot points
PlotEffCost <- function(from,
                        to ,
                        AllCost, ##all criterion up to now (all costfunction)
                        UserCost,
                        title1,
                        DesignType,
                        plot_main = TRUE,
                        ...){

  ## plot the evolution per each iteration
  # from: start iteration
  # to: last iteration
  # AllCost: all the cost values that should be plotted
  # title1: what title do u want? for example ICA, FWA or multiple
  # DesignType: 'standardized' or 'minimax'or 'locally', standardized is not available anymore
  # UserCost: the cost value that user has. if not NULL, then the relative efficiency is plotted.
  cex.main <- .9
  prec <- 8
  if (is.null(UserCost)){
    ##cost function plot
    if (plot_main)
      main1 <- paste(title1, ": ", round( AllCost[length(AllCost)], prec), paste = "") else
        main1 = NULL

      plot(x = from:to, y = AllCost ,
           xlim = c(from, to),  ylab = "Imperialist Cost", type = "s",
           main = main1,
           cex.main = cex.main,
           xaxt = "n",...)

  }else{
    ##Efficiency
    Efficiency <- switch(DesignType, "minimax" = UserCost/AllCost, "standardized" = AllCost/UserCost, "locally" =  UserCost/AllCost)
    if (plot_main)
      main1 <- paste(title1, ": ", round(Efficiency[length(Efficiency)], prec), paste = "") else
        main1 = NULL

      plot(x = from:to, y = Efficiency ,
           xlim = c(from, to),  ylab = "Efficiency", type = "s",
           main = main1,
           cex.main = cex.main,
           xaxt = "n",...)

  }
  ## here we plot the axis to control everything.
  if (to < 5)
    axis(side = 1, at = c(from:to),
         labels = c(from:to)) else{
           pretty_plot <- pretty(from:to)
           axis(side = 1, at = pretty_plot,
                labels = pretty_plot)
         }
}

######################################################################################################*
######################################################################################################*
#' @importFrom stats optim
optim2 <- function(fn, ..., lower , upper, control = list(factr = sqrt(.Machine$double.eps)), n_seg){
  # like optim, but with n_seg. divides the uppper-lower into n_seg intervals and starts to find the local optima
  # with respect to each endpoints.
  # The outputs is a list of minima (cost values is the last column) and the number of function evaluations.
  if(!is.function(fn))
    stop("'fn' must be a function.")
  fn1 <- function(arg)fn(arg, ...)

  u <- (upper - lower)/n_seg

  partition <- t(sapply(X = (0):(n_seg),
                        FUN = function(i) lower + i * u ))

  if(length(lower) == 1)
    partition <- t(partition)
  p0 <- do.call(`expand.grid`,as.data.frame(partition))
  result <- sapply(X = 1:dim(p0)[1],
                   FUN = function(i)optim(par = p0[i, ], fn = fn1,
                                          lower = lower, upper = upper,
                                          method = "L-BFGS-B",
                                          control = list(factr = control$factr)))

  minima <- result[1:2, 1:dim(p0)[1]]

  minima <- matrix(unlist(minima, use.names=FALSE), ncol = length(lower) + 1, byrow = TRUE)
  minima <- minima[!duplicated(round(minima[, -dim(minima)[2]], 3)), , drop = FALSE]
  counts <- sum(sapply(result[3,], "[[", 1))

  return(list(minima = minima, counts = counts))
}
######################################################################################################*
######################################################################################################*
construct_pen <- function(const, npred, k){
  dim_ind <- list()
  count1 <- 1
  count2 <- k
  for (j in 1:npred){
    dim_ind[[j]] <- count1:count2
    count1 <- count2 + 1
    count2 <- count2 + k
  }

  temp1 <- c()
  for (i in 1:nrow(const$ui)){
    temp2 <- c()
    for(l in 1:npred){
      temp2[l] <- paste(const$ui[i, l], " * x[c(", paste(dim_ind[[l]], collapse = ","), ")]", sep = "", collapse = " + ")
    }
    temp3 <- paste(temp2, collapse = " + ")
    temp1[i] <- paste("pmin(0, (", paste(const$ci[i], temp3, sep = " + "), "))^2", sep = "")
  }
  pen_char <- paste("sum(", paste(const$coef, "*", temp1, collapse = "+"), ")")
  pen <- NA ## to define the variable in the global environment and avoid  check Note
  pen_func <- paste("pen <- function(x)return(", pen_char, ")", sep = "")
  eval(parse(text = pen_func))
  return(pen)
}
######################################################################################################*
######################################################################################################*
is.formula <- function(x){
  inherits(x, "formula")
}
######################################################################################################*
######################################################################################################*
#' @importFrom methods formalArgs
check_common_args <- function(fimfunc, formula, predvars, parvars, family,
                              lx, ux, iter, k, paramvectorized, prior, x,
                              user_crtfunc, user_sensfunc){

  # adding to prevent errors @seongho on 06192020
  mu = NULL
  grad = NULL

  ### it checks the formula, k , lx, ux and iter and returns the fisher information matrix
  if (is.null(fimfunc) & missing(formula))
    stop("either 'fimfunc' or 'formula' must be given")
  if (!is.null(fimfunc) & !missing(formula))
    stop("one of 'fimfunc' and 'formula' must be given")


  if (is.null(fimfunc)){
    if (is.null(family))
      stop("Bug: 'family' must not be NULL in this function!")
    if (is.character(family))
      family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family))
      family <- family()
    if (is.null(family$family)) {
      print(family)
      stop("'family' not recognized")
    }
    if (!all(is.character(parvars)))
      stop("'paravasr' must be character")
    if (!all(is.character(predvars)))
      stop("'predvars' must be character")
    if (!is.formula(formula))
      stop("'formula' must be of formula type")

    ### desvar: from user matrix!
    ## checking the formula and predvars and parvars
    # remove spaces fom predvars
    predvars <- gsub(" ", "", predvars, fixed = TRUE)
    allvars <- all.vars(formula)
    parvars_list <- strsplit(parvars, "=")
    fixedpars <- sapply(1:length(parvars_list), FUN = function(j)parvars_list[[j]][2])
    fixedpars <- gsub(" ", "", fixedpars, fixed = TRUE)
    parvars <- sapply(parvars_list, '[[', 1)
    parvars <- gsub(" ", "", parvars, fixed = TRUE)
    num_unknown_param <- sum(is.na(fixedpars))

    ## cheking parvars
    if (!all(parvars %in% allvars)){
      parvars_notinformula <- paste(parvars[which(!parvars %in% allvars)], collapse = ", ")
      stop(parvars_notinformula, " is (are) not given in ", formula)
    }
    if (length(unique(parvars)) != length(parvars))
      stop(paste(parvars[duplicated(parvars)], collapse = ", "),  " is replicated in 'parvars'")

    ## cheking predvars
    if (!all(predvars %in% allvars)){
      predvars_notinformula <- paste(predvars[which(!predvars %in% allvars)], collapse = ", ")
      stop(predvars_notinformula, " is (are) not given in ", formula)
    }
    if (length(unique(predvars)) != length(predvars))
      stop(paste(predvars[duplicated(predvars)], collapse = ", "),  " is replicated in 'predvars'")

    ## checking if parvars and predvars have interaction
    if (any(parvars %in% predvars))
      stop(paste(parvars[which((predvars %in% parvars))], sep = ", "), " is (are) repeated in both 'parvars' and 'predvars'")
    mu <- create_mu(formula = formula, predvars = predvars, parvars = parvars,  paramvectorized = paramvectorized, fixedpars = fixedpars)
    grad <- create_grad(formula = formula, predvars = predvars, parvars = parvars,  paramvectorized = paramvectorized, fixedpars = fixedpars)


    mu_sens <- create_mu(formula = formula, predvars = predvars, parvars = parvars,  paramvectorized = FALSE, fixedpars = fixedpars)
    grad_sens <- create_grad(formula = formula, predvars = predvars, parvars = parvars, paramvectorized = FALSE, fixedpars = fixedpars)


    fimfunc_formula <- function(x, w, param){
      fim(x = x, w =w, param = param, grad = grad, mu = mu, family = family, paramvectorized = paramvectorized)
    }

    ## becasue we use it for sensitivity check we don to vectorizedit with respect to the parameters but the design
    ## we do not have a set of parameters here

    fimfunc_sens_formula <- function(x, w, param){
      fim(x = x, w =w, param = param, grad = grad_sens, mu = mu_sens, family = family, paramvectorized = FALSE)
    }
  }else{


    if (!is.function(fimfunc))
      stop(" \"fimfunc\" should be a \"function\"")
    if (!all(formalArgs(fimfunc) %in% c("x", "w", "param")))
      stop("arguments of 'fimfunc' must include 'x', 'w' and 'param'")
    fimfunc_formula <- NULL
    fimfunc_sens_formula <- NULL
    num_unknown_param <- NA


  }
  if (missing(lx))
    stop("\"lx\" is missing")
  if (missing(ux))
    stop("\"ux\" is missing")
  if (length(lx) != length(ux))
    stop("Length of \"lx\" is not equal to length of \"ux\"")
  if (!all(lx <= ux))
    stop("'lx' must be less than 'ux'")
  if (is.null(fimfunc))
    if (length(predvars) != length(lx))
      stop("length of 'lx' is not equal to the number of design parameters")
  if (missing(k))
    stop("\"k\" is missing")
  if (!is.numeric(k) || (k %% 1) != 0 || k <= 0)
    stop("\"k\" must be a positive integer number")
  if (missing(iter))
    stop("\"iter\" is missing")
  if (!is.numeric(iter) || (iter %% 1) != 0 || iter <= 0)
    stop("\"iter\" must be a positive integer number")
  if (!is.null(prior)){
    if (class(prior) != "cprior")
      stop("'prior' must be a list of class 'cprior'. See ?normal")
    if (!is.list(prior))
      stop("'prior' must be a list")
    if (is.null(prior$lower) || is.null(prior$upper) || is.null(prior$fn) || is.null(prior$npar))
      stop("please complete the 'prior' list. 'min', 'max', 'npar' and 'fn' are required in a list")
    if (length(prior$lower) != length(prior$upper))
      stop("length of min and max in 'prior' is not equal")
    if (!all(prior$lower < prior$upper))
      stop("min must be less than max in 'prior'")
    fixedpars <- NULL
  }
  if (!is.null(x)){
    if (is.null(k))
      if (length(x)/k != length(lx))
        stop("Length of 'lx' is not equal to 'length(x)/k'. Check 'k', 'x', 'lx' and 'ux' to match.")
  }

  ###############*
  # user ----
  ###############*
  if (!is.null(user_crtfunc)){
    if (!is.function(user_crtfunc))
      stop("'user_crtfunc' must be a 'function'")
    if (!all(c("x", "w", "fimfunc") %in% formalArgs(user_crtfunc)))
      stop("'user_crtfunc' must have arguments 'x', 'w' and 'fimfunc'")
    if (!is.null(fimfunc)){
      if (!any(formalArgs(user_crtfunc) == "param"))
        stop("'user_crtfunc' must have an argument named 'param'")
      user_crtfunc2 <- user_crtfunc
    }else{
      if (!all(c(parvars) %in% formalArgs(user_crtfunc)))
        stop("'crtfunc' must additionally have arguments, ", paste(parvars, collapse = ", "))
      user_crtfunc2 <- create_user_crtfunc(parvars = parvars, user_crtfunc = user_crtfunc, paramvectorized = paramvectorized)
      fimfunc_sens_formula <- fimfunc_formula <- create_user_fim(parvars = parvars, grad = grad, mu = mu, family = family, paramvectorized = paramvectorized)
    }
  }else
    user_crtfunc2 <- NULL


  if (!is.null(user_sensfunc)){

    if (is.null(user_crtfunc))
      warning("The sensitivity function for the criterion is already implemented. Forgot to specify your 'crtfunc'?")
    if (!is.function(user_sensfunc))
      stop("'user_sensfunc' must be a 'function'")
    if (!all(c("xi_x", "x", "w", "fimfunc") %in% formalArgs(user_sensfunc)))
      stop("'user_sensfunc' must have arguments 'xi_x', 'x', 'w' and 'fimfunc'")
    if (!is.null(fimfunc)){
      if (!any(formalArgs(user_sensfunc) == "param"))
        stop("'user_sensfunc' must have an argument named 'param'")
      user_sensfunc2 <- user_sensfunc
    }else{
      if (!all(c(parvars) %in% formalArgs(user_sensfunc)))
        stop("'sensfunc' must additionally have arguments, ", paste(parvars, collapse = ", "))
      user_sensfunc2 <- create_user_sensfunc(parvars = parvars, user_sensfunc = user_sensfunc, paramvectorized = paramvectorized)
    }
  }else
    user_sensfunc2 <- NULL

  return(list(fimfunc_formula = fimfunc_formula,
              fimfunc_sens_formula = fimfunc_sens_formula,
              num_unknown_param = num_unknown_param,
              user_crtfunc = user_crtfunc2, user_sensfunc = user_sensfunc2
              # SK@03052020
              ,mu = mu
              ,grad = grad
              ))
}

######################################################################################################*
######################################################################################################*
return_ld_ud <- function(sym, equal_weight, k, npred, lx, ux, is.only.w){
  # create the upper bound and lower bound for and the dimension of decision variables
  #x_length can be
  if (!is.only.w){
    if (sym) {
      ## if the number of design points be odd then the middle point of the design shouyld be the symmetric point!
      ## the number of design point can be one less then w for odd k
      if (k %% 2 == 0) {
        w_length <-  k / 2
        x_length <-  k / 2
      }else{
        x_length <-  floor(k / 2)
        w_length <- floor(k / 2) + 1
      }
    }else{
      x_length  <-  k * npred
      w_length <-  k
    }
    # so if sym == TRUE then we have two possiblity:
    # the number of design points is odd, then the x is one element less than w and in the missed point is sym_point
    # but if number of design be even then the length of x is equal to length of w
    ld <- c()
    ud <- c()
    for (i in 1:npred) {
      ld <- c(ld, rep(lx[i], x_length / npred))
      ud <- c(ud, rep(ux[i], x_length / npred))
    }
    #now we should add the wights!
    if (!equal_weight) {
      ld <- c(ld, rep(0, w_length))
      ud <- c(ud, rep(1, w_length))
    }
  }else{
    # we dont use lx, ux here!
    w_length <-  k
    ld <- rep(0, w_length)
    ud <- rep(1, w_length)
  }
  return(list(ld = ld, ud = ud))
}
######################################################################################################*
######################################################################################################*
check_initial <- function(initial, ld, ud){
  ## dealing with initial countries if set by user
  if (!is.null(initial)) {
    # convert the initial to matrix if it is a vector!
    if (!is.matrix(initial))
      initial <- t(as.matrix(initial))
    initial <- round(initial, 8) #because it may produce strange error when cheking the lower upper bound!!!
    # check if the length is true!
    if (!is.na(initial) && dim(initial)[2] != length(ld))
      stop("The number of columns of 'initial' does not match with length of countries and should be", length(ld))
    #  we use round to protect the function from strange behaviour
    initial_out <-  sapply(1:dim(initial)[1], function(j) any(round(initial[j,], 5) > round(ud, 5)) || any(round(initial[j,], 5) < round(ld, 5)))
    if (any(initial_out))
      stop("The initial vaule(s) in row ", paste(which(initial_out), collapse = ", "), " are (is) out of bound.")
  }
  return(initial)
}
######################################################################################################*
######################################################################################################*
print_xw_char <- function(x, w,  npred, is.only.w, equal_weight){
  # if (npred == 1)
  #   return(x)
  if (!is.only.w){
    k <- length(x)/npred
    if ( k != length(w))
      stop("'k' is not equal to the length of 'w'")
    x <- sprintf("%.5f", round(x,5))
    x <- format(x, width = max(nchar(x)))
    if (npred != 1){
      brackets_left  <- brackets_right <- rep("|", npred)
      brackets_left[1] <- "/"
      brackets_left[npred] <- "\\"
      brackets_right[1] <- "\\"
      brackets_right[npred] <- "/"
    }else
      brackets_left  <- brackets_right <- rep("", npred)

    x_mat <- matrix(x, nrow = npred, byrow = TRUE)
    x_char <- sapply(X = 1:ncol(x_mat), function(j)paste(brackets_left, x_mat[, j], brackets_right, sep = ""))
    # adding a row points
    point_char <- format(paste("Points", 1:k, sep = ""), width = max(nchar(x_char)), justify = "centre")
    # adding the weights
    if (!equal_weight){
      w <- sprintf("%.3f", round(w,3))
      w <- format(w, width = max(nchar(x_char)), justify = "centre")
      weight_char <- format(paste("Weights", 1:k, sep = ""), width = max(nchar(x_char)), justify = "centre")
      x_char <- rbind(point_char, x_char, weight_char, w)
    } else
      x_char <- rbind(point_char, x_char)

    x_char <- sapply(X = 1:nrow(x_char), function(i)paste(x_char[i, ], collapse = " "))
    outchar <- paste(x_char, collapse = "\n ")
  }else{
    w <- sprintf("%.3f", round(w,3))
    outchar <- paste("Weights:", paste(w, collapse = " "))
  }

  return(outchar)

}

######################################################################################################*
######################################################################################################*
# plot_sens <- function(lower, upper, chi, chi_deriv, x, x_deriv){
#
#   # chi design space as grid chi
#   # derivative values for chi
#   # x design points
#   # x_deriv derivative values at design points
#   plot(chi, chi_deriv, type = "l",
#        col = "blue",   xlab = "Design Interval",
#        ylab = "Sensitivity Function", xaxt = "n", main = "Sensitivity Plot")
#   abline(h = 0, v = c(x) ,col = "grey", lty = 3)
#
#   points(x = x,  y = x_deriv, col = "red" ,pch = 16, cex = 1)
#
#   axis(side = 1, at = c(lower, chi, upper, 0),
#        labels = round(c(lower, chi, upper, 0), 4))
# }

######################################################################################################*
######################################################################################################*
create_FIM_LLTM <- function(Q, normalization){
  wlambda <- function(x, w, param, q){
    # x is the design points that are the values for the ability parameters
    #qparam <- apply(q * param, 2, sum)
    qparam <- -sum(q * param)
    # argument of lambda for each param vector
    arg_lambda <- sapply(1:length(qparam), function(k)sum(w * exp(qparam[k] + x)/(1 + exp(qparam[k] + x))^2))
  }

  fim_LLTM <- function(x, w, param){
    Qmat <- Q
    nitems <- nrow(Qmat)
    if (!is.null(normalization))
      param <- c(normalization, param)
    lmat <- lapply(1:nitems, FUN = function(j)wlambda(q = Qmat[j, ],x = x, w = w, param = param) * Qmat[j, ] %*% t(Qmat[j, ]))
    lmat <- apply(simplify2array(lmat), c(1, 2), sum)
    return(lmat)
  }
  fimfunc <- fim_LLTM
  return(fimfunc = fim_LLTM)
}


#############################################################################################################*
#############################################################################################################*


create_user_crtfunc <- function(parvars, user_crtfunc, paramvectorized = FALSE){
  npar <- length(parvars)
  crtfunc_formula <- NA ## to define the variable in the global environment and avoid R CMD check Note
  crtfunc_char <- "crtfunc_formula <- function(x, w, fimfunc, param)\n{\n  "
  ### for parameters
  if (!paramvectorized){
    crtfunc_char <- paste(crtfunc_char,  parvars[1], " <- param[1]", " \n  ", sep = "")
    if (npar > 1) {
      for (j in 2:npar) {
        crtfunc_char <- paste(crtfunc_char, parvars[j], " <- param[", j, "]", " \n  ", sep = "")
      }
    }
  }else{
    crtfunc_char <- paste(crtfunc_char,  parvars[1], " <- param[, 1]", " \n  ", sep = "")
    if (npar > 1) {
      for (j in 2:npar) {
        crtfunc_char <- paste(crtfunc_char, parvars[j], " <- param[, ", j, "]", " \n  ", sep = "")
      }
    }
  }
  crtfunc_char_func <- paste("out <-  user_crtfunc(", paste(parvars, parvars, sep = " = ", collapse = ", "), ",x = x, w = w, fimfunc = fimfunc)", sep = "")
  crtfunc_char <- paste( crtfunc_char,  crtfunc_char_func, "\n  return(out)\n}",  sep = "")
  #cat(crtfunc_char)
  eval(parse(text =  crtfunc_char))
  return( crtfunc_formula)
}

#############################################################################################################*
#############################################################################################################*
create_user_sensfunc <- function(parvars, user_sensfunc, paramvectorized = FALSE){
  npar <- length(parvars)
  sensfunc_formula <- NA ## to define the variable in the global environment and avoid R CMD check Note
  sensfunc_char <- "sensfunc_formula <- function(xi_x, x, w, fimfunc, param)\n{\n  "
  ### for parameters
  if (!paramvectorized){
    sensfunc_char <- paste(sensfunc_char,  parvars[1], " <- param[1]", " \n  ", sep = "")
    if (npar > 1) {
      for (j in 2:npar) {
        sensfunc_char <- paste(sensfunc_char, parvars[j], " <- param[", j, "]", " \n  ", sep = "")
      }
    }
  }else{
    sensfunc_char <- paste(sensfunc_char,  parvars[1], " <- param[, 1]", " \n  ", sep = "")
    if (npar > 1) {
      for (j in 2:npar) {
        sensfunc_char <- paste(sensfunc_char, parvars[j], " <- param[, ", j, "]", " \n  ", sep = "")
      }
    }
  }
  sensfunc_char_func <- paste("out <-  user_sensfunc(", paste(parvars, parvars, sep = " = ", collapse = ", "), ",xi_x = xi_x, x = x, w = w, fimfunc = fimfunc)", sep = "")
  sensfunc_char <- paste( sensfunc_char,  sensfunc_char_func, "\n  return(out)\n}",  sep = "")
  eval(parse(text =  sensfunc_char))

  return(sensfunc_formula)
}
#############################################################################################################*
#############################################################################################################*
create_user_fim <- function(parvars, grad, mu, family, paramvectorized = FALSE){
  # I need the arguments grad, mu and family to be present
  # to be recognised in the eval function

  npar <- length(parvars)
  ## to define the variables in the global environment and avoid R CMD check Note
  fimfunc_formula_original <- NA
  if (npar >=2)
    pars_char <- paste(parvars, collapse = ", ") else
      pars_char <- parvars
  if (paramvectorized) # for bayesian designs every parameter will be in one row
    pars_char_vec <- paste0("param <- cbind(", paste(pars_char, collapse = ","), ")") else
      pars_char_vec <- paste0("param <- c(", paste(pars_char, collapse = ","), ")")
  fim_char <- paste0("fimfunc_formula_original <- function(",
                     pars_char,
                     ", x, w)\n{\n",
                     pars_char_vec, "\n",
                     "out <-  fim(x = x, w = w, param = param, grad = grad, mu = mu, family = family, paramvectorized = paramvectorized)\nreturn(out)\n}")
  eval(parse(text =  fim_char))
  return(fimfunc_formula_original)
}
#############################################################################################################*
#############################################################################################################*
# create_FIM_paramvectorizes <- function(parvars, FIM){
#   npar <- length(parvars)
#   FIM_formula_parallel <- NA ## to define the variable in the global environment and avoid R CMD check Note
#   FIM_char <- "FIM_formula_parallel <- function(x, w, param)\n{\n  "
#   ### for parameters
#   FIM_char <- paste(FIM_char,  parvars[1], " <- param[, 1]", " \n  ", sep = "")
#   if (npar > 1) {
#     for (j in 2:npar) {
#       FIM_char <- paste(FIM_char, parvars[j], " <- param[, ", j, "]", " \n  ", sep = "")
#     }
#
#   }
#   FIM_char_func <- paste("out <-  FIM(", paste(parvars, parvars, sep = " = ", collapse = ", "), ", x = x, w = w)", sep = "")
#   FIM_char <- paste(FIM_char,  FIM_char_func, "\n  return(out)\n}",  sep = "")
#   #cat(crtfunc_char)
#   eval(parse(text =  FIM_char))
#   return(FIM_formula_parallel)
# }
