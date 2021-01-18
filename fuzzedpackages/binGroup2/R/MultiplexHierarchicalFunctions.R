###############################################################################
# This file contains functions to calculate the expected testing expenditure
# of hierarchical testing, and the individual specific measures of accuracy
# Pooling Sensitivity, Pooling Specificity, Pooling Positive Predicitive Value,
# and Pooling Negative Predicitive Value.
###############################################################################

# Begin Exp_tests.R functions

# Find E(T) for hierarchical group testing with multiplex assays (K=2) in
# heterogeneous populations
# Author: Chris Bilder

###############################################################################
# Testing error calculations: 1 - P(G_{sj} = (0,0) | G~_{sj} = g~_{sj})
#   Order for diseases 1 and 2 in the conditioning part
#   1   2
#   0   0
#   0   1
#   1   0
#   1   1

test.error2 <- function(Se, Sp, s, k, kbar) {
  # The Product k = 1, ... , K is done because there are both the k and kbar 
  #   diseases in kronecker()
  1 - kronecker(c(Sp[k,s], (1 - Se[k,s])), c(Sp[kbar,s], (1 - Se[kbar,s])))
}



###############################################################################
# Functions help to find the correct indices. For K > 2, it may be easier to 
#   use matrix multiplication like with MRCV research.
# Could find K within functions here but more efficient to find only once in 
#   the function which calls these.

# Row number in joint.p is for an infection k positive and all other infections 
#   are negatives. This is the same function as row.for.k() in PSe_PSp.R - If 
#   changes are made to it here, changes may need to be made for Exp_tests.R
row.for.k <- function(k, K, ordering) {
  row.numbers <- 1:2^K   # Row numbers in ordering
  row.numbers[ordering[,k] == 1  & rowSums(ordering) == 1]
}



# Row in joint.p where disease kbar is 0's
#   This is the same function as all.zero.kbar() in PSe_PSp.R - If changes are 
#   made to it here changes may need to be made for Exp_tests.R
all.zero.kbar <- function(k, K, ordering) {
  row.numbers <- 1:2^K   # Row numbers in ordering
  row.numbers[ordering[,-k] == 0]
}



# Row in joint.p matrix disease k is 1
#   This is the same function as all.one.k() in PSe_PSp.R - If changes are made 
#   to it here changes may need to be made for Exp_tests.R
all.one.k <- function(k, K, ordering) {
  row.numbers <- 1:2^K   # Row numbers in ordering
  row.numbers[ordering[,k] == 1]
}



# Row in joint.p where disease k is 0's
all.zero.k <- function(k, K, ordering) {
  row.numbers <- 1:2^K   # Row numbers in ordering
  row.numbers[ordering[,k] == 0]
}



###############################################################################
# Transition probability calculations: P(G~_{Sj}^{(s'+1)} = ___ | G~_{Sj}^{(s')} = ___)
#   Much of this programming took advantage of code first used for PSp 
#   calculations in PSePSpAllStages() - k, kbar, and K are included in a 
#   general manner here for possible future generalizations and/or PSe and PSp 
#   new programming

trans.prob.calc <- function(s, group.numb, group.mem, joint.p, k = 1, K = 2, 
                            ordering) {

  # P(G~_{Sj}^{(s'+1)} = (0,0) | G~_{Sj}^{(s')} = (0,0))
  prob00to00 <- 1
  # P(G~_{Sj}^{(s'+1)} = (0,1) | G~_{Sj}^{(s')} = (0,0))
  prob00to01 <- 0
  prob00to10 <- 0
  prob00to11 <- 0

  index.all.zeros <- rowSums(ordering) == 0
  index.all.zero.k <- all.zero.k(k = k, K = K, ordering = ordering)
  index.all.zero.kbar <- all.zero.kbar(k = k, K = K, ordering = ordering)

  # All people in group at stage s for group j
  group.numb.prev.stage <- group.mem[s, group.numb == group.mem[s+1,] & !is.na(group.mem[s+1,])][1]
  index.stage.s <- group.mem[s,] == group.numb.prev.stage & !is.na(group.mem[s,])
  # All people in group at stage s+1 for group j
  index.stage.splus1 <- group.mem[s+1,] == group.numb & !is.na(group.mem[s+1,])

  # prod(p_b00, b in B_{Sj}^{(s')} where j is group containing individual a
  prod.pb00.s <- prod(as.matrix(joint.p[index.all.zeros, 
                                        index.stage.s]))
  # prod(p_b00, b in B_{Sj}^{(s'+1)})
  prod.pb00.splus1 <- prod(as.matrix(joint.p[index.all.zeros, 
                                             index.stage.splus1]))
  # prod(p_b0+, b in B_{Sj}^{(s')}
  prod.pb0plus.s <- prod(colSums(as.matrix(joint.p[index.all.zero.k, 
                                                   index.stage.s])))
  # prod(p_b0+, b in B_{Sj}^{(s'+1)}
  prod.pb0plus.splus1 <- prod(colSums(as.matrix(joint.p[index.all.zero.k, 
                                                        index.stage.splus1])))
  # prod(p_b+0, b in B_{Sj}^{(s')}
  prod.pbplus0.s <- prod(colSums(as.matrix(joint.p[index.all.zero.kbar, 
                                                   index.stage.s])))
  # prod(p_b+0, b in B_{Sj}^{(s'+1)}
  prod.pbplus0.splus1 <- prod(colSums(as.matrix(joint.p[index.all.zero.kbar, 
                                                        index.stage.splus1])))

  # Who was in group at stage s that did not make it into group at stage s + 1?
  bar.group.mem <- group.mem[s,] == group.numb.prev.stage & !(group.numb == group.mem[s+1,]) & !is.na(group.mem[s+1,])

  # prod(p_b0+, b in Bbar_{Sj}^{(s'+1)}
  prod.pb0plus.bar.splus1 <- prod(colSums(as.matrix(joint.p[index.all.zero.k, 
                                                            bar.group.mem])))
  # prod(p_b+0, b in Bbar_{Sj}^{(s'+1)}
  prod.pbplus0.bar.splus1 <- prod(colSums(as.matrix(joint.p[index.all.zero.kbar, 
                                                            bar.group.mem])))


  # P(G~_{Sj}^{(s'+1)} = (0,0) | G~_{Sj}^{(s')} = (0,1))
  prob01to00 <- (prod.pb00.splus1 * prod.pb0plus.bar.splus1 - prod.pb00.s) / (prod.pb0plus.s - prod.pb00.s)
  # P(G~_{Sj}^{(s'+1)} = (0,1) | G~_{Sj}^{(s')} = (0,1))
  prob01to01 <- (prod.pb0plus.s - prod.pb00.splus1 * prod.pb0plus.bar.splus1) / (prod.pb0plus.s - prod.pb00.s)
  # prob01to00 + prob01to01 = 1 as it should, could use 1 - prob01to00 then
  prob01to10 <- 0
  prob01to11 <- 0

  prob10to00 <- (prod.pb00.splus1 * prod.pbplus0.bar.splus1 - prod.pb00.s) / (prod.pbplus0.s - prod.pb00.s)
  prob10to01 <- 0
  prob10to10 <- (prod.pbplus0.s - prod.pb00.splus1 * prod.pbplus0.bar.splus1) / (prod.pbplus0.s - prod.pb00.s)
  prob10to11 <- 0

  denominator <- 1 - prod.pbplus0.s - prod.pb0plus.s + prod.pb00.s
  prob11to00 <- (prod.pb00.splus1 - prod.pb00.splus1 * prod.pbplus0.bar.splus1 -
                   prod.pb00.splus1 * prod.pb0plus.bar.splus1 +  prod.pb00.s) / denominator
  prob11to01 <- (prod.pb0plus.splus1 - prod.pb00.splus1 - prod.pb0plus.s +
                   prod.pb00.splus1 * prod.pb0plus.bar.splus1) / denominator
  prob11to10 <- (prod.pbplus0.splus1 - prod.pb00.splus1 - prod.pbplus0.s +
                   prod.pb00.splus1 * prod.pbplus0.bar.splus1 ) / denominator
  prob11to11 <- (1 - prod.pbplus0.splus1 - prod.pb0plus.splus1 + prod.pb00.splus1) / denominator
  # prob11to00 + prob11to01 + prob11to10 + prob11to11   # = 1 as it should

  c(prob00to00, prob00to01, prob00to10, prob00to11,
    prob01to00, prob01to01, prob01to10, prob01to11,
    prob10to00, prob10to01, prob10to10, prob10to11,
    prob11to00, prob11to01, prob11to10, prob11to11)
}



###############################################################################
# Main function used to calculate E(T)
#   NOTE: Lists are used a lot in this function to store information so that I 
#         can dynamically add additional information depending on the stage 
#         number and group membership matrix. A matrix storage method would not 
#         be as flexible. Also, lists within lists are used for a 2D list - 
#         again, this is very flexible.
#   k, kbar, and K are included for possible future generalizations and/or PSe 
#         and PSp new programming

# Brianna Hitt - 03-06-19
# Function edited - removed the K argument and define K in the first lines 
#   of the function

ET.all.stages.new <- function(joint.p, group.mem, Se, Sp, ordering) {

  k <- 1    # k represents the disease of interest

  # Number of diseases
  K <- 2  # K represents the total number of diseases
  infections <- 1:K
  kbar <- infections[-k]   # Infections other than k
  # the above values (k, kbar, K) are set here in anticipation of
  #   future generalizations of this function; until the
  #   generalization is made, the values above should not be changed
  # when the generalization is made, substitute K <- ncol(ordering)
  #   for K <- 2
  # until the generalization is made, kbar <- 2

  # Number of stages
  S <- nrow(group.mem)

  # Number of groups to be tested at each stage
  c <- numeric(S)
  c[1] <- 1
  # Number of subgroups to be formed when a positive group occurs for stage s 
  #   group j
  m.sj <- list()
  # Size of group j in stage s
  I.sj <- list()
  I.sj[[1]] <- ncol(group.mem)

  # Find c_s, I_sj, m_sj
  for (s in 2:S) {

    c[s] <- length(unique(group.mem[s,])) - anyNA(group.mem[s,])
    # Find I.sj - probably could combine with the next loop, one problem is 
    #   getting group sizes for last stage due to 1:c[s-1] index
    group.size <- numeric(c[s])
    for (group.numb in 1:c[s]) {
      # Since NA's get included in the vector, I needed to include the & part here
      group.size[group.numb] <- length(group.mem[s, group.mem[s,] == group.numb & !is.na(group.mem[s,])])
    }
    I.sj[[s]] <- group.size

    # Find m.sj - NOTE: m.sj = 0 is put into the list - if do not want this, 
    #   use an if (I.sj[[s]][i] != 1) {} like code
    save.m <- numeric(c[s-1])
    for (group.numb in 1:(c[s-1])) {
      save.m[group.numb] <- length(unique(group.mem[s, group.mem[s-1,] == group.numb & !is.na(group.mem[s,])] ))
    }
    m.sj[[s-1]] <- save.m
  }


  #################
  # P(G_11+ > 0)

  prob <- list()
  prob[[1]]<-list()

  # P(G~_11 = (0,0))
  index.all.zeros <- rowSums(ordering) == 0
  theta00 <- prod(as.matrix(joint.p[index.all.zeros,]))

  # P(G~_11 = (0,1))
  index.all.zero.k <- all.zero.k(k = k, K = K, ordering = ordering)
  theta01 <- prod(colSums(as.matrix(joint.p[index.all.zero.k,]))) - theta00

  # P(G~_11 = (1,0))
  index.all.zero.kbar <- all.zero.kbar(k = k, K = K, ordering = ordering)
  theta10 <- prod(colSums(as.matrix(joint.p[index.all.zero.kbar,]))) - theta00

  # P(G~_11 = (1,1))
  theta11 <- 1 - theta01 - theta10 - theta00

  theta <- c(theta00, theta01, theta10, theta11)
  # If one were to sum (testing error part) * P(G~_11 = g~_11), the result 
  #   would be P(G_11+ > 0)
  prob[[1]][[1]] <- test.error2(Se = Se, Sp = Sp, s = 1, k = k, kbar = kbar) * theta


  #################
  # P(G_11+ > 0, ..., G_sj+ > 0) = P(G_{sj+}^(1) > 0, ..., G_{sj+}^(s) > 0)

  joint.prob <- list()
  joint.prob[[1]] <- list()
  joint.prob[[1]][[1]] <- prob[[1]][[1]]

  ET.sj <- list()
  ET.sj[[1]] <- list()
  ET.sj[[1]][[1]] <- m.sj[[1]][[1]] * sum(joint.prob[[1]][[1]])

  if (S > 2) {

    for (s in 2:(S-1)) {

      prob[[s]] <- list()
      joint.prob[[s]] <- list()
      ET.sj[[s]] <- list()

      for (j in 1:c[s]) {

        if (I.sj[[s]][j] != 1) {

          # Testing error and conditional probability part for stage s'+1
          test.error.part <- rep(test.error2(Se = Se, Sp = Sp, s = s, k = k, 
                                             kbar = kbar), times = 4)
          prob[[s]][[j]] <- test.error.part * trans.prob.calc(s = s - 1, 
                                                              group.numb = j, group.mem = group.mem,
                                                              joint.p = joint.p, K = K, ordering = ordering)

          group.numb.prev.stage <- group.mem[s-1, j == group.mem[s,] & !is.na(group.mem[s,])][1]

          if (s == 2) {
            # Testing error and conditional probability part for stage 1 to 2 
            #   multiplied by the testing error and stage 1 part. The 
            #   elementwise multiplication is properly aligned as
            #  g111 g112    g2j1 g2j2
            #     0    0       0    0
            #     0    0       0    1
            #     0    0       1    0
            #     0    0       1    1
            #     0    1       0    0
            #     0    1       0    1
            # ...
            joint.prob[[s]][[j]] <- rep(joint.prob[[s-1]][[group.numb.prev.stage]], 
                                        each = 4) * prob[[s]][[j]]
            # Reason for rep( , each = 4) when s = 2: P(G_11 = g_11) has only 
            #   four values and the transitions have 16 - writing out E(T) for 
            #   S = 3 helps to see the need for this special case
          } else
          {
            # Components of P(G_{sj+}^(1) > 0, ..., G_{sj+}^(s) > 0, G~_{sj+}^(1) > 0, ..., G~_{sj+}^(s) > 0)
            joint.prob[[s]][[j]] <- rep(joint.prob[[s-1]][[group.numb.prev.stage]], 
                                        each = 4) * rep(prob[[s]][[j]], times = 4)   
            # Reason for the rep() showing an example:
            #  joint.prob[[s-1]]      prob[[s]][[j]]
            #  G11	  G21 			      G21	G31
            #   1	2   1 2	            1 2 1 2
            #   0	0		0	0		          0	0	0	0
            #   						          0	0	0	1
            #   					          	0	0	1	0
            #   						          0	0	1	1
            #   0	0		0	1		          0	1	0	0
            #   						          0	1	0	1
            #   						          0	1	1	0
            #   						          0	1	1	1
            #   0	0		1	0		          1	0	0	0
            #   						          1	0	0	1
            #   						          1	0	1	0
            #   						          1	0	1	1
            #   0	0		1	1		          1	1	0	0
            #   						          1	1	0	1
            #   						          1	1	1	0
            #   						          1	1	1	1
            #   0	1		0	0		          0	0	0	0
            #   						          0	0	0	1
          }

          # print(m.sj[[s]][[j]] * sum(joint.prob[[s]][[j]]))
          #E(T_sj)
          ET.sj[[s]][[j]] <- m.sj[[s]][[j]] * sum(joint.prob[[s]][[j]])  

        }
      }
    }
  }


  ExpT <- 1 + sum(unlist(ET.sj))

  list(Expt = ExpT, c = c, I.sj = I.sj, m.sj = m.sj)
}



# Begin PSe_PSp.R functions

###############################################################################
# Calculate PSe and PSp for any number of stages and K=2
#
# The calculations for PSe and PSp are performed by breaking up the large
#   summation for each of them into parts. These parts then are simply put into
#   long vectors and elementwise multiplication is performed. Elements in
#   the resulting vector are then summed to produce the desired results.
#
# Notes on understanding the code
# 1) The individual "i" index of the paper is denoted by an "a" here instead.
# 2) The object S.all is "S" in the paper. The object S is "L" in the paper.
#
# Author: Chris Bilder
###############################################################################
# Additional functions needed for PSeAllStages()

# Used in testing error part
#   Calculates 1 - P(G_{Lj+}^{(s')} = 0 | G~_{Ljkbar}^{(s')} = g~_{Ljkbar}^{(s')}, Y~_ak = 1)
#   g~_{Ljkbar}^{(s') = 0 part comes before 1 part
test.error <- function(Se, Sp, s, k, kbar) {
  1 - (1 - Se[k,s]) * c(Sp[kbar,s], (1 - Se[kbar,s]))
}



# Used in transition probability part
# For k = 1, prod(p_b+0) over b <> a in stage s for group of interest

# Brianna Hitt - 03-06-19
# Function edited - removed the K argument and define K and infections in the 
#   first lines of the function

lambda.calc <- function(joint.p, group.mem, s, a, k, ordering) {

  # number of diseases - this is here for future generalizations
  K <- ncol(ordering)
  infections <- 1:K

  # I.11 is I because this function only called for S > 2
  I.11 <- ncol(group.mem)   

  # Included "& !is.na(group.mem[s,])" to deal with NAs for algorithms 
  #   where it is not possible for everyone to be tested in all stages
  # All people in group except individual a
  index.p.part.no.a <- group.mem[s,] == group.mem[s,a] & !is.na(group.mem[s,]) & (1:I.11 != a)  
  index.pb.plus.AllZeros <- all.zero.kbar(k = k, K = K, ordering = ordering)

  # prod.pb.plus.AllZeros.SubGroup
  #  Added as.matrix() for the case when subgroup size was only 2 leading to
  #    joint.p being only one individual (a vector) and colSums() will not work for it
  prod(colSums(as.matrix(joint.p[index.pb.plus.AllZeros, index.p.part.no.a])))
}



###############################################################################
# Called by PSeAllStages() to calculate PSe
#   Calculate PSe for K = 2 infections; There are some coding conventions used 
#   here to help generalize to K > 2 in future

# Brianna Hitt - 03-06-19
# Function edited - removed the K argument and define K and infections in the 
#   first lines of the function

PSeAllStages <- function(joint.p, group.mem, a, k, Se, Sp, ordering) {

  # number of diseases - this is here for future generalizations
  K <- ncol(ordering)
  infections <- 1:K

  # Infections other than k
  kbar <- infections[-k]  
  # Number of total stages for algorithm
  S.all <- nrow(group.mem)    
  # Number of stages for individual a where assume group.mem is coded 
  #   correctly - NAs only appear at end
  S <- sum(!is.na(group.mem[,a])) 


  #############################################################################
  #Testing error calculations excluding the leading S_e:Ljk

  # Stage 1
  test.error.part <- test.error(Se = Se, Sp = Sp, s = 1, k = k, kbar = kbar)  
  if (S > 2) {
    for(s in 2:(S-1)) {
      test.error.part.newstage <- test.error(Se = Se, Sp = Sp, s = s, k = k, 
                                             kbar = kbar)
      test.error.part <- kronecker(test.error.part, test.error.part.newstage)
    }
  }
  # Order for four stages: (stage 1, stage2, stage3) = (0,0,0), (0,0,1), 
  #                                                    (0,1,0), (0,1,1),
  #                                                    (1,0,0), (1,0,1), 
  #                                                    (1,1,0), (1,1,1)
  # Could also use expand.grid() and take prod() of rows


  ########################################################
  # Calculations needed for P(G~_{LjKBAR}^{(1)} = g~_{LjKBAR}^{(1)} | Y~_ak = 1)
  #   = P(G~_11k = g~_11k | Y~_ak = 1)
  # Note by including the ordering object here, this is how I can look at 
  #   disease 2 when needed. Thus, disease 2 essentially becomes disease 1.

  index.pa.1andAllZeros <- row.for.k(k = k, K = K, ordering = ordering)
  # For k = 1, p_a10
  pa.1andAllZeros <- joint.p[index.pa.1andAllZeros, a]
  # The "AllZeros" terminology is in anticipation of generalizing for K > 2. 
  #   For now, it is just one 0.

  index.pa.1andPlus <- all.one.k(k = k, K = K, ordering = ordering)
  # For k = 1, p_a1+
  pa.1andPlus <- sum(joint.p[index.pa.1andPlus, a])

  index.pb.plus.AllZeros <- all.zero.kbar(k = k, K = K, ordering = ordering)
  # For k = 1, prod(p_b+0) over b <> a
  prod.pb.plus.AllZeros <- prod(colSums(joint.p[index.pb.plus.AllZeros,-a]))

  # P(G~_11kbar = 0 | Y~_ak = 1)
  end.prob <- pa.1andAllZeros * prod.pb.plus.AllZeros / pa.1andPlus
  end.prob.vec <- rep(x = c(end.prob, (1 - end.prob)), each = 2^(S-2))
  # Order for four stages: (stage 1, stage2, stage3) = (0,0,0), (0,0,1), 
  #                                                    (0,1,0), (0,1,1),
  #                                                    (1,0,0), (1,0,1), 
  #                                                    (1,1,0), (1,1,1)


  ########################################################
  # Transition probability calculations
  #   P(G~_{Ljkbar}^{(s'+1)} = g~_{LjKBAR}^{(s'+1)} | G~_{Ljkbar}^{(s')} = g~_{LjKBAR}^{(s')}, Y~_ak = 1)

  # Need this here for S = 2 so that calculations work when execute
  nu <- 1  
  #   save.stage <- test.error.part*nu*end.prob.vec
  # Use a list storage format so that it can be dynamically extended to a new 
  #   size depending on S
  nu.trans.prob <- list()  
  # For diagnostic purposes only, not printed or saved outside of function
  # nu.trans.label <- list()  
  if (S > 2) {
    for (s in 1:(S-2)) {

      # For k = 1, prod(p_b+0) over b <> a in stage s' + 1 for group of interest
      lambda.minus.low <- lambda.calc(joint.p = joint.p, group.mem = group.mem, 
                                      s = s + 1, a = a,
                                      k = k, ordering = ordering)

      # For k = 1, prod(p_b+0) over b <> a in stage s' for group of interest
      lambda.minus.high <- lambda.calc(joint.p = joint.p, group.mem = group.mem, 
                                       s = s, a = a,
                                       k = k, ordering = ordering)

      # nu notation: nu.(stage s').(stage s'+1)

      # P(G~_{Ljkbar}^{(s'+1)} = 0 | G~_{Ljkbar}^{(s')} = 0, Y~_ak = 1)
      nu.0.0 <- 1

      # P(G~_{Ljkbar}^{(s'+1)} = 0 | G~_{Ljkbar}^{(s')} = 1, Y~_ak = 1)
      nu.1.0 <- pa.1andAllZeros * (lambda.minus.low - lambda.minus.high) /
        (pa.1andPlus - pa.1andAllZeros * lambda.minus.high)

      # P(G~_{Ljkbar}^{(s'+1)} = 1 | G~_{Ljkbar}^{(s')} = 0, Y~_ak = 1)
      nu.0.1 <- 0

      # P(G~_{Ljkbar}^{(s'+1)} = 1 | G~_{Ljkbar}^{(s')} = 1, Y~_ak = 1)
      nu.1.1 <- (pa.1andPlus - pa.1andAllZeros * lambda.minus.low) /
        (pa.1andPlus - pa.1andAllZeros * lambda.minus.high)

      nu.trans.prob[[s]] <- c(nu.0.0, nu.0.1, nu.1.0, nu.1.1)
      # Diagnostic purposes
      # nu.trans.label[[s]] <- paste("Transition from", s , "to", s+1)  
    }

    nu <- nu.trans.prob[[1]]
    if (S > 3) {
      for (s in 2:(S-2)) {
        nu <- rep(x = nu.trans.prob[[s]], times = 2) * rep(x = nu, each = 2)
      }
    }
  }
  # Example of nu ordering for S = 4
  #                  g1  g2  g3  nu.g2->g3  nu.g1->g2
  #                   0   0   0     0->0       0->0
  #                   0   0   1     0->1       0->0
  #                   0   1   0     1->0       0->1
  #                   0   1   1     1->1       0->1
  #                   1   0   0     0->0       1->0
  #                   1   0   1     0->1       1->0
  #                   1   1   0     1->0       1->1
  #                   1   1   1     1->1       1->1
  # In the PSe expression, the first summation is over g1, the second summation 
  #   is over g2, ..., so this controls the ordering of terms in the nu vector

  save.stage <- test.error.part * nu * end.prob.vec
  # All vectors are in the correct order for elementise multiplication

  Se[k,S] * sum(save.stage)
}

###############################################################################
# Calculate PSp and call PSeAllStages() to calculate PSe; all calculations are 
#   for K = 2
# The computational approach here is to find all possible terms in the large 
#   summations for PSp and PSe and then simply sum these terms up. For example, 
#   consider the case of finding P(Y_ik = 1) and S = 3. We need to sum over 
#   g~_{3j}^{(1)}, g~_{3j}^{(2)}, and g~_3jk. There are a total of 4 * 4 * 2 = 32 
#   separate terms in the summation. At intermediate steps, these terms are put
#   in the order of
#
#   g~_{3jk}^(1) g~_{3jkbar}^(1)    g~_{3jk}^(2) g~_{3jkbar}^(2)    g~_3jk
#              0               0               0               0         0
#              0               0               0               0         1
#              0               0               0               1         0
#              0               0               0               1         1
#              0               0               1               0         0
#              0               0               1               0         1
#              0               0               1               1         0
#              0               0               1               1         1
#
#              0               1               0               0         0
#              0               1               0               0         1
#
#...
#
#              1               1               1               1         0
#              1               1               1               1         1
#
# Once all of the terms are found, they are summed to find P(Y_ik = 1).


# Calculate PSe and PSp for any number of stages and K=2
#
# The calculations for PSe and PSp are performed by breaking up the large
#   summation for each of them into parts. These parts then are simply put into
#   long vectors and elementwise multiplication is performed. Elements in
#   the resulting vector are then summed to produce the desired results.
#
# Notes on understanding the code
# 1) The individual "i" index of the paper is denoted by an "a" here instead.
# 2) The object S.all is "S" in the paper. The object S is "L" in the paper.

# Brianna Hitt - 03-06-19
# Function edited - removed the K argument and define K and infections in the 
#   first lines of the function

PSePSpAllStages <- function(joint.p, group.mem, a, k, Se, Sp, ordering) {

  # number of diseases - this is here for future generalizations
  K <- ncol(ordering)
  infections <- 1:K

  kbar <- infections[-k]   # Infections other than k
  S.all <- nrow(group.mem) # Number of total stages for algorithm
  # Number of stages for individual a where assume group.mem is coded 
  #   correctly - NAs only appear at end
  S <- sum(!is.na(group.mem[,a])) 
  I.11 <- ncol(group.mem)  # I=I.11 for S > 2 and I for S = 2



  ################################
  # Testing error calculations for P(Y_ik = 1). This will compute all needed 1 
  #   - Prod((1 - Se)^g~ * Sp^(1-g~)) values and put into the correct order
  # Stage 1
  test.error.part <- test.error2(Se = Se, Sp = Sp, s = 1, k = k, kbar = kbar) 
  # Stages 2 to S - 1
  if (S > 2) {
    for(s in 2:(S-1)) {
      test.error.part.newstage <- test.error2(Se = Se, Sp = Sp, s = s, k = k, 
                                              kbar = kbar)
      test.error.part <- kronecker(test.error.part, test.error.part.newstage)
    }
  }
  # Last stage
  test.error.part <- kronecker(test.error.part, c((1-Sp[k,S]), Se[k,S]))
  #Example ordering when S = 2 for G11 and g2ak
  #  g1k g1kbar g2ak
  #   0    0     0
  #   0    0     1
  #   0    1     0
  #   0    1     1
  #   1    0     0
  #   1    0     1
  #   1    1     0
  #   1    1     1



  ###################################
  # theta_g~_11 = P(G~_{Lj}^{(1)} = g~_{Lj}^{(1)}) = P(G~_11 = g~_11) -
  # Note by including the ordering object here, this is how I can look at 
  #   disease 2 when needed. Thus, disease 2 essentially becomes disease 1.
  # k is the first infection and kbar is the second infection given in my 
  #   code labels below
  index.all.zeros <- rowSums(ordering) == 0
  theta00 <- prod(as.matrix(joint.p[index.all.zeros,]))

  index.all.zero.k <- all.zero.k(k = k, K = K, ordering = ordering)
  theta01 <- prod(colSums(as.matrix(joint.p[index.all.zero.k,]))) - theta00

  index.all.zero.kbar <- all.zero.kbar(k = k, K = K, ordering = ordering)
  theta10 <- prod(colSums(as.matrix(joint.p[index.all.zero.kbar,]))) - theta00

  theta11 <- 1 - theta01 - theta10 - theta00

  theta <- rep(x = c(theta00, theta01, theta10, theta11), each = 2 * 4^(S-2))
  # Reason for "each = 2 * 4^(S-2)": These four items need to be repeated each 
  #   2 * 4^(S-2) times. In the case of S = 3, this means each is repeated 8 
  #   times. The comment before the main function shows that each set of g1k 
  #   and g1kbar is repeated 8 times.



  #############################################################################
  # Transition from stage S-1 to stage S, i.e., P(G~_Sjk = g~_Sjk | G~_{Sj}^{(S-1)} = g~_{Sj}^{(S-1)})
  #   Example notation for S = 2: prob01to0 means for g1k = 0, g1kbar = 1, g2ak = 0
  #   i.e., prob k kbar to k
  prob00to0 <- 1
  prob00to1 <- 0

  # All people in group containing individual a except individual a
  index.p.part.no.a <- group.mem[S-1,] == group.mem[S-1,a] & !is.na(group.mem[S-1,]) & (1:I.11 != a)
  # All people in group containing individual a
  index.p.part.a <- group.mem[S-1,] == group.mem[S-1,a] & !is.na(group.mem[S-1,])

  #If k = 1, this would be Prod(p_b+0, b <> a).
  prod.p.bplus0.no.a <- prod(colSums(as.matrix(joint.p[index.all.zero.kbar, 
                                                       index.p.part.no.a])))

  #If k = 1, Prod(p_b+0)
  prod.p.bplus0 <- prod(colSums(as.matrix(joint.p[index.all.zero.kbar, 
                                                  index.p.part.a])))

  # Prod(p_b00)
  prod.pb00 <- prod(as.matrix(joint.p[index.all.zeros, index.p.part.a]))

  #p_a00
  pa00 <- joint.p[index.all.zeros, a]

  # For k = 1, p_a0+
  pa.0andPlus <- sum(joint.p[index.all.zero.k, a])

  #If k = 1, Prod(p_b0+)
  prod.p.b0plus <- prod(colSums(as.matrix(joint.p[index.all.zero.k, 
                                                  index.p.part.a])))

  # P(G~_Sjk = 0 | G~_{Sj}^{(S-1)} = (0,1)) with k = 1, kbar = 2
  #   If want k = 2, the code takes this into account with how the ordering 
  #   object is used to re-index what is considered k and kbar 
  #   (i.e., "2" becomes "1")
  prob01to0 <- 1

  # P(G~_Sjk = 1 | G~_{Sj}^{(S-1)} = (0,1)) with k = 1, kbar = 2
  prob01to1 <- 0

  # P(G~_Sjk = 0 | G~_{Sj}^{(S-1)} = (1,0)) with k = 1, kbar = 2
  # k=1, S=2: P(G~_2a1 = 0| G~_11=(1,0))
  prob10to0 <- (pa00 * prod.p.bplus0.no.a - prod.pb00)/(prod.p.bplus0 - prod.pb00)

  # P(G~_Sjk = 1 | G~_{Sj}^{(S-1)} = (1,0)) with k = 1, kbar = 2
  # Could do compliment of above
  prob10to1 <- (prod.p.bplus0 - pa00 * prod.p.bplus0.no.a)/(prod.p.bplus0 - prod.pb00)

  # P(G~_Sjk = 0 | G~_{Sj}^{(S-1)} = (1,1)) with k = 1, kbar = 2
  prob11to0 <- (pa.0andPlus - pa00 * prod.p.bplus0.no.a - prod.p.b0plus + prod.pb00) /
    (1 - prod.p.bplus0 - prod.p.b0plus + prod.pb00)

  # P(G~_Sjk = 1 | G~_{Sj}^{(S-1)} = (1,1)) with k = 1, kbar = 2
  prob11to1 <- (1 - pa.0andPlus - prod.p.bplus0 + pa00 * prod.p.bplus0.no.a) /
    (1 - prod.p.bplus0 - prod.p.b0plus + prod.pb00)

  #Example ordering for G11 and g2ak
  #  g1k g1kbar g2ak
  #   0    0     0
  #   0    0     1
  #   0    1     0
  #   0    1     1
  #   1    0     0
  #   1    0     1
  #   1    1     0
  #   1    1     1
  prob.Sminus1.k.kbar.to.Sk <- rep(x = c(prob00to0, prob00to1, prob01to0, 
                                         prob01to1, prob10to0, prob10to1, 
                                         prob11to0, prob11to1), times = 4^(S-2))
  # "times = 4^(S-2)" corresponds to: 1) There are 2^K = 2^2 = 4 possible 
  #    values for each G~_sj vector. Because these probabilties correspond 
  #    to stages S and S-1, the probabilities need to be repeated over the 
  #    first S - 2 stages. Thus, these probabilities correspond to the last 
  #    3 columns of the example ordering given in the comments for this function.
  #
  #Example ordering for Gstage1, Gstage2, and g3ak
  #  g1k g1kbar g2k g2kbar g3ak
  #   0    0     0     0    0
  #   0    0     0     0    1
  #   0    0     0     1    0
  #   0    0     0     1    1
  #   0    0     1     0    0
  #   0    0     1     0    1
  #   0    0     1     1    0
  #   0    0     1     1    1

  #   0    1     0     0    0
  #   0    1     0     0    1
  # ...



  ###################################

  #Transition from stage s-1 to stage s
  #  P(G~_{Sj}^{(s'+1)} = g~_{Sj}^{(s'+1)} | G~_{Sj}^{(s')} = g~_{Sj}^{(s')})
  # For example, with S = 3, we need to find P(G~_{3j}^(2) = g~_{3j}^(2) | G~_11 = g_11)
  # Notation for S = 3: prob00to00 means for g1k = 0, g1kbar = 0 to g2ak = 0, g2akbar = 0
  #    prob k kbar at stage s' to k kbar at stage s'+1

  trans.prob <- 1
  trans.prob.list <- list()
  if (S > 2) {
    for (s in 1:(S-2)) {

      # For K > 2, it may be better to put all probabilties in a matrix or 
      #   K-dimensional array

      # P(G~_{Sj}^{(s'+1)} = (0,0) | G~_{Sj}^{(s')} = (0,0))
      prob00to00 <- 1
      # P(G~_{Sj}^{(s'+1)} = (0,1) | G~_{Sj}^{(s')} = (0,0))
      prob00to01 <- 0
      prob00to10 <- 0
      prob00to11 <- 0

      # All people in group at stage s for individual a
      index.stage.s <- group.mem[s,] == group.mem[s,a] & !is.na(group.mem[s,])
      # All people in group at stage s+1 for individual a
      index.stage.splus1 <- group.mem[s+1,] == group.mem[s+1,a] & !is.na(group.mem[s+1,])

      # S=3, prod(p_b00, b in B_{Sj}^{(s')} where j is group containing individual a
      prod.pb00.s <- prod(as.matrix(joint.p[index.all.zeros, index.stage.s]))
      # S=3, prod(p_b00, b in B_{Sj}^{(s'+1)})
      prod.pb00.splus1 <- prod(as.matrix(joint.p[index.all.zeros, 
                                                 index.stage.splus1]))
      # S=3, prod(p_b0+, b in B_{Sj}^{(s')}
      prod.pb0plus.s <- prod(colSums(as.matrix(joint.p[index.all.zero.k, 
                                                       index.stage.s])))
      # S=3, prod(p_b0+, b in B_{Sj}^{(s'+1)}
      prod.pb0plus.splus1 <- prod(colSums(as.matrix(joint.p[index.all.zero.k, 
                                                            index.stage.splus1])))
      # S=3, prod(p_b+0, b in B_{Sj}^{(s')}
      prod.pbplus0.s <- prod(colSums(as.matrix(joint.p[index.all.zero.kbar, 
                                                       index.stage.s])))
      # S=3, prod(p_b+0, b in B_{Sj}^{(s'+1)}
      prod.pbplus0.splus1 <- prod(colSums(as.matrix(joint.p[index.all.zero.kbar, 
                                                            index.stage.splus1])))

      #Who was in group at stage s that did not make it into stage s + 1 group?
      group.numb <- group.mem[s+1,a]  # Group number of individual a stage s + 1
      group.numb.prev.stage <- group.mem[s,a]
      bar.group.mem <- group.mem[s,] == group.numb.prev.stage & !(group.numb == group.mem[s+1,]) & !is.na(group.mem[s+1,])
      # prod(p_b0+, b in Bbar_{Sj}^{(s'+1)}
      prod.pb0plus.bar.splus1 <- prod(colSums(as.matrix(joint.p[index.all.zero.k, 
                                                                bar.group.mem])))
      # prod(p_b+0, b in Bbar_{Sj}^{(s'+1)}
      prod.pbplus0.bar.splus1 <- prod(colSums(as.matrix(joint.p[index.all.zero.kbar, 
                                                                bar.group.mem])))

      # P(G~_{Sj}^{(s'+1)} = (0,0) | G~_{Sj}^{(s')} = (0,1))
      prob01to00 <- (prod.pb00.splus1 * prod.pb0plus.bar.splus1 - prod.pb00.s) / (prod.pb0plus.s - prod.pb00.s)
      # P(G~_{Sj}^{(s'+1)} = (0,1) | G~_{Sj}^{(s')} = (0,1))
      prob01to01 <- (prod.pb0plus.s - prod.pb00.splus1 * prod.pb0plus.bar.splus1) / (prod.pb0plus.s - prod.pb00.s)
      # prob01to00 + prob01to01 = 1 as it should - could use 1 - prob01to00 then
      prob01to10 <- 0
      prob01to11 <- 0

      prob10to00 <- (prod.pb00.splus1 * prod.pbplus0.bar.splus1 - prod.pb00.s) / (prod.pbplus0.s - prod.pb00.s)
      prob10to01 <- 0
      prob10to10 <- (prod.pbplus0.s - prod.pb00.splus1 * prod.pbplus0.bar.splus1) / (prod.pbplus0.s - prod.pb00.s)
      prob10to11 <- 0

      denominator <- 1 - prod.pbplus0.s - prod.pb0plus.s + prod.pb00.s
      prob11to00 <- (prod.pb00.splus1 - prod.pb00.splus1 * prod.pbplus0.bar.splus1 -
                       prod.pb00.splus1 * prod.pb0plus.bar.splus1 +  prod.pb00.s) / denominator
      prob11to01 <- (prod.pb0plus.splus1 - prod.pb00.splus1 - prod.pb0plus.s +
                       prod.pb00.splus1 * prod.pb0plus.bar.splus1) / denominator
      prob11to10 <- (prod.pbplus0.splus1 - prod.pb00.splus1 - prod.pbplus0.s +
                       prod.pb00.splus1 * prod.pbplus0.bar.splus1 ) / denominator
      prob11to11 <- (1 - prod.pbplus0.splus1 - prod.pb0plus.splus1 + prod.pb00.splus1) / denominator
      # prob11to00 + prob11to01 + prob11to10 + prob11to11   # = 1 as it should

      trans.prob.list[[s]] <- c(prob00to00, prob00to01, prob00to10, prob00to11,
                                prob01to00, prob01to01, prob01to10, prob01to11,
                                prob10to00, prob10to01, prob10to10, prob10to11,
                                prob11to00, prob11to01, prob11to10, prob11to11)
    }

    trans.prob <- rep(x = trans.prob.list[[1]], each = 2 * 4^(S-3))
    # This is for stage 1 to stage 2. There are S - 3 transitions remaining that 
    #   each have 4 possible values.
    # At stage S, there are 2 possible values for g_Sak.
    if (S > 3) {
      for (s in 2:(S-2)) {
        trans.prob <- rep(x = rep(x = trans.prob.list[[s]], each = 2 * 4^(S-s-2)), times = 4^(s-1)) * trans.prob
        # For rep(x = trans.prob.list[[s]], each = 2 * 4^(S-s-2)): There are 
        #   (S-s-2) stages remaining that each have 4 possible values.
        # For rep(x = rep(... ), times = 4^(s-1)): There are 4^(s-1) stages 
        #   prior to this.
        # Example with S = 4: The "each = 2 * 4^(S-s-2)" = 2 due to two 
        #   possible values of g4ak. Also, the "times = 4^(s-1)" is due to the 
        #   four possible values of g~_1 = (g~_11, g~_12).
      }
    }
  }
  # Example ordering for S = 3
  #  g1k g1kbar g2k g2kbar g3ak
  #   0    0     0     0    0
  #   0    0     0     0    1  each = 2 - notice duplicates for each value of g3ak
  #   0    0     0     1    0
  #   0    0     0     1    1
  #   0    0     1     0    0
  #   0    0     1     0    1
  #   0    0     1     1    0
  #   0    0     1     1    1

  #   0    1     0     0    0
  #   0    1     0     0    1
  # ...

  # Need to include new transition probabilities in the sum
  probYak <- sum(test.error.part * prob.Sminus1.k.kbar.to.Sk * trans.prob * theta)

  ##################################
  # Calculate PSp
  PSe <- PSeAllStages(joint.p = joint.p, group.mem = group.mem, a = a, k = k,
                      Se = Se, Sp = Sp, ordering = ordering)

  # From PSeAllStages() - could have this function return a list instead and 
  #   grab these values out of it
  index.pa.1andPlus <- all.one.k(k = k, K = K, ordering = ordering)
  pa.1andPlus <- sum(joint.p[index.pa.1andPlus, a])    # For k = 1, p_a1+

  PSp <- 1 - (probYak - PSe*pa.1andPlus)/(1-pa.1andPlus)

  PPPV <- pa.1andPlus * PSe / (pa.1andPlus * PSe + (1 - pa.1andPlus) * (1 - PSp))
  PNPV <- (1 - pa.1andPlus) * PSp / ((1 - pa.1andPlus) * PSp + pa.1andPlus * (1 - PSe))

  list(PSe = PSe, PSp = PSp, PPPV = PPPV, PNPV = PNPV)
}
###############################################################################




# Start new functions for hierarchical testing with two diseases
###############################################################################
# functions from Chris's Section 4 materials for Bilder et al. (2018)

# Create group membership matrix using information from parts() in partition 
#   package

# Create a group membership matrix M that can be used with E(T)
# one.config = one column from parts()
# Assumes groups with only one individual are not until end of vector used to 
#   denote a stage

# Brianna Hitt - 04/05/2019
# add an IF statement to create a group membership matrix for a configuration 
#   of all 1's

group.membership.mat <- function(one.config) {
  
  I <- sum(one.config)
  
  # First group of second stage
  save.2nd.stage <- rep(x = 1, times = one.config[1])
  
  # Remaining groups of second stage
  for(i in 2:length(one.config[one.config != 0])) {
    save.2nd.stage <- c(save.2nd.stage, rep(i, times = one.config[i]))
  }
  # save.2nd.stage
  
  # Find 3rd stage for configurations of all 1s
  if(isTRUE(all.equal(one.config, rep(x=1, times=I)))){
    save.third.stage <- rep(x=NA, times=I)
  } else{
    # Find 3rd stage for all other configurations
    counts <- table(save.2nd.stage)  # Number per group stage 2
    group.num <- as.numeric(names(counts))  # Numerical group names
    # Just group numbers with more than one individual
    group.grtr1 <- group.num[counts > 1]  
    # Groups with more than one individual
    # save.2nd.stage <= max(group.grtr1)  
    # Number of individuals tested in stage 3
    ind.test <- sum(save.2nd.stage <= max(group.grtr1)) 
    save.third.stage <- rep(x = NA, times = I)  # Stage 3 filled with NAs
    # Replace NAs with group numbers
    save.third.stage[1:ind.test] <- 1:ind.test   
    # save.third.stage
  }
  
  matrix(data = c(rep(x = 1, times = I), save.2nd.stage, save.third.stage),
         nrow = 3, ncol = I, byrow = TRUE, 
         dimnames = list(Stage = 1:3, Individual = 1:I))
}

###############################################################################
# Find joint probabilities

# Simulate P
joint.p.create <- function(alpha, I = 10) {
  p.unordered <- t(rdirichlet(n = I, shape = alpha))
  # Assumes p00 is given first
  p.ordered <- p.unordered[,order(1 - p.unordered[1,])]
  p.ordered
}

# Create row labels for the 0-1 combinations used with joint.p matrix
joint.p.row.labels <- function(ordering) {
  save.it <- character(length = nrow(ordering))
  for (i in 1:nrow(ordering)) {
    save.it[i] <- paste(ordering[i,], collapse="")
  }
  save.it
}




# wrapper functions for ET.all.stages.new and PSePSpAllStages

# write a wrapper function for ET.all.stages.new that calculates ET for
#   informative two-stage testing - do this by calling ET.all.stages.new
#   two separate times - once for groups that are tested in both stages
#   and once for groups that are tested only in the first stage of the algorithm

ET.all.algs <- function(joint.p, group.mem, Se, Sp, ordering){
  
  # find the number of groups in the 1st stage of testing
  max.groups <- max(group.mem[1,])
  
  ET.save <- numeric(max.groups)
  for(group.num in 1:max.groups){
    
    # account for individuals tested individually 
    if(sum(group.mem[1,]==group.num) > 1){
      # ET.all.stages.new() may not work properly when the stage 1 group number 
      #   is different from 1.
      # Thus, a quick fix is to change the group number to 1 when needed.
      one.group <- group.mem[, group.mem[1,]==group.num]
      one.group[1,] <- rep(x=1, times=ncol(one.group))
      
      ET.part <- ET.all.stages.new(joint.p=joint.p[, group.mem[1,]==group.num], 
                                   group.mem=one.group, Se=Se, Sp=Sp, 
                                   ordering=ordering)
      ET.save[group.num] <- ET.part$Expt
    } else{
      ET.save[group.num] <- 1
    }
    ET <- sum(ET.save)
  }
  
  # find c, I.sj, and m.sj for the overall configuration
  config.info <- ET.all.stages.new(joint.p=joint.p, group.mem=group.mem, Se=Se, 
                                   Sp=Sp, ordering=ordering)[2:4]
  
  c("Expt"=ET, config.info)
}


# write a wrapper function for ET.all.stages.new that calculates ET for
#   informative two-stage testing - do this by calling ET.all.stages.new
#   two separate times - once for groups that are tested in both stages
#   and once for groups that are tested only in the first stage of the algorithm

# write a wrapper function for PSePSpAllStages that calculates accuracy measures
#   for all algorithms - want to allow informative two-stage testing where 
#   individuals may not be tested in every stage of the algorithm - do this by 
#   calling PSePSpAllStages only for individuals who are tested twice - assign 
#   values for individuals who are tested only once
# also want to calculate measures for both diseases (k=1 and k=2) and for 
#   as many individuals as requested by the user

PSePSp.all.algs <- function(joint.p=joint.p, group.mem=group.mem, Se=Se, Sp=Sp, 
                            a=1, ordering=ordering){
  index.p.1plus <- all.one.k(k=1, K=2, ordering=ordering)
  p.1plus <- colSums(joint.p[index.p.1plus,])    
  
  index.p.plus1 <- all.one.k(k=2, K=2, ordering=ordering)
  p.plus1 <- colSums(joint.p[index.p.plus1,])
  
  num.ind <- ncol(group.mem)
  accuracy.dis1 <- matrix(data=NA, nrow=num.ind, ncol=5)
  accuracy.dis2 <- matrix(data=NA, nrow=num.ind, ncol=5)
  count <- 1
  
  for(i in 1:num.ind){
    # this condition, for informative two-stage testing only
    if(is.na(group.mem[2,i])){
      pa.1plus <- p.1plus[i]    # For k = 1, p_a1+
      PSe1 <- Se[1,1]
      PSp1 <- Sp[1,1]
      PPPV1 <- (pa.1plus*PSe1)/((pa.1plus*PSe1) + ((1-pa.1plus)*(1-PSp1)))
      PNPV1 <- ((1-pa.1plus)*PSp1)/(((1-pa.1plus)*PSp1) + (pa.1plus*(1-PSe1)))
      
      pa.plus1 <- p.plus1[i]    # For k = 2, p_a+1
      PSe2 <- Se[2,1]
      PSp2 <- Sp[2,1]
      PPPV2 <- (pa.plus1*PSe2)/((pa.plus1*PSe2) + ((1-pa.plus1)*(1-PSp2)))
      PNPV2 <- ((1-pa.plus1)*PSp2)/(((1-pa.plus1)*PSp2) + (pa.plus1*(1-PSe2)))
      
      accuracy.dis1[count,] <- c(i, PSe1, PSp1, PPPV1, PNPV1)
      accuracy.dis2[count,] <- c(i, PSe2, PSp2, PPPV2, PNPV2)
    } else{
      res.dis1 <- PSePSpAllStages(joint.p=joint.p, group.mem=group.mem,
                                  a=i, k=1, Se=Se, Sp=Sp, ordering=ordering)
      res.dis2 <- PSePSpAllStages(joint.p=joint.p, group.mem=group.mem,
                                  a=i, k=2, Se=Se, Sp=Sp, ordering=ordering)
      
      accuracy.dis1[count,] <- c(i, res.dis1$PSe, res.dis1$PSp, res.dis1$PPPV, 
                                 res.dis1$PNPV)
      accuracy.dis2[count,] <- c(i, res.dis2$PSe, res.dis2$PSp, res.dis2$PPPV, 
                                 res.dis2$PNPV)
    }
    count <- count + 1
  }
  colnames(accuracy.dis1) <- c("Individual", "PSe", "PSp", "PPPV", "PNPV")
  PSe1.vec <- accuracy.dis1[,2]
  PSp1.vec <- accuracy.dis1[,3]
  overall.PSe1 <- sum(p.1plus*PSe1.vec)/sum(p.1plus)
  overall.PSp1 <- sum((1-p.1plus)*PSp1.vec)/sum(1-p.1plus)
  overall.PPPV1 <- sum(p.1plus*PSe1.vec)/sum((p.1plus*PSe1.vec) + ((1-p.1plus)*(1-PSp1.vec)))
  overall.PNPV1 <- sum((1-p.1plus)*PSp1.vec)/sum(((1-p.1plus)*PSp1.vec) + (p.1plus*(1-PSe1.vec)))
  
  colnames(accuracy.dis2) <- c("Individual", "PSe", "PSp", "PPPV", "PNPV")
  PSe2.vec <- accuracy.dis2[,2]
  PSp2.vec <- accuracy.dis2[,3]
  overall.PSe2 <- sum(p.plus1*PSe2.vec)/sum(p.plus1)
  overall.PSp2 <- sum((1-p.plus1)*PSp2.vec)/sum(1-p.plus1)
  overall.PPPV2 <- sum(p.plus1*PSe2.vec)/sum((p.plus1*PSe2.vec) + ((1-p.plus1)*(1-PSp2.vec)))
  overall.PNPV2 <- sum((1-p.plus1)*PSp2.vec)/sum(((1-p.plus1)*PSp2.vec) + (p.plus1*(1-PSe2.vec)))
  
  accuracy.dis1 <- accuracy.dis1[order(accuracy.dis1[,1]),]
  accuracy.dis2 <- accuracy.dis2[order(accuracy.dis2[,1]),]
  
  list("Disease1"=list("Individual"=accuracy.dis1[a,],
                       "Overall"=matrix(data=c(overall.PSe1, overall.PSp1, 
                                               overall.PPPV1, overall.PNPV1), 
                                        nrow=1, ncol=4, 
                                        dimnames=list("", c("PSe", "PSp", "PPPV", "PNPV")))),
       "Disease2"=list("Individual"=accuracy.dis2[a,], 
                       "Overall"=matrix(data=c(overall.PSe2, overall.PSp2, 
                                               overall.PPPV2, overall.PNPV2),
                                        nrow=1, ncol=4, 
                                        dimnames=list("", c("PSe", "PSp", "PPPV", "PNPV")))))
}




# Start TwoDisease.Hierarchical() function
###############################################################################
# Find E(T) and accuracy measures for hierarchical group testing with
#   multiplex assays (K=2) in heterogeneous populations
# This is a wrapper function that calls both ET.all.algs and PSePSp.all.algs
# Author: Brianna Hitt
# Date: 03-13-2019

###############################################################################
# updated 04-13-19 with new wrapper functions (rather than using 
#   ET.all.stages.new and PSePSpAllStages directly)

TwoDisease.Hierarchical <- function(joint.p, group.mem, Se, Sp, ordering, a=1, 
                                    accuracy=TRUE){
  
  # save results from ET.all.stages.new() and PSePSpAllStages()
  ET.results <- ET.all.algs(joint.p=joint.p, group.mem=group.mem,
                            Se=Se, Sp=Sp, ordering=ordering)
  
  # print results
  if(accuracy==TRUE){
    acc.results <- PSePSp.all.algs(joint.p=joint.p, group.mem=group.mem,
                                   Se=Se, Sp=Sp, ordering=ordering, a=a)
    
    c(ET.results, acc.results)
    # results <- list(ExpT=ET.results$Expt, c=ET.results$c,
    #                 I.sj=ET.results$I.sj, m.sj=ET.results$m.sj,
    #                 Accuracy=Acc.results)
  } else{
    ET.results
    # results <- list(ExpT=ET.results$Expt, c=ET.results$c,
    #                 I.sj=ET.results$I.sj, m.sj=ET.results$m.sj)
  }
}

#
