########################################################################################
##  Copyright (C) 2017,  Constantinos Tsirogiannis.  Email: tsirogiannis.c@gmail.com  ##
##                                                                                    ##
##  This file is part of CNull.                                                       ##
##                                                                                    ##
##  CNull is free software: you can redistribute it and/or modify                     ##
##  it under the terms of the GNU General Public License as published by              ##
##  the Free Software Foundation, either version 3 of the License, or                 ##
##  (at your option) any later version.                                               ##
##                                                                                    ##
##  CNull is distributed in the hope that it will be useful,                          ##
##  but WITHOUT ANY WARRANTY; without even the implied warranty of                    ##
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                     ##
##  GNU General Public License for more details.                                      ##
##                                                                                    ##
##  You should have received a copy of the GNU General Public License                 ##
##  along with CNull.  If not, see <http://www.gnu.org/licenses/>                     ##
########################################################################################

####################################
####################################
### Permutation sampling methods ###
####################################
####################################

#############################################################
#############################################################
################# Alpha diversity functions #################
#############################################################
#############################################################

permutation.communities.a = function(matrix,reps=1000)
{
  samples = communities_permutation_sampling_alpha(as.matrix(matrix), reps)

  if(nrow(samples)==0)
    return()
 
  colnames(samples) = colnames(matrix)  
  return (samples)

} # permutation.communities.a(...)


permutation.random.values.a = function(matrix,f,args,reps=1000)
{
  samples = communities_permutation_sampling_alpha(as.matrix(matrix), reps)

  if(nrow(samples)==0)
    return()

  colnames(samples) = colnames(matrix)

  vals = f(samples,args)

  return (vals)

} # permutation.random.values.a(...)

permutation.moments.a = function(matrix,f,args,reps=1000)
{
  samples = communities_permutation_sampling_alpha(as.matrix(matrix), reps)

  if(nrow(samples)==0)
    return()

  colnames(samples) = colnames(matrix)

  vals = f(samples,args)
  mn = mean(vals)
  vr = var(vals)

  return (list(mean=mn,variance=vr))

} # permutation.moments.a(...)


permutation.pvalues.a = function(matrix,f,args,reps=1000)
{
  samples = communities_permutation_sampling_alpha(as.matrix(matrix), reps)

  if(nrow(samples)==0)
    return()

  colnames(samples) = colnames(matrix)

  original.vals = f(matrix,args)
  randomized.vals = f(samples,args)
  
  pvals = compute_pvalues(original.vals,randomized.vals)

  return (pvals)

} # permutation.pvals.a(...)



########################################################
########################################################
## Beta diversity functions, single matrix, all pairs ##
########################################################
########################################################

permutation.communities.b = function(matrix,reps=1000)
{
  samples = communities_permutation_sampling_beta_interleaved_matrices(as.matrix(matrix), reps)

  if(nrow(samples)==0)
    return()

  colnames(samples) = colnames(matrix)

  return (samples)

} # permutation.communities.b(matrix,...)

permutation.random.values.b = function(matrix,f,args,reps=1000)
{
  samples = communities_permutation_sampling_beta(as.matrix(matrix), reps)
  
  samples.a = samples[[1]]
  samples.b = samples[[2]]

  if(nrow(samples.a)==0)
    return()

  vals = vector(mode="numeric", length=reps)

  for(i in 1:reps)
  { 
    tmp.matrix = rbind(samples.a[i,],samples.b[i,])  
    colnames(tmp.matrix) = colnames(matrix)
    res = f(tmp.matrix,args)
    vals[i] = res[1,1] 
  }

  return (vals)

} # permutation.random.values.b(matrix, ...)

permutation.moments.b = function(matrix,f,args,reps=1000)
{
  samples = communities_permutation_sampling_beta(as.matrix(matrix), reps)
  
  samples.a = samples[[1]]
  samples.b = samples[[2]]

  if(nrow(samples.a)==0)
    return()

  vals = vector(mode="numeric", length=reps)

  for(i in 1:reps)
  { 
    tmp.matrix = rbind(samples.a[i,],samples.b[i,])  
    colnames(tmp.matrix) = colnames(matrix)
    res = f(tmp.matrix,args)
    vals[i] = res[1,1] 
  }

  mn = mean(vals)
  vr = var(vals)

  return (list(mean=mn,variance=vr))

} # permutation.moments.b(matrix, ...)

permutation.pvalues.b = function(matrix, f, args, observed.vals, reps=1000)
{
  samples = communities_permutation_sampling_beta(as.matrix(matrix), reps)
  
  samples.a = samples[[1]]
  samples.b = samples[[2]]

  if(nrow(samples.a)==0)
    return()
  
  randomized.vals = vector(mode="numeric", length=reps)

  for(i in 1:reps)
  { 
    tmp.matrix = rbind(samples.a[i,],samples.b[i,])  
    colnames(tmp.matrix) = colnames(matrix)
    res = f(tmp.matrix,args)
    randomized.vals[i] = res[1,1] 
  }

  pvals = compute_pvalues(observed.vals,randomized.vals)

  return (pvals)

} # permutation.pvals.b(matrix, ...)

###################################
###################################
### Unit-based sampling methods ###
###################################
###################################

#############################################################
#############################################################
################# Alpha diversity functions #################
#############################################################
#############################################################

individual.based.communities.a = function(matrix,reps=1000)
{
  samples = communities_individual_based_sampling_alpha(as.matrix(matrix), reps) 

  if(nrow(samples)==0)
    return()

  colnames(samples) = colnames(matrix)  
  return (samples)

} # individual.based.communities.a(...)


individual.based.moments.a = function(matrix,f,args,reps=1000)
{
  samples = communities_individual_based_sampling_alpha(as.matrix(matrix), reps)

  if(nrow(samples)==0)
    return()

  colnames(samples) = colnames(matrix)

  vals = f(samples,args)
  mn = mean(vals)
  vr = var(vals)

  return (list(mean=mn,variance=vr))

} # individual.based.moments.a(...)

individual.based.random.values.a = function(matrix,f,args,reps=1000)
{
  samples = communities_individual_based_sampling_alpha(as.matrix(matrix), reps)

  if(nrow(samples)==0)
    return()

  colnames(samples) = colnames(matrix)

  vals = f(samples,args)

  return (vals)

} # individual.based.random.values.a(...)

individual.based.pvalues.a = function(matrix,f,args,reps=1000)
{
  samples = communities_individual_based_sampling_alpha(as.matrix(matrix), reps)

  if(nrow(samples)==0)
    return()

  colnames(samples) = colnames(matrix)

  original.vals = f(matrix,args)
  randomized.vals = f(samples,args)
  
  pvals = compute_pvalues(original.vals,randomized.vals)

  return (pvals)

} # individual.based.pvals.a(...)

########################################################
########################################################
## Beta diversity functions, single matrix, all pairs ##
########################################################
########################################################

individual.based.communities.b = function(matrix,reps=1000)
{
  samples = communities_individual_based_sampling_beta_interleaved_matrices(as.matrix(matrix), reps)

  if(nrow(samples)==0)
    return()

  colnames(samples) = colnames(matrix)

  return (samples)

} # individual.based.communities.b(matrix,...)

individual.based.random.values.b 