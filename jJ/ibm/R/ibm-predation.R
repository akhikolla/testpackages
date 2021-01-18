
localPredationInteractions = function(pop, R, N=NULL, P=NULL) {
  pred = localIntegration(x=.getMatrix(x=pop$N, size=N), 
                          y=.getMatrix(x=pop$P, size=P), 
                          R1=R$N, R2=R$P)
  NP = pred$x 	# number of predators near to each prey
  PN = pred$y		# number of preys near to each predator
  return(list(NP=NP, PN=PN))
}
