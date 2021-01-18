################
## All Generics
################

setGeneric("marg", function(object) standardGeneric("marg"))
setGeneric("cond", function(object) standardGeneric("cond"))
setGeneric("mcor", function(object) standardGeneric("mcor"))
setGeneric("pcor", function(object) standardGeneric("pcor"))
setGeneric("postExpSigma", function(object, vars.method="eb") standardGeneric("postExpSigma"))
setGeneric("postExpOmega", function(object, vars.method="eb") standardGeneric("postExpOmega"))
setGeneric("plotML", function(object, ...) standardGeneric("plotML"))
setGeneric("plotCor", function(object, type = object@type, order = 'original', by = "marginal") standardGeneric("plotCor"))
setGeneric("plotAdj", function(object, type="conditional", order = 'original') standardGeneric("plotAdj"))
setGeneric("bgraph", function(object) standardGeneric("bgraph"))
setGeneric("ugraph", function(object) standardGeneric("ugraph"))
