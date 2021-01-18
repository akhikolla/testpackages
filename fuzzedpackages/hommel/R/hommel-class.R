setClass("hommel",
         representation(
           p = "numeric",               #(m) store vector of p-values (possibly with names)
           jumpalpha = "numeric",       #(m+1) stores vector of alphas where h(alpha) jumps
           sorter = "integer",          #(m) permutation to sort p
           adjusted = "numeric",        #(m) stores adjusted p-values 
           simesfactor = "numeric",     #(m+1) stores simesfactor (start at 0)
           simes = "logical"            # stores whether adjusted p-values ar calculated based on a Simes test (if TRUE) or on Hommel's test from 1983 (if FALSE)
         )
)

setMethod("p.adjust", "hommel", function(p, method, n) {
  if(!missing(method)) {
    if (p@simes)
      method = "Hommel's method (Simes inequality assumed)"
    else
      method = "Hommel's method (Simes inequality not assumed)"
    stop(paste("Method argument not used. Method is", method))
  }
  if (!missing(n)) {
    stop("argument n not used.")
  }
  p@adjusted
})

setMethod("show", "hommel", function(object) {
  cat("A hommel object for", length(object@p), "hypotheses.\n")
  if (object@simes)
    cat("Simes inequality is assumed.\n")
  else
    cat("Simes inequality is not assumed.\n")
  cat("Use p.adjust(), discoveries() or localtest() to access this object.\n")
})


setMethod("summary", "hommel", function(object, alpha = 0.05, ...) {
  show(object)
  
  cat("\n")
  
  disc <- discoveries(object, alpha=alpha)
  
  cat("With", 1-alpha, "confidence: at least", sum(disc), "discoveries.\n")
  cat(sum(object@adjusted <= alpha), " hypotheses with adjusted p-values below ", alpha, ".\n", sep="")
})

