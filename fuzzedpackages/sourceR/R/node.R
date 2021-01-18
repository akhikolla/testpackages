# Testing reference classes

#
#' Node
#'
#' This is a base class representing a node in a DAG. Is not intended to be used by a regular user.  Developers only here!
#'
#' @docType class
#' @name Node
#' @importFrom R6 R6Class
#' @keywords DAG node
#' @return Object of \code{\link{Node}}
#' @format Object of \code{\link{R6Class}} with methods for constructing a DAG.
#' @field parents a list of parent nodes
#' @field children a list of child nodes
#' @field name a tag name applied to the node
#' @section Methods:
#' \describe{
#'   \item{\code{new(parents = list(), children = list(), name)}}{creates a new \link{Node} with parent nodes, child nodes, and a name.}
#'   \item{\code{logDensity()}}{calculate the log probability density/mass function evaluated at the current node value.}
#'   \item{\code{addChild(node)}}{add \code{node} as a child.  Returns \code{node}.}
#'   \item{\code{addParent(node)}}{add \code{node} as a parent.  Returns \code{node}.}
#'   \item{\code{removeParent(name)}}{remove the parent node named \code{name}.  Returns \code{node}.}
#'   \item{\code{removeChild(name)}}{remove the child node named \code{name}.  Returns \code{node}.}
#'   }
Node <- R6::R6Class(
  "Node",
  public = list(
    name = NA,
    parents = NA,
    children = NA,
    initialize = function(parents = list(), children = list(), name) {
      self$name <- name
      self$parents <- parents
      self$children <- children
    },
    logDensity = function() {
      "Return the log probability density|mass function"
      1
    },
    addChild = function(node, name) {
      if(missing(name)) {
        name = node$name
      }
      self$children[[name]] <- node
      node$parents[[self$name]] <- self
      node
    },
    addParent = function(node, name) {
      if(missing(name)) name <- node$name
      self$parents[[name]] <- node
      node$children[[self$name]] <- self
      node
    },
    removeChild = function(name) {
      node <- self$children[[name]]
      self$children[[name]] <- NULL
      node
    },
    removeParent = function(name) {
      node <- self$parents[[name]]
      self$parents[[name]] <- NULL
      node
    }
  )
)



#' StochasticNode
#'
#' Represents a stochastic node in a DAG
#'
#' Derived from \link{Node}, please see base class documentation.
#'
#' @docType class
#' @name StochasticNode
#'
#' @importFrom R6 R6Class
#' @export
#' @keywords DAG node
#' @return Object of \code{\link{StochasticNode}}
#' @format Object of \code{\link{R6Class}} with methods for constructing a DAG.
#'
#' @field data contains the node's data
#' @section Methods:
#' \describe{
#'   \item{\code{logPosterior()}}{return the value of the log posterior distribution of the node.}
#'   \item{\code{getData()}}{returns the node's data.}
#'   }
StochasticNode <- R6::R6Class(
  "StochasticNode",
  inherit = Node,
  public = list(
    data = NA,
    logDensity = function()
      1,
    logPosterior = function()
      self$logDensity() + sum(sapply(self$children, function(node)
        node$logDensity())),
    getData = function()
      self$data
  )
)


#' DataNode
#'
#' Represents a static data node in a DAG.
#'
#' Derived from \link{Node}, please see base class documentation.
#'
#' @docType class
#' @name DataNode
#'
#' @importFrom R6 R6Class
#' @export
#' @keywords DAG node
#' @return Object of \code{\link{DataNode}}
#' @format Object of \code{\link{R6Class}} with methods for constructing a DAG.
#' @field data the data
#'
#' @section Methods:
#' \describe{
#'   \item{\code{getData()}}{returns the node's data.}
#'   }
DataNode <- R6::R6Class("DataNode",
                    inherit = Node,
                    public = list(
                      data = NA,
                      initialize = function(data, name) {
                        super$initialize(name=name)
                        self$data <- data
                      },
                      getData = function()
                        self$data
                    )
)



#' FormulaNode
#'
#' Represents a formula node in a DAG.  Inherit from this node to specify some kind of
#' formula within the DAG, e.g. a linear predictor and/or link function.  Override the
#' FormulaNode$getData() method to apply your own function.
#'
#' Derived from \link{Node}, please see base class documentation.
#'
#' @docType class
#' @name FormulaNode
#'
#' @importFrom R6 R6Class
#' @export
#' @keywords DAG node
#' @return Object of \code{\link{FormulaNode}}
#' @format Object of \code{\link{R6Class}} with methods for constructing a DAG.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{getData()}}{returns the node's transformed data.}
#'   }
FormulaNode <- R6::R6Class(
  "FormulaNode",
  inherit = Node,
  public = list(
    getData = function() {
      sapply(self$children, function(child)
        child$getData())
    },
    logDensity = function()
      sapply(self$children, function(node)
        node$logDensity())
  )
)



#' PoissonNode
#'
#' Represents a Poisson distribution node in a DAG
#'
#' Derived from \link{StochasticNode}, please see base class documentation.
#'
#' @docType class
#' @name PoissonNode
#'
#' @importFrom R6 R6Class
#' @export
#' @keywords DAG node
#' @return Object of \code{\link{PoissonNode}}
#' @format Object of \code{\link{R6Class}} with methods for constructing a DAG.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new(data, lambda, offset)}}{create a PoissonNode, with mean \link{Node} \code{lambda}, and offset \link{Node} \code{offset}.}
#'   }
PoissonNode <- R6::R6Class(
  "PoissonNode",
  inherit = StochasticNode,
  public = list(
    initialize = function(data, lambda = NULL, offset = NULL, name) {
      super$initialize(name = name)
      self$data <- data
      if (!is.null(lambda))
        self$addParent(lambda, 'lambda')
      if (!is.null(offset))
        self$addParent(offset, 'offset')
    },
    logDensity = function() {
      lambda <- as.data.frame(sapply(self$parents,
                                     function(parent) parent$getData()))
      lambda <- apply(lambda, 1, prod)
      sum(dpois(self$data, lambda, log = T))
    }
  )
)


#' GammaNode
#'
#' Represents a Gamma distribution node in a DAG.  Requires parent nodes for shape and rate respectively as
#' specified in \link{dgamma}.
#'
#' Derived from \link{StochasticNode}, please see base class documentation.
#'
#' @docType class
#' @name GammaNode
#'
#' @importFrom R6 R6Class
#' @export
#' @keywords DAG node
#' @return Object of \code{\link{GammaNode}}
#' @format Object of \code{\link{R6Class}} with methods for constructing a DAG.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new(data, shape=1, rate=1)}}{Create a Gamma node with data \code{data} and Nodes
#'    \code{shape} and \code{rate} as specified in \link{dgamma}.}
#'   }
GammaNode <- R6::R6Class(
  "GammaNode",
  inherit = StochasticNode,
  public = list(
    initialize = function(data, shape, rate, name) {
      super$initialize(name = name)
      self$data <- data
      self$addParent(shape, 'shape')
      self$addParent(rate, 'rate')
      #self$parents <- list(shape=shape, rate=rate)
    },
    logDensity = function()
      dgamma(self$data, self$parents$shape$getData(),
             self$parents$rate$getData(), log = T)
  )
)


#' DirichletNode
#'
#' Represents a d-dimensional Dirichlet distribution node in a DAG.
#'
#' Derived from \link{StochasticNode}, please see base class documentation.
#'
#' @docType class
#' @name DirichletNode
#'
#' @importFrom R6 R6Class
#' @export
#' @keywords DAG node
#' @return Object of \code{\link{DirichletNode}}
#' @format Object of \code{\link{R6Class}} with methods for constructing a DAG.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new(data, alpha)}}{Create a DirichletNode with data vector \code{data} (length > 1) and parameter vector
#'    \code{alpha}}.
#'   }
DirichletNode <- R6::R6Class(
  "DirichletNode",
  inherit = StochasticNode,
  public = list(
    initialize = function(data, alpha, name) {
      super$initialize(name = name)
      self$data <- data
      self$addParent(alpha, name = 'alpha')
    },
    logDensity = function() {
      lgamma(sum(self$parents$alpha$getData())) -
        sum(lgamma(self$parents$alpha$getData())) +
        sum((self$parents$alpha$getData() - 1) * log(self$data))
    }
  )
)


DirichletNode2 <- R6::R6Class(
  "DirichletNode2",
  inherit = StochasticNode,
  public = list(
    initialize = function(data, alpha, name) {
      super$initialize(name = name)
      self$data <- data
      self$addParent(alpha, name='alpha')
    },
    logDensity = function() {
      sum(dgamma(self$data, shape=self$parents$alpha$getData(), rate=1, log=T))
    },
    getData = function() {
      self$data / sum(self$data)
    }
  )
)


#' DirichletProcessNode
#'
#' Represents a Dirichlet process as a single node in a DAG.
#'
#' Derived from \link{StochasticNode}, please see base class documentation.
#'
#' @docType class
#' @name DirichletProcessNode
#'
#' @importFrom R6 R6Class
#' @export
#' @keywords DAG node
#' @return Object of \code{\link{DirichletProcessNode}}
#' @format Object of \code{\link{R6Class}} with methods for constructing a DAG.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new(theta, s, alpha, base, ...)}}{Create a DirichletProcessNode with value vector
#'   \code{theta} (length > 1), initial grouping vector \code{s}, concentration parameter \code{alpha}, and
#'   base distribution \code{base}.  Base should be a distribution function (dnorm, dgamma, etc) whose parameters
#'   are specified in \code{...}.}
#'   }
DirichletProcessNode <- R6::R6Class( # TODO: Make this accept a generic base distribution, not just Gamma
  "DirichletProcessNode",
  inherit = StochasticNode,
  public = list(
    theta = NA,
    s = NA,
    base = NA,
    conc = NA,
    baseShape = NA,
    baseRate = NA,
    idBucket = NA,
    initialize = function(theta, s, alpha, base, shape, rate, name) {
      super$initialize(name = name)
      self$conc <- alpha
      self$base <- dgamma
      self$baseShape <- shape
      self$baseRate <- rate
      self$idBucket <- queue()
      for(i in 1:length(s)) enqueue(self$idBucket, as.character(i))
      keys = replicate(length(theta), dequeue(self$idBucket))
      self$theta <- HashTable(keys, theta)
      self$s <- keys[s]
    },
    getData = function()
      self$theta$find(self$s),
    getDensity = function(i)
      sum(self$base(self$theta$find(self$s[i]),shape=self$baseShape, rate=self$baseRate ,log = T))
  )
)

