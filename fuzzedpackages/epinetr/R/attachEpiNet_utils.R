# Generate an EpiNet object stochastically
buildNetwork <- function(n, k, m, additive, scaleFree, pop) {
  retval <- constructEpiNet(n, k)

  degree <- rep(0, n)
  network <- NULL
  tempk <- sort(k)

  if (length(additive) == 1) {

    # Use a number of randomly selected QTLs to be additive
    if (!is.numeric(additive) || additive >= n || additive < 0 || additive %% 1 !=
      0) {
      stop("additive must be a non-negative integer less than the number of QTL")
    }

    additive <- sample(n, additive)
  } else if (length(additive) > 1) {

    # Verify that additive IDs are valid
    if (!all(additive %in% pop$map$SNP)) {
      stop("Invalid additive QTL IDs")
    }

    # Translate IDs to indices
    additive <- match(additive, pop$map$SNP)

    # Verify that all indices are valid QTLs
    if (!all(additive %in% pop$qtl)) {
      stop("Invalid additive QTL IDs")
    }

    # Translate indices to QTLs
    additive <- match(additive, pop$qtl)
  } else {
    stop("additive vector cannot have length 0")
  }

  for (i in tempk) {
    incmat <- rn(degree, k = i, m = m, additive = additive, scaleFree = scaleFree)
    if (!is.null(incmat)) {
      network <- cbind(network, incmat)
      degree <- degree + rowSums(incmat)
    }
  }

  # retval$Incidence = apply(network, 2, pow3vec)
  retval$Incidence <- network
  retval$Seeds <- sample((2^20):(2^30), ncol(retval$Incidence))

  return(retval)
}


# Build an EpiNet object based on a user-supplied incidence matrix
userNetwork <- function(incmat) {
  # Ensure the matrix is composed of numbers
  if (!is.numeric(incmat[1])) {
    stop("Incidence matrix must be numeric")
  }

  net <- list()
  class(net) <- "EpiNet"

  # Limit values to 0 and 1
  dimensions <- dim(incmat)
  incmat <- round(abs(incmat))
  incmat[incmat > 0] <- 1
  incmat <- matrix(as.integer(incmat), nrow = dimensions[1], ncol = dimensions[2])
  net$Incidence <- incmat

  # Get the seeds for the interactions
  net$Seeds <- sample((2^20):(2^30), dimensions[2])

  return(net)
}


# Construct a network incidence matrix
rn <- function(degree, m = 1, k = 2, additive = 0, scaleFree = FALSE) {
  if (!is.numeric(degree) || any(degree %% 1 != 0) || any(degree < 0) ||
    length(degree) < 2) {
    stop("degree must be a vector of non-negative integers of length greater than 1")
  }

  # Determine number of nodes
  n <- length(degree)

  if (!is.numeric(m) || m %% 1 != 0 || m < 1 || length(m) != 1) {
    stop("m must be a single positive integer")
  }

  if (!is.numeric(k) || k %% 1 != 0 || k < 2 || k > n || length(k) != 1) {
    stop("k must be an integer between 2 and the number of QTL inclusive")
  }

  # Start with enough unconnected nodes so that we can form an initial
  # connected graph, given parameters m and k
  initnodes <- 1
  while (choose(initnodes, k - 1) < m) initnodes <- initnodes + 1

  if (length(additive) > n - initnodes - 1) {
    stop("additive is too large, given m and k")
  }

  # Shuffle the nodes
  nodes <- sample(n, n)

  # Remove nodes considered to be purely additive
  nodes <- nodes[!(nodes %in% additive)]

  # Empty incidence matrix
  incmat <- matrix(0, n, (length(nodes) - initnodes) * m)

  # Add nodes to the graph
  column <- 1
  for (i in (initnodes + 1):(length(nodes))) {
    if (i == 2) {
      # It's simple enough to simply connect the first node to the second if
      # there are only two; this gets around combn taking a vector x of
      # length 1 as 1:x.
      combinations <- nodes[1]
    } else {
      # Find all the ways the new node can connect to the previous nodes.
      if (m == 1) {
        combinations <- sample(nodes[1:(i - 1)], k - 1)
      } else {
        combinations <- t(RcppAlgos::comboGeneral(nodes[1:(i - 1)], k - 1))

        # Get probabilities from degrees if scale-free
        pp <- apply(combinations, 2, function(x) sum(degree[x]))
        if (sum(pp) > 0 && scaleFree) {
          pp <- pp / sum(pp)

          # Sample from the different ways of connecting to the previous nodes
          combinations <- combinations[, sample(ncol(combinations), m,
                                                prob = pp
          )]
        } else {
          combinations <- combinations[, sample(ncol(combinations), m)]
        }
      }
    }

    for (j in 1:m) {
      # Add the new node
      incmat[nodes[i], column] <- 1
      if (!scaleFree)
        degree[nodes[i]] <- degree[nodes[i]] + 1

      # Connect it to the previous nodes
      if (is.null(dim(combinations)) && m == 1) {
        incmat[combinations, column] <- 1
        if (!scaleFree)
          degree[combinations] <- degree[combinations] + 1
      } else if (is.null(dim(combinations)) && k == 2) {
        incmat[combinations[j], column] <- 1
        if (!scaleFree)
          degree[combinations[j]] <- degree[combinations[j]] + 1
      } else {
        incmat[combinations[, j], column] <- 1
        if (!scaleFree)
          degree[combinations[, j]] <- degree[combinations[, j]] +
            1
      }
      column <- column + 1
    }
  }

  return(incmat)
}


# Construct an EpiNet object
constructEpiNet <- function(n, k) {
  if (!is.numeric(n)) {
    stop("n must be single integer greater than 2")
  }

  if (length(n) != 1 | n < 2 | n %% 1 != 0) {
    stop("n must be single integer greater than 2")
  }

  if (!is.numeric(k)) {
    stop("k must be a integer vector with all values between 2 and n inclusive")
  }

  if (length(k) < 1 | any(k %% 1 != 0) | any(k < 2) | any(k > n)) {
    stop("k must be a integer vector with all values between 2 and n inclusive")
  }

  net <- list()
  class(net) <- "EpiNet"

  return(net)
}
