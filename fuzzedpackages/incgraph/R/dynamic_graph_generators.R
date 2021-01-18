#' Generate a dynamic network
#'
#' @usage
#' generate.dynamic.network(
#'   model, amnt.nodes, amnt.edges, amnt.operations, trace = T, ...)
#'
#' @param model The network model with which to generate the network; \code{"BA"} for Barabási–Albert, \code{"ER"} for Erdős–Rényi, or \code{"GEO"} for geometric
#' @param amnt.nodes the number of nodes in the network at any given type
#' @param amnt.edges the number of edges in the network at any given type
#' @param amnt.operations the number of edge additions/deletions to generate
#' @param trace will print output text if \code{TRUE}
#' @param amnt.dimensions (only GEO) the number of dimensions in which to operate
#' @param offset.exponent (only BA) the offset exponent for the weighted sampling
#' @param ... extra parameters to pass to a specific network generator
#'
#' @return A list containing the starting network \code{network} and the dynamic operations peformed on it \code{operations}.
#' @export
#'
#' @rdname generate.dynamic.network
#'
#' @examples
#' # dyn.net.ba <- generate.dynamic.network("BA", 300, 300, 1000)
#' dyn.net.er <- generate.dynamic.network("ER", 300, 300, 1000)
#' dyn.net.geo <- generate.dynamic.network("GEO", 300, 300, 1000)
generate.dynamic.network <- function(model, amnt.nodes, amnt.edges, amnt.operations, trace = T, ...) {
  generate.dynnet <- switch(
    model,
    GEO = generate.geometric,
    BA = generate.barabasialbert,
    ER = generate.erdosrenyi
  )
  generate.dynnet(amnt.nodes, amnt.edges, amnt.operations, trace = trace, ...)
}

#' @rdname generate.dynamic.network
#' @importFrom dplyr bind_rows
#' @importFrom utils head
#' @importFrom stats runif quantile
#' @export
generate.geometric <- function(amnt.nodes, amnt.edges, amnt.operations, amnt.dimensions = 3, trace = T) {
  if (amnt.nodes > 24000) stop("This function currently does not support amnt.nodes > 24000")
  initialisation.t <- system.time({
    locations <- matrix(stats::runif(amnt.nodes*amnt.dimensions), ncol=amnt.dimensions)

    idx <- matrix(sample.int(nrow(locations), 2 * 1000000, replace = T), ncol = 2)
    idx <- idx[idx[,1] != idx[,2], ]
    vec <- sqrt(rowSums((locations[idx[,1],]-locations[idx[,2],])^2))

    about.cutoff.dist <- stats::quantile(vec, amnt.edges / (amnt.nodes * (amnt.nodes-1)))
    names(about.cutoff.dist) <- NULL
    num.bins <- floor(1/about.cutoff.dist) # lower 1 to .75 or .5 ifneedbe
    bins <- ceiling(locations * num.bins)-1
    row.mult <- num.bins^(seq_len(amnt.dimensions)-1)
    bin.ids <- apply(bins, 1, function(z) sum(z * row.mult))

    bin.df <- do.call(expand.grid, lapply(seq_len(amnt.dimensions), function(x) seq_len(num.bins)-1))
    bin.df$id <- apply(bin.df, 1, function(z) sum(z * row.mult))
    num.different.bins <- nrow(bin.df)
    check.bins <-
      lapply(seq_len(num.different.bins)-1, function(i) {
        neighbour.bins <- unlist(sapply(seq(0, amnt.dimensions-1), function(z) c(i + num.bins^z, i - num.bins^z), simplify = F))
        neighbour.bins <- neighbour.bins[0 <= neighbour.bins & neighbour.bins < num.different.bins]
        c(i, neighbour.bins)
      })

    binned.idx <- lapply(seq(0, num.different.bins-1), function(b) which(bin.ids == b))

    dists.between.bins <- lapply(seq_along(binned.idx), function(bin.i) {
      df <- expand.grid(
        i = binned.idx[[bin.i]],
        j = unlist(binned.idx[check.bins[[bin.i]]+1])
      )
      df <- df[df$i < df$j,,drop=F]
      df$dist <- sqrt(rowSums((locations[df$i,,drop=F] - locations[df$j,,drop=F])^2))
      df
    })

    begin.df <- dplyr::bind_rows(dists.between.bins)
    ord <- utils::head(order(begin.df$dist), amnt.edges)
    begin.df <- begin.df[ord,1:2,drop=F]
    current.df <- begin.df
    operations.df <- data.frame(type = factor(NULL, levels = c("ADD", "REM")), i = numeric(), j = numeric())
  })

  operations.t <- system.time({
    while (nrow(operations.df) < amnt.operations) {
      if (trace) cat("Ops progress: ", nrow(operations.df), "/", amnt.operations, "\n", sep="")
      node <- ceiling(runif(1) * amnt.nodes)
      locations[node,] <- runif(amnt.dimensions)

      prev.bin <- bin.ids[[node]]
      bins[node,] <- ceiling(locations[node,] * num.bins)-1
      new.bin <- sum(bins[node,] * row.mult)
      bin.ids[[node]] <- new.bin

      rem.bins <- check.bins[[prev.bin+1]]+1
      add.bins <- check.bins[[new.bin+1]]+1

      recalculate.bins <- unique(c(add.bins, rem.bins))
      dists.between.bins[recalculate.bins+1] <- lapply(recalculate.bins, function(bin.i) {
        df <- expand.grid(
          i = binned.idx[[bin.i]],
          j = unlist(binned.idx[check.bins[[bin.i]]+1])
        )
        df <- df[df$i < df$j,,drop=F]
        df$dist <- sqrt(rowSums((locations[df$i,,drop=F] - locations[df$j,,drop=F])^2))
        df
      })

      new.df <- dplyr::bind_rows(dists.between.bins)
      ord <- utils::head(order(new.df$dist), amnt.edges)
      new.df <- new.df[ord,1:2,drop=F]

      tmp.new.df <- new.df
      tmp.new.df$ADD <- TRUE
      tmp.cur.df <- current.df
      tmp.cur.df$REM <- TRUE

      new.operations <- dplyr::bind_rows(tmp.new.df, tmp.cur.df)
      new.operations <- new.operations[is.na(new.operations$ADD) != is.na(new.operations$REM),,drop=F]
      new.operations$type <- factor(ifelse(is.na(new.operations$ADD), "REM", "ADD"), levels = c("ADD", "REM"))
      new.operations <- new.operations[sample.int(nrow(new.operations)),c("i", "j", "type")]

      operations.df <- dplyr::bind_rows(operations.df, new.operations)
      current.df <- new.df
    }

    operations.df <- utils::head(operations.df, amnt.operations)
  })

  timings <- data.frame(initialisation = initialisation.t[["elapsed"]], operations = operations.t[["elapsed"]])

  list(network = begin.df, operations = operations.df, timings = timings)
}

#' @rdname generate.dynamic.network
#' @export
generate.barabasialbert <- function(amnt.nodes, amnt.edges, amnt.operations, offset.exponent=1, trace=T) {
  if (round(amnt.edges / amnt.nodes) != amnt.edges / amnt.nodes) stop("the amount of edges needs to be a multiple of the amount of nodes.")
  initialisation.t <-
    system.time({
      m <- amnt.edges / amnt.nodes

      init.net <- function() {
        degree <- rep(0, amnt.nodes)
        neighbours <- lapply(seq_len(amnt.nodes), function(i) numeric(0))
        back.neighbours <- lapply(seq_len(amnt.nodes), function(i) numeric(0))
        operations <- data.frame(type=character(), i=numeric(), j=numeric())
        list(degree=degree, neighbours=neighbours, back.neighbours=back.neighbours, operations=operations)
      }

      add.edge <- function(net, i, j) {
        net$degree[c(i,j)] <- net$degree[c(i,j)] + 1
        net$neighbours[[i]] <- c(net$neighbours[[i]], j)
        net$back.neighbours[[j]] <- c(net$back.neighbours[[j]], i)
        net$operations <- rbind(net$operations, data.frame(op="ADD", i=i, j=j))
        net
      }

      remove.edge <- function(net, i, j) {
        net$degree[c(i,j)] <- net$degree[c(i,j)] - 1
        n <- net$neighbours[[i]]
        n <- n[n != j]
        net$neighbours[[i]] <- n
        m <- net$back.neighbours[[j]]
        m <- m[m != i]
        net$back.neighbours[[j]] <- m
        net$operations <- rbind(net$operations, data.frame(op="REM", i=i, j=j))
        net
      }

      remove.node <- function(net, i) {
        n <- net$neighbours[[i]]
        for (j in n) {
          net <- remove.edge(net, i, j)
        }
        net
      }

      reset.operations <- function(net) {
        net$operations <- data.frame(type=character(), i=numeric(), j=numeric())
        net
      }

      add.node.with.random.edges <- function(net, i) {
        bn <- net$back.neighbours[[i]]
        poss.neighs <- seq_len(i-1) # only smaller than i
        if (length(bn) > 0) {
          poss.neighs <- poss.neighs[-bn]
        }
        w <- net$degree[poss.neighs]
        w <- (w / sum(w))
        w <- w ^ offset.exponent
        neighs <- sample(poss.neighs, size=m, replace=F, prob=w)
        for (j in neighs) {
          net <- add.edge(net, i, j)
        }
        net
      }

      net <- init.net()

      init.m <- m#max(m, 5)
      # initialise network with m nodes
      for (i in seq_len(init.m)+1) {
        for (j in seq_len(i-1)) {
          net <- add.edge(net, i, j)
        }
      }

      # add nodes
      for (i in (init.m+2):amnt.nodes) {
        net <- add.node.with.random.edges(net, i)
      }

      final.network <- bind_rows(lapply(seq_len(amnt.nodes), function(i) data.frame(i=rep(i, length(net$neighbours[[i]])), j=net$neighbours[[i]])))
      net <- reset.operations(net)
    })
  operations.t <-
    system.time({
      while (nrow(net$operations) < amnt.operations) {
        if (trace)  cat("Ops progress: ", nrow(net$operations), "/", amnt.operations, "\n", sep="")
        i <- sample((init.m+2):amnt.nodes, 1)
        net <- remove.node(net, i)
        net <- add.node.with.random.edges(net, i)
      }

      net$operations <- net$operations[seq_len(amnt.operations),]
    })

  timings<-data.frame(initialisation=initialisation.t[["elapsed"]], operations=operations.t[["elapsed"]])

  list(network=final.network, operations=net$operations, timings=timings)
}

#' @rdname generate.dynamic.network
#' @importFrom dplyr bind_rows
#' @export
generate.erdosrenyi <- function(amnt.nodes, amnt.edges, amnt.operations, trace = T) {
  if (amnt.edges > amnt.nodes * (amnt.nodes - 1) / 2)
    stop("amnt.edges should be less than or equal to amnt.nodes * (amnt.nodes-1) / 2")

  initialisation.t <-
    system.time({
      v <- sapply(1:(amnt.nodes-1), function(i) i * (i+1) / 2)
      max.x <- max(v)
      xs <- sample(1:max.x, amnt.edges, replace=F)

      edges <- dplyr::bind_rows(lapply(xs, function(x) {
        i <- min(which(v >= x))+1
        j <- x - ifelse(i > 2, v[[i-2]], 0)
        data.frame(i=i, j=j)
      }))

      original.network <- edges
    })

  operations.t <-
    system.time({
      operations <- data.frame(type=character(), i=numeric(), j=numeric())
      while (nrow(operations) < amnt.operations) {
        rem.y <- sample(1:length(xs), 1)
        operations <- rbind(operations, data.frame(op="REM", edges[rem.y,]))
        edges <- edges[-rem.y,]
        xs <- xs[-rem.y]

        add.y <- sample((1:max.x)[-xs], 1)
        i <- min(which(v >= add.y))+1
        j <- add.y - ifelse(i > 2, v[[i-2]], 0)
        edge <- data.frame(i=i, j=j)
        edges <- rbind(edges, edge)
        xs <- c(xs, add.y)
        operations <- rbind(operations, data.frame(op="ADD", edge))
      }
    })
  timings <- data.frame(initialisation = initialisation.t[["elapsed"]], operations = operations.t[["elapsed"]])

  list(network = original.network, operations = operations, timings = timings)
}
