# Functions for converting ape-format phyogenies to
# 'peth': a format better suited for time-forward simulation.

## Convert a tree in ape format to one in new format (return a new 'peth' tree object). Only works for ntip > 4

ape2peth <- function(tree) {
    edge <- tree$edge
    edge.length <- tree$edge.length
    nedge <- length(edge[,1])
    ntip <- Ntip(tree)

    # Set up ape-to-peth map and vector to keep track of running branches.
    running <- rep(FALSE, 2 * ntip - 2)
    map <- matrix(ncol=2, nrow=nedge + 1)

    # Assign root
    map[1,1] <- root <- edge[1,1]
    map[1,2] <- 1

    # splitting gives the NODE, as labelled in edge, which is speciating
    splitting <- root
	edge_counter <- 1

    # Keep track of how many nodes we've already assigned.
    node_counter <- 1

    # Produce vector of nodes to split for future simulations (peth numbering) and Produce vector of times between speciation events
    nts <- tts <- c()

    # Loop over speciation events
    for(i in seq(2, 2*(ntip-1), by=2)) {
        
        # 'splitting' is the splitting node in ape format
        # 'peth_split' is the splitting node in peth format
        peth_split <- map[ which(map[,1] == splitting), 2 ][1]
        nts <- c(nts, peth_split)

        # The branches coming from this node are now running
        running[ which(edge[,1]==splitting) ] <- TRUE
        running[ which(edge[,2]==splitting) ] <- FALSE

        # We're at a speciation event. Find out which nodes come from it.
        new_edges <- which(edge[,1] == splitting)
        new_nodes <- edge[new_edges, 2]

        # Get which of the new nodes comes from the shortest branch
        # We want this one to inherit the label-number of its parent.
        new_node1 <- new_nodes[order(edge.length[new_edges])[1]]
        new_node2 <- new_nodes[order(edge.length[new_edges])[2]]

        # First new node keeps ancestor number (splitting node), 
        # other new node takes node_counter+1
        map[i,1] <- new_node1
        map[i,2] <- peth_split
        map[i+1,1] <- new_node2
        map[i+1,2] <- node_counter + 1

        node_counter <- node_counter + 1

        # time_to_split = time to next split = shortest branch running.
        time_to_split <- min(edge.length[running])
        tts <- c(tts, time_to_split)

        # Subtract time to next split from running branches.
        edge.length[running] <- edge.length[running] - time_to_split

        # Which node splits next? 
        # Of the running branches, take the shortest one and then the splitting node 
        # is the second elemtn of edge for that branch (i.e. the child)
        el_temp <- edge.length
        el_temp[!running] <- 9e99
        splitting <- edge[ which.min(el_temp),2 ]
    }
    
    # Need to know what order simulated data will be in, relative to ape's expectatons
    ordered_map <- order(map[,1])
    want <- ordered_map[1:ntip]
    dord <- rev(map[want,2])

    # Need to know what order simulated data will be in, relative to ape's expectatons. This is given by working from bottom up through map
    dord <-  1:ntip
    for(i in 1:ntip) {
        want <- which(map[,1]==i)
        dord[i] <- map[want , 2]
    }

    ptree <- list(map, dord, nts, tts)
    class(ptree) = 'pethtree'
    names(ptree) = c('map', 'data_order', 'splitting_nodes', 'times')
    return(ptree)
}

## Convert ape to peth, or pass on peth tree, or return error if incorrect

checktree <- function(phy) {
    if(class(phy) == "phylo") {
        phy = ape2peth(phy)
    } else if(class(phy) != "pethtree") {
        stop("Tree incorrectly formatted.")
    }
    phy
}

# Convert time-matrices of when species become sympatric/allopatric to vectors
# If matrix is NA then all elements are set to constant number n

vectortime <- function(sympatric, degree, n.tips) {
  if (all(is.na(sympatric))) {
    sympatric.out <- replicate(n.tips ^ 2, degree)
  } else {
    sympatric.out <- as.vector(sympatric)
  }
  sympatric.out
}

## Reorder a dataset to match the order of tips in an ape-format tree

reorder_data <- function(tree, y, ntraits) {
    data.re.ord.out <- matrix(ncol=ntraits, nrow=length(tree$data_order))
    for (i in 1:ntraits) {
        data.re.ord.out[,i] <- y[seq(i, length(y), by=ntraits)]
        data.re.ord.out[,i] <- data.re.ord.out[ ,i][tree$data_order]
    }
    data.re.ord.out
}

## Generate a matrix of the times at which lineages come into sympatry if each new lineage starts interacting after a delay

symp_matrix <- function(tree, delay=0) {
    
    tree.check = checktree(tree)
    ntip <- length(tree.check$data_order)
    sympatry.out <- matrix(nrow=ntip, ncol=ntip, 0)
    if(class(tree) == "phylo") rownames(s) <- colnames(s) <- tree$tip.label

    # find ages (time from root) of tip lineages
    age <- 1:ntip
    age[1] <- age[2] <- 0
    for(i in 1:ntip) age[i+2] = sum(t$times[1:i])

    # delay measured in terms of mean time between speciation events.
    time_between_splits <- mean(abs(diff(tree.check$times)))
    delay <- delay * time_between_splits

    # apply starttimes and delay to symp-matrix s
    for (i in 1:ntip) 
    {
        for (j in 1:ntip)
        {
            sympatry.out[i,j] = max(c(age[i], age[j])) + delay
        }
    }
    return(sympatry.out)
}

## Get summary statistics for a given tree and dataset

summary_stats <- function(phy, y, est.blomberg.k=FALSE) {
    difs <- as.matrix(dist(y))            # Euclidian distance
    difs[which(difs == 0)] <- NA               # Ignore matrix diagonal
    gap <- apply(difs, 1, min, na.rm=T)  
    y <- as.matrix(y)
    ntraits <- length(y[1,])
    output <- list(mean.gap=mean(gap), sd.gap=sd(gap))

    # Summary statistics: mean and sd of gaps between neighbours. Plus Blomberg's K optionally.
    if(est.blomberg.k) {
        k.est <- sapply(1:ntraits, function(k.int) {
        		y.nit <- as.matrix(y[, k.int])
        		rownames(y.nit) <- phy$tip.label
        		blomberg.k(y=y.nit, phy=phy)
        		}
        	)
        k.est <- mean(k.est)
        output$blom.k <- k.est
    }
    return(output)
  }


   