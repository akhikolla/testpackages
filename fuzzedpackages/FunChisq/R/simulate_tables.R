# simulate_tables : Generates simulated contingency tables (with or without noise) of a given "type".
#
# Parameters details
# n              : Expected sample size, must be at least nrow for "functional", "many.to.one", "discontinuous"
#                  and nrow * ncol for "independent" and "dependent.non.functional" [DEFAULT = 100].
# nrow           : Expected number of rows in the table, must be greater than 2 (greater than 3 for type = "many.to.one") [DEFAULT = 3].
# ncol           : Expected number of columns in the table , must be greater than 2 [DEFAULT = 3].
# types include  : functional: y is a function of x but x may or may not be a function of y [DEFAULT] ,
#                  many.to.one: y is a function of x and x is not a function of y,
#                  dependent.non.functional: Non functional table with statistical dependency,
#                : discontinuous: y is a duscontinuous function of x but x may or may not be function of y,
#                  independent: Independent tables (null population).
# n.tables       : number of tables to be generated.[DEFAULT = 1].
# row.marginal   : A vector of row probabilities [DEFAULT = Equally likely].
# col.marginal   : A vector of column probabilities (Only required for type = "independent") [DEFAULT = Equally likely].

# Created by     : Ruby Sharma
# Date           : October 16 2016
#
# Modified       : Ruby Sharma. April 27 2017
# Version        : 0.0.2
# Updates        : Added pattern type "discontinuous"
#                : Added noise.model parameter
#
# Modified       : Ruby Sharma December 2 2018
# Version        : 0.0.3
# Updates        : Allowed small sample size for independent tables
#                : Modified the coding for independent tables

# Modified       : MS December 4, 2018
# Version        : 0.0.4
# Updates        : Introduced simulate_independent_tables() from the
#                : preivous code, where noise is applied along both
#                : rows and columns (previously only along the rows).
#                : Modified prelim.check()

# Modified       : Ruby Sharma April 20, 2020
# Version        : 0.0.5
# Updates        : parameter margin has been added to apply noise
#                : along the rows or columns or both ways.


# Modified       : Ruby Sharma May 26, 2020
# Version        : 0.0.6
# Updates        : All changes are for functional tables.
#                : Allowed user to provide col.marginal,
#                : (previously only row.marginal)
#                : Name of functional.table is changed to disc.functional.pattern.
#                : Pattern table is not dependent on row and column marginals.
#                : (previously depdendent on row marginals).
#                : If samples are distributed according to column marginal then, samples
#                : within the columns are distributed accroding to row marginals.
#                : Added a new method decide.marginal() : if both marginals are provided,
#                : randomly select one with probability of 0.5 or none of them provided.
#                : Default row and column marginals are NULL.
#                : Added a new method sample.in.cols() to distribute samples according to
#                : each columns accrording to given column marginals.
#                : Added a new method dis.samples() to distribute samples within a column
#                : accroding to row marginals.

simulate_tables <- function(n = 100,
                            nrow = 3,
                            ncol = 3,
                            type = c(
                              "functional",
                              "many.to.one",
                              "discontinuous",
                              "independent",
                              "dependent.non.functional"
                            ),
                            n.tables = 1,
                            row.marginal = NULL,
                            col.marginal = NULL,
                            noise = 0.0,
                            noise.model = c("house", "candle"),
                            margin = 0)
{
  type <- match.arg(type)
  noise.model <- match.arg(noise.model)

  prelim.check(nrow, ncol, n, row.marginal, col.marginal,
               n.tables)

  if (ncol < 2 || nrow < 2)
    stop("ERROR: Numbers of rows and columns must be >= 2!\n")

  if (type != "independent")
  {
    #Intialization
    pattern.list <- list()
    sample.list <- list()
    noise.list <- list()
    p.value.list <- list()

    for (i in seq(n.tables))
    {
      alltables = table.generate(nrow,
                                 ncol,
                                 type,
                                 n,
                                 row.marginal,
                                 col.marginal,
                                 noise,
                                 noise.model,
                                 margin)
      pattern.list[[i]] = alltables$pattern.table
      sample.list[[i]] = alltables$sampled.table
      noise.list[[i]] = alltables$noise.table
      p.value.list[[i]] = alltables$p.value
    }

    tbls <-
      list(
        pattern.list = pattern.list,
        sample.list = sample.list,
        noise.list = noise.list,
        pvalue.list = p.value.list
      )

  } else {
    tbls <- simulate_independent_tables(n,
                                        nrow,
                                        ncol,
                                        n.tables,
                                        row.marginal,
                                        col.marginal,
                                        noise,
                                        noise.model,
                                        margin)
  }

  # Return a list of pattern table, sampled contingency table
  #   and noise table:
  return(tbls)
}

simulate_independent_tables <- function(n,
                                        nrow,
                                        ncol,
                                        n.tables,
                                        p.row.marginal,
                                        p.col.marginal,
                                        noise,
                                        noise.model,
                                        margin)
{
  if (is.null(p.row.marginal))
  {
    p.row.marginal = rep(1 / nrow, nrow)
  }

  if (is.null(p.col.marginal))
  {
    p.col.marginal = rep(1 / ncol, ncol)
  }

  if (n < 0 && (n %% 1 != 0)) {
    stop("ERROR: n must be a positive integer!\n")
  }

  prob.table <- p.row.marginal %*% t(p.col.marginal)

  pattern.table <- matrix(0, nrow = nrow, ncol = ncol)
  pattern.table[prob.table > 0] <- 1

  pattern.list <- lapply(seq(n.tables), function(i) {
    return(pattern.table)
  })

  sampled.V <- rmultinom(n.tables, n, prob.table)
  sample.list <- lapply(as.data.frame(sampled.V), function(v) {
    mat <- matrix(v,
                  nrow = nrow,
                  ncol = ncol,
                  byrow = FALSE)
    return(mat)
  })


  # apply noise along both row and columns
  noise.list <- add.noise(
    tables = sample.list,
    u = noise,
    noise.model = noise.model,
    margin = margin
  )

  p.value.list = lapply(1:n.tables, function(k) {
    return(chisq.test.pval(sample.list[[k]]))
  })

  return(
    list(
      pattern.list = pattern.list,
      sample.list = sample.list,
      noise.list = noise.list,
      pvalue.list = p.value.list
    )
  )
}

#generating pattern table, sampled contingency table and noise table
table.generate = function(nrow,
                          ncol,
                          type,
                          n,
                          row.marginal,
                          col.marginal,
                          noise,
                          noise.model,
                          margin)
{
  allmar = decide.marginal(nrow, ncol, row.marginal, col.marginal)
  row.marginal = allmar$Xmarg
  col.marginal = allmar$Ymarg
  mar.type = allmar$Martype

  if (type == "dependent.non.functional") {
     if (n < (nrow * ncol))
        stop(paste(
        "For dependent.non.functional, n must be greater than or equal to",
        (nrow * ncol)
      ))
      col.marginal = rep(1/ncol, ncol)
      sam.val.row = sample.in.rows(n, row.marginal, type, ncol)
      pattern.table = nonfunctional.table(nrow, ncol, row.marginal, sam.val.row)
      prob.table = non.functional.prob(nrow, ncol, pattern.table)
      sampled.table = dis.sample.prob(nrow, ncol, sam.val.row, prob.table, col.marginal)

      sampled.table = is_dependent(n,
                                   nrow,
                                   ncol,
                                   row.marginal,
                                   col.marginal,
                                   sampled.table,
                                   sam.val.row)

      p.val = chisq.test.pval(sampled.table)

  } else if (type == "many.to.one") {
      if (nrow < 3)
        stop("For many.to.one, number of rows must be at least be 3!")

      if (length(which(row.marginal != 0)) < 3)
        stop("For many.to.one, at least three non-zero row probabilities are expected!")

      if (length(which(col.marginal != 0)) < 2)
        stop("For many.to.one, at least two non-zero column probabilities are expected!")

      if (n < nrow)
        stop(paste("n must be greater than or equal to", nrow))

      pattern.table = many.to.one.table(nrow, ncol)

      if (mar.type == "col") {
        sam.val.col = sample.in.cols(n, pattern.table, col.marginal,  ncol)
        sampled.table = dis.samples(n, pattern.table, sam.val.col, row.marginal, nrow, ncol)
      } else {
        sam.val.row = sample.in.rows(n, row.marginal, type, ncol)
        sampled.table = dis.sample.prob(nrow, ncol, sam.val.row, pattern.table, col.marginal)
      }

    p.val = FunChisq::fun.chisq.test(sampled.table)$p.value

  } else if (type == "discontinuous") {
      if (n < nrow)
        stop(paste("n must be greater than or equal to", nrow))

      pattern.table = discontinuous.table(nrow, ncol)
      if (mar.type == "col") {
        sam.val.col = sample.in.cols(n, pattern.table, col.marginal,  ncol)
        sampled.table = dis.samples(n, pattern.table, sam.val.col, row.marginal, nrow, ncol)
      } else {
         sam.val.row = sample.in.rows(n, row.marginal, type, ncol)
         sampled.table = dis.sample.prob(nrow, ncol, sam.val.row, pattern.table, col.marginal)
      }

      p.val = FunChisq::fun.chisq.test(sampled.table)$p.value

  } else {
      if (n < nrow)
        stop(paste("n must be greater than or equal to", nrow))

      pattern.table = disc.functional.pattern(nrow, ncol)

      if (mar.type == "col") {
        sam.val.col = sample.in.cols(n, pattern.table, col.marginal,  ncol)
        sampled.table = dis.samples(n, pattern.table, sam.val.col, row.marginal, nrow, ncol)

      } else {
        sam.val.row = sample.in.rows(n, row.marginal, type, ncol)
        sampled.table = dis.sample.prob(nrow, ncol, sam.val.row, pattern.table, col.marginal)
      }

      p.val = FunChisq::fun.chisq.test(sampled.table)$p.value

  }

  noise.table = add.noise(
    tables = sampled.table,
    u = noise,
    noise.model = noise.model,
    margin = margin
  )

  # return singular tables
  list(
    pattern.table = pattern.table,
    sampled.table = sampled.table,
    noise.table = noise.table,
    p.value = p.val
  )
}

# distributing samples across rows guided by row probabilities
sample.in.rows = function(n, row.marginal, type, ncol)
{
  non.zero.rows.length = length(which(row.marginal != 0))
  if (type != "independent") {
    n = n - non.zero.rows.length
    sam.val = rmultinom(1, n, row.marginal)
    ind = which(row.marginal != 0)
    sam.val[ind] = sam.val[ind] + 1
  } else {
    n = n - (non.zero.rows.length * ncol)
    sam.val = rmultinom(1, n, row.marginal)
    ind = which(row.marginal != 0)
    sam.val[ind] = sam.val[ind] + ncol
  }
  return(sam.val)
}

sample.in.cols = function(n, pattern.table, col.marginal,  ncol)
{
  zeroCol = which(colSums(pattern.table) == 0)
  col.marginal[zeroCol] = 0
  colSamp = rmultinom(1, n, col.marginal)
  return(colSamp)
}

# Distribute samples according to column distribution
dis.samples = function(n,
                       pattern.table,
                       sam.val.col,
                       row.marginal,
                       nrow,
                       ncol)
{
  sampled.table = matrix(0, nrow = nrow, ncol = ncol)

  for (i in 1:ncol)
  {
    if (sam.val.col[i] != 0) {
      non.zero.rows = which(pattern.table[, i] == 1)
      prob = row.marginal[non.zero.rows]
      sampled.table[non.zero.rows, i] = rmultinom(1, sam.val.col[i], prob)
    }
  }
  return(sampled.table)
}


# distribute samples according to row marginals
dis.sample.prob = function(nrow, ncol, sam.val.row, table, col.marginal)
{
  sampled.table = matrix(0, nrow = nrow, ncol = ncol)
  for (i in 1:nrow)
  {
    # determine if all columns in the ith row of the supplied table are zero
    all.zero.column = all(table[i, ] == 0)
    if (!all.zero.column) {
      non.zero.columns = which(table[i, ] != 0)
      sampled.table[i,non.zero.columns] = rmultinom(1, sam.val.row[i], col.marginal[non.zero.columns])

      #size = sam.val.row[i] - length(non.zero.columns)
      #sam.val.cell = rmultinom(1, size, table[i, ])
      #sam.val.cell[non.zero.columns] = sam.val.cell[non.zero.columns] + 1
      #sampled.table[i, ] = sam.val.cell
    }
  }
  return(sampled.table)
}

# generating pattern table for "functional"
disc.functional.pattern = function(nrow, ncol)
{
  pattern.table = matrix(0, ncol = ncol, nrow = nrow)

  for (i in 1:nrow)
  {
    index =  sample(1:ncol, 1)
    pattern.table[i, index] = 1

  }
  # check for constant functions
  pattern.table = not_constant(ncol, pattern.table)
  return(pattern.table)
}

# generating pattern table for "many.to.one"
many.to.one.table = function(nrow, ncol)
{
  pattern.table = matrix(0, ncol = ncol, nrow = nrow)
  # get the functional table
  pattern.table = disc.functional.pattern(nrow, ncol)
  # check for non-monotonicity
  pattern.table = is_many.to.one(pattern.table)

  return(pattern.table)
}

# generating pattern table for "dependent.non.functional"
nonfunctional.table = function(nrow, ncol, row.marginal, sam.val.row)
{
  pattern.table = matrix(0, ncol = ncol, nrow = nrow)

  for (i in 1:nrow)
  {
    if (sam.val.row[i] != 0) {
      if (sam.val.row[i] < ncol) {
        index = sample(1:ncol, sample(1:sam.val.row[i], 1))
      } else {
        index = sample(1:ncol, sample(1:ncol, 1))
      }

      pattern.table[i, index] = 1
    }
  }
  pattern.table = make.non.functional(pattern.table, ncol, sam.val.row, nonfunc = "notf.x")
  pattern.table = make.non.functional(t(pattern.table), ncol, sam.val.row, nonfunc = "notf.y")
  pattern.table = t(pattern.table)
  return(pattern.table)
}


# generating pattern table for "discontinuous function"
discontinuous.table = function(nrow, ncol)
{
  pattern.table = matrix(0, ncol = ncol, nrow = nrow)

  sample.from = seq(ncol(pattern.table))

  for (i in 1:nrow)
  {
    if (length(sample.from) == 1) {
      index = sample.from
    } else {
      index = sample(sample.from, 1)
    }
    sample.from = seq(ncol(pattern.table))
    prev.col.ind = index
    pattern.table[i, index] = 1
    sample.from = sample.from[-which(sample.from == prev.col.ind)]
 }

  return(pattern.table)
}


# introducting atleast two entries in one row for nonfunctional table
make.non.functional = function(pattern.table, ncol, sam.val.row, nonfunc)
{
  indexes = non.zero.index(pattern.table)
  rows = indexes$rows
  cols = indexes$cols


  if (nonfunc == "notf.x") {
    if (anyDuplicated(rows) == 0) {
      only.row = which(sam.val.row > 1)
      chng.from = c(only.row)
      if (length(chng.from) == 1) {
        chng.row.index = chng.from
      } else {
        chng.row.index = sample(chng.from, 1)
      }
      prev.col.index = cols[which(rows == chng.row.index)]
      chng.col.index = sample.sec.col.ind(prev.col.index, ncol)
      pattern.table[chng.row.index, chng.col.index] = 1
    }
  }
  if (nonfunc == "notf.y") {
    nr = nrow(pattern.table)
    nc = ncol(pattern.table)
    wth.more.samp <- c()
    wth.less.samp <- c()
    for (i in 1:nc)
    {
      sample.lim <- length(cols[cols[] == i])
      if (sam.val.row[i] > sample.lim) {
        wth.more.samp = c(i, wth.more.samp)
      } else {
        wth.less.samp = c(i, wth.less.samp)
      }

    }
    if (anyDuplicated(rows) == 0) {
      chng.from = rows
      if (length(chng.from) == 1) {
        chng.row.index = chng.from
      } else {
        chng.row.index = sample(chng.from, 1)
      }

      prev.col.index = cols[which(rows == chng.row.index)]
      chng.col.index = sample.sec.col.ind_notf.y(prev.col.index, nc, wth.more.samp)
      if (chng.col.index == 0) {
        if (length(rows[cols[] == wth.less.samp]) == 1) {
          chng.row.index = rows[cols[] == wth.less.samp]
        } else {
          chng.row.index = sample(rows[cols[] == wth.less.samp], 1)
        }
        chng.col.index = wth.more.samp
      }
      pattern.table[chng.row.index, chng.col.index] = 1
    }
  }

  return(pattern.table)
}


# generating probability table for "independent"
indep.prob.table = function(nrow, ncol, row.marginal, col.marginal)
{
  prob.table = matrix(0, ncol = ncol, nrow = nrow)

  # multiplying row probability and column probability in in each cell of prob.table
  prob.table = row.marginal %*% t(col.marginal)
  return(prob.table)
}

# generating probability table for "dependent.non.functional"
non.functional.prob = function(nrow, ncol, pattern.table)
{
  prob.table = matrix(0, ncol = ncol, nrow = nrow)

  indexes = non.zero.index(pattern.table)
  rows = indexes$rows
  cols = indexes$cols

  for (i in 1:nrow) {
    if (is.element(i, rows)) {
      col.ind = which(rows == i)
      freq.no = length(col.ind)
      prob.ele.row = 1 / freq.no
      prob.table[i, cols[col.ind]] = prob.ele.row
    }
  }
  return(prob.table)
}


# check for constant function, if found then change one column index
not_constant = function(ncol, pattern.table)
{
  indexes = non.zero.index(pattern.table)
  rows = indexes$rows
  cols = indexes$cols

  if (length(unique(cols)) == 1) {
    index = cols[1]
    chng.col.index = sample.sec.col.ind(index, ncol)
    chng.row.index = sample(rows, 1)
    pattern.table[chng.row.index, index] = 0
    pattern.table[chng.row.index, chng.col.index] = 1
  }

  return(pattern.table)
}


# check monotonicity, if found then atleast two rows would share samples in the same column
is_many.to.one = function(pattern.table)
{
  indexes = non.zero.index(pattern.table)
  rows = indexes$rows
  cols = indexes$cols

  if (length(unique(cols)) == length(rows)) {
    id = sample(rows, 2)
    pattern.table[id[2], cols[which(rows == id[2])]] = 0
    pattern.table[id[2], cols[which(rows == id[1])]] = 1
  }

  return(pattern.table)
}

#check whether the dependent.non.functional table is dependent, if not make dependent
is_dependent = function(n,
                        nrow,
                        ncol,
                        row.marginal,
                        col.marginal,
                        sampled.table,
                        sam.val.row)
{
  #sam.val.indep =  sample.in.rows(n, row.marginal, type = "independent", ncol)
  #expec.prob.table = matrix(0, ncol = ncol, nrow = nrow)
  #expec.prob.table = indep.prob.table(nrow, ncol, row.marginal, col.marginal)
  #indep.sampled.table = dis.sample.prob(nrow, ncol, sam.val.indep, expec.prob.table, col.marginal)
  indep.sampled.table = simulate_independent_tables(n, nrow, ncol, n.tables = 1, p.row.marginal = row.marginal,
                              p.col.marginal = col.marginal,
                              noise.model = "house", noise = 0.0, margin = 0)$sample.list[[1]]


  difference = indep.sampled.table - sampled.table
  if (all(difference == 0)) {
    indexes = non.zero.index(sampled.table)
    rows = indexes$rows
    cols = indexes$cols
    sel.row <- sample(rows, 1)
    sel.from <- cols[which(rows == sel.row)]

    if (length(sel.from) == 1) {
      sel.colj1 = sel.from
    } else {
      sel.colj1 <- sample(sel.from, 1)
    }

    sel.from <- 1:ncol
    sel.colj2 <- sample(sel.from[-which(sel.from == sel.colj1)], 1)
    sampled.table[sel.row, sel.colj2] <-
      sampled.table[sel.row, sel.colj2] + sampled.table[sel.row, sel.colj1]
    sampled.table[sel.row, sel.colj1] <- 0
  }

  return(sampled.table)

}


# sorting row and column indexes on the basis of row
sort.index = function(rows, cols)
{
  row.col.ind = matrix(nrow = length(rows), ncol = 2)
  row.col.ind[, 1] = rows
  row.col.ind[, 2] = cols
  row.col.ind = row.col.ind[order(row.col.ind[, 1], row.col.ind[, 2]), ]
  rows = row.col.ind[, 1]
  cols = row.col.ind[, 2]

  list(row.in = rows, col.in = cols)
}


# retrieving non zero indexes from the table
non.zero.index = function(table)
{
  rows = row(table)[which(!table == 0)]
  cols = col(table)[which(!table == 0)]
  sorted = sort.index(rows, cols)
  rows = sorted$row.in
  cols = sorted$col.in

  list(rows = rows, cols = cols)
}

# sampling secondary column index
sample.sec.col.ind = function(index, ncol)
{
  if (index == 1)
    vec.to.sample = c(2:ncol)

  if (index == ncol)
    vec.to.sample = c(1:(ncol - 1))

  if (index > 1 && index < ncol)
    vec.to.sample = c(1:(index - 1), (index + 1):ncol)

  if (length(vec.to.sample) == 1) {
    chng.col.index = vec.to.sample
  } else {
    chng.col.index = sample(vec.to.sample, 1)
  }

  return(chng.col.index)
}

# Parameters checking
prelim.check <- function(nrow,
                         ncol,
                         n,
                         row.marginal,
                         col.marginal,
                         n.tables)
{
  if (class(nrow) != "numeric" && class(nrow) != "integer")
    stop("ERROR: nrow must be numeric!\n")

  if (class(ncol) != "numeric" && class(ncol) != "integer")
    stop("ERROR: ncol must be numeric!\n")

  if (class(n) != "numeric" && class(n) != "integer")
    stop("ERROR: n must be numeric!\n")

  if (!is.null(row.marginal)) {
    if (class(row.marginal) != "numeric" &&
        class(row.marginal) != "integer")
      stop("ERROR: row.marginal must be numeric!\n")

    if (any(row.marginal < 0))
      stop("ERROR: row.marginal must be >= 0!\n")


    if (length(row.marginal) < nrow)
      stop("ERROR: Row marginal probabilites for all rows expected!\n")

    if (sum(row.marginal != 0) < 2)
      stop("ERROR: Two or more non-zero row probabilites expected!\n")

  }
  if (!is.null(col.marginal)) {
    if (class(col.marginal) != "numeric" &&
        class(col.marginal) != "integer")
      stop("ERROR: col.marginal must be numeric!\n")

    if (any(col.marginal < 0))
      stop("ERROR: col.marginal must be >= 0!\n")

    if (length(col.marginal) < ncol)
      stop("ERROR: Column marginal probabilites for all columns expected!\n")

    if (sum(col.marginal != 0) < 2)
      stop("ERROR: Two or more non-zero column probabilites expected!\n")

  }

  if (n.tables <= 0 || (n.tables %% 1 != 0))
    stop("ERROR: n.tables must be a positive integer!\n")

}

chisq.test.pval <- function(table)
{
  if (length(table[table != 0]) == 1) {
    pval = NA
  } else {
      # identify non-zero rows and non-zero columns
      non.zero.rows <-
        apply(table, 1, function(row) {
        0 != sum(abs(row))
      })
      non.zero.cols <-
        apply(table, 2, function(col) {
        0 != sum(abs(col))
      })

    # perform Pearson chi-square test
    chisq <- suppressWarnings(chisq.test(table[non.zero.rows, non.zero.cols])$statistic)

    # compute p-value using the orgional table size
    pval <- pchisq(chisq, prod(dim(table) - 1), lower.tail = FALSE)
  }

  names(pval) <- NULL
  return(pval)
}

#sampling second column for making x!=f(y) for dependent.non.functional
sample.sec.col.ind_notf.y <- function(index , ncol, only.col)
{
  if (index == 1)
    vec.to.sample = c(2:ncol)

  if (index == ncol)
    vec.to.sample = c(1:(ncol - 1))

  if (index > 1 && index < ncol)
    vec.to.sample = c(1:(index - 1), (index + 1):ncol)

  vec.to.sample =  vec.to.sample[vec.to.sample %in% only.col]

  if (length(vec.to.sample) == 0) {
    chng.col.index = 0
  } else if (length(vec.to.sample) == 1) {
    chng.col.index = vec.to.sample
  } else {
    chng.col.index = sample(vec.to.sample, 1)
  }
  return(chng.col.index)
}

# Decide the marginal type depending on the user input
# If both marginals are null, assign it uniform probabilities
# and sample one of them for creating the desired table
# If both marginals are provided, select one to create the table
# If one of the marginal is provided, use that one to create the table.
decide.marginal <- function(nrow, ncol, row.marginal, col.marginal)
{
  if (is.null(row.marginal) && is.null(col.marginal)) {
    row.marginal = rep(1 / nrow, nrow)
    col.marginal = rep(1 / ncol, ncol)
  }

  if (!is.null(row.marginal) && !is.null(col.marginal)) {
    martype = sample(c("row", "col"), 1)
  }

  if (is.null(row.marginal) && !is.null(col.marginal)) {
    row.marginal = rep(1 / nrow, nrow)
    martype = "col"
  } else if (!is.null(row.marginal) && is.null(col.marginal)) {
    col.marginal = rep(1 / ncol, ncol)
    martype = "row"
  }

  list(Xmarg = row.marginal,
       Ymarg = col.marginal,
       Martype = martype)
}
