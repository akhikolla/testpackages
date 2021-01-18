
checkDataSet <- function(X, fracNA = 0.5, numDiscrete = 3, precScale = 1e-12, 
                        silent = FALSE, cleanNAfirst = "automatic") {
  # This function checks the dataset X, and sets aside certain
  # columns and rows that do not satisfy the conditions.
  #
  # fracNA      : only consider columns and rows with fewer NAs than this.
  # numDiscrete : a column that takes on numDiscrete or fewer values
  #               will be considered discrete and not used in the analysis.
  # precScale   : only consider columns whose scale is > precScale.
  #               Here scale is measured by the median absolute deviation.
  # cleanNAfirst: takes "rows", "columns" or "automatic"
  
  
  if (is.null(fracNA)) {
    fracNA <- 0.5
  }
  if (is.null(numDiscrete)) {
    numDiscrete <- 3
  }
  if (is.null(precScale)) {
    precScale <- 1e-12
  }
  if (is.null(silent)) {
    silent <- FALSE
  }
  if (is.null(cleanNAfirst)) {
    cleanNAfirst <- "automatic"
  }
  
  wnq <- function(string, qwrite = 1) { # auxiliary function
    # writes a line without quotes
    if (qwrite == 1) write(noquote(string), file = "", ncolumns = 100)
  }
  
  pnq <- function(string, qwrite = 1) { # auxiliary function
    # prints a line without quotes
    if (qwrite == 1) print(noquote(string))
  } 
  
  if (!is.data.frame(X) & !is.matrix(X)) {
    stop("The input data must be a matrix or a data frame") }    
  
  n <- nrow(X)
  if (n < 3) stop(" The input data must have at least 3 rows (cases)")  
  d <- ncol(X)
  # if (d < 2) stop(" The input data must have at least 2 columns (variables)") 
  
  if (!silent) {
    wnq(" ")
    wnq(paste(" The input data has ", n, " rows and ",
              d, " columns.", sep = ""))
  }
  # Add column names and row names if these were not given:
  if (is.matrix(X)) { X <- data.frame(X) } 
  # This keeps the row names and column names if they exist, else creates
  # them as 1, 2, 3, ... for rows and V1, V2, V3,... for columns.
  
  remX <- X # remX will be the remaining part of the data
  colInAnalysis <- sapply(remX,is.numeric)
  numgoodcol <- sum(colInAnalysis)
  vecNotNumeric <- (colInAnalysis == FALSE)
  numbadcol <- sum(vecNotNumeric) # can be 0
  namesNotNumeric <- NULL  
  if (numbadcol > 0) {
    if (!silent) {
      wnq(" ")
      wnq(paste(" The input data contained ", numbadcol,
                " non-numeric columns (variables).", sep = ""))
      wnq(" Their column names are:")
      wnq(" ")
    }
    namesNotNumeric <- colnames(remX)[vecNotNumeric]
    if (!silent) {
      pnq(namesNotNumeric)   
      wnq(" ")    
    }
    if (numgoodcol > 1) {
      if (!silent) {
        wnq(" These columns will be ignored in the analysis.")
        wnq(paste(" We continue with the remaining ", numgoodcol,
                  " numeric columns:", sep = "")) 
      }
      remX <- remX[colInAnalysis]
    } else { 
      if (numgoodcol == 0) stop(" No columns remain, so we stop.")
      if (numgoodcol == 1) stop(" Only 1 column remains, so we stop.")
    }    
    if (!silent) {
      wnq(" ")
      pnq(names(which(colInAnalysis)))
    }
  }
  # Turn data into a matrix
  remX <- data.matrix(remX, rownames.force = TRUE)
  # This keeps the row names and remaining column names.
  
  # 2. Deselect column(s) containing only the case number
  #
  caseNumber <- (1:nrow(remX))
  distFromCaseNumber <- function(colj, caseNumber) {
    # assumes colj is a vector
    mean(abs(colj - caseNumber), na.rm = FALSE) # can be NA
  }
  dists <- apply(remX, 2, distFromCaseNumber, caseNumber)
  dists[!is.finite(dists)] <- 1 # takes out NA, NaN, Inf, -Inf
  vecbadcol  <- (dists == 0)
  numbadcol  <- sum(vecbadcol) # can be 0
  goodcol    <- (vecbadcol == FALSE)
  numgoodcol <- sum(goodcol)
  namesCaseNumber <- NULL
  if (numbadcol > 0) {
    if (!silent) {
      wnq(" ")
      wnq(paste(" The data contained ", numbadcol, " columns that were",
                " identical to the case number", sep = ""))
      wnq(" (number of the row).")
      wnq(" Their column names are:")
      wnq(" ")
    }
    namesCaseNumber <- colnames(remX)[vecbadcol]
    if (!silent) {
      pnq(namesCaseNumber)
      wnq(" ")
    }
    if (numgoodcol > 1) {
      if (!silent) {
        wnq(" These columns will be ignored in the analysis.")
        wnq(paste(" We continue with the remaining ", numgoodcol,
                  " columns:", sep = "")) 
      }
      remX <- remX[, goodcol, drop = FALSE]
    } else { 
      if (numgoodcol == 0) stop(" No columns remain, so we stop.")
      if (numgoodcol == 1) stop(" Only 1 column remains, so we stop.")
    }    
    # Update array colInAnalysis using goodcol:
    colInAnalysis[colInAnalysis == TRUE] <- goodcol # overwrites correctly
    if (!silent) {
      wnq(" ")
      pnq(names(which(colInAnalysis)))
    }
  }
  
  
  if (cleanNAfirst == "automatic") {
    if (dim(remX)[2] >=  5 * dim(remX)[1]) { # more variables than observations
      cleanNAfirst = "columns"
    } else {
      cleanNAfirst = "rows"
    }
  }
  
  if (cleanNAfirst == "columns") {
    # 3. Deselect variables with over fracNA% of missing values
    #    (e.g. fracNA=0.20). Then update the vector colInAnalysis.
    #
    remX[!is.finite(remX)] <- NA # sets NA, NaN, Inf, -Inf all to NA
    acceptNA <- nrow(remX) * fracNA
    NAcounts   <- colSums(is.na(remX))
    goodcol    <- (NAcounts <= acceptNA)
    numgoodcol <- sum(goodcol)
    vecNAcol   <- (goodcol == FALSE)
    numNAcol   <- sum(vecNAcol)
    namesNAcol <- NULL
    if (numNAcol > 0) {
      if (!silent) {
        wnq(" ")
        wnq(paste(" The data contained ", numNAcol, " columns with over ",
                  round(100 * fracNA, 2), "% of NAs.", sep = ""))
        wnq(" Their column names are:")
        wnq(" ")    
      }
      namesNAcol <- colnames(remX)[vecNAcol]
      if (!silent) {
        pnq(namesNAcol)    
        wnq(" ")    
      }
      if (numgoodcol > 1) {
        if (!silent) {
          wnq(" These columns will be ignored in the analysis.")
          wnq(paste(" We continue with the remaining ", numgoodcol,
                    " columns:", sep = ""))
        }
        remX <- remX[, goodcol, drop = FALSE]
      } else { 
        if (numgoodcol == 0) stop(" No columns remain, so we stop.")
        if (numgoodcol == 1) stop(" Only 1 column remains, so we stop.")
      }    
      colInAnalysis[colInAnalysis == TRUE] <- goodcol
      if (!silent) {
        wnq(" ")
        pnq(names(which(colInAnalysis)))
      }
    }
    
    
    # 4. Deselect rows with too many NAs.
    #    Create the vector rowInAnalysis.
    #
    acceptNA   <- ncol(remX) * fracNA
    NAcounts   <- rowSums(is.na(remX))
    goodrow    <- (NAcounts <= acceptNA)
    numgoodrow <- sum(goodrow)
    vecNArow   <- (goodrow == FALSE)
    numNArow   <- sum(vecNArow)
    rowInAnalysis <- goodrow # in case we need to remove more rows later.
    namesNArow <- NULL
    if (numNArow > 0) {
      if (!silent) {
        wnq(" ")
        wnq(paste(" The data contained ", numNArow, " rows with over ",
                  round(100 * fracNA, 2), "% of NAs.", sep = ""))
        wnq(" Their row names are:")
        wnq(" ")
      }
      namesNArow <- rownames(remX)[vecNArow]
      if (!silent) {
        pnq(namesNArow)    
        wnq(" ")
      }
      if (numgoodrow > 2) {
        if (!silent) {
          wnq(" These rows will be ignored in the analysis.")
          wnq(paste(" We continue with the remaining ", numgoodrow,
                    " rows:", sep = "")) 
        }
        remX <- remX[goodrow, , drop = FALSE]
      } else { 
        if (numgoodrow == 0) stop(" No rows remain, so we stop.")
        if (numgoodrow == 1) stop(" Only 1 row remains, so we stop.")
        if (numgoodrow == 2) stop(" Only 2 rows remain, so we stop.")      
      }
      if (!silent) {
        wnq(" ")    
        pnq(names(which(rowInAnalysis)))
      }
    }
    
  } else {
    
    # 3. Deselect rows with too many NAs.
    #    Create the vector rowInAnalysis.
    remX[!is.finite(remX)] <- NA # sets NA, NaN, Inf, -Inf all to NA
    acceptNA   <- ncol(remX) * fracNA
    NAcounts   <- rowSums(is.na(remX))
    goodrow    <- (NAcounts <= acceptNA)
    numgoodrow <- sum(goodrow)
    vecNArow   <- (goodrow == FALSE)
    numNArow   <- sum(vecNArow)
    rowInAnalysis <- goodrow # in case we need to remove more rows later.
    namesNArow <- NULL
    if (numNArow > 0) {
      if (!silent) {
        wnq(" ")
        wnq(paste(" The data contained ", numNArow, " rows with over ",
                  round(100 * fracNA, 2), "% of NAs.", sep = ""))
        wnq(" Their row names are:")
        wnq(" ")
      }
      namesNArow <- rownames(remX)[vecNArow]
      if (!silent) {
        pnq(namesNArow)    
        wnq(" ")
      }
      if (numgoodrow > 2) {
        if (!silent) {
          wnq(" These rows will be ignored in the analysis.")
          wnq(paste(" We continue with the remaining ", numgoodrow,
                    " rows:", sep = "")) 
        }
        remX <- remX[goodrow, , drop = FALSE]
      } else { 
        if (numgoodrow == 0) stop(" No rows remain, so we stop.")
        if (numgoodrow == 1) stop(" Only 1 row remains, so we stop.")
        if (numgoodrow == 2) stop(" Only 2 rows remain, so we stop.")      
      }
      if (!silent) {
        wnq(" ")    
        pnq(names(which(rowInAnalysis)))
      }
    }
    
    # 4. Deselect variables with over fracNA% of missing values
    #    (e.g. fracNA=0.20). Then update the vector colInAnalysis.
    #
    acceptNA <- nrow(remX) * fracNA
    NAcounts   <- colSums(is.na(remX))
    goodcol    <- (NAcounts <= acceptNA)
    numgoodcol <- sum(goodcol)
    vecNAcol   <- (goodcol == FALSE)
    numNAcol   <- sum(vecNAcol)
    namesNAcol <- NULL
    if (numNAcol > 0) {
      if (!silent) {
        wnq(" ")
        wnq(paste(" The data contained ", numNAcol, " columns with over ",
                  round(100 * fracNA, 2), "% of NAs.", sep = ""))
        wnq(" Their column names are:")
        wnq(" ")    
      }
      namesNAcol <- colnames(remX)[vecNAcol]
      if (!silent) {
        pnq(namesNAcol)    
        wnq(" ")    
      }
      if (numgoodcol > 1) {
        if (!silent) {
          wnq(" These columns will be ignored in the analysis.")
          wnq(paste(" We continue with the remaining ", numgoodcol,
                    " columns:", sep = ""))
        }
        remX <- remX[, goodcol, drop = FALSE]
      } else { 
        if (numgoodcol == 0) stop(" No columns remain, so we stop.")
        if (numgoodcol == 1) stop(" Only 1 column remains, so we stop.")
      }    
      colInAnalysis[colInAnalysis == TRUE] <- goodcol
      if (!silent) {
        wnq(" ")
        pnq(names(which(colInAnalysis)))
      }
    }
  }
  
  
  # 5. Deselect discrete variables, loosely defined as variables that
  #    take on numDiscrete or fewer values, such as binary variables.
  #
  countValues <- function(colj) {
    # assumes colj is a vector
    # length(unique(colj))    # counts NA as a value    
    sum(!is.na(unique(colj))) # only counts non-NAs
  }
  valueCount  <- apply(remX, 2, countValues)
  goodcol     <- (valueCount > numDiscrete)
  numgoodcol  <- sum(goodcol)
  vecbadcol   <- (goodcol == F)
  numbadcol   <- sum(vecbadcol) # can be 0
  namesDiscrete <- NULL
  if (numbadcol > 0) {
    if (!silent) {
      wnq(" ")
      wnq(paste(" The data contained ", numbadcol, " discrete columns with ",
                numDiscrete, " or fewer values.", sep = ""))
      wnq(" Their column names are:")
      wnq(" ")
    }
    namesDiscrete <- colnames(remX)[vecbadcol]
    if (!silent) {
      pnq(namesDiscrete)    
      wnq(" ")
    }
    if (numgoodcol > 1) {
      if (!silent) {
        wnq(" These columns will be ignored in the analysis.")
        wnq(paste(" We continue with the remaining ", numgoodcol,
                  " columns:", sep = "")) 
      }
      remX <- remX[, goodcol, drop = FALSE]
    } else { 
      if (numgoodcol == 0) stop(" No columns remain, so we stop.")
      if (numgoodcol == 1) stop(" Only 1 column remains, so we stop.")
    }
    colInAnalysis[colInAnalysis == TRUE] <- goodcol
    if (!silent) {
      wnq(" ")
      pnq(names(which(colInAnalysis)))
    }
  }
  
  # 6. Deselect columns for which the median absolute deviation is
  #    zero. This is equivalent to saying that 50% or more of its
  #    values are equal.
  #
  colScale    <- apply(remX, 2, mad, na.rm = TRUE)
  goodcol     <- (colScale > precScale)
  numgoodcol  <- sum(goodcol)
  vecbadcol   <- (goodcol == FALSE)
  numbadcol   <- sum(vecbadcol) # can be 0
  namesZeroScale <- NULL
  if (numbadcol > 0) {
    if (!silent) {
      wnq(" ")
      wnq(paste(" The data contained ", numbadcol, " columns with zero",
                " or tiny median absolute deviation.", sep = ""))
      wnq(" Their column names are:")
      wnq(" ")
    }
    namesZeroScale <- colnames(remX)[vecbadcol]
    if (!silent) {
      pnq(namesZeroScale)
      wnq(" ")
    }
    if (numgoodcol > 1) {
      if (!silent) {
        wnq(" These columns will be ignored in the analysis.")
        wnq(paste(" We continue with the remaining ", numgoodcol,
                  " columns:", sep = "")) 
      }
      remX <- remX[, goodcol, drop = FALSE]
    } else { 
      if (numgoodcol == 0) stop(" No columns remain, so we stop.")
      if (numgoodcol == 1) stop(" Only 1 column remains, so we stop.")
    }    
    colInAnalysis[colInAnalysis == TRUE] <- goodcol
    if (!silent) {
      wnq(" ")
      pnq(names(which(colInAnalysis)))
    }
  }
  
  # check whether we have reduced the size of X
  if (nrow(remX) < n | ncol(remX) < d) {
      wnq(" ")
      wnq(paste(" The final data set we will analyze has ",
                nrow(remX), " rows and ", ncol(remX),
                " columns.", sep = ""))
      wnq(" ")
  }
  
  
  return(list(colInAnalysis = which(colInAnalysis),
              rowInAnalysis = which(rowInAnalysis),
              namesNotNumeric = namesNotNumeric,
              namesCaseNumber = namesCaseNumber,
              namesNAcol = namesNAcol,
              namesNArow = namesNArow,
              namesDiscrete = namesDiscrete,
              namesZeroScale = namesZeroScale,
              remX = remX))
}

