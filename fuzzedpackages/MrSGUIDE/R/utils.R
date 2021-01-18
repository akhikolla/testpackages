.charLevels <- function(x) {
    tmp = as.character(x)
    tmp[which(is.na(x))] = "12345679"
    sort(unique(tmp))
}

.dictRename <- function(dictList, cvarName = NULL) {
    if (!is.null(cvarName) && is.null(names(dictList))) names(dictList) = cvarName
    .naTchar <- function(x) {
        if (any(is.na(x))) {
            ind = which(is.na(x))
            x[ind] = "12345679"
        }
        return(x)
    }

    lapply(dictList, .naTchar)
}


.getAlpha <- function(bootAlpha) {
    res <- c(0.025, 0.025, 0.025)
    if (!is.null(bootAlpha)) {
        alphaVec <- bootAlpha[, 1]
        ind <- bootAlpha[, 2]
        ove <- bootAlpha[, 3]
        alpha <- max(alphaVec) * 2
        .helper <- function(gamma, alphaVec, alpha = 0.05) {
            p <- which(gamma < 1 - alpha)
            if (p[1] == 1) return(alpha)
            else {
                p = p[1]
                f = (gamma[p - 1] - 1 + alpha) / (gamma[p - 1] - gamma[p])
                alpha = (1 - f) * alphaVec[p - 1] + f * alphaVec[p]
            }
            return(alpha)
        }
        indAlpha = .helper(ind, alphaVec, alpha)
        oveAlpha = .helper(ove, alphaVec, alpha)
        res <- c(alpha / 2, indAlpha, oveAlpha)
    }
}

#' Check the condition and print message
#' @param logicstatus whehter it is true
#' @param message message want to print
#' @noRd
.check <- function(logicstatus, message) {
    if (!logicstatus)
        stop(message, call. = FALSE)
}

.check_file <- function(file) {
    stopifnot(!is.na(file))
    stopifnot(!is.null(file))
}

.data_process <- function(dataframe, role) {
    .check(NCOL(dataframe) == length(role),
           "column number of data frame is not same as role length.")
    varName = colnames(dataframe)

    dr = which(role == "d")

    .check(length(dr) > 0, "No dependent variable in role.")

    rr = which(role == "r")

    .check(length(rr) > 0, "No treatment variable in role.")
    .check(length(rr) < 2, paste("Current version can only deal with one treatment variable. Find ", length(rr), "in role."))

    sr = which(role %in% c("n", "s"))
    fr = which(role %in% c("n", "f"))
    hr = which(role %in% c("h"))

    nr = sort(union(union(sr, fr), hr))
    splitIndex = which(nr %in% sr) - 1
    fitIndex = which(nr %in% fr) - 1
    holdIndex = which(nr %in% hr) - 1

    cr = which(role %in% "c")
    xr = which(role == "x")

    .check(length(sr) + length(cr) > 0, "No split variable in role.")

    for (i in c(sr, fr, hr)) {
        .check(!class(dataframe[, i]) %in% c("character", "factor"),
               paste(varName[i], "seems a categorical variable? please change the role vector") )
    }

    for (i in cr) {
        .check(class(dataframe[, i]) != "numeric",
               paste(varName[i], "seems not a categorical variable. Please change the role vector"))
    }

    for (i in dr) {
        if (class(dataframe[, i]) %in% c("character", "factor"))
            warning(paste(varName[i], "seems a categorical variable? Will force to treat as numerical variable"))
    }

    Yori = as.matrix(dataframe[, dr])
    missInd <- !is.na(Yori)

    if (NCOL(missInd) > 1) {
        non_miss <- as.logical(apply(missInd, 1, prod))
    } else {
        non_miss <- missInd
    }

    if (sum(non_miss) < NROW(Yori)) {
        warning(paste(" Dependent variables have missing value.\n Only use non-missing record for subgroup identification.\n Total is: ", sum(non_miss), "Original: ", NROW(Yori)))
    }

    Y = as.matrix(Yori[non_miss, ])
    cLevels = lapply(dataframe[non_miss,][cr], .charLevels)
    cXL = characterDict(dataframe[non_miss, cr], cLevels)
    nX = dataFramToNumeric(dataframe[non_miss, nr])

    tLevels = list(.charLevels(dataframe[non_miss, rr]))
    TrtL = characterDict(dataframe[non_miss, rr], tLevels)
    numVarName = varName[nr]
    catVarName = varName[cr]
    newVar = c(numVarName, catVarName)

    res = list(
        Y = Y,
        cLevels = cLevels,
        cXL = cXL,
        nX = nX,
        tLevels = tLevels,
        TrtL = TrtL,
        numVarName = numVarName,
        catVarName = catVarName,
        newVar = newVar,
        splitIndex = splitIndex,
        holdIndex = holdIndex,
        fitIndex = fitIndex,
        non_miss = non_miss,
        varName = varName
    )
}


#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`
