.yamlpretty <- function(node, clevels, nodeID) {
    if (node$Type == 'Terminal') {
        node$Size = sum(nodeID == node$ID)
        return(node)
    } else {
        if (node$Role != 'num') {
            node$ThreshSet = clevels[[node$SplitVar]][as.integer(node$ThreshSet)]
        }
        node$Left = .yamlpretty(node$Left, clevels, nodeID)
        node$Right = .yamlpretty(node$Right, clevels, nodeID)
        return(node)
    }
}

#' print node information
#'
#' @param node node object
#' @param depth current depth
#' @param digits digit print
#' @param long default=TRUE useless
#' @param yName outcomes names
#' @param trtName treatment name
#' @param tlevels treatment levels
#' @param clevels categorical variables' levels
#' @param ... pass to cat
#'
#' @noRd
#'
#'
print_node <- function(node, depth = 0, digits = 3,
                       long = TRUE, yName, trtName,
                       tlevels, clevels, ...) {
    if (node$Type == 'Terminal') {
        cat(rep(' ', depth), 'ID: ', node$ID, ', Size: ', node$Size, ' [Terminal]\n', sep = '', ...)
        if (long) {
            cat(rep(' ', depth), 'Outcome Models: \n', sep = '', ...)
            for(i in seq_along(yName)) {
                cat(rep(' ', depth + 4), yName[i], rep(' ', 8), 'Est', rep(' ', 8), 'SE\n', sep = '', ...)
                # cat(rep(' ', depth + 4), format(c(yName[i], 'Estimate', 'SE'), justify = 'centre'), '\n', sep = '', ...)
                xvar <- node$FitIndex[[i]]
                #coefT <- matrix(NA, length(tlevels), 2)
                #rownames(coefT) <- paste0(trtName, '.', tlevels)
                #colnames(coefT) <- c('Estimate', 'Std.Err')
                if (length(xvar) > 0) {
                    for (j in seq_along(xvar)) {
                        cat(rep(' ', depth + 4), xvar[j], '    ',
                            round(node$Parms[[i]][j], digits), '\n', sep = '', ...)
                    }
                }

                for (j in seq_along(tlevels)) {
                    #coefT[j, 1] <- node$Trts[[i]][j]
                    #coefT[j, 2] <- node$SEs[[i]][j]
                    cat(rep(' ', depth + 4), paste0(trtName,'.', tlevels[j]), '    ',
                        round(node$Trts[[i]][j], digits), '    ',
                        round(node$SEs[[i]][j], digits),'\n', sep = '', ...)
                    # cat(rep(' ', depth + 4), format(c(paste0(trtName,'.', tlevels[j]),
                    #                                   round(node$Trts[[i]][j], digits), round(node$SEs[[i]][j], digits)), justify="centre"), '\n', sep = '', ...)
                }
                cat(rep(' ', depth), rep('- ', 14), '\n', sep = '', ...)
            }
        }

    } else {
        cat(rep(' ', depth), 'ID: ', node$ID, ', ', sep = '', ...)
        if (node$Role == 'num') {
            if (node$MisDirection != 'A') {
                cat(node$SplitVar, ' <=', ifelse(node$MisDirection == 'L', '* ', ' '), round(node$Threshold, digits), '\n', sep = '', ...)
            } else {
                cat(node$SplitVar, ' = NA\n', sep = '', ...)
            }
        } else {
            if (node$MisDirection != 'A') {
                cat(node$SplitVar, ' = { ', paste0(node$ThreshSet, collapse = ', '), ifelse(node$MisDirection == 'L', ', NA', ''), ' }\n', sep = '', ...)
            } else {
                cat(node$SplitVar, ' = NA\n', sep = '', ...)
            }

        }
        print_node(node$Left, depth + 4, digits, long, yName, trtName, tlevels, clevels, ...)

        cat(rep(' ', depth), 'ID: ', node$ID, ', ' , sep = '', ...)
        if (node$Role == 'num') {
            if (node$MisDirection != 'A') {
                cat(node$SplitVar, ' >', ifelse(node$MisDirection == 'L', ' ', '* '), round(node$Threshold, digits), '\n', sep = '', ...)
            } else {
                cat(node$SplitVar, ' != NA\n', sep = '', ...)
            }
        } else {
            if (node$MisDirection != 'A') {
                varLevel = clevels[[node$SplitVar]]
                cat(node$SplitVar, ' = { ', paste0(varLevel[which(!varLevel %in% node$ThreshSet)], collapse = ', '), ifelse(node$MisDirection == 'L', '', ', NA'), ' }\n', sep = '', ...)
            } else {
                cat(node$SplitVar, ' != NA\n', sep = '', ...)
            }
        }
        print_node(node$Right, depth + 4, digits, long, yName, trtName, tlevels, clevels, ...)
    }
}

#' Print fitted regression tree
#'
#' @param mrsobj MrSGUIDE object
#' @param digits digits pass to coefficient
#' @param details whether to print fitting details
#' @param ... parameter pass to \code{print_node}
#'
#' @return print tree information into console
#'
#' @examples
#' library(MrSGUIDE)
#' set.seed(1234)
#'
#' N = 200
#' np = 3
#'
#' numX <- matrix(rnorm(N * np), N, np) ## numerical features
#' gender <- sample(c('Male', 'Female'), N, replace = TRUE)
#' country <- sample(c('US', 'UK', 'China', 'Japan'), N, replace = TRUE)
#'
#' z <- sample(c(0, 1), N, replace = TRUE) # Binary treatment assignment
#'
#' y1 <- numX[, 1] + 1 * z * (gender == 'Female') + rnorm(N)
#' y2 <- numX[, 2] + 2 * z * (gender == 'Female') + rnorm(N)
#'
#' train <- data.frame(numX, gender, country, z, y1, y2)
#' role <- c(rep('n', 3), 'c', 'c', 'r', 'd', 'd')
#'
#' mrsobj <- MrSFit(dataframe = train, role = role)
#' printTree(mrsobj, digits = 2, details=TRUE)
#' printTree(mrsobj, digits = 2, details=FALSE)
#'
#' @export
printTree <- function(mrsobj, digits = 3, details = TRUE, ...) {
    print_node(mrsobj$treeRes, depth = 0, digits, details,
               mrsobj$ynames, mrsobj$trtname, mrsobj$tLevels[[1]], mrsobj$cLevels, ...)
}

#' Write Latex file for GUIDE regression Tree
#'
#' @param mrsobj MrSGUIDE object
#' @param file latex filename
#' @param digits digits pass to coefficient
#' @param ... parameters pass to cat function
#'
#' @return write txt file into disk
#'
#' @examples
#'
#' library(MrSGUIDE)
#' set.seed(1234)
#'
#' N = 200
#' np = 3
#'
#' numX <- matrix(rnorm(N * np), N, np) ## numerical features
#' gender <- sample(c('Male', 'Female'), N, replace = TRUE)
#' country <- sample(c('US', 'UK', 'China', 'Japan'), N, replace = TRUE)
#'
#' z <- sample(c(0, 1), N, replace = TRUE) # Binary treatment assignment
#'
#' y1 <- numX[, 1] + 1 * z * (gender == 'Female') + rnorm(N)
#' y2 <- numX[, 2] + 2 * z * (gender == 'Female') + rnorm(N)
#'
#' train <- data.frame(numX, gender, country, z, y1, y2)
#' role <- c(rep('n', 3), 'c', 'c', 'r', 'd', 'd')
#'
#' mrsobj <- MrSFit(dataframe = train, role = role)
#' \dontshow{.old_wd <- setwd(tempdir())}
#' writeTex(mrsobj, 'test.tex')
#' \dontshow{setwd(.old_wd)}
#'
#' @export
writeTex <- function(mrsobj, file, digits = 3, ...) {
    cat('\\documentclass[12pt]{article}\n', file = file, sep = "", ...)
    cat(' %File creation date:', as.character(Sys.time()), '\n', file = file, append = TRUE,sep = "", ...)
    cat('\\usepackage{pstricks,pst-node,pst-tree}\n', file = file, append = TRUE,sep = "", ...)
    cat('\\usepackage{geometry}
\\usepackage{lscape}
\\pagestyle{empty}
\\begin{document}
%\\begin{landscape}
\\begin{center}
\\psset{linecolor=black,tnsep=1pt,tndepth=0cm,tnheight=0cm,treesep=.8cm,levelsep=50pt,radius=10pt}\n',  file = file, append = TRUE, sep = "", ...)
    treatNode <- .getTrt(mrsobj$nodeMap, mrsobj$ynames, mrsobj$tLevels[[1]])
    .writetex(mrsobj$treeRes, file, treatNode, digits, ...)
    cat("\\end{center}
%\\end{landscape}
\\end{document}\n", file = file, append = TRUE, sep = "", ...)
}

.writetex <- function(node, texfile, treatNode, depth = 0, digits = 3, ...) {
    if (node$Type == 'Terminal') {
        trt = treatNode[treatNode$Node == node$ID, 'Estimate']
        fillcolor = ifelse(mean(trt) > 0, "red", "yellow")
        cat(rep(' ', depth), "\\Tcircle[fillcolor=",fillcolor,",fillstyle=solid]{ ", node$ID, " }~{\\makebox[0pt][c]{\\em ",
            node$Size, " }}\n", file = texfile, append = TRUE, sep = "", ...)
        return()
    } else {
        cat(rep(' ', depth), "\\pstree[treemode=D]{\\Tcircle{ ", node$ID,
            " }~[tnpos=l]{\\shortstack[r]{\\texttt{\\detokenize{",
            node$SplitVar, "}}\\\\", file = texfile, append = TRUE, sep = "", ...)

        if (node$Role == 'num') {
            if (node$MisDirection != 'A') {
                cat('$\\leq', ifelse(node$MisDirection == 'L', '_*$', '$'),
                    round(node$Threshold, digits), '}}\n', sep = '', file = texfile, append = TRUE, ...)
            } else {
                cat('$=$', 'NA}}\n', sep = '', file = texfile, append = TRUE, ...)
            }
        } else {
            if (node$MisDirection != 'A') {
                cat('$\\in$ \\{ ', paste0(node$ThreshSet, collapse = ', '), ifelse(node$MisDirection == 'L', ', NA', ''), '\\}}}\n', sep = '', file = texfile, append = TRUE, ...)
            } else {
                cat('$=$', 'NA}}\n', sep = '', file = texfile, append = TRUE,...)
            }
        }
        cat(rep(' ', depth), "}{\n", file = texfile, sep = '', append = TRUE, ...)
        .writetex(node$Left, texfile, treatNode, depth = depth + 4,  digits = 3, ...)
        .writetex(node$Right, texfile, treatNode, depth = depth + 4, digits = 3, ...)
        cat(rep(' ', depth), "}\n", file = texfile, sep = '', append = TRUE, ...)
    }
}
