

calcOverlap <- function(conditions, which.conditions = NULL)
{
    if (is.null(dim(conditions)))
    {
        return(sum(drop(conditions)))
    }
    if (is.null(which.conditions))
    {
        which.conditions <- 1:ncol(conditions)
    }
    if (any(which.conditions > ncol(conditions)))
    {
        stop("too many conditions requested")
    }
    ## looks redundant but this handles factors
    sum(apply(1 *(conditions[,which.conditions, drop = FALSE] == 1), 1, prod) == 1)
}


#' plotting function for venn diagrams of overlapping conditions
#'
#' @param conditions condition matrix such as the one given to vennLasso() function. It can have up to 5 conditions
#' @param condition.names names of the conditions (equal to the number of columns of conditions)
#' @param fill.colors vector of colors for plotting. Set fill.colors = NULL for no colors
#' @param lty standard 'lty' graphical parameter for line type around circles. default is no lines
#' @param ... other graphical parameters for the plot
#' @import VennDiagram
#' @importFrom grid grid.newpage
#' @export
#' @examples
#' library(Matrix)
#' 
#' set.seed(123)
#' n.obs <- 200
#' n.vars <- 50
#'
#' true.beta.mat <- array(NA, dim = c(3, n.vars))
#' true.beta.mat[1,] <- c(-0.5, -1, 0, 0, 2, rep(0, n.vars - 5))
#' true.beta.mat[2,] <- c(0.5, 0.5, -0.5, -0.5, 1, -1, rep(0, n.vars - 6))
#' true.beta.mat[3,] <- c(0, 0, 1, 1, -1, rep(0, n.vars - 5))
#' rownames(true.beta.mat) <- c("1,0", "1,1", "0,1")
#' true.beta <- as.vector(t(true.beta.mat))
#'
#' x.sub1 <- matrix(rnorm(n.obs * n.vars), n.obs, n.vars)
#' x.sub2 <- matrix(rnorm(n.obs * n.vars), n.obs, n.vars)
#' x.sub3 <- matrix(rnorm(n.obs * n.vars), n.obs, n.vars)
#'
#' x <- as.matrix(rbind(x.sub1, x.sub2, x.sub3))
#'
#' conditions <- as.matrix(cbind(c(rep(1, 2 * n.obs), rep(0, n.obs)),
#'                               c(rep(0, n.obs), rep(1, 2 * n.obs))))
#'
#' y <- rnorm(n.obs * 3, sd = 3) + drop(as.matrix(bdiag(x.sub1, x.sub2, x.sub3)) %*% true.beta)
#'
#' fit <- vennLasso(x = x, y = y, groups = conditions)
#'
#' vennobj <- plotVenn(conditions)
#'
plotVenn <- function(conditions,
                     condition.names = NULL,
                     lty = "blank",
                     fill.colors = c("royalblue1", "goldenrod1", "mediumvioletred", "turquoise3", "firebrick1"),
                     ...) {
    grid.newpage()

    ## really boring case. why even plot this?
    if (ncol(conditions) == 1)
    {
        conditions <- drop(conditions)
        out <- draw.single.venn(sum(conditions),
                                fill = fill.colors[1],
                                lty = lty,
                                ...)
    }

    if (!is.null(fill.colors))
    {
        if (length(fill.colors) < ncol(conditions))
        {
            # add more colors if needed
            fill.colors <- c(fill.colors, c("royalblue1", "goldenrod1", "mediumvioletred", "turquoise3", "firebrick1"))
        }
    }

    if (is.null(condition.names))
    {
        condition.names <- colnames(conditions)
        if (is.null(condition.names))
        {
            condition.names <- paste("Condition:", 1:ncol(conditions))
        }
    } else
    {
        if (length(condition.names) != ncol(conditions))
        {
            stop("Must provide number of condition names
                 equal to the number of columns in conditions matrix")
        }
    }


    if (ncol(conditions) == 2)
    {
        out <- draw.pairwise.venn(calcOverlap(conditions, 1),
                                  calcOverlap(conditions, 2),
                                  calcOverlap(conditions, c(1, 2)),
                                  category = condition.names,
                                  fill     = fill.colors[1:2],
                                  lty      = lty,
                                  ...)
    }
    if (ncol(conditions) == 3)
    {
        out <- draw.triple.venn(calcOverlap(conditions, c(1)),
                                calcOverlap(conditions, c(2)),
                                calcOverlap(conditions, c(3)),
                                calcOverlap(conditions, c(1, 2)),
                                calcOverlap(conditions, c(2, 3)),
                                calcOverlap(conditions, c(1, 3)),
                                calcOverlap(conditions, c(1, 2, 3)),
                                category = condition.names,
                                fill     = fill.colors[1:3],
                                lty      = lty,
                                ...)
    }
    if (ncol(conditions) == 4)
    {
        out <- draw.quad.venn(calcOverlap(conditions, c(1)),
                              calcOverlap(conditions, c(2)),
                              calcOverlap(conditions, c(3)),
                              calcOverlap(conditions, c(4)),
                              calcOverlap(conditions, c(1, 2)),
                              calcOverlap(conditions, c(1, 3)),
                              calcOverlap(conditions, c(1, 4)),
                              calcOverlap(conditions, c(2, 3)),
                              calcOverlap(conditions, c(2, 4)),
                              calcOverlap(conditions, c(3, 4)),
                              calcOverlap(conditions, c(1, 2, 3)),
                              calcOverlap(conditions, c(1, 2, 4)),
                              calcOverlap(conditions, c(1, 3, 4)),
                              calcOverlap(conditions, c(2, 3, 4)),
                              calcOverlap(conditions, c(1, 2, 3, 4)),
                              category = condition.names,
                              fill     = fill.colors[1:4],
                              lty      = lty,
                              ...)
    }
    if (ncol(conditions) == 5)
    {
        out <- draw.quintuple.venn(
            calcOverlap(conditions, c(1)),
            calcOverlap(conditions, c(2)),
            calcOverlap(conditions, c(3)),
            calcOverlap(conditions, c(4)),
            calcOverlap(conditions, c(5)),
            calcOverlap(conditions, c(1, 2)),
            calcOverlap(conditions, c(1, 3)),
            calcOverlap(conditions, c(1, 4)),
            calcOverlap(conditions, c(1, 5)),
            calcOverlap(conditions, c(2, 3)),
            calcOverlap(conditions, c(2, 4)),
            calcOverlap(conditions, c(2, 5)),
            calcOverlap(conditions, c(3, 4)),
            calcOverlap(conditions, c(3, 5)),
            calcOverlap(conditions, c(4, 5)),
            calcOverlap(conditions, c(1, 2, 3)),
            calcOverlap(conditions, c(1, 2, 4)),
            calcOverlap(conditions, c(1, 2, 5)),
            calcOverlap(conditions, c(1, 3, 4)),
            calcOverlap(conditions, c(1, 3, 5)),
            calcOverlap(conditions, c(1, 4, 5)),
            calcOverlap(conditions, c(2, 3, 4)),
            calcOverlap(conditions, c(2, 3, 5)),
            calcOverlap(conditions, c(2, 4, 5)),
            calcOverlap(conditions, c(3, 4, 5)),
            calcOverlap(conditions, c(1, 2, 3, 4)),
            calcOverlap(conditions, c(1, 2, 3, 5)),
            calcOverlap(conditions, c(1, 2, 4, 5)),
            calcOverlap(conditions, c(1, 3, 4, 5)),
            calcOverlap(conditions, c(2, 3, 4, 5)),
            calcOverlap(conditions, c(1, 2, 3, 4, 5)),
            category = condition.names,
            fill     = fill.colors[1:5],
            lty      = lty,
            ...
        )
    }
    
    if (!exists("out"))
    {
        return("Too many conditions; cannot plot venn diagram")
    } else 
    {
        return(invisible(out))
    }
}






