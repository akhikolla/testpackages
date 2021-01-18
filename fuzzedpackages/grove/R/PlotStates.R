#' @title Function to plot the hidden states
#'
#' @description This function plots on a tree the state of each latent variables.
#'
#' @param grove.obj Output from function \code{FAnova}.
#' @param block Which block to plot.
#' @param legend If \code{TRUE}, show legend. 
#' @param main Main title.
#' @param prior If \code{TRUE}, plot prior state probabilities. If \code{FALE},
#' plot posterior state probabilities. 
#' 
#' @return A plot.
#' @export
#' @examples \dontrun{
#' data <- GenerateSyntheticAnova(st.dev = 5, n.replicates = 5)
#' W <- DWT(data$noisy.Y)
#' X <- data$X
#' ans <- FAnova(W, X, ~ 1 + factorA + factorB)
#' PlotStates(ans)
#' PlotStates(ans, block = "factorA")
#' PlotStates(ans, block = "factorB")}

PlotStates <- function(grove.obj, 
                       block = "Intercept", 
                       legend = FALSE, 
                       main = NULL, 
                       prior = FALSE) {
  
  if (class(grove.obj) != "grove")  {
    stop("Input should be a grove class object")
  }
  
  temp <- gsub(" + ", 
               "+", 
               as.character(grove.obj$data$formula)[2], fixed = TRUE)
  temp <- strsplit(temp, split = "+", fixed = TRUE)[[1]]
  if(block == "Intercept")
    index <- 1
  else  {
    index <- which(temp == block) 
    if(length(index) == 0) {
      stop("No block with such a name")
    }    
  }

  p <- grove.obj$data$p
  pos <- rep(0, 2 ^ length(p))
  for (i in 0:(2 ^ (length(p)) - 1))  {
    pos[i + 1] = as.integer(intToBits(i))[index]
  }
  states <- which(pos == 1)

  if (prior) {
    S <- grove.obj$samples$Sprior
  } else {
    S <- grove.obj$samples$S
  }
  J <- length(S)
  colori <- colorRampPalette(c("gray85", "darkblue"))(100)  
  idx <- 2 ^ (0:(J - 1))
  ylim <- c(ifelse(legend, -2, 0), J + 2)
  plot(1, 
       type = "n", 
       axes = F, 
       xlab = "",
       ylab = "scale", 
       xlim = c(0, 1),
       ylim = ylim)
     
  if (is.null(main)) {
    title(main = block)
  } else {
    title(main = main)
  } 
  i <- 0
  for (j in 0 : (J - 1))  {
     
    if (j == 0) {
      col <- colori[rev(ceiling(sum(S[[j + 1]][states]) * 99))]
    } else {
      col <- colori[rev(ceiling(colSums(S[[j + 1]][states, ]) * 99))]
    }
    
    points(seq(0, 1, by = 1 / (2 ^ j + 1))[2 : (2 ^ j + 1)],
           rep(J - j, 2 ^ j), 
           col = col, 
           pch = 15, 
           cex = 1.5)
    
    if (j == 0) {
      i <- i + 2 ^ j 
    }
  }
  axis(2, at = 1 : J, labels = rev(0 : (J - 1)), las = 2)
  
  if (legend) {
    x <- seq(0.2, 0.95, length=11)
    labels <- seq(1, 100, length = 11)
    points(x, 
           rep(-0.5, 11), 
           col = colori[round(labels)], 
           pch = 15, 
           cex = 1.5)
    text(x, -0.7, labels = labels, pos = 1)
    rect(0.0, -2.2, 1, 0.1)
    if (prior) {
      if (block == "Intercept") {
        text(0.1, -1.1 , expression(P(S[jk]==1)))
      } else {
        text(0.1, -1.1 , expression(P(R[jk]==1)))
      }
    } else {
      if (block == "Intercept") {
        text(0.1, -1.1 , expression(P(S[jk]==1~'|'~D)))
      } else {
        text(0.1, -1.1 , expression(P(R[jk]==1~'|'~D)))       
      }
    }    
  }
}

