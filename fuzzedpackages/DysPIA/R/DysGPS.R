#' @title DysGPS: Calculates Dysregulated gene pair score (DysGPS) for each gene pair
#'
#' @description Calculates Dysregulated gene pair score (DysGPS) for each gene pair. Two-sample Welch's T test of gene pairs between case and control samples. The package 'DysPIAData' including the background data is needed to be loaded.
#'
#' @param dataset Matrix of gene expression values (rownames are genes, columnnames are samples).
#' @param class.labels Vector of category labels.
#' @param controlcharacter Charactor of control group in the class labels.
#' @param casecharacter Charactor of case group in the class labels.
#' @param background Matrix of the gene pairs' background. The default is `combined_background`,
#'                        which includes real pathway gene pairs and randomly producted gene pairs. The 'combined_background' was incluede in 'DysPIAData'.
#' @return A vector of DysGPS for each gene pair.
#' @import DysPIAData
#' @importFrom stats na.omit t.test var
#' @examples
#' data(gene_expression_p53, class.labels_p53,sample_background)
#' DysGPS_sample<-DysGPS(gene_expression_p53, class.labels_p53,
#'  "WT", "MUT", sample_background)
#'
DysGPS <- function(dataset, class.labels, controlcharacter, casecharacter,
                          background = combined_background) {

  caseloca <- which(class.labels == casecharacter)
  controlloca <- which(class.labels == controlcharacter)
  caseSet <- dataset[, caseloca]
  contrSet <- dataset[, controlloca]

  data_var <- apply(caseSet, 1, var) * apply(contrSet, 1, var)
  if (sum(data_var == 0) > 0){
    dataset <- dataset[-which(data_var == 0), ]
  }
  rm(caseSet, contrSet)

  # T test of gene pairs between case and control samples
  T_test <- function(t, c1, c2) {
    result <- t.test(t[c1], t[c2])
    return(result$statistic)
  }

  score <- as.numeric()
  group_num <- ceiling(nrow(background)/100000)

  for (i in 1:group_num){
    background_i <- background[(100000*(i-1)+1):min(100000*i, nrow(background)),]
    location <- matrix(0, nrow(background_i), 2)
    location[, 1] <- match(background_i[, 1], rownames(dataset))
    location[, 2] <- match(background_i[, 2], rownames(dataset))
    location <- na.omit(location)
    caseExp.1 <- dataset[location[, 1], caseloca]
    caseExp.2 <- dataset[location[, 2], caseloca]
    contrExp.1 <- dataset[location[, 1], controlloca]
    contrExp.2 <- dataset[location[, 2], controlloca]

    # Normalize gene expression data of case samples
    x1 <- scale(t(caseExp.1))
    x2 <- scale(t(caseExp.2))
    rm(caseExp.1, caseExp.2)
    x <- x1*x2
    colnames(x) <- paste0(colnames(x1), "|", colnames(x2))
    rm(x1, x2)
    y1 <- scale(t(contrExp.1))
    y2 <- scale(t(contrExp.2))
    rm(contrExp.1, contrExp.2)
    y<-y1*y2
    rm(y1, y2)
    colnames(y) <- colnames(x)

    z <- t(rbind(x, y))
    case_number <- nrow(x)
    contr_number <- nrow(y)
    rm(x, y)

    score_i <- apply(z, 1, T_test, c1=1:case_number,
                     c2=(case_number+1):(case_number+contr_number))

    rm(z)
    score <- c(score, score_i)
  }
    return(score)
}
