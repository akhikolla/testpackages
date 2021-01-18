# Wrapper function:

SILGGM <- function(x, method = NULL, lambda = NULL, global = FALSE, alpha = NULL, ndelta = NULL, true_graph = NULL, cytoscape_format = FALSE, csv_save = FALSE, directory = NULL){
  result <- SILGGMCpp(x, method, lambda, global, alpha, ndelta, true_graph)

  if(cytoscape_format == FALSE){
    return(result)
  }
  else{
    cytoscape_format_table = list()

    if(is.null(method) || method == "D-S_NW_SL"){
        pair <- melt(upper.tri(result$precision))
        pair <- pair[pair[,3] %in% TRUE,c(1,2)]
        row.names(pair) <- NULL
        precision <- result$precision[upper.tri(result$precision)]
        partialCor <- result$partialCor[upper.tri(result$partialCor)]
        CI_low_precision <- result$CI_low_precision[upper.tri(result$CI_low_precision)]
        CI_high_precision <- result$CI_high_precision[upper.tri(result$CI_high_precision)]
        z_score_precision <- result$z_score_precision[upper.tri(result$z_score_precision)]
        p_precision <- result$p_precision[upper.tri(result$p_precision)]
        if(global == FALSE){
          corr_table <- cbind(pair, precision, partialCor, CI_low_precision, CI_high_precision, z_score_precision, p_precision)
          colnames(corr_table) <- c("gene1","gene2","precision","partialCor","CI_low_precision","CI_high_precision","z_score_precision","p_precision")
          if(!is.null(colnames(x))){
            gene_symbol <- colnames(x)
            corr_table[,1] <- gene_symbol[as.numeric(corr_table[,1])]
            corr_table[,2] <- gene_symbol[as.numeric(corr_table[,2])]
          }
          cytoscape_format_table$cytoscape_format_table <- corr_table
          if(csv_save == TRUE){
            cat("Save the Cytoscape format table. \n")
            if(is.null(directory)){
              directory <- tempdir()
              cat("The directory is not specified, so we use a temporary directory from the R session. \n")
              cat("The directory is ", directory, ".\n",sep = "")
            }
            else{
              cat("The directory is specified and it is ", directory, ".\n",sep = "")
            }
            write.csv(corr_table,file=paste(directory, "Cytoscape_D-S_NW_SL.csv", sep = "/"),row.names = FALSE)
          }
          return(cytoscape_format_table)
        }
        else{
          if(is.null(true_graph)){
            cytoscape_format_table$threshold <- result$threshold
            cytoscape_format_table$FDR <- result$FDR
            global_decision <- NULL
            if(is.null(alpha)){
              alpha_length <- 2
            }
            else{
              alpha_length <- length(alpha)
            }
            for(i in 1:alpha_length){
              global_decision1 <- result$global_decision[[i]][upper.tri(result$global_decision[[i]])]
              global_decision <- cbind(global_decision, global_decision1)
            }
            if(is.null(alpha)){
              colnames(global_decision) <- c("global_decision_0.05","global_decision_0.1")
            }
            else{
              for(i in 1:alpha_length){
                colnames(global_decision)[i] <- paste("global_decision",alpha[i],sep="_")
              }
            }
            corr_table <- cbind(pair, precision, partialCor, CI_low_precision, CI_high_precision, z_score_precision, p_precision, global_decision)
            colnames(corr_table) <- c("gene1","gene2","precision","partialCor","CI_low_precision","CI_high_precision","z_score_precision","p_precision",colnames(global_decision))
            if(!is.null(colnames(x))){
              gene_symbol <- colnames(x)
              corr_table[,1] <- gene_symbol[as.numeric(corr_table[,1])]
              corr_table[,2] <- gene_symbol[as.numeric(corr_table[,2])]
            }
            cytoscape_format_table$cytoscape_format_table <- corr_table
            if(csv_save == TRUE){
              cat("Save the Cytoscape format table. \n")
              if(is.null(directory)){
                directory <- tempdir()
                cat("The directory is not specified, so we use a temporary directory from the R session. \n")
                cat("The directory is ", directory, ".\n",sep = "")
              }
              else{
                cat("The directory is specified and it is ", directory, ".\n",sep = "")
              }
              write.csv(corr_table,file=paste(directory, "Cytoscape_D-S_NW_SL.csv", sep = "/"),row.names = FALSE)
            }
            return(cytoscape_format_table)
          }
          else{
            cytoscape_format_table$threshold <- result$threshold
            cytoscape_format_table$FDR <- result$FDR
            cytoscape_format_table$power <- result$power
            global_decision <- NULL
            if(is.null(alpha)){
              alpha_length <- 2
            }
            else{
              alpha_length <- length(alpha)
            }
            for(i in 1:alpha_length){
              global_decision1 <- result$global_decision[[i]][upper.tri(result$global_decision[[i]])]
              global_decision <- cbind(global_decision, global_decision1)
            }
            if(is.null(alpha)){
              colnames(global_decision) <- c("global_decision_0.05","global_decision_0.1")
            }
            else{
              for(i in 1:alpha_length){
                colnames(global_decision)[i] <- paste("global_decision",alpha[i],sep="_")
              }
            }
            corr_table <- cbind(pair, precision, partialCor, CI_low_precision, CI_high_precision, z_score_precision, p_precision, global_decision)
            colnames(corr_table) <- c("gene1","gene2","precision","partialCor","CI_low_precision","CI_high_precision","z_score_precision","p_precision",colnames(global_decision))
            if(!is.null(colnames(x))){
              gene_symbol <- colnames(x)
              corr_table[,1] <- gene_symbol[as.numeric(corr_table[,1])]
              corr_table[,2] <- gene_symbol[as.numeric(corr_table[,2])]
            }
            cytoscape_format_table$cytoscape_format_table <- corr_table
            if(csv_save == TRUE){
              cat("Save the Cytoscape format table. \n")
              if(is.null(directory)){
                directory <- tempdir()
                cat("The directory is not specified, so we use a temporary directory from the R session. \n")
                cat("The directory is ", directory, ".\n",sep = "")
              }
              else{
                cat("The directory is specified and it is ", directory, ".\n",sep = "")
              }
              write.csv(corr_table,file=paste(directory, "Cytoscape_D-S_NW_SL.csv", sep = "/"),row.names = FALSE)
            }
            return(cytoscape_format_table)
          }
    }
    }

    if(method == "D-S_GL"){
      pair <- melt(upper.tri(result$precision))
      pair <- pair[pair[,3] %in% TRUE,c(1,2)]
      row.names(pair) <- NULL
      precision <- result$precision[upper.tri(result$precision)]
      partialCor <- result$partialCor[upper.tri(result$partialCor)]
      CI_low_precision <- result$CI_low_precision[upper.tri(result$CI_low_precision)]
      CI_high_precision <- result$CI_high_precision[upper.tri(result$CI_high_precision)]
      z_score_precision <- result$z_score_precision[upper.tri(result$z_score_precision)]
      p_precision <- result$p_precision[upper.tri(result$p_precision)]
      if(global == FALSE){
        corr_table <- cbind(pair, precision, partialCor, CI_low_precision, CI_high_precision, z_score_precision, p_precision)
        colnames(corr_table) <- c("gene1","gene2","precision","partialCor","CI_low_precision","CI_high_precision","z_score_precision","p_precision")
        if(!is.null(colnames(x))){
          gene_symbol <- colnames(x)
          corr_table[,1] <- gene_symbol[as.numeric(corr_table[,1])]
          corr_table[,2] <- gene_symbol[as.numeric(corr_table[,2])]
        }
        cytoscape_format_table$cytoscape_format_table <- corr_table
        if(csv_save == TRUE){
          cat("Save the Cytoscape format table. \n")
          if(is.null(directory)){
            directory <- tempdir()
            cat("The directory is not specified, so we use a temporary directory from the R session. \n")
            cat("The directory is ", directory, ".\n",sep = "")
          }
          else{
            cat("The directory is specified and it is ", directory, ".\n",sep = "")
          }
          write.csv(corr_table,file=paste(directory, "Cytoscape_D-S_GL.csv", sep = "/"),row.names = FALSE)
        }
        return(cytoscape_format_table)
      }
      else{
        if(is.null(true_graph)){
          cytoscape_format_table$threshold <- result$threshold
          cytoscape_format_table$FDR <- result$FDR
          global_decision <- NULL
          if(is.null(alpha)){
            alpha_length <- 2
          }
          else{
            alpha_length <- length(alpha)
          }
          for(i in 1:alpha_length){
            global_decision1 <- result$global_decision[[i]][upper.tri(result$global_decision[[i]])]
            global_decision <- cbind(global_decision, global_decision1)
          }
          if(is.null(alpha)){
            colnames(global_decision) <- c("global_decision_0.05","global_decision_0.1")
          }
          else{
            for(i in 1:alpha_length){
              colnames(global_decision)[i] <- paste("global_decision",alpha[i],sep="_")
            }
          }
          corr_table <- cbind(pair, precision, partialCor, CI_low_precision, CI_high_precision, z_score_precision, p_precision, global_decision)
          colnames(corr_table) <- c("gene1","gene2","precision","partialCor","CI_low_precision","CI_high_precision","z_score_precision","p_precision",colnames(global_decision))
          if(!is.null(colnames(x))){
            gene_symbol <- colnames(x)
            corr_table[,1] <- gene_symbol[as.numeric(corr_table[,1])]
            corr_table[,2] <- gene_symbol[as.numeric(corr_table[,2])]
          }
          cytoscape_format_table$cytoscape_format_table <- corr_table
          if(csv_save == TRUE){
            cat("Save the Cytoscape format table. \n")
            if(is.null(directory)){
              directory <- tempdir()
              cat("The directory is not specified, so we use a temporary directory from the R session. \n")
              cat("The directory is ", directory, ".\n",sep = "")
            }
            else{
              cat("The directory is specified and it is ", directory, ".\n",sep = "")
            }
            write.csv(corr_table,file=paste(directory, "Cytoscape_D-S_GL.csv", sep = "/"),row.names = FALSE)
          }
          return(cytoscape_format_table)
        }
        else{
          cytoscape_format_table$threshold <- result$threshold
          cytoscape_format_table$FDR <- result$FDR
          cytoscape_format_table$power <- result$power
          global_decision <- NULL
          if(is.null(alpha)){
            alpha_length <- 2
          }
          else{
            alpha_length <- length(alpha)
          }
          for(i in 1:alpha_length){
            global_decision1 <- result$global_decision[[i]][upper.tri(result$global_decision[[i]])]
            global_decision <- cbind(global_decision, global_decision1)
          }
          if(is.null(alpha)){
            colnames(global_decision) <- c("global_decision_0.05","global_decision_0.1")
          }
          else{
            for(i in 1:alpha_length){
              colnames(global_decision)[i] <- paste("global_decision",alpha[i],sep="_")
            }
          }
          corr_table <- cbind(pair, precision, partialCor, CI_low_precision, CI_high_precision, z_score_precision, p_precision, global_decision)
          colnames(corr_table) <- c("gene1","gene2","precision","partialCor","CI_low_precision","CI_high_precision","z_score_precision","p_precision",colnames(global_decision))
          if(!is.null(colnames(x))){
            gene_symbol <- colnames(x)
            corr_table[,1] <- gene_symbol[as.numeric(corr_table[,1])]
            corr_table[,2] <- gene_symbol[as.numeric(corr_table[,2])]
          }
          cytoscape_format_table$cytoscape_format_table <- corr_table
          if(csv_save == TRUE){
            cat("Save the Cytoscape format table. \n")
            if(is.null(directory)){
              directory <- tempdir()
              cat("The directory is not specified, so we use a temporary directory from the R session. \n")
              cat("The directory is ", directory, ".\n",sep = "")
            }
            else{
              cat("The directory is specified and it is ", directory, ".\n",sep = "")
            }
            write.csv(corr_table,file=paste(directory, "Cytoscape_D-S_GL.csv", sep = "/"),row.names = FALSE)
          }
          return(cytoscape_format_table)
        }
      }
    }

    if(method == "GFC_SL"){
      pair <- melt(upper.tri(result$T_stat))
      pair <- pair[pair[,3] %in% TRUE,c(1,2)]
      row.names(pair) <- NULL
      T_stat <- result$T_stat[upper.tri(result$T_stat)]
      if(is.null(true_graph)){
        cytoscape_format_table$threshold <- result$threshold
        cytoscape_format_table$FDR <- result$FDR
        global_decision <- NULL
        if(is.null(alpha)){
          alpha_length <- 2
        }
        else{
          alpha_length <- length(alpha)
        }
        for(i in 1:alpha_length){
          global_decision1 <- result$global_decision[[i]][upper.tri(result$global_decision[[i]])]
          global_decision <- cbind(global_decision, global_decision1)
        }
        if(is.null(alpha)){
          colnames(global_decision) <- c("global_decision_0.05","global_decision_0.1")
        }
        else{
          for(i in 1:alpha_length){
            colnames(global_decision)[i] <- paste("global_decision",alpha[i],sep="_")
          }
        }
        corr_table <- cbind(pair, T_stat, global_decision)
        colnames(corr_table) <- c("gene1","gene2","test_statistic",colnames(global_decision))
        if(!is.null(colnames(x))){
          gene_symbol <- colnames(x)
          corr_table[,1] <- gene_symbol[as.numeric(corr_table[,1])]
          corr_table[,2] <- gene_symbol[as.numeric(corr_table[,2])]
        }
        cytoscape_format_table$cytoscape_format_table <- corr_table
        if(csv_save == TRUE){
          cat("Save the Cytoscape format table. \n")
          if(is.null(directory)){
            directory <- tempdir()
            cat("The directory is not specified, so we use a temporary directory from the R session. \n")
            cat("The directory is ", directory, ".\n",sep = "")
          }
          else{
            cat("The directory is specified and it is ", directory, ".\n",sep = "")
          }
          write.csv(corr_table,file=paste(directory, "Cytoscape_GFC_SL.csv", sep = "/"),row.names = FALSE)
        }
        return(cytoscape_format_table)
      }
      else{
        cytoscape_format_table$threshold <- result$threshold
        cytoscape_format_table$FDR <- result$FDR
        cytoscape_format_table$power <- result$power
        global_decision <- NULL
        if(is.null(alpha)){
          alpha_length <- 2
        }
        else{
          alpha_length <- length(alpha)
        }
        for(i in 1:alpha_length){
          global_decision1 <- result$global_decision[[i]][upper.tri(result$global_decision[[i]])]
          global_decision <- cbind(global_decision, global_decision1)
        }
        if(is.null(alpha)){
          colnames(global_decision) <- c("global_decision_0.05","global_decision_0.1")
        }
        else{
          for(i in 1:alpha_length){
            colnames(global_decision)[i] <- paste("global_decision",alpha[i],sep="_")
          }
        }
        corr_table <- cbind(pair, T_stat, global_decision)
        colnames(corr_table) <- c("gene1","gene2","test_statistic",colnames(global_decision))
        if(!is.null(colnames(x))){
          gene_symbol <- colnames(x)
          corr_table[,1] <- gene_symbol[as.numeric(corr_table[,1])]
          corr_table[,2] <- gene_symbol[as.numeric(corr_table[,2])]
        }
        cytoscape_format_table$cytoscape_format_table <- corr_table
        if(csv_save == TRUE){
          cat("Save the Cytoscape format table. \n")
          if(is.null(directory)){
            directory <- tempdir()
            cat("The directory is not specified, so we use a temporary directory from the R session. \n")
            cat("The directory is ", directory, ".\n",sep = "")
          }
          else{
            cat("The directory is specified and it is ", directory, ".\n",sep = "")
          }
          write.csv(corr_table,file=paste(directory, "Cytoscape_GFC_SL.csv", sep = "/"),row.names = FALSE)
        }
        return(cytoscape_format_table)
      }
    }

    if(method == "GFC_L"){
      pair <- melt(upper.tri(result$T_stat))
      pair <- pair[pair[,3] %in% TRUE,c(1,2)]
      row.names(pair) <- NULL
      T_stat <- result$T_stat[upper.tri(result$T_stat)]
      if(is.null(true_graph)){
        cytoscape_format_table$threshold <- result$threshold
        cytoscape_format_table$FDR <- result$FDR
        global_decision <- NULL
        if(is.null(alpha)){
          alpha_length <- 2
        }
        else{
          alpha_length <- length(alpha)
        }
        for(i in 1:alpha_length){
          global_decision1 <- result$global_decision[[i]][upper.tri(result$global_decision[[i]])]
          global_decision <- cbind(global_decision, global_decision1)
        }
        if(is.null(alpha)){
          colnames(global_decision) <- c("global_decision_0.05","global_decision_0.1")
        }
        else{
          for(i in 1:alpha_length){
            colnames(global_decision)[i] <- paste("global_decision",alpha[i],sep="_")
          }
        }
        corr_table <- cbind(pair, T_stat, global_decision)
        colnames(corr_table) <- c("gene1","gene2","test_statistic",colnames(global_decision))
        if(!is.null(colnames(x))){
          gene_symbol <- colnames(x)
          corr_table[,1] <- gene_symbol[as.numeric(corr_table[,1])]
          corr_table[,2] <- gene_symbol[as.numeric(corr_table[,2])]
        }
        cytoscape_format_table$cytoscape_format_table <- corr_table
        if(csv_save == TRUE){
          cat("Save the Cytoscape format table. \n")
          if(is.null(directory)){
            directory <- tempdir()
            cat("The directory is not specified, so we use a temporary directory from the R session. \n")
            cat("The directory is ", directory, ".\n",sep = "")
          }
          else{
            cat("The directory is specified and it is ", directory, ".\n",sep = "")
          }
          write.csv(corr_table,file=paste(directory, "Cytoscape_GFC_L.csv", sep = "/"),row.names = FALSE)
        }
        return(cytoscape_format_table)
      }
      else{
        cytoscape_format_table$threshold <- result$threshold
        cytoscape_format_table$FDR <- result$FDR
        cytoscape_format_table$power <- result$power
        global_decision <- NULL
        if(is.null(alpha)){
          alpha_length <- 2
        }
        else{
          alpha_length <- length(alpha)
        }
        for(i in 1:alpha_length){
          global_decision1 <- result$global_decision[[i]][upper.tri(result$global_decision[[i]])]
          global_decision <- cbind(global_decision, global_decision1)
        }
        if(is.null(alpha)){
          colnames(global_decision) <- c("global_decision_0.05","global_decision_0.1")
        }
        else{
          for(i in 1:alpha_length){
            colnames(global_decision)[i] <- paste("global_decision",alpha[i],sep="_")
          }
        }
        corr_table <- cbind(pair, T_stat, global_decision)
        colnames(corr_table) <- c("gene1","gene2","test_statistic",colnames(global_decision))
        if(!is.null(colnames(x))){
          gene_symbol <- colnames(x)
          corr_table[,1] <- gene_symbol[as.numeric(corr_table[,1])]
          corr_table[,2] <- gene_symbol[as.numeric(corr_table[,2])]
        }
        cytoscape_format_table$cytoscape_format_table <- corr_table
        if(csv_save == TRUE){
          cat("Save the Cytoscape format table. \n")
          if(is.null(directory)){
            directory <- tempdir()
            cat("The directory is not specified, so we use a temporary directory from the R session. \n")
            cat("The directory is ", directory, ".\n",sep = "")
          }
          else{
            cat("The directory is specified and it is ", directory, ".\n",sep = "")
          }
          write.csv(corr_table,file=paste(directory, "Cytoscape_GFC_L.csv", sep = "/"),row.names = FALSE)
        }
        return(cytoscape_format_table)
      }
    }

    if(method == "B_NW_SL"){
      pair <- melt(upper.tri(result$precision))
      pair <- pair[pair[,3] %in% TRUE,c(1,2)]
      row.names(pair) <- NULL
      precision <- result$precision[upper.tri(result$precision)]
      partialCor <- result$partialCor[upper.tri(result$partialCor)]
      CI_low_precision <- result$CI_low_precision[upper.tri(result$CI_low_precision)]
      CI_high_precision <- result$CI_high_precision[upper.tri(result$CI_high_precision)]
      CI_low_partialCor <- result$CI_low_partialCor[upper.tri(result$CI_low_partialCor)]
      CI_high_partialCor <- result$CI_high_partialCor[upper.tri(result$CI_high_partialCor)]
      z_score_precision <- result$z_score_precision[upper.tri(result$z_score_precision)]
      z_score_partialCor <- result$z_score_partialCor[upper.tri(result$z_score_partialCor)]
      p_precision <- result$p_precision[upper.tri(result$p_precision)]
      p_partialCor <- result$p_partialCor[upper.tri(result$p_partialCor)]
      if(global == FALSE){
        corr_table <- cbind(pair, precision, partialCor, CI_low_precision, CI_high_precision, CI_low_partialCor, CI_high_partialCor, z_score_precision, z_score_partialCor, p_precision, p_partialCor)
        colnames(corr_table) <- c("gene1","gene2","precision","partialCor","CI_low_precision","CI_high_precision", "CI_low_partialCor", "CI_high_partialCor", "z_score_precision", "z_score_partialCor", "p_precision", "p_partialCor")
        if(!is.null(colnames(x))){
          gene_symbol <- colnames(x)
          corr_table[,1] <- gene_symbol[as.numeric(corr_table[,1])]
          corr_table[,2] <- gene_symbol[as.numeric(corr_table[,2])]
        }
        cytoscape_format_table$cytoscape_format_table <- corr_table
        if(csv_save == TRUE){
          cat("Save the Cytoscape format table. \n")
          if(is.null(directory)){
            directory <- tempdir()
            cat("The directory is not specified, so we use a temporary directory from the R session. \n")
            cat("The directory is ", directory, ".\n",sep = "")
          }
          else{
            cat("The directory is specified and it is ", directory, ".\n",sep = "")
          }
          write.csv(corr_table,file=paste(directory, "Cytoscape_B_NW_SL.csv", sep = "/"),row.names = FALSE)
        }
        return(cytoscape_format_table)
      }
      else{
        if(is.null(true_graph)){
          cytoscape_format_table$threshold <- result$threshold
          cytoscape_format_table$FDR <- result$FDR
          global_decision <- NULL
          if(is.null(alpha)){
            alpha_length <- 2
          }
          else{
            alpha_length <- length(alpha)
          }
          for(i in 1:alpha_length){
            global_decision1 <- result$global_decision[[i]][upper.tri(result$global_decision[[i]])]
            global_decision <- cbind(global_decision, global_decision1)
          }
          if(is.null(alpha)){
            colnames(global_decision) <- c("global_decision_0.05","global_decision_0.1")
          }
          else{
            for(i in 1:alpha_length){
              colnames(global_decision)[i] <- paste("global_decision",alpha[i],sep="_")
            }
          }
          corr_table <- cbind(pair, precision, partialCor, CI_low_precision, CI_high_precision, CI_low_partialCor, CI_high_partialCor, z_score_precision, z_score_partialCor, p_precision, p_partialCor, global_decision)
          colnames(corr_table) <- c("gene1","gene2","precision","partialCor","CI_low_precision","CI_high_precision", "CI_low_partialCor", "CI_high_partialCor", "z_score_precision", "z_score_partialCor", "p_precision", "p_partialCor", colnames(global_decision))
          if(!is.null(colnames(x))){
            gene_symbol <- colnames(x)
            corr_table[,1] <- gene_symbol[as.numeric(corr_table[,1])]
            corr_table[,2] <- gene_symbol[as.numeric(corr_table[,2])]
          }
          cytoscape_format_table$cytoscape_format_table <- corr_table
          if(csv_save == TRUE){
            cat("Save the Cytoscape format table. \n")
            if(is.null(directory)){
              directory <- tempdir()
              cat("The directory is not specified, so we use a temporary directory from the R session. \n")
              cat("The directory is ", directory, ".\n",sep = "")
            }
            else{
              cat("The directory is specified and it is ", directory, ".\n",sep = "")
            }
            write.csv(corr_table,file=paste(directory, "Cytoscape_B_NW_SL.csv", sep = "/"),row.names = FALSE)
          }
          return(cytoscape_format_table)
        }
        else{
          cytoscape_format_table$threshold <- result$threshold
          cytoscape_format_table$FDR <- result$FDR
          cytoscape_format_table$power <- result$power
          global_decision <- NULL
          if(is.null(alpha)){
            alpha_length <- 2
          }
          else{
            alpha_length <- length(alpha)
          }
          for(i in 1:alpha_length){
            global_decision1 <- result$global_decision[[i]][upper.tri(result$global_decision[[i]])]
            global_decision <- cbind(global_decision, global_decision1)
          }
          if(is.null(alpha)){
            colnames(global_decision) <- c("global_decision_0.05","global_decision_0.1")
          }
          else{
            for(i in 1:alpha_length){
              colnames(global_decision)[i] <- paste("global_decision",alpha[i],sep="_")
            }
          }
          corr_table <- cbind(pair, precision, partialCor, CI_low_precision, CI_high_precision, CI_low_partialCor, CI_high_partialCor, z_score_precision, z_score_partialCor, p_precision, p_partialCor, global_decision)
          colnames(corr_table) <- c("gene1","gene2","precision","partialCor","CI_low_precision","CI_high_precision", "CI_low_partialCor", "CI_high_partialCor", "z_score_precision", "z_score_partialCor", "p_precision", "p_partialCor", colnames(global_decision))
          if(!is.null(colnames(x))){
            gene_symbol <- colnames(x)
            corr_table[,1] <- gene_symbol[as.numeric(corr_table[,1])]
            corr_table[,2] <- gene_symbol[as.numeric(corr_table[,2])]
          }
          cytoscape_format_table$cytoscape_format_table <- corr_table
          if(csv_save == TRUE){
            cat("Save the Cytoscape format table. \n")
            if(is.null(directory)){
              directory <- tempdir()
              cat("The directory is not specified, so we use a temporary directory from the R session. \n")
              cat("The directory is ", directory, ".\n",sep = "")
            }
            else{
              cat("The directory is specified and it is ", directory, ".\n",sep = "")
            }
            write.csv(corr_table,file=paste(directory, "Cytoscape_B_NW_SL.csv", sep = "/"),row.names = FALSE)
          }
          return(cytoscape_format_table)
        }
      }
    }

  }
}

