calculate_p_value <- function(mu,std,val)
{
  if (val<=0 && std<=0 && mu<=0)
  {
    z <- 0;
  }
  else
  {
    z <- ((val-mu)/std);
  }
  p_val <- 2*pnorm(-abs(z));
  return(p_val)
}

#Fast implementation removing the node with least contribution to CGHD
differential_subnetwork_analysis_closedform <- function(ghd_val,mu_perm,p,matrixA,matrixB,threshold){
  df <- NULL
  listids <- NULL
  original_rownames <- row.names(matrixA)
  count=1
  i = 0;
  while(p>4){
    delta_prev <- ghd_val-mu_perm;
    if (is.null(listids)){
        newoutputlist <- NULL
        N <- nrow(matrixA)
        global_A <- (1/(N*(N-1)))*(sum(matrixA)-sum(diag(matrixA)))
        global_B <- (1/(N*(N-1)))*(sum(matrixB)-sum(diag(matrixB)))
        A <- matrixA - global_A
        B <- matrixB - global_B
        
        contri_rem_ghd_node <- ghd_val - (1/(N*(N-1)))*(colSums(A^2) - diag(A^2) + colSums(B^2) - diag(B^2) - 2*colSums(A*B) + 2*diag(A*B));
        contri_rem_mu <- mu_perm -((1/(N*(N-1)))*(colSums(matrixA^2) - diag(matrixA^2) + colSums(matrixB^2) - diag(matrixB^2))
                             - (2/((N^2)*((N-1)^2)))*((colSums(matrixA)-diag(matrixA))*sum(matrixB)+(colSums(matrixB)-diag(matrixB))*sum(matrixA) - (colSums(matrixA)-diag(matrixA))*(colSums(matrixB)-diag(matrixB))))
        
        rm(A)
        rm(B)
        gc()
        outputlist <- (ghd_val - contri_rem_ghd_node)-(mu_perm-contri_rem_mu);
        
    }
    else if (p_val>threshold) {
      a1 <- matrixA[-listids,-listids]
      b1 <- matrixB[-listids,-listids]
      newoutputlist <- NULL;
        newoutputlist <- foreach (i=1:p, .packages='Matrix', .combine=rbind, .inorder = TRUE) %dopar% 
        {
          gv_x <- GHD_Fast(a1[-i,-i],b1[-i,-i])
          mu_x <- MU_Fast(a1[-i,-i],b1[-i,-i])
          delta_x <- gv_x -mu_x
          tupele <- cbind(i,delta_x)
        }
      rm(a1)
      rm(b1)
      gc()
    }
    if (is.null(listids))
    {
      vector_v <- order(outputlist,decreasing=FALSE);
      find_max_id <- vector_v[1];
      vector_v <- vector_v[-1];
      dim_name <- row.names(matrixA)[find_max_id]
    }
    else if (is.null(newoutputlist)){
      find_max_id <- vector_v[1]
      vector_v <- vector_v[-1];
      dim_name <- row.names(matrixA)[find_max_id]
    }
    else if (!is.null(newoutputlist))
    {
      #if (length(newoutputlist)>2){
      find_max_id <- newoutputlist[which.max(newoutputlist[,2]),1];
      find_max_id <- as.numeric(find_max_id);
      newoutputlist <- newoutputlist[-find_max_id,];
      #  id <- which(newoutputlist[,1]==find_max_id);
      #  newoutputlist <- newoutputlist[-id,];
      #}else{
      #  find_max_id <- newoutputlist[1];
      #  id <- 1;
      #  newoutputlist <- NULL;
      #}
      dim_name <- original_rownames[-listids][find_max_id] ;
    }
    actual_id <- which(original_rownames==dim_name)
    listids <- c(listids,actual_id)
    p <- p-1;
    
    #Calculate p-value for remaining nodes
    ghd_val <- GHD_Fast(matrixA[-listids,-listids],matrixB[-listids,-listids]);
    mu_perm <- MU_Fast(matrixA[-listids,-listids],matrixB[-listids,-listids])
    std_perm <- STD_Fast(matrixA[-listids,-listids],matrixB[-listids,-listids])
    std_perm <- sqrt(abs(std_perm)/((p^3)*(p-1)^3))
    p_val <- calculate_p_value(mu_perm,std_perm,ghd_val);
    #Contribution per iteration
    #print(paste(count,threshold,dim_name,p_val,ghd_val,mu_perm,std_perm))
    count=count+1
    df <- rbind(df,c(actual_id,dim_name,p_val,ghd_val,mu_perm,std_perm))
    gc()
  }
  p_val_list <- df[,3];
  V7 <- p.adjust(p_val_list,method="holm");
  df <- cbind(df,V7);
  df <- as.data.frame(df)
  return(df)
}


#Parallel Implementation for Original 
differential_subnetwork_analysis_original <- function(ghd_val,mu_perm,p,matrixA,matrixB,threshold){
  df <- NULL
  listids <- NULL
  original_rownames <- row.names(matrixA)
  count=1
  i <- 0;
  while(p>4){
    delta_prev <- ghd_val-mu_perm;
    if (is.null(listids)){
        outputlist <- foreach (i=1:p, .packages='Matrix', .combine=rbind, .inorder=TRUE) %dopar% 
        {
          gv_x <- GHD_Fast(matrixA[-i,-i],matrixB[-i,-i])
          mu_x <- MU_Fast(matrixA[-i,-i],matrixB[-i,-i])
          delta_x <- gv_x -mu_x
          tupele <- cbind(i,delta_x)
        }
    }
    else{
      a1 <- matrixA[-listids,-listids]
      b1 <- matrixB[-listids,-listids]
        outputlist <- foreach (i=1:p, .packages='Matrix', .combine=rbind, .inorder = TRUE) %dopar% 
        {
          gv_x <- GHD_Fast(a1[-i,-i],b1[-i,-i])
          mu_x <- MU_Fast(a1[-i,-i],b1[-i,-i])
          delta_x <- gv_x -mu_x
          tupele <- cbind(i,delta_x)
        }
      rm(a1)
      rm(b1)
      gc()
    }
    find_max_id <- outputlist[which.max(outputlist[,2]),1]
    if (is.null(listids))
    {
      dim_name <- row.names(matrixA)[find_max_id]
    }
    else{
      dim_name <- row.names(matrixA[-listids,-listids])[find_max_id]
    }
    actual_id <- which(original_rownames==dim_name)
    listids <- c(listids,actual_id)
    p <- p-1;
    
    #Calculate p-value for remaining nodes
    small_delta_max <- max(outputlist)
    ghd_val <- GHD_Fast(matrixA[-listids,-listids],matrixB[-listids,-listids]);
    mu_perm <- MU_Fast(matrixA[-listids,-listids],matrixB[-listids,-listids]);
    std_perm <- STD_Fast(matrixA[-listids,-listids],matrixB[-listids,-listids])
    std_perm <- sqrt(abs(std_perm)/((p^3)*(p-1)^3))
    p_val <- calculate_p_value(mu_perm,std_perm,ghd_val);
    #To keep track of iteration number
    #print(paste(count,dim_name,p_val,ghd_val,mu_perm,std_perm))
    count=count+1
    df <- rbind(df,c(actual_id,dim_name,p_val,ghd_val,mu_perm,std_perm))
    gc()
  }
  p_val_list <- df[,3];
  p_val_new_list <- p.adjust(p_val_list,method="holm");
  V7 <- p_val_new_list;
  df <- cbind(df,V7);
  df <- as.data.frame(df)
  return(df)
}


#Fast implementation of Fast Approximation
differential_subnetwork_analysis_fastapprox <- function(ghd_val,mu_perm,p,matrixA,matrixB,threshold){
  df <- NULL
  listids <- NULL
  original_rownames <- row.names(matrixA)
  count=1
  i <- 0;
  p_val <- 0;
  while(p>4){
    delta_prev <- ghd_val-mu_perm;
    if (is.null(listids)){
      newoutputlist <- NULL
        outputlist <- foreach (i=1:p, .packages='Matrix', .combine=rbind, .inorder = TRUE) %dopar% 
        {
          gv_x <- GHD_Fast(matrixA[-i,-i],matrixB[-i,-i])
          mu_x <- MU_Fast(matrixA[-i,-i],matrixB[-i,-i])
          delta_x <- gv_x - mu_x
          tupele <- cbind(i,delta_x)
        }
    }
    else if (p_val>threshold){
      a1 <- matrixA[-listids,-listids]
      b1 <- matrixB[-listids,-listids]
          newoutputlist <- foreach (i=1:p, .packages='Matrix', .combine=rbind, .inorder = TRUE) %dopar% 
          {
            gv_x <- GHD_Fast(a1[-i,-i],b1[-i,-i])
            mu_x <- MU_Fast(a1[-i,-i],b1[-i,-i])
            delta_x <- (gv_x-mu_x);
            tupele <- cbind(i,delta_x)
          }
      rm(a1)
      rm(b1)
      gc()
    }
    if (is.null(listids))
    {
      vector_v <- order(outputlist[,2],decreasing=TRUE);
      find_max_id <- vector_v[1];
      vector_v <- vector_v[-1];
      dim_name <- row.names(matrixA)[find_max_id]
    }
    else if (is.null(newoutputlist)){
      find_max_id <- vector_v[1]
      vector_v <- vector_v[-1];
      dim_name <- row.names(matrixA)[find_max_id]
    }
    else if (!is.null(newoutputlist))
    {
        find_max_id <- newoutputlist[which.max(newoutputlist[,2]),1];
        newoutputlist <- newoutputlist[-find_max_id,]
        dim_name <- original_rownames[-listids][find_max_id] ;
    }
    actual_id <- which(original_rownames==dim_name)
    listids <- c(listids,actual_id)
    p <- p-1;
    
    #Calculate p-value for remaining nodes
    ghd_val <- GHD_Fast(matrixA[-listids,-listids],matrixB[-listids,-listids]);
    mu_perm <- MU_Fast(matrixA[-listids,-listids],matrixB[-listids,-listids]);
    std_perm <- STD_Fast(matrixA[-listids,-listids],matrixB[-listids,-listids])
    std_perm <- sqrt(abs(std_perm)/((p^3)*(p-1)^3))
    gc()
    p_val <- calculate_p_value(mu_perm,std_perm,ghd_val);
    #To keep track of iteration number
    #print(paste(count,threshold,dim_name,p_val,ghd_val,mu_perm,std_perm))
    count=count+1
    df <- rbind(df,c(actual_id,dim_name,p_val,ghd_val,mu_perm,std_perm))
    gc()
  }
  p_val_list <- df[,3];
  p_val_new_list <- p.adjust(p_val_list,method="holm");
  V7 <- p_val_new_list;
  df <- cbind(df,V7);
  df <- as.data.frame(df)
  return(df)
}

diffnet <- function(g_A = sample_grg(6,0.15,torus=TRUE,coords=TRUE), g_B = permute(g_A,c(sample(5),6)),
                    p = 6, threshold = 1e-50, approach = "closed-form"){
  
  p_value_vector_approach <- rep(0,p)
  
  #Get the adjaceny matrices from 
  actual_nodes <- c(1:p)
  adjacencyA <- get.adjacency(g_A);
  adjacencyB <- get.adjacency(g_B);
  if (is.null(rownames(adjacencyA)))
  {
    rownames(adjacencyA) <- paste("N",c(1:p),sep="")
  }
  if (is.null(rownames(adjacencyB)))
  {
    rownames(adjacencyB) <- paste("N",c(1:p),sep="")
  }
    
  cosine_sim_A <- cosSparse(t(adjacencyA),Y=t(adjacencyA),norm = norm2);
  cosine_sim_B <- cosSparse(t(adjacencyB),Y=t(adjacencyB),norm = norm2)
    
  ghd_val <- GHD_Fast(cosine_sim_A,cosine_sim_B)
  mu_permutation <- MU_Fast(cosine_sim_A,cosine_sim_B)
    
  if (approach=="closed-form")
  {
    #Closed-Form approach
    df_approach <- differential_subnetwork_analysis_closedform(ghd_val,mu_permutation,p,cosine_sim_A,cosine_sim_B,threshold)
  }
  else if (approach == "original")
  {
    #Original dGHD approach
    df_approach <- differential_subnetwork_analysis_original(ghd_val,mu_permutation,p,cosine_sim_A, cosine_sim_B, threshold)
  }
  else{
    #Fast-approximation approach
    df_approach <- differential_subnetwork_analysis_fastapprox(ghd_val,mu_permutation,p,cosine_sim_A,cosine_sim_B, threshold)
  }
  
  #Perform Ordering for ROC Analysis
  pval <- as.numeric(as.character(df_approach[,3]))
  total_nodes <- as.integer(as.character(df_approach[,1]))
  nodes_not_list <- setdiff(actual_nodes,total_nodes)
  total_nodes <- c(total_nodes,nodes_not_list)
  outputpval <- c(pval,t(rep(1,length(nodes_not_list))))
  
  indices <- order(total_nodes)
  p_value_vector_approach <- outputpval[indices];
  
  return(p_value_vector_approach)
}
