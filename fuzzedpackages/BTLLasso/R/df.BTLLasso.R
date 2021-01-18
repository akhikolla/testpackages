df.BTLLasso <- function(coefs, design, m){

  df <- c()

  for(l in 1:nrow(coefs)){
    df.l <- design$n.theta
    coefs.l <- coefs[l,]
    start <- design$n.theta+1
    
    if(design$n.order>0){
      end <- start+design$n.order-1
      xhelp <- coefs.l[start:end]
      df.l <- df.l + length(unique(xhelp[xhelp!=0]))
      start <- end+1
    }
    
    if(design$n.intercepts>0){
      end <- start+design$n.intercepts
      xhelp <- coefs.l[start:end]
      df.l <- df.l + length(unique(xhelp)) - 1
      start <- end+1
    }
    
    if(design$p.X>0){
      for(ll in 1:design$p.X){
        end <- start+m-1
        xhelp <- coefs.l[start:end]
        df.l <- df.l + length(unique(xhelp)) - 1
        start <- end+1
      }
      
    }
      
      if(design$p.Z1>0){
        for(ll in 1:design$p.Z1){
          end <- start+m-1
          xhelp <- coefs.l[start:end]
          df.l <- df.l + length(unique(xhelp[xhelp!=0]))
          start <- end+1
        }
      }
      
    if(design$p.Z2>0){
        end <- start+design$p.Z2-1
        xhelp <- coefs.l[start:end]
        df.l <- df.l + length(unique(xhelp[xhelp!=0]))
        start <- end+1
    }
    df[l] <- df.l
  }
  df
}