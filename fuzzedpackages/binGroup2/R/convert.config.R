# Start convert.config() function
###############################################################################
# Brianna Hitt
# Support function for OTC functions for one and two diseases
# converts configurations to a usable format for OTC functions

# updated convert.config() allows two disease functions
# updated 04-08-18 to allow for two disease assays
# function to convert the configurations to a single column

convert.config <- function(algorithm, results, diseases=1, old.label="pool.sz", 
                           new.label="pool.szs", sep=","){
  new.results <- cbind(results, NA)
  colnames(new.results)[dim(new.results)[2]] <- new.label
  
  index.pools <- which(colnames(new.results)==old.label)[1]
  
  if(algorithm %in% c("D2")){
    if(diseases==1){
      final <- data.frame(I=as.numeric(new.results[,which(colnames(new.results)=="I")]),
                          ET=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="ET")]), 4), nsmall=4)),
                          value=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="value")]), 4), nsmall=4)),
                          PSe=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="PSe")]), 4), nsmall=4)),
                          PSp=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="PSp")]), 4), nsmall=4)),
                          PPPV=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="PPPV")]), 4), nsmall=4)),
                          PNPV=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="PNPV")]), 4), nsmall=4)),
                          row.names=NULL)
    } else if(diseases==2){
      final <- data.frame(I=as.numeric(new.results[,which(colnames(new.results)=="I")]),
                          ET=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="ET")]), 4), nsmall=4)),
                          value=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="value")]), 4), nsmall=4)),
                          PSe1=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="PSe1")]), 4), nsmall=4)),
                          PSp1=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="PSp1")]), 4), nsmall=4)),
                          PPPV1=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="PPPV1")]), 4), nsmall=4)),
                          PNPV1=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="PNPV1")]), 4), nsmall=4)),
                          PSe2=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="PSe2")]), 4), nsmall=4)),
                          PSp2=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="PSp2")]), 4), nsmall=4)),
                          PPPV2=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="PPPV2")]), 4), nsmall=4)),
                          PNPV2=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="PNPV2")]), 4), nsmall=4)),
                          row.names=NULL)
    }
  } else if(algorithm %in% c("A2", "IA2", "A2M")){
    if(diseases==1){
      final <- data.frame(I=as.numeric(new.results[,which(colnames(new.results)=="I")]),
                          N=as.numeric(new.results[,which(colnames(new.results)=="N")]),
                          ET=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="ET")]), 4), nsmall=4)),
                          value=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="value")]), 4), nsmall=4)),
                          PSe=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="PSe")]), 4), nsmall=4)),
                          PSp=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="PSp")]), 4), nsmall=4)),
                          PPPV=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="PPPV")]), 4), nsmall=4)),
                          PNPV=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="PNPV")]), 4), nsmall=4)),
                          row.names=NULL)
    } else if(diseases==2){
      final <- data.frame(I=as.numeric(new.results[,which(colnames(new.results)=="I")]),
                          N=as.numeric(new.results[,which(colnames(new.results)=="N")]),
                          ET=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="ET")]), 4), nsmall=4)),
                          value=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="value")]), 4), nsmall=4)),
                          PSe1=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="PSe1")]), 4), nsmall=4)),
                          PSp1=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="PSp1")]), 4), nsmall=4)),
                          PPPV1=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="PPPV1")]), 4), nsmall=4)),
                          PNPV1=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="PNPV1")]), 4), nsmall=4)),
                          PSe2=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="PSe2")]), 4), nsmall=4)),
                          PSp2=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="PSp2")]), 4), nsmall=4)),
                          PPPV2=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="PPPV2")]), 4), nsmall=4)),
                          PNPV2=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="PNPV2")]), 4), nsmall=4)),
                          row.names=NULL)
    }
  } else if(algorithm %in% c("ID2")){
    for(i in 1:dim(new.results)[1]){
      
      config <- paste(new.results[i,index.pools], new.results[i,(index.pools+1)], 
                      sep=sep)
      for(j in (index.pools+2):(dim(new.results)[2]-1)){
        if(new.results[i,j]==0){
          config <- config
        } else{
          config <- paste(config, new.results[i,j], sep=sep)
        }
      }
      new.results[i,dim(new.results)[2]] <- config
    }
    if(diseases==1){
      final <- data.frame(N=as.numeric(new.results[,which(colnames(new.results)=="N")]),
                          config=new.results[,which(colnames(new.results)==new.label)],
                          ET=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="ET")]), 4), nsmall=4)),
                          value=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="value")]), 4), nsmall=4)),
                          PSe=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="PSe")]), 4), nsmall=4)),
                          PSp=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="PSp")]), 4), nsmall=4)),
                          PPPV=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="PPPV")]), 4), nsmall=4)),
                          PNPV=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="PNPV")]), 4), nsmall=4)),
                          row.names=NULL)
    } else if(diseases==2){
      final <- data.frame(N=as.numeric(new.results[,which(colnames(new.results)=="N")]),
                          config=new.results[,which(colnames(new.results)==new.label)],
                          ET=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="ET")]), 4), nsmall=4)),
                          value=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="value")]), 4), nsmall=4)),
                          PSe1=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="PSe1")]), 4), nsmall=4)),
                          PSp1=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="PSp1")]), 4), nsmall=4)),
                          PPPV1=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="PPPV1")]), 4), nsmall=4)),
                          PNPV1=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="PNPV1")]), 4), nsmall=4)),
                          PSe2=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="PSe2")]), 4), nsmall=4)),
                          PSp2=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="PSp2")]), 4), nsmall=4)),
                          PPPV2=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="PPPV2")]), 4), nsmall=4)),
                          PNPV2=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="PNPV2")]), 4), nsmall=4)),
                          row.names=NULL)
    }
  } else if(algorithm %in% c("D3", "ID3")){
    for(i in 1:dim(new.results)[1]){
      
      config <- paste(new.results[i,index.pools], new.results[i,(index.pools+1)], 
                      sep=sep)
      for(j in (index.pools+2):(dim(new.results)[2]-1)){
        if(new.results[i,j]==0){
          config <- config
        } else{
          config <- paste(config, new.results[i,j], sep=sep)
        }
      }
      new.results[i,dim(new.results)[2]] <- config
    }
    if(diseases==1){
      final <- data.frame(I=as.numeric(new.results[,which(colnames(new.results)=="I")]),
                          config=new.results[,which(colnames(new.results)==new.label)],
                          ET=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="ET")]), 4), nsmall=4)),
                          value=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="value")]), 4), nsmall=4)),
                          PSe=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="PSe")]), 4), nsmall=4)),
                          PSp=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="PSp")]), 4), nsmall=4)),
                          PPPV=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="PPPV")]), 4), nsmall=4)),
                          PNPV=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="PNPV")]), 4), nsmall=4)),
                          row.names=NULL)
    } else if(diseases==2){
      final <- data.frame(I=as.numeric(new.results[,which(colnames(new.results)=="I")]),
                          config=new.results[,which(colnames(new.results)==new.label)],
                          ET=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="ET")]), 4), nsmall=4)),
                          value=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="value")]), 4), nsmall=4)),
                          PSe1=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="PSe1")]), 4), nsmall=4)),
                          PSp1=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="PSp1")]), 4), nsmall=4)),
                          PPPV1=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="PPPV1")]), 4), nsmall=4)),
                          PNPV1=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="PNPV1")]), 4), nsmall=4)),
                          PSe2=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="PSe2")]), 4), nsmall=4)),
                          PSp2=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="PSp2")]), 4), nsmall=4)),
                          PPPV2=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="PPPV2")]), 4), nsmall=4)),
                          PNPV2=as.numeric(format(round(as.numeric(new.results[,which(colnames(new.results)=="PNPV2")]), 4), nsmall=4)),
                          row.names=NULL)
    }
  }
  final
}

#