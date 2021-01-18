process.data<-function(data,indexed=FALSE,mad=FALSE)
{
    col.names<-colnames(data)
    if(!is.data.frame(data))
    {
        data<-data.frame(data)
    }
    if(!indexed && is.null(col.names))
    {
        data<-cbind(data.frame(seq(1:nrow(data))),data)
        colnames(data)<-c("k",paste("V.",1:(ncol(data)-1),sep=""))               
    }
    if(!indexed && !is.null(col.names))
    {
        data<-cbind(data.frame("k"=seq(1:nrow(data))),data)
    }
    if(indexed && is.null(col.names))
    {
        colnames(data)<-c("k",paste("V.",1:(ncol(data)-1),sep=""))               
    }
    # apply mad if required
    if(mad)
    {
        usr.names<-names(data)
        data<-cbind(data[,1],data.frame(Map(function(i) data[,i]*sqrt(2)/mad(diff(data[,i])),2:ncol(data))))
        names(data)<-usr.names
    }
    return(data)
}

