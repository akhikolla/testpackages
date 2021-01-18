## This work is licensed under the Creative Commons
## Attribution-ShareAlike 3.0 Unported License. To view a copy of this
## license, visit http://creativecommons.org/licenses/by-sa/3.0/ or send
## a letter to Creative Commons, 171 Second Street, Suite 300, San
## Francisco, California, 94105, USA.
##
## (c) Copyright Barry Rowlingson 2010
## B.Rowlingson@lancaster.ac.uk

queue = function(statistics=FALSE){
###
### queue implmentation (c) Barry Rowlingson 2010
###
### creates a queue object. the 'statistics' flag adds basic logging
###
  e=new.env()
  q=list()
  assign("q",q,envir=e)
  assign("stats",NULL,envir=e)
  if(statistics){
    assign("stats",list(maxsize=0,on=0,off=0,minsize=0),envir=e)
  }
  class(e)=c("queue","environment")
  e
}

getStats = function(eq){
### get the statistics from the queue
  return(get("stats",eq))
}

resetStats=function(eq){
### reset the statistics from the queue
  q=get("q",envir=eq)
  m=list(maxsize=length(q),minsize=length(q),on=0,off=0)
  assign("stats",m,envir=eq)
  return(m)
}

### S3 Generic Functions
enqueue=function(eq,v){UseMethod("enqueue")}
dequeue=function(eq,v){UseMethod("dequeue")}

enqueue.queue=function(eq,v){
###
### S3 Method for adding a value 'v' to the queue
###
  ## add the value to the list
  q=c(v,get("q",envir=eq))
  ## stick the list back in the environment
  assign("q",q,envir=eq)

  ## process the statistics if they're there:
  if(length({m=getStats(eq)})>0){
    m$on=m$on+1
    if(length(q)>m$maxsize){
      m$maxsize=length(q)
    }
    assign("stats",m,envir=eq)
  }

  ## may as well return something, so here's the value back
  v
}


dequeue.queue=function(eq){
###
### S3 Method for taking something off the list
###
  ## get the queue, check if anything on it:
  q=get("q",envir=eq)
  if(length(q)==0){
    stop("Attempt to take element from empty queue")
  }

  ## take the last value
  v=q[[length(q)]]

  ## take the last value off:
  if(length(q)==1){
    assign("q",list(),eq)
  }else{
    assign("q",q[1:(length(q)-1)],eq)
  }
  ## update stats
  if(length({m=getStats(eq)})>0){
    m$off=m$off+1
    if(length(get("q",envir=eq))<m$minsize){
      m$minsize=length(q)-1
    }
    assign("stats",m,envir=eq)
  }
  return(v)
}

### print method
print.queue=function(x,...){
  print(get("q",envir=x))
}

summary.queue=function(object,...){
  cat("Queue object with ",length(get("q",envir=object))," elements\n")
  if(length({m=getStats(object)})>0){
    cat(m$on,"added,",m$off,"removed. min/max size=",m$minsize,"/",m$maxsize,"\n")
  }

}


