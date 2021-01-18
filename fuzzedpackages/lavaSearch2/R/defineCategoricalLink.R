### defineCategoricalLink.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 26 2017 (16:35) 
## Version: 
## last-updated: aug  6 2018 (15:32) 
##           By: Brice Ozenne
##     Update #: 156
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * documentation - defineCategoricalLink
#' @title Identify Categorical Links in LVM
#' @description Identify categorical links in latent variable models.
#' @name defineCategoricalLink
#' 
#' @param object a \code{lvm} model.
#' @param link [character] the links to be analyzed. If \code{NULL}, all the coefficients from the lvm model are used instead.
#' @param data [data.frame] the dataset that will be used to fit the model. If \code{NULL}, a simulated data will be generated from the model.
#'
#' @return a \code{data.frame} with a description of each link in rows. \cr
#' The column factitious identify the links that will be replaced with other links
#' (e.g. "Y1~X1" becomes "Y1~X1b" and "Y1~X1c"). \cr
#' 
#' @examples
#' \dontrun{
#' defineCategoricalLink <- lavaSearch2:::defineCategoricalLink
#' defineCategoricalLink.lvm <- lavaSearch2:::defineCategoricalLink.lvm
#' 
#' ## linear model
#' m <- lvm(Y1~X1+X2,Y2~X1+X3)
#' categorical(m, K = 3) <- "X1"
#' try(defineCategoricalLink(m)) # error
#'
#' categorical(m, K = 3, labels = 1:3) <- "X1"
#' defineCategoricalLink(m)
#' defineCategoricalLink(m, "Y~X1")
#' defineCategoricalLink(m, "X1:0|1")
#' defineCategoricalLink(m, "X1:1|2")
#' defineCategoricalLink(m, c("X1:0|1", "X1:1|2"))
#' defineCategoricalLink(m, c("Y~X1","Y~X2"))
#' defineCategoricalLink(m, c("Y~X2","Y~X1"))
#'
#' ## latent variable model
#' m <- lvm()
#' regression(m) <- c(y1,y2,y3)~u
#' regression(m) <- u~x1+x2
#' latent(m) <- ~u
#' covariance(m) <- y1~y2
#' categorical(m, labels = as.character(1:3)) <- "X1"
#'
#' defineCategoricalLink(m)
#'}
#' 
#' @concept setter
#' @keywords internal 
`defineCategoricalLink` <-
  function(object, link, data) UseMethod("defineCategoricalLink")


## * defineCategoricalLink.lvm
#' @rdname defineCategoricalLink
defineCategoricalLink.lvm <- function(object, link = NULL, data = NULL){

    ### ** normalize arguments
    if(is.null(link)){
        link <- stats::coef(object)
    }
    if(is.null(data)){
        data <- lava::sim(object, 1e2)
    }
    
    ### ** identify the type of regression variable (continuous or categorical)
    index.cat <- which(link %in% unlist(object$attributes$ordinalparname))
    index.Ncat <- setdiff(1:length(link), index.cat)
    link.Ncat <- setdiff(link[index.Ncat], names(object$attributes$ordinalparname))

    ### ** caracterize links involving categorical variables    
    if(length(index.cat)>0){
        link.cat <- link[index.cat]
        xCAT <- lava_categorical2dummy(object, data)$x

        ## *** find exogenous variable
        X.name.cat <- sapply(link.cat, function(iL){
            test <- unlist(lapply(object$attributes$ordinalparname, function(vec.coef){iL %in% vec.coef}))
            return(names(object$attributes$ordinalparname)[test])
        })
        UX.name.cat <- unique(X.name.cat)
    
        ## *** find the level of the exogenous variable
        X.level.cat <- unlist(lapply(UX.name.cat, function(iL){ 
            if(iL %in% names(xCAT$attributes$labels)){
                labels <- xCAT$attributes$labels[[iL]]
                index.label <- which(object$attributes$ordinalparname[[iL]] %in% link.cat)                
                return(labels[1+index.label])
            }else {
                stop("Categorical variables must have labels. Specify argument \'labels\' when calling categorical. \n")
            }            
        }))

        ## *** find endogenous variable
        M.link <- xCAT$M[paste0(X.name.cat,X.level.cat),,drop = FALSE]
        M.link <- cbind(M.link, as.numeric(rowSums(M.link)==0))

        convertion.back <- stats::setNames(X.name.cat,paste0(X.name.cat,X.level.cat))
                
        indexLink <- which(M.link==1, arr.ind = TRUE)
        Y.name.allcat <- colnames(M.link)[indexLink[,"col"]]
        X.name.allcat <- as.character(convertion.back[rownames(M.link)][indexLink[,"row"]])
        
        ## *** characterize all links
        Xcat.name.allcat <- rownames(M.link)[indexLink[,"row"]]
        X.level.allcat <- as.character(X.level.cat[indexLink[,"row"]])
        external.link.allcat <- link[index.cat[indexLink[,"row"]]]
        original.link.allcat <- paste0(Y.name.allcat, lava.options()$symbol[1], X.name.allcat)
        original.link.allcat[Y.name.allcat == ""] <- gsub("~","",original.link.allcat[Y.name.allcat == ""])
        cat.link.allcat <- paste0(Y.name.allcat, lava.options()$symbol[1], Xcat.name.allcat)
        cat.link.allcat[Y.name.allcat == ""] <- gsub("~","",cat.link.allcat[Y.name.allcat == ""])
        
        df.cat <- data.frame(link = cat.link.allcat,
                             endogenous = Y.name.allcat,
                             exogenous = X.name.allcat,
                             type = "categorical",
                             factitious = FALSE,
                             level = X.level.allcat,
                             originalLink = original.link.allcat,
                             externalLink = external.link.allcat,
                             stringsAsFactors = FALSE)       

    }else{
            df.cat <- NULL
    }

### ** caracterize links involving continuous variables    
    if(length(index.Ncat)>0){

        var.tempo <- initVarLinks(link.Ncat)
        Y.name.Ncat <- var.tempo$var1
        X.name.Ncat <- var.tempo$var2
        test.factitious <- X.name.Ncat %in% names(object$attributes$ordinalparname)

        df.Ncat <- data.frame(link = link.Ncat,
                              endogenous = Y.name.Ncat,
                              exogenous = X.name.Ncat,
                              type = "continuous",
                              factitious = test.factitious,
                              level = NA,
                              originalLink = link.Ncat,
                              externalLink = NA,
                              stringsAsFactors = FALSE)       
    }else{
        df.Ncat <- NULL
    }
    
    ### ** export
    out <- rbind(df.Ncat,df.cat)
    return(out)
}

#----------------------------------------------------------------------
### defineCategoricalLink.R ends here
