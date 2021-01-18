#' @title Visualisation of multiple locus association mapping results
#' @description    A interactive plotting function that provides additional information on the significant 
#'     marker-trait associations found by \code{\link{AM}}
#' @param  AMobj  the (list) object obtained from running \code{\link{AM}}. 
#' @param  itnum  the iteration number of the model building process whose results are to be viewed.
#' @param  chr  either "All" for all chromosomes or the label of a specific chromosome.
#' @param  type either "Manhattan" or  "Score"  for a manhattan plot or a plot of the 
#'    score statistics across SNPs. 
#' @param interactive boolean parameter. When \code{TRUE}, an interactive plot is generated. When 
#'         \code{FALSE}, a ggplot object is returned which can be saved as an image to file 
#          or viewed on-screen.  
#'
#' @details
#'  A function useful for viewing the strength of association across the whole genome and 
#'  how this association changes as the model is built. 
#'
#'  The  score statistics (\code{type="Score"}) or p-values  of the score statistics
#'  (\code{type="Manhattan"}) are plotted against the location of the SNPs.  The orange
#'  vertical  lines 
#'  denote the location of the SNPs already found by \code{\link{AM}}. 
#'  The red vertical line is the location of the SNP in strongest association with the trait 
#'  at that iteration number. 
#'
#'  The vertical lines are numbered according to the order in which the snp-trait associations were found by the model. 
#'
#'  A single chromosome or all (\code{chr="All"}) chromosomes can be viewed.  
#'
#'  By setting \code{itnum} to different values, how the score statistics or p-values 
#'  increase/decrease 
#'  over the model building process can be observed. 
#'
#' 
#' @examples
#'  \dontrun{
#'   # Since the following code takes longer than 5 seconds to run, it has been tagged as dontrun. 
#'   # However, the code can be run by the user. 
#'   #
#'
#'   #---------------
#'   # read the map 
#'   #---------------
#'   #
#'   # File is a plain space separated text file with the first row 
#'   # the column headings
#'   complete.name <- system.file('extdata', 'map.txt', 
#'                                    package='Eagle')
#'   map_obj <- ReadMap(filename=complete.name) 
#'
#'  # to look at the first few rows of the map file
#'  head(map_obj)
#'
#'   # read marker data
#'   complete.name <- system.file('extdata', 'geno.ped', 
#'                                      package='Eagle')
#'   geno_obj <- ReadMarker(filename=complete.name,  type='PLINK', availmemGb=8) 
#'  
#'   # read phenotype data
#'   complete.name <- system.file('extdata', 'pheno.txt', package='Eagle')
#'   
#'   pheno_obj <- ReadPheno(filename=complete.name)
#'            
#'   # Perform multiple-locus genome-wide association mapping 
#'   res <- AM(trait = 'y',
#'                            fformula=c("cov1 + cov2"),
#'                            map = map_obj,
#'                            pheno = pheno_obj,
#'                            geno = geno_obj)
#'
#'  # Plotting the p-values from the first iteration of the module building process. 
#'  # You can see why Eagle has identified the SNP that is has. 
#'   PlotAM(AMobj=res, itnum=1)
#'
#'
#'  # Plotting the results from the final step of the model building process
#'  # By accounting for the effect of SNP in strong association with the trait, the 
#'  # strength of association changes across the genome. 
#'   PlotAM(AMobj=res, itnum=3)
#'
#'
#'  # Suppose you want to save the above plot to a jpeg file called myplot.jpg
#'  jpeg("./myplot.jpg", width=1200, height=800)
#'  PlotAM(AMobj=res, itnum=3, interactive=FALSE)
#'  dev.off()
#'
#'
#'  }
#'
#' 
#' 
#' @seealso \code{\link{AM}}
#'
PlotAM <- function(AMobj=NULL, itnum=1, chr="All", type="Manhattan", interactive=TRUE )
{


 if(is.null(AMobj)){
    message(" PlotAM function requires AMobj object to be specified. This object is obtained by running AM().")
    return(NULL)
    }
 if(!is.list(AMobj)){
    message(" PlotAM function requires AMobj object to be a list object.")
    return(NULL)
   }

if (!is.integer(itnum)){

 if(itnum < 1 ||  (itnum > length(AMobj$outlierstat)  )){
    message(" The parameter itnum must have an integer value between ", 1, " and ",
         length(AMobj$outlierstat))
    return(NULL)
  }
}

 if (!(chr=="All")  &&  !any(chr %in%  unique(AMobj$map[,2])) ) {
    message(" The parameter chr does not match any of the chromosome labels. ")
    message(" The allowable chromosome labels are ", c("All", cat(unique(AMobj$map[,2]))))
    return(NULL)
 }



# we do not have a map
xindx <- 1:length( AMobj$outlierstat[[itnum]] )
xvals <- xindx

yvals <- AMobj$outlierstat[[itnum]]
isit  <- IsItBigger(vals=AMobj$outlierstat, itnum=itnum )
bigger <- isit$bigger

chrm <- rep(1, length(xindx))
pos  <- xvals

# map exisits
if(!is.null(AMobj$map)){
  # plotting all the chromosomes - more difficult
  # reordering based on chrm then map position
  oindx <- order(AMobj$map[,2], AMobj$map[, ncol(AMobj$map)])   # **** VERY IMPORTANT OBJECT 
  yvals <- AMobj$outlierstat[[as.numeric(itnum)]][oindx]  ## reordering yvals
  mapordered <- AMobj$map[oindx,]
  chrm <- mapordered[,2]

  # Cumulative map position (genome-wide position)
  pos <- rep(NA, nrow(mapordered))
  t <- 0
  for(ii in 1:nrow(mapordered)){
    t <- t + mapordered[ii,ncol(mapordered)]
    pos[ii] <- t
  }






  if(length(AMobj$Indx) > 1){
    # AMobj$Indx has column positions of SNPs in geno file but if map order has changed,
    #   then this Indx also need to change accordingly
    Indx <- rep(NA, length(AMobj$Mrk))  # contained reordered AMobj$Indx 
    NewMrkPos <- rep(NA, length(AMobj$Mrk)) # new map position of sig SNP based on mapordered
    for(ii in 2:length(AMobj$Mrk)){
      Indx[ii] <- which(mapordered[,1]==AMobj$Mrk[ii])
      NewMrkPos[ii] <- pos[Indx[ii]]  ## map position
    }
    Indx <- Indx[!is.na(Indx)]
    NewMrkPos <- NewMrkPos[!is.na(NewMrkPos)]
  }





  if( as.numeric(itnum)  > 1){
    bigger <- rep(""    , length(yvals) )
    percentagechange <- rep(0, length(yvals) )


    a <-  AMobj$outlierstat[[as.numeric(itnum)]][oindx]
    b <-  AMobj$outlierstat[[as.numeric(itnum) - 1 ]][oindx]

    indx <- which(  a >  b )
    bigger[indx] <- "Increased value"
    percentagechange[indx] <-  (( b - a)) [indx]

    indx <- which(  a <=  b )
    bigger[indx] <- "Decreased value"
    percentagechange[indx] <-  (( a - b)) [indx]
      

  }  ## end if( as.numeric(itnum)  > 1)



}  ##  if(!is.null(map))

xlabel <- "Map Position (bp)"
if(is.null(AMobj$map ))
   xlabel <- "Column Position of SNP"

ylabel <- "Score Statistic"
if(type=="Manhattan")
   ylabel <- "-log10(p value)"


 # addition on SNP-trait positions on map


   # place on -lgo10 scale if manhattan selected
   if(type=="Manhattan"){
       yvals[is.nan(yvals)] <- 0
       yvals[yvals < 0] <- 0  ## rounding error - very close to 0 when negative
       ts <- sqrt(yvals)
       pval <- 1 - pnorm(ts)
       logp <- -1*log10(pval)
       yvals <- logp
   }


   # create data frame for plotting 
   df <- data.frame(xvals=pos, yvals=yvals, chrm=chrm )


   # subset df on chromosome if chr!="All"
   if(chr!="All"){
       indx <- which(df$chrm==chr)
       df <- df[indx,]
       if (itnum >1){
          bigger <- bigger[indx]
          percentagechange <- percentagechange[indx]
       }
       # df <- subset(df, chrm==chr)
   }


  geomX <- NULL
  if(length(AMobj$Indx) > 1){
    Labels <- 1:length(Indx)
    indx <-  which( NewMrkPos %in% df$xvals )
    if (length(indx) > 0 ){
     geomX <- NewMrkPos[indx] 
     geomLabels <- Labels[indx ]
    }
   }


   if(itnum==1){
       p  <- ggplot(data=df, aes(x=xvals, y=yvals )) + geom_point()

   } else {
       # making the points easier to view 
       percentagechange <- abs(percentagechange)  ## getting rid of negs as we don't need to worry about
                                                  ## the sign 
       sd <- sqrt(var(percentagechange))
       indx <- which( 2.0 *sd <  percentagechange)
       percentagechange[indx]  <- 2.0*sd  ## setting outliers to 3

       p  <- ggplot(data=df, aes(x=xvals, y=yvals , color=bigger, size=percentagechange ) )  + 
           geom_point() +  scale_color_manual(values=c("#6600CC","#4C9900" ))
       p <- p + scale_size(percentagechange, range=c(0.4  ,2 ), guide="none") 
 
  

   }



   p <- p + theme_hc()
   p <- p + ylab(ylabel) + xlab(xlabel)
   p <- p +  theme(legend.title=element_blank())  ## no legend title
   p <- p + theme(legend.position="right")

   if(!is.null(geomX)){
      if (itnum >1){
         # draw vertical line for found snp-trait associations
         for(ii in geomLabels) {
            if (ii < itnum){
               yadj <- sample(seq(0.5,0.9,0.1), 1)
               p <- p + geom_vline(xintercept = geomX[which(ii==geomLabels)], 
                        linetype="solid", color="#FFE4B5", size=0.75)
               p <- p + annotate("text", size=6, label=ii,
                       x=(geomX[which(ii==geomLabels)]  - ( diff(range(df$xvals))*0.00) )   ,
                        y = max(df$yvals)*yadj )
             } ## if ( ii < itnum)

             if (ii == itnum){
                  p <- p + geom_vline(xintercept = geomX[which(ii==geomLabels)] , 
                             linetype="solid", color="red", size=0.75)
             }
           }  ## end for ii
      } else { 

         ii  <- 1
         p <- p + geom_vline(xintercept = geomX[which(ii==geomLabels)], 
                   linetype="solid", color="red", size=0.75)

      } ## end if itnum 
    #  p <- p + scale_size(guide='none')
   } ## if !is.null
   p <- p + theme( axis.text.x = element_text(size=16), axis.text.y = element_text(size=16),
                           axis.title.x=element_text(size=16), axis.title.y=element_text(size=16),
                           legend.text=element_text(size=14) )
   p <- p  + guides(colour = guide_legend(override.aes = list(size=8)))  ## changing point size in legend

   p <- p + theme(legend.position = 'bottom', legend.spacing.x = unit(0.5, 'cm'))

           
   if(chr=="All"){
     # entire genome
     txt2 <- "across all chromosomes"
     if(type=="Manhattan"){
       txt1 <- " -log p value of the score statistic"
       txt3 <- "-log p value"
     } else {
       txt1 <- "score statistic"
       txt3 <- txt1
     } ## inner else
   } else{
     # chrm selected
     txt2 <- paste("on chromosome", chr)
     if(type=="Manhattan"){
       txt1 <- " -log p value of the score statistic"
       txt3 <- "-log p value"
     } else {
       txt1 <- "score statistic"
       txt3 <- txt1
     } ## inner else
   }  ## end outer else



   # adjust  Y axis scale
   p <- p + ylim(0, (max(df$yvals)*1.2) )


   # increase point size in legend
   p <- p + guides(colour = guide_legend(override.aes = list(size=5)))
 
 if(interactive){
  # launch interactive graphics
  ggplotly(p)  %>% layout(legend = list( orientation = "h", x=0, y=-0.4, itemsizing="constant"  ))
  } else {
    return(p)
  }


}  ## end function PlotAM ... 


IsItBigger <- function(vals, itnum, xindx=NULL){
   ## calculates percentage change in size (increase or decrease)
   bigger <- NULL
   percentagechange <- NULL
   if(as.numeric(itnum) > 1){
       # entire genome
       bigger <- rep("" , length(vals[[as.numeric(itnum)]]) ) ## initialize 
       percentagechange <- rep(0, length(vals[[as.numeric(itnum)]]) ) ## initialize
       a <-  vals[[as.numeric(itnum)]]
       b <- vals[[as.numeric(itnum) - 1 ]]
       # a > b       
       indx <- which( a  >  b )
       bigger[indx] <- "Increased value"
       percentagechange[indx] <- ( (a  -    b ) )[indx]
       # a < b
       indx <- which( a <=  b  )
       bigger[indx] <- "Decreased value"
       percentagechange[indx] <- ( (b  -    a ) )[indx]

       if(!is.null(xindx)){
          # reduce to chromosome
          bigger <- bigger[xindx]
          percentagechange <- percentagechange[xindx]
       }
     }
     res <- list(bigger=bigger, percentagechange=percentagechange)
     return(res)
  }



