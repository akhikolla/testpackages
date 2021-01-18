cellMap = function(D, R, indcells = NULL, indrows = NULL,
                   standOD=NULL,showVals=NULL,rowlabels="",
                   columnlabels="",mTitle="", rowtitle="",
                   columntitle="",showrows=NULL, showcolumns=NULL,
                   nrowsinblock=1, ncolumnsinblock=1,autolabel=TRUE,
                   columnangle=90,sizetitles=1.1,adjustrowlabels=1,
                   adjustcolumnlabels=1, colContrast=1,outlyingGrad=TRUE,
                   darkestColor = sqrt(qchisq(0.999,1)), 
                   drawCircles = TRUE) {
  # Draws a cellmap, possibly of a subset of rows and columns of the data,
  # and possibly combining cells into blocks. 
  # The inputs are:
  #
  # D            : the data matrix (required input argument)
  # R            : matrix of cell residuals (required input argument)
  # indcells     : indices of outlying cells (required input argument)
  # indrows      : indices of outlying rows (required input argument)
  # standOD      : standardized Orthogonal Distance of each row
  # showVals     : whether to show the entries of "D" or "R" in the map. 
  #                Default is "D".
  # columnlabels      : labels for the x-axis
  # rowlabels      : labels for the y-axis
  # mTitle       : main title of the cellMap  
  # columntitle       : title for the x-axis
  # rowtitle       : title for the y-axis
  # showcolumns   : indices of the cells that will be shown, in the x direction
  # showrows   : indices of the cells that will be shown, in the y direction  
  # ncolumnsinblock   : size of combination blocks in the x direction
  # nrowsinblock   : size of combination blocks in the y direction
  # autolabel    : automatically combines labels of cells in blocks.  
  #                If FALSE, you must provide the final columnlabels and/or rowlabels.
  # columnangle       : angle of the labels on the x-axis
  # sizetitles       : size of title for x-axis and y-axis
  # adjustcolumnlabels : adjust x-labels: 0=left, 0.5=centered, 1=right
  # adjustrowlabels : adjust y-labels: 0=left, 0.5=centered, 1=right
  # colContrast  : adjust color contrast
  # outlyingGrad : use gradient colors for outlyingness
  # drawCircles: whether or not to draw the circles indicating casewise outliers
  
  funcSqueeze = function(Xin,n,d,ncolumnsinblock,nrowsinblock,colContrast) {
    # function to combine cells into blocks
    Xblock = matrix(0,nrow=n,ncol=d)
    Xblockgrad = matrix(0,nrow=n,ncol=d)
    for (i in 1:n) {
      for (j in 1:d) {
        Xsel = Xin[(1+((i-1)*nrowsinblock)):(i*nrowsinblock),
                   (1+((j-1)*ncolumnsinblock)):(j*ncolumnsinblock)]
        seltable = tabulate(Xsel,nbins=4) #changed
        # cnt0 = (nrowsinblock*ncolumnsinblock) - sum(seltable) #changed
        if (sum(seltable) > 0) {
          indmax = which(seltable==max(seltable))[1]
          cntmax = seltable[indmax]
          gradmax = (cntmax / (ncolumnsinblock*nrowsinblock))^(1/colContrast)
        } else {
          indmax = 0
          gradmax = 1
        }
        Xblock[i,j] = indmax
        Xblockgrad[i,j] = gradmax
      }
    }
    return(list(X=Xblock,Xgrad=Xblockgrad))
  }
  
  variable <- rownr <- rescaleoffset <- x <- y <- NULL
  
  type = "cell" # forced
  
  n = nrow(R)
  d = ncol(R)
  
  #check input arguments 
  blockMap = FALSE
  if (ncolumnsinblock > 1 | nrowsinblock > 1 ){
    blockMap = TRUE
    if (ncolumnsinblock > d) stop('Input argument ncolumnsinblock cannot be larger than d')
    if (nrowsinblock > n) stop('Input argument nrowsinblock cannot be larger than n')
    if (!is.null(showVals)) warning('The option showVals=D or showVals=R cannot
                                    be combined with ncolumnsinblock or nrowsinblock greater than 1,
                                    so showVals is set to NULL here.')
    showVals = NULL
  }
  
  if(!blockMap){
    if (!all(dim(R) == dim(D))) stop('Dimensions of D and R must match')
  }
  
  if(!(blockMap & autolabel==FALSE)){
    if (length(columnlabels) > 0 & length(columnlabels)!= d) {
      stop(paste('Number of columnlabels does not match d = ',d,sep=""))}
    if (length(rowlabels) > 0 & length(rowlabels)!= n) {
      stop(paste('Number of rowlabels does not match n = ',n,sep=""))}
  }
  
  if (!is.null(showVals)) {
    if (!showVals %in% c("D", "R")) {
      stop(paste("Invalid \"showVals\" argument. Should be one of: NULL, \"D\", \"R\""))
    }
  }
  
  if(is.null(indcells)) indcells = which(abs(R) > sqrt(qchisq(0.99,1)))
  
  if(!(is.null(showcolumns) & is.null(showrows))){
    # here we extract the rows and columns that will be shown.
    if(is.null(showcolumns)) { showcolumns = 1:d } else {
      if(!(all(showcolumns %in% 1:d))) stop(" showcolumns goes out of bounds")}
    if(is.null(showrows)) { showrows = 1:n } else {
      if(!(all(showrows %in% 1:n))) stop(" showrows goes out of bounds")}
    
    tempMat = matrix(0,n,d)
    tempMat[indcells] = 1
    tempMat = tempMat[showrows,showcolumns]
    indcells = which(tempMat == 1)
    
    tempVec = rep(0,n)
    tempVec[indrows] = 1
    tempVec = tempVec[showrows]
    indrows = which(tempVec == 1) # also works if indrows was empty
    rm(tempMat,tempVec)
    
    if(!blockMap) D = D[showrows,showcolumns] 
    R = R[showrows,showcolumns]
    if (!(blockMap & autolabel==FALSE)) columnlabels = columnlabels[showcolumns]
    if (!(blockMap & autolabel==FALSE)) rowlabels = rowlabels[showrows]
    n = nrow(R)
    d = ncol(R)
    if(!is.null(standOD)) standOD = standOD[showrows] 
  }  
  
  # force outlyingGrad when type is residual
  if (type=="residual") outlyingGrad=1
  
  # create the matrix which indicates the color of each cell
  # 0=yellow, 1=blue, 2=red, 3=black, 4=white
  X = matrix(0,n,d)
  # X[indrows,] = 3
  Xrow = matrix(0,n,1)
  Xrow[indrows,1] = 3
  
  if (type=="cell" | blockMap){
    pcells = indcells[indcells %in% which(R>=0)]
    ncells = indcells[indcells %in% which(R<0)]
  } else { #residual
    pcells = which(R>=0)
    ncells = which(R<0)
  }
  X[ncells] = 1
  X[pcells] = 2
  X[is.na(R)] = 4 #changed
  
  if (blockMap) { # in this case the input D will be ignored  
    n = floor(n/nrowsinblock)
    d = floor(d/ncolumnsinblock)
    
    # create new X and Xgrad
    result = funcSqueeze(X,n,d,ncolumnsinblock,nrowsinblock,colContrast)
    X = result$X
    Xgrad = result$Xgrad
    result = funcSqueeze(Xrow,n,1,1,nrowsinblock,colContrast)
    Xrowgrad = result$Xgrad
    Xrowgrad[result$X==0] = 0
    
    # if autolabel=F, labels{x,y} will be used for the blocks.
    if (autolabel==TRUE) { # automatically combine labels for blocks
      
      if (ncolumnsinblock>1 & length(columnlabels)>0) {
        labx = columnlabels
        columnlabels = rep(0,d)
        for(ind in 1:d) {
          columnlabels[ind] = paste(labx[(1+((ind-1)*ncolumnsinblock))],"-",
                               labx[(ind*ncolumnsinblock)],sep="")
        }
      }
      
      if (nrowsinblock>1 & length(rowlabels)>0) {
        laby = rowlabels
        rowlabels = rep(0,n)
        for(ind in 1:n) {
          rowlabels[ind] = paste(laby[(1+((ind-1)*nrowsinblock))],"-",
                               laby[(ind*nrowsinblock) ])
        }
      }
    } else {
      if (length(columnlabels) > 0 & length(columnlabels)!= d) {
        stop(paste(' autolabel=FALSE and number of columnlabels is ',
                   length(columnlabels),' but should be ',d,sep=""))}
      if (length(rowlabels) > 0 & length(rowlabels)!= n) {
        stop(paste(' autolabel=FALSE and number of rowlabels is ',
                   length(rowlabels),' but should be ',n,sep=""))}      
    }
    
    # Melt data matrices for cellMap
    Xdf = data.frame(cbind(seq(1,n,1),X))
    colnames(Xdf) = c("rownr",seq(1,d,1))
    rownames(Xdf) = NULL
    Xdf$rownr = with(Xdf, reorder(rownr, seq(n,1,-1)))
    mX = melt(Xdf,id.var="rownr", value.name = "CatNr") 
    
    Xgraddf = data.frame(cbind(seq(1,n,1),Xgrad))
    colnames(Xgraddf) = c("rownr",seq(1,d,1))
    rownames(Xgraddf) = NULL
    Xgraddf$rownr = with(Xgraddf, reorder(rownr, seq(n,1,-1)))
    mXgrad = melt(Xgraddf,id.var="rownr", value.name = "grad") 
    
    # Combine melted data
    mX$grad = mXgrad$grad
    mX$rescaleoffset = mXgrad$grad + 10*mX$CatNr
    
    mXrow = data.frame(rownr=1:n,rescaleoffset=Xrowgrad+ 10*3)
    
    scalerange = c(0,1)
    gradientends = scalerange + rep(c(0,10,20,30,40), each=2)  #changed
    if (type=="cell") colorends = c("yellow", "yellow", "yellow", "blue",
                                    "yellow", "red","white", "black", 
                                    "yellow", "white")
    if (type=="residual") colorends = c("white", "white", "white", "blue", 
                                        "white", "red","white", "black",
                                        "white", "white")
  } else { # no blockMap
    
    # Melt data matrices for cellMap
    Ddf = data.frame(cbind(seq(1,n,1),D))
    colnames(Ddf) = c("rownr",seq(1,d,1))
    rownames(Ddf) = NULL
    Ddf$rownr = with(Ddf, reorder(rownr, seq(n,1,-1)))
    mD = melt(Ddf,id.var="rownr") 
    
    Rdf = data.frame(cbind(seq(1,n,1),R))
    colnames(Rdf) = c("rownr",seq(1,d,1))
    rownames(Rdf) = NULL
    Rdf$rownr = with(Rdf, reorder(rownr, seq(n,1,-1)))
    mR = melt(Rdf,id.var="rownr") 
    
    Xdf = data.frame(cbind(seq(1,n,1),X))
    colnames(Xdf) = c("rownr",seq(1,d,1))
    rownames(Xdf) = NULL
    Xdf$rownr = with(Xdf, reorder(rownr, seq(n,1,-1)))
    mX = melt(Xdf,id.var="rownr", value.name = "CatNr") 
    
    if (!is.null(showVals)) {
      # Combine melted data
      if (showVals=="D") mX$data = mD$value
      if (showVals=="R") mX$data = mR$value
    }
    
    if (!outlyingGrad) {
      mX$rescaleoffset = 10*mX$CatNr
      
      scalerange = c(0,1)
      gradientends = scalerange + rep(c(0,10,20,30,40), each=2)
      gradientends
      colorends = c("yellow", "yellow", "blue", "blue", "red", "red",
                    "white", "black","white", "white")
    } else {
      Xgrad = matrix(NA,n,d)
      if (type=="cell") {
        Xgrad[indcells] = abs(R[indcells]) 
        limL = sqrt(qchisq(0.9,1))
      } else {
        Xgrad = abs(R) 
        limL = 0
      }
      
      limH = darkestColor
      Xgrad[Xgrad>limH] = limH
      Xgrad = ( (Xgrad -limL) / (limH - limL) )^colContrast
      Xgrad[is.na(Xgrad)] = 0
      
      Xgraddf = data.frame(cbind(seq(1,n,1),Xgrad))
      colnames(Xgraddf) = c("rownr",seq(1,d,1))
      rownames(Xgraddf) = NULL
      Xgraddf$rownr = with(Xgraddf, reorder(rownr, seq(n,1,-1)))
      mXgrad = melt(Xgraddf,id.var="rownr", value.name = "grad") 
      
      mX$grad = mXgrad$grad
      mX$rescaleoffset = mXgrad$grad + 10*mX$CatNr
      scalerange = c(0,1)
      gradientends = scalerange + rep(c(0,10,20,30,40), each=2)
      if (type=="cell") colorends = c("yellow", "yellow", "yellow", "blue",
                                      "yellow", "red", "white", "black",
                                      "white", "white")  
      if (type=="residual") colorends = c("white", "white", "white", "blue",
                                          "white", "red", "white", "black",
                                          "white", "white")  
    }
    
    tempVec = rep(0,n)
    tempVec[indrows] = 1
    mXrow = data.frame(rownr=1:n,rescaleoffset=40-(10*tempVec) )
    rm(tempVec)
    if (is.null(standOD)) {
      mXrow$rescaleoffset[indrows] = mXrow$rescaleoffset[indrows] + 1
    } else {
      limL = 1
      limH = 3
      standOD[standOD>limH] = limH
      standOD = ( (standOD - limL) / (limH - limL) )^colContrast
      mXrow$rescaleoffset[indrows] = mXrow$rescaleoffset[indrows] + 
        standOD[indrows]
    }
  }
  
  rowlabels = rev(rowlabels) # will be plotted from bottom to top
  base_size = 10
  columnlabels = c(columnlabels,"","")
  
  circleFun = function(centerx, centery, r, npoints) {
    tt = seq(0, 2 * pi, length.out = npoints)
    xx = centerx + r * cos(tt)
    yy = centery + r * sin(tt)
    return(c(xx, yy))
  }
  if (drawCircles) {
    centerx = d + 1
    centery = n:1
    radius = 0.4
    npoints = 100
    circlePoints = mapply(circleFun, centerx, centery, radius,
                          npoints)
    positions = data.frame(rownr = rep(1:n, each = npoints),
                           x = c(circlePoints[1:npoints, ]), 
                           y = c(circlePoints[(npoints +1):(2*npoints), ]))
    datapoly = merge(mXrow, positions, by = c("rownr"))
  }
  
  ggp = ggplot(data = mX, aes(variable, rownr)) + {
    if (blockMap) 
      geom_tile(aes(fill = rescale(rescaleoffset, 
                                   from = range(gradientends))), 
                color = "white")
  } + {
    if (!blockMap & outlyingGrad) 
      geom_tile(aes(fill = rescale(rescaleoffset, 
                                   from = range(gradientends))), 
                color = "white")
  } + {
    if (!blockMap & !outlyingGrad) 
      geom_tile(aes(fill = rescale(rescaleoffset, 
                                   from = range(gradientends))), 
                colour = "white")
  } + {
    if (drawCircles)
      geom_polygon(data = datapoly,
                   aes(x = x, y = y,
                       fill = rescale(rescaleoffset, from = range(gradientends)),
                       group = rownr), colour = "black") } + 
    scale_fill_gradientn(colours = colorends, 
                         values = rescale(gradientends), 
                         rescaler = function(x, ...) x, oob = scales::squish) + 
    ggtitle(mTitle) + coord_fixed() +
    theme_classic(base_size = base_size * 1) +
    labs(x = columntitle, y = rowtitle) + 
    scale_x_discrete(expand = c(0, 0), 
                     limits = as.factor(seq(1, d + 2, 1)), 
                     labels = columnlabels) + 
    scale_y_discrete(expand = c(0, 0), labels = rowlabels) + 
    theme(legend.position = "none", axis.ticks = element_blank(), 
          plot.title = element_text(size = base_size * 2, hjust = 0.5, 
                                    vjust = 1, face = "bold"), 
          axis.text.x = element_text(size = base_size * 1.8, 
                                     angle = columnangle, hjust = adjustcolumnlabels,
                                     vjust = 0.5, colour = "black"), 
          axis.text.y = element_text(size = base_size * 1.8, 
                                     angle = 0, hjust = adjustrowlabels, colour = "black"), 
          axis.title.x = element_text(colour = "black", size = 
                                        base_size * sizetitles, vjust = 1), 
          axis.title.y = element_text(colour = "black", size = 
                                        base_size * sizetitles, vjust = 0), 
          axis.line.x = element_blank(),
          panel.border = element_blank()) + 
    annotate(geom = "segment", x = 0.5, xend = d + 0.5, 
             y = 0.5, yend = 0.5) + 
    annotate(geom = "segment", x = 0.5, xend = d + 0.5,
             y = n + 0.5, yend = n + 0.5) +
    annotate(geom = "segment", x = d + 0.5, xend = d + 0.5, 
             y = 0.5, yend = n + 0.5)
  if (!is.null(showVals)) {
    txtcol = mX$CatNr
    txtcol[txtcol==0] = "black"
    txtcol[txtcol==1] = "white"
    txtcol[txtcol==2] = "white"
    # txtcol[txtcol==3] = "white"
    if (type=="residual") {
      txtcol[]="black"
      txtcol[mXgrad$grad>0.5] = "white"
    }
    txtcol[txtcol==4] = "black"
    # the following line places `NA' in missing cells
    # the size of the symbols is also specified here
    ggp = ggp + geom_text(aes(label = ifelse(is.na(data),
                                             sprintf("%1.0f",data),
                                             round(data,1))),size = base_size*0.5, 
                          colour=txtcol, na.rm = TRUE) 
  } 
  return(ggp)
} # ends cellMap