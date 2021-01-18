MultiManhattan<-function(ResultIntermediate,ResultFinal,mar=c(2.9,2.8,0.7,2.8),LabDistance=1.5,
                         ScaleDistance=0.4,LabelSize=0.8,ScaleSize=0.7,AxisLwd=5,
                         TckLength=-0.03,LogTimes=2,LODTimes=1.2,lodline=3,
                         dirplot=getwd(),PlotFormat="tiff",width=28000, height=7000, 
                         pointsize = 60,res=600,MarkGene=FALSE,Pos_x=NULL,Pos_y=NULL,
                         GeneName=NULL,GeneNameColour=NULL,...){
  
  ###########Data process#################
  ###########intermediate result
  ##################################
  if(is.character(ResultIntermediate)==TRUE&is.character(ResultFinal)==TRUE){
  ResultIntermediate<-as.matrix(read.csv(ResultIntermediate,header = F))
  ResultFinal<-as.matrix(read.csv(ResultFinal,header = F))
  }else{
  ResultIntermediate<-as.matrix(ResultIntermediate)
  ResultFinal<-as.matrix(ResultFinal)  
  }
  ######################################
  cona<-ResultIntermediate[1,]
  ResultIntermediate<-ResultIntermediate[-1,]
  ResultFinal<-matrix(ResultFinal[-1,],,14)
  logp_4method<-as.matrix(ResultIntermediate[,which(substr(cona,3,7)=="log10")])
  logp_4method<-apply(logp_4method,2,as.numeric)
  p_4method<-10^-logp_4method
  p_median<-apply(p_4method,1,median)
  locsub<-which(p_median==0)
  pmin<-min(p_median[p_median!=0])
  subvalue<-10^(1.1*log10(pmin))
  p_median[locsub]<-subvalue
  data_p<-as.matrix(p_median)
  data_num<-as.matrix(seq(1:length(p_median)))
  data_chr<-as.matrix(ResultIntermediate[,4])
  data_pos<-as.matrix(ResultIntermediate[,5])
  manresult<-cbind(data_chr,data_pos,data_p,data_num)
  manresult<-apply(manresult,2,as.numeric)
  colnames(manresult)<-c("Chromosome","BPnumber","P-value","SNPname")
  manresult<-as.data.frame(manresult)
  #######final result##################
  data_fin_method<-unique(ResultFinal[,3])
  data_fin_method_length<-1:length(unique(ResultFinal[,3]))
  for(r in 1:length(unique(ResultFinal[,3]))){
    ResultFinal[which(ResultFinal[,3]==data_fin_method[r]),3]<-r
  }
  data_fin_mark<-matrix(ResultFinal[,c(5,6,8,3)],,4)
  data_fin_mark<-matrix(apply(data_fin_mark,2,as.numeric),,4)
  data_fin_mark_chr<-matrix(data_fin_mark[order(data_fin_mark[,1]),],,4)
  data_fin_mark_order<-numeric()
  for(i in c(unique(data_fin_mark_chr[,1]))){
    data_fin_mark_erery_chr<-matrix(data_fin_mark_chr[which(data_fin_mark_chr[,1]==i),],,4)
    data_fin_mark_pos<-matrix(data_fin_mark_erery_chr[order(data_fin_mark_erery_chr[,2]),],,4)
    all_pos<-unique(data_fin_mark_pos[,2])
    all_pos_maxlod<-numeric()
    for(ii in 1:length(all_pos)){
      all_pos_every<-matrix(data_fin_mark_pos[which(data_fin_mark_pos[,2]==all_pos[ii]),],,4)
      lod_me<-median(all_pos_every[,3])
      all_pos_every_median<-c(all_pos_every[1,1:2],lod_me,all_pos_every[1,4])
      if(nrow(all_pos_every)>=2){
        all_pos_every_median<-c(all_pos_every[1,1:2],lod_me,max(data_fin_mark[,4])+1)
      }
      all_pos_maxlod<-rbind(all_pos_maxlod,all_pos_every_median)
    }
    data_fin_mark_order<-rbind(data_fin_mark_order,all_pos_maxlod)
  }
  snpOfInterest<-numeric()
  for(i in c(unique(data_fin_mark_order[,1]))){
    manresult_chr<-manresult[which(manresult[,1]==i),]
    data_fin_mark_order_chr<-matrix(data_fin_mark_order[which(data_fin_mark_order[,1]==i),],,4)
    mark_loc<-manresult_chr[which(manresult_chr[,2]%in%data_fin_mark_order_chr[,2]),4]
    snpOfInterest<-c(snpOfInterest,mark_loc) 
  }
  bpnumber <- numeric()
  chrnum <- unique(manresult[,1])
  for(i in 1:length(chrnum))
  {
    bpnumber <- rbind(bpnumber,as.matrix(c(1:length(which(manresult[,1]==chrnum[i])))))
  }
  manresult2<-cbind(manresult[,1],bpnumber,manresult[,3:4])
  colnames(manresult2)<-c("Chromosome","BPnumber","P-value","SNPname")
  ##########prepare for data#############################
  x<-manresult2;col=c("lightgreen","lightskyblue");logp=TRUE
  chr = "Chromosome";bp ="BPnumber";p ="P-value";snp="SNPname";
  highlight<-snpOfInterest
  CHR=BP=P=index=NULL
  d=data.frame(CHR=x[[chr]], BP=x[[bp]], P=x[[p]])
  if (!is.null(x[[snp]])) d=transform(d, SNP=x[[snp]])
  d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
  d <- d[order(d$CHR, d$BP), ]
  if (logp) {
    d$logp <- -log10(d$P)
  } else {
    d$logp <- d$P
  }
  d$pos=NA
  d$index=NA
  ind = 0
  for (i in unique(d$CHR)){
    ind = ind + 1
    d[d$CHR==i,]$index = ind
  }
  
  nchr = length(unique(d$CHR))
  if (nchr==1) { ## For a single chromosome
    ## Uncomment the next two linex to plot single chr results in Mb
    #options(scipen=999)
    #d$pos=d$BP/1e6
    d$pos=d$BP
    ticks=floor(length(d$pos))/2+1
    xlabel = paste('Chromosome',unique(d$CHR),'position')
    labs = ticks
  } else { ## For multiple chromosomes
    lastbase=0
    ticks=NULL
    for (i in unique(d$index)) {
      if (i==1) {
        d[d$index==i, ]$pos=d[d$index==i, ]$BP
      } else {
        lastbase=lastbase+tail(subset(d,index==i-1)$BP, 1)
        d[d$index==i, ]$pos=d[d$index==i, ]$BP+lastbase
      }
      # Old way: assumes SNPs evenly distributed
      # ticks=c(ticks, d[d$index==i, ]$pos[floor(length(d[d$index==i, ]$pos)/2)+1])
      # New way: doesn't make that assumption
      ticks = c(ticks, (min(d[d$index == i,]$pos) + max(d[d$index == i,]$pos))/2 + 1)
    }
    xlabel = 'Chromosomes'
    #labs = append(unique(d$CHR),'') ## I forgot what this was here for... if seems to work, remove.
    labs <- unique(d$CHR)
  }
  xmax = ceiling(max(d$pos) * 1.03)
  xmin = floor(max(d$pos) * -0.03)
  ########draw plot#######################
  if(PlotFormat=="png"){
    png(paste(dirplot,"/","Manhattan plot.png",sep=""),width=width, height=height, units= "px", pointsize = pointsize,res=res)
  }else if(PlotFormat=="tiff"){
    tiff(paste(dirplot,"/","Manhattan plot.tiff",sep=""),width=width, height=height, units= "px", pointsize = pointsize,res=res,compression="lzw")
  }else if(PlotFormat=="jpeg"){
    jpeg(paste(dirplot,"/","Manhattan plot.jpeg",sep=""),width=width, height=height, units= "px", pointsize = pointsize,res=res)
  }else if(PlotFormat=="pdf"){
    pdf(paste(dirplot,"/","Manhattan plot.pdf",sep=""),width=width,height=height,pointsize=pointsize)
  }
  
  par(mar=mar)
  def_args <- list(xaxt='n',yaxt="n",bty='n', xaxs='i', yaxs='i', las=1, pch=20,
                   xlim=c(xmin,xmax), ylim=c(0,LogTimes*max(d$logp)),
                   xlab=xlabel,ylab="",mgp=c(LabDistance,0,0),cex.lab=LabelSize)
  
  dotargs <- list(NULL)
  do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))
  #axis(1, at=ticks, labels=labs,lwd=CoorLwd,tck=TckLen,mgp=c(2.5,HorTckDis,0.5),cex.axis=TckLwd)
  axis(1, at=ticks, labels=labs,lwd=AxisLwd,tck=TckLength,mgp=c(LabDistance+1.2,ScaleDistance-0.2,0.5),cex.axis=ScaleSize,...)
  
  suppressWarnings(axis(2, at=seq(0,LogTimes*max(d$logp),ceiling(LogTimes*max(d$logp)/5)),lwd=AxisLwd,tck=TckLength,mgp=c(2.2,ScaleDistance,0),cex.axis=ScaleSize,...))
  mtext(expression(-log[10]('P-value')),side=2,line=LabDistance,cex=LabelSize,font=1)
  # Create a vector of alternatiting colors
  col=rep(col, max(d$CHR))
  # Add points to the plot
  if (nchr==1) {
    with(d, points(pos, logp, pch=20, col=col[1],...))
  } else {
    # if multiple chromosomes, need to alternate colors and increase the color index (icol) each chr.
    icol=1
    for (i in unique(d$index)) {
      with(d[d$index==unique(d$index)[i], ], points(pos, logp, col=col[icol], pch=20,...))
      icol=icol+1
    }
  }
  d.highlight=d[which(d$SNP %in% highlight), ]
  highlight_LOD<-as.numeric(data_fin_mark_order[,3])
  d.highlight<-as.data.frame(cbind(d.highlight,highlight_LOD))
  
  ################################
  par(new=T)
  def_args <- list(xaxt='n', yaxt='n',bty='n', xaxs='i', yaxs='i', las=1, pch=20,
                   xlim=c(xmin,xmax), ylim=c(0,LODTimes*max(highlight_LOD)),xlab="",ylab="")
  dotargs <- list(NULL)
  do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))
  suppressWarnings(axis(4,mgp=c(1.4,ScaleDistance,0),at=seq(0,LODTimes*max(highlight_LOD),ceiling(LODTimes*max(highlight_LOD)/5)),col="magenta",col.ticks="magenta",col.axis="magenta",lwd=AxisLwd,tck=TckLength,cex.axis=ScaleSize,...))
  mtext("LOD score",side=4,line=LabDistance,cex=LabelSize,font=1,col="magenta")
  abline(h=lodline,col="gray25",lty=2,lwd=2)
  peach_colors<-c("magenta","deepskyblue2")
  col_pos<-list(NULL)
  method_num<-sort(unique(data_fin_mark_order[,4]))
  
  if(max(unique(ResultFinal[,3]))<max(unique(data_fin_mark_order[,4]))){
    col_pos[[1]]<-which(data_fin_mark_order[,4]==max(method_num))
    col_pos[[2]]<-which(data_fin_mark_order[,4]!=max(method_num))
  }else{
    if(length(unique(ResultFinal[,3]))==1){
      col_pos[[1]]<-which(data_fin_mark_order[,4]==max(method_num))
    }else{
      col_pos[[1]]<-1:nrow(data_fin_mark_order)
      
    }
  }
  if(length(col_pos)>1&&length(col_pos[[2]])!=0){
    with(d.highlight, points(pos[col_pos[[2]]], highlight_LOD[col_pos[[2]]], col=peach_colors[2], pch=20,...))
    with(d.highlight, points(pos[col_pos[[2]]], highlight_LOD[col_pos[[2]]], col=peach_colors[2], type="h",lty=2,...))
    with(d.highlight, points(pos[col_pos[[1]]], highlight_LOD[col_pos[[1]]], col=peach_colors[1], pch=20,...))
    with(d.highlight, points(pos[col_pos[[1]]], highlight_LOD[col_pos[[1]]], col=peach_colors[1], type="h",lty=2,...))
  }else{
    with(d.highlight, points(pos[col_pos[[1]]], highlight_LOD[col_pos[[1]]], col=peach_colors[1], pch=20,...))
    with(d.highlight, points(pos[col_pos[[1]]], highlight_LOD[col_pos[[1]]], col=peach_colors[1], type="h",lty=2,...))
  }
  
  if(MarkGene==TRUE){
    text(c(Pos_x),c(Pos_y),c(GeneName),font=3,cex=LabelSize-0.3,col=GeneNameColour,...)
    lod_information<-cbind(d.highlight[,6],d.highlight[,8])
    colnames(lod_information)<-c("Genome position (Pos_x)","LOD (Pos_y)")
    write.table(lod_information,paste(dirplot,"/","Reference information to mark gene.csv",sep=""),sep=",",row.names = F)
  }
   dev.off()
}





