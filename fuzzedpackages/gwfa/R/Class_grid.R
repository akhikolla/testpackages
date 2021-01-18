#Carroyage est le nom de la classe definie
#les slots sont des variables typees
setClass("Grid",
         slots = list(cell_size = "numeric",bandwith= "numeric",radius="vector",q="vector",sample_size="numeric"),
         contains = "data.frame")


#fonction generique = utilisable pour plusieurs classes differentes 
setGeneric(
  #creation de fond de carte dont la projection a le code epsg defini
  name="creation_fond_de_carte",
  def=function(object,epsg){standardGeneric("creation_fond_de_carte")}
)

#on definit plus specifiquement la fonction creation_fond_de_carte pour un objet carroyage
setMethod(
  f="creation_fond_de_carte",
  signature="Grid",
  definition= function(object,epsg)
  {
    #fonction de creation de carreaux de cote taille_carreau
    carreau<-function(x,y,taille_carreau) {(Polygons(list(Polygon(cbind(c(0,taille_carreau,taille_carreau,0,0)+(x-(taille_carreau)/2),c(0,0,taille_carreau,taille_carreau,0)+(y-(taille_carreau)/2)))),paste(x,y,sep="_")))}
    
    #on applique la fonction carreau
    #MoreArgs contient les parametres supplementaires de la fonction carreau
    grille=mapply(carreau,object[,"x"],object[,"y"],MoreArgs=list(taille_carreau=object@cell_size))
    
    #on cree un spatial polygon avec un code epsg de projection defini
    grille_spat = SpatialPolygons((grille),proj4string=CRS(paste("+init=epsg:",epsg,sep="")))
    df=data.frame(object@.Data);names(df)=names(object)
    data=data.frame(ID=paste(df[,"x"],df[,"y"],sep="_"),df)
    
    #un SpatialPolygonsDataFrame est un SpatialPolygon auquel on attache une table d attributs
    return(SpatialPolygonsDataFrame(grille_spat, data, match.ID = "ID"))
  }
)

setMethod(
  `[`,
  signature=signature(x="Grid"),
  function(x, ...){
    # Save the original
    storedtdt <- x
    # Use the fact that x is a subclass to "data.frame"
    Nargs <- nargs()
    hasdrop <- "drop" %in% names(sys.call())
    if(Nargs==2) {
      tmpDF <- `[.data.frame`(x, i=TRUE, j=i, ..., drop=FALSE)
    } else if((Nargs==3 && hasdrop)) {
      tmpDF <- `[.data.frame`(x, i=TRUE, j=i, ..., drop)
    } else if(hasdrop) {
      tmpDF <- `[.data.frame`(x, i, j, ..., drop)
    } else {
      tmpDF <- `[.data.frame`(x, i, j, ...)
    }
    # Reintegrate the results
    if (inherits(x=tmpDF, what="data.frame")){
      for(sName in names(getSlots("data.frame"))){
        slot(storedtdt, sName) <- slot(tmpDF, sName)
      }
      return(storedtdt)
    } else {
      return(tmpDF)
    }
  })

setMethod(
  `[<-`,
  signature=signature(x="Grid"),
  function(x, ..., value){
    # Save the original
    storedtdt <- x
    # Use the fact that x is a subclass to "data.frame"
    Nargs <- nargs()
    if (any(!names(sys.call()) %in% c("", "i", "j", "value"))) {
      stop("extra arguments are not allowed")
    }
    tmpDF <- data.frame(x)
    if(Nargs==3) {
      if (missing(i)) i <- j
      tmpDF[i] <- value
    } else if(Nargs==4) {
      tmpDF[i, j] <- value
    }
    # Reintegrate the results
    for(sName in names(getSlots("data.frame"))){
      slot(storedtdt, sName) <- slot(tmpDF, sName)
    }
    return(storedtdt)
  })


#fonction pour transformer un lissage en un carroyage (fond de carte)



grid_to_spdf<-function(df,epsg,cell_size=NULL,bandwith=NULL,radius=NULL,q=NULL,sample_size=NULL)
{
  if (is.null(cell_size)){cell_size=df@cell_size}
  if (is.null(bandwith)){bandwith=df@bandwith}
  if (is.null(radius)){radius=df@radius}
  if (is.null(q)){q=df@q}
  if (is.null(sample_size)){sample_size=df@sample_size}
  
  xcol="x";ycol="y"
  car=df[,c(xcol,ycol,(liste_var=names(df)[!(names(df)%in% c(xcol,ycol))]))]
  names(car)=c("x","y",liste_var)
  return(creation_fond_de_carte(new(Class = "Grid",car,cell_size=cell_size,bandwith=bandwith,radius=radius,sample_size=sample_size),epsg))
}

#constructeur grand public

gwfa<-function(points,q=0,radius,bandwith,sample_size,cell_size)
{
  xmin=floor(min(points$x)/cell_size)*cell_size
  xmax=ceiling(max(points$x)/cell_size)*cell_size
  ymin=floor(min(points$y)/cell_size)*cell_size
  ymax=ceiling(max(points$y)/cell_size)*cell_size
  
  coord=merge(data.frame(x=seq(xmin,xmax,cell_size)),data.frame(y=seq(ymin,ymax,cell_size) ) )
  
  #coord=coord[1:50,]
  res_gwfa=gwfa_c(points[,1],points[,2],coord[,1],coord[,2],bandwith,sample_size,radius,q)
  res=data.frame(coord,res_gwfa)
  
  k=1
  namesQR=list()
  for (i in 1:length(q))
  {
    for (j in 1:length(radius))
    {
      namesQR[k]=paste0("M_Q_",q[i],"_R_",round(radius[j]))
      k=k+1
    }
  }
  namesQR=do.call(c,namesQR)
  
  names(res)=c("x","y","count",namesQR)
  
  res=new(Class = "Grid",res,cell_size=cell_size,bandwith=bandwith,radius=radius,q=q,sample_size=sample_size)
  
  return(res)  
}


