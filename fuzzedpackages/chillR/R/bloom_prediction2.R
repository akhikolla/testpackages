#' Bloom prediction from chilling and forcing requirements, assumed to be
#' fulfilled strictly in sequence - version 2
#' 
#' This is a pretty rudimentary function to predict phenological dates from
#' chilling and forcing requirements and hourly chilling and forcing data. Note
#' that there are enormous uncertainties in these predictions, which are hardly
#' ever acknowledged. So please use this function with caution.
#' 
#' This function is an update to the bloom_prediction function, which was
#' quite slow and didn't allow testing multiple chilling and heat requirements.
#' In this updated version, chilling and heat requirements can be supplied
#' as vectors, which are interpreted in sequence, with each pair of Chill_req
#' and Heat_req values matched according to their position in the vectors.
#' Through the permutations argument, it is also possible to compute stage
#' occurrence dates for all possible combinations of the requirements
#' specified by the Chill_req and Heat_req vectors.
#' 
#' The model allows specifying any numeric column as the chill and heat columns,
#' indicated by the Chill_model and Heat_model parameters.
#' 
#' @param HourChillTable a data frame resulting from the chilling_hourtable
#' function.
#' @param Chill_req numeric vector indicating one or multiple chilling
#' requirements of the particular growth stage (in the unit specified by
#' "Chill_model")
#' @param Heat_req numeric vector indicating one or multiple heat requirements
#' of the particular growth stage (in Growing Degree Hours)
#' @param permutations boolean parameter indicating whether all possible
#' combinations of the supplied chilling and heat requirements should be used.
#' Defaults to FALSE, which means that the function matches chilling and heat
#' requirements according to their positions in the Chill_req and Heat_req
#' vectors and only predicts stage occurrence dates for these combinations.
#' @param Chill_model character string specifying the chill model to use. This
#' has to correspond to the name of the column in HourChillTable that contains
#' the chill accumulation (default is "Chill_Portions" for units of the Dynamic
#' Model).
#' @param Heat_model character string specifying the heat model to use. This
#' has to correspond to the name of the column in HourChillTable that contains
#' the heat accumulation (e.g "GDH").
#' @param Start_JDay numeric parameter indicating the day when chill
#' accumulation is supposed to start
#' @param infocol a vector of length length(Chill_req) which contains additional
#' information for each element of the vector. This is preserved and included
#' in the output. This only works when permutation=FALSE, and is meant to
#' facilitate recognition of particular phenological events in the output. 
#' @return data frame containing the predicted Julian dates of chilling requirement
#' fulfillment and timing of the phenological stage. Columns are Season, Creq, Hreq,
#' Creq_full (day when the chilling requirement is fulfilled) and Pheno_date
#' (the predicted date of the phenological event).
#' @author Eike Luedeling
#' @references Model references:
#' 
#' Dynamic Model:
#' 
#' Erez A, Fishman S, Linsley-Noakes GC, Allan P (1990) The dynamic model for
#' rest completion in peach buds. Acta Hortic 276, 165-174
#' 
#' Fishman S, Erez A, Couvillon GA (1987a) The temperature dependence of
#' dormancy breaking in plants - computer simulation of processes studied under
#' controlled temperatures. J Theor Biol 126(3), 309-321
#' 
#' Fishman S, Erez A, Couvillon GA (1987b) The temperature dependence of
#' dormancy breaking in plants - mathematical analysis of a two-step model
#' involving a cooperative transition. J Theor Biol 124(4), 473-483
#' 
#' Growing Degree Hours:
#' 
#' Anderson JL, Richardson EA, Kesner CD (1986) Validation of chill unit and
#' flower bud phenology models for 'Montmorency' sour cherry. Acta Hortic 184,
#' 71-78
#' 
#' 
#' @keywords bloom prediction
#' @examples
#' 
#' 
#' hourtemps<-stack_hourly_temps(fix_weather(KA_weather[which(KA_weather$Year>2008),]),latitude=50.4)
#' 
#' CT<-chilling_hourtable(hourtemps,Start_JDay=305)
#' 
#' bloom_prediction2(CT,c(30,40,50),c(1000,1500,2000))
#' bloom_prediction2(CT,c(30,40,50),c(1000,1500,2000),permutations=TRUE)
#' 
#' 
#' @export bloom_prediction2
bloom_prediction2 <-
  function (HourChillTable,Chill_req,Heat_req,permutations=FALSE,
            Chill_model="Chill_Portions",Heat_model="GDH",Start_JDay=305,infocol=NULL) 
  {
    
    if(is.null(HourChillTable)) stop("no HourChillTable provided.",call. = FALSE)
    if(!is.data.frame(HourChillTable)) stop("HourChillTable is not a data.frame.",call. = FALSE)
    if(is.null(Chill_model)) stop("no Chill_model provided.",call. = FALSE)
    if(!Chill_model %in% colnames(HourChillTable)) stop("HourChillTable doesn't contain a column named ",Chill_model,call. = FALSE)
    if(!Heat_model %in% colnames(HourChillTable)) stop("HourChillTable doesn't contain a column named ",Heat_model,call. = FALSE)
    if(!is.numeric(HourChillTable[,Chill_model])) stop("column ",Chill_model," not numeric.",call. = FALSE)
    if(!is.numeric(HourChillTable[,Heat_model])) stop("column ",Heat_model," not numeric.",call. = FALSE)
    if(is.null(Chill_req)) stop("no Chill_req provided.",call. = FALSE)
    if(!is.numeric(Chill_req)) stop("Chill_req not numeric.",call. = FALSE)
    if(is.null(Heat_req)) stop("no Heat_req provided.",call. = FALSE)
    if(!is.numeric(Heat_req)) stop("Heat_req not numeric.",call. = FALSE)
    if(is.null(Start_JDay)) stop("no Start_JDay provided.",call. = FALSE)
    if(Start_JDay>366) stop("Start_JDay can't be greater than 366",call. = FALSE)
    if(Start_JDay<1) stop("Start_JDay can't be less than 1",call. = FALSE)
    if(!"JDay" %in% colnames(HourChillTable))
      HourChillTable<-make_JDay(HourChillTable)
    if(!"Season" %in% colnames(HourChillTable))
    {HourChillTable[which(HourChillTable$JDay>=Start_JDay),"Season"]<-
      HourChillTable[which(HourChillTable$JDay>=Start_JDay),"Year"]
    HourChillTable[which(HourChillTable$JDay<Start_JDay),"Season"]<-
      HourChillTable[which(HourChillTable$JDay<Start_JDay),"Year"]-1}
    
    cc<-HourChillTable[,Chill_model]
    hh<-HourChillTable[,Heat_model]
    sea<-HourChillTable$Season
    stdd<-HourChillTable$JDay
    for (s in unique(sea))
      cc[which(sea==s)]<-cc[which(sea==s)]-cc[which(sea==s&stdd==round(Start_JDay))][1]
    
    for (s in unique(sea))
      hh[which(sea==s)]<-hh[which(sea==s)]-hh[which(sea==s&stdd==round(Start_JDay))][1]
    
    if(permutations)
      reqs<-expand.grid(unique(sea),Chill_req,Heat_req) else
      {if(!length(Chill_req)==length(Heat_req))
        stop("Chill_req and Heat_req are of different length.",call. = FALSE)
        reqs<-data.frame(Chill_req,Heat_req)
        reqs<-data.frame(Season=rep(unique(sea),each=nrow(reqs)),reqs)}
    colnames(reqs)<-c("Season","Creq","Hreq")
    if(!permutations & !is.null(infocol))
      if(length(infocol)==length(Chill_req)) reqs<-cbind(infocol,reqs)
    
    results <- data.frame()
    chill <- cc
    chill2 <- chill[c(2:length(chill), 1)]
    heat<-hh
    heat2 <- heat[c(2:length(heat), 1)]
    seas<-HourChillTable$Season
    
    unireq<-unique(reqs[,c("Creq","Season")])
    
    unireq[,"Chill_comp"]<-as.numeric(sapply(1:nrow(unireq),function(x)
      which(chill<unireq$Creq[x]&chill2>=unireq$Creq[x]&seas==unireq$Season[x])[1]))
    unireq[which(!is.na(unireq$Chill_comp)),"Heat_on_CR"]<-
      as.numeric(sapply(unireq$Chill_comp[which(!is.na(unireq$Chill_comp))],function(x)
        HourChillTable[x,Heat_model]))
    unireq[which(!is.na(unireq$Chill_comp)),"Chill_comp"]<-
      as.numeric(sapply(unireq$Chill_comp[which(!is.na(unireq$Chill_comp))],function(x)
        HourChillTable$JDay[x]))
    
    reqs<-merge(reqs,unireq,by = c("Creq","Season"))
    reqs[,"Heat_on_stage"]<-reqs$Hreq+reqs$Heat_on_CR
    
    heatstage<-reqs$Heat_on_stage
    reqseas<-reqs$Season
    reqs[,"Heat_comp"]<-as.numeric(sapply(1:nrow(reqs),function(x)
      which(heat<heatstage[x]&heat2>=heatstage[x]&seas==reqseas[x])[1]))
    reqs[which(!is.na(reqs$Heat_comp)),"Pheno_date"]<-
      as.numeric(sapply(reqs$Heat_comp[which(!is.na(reqs$Heat_comp))],function(x)
        HourChillTable$JDay[x]))
    
    if("infocol" %in% colnames(reqs))
      reqs<-reqs[,c("infocol","Season","Creq","Hreq","Chill_comp","Pheno_date")] else
        reqs<-reqs[,c("Season","Creq","Hreq","Chill_comp","Pheno_date")]
    
    return(reqs)
  }
