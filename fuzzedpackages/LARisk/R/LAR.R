#'Estimate LAR with Single ID
#'
#'\code{LAR} is used to estimate lifetime attributable radiation-related cancer risk for data with single ID.
#'
#'@param data data frame containing 'id', 'sex', 'birth', 'dosedist', 'dose1', 'dose2', 'dose3', 'site', 'exposure_rate'. See 'Details'.
#'@param weight_site vector containing the name of cancer sites to give weights.
#'@param weight_value numeric vector containing the value between 0 and 1 which is a weight on ERR model.
#'@param current number of current year. default is year of the system time.
#'@param sim number of iteration of simulation.
#'@param seed number of seed.
#'@param basepy number of base person year
#'@param DDREF logical. Whether to apply the dose and dose-rate effectiveness factor.
#'@param excel logical. Whether to extract the result as csv file.
#'@param filename a string naming the file to save (.csv file).
#'@param ci confidence level of the confidence interval.
#'@param changedata logical, whether to change the data of lifetime table and incidence rate.
#'@param dbaseline a path or data frame of the new lifetime table.
#'@param dincidence a path or data frame of the new incidence rate table.
#'@param rounddigit the number of decimal points to print.
#'
#'@details The data to be put in \code{LAR} should include some prerequisite information, which includes id, sex and birth of people(or person), distribution of dose, doses of exposed radiation, sites where exposed, and exposure rate.
#'Also, the variable names should be written as expressed.
#'The maximum age in \code{LAR} is set as 100 years old. If the data contains a birth year which makes attained age(= current - birth) over 100, the result has no useful value.
#'For some variables, there is a fixed format. \code{sex} can have the component 'male' or 'female'.
#'\code{dosedist} can have the component 'fixedvalue', 'normal', 'lognormal', 'triangular', 'logtriangular', 'uniform', 'loguniform'.
#'\code{site} can have the component 'stomach', 'colon', 'liver', 'lung', 'breast', 'ovary', 'uterus', 'prostate', 'bladder', 'brain/cns',
#''thyroid', 'remainder', 'oral', 'oesophagus', 'rectum', 'gallbladder', 'pancreas', 'kidney', 'leukemia'.
#'\code{exposure_rate} can have the component 'acute' or 'chronic'.
#'
#'@return \code{LAR}
#'@return     Cancer incidence probability per 100,000 persons to radiation exposure for their lifetime after exposed  year.
#'@return \code{LBR}
#'@return     Lifetime baseline risk. Cumulative baseline probability of having cancer over the maximum lifetime without radiation exposure after  exposed year.
#'@return \code{LFR}
#'@return     Lifetime fractional risk. Ratio LAR/LBR.
#'@return \code{Future_LAR}
#'@return     LAR after current year.
#'@return \code{BFR}
#'@return     Baseline future risk. Cumulative baseline probability of having cancer over the maximum lifetime without radiation exposure after  current year.
#'@return \code{TFR}
#'@return     Total future risk. Future LAR + BFR
#'
#'@references Berrington de Gonzalez, A., Iulian Apostoaei, A., Veiga, L., Rajaraman, P., Thomas, B., Owen Hoffman, F., Gilbert, E. and Land, C. (2012). RadRAT: a radiation risk assessment tool for lifetime cancer risk projection. \emph{Journal of Radiological Protection}, \bold{32(3)}, pp.205-222.
#'@references National Research Council (NRC) and Committee to Assess Health Risks from Exposure to Low Levels of Ionizing Radiation (2005) \emph{Health Risks from Exposure to Low Levels of Ionizing Radiation: BEIR VII Phase 2} (Washington, DC: National Academy of Sciences)
#'
#'@examples
#'data<-data.frame(id='a001', birth=1970, exposure=1980, dosedist='fixedvalue',
#'                 dose1=10, dose2=0, dose3=0, sex='male', site='colon',
#'                 exposure_rate='acute')
#'
#'LAR(data)
#'@export
#'@importFrom Rcpp evalCpp
#'@importFrom stats quantile
#'@importFrom utils write.csv
#'@useDynLib LARisk, .registration = TRUE

LAR <- function(data=data, weight_site="no", weight_value=0, current=as.numeric(substr(Sys.time(),1,4)), sim=300,
                seed=99, basepy=100000, DDREF=TRUE, excel=FALSE, filename=NULL, ci=0.9,
                changedata=FALSE, dbaseline=0, dincidence=0, rounddigit=4){
  set.seed(seed)

  SUMMAT <- F_SUMMAT <- LARMAT_leukemia <- F_LARMAT_leukemia <- c()

  data <- data[data$dose1!=0,]  # delete dose1=0

  data$sex <- tolower(data$sex)
  data$site <- tolower(data$site)
  data$exposure_rate <- tolower(data$exposure_rate)

  ## data check ####
  if(all(excel,is.null(filename)))
    stop("'filename' must be specified")

  if(changedata==TRUE){
    if(length(dbaseline)==1 | length(dincidence)==1)
      stop("Put the data")
    if(dim(dbaseline)[1]!=101 | dim(dbaseline)[2]!=3)
      stop("Put the baseline data in the correct format")
    if(dim(dincidence)[1]!=1919 | dim(dincidence)[2]!=4)
      stop("Put the incidence data in the correct format")
  }

  if(!all(data$sex %in% c("male", "female")))
    stop("'sex' hase invalid component.")

  if(!all(data$site %in% c("stomach","colon","liver","lung","breast","ovary","uterus","prostate","bladder","brain/cns",
                          "thyroid","remainder","oral","oesophagus","rectum","gallbladder","pancreas","kidney","leukemia")))
    stop("'site' has invalid component.")

  if(!all(data$exposure_rate %in% c("acute","chronic")))
    stop("'exposure_rate' has invalid component.")

  if(!all(data$dosedist %in% c( 'fixedvalue', 'normal', 'lognormal', 'triangular', 'logtriangular', 'uniform', 'loguniform')))
    stop("'dosedist' has invalid component.")


  if(any((data$sex=="male")&(data$site=="uterus" | data$site=="breast"| data$site=="ovary")))
    stop("'male' cannot calculate about uterus, breast or ovary")

  if(any(data$sex=="female"&data$site=="prostate"))
    stop("'female' cannot calculate about prostate")

  if(length(unique(data$birth))!=1)
    stop("Pairs of 'birth' is inconsistent")

  if(!all(data$birth<=data$exposure))
    stop("'exposure' values are improper")


  ## changedata
  if(changedata==FALSE){
    changedata1<-2
    dbaseline<-0
    dincidence<-0
  }else{
    changedata1<-1
    dbaseline<-c(dbaseline$Prob_d_m,dbaseline$Prob_d_m)
    dincidence<-c(dincidence$Rate_m,dincidence$Rate_f)
  }


  ## get the results from C ####
  for(k in 1:dim(data)[1]){

    ## input data ####
    birth<-data$birth[k]
    exposure<-data$exposure[k]
    dosedist<-tolower(data$dosedist[k])
    doseinfo<-as.numeric(c(data$dose1[k],data$dose2[k],data$dose3[k]))
    sex<-tolower(data$sex[k])
    site<-tolower(data$site[k])
    exposure_rate<-tolower(data$exposure_rate[k])
    age<-current-birth
    exposeage<-exposure-birth

    #sex1
    if(sex=="male"){
      sex1<-1
    }else{
      sex1<-2
    }

    #exposure_rate1
    if(exposure_rate=="acute"){
      exposure_rate1<-1
    }else{
      exposure_rate1<-2
    }

    #dosedist1
    switch(dosedist,
           'fixedvalue'={ dosedist1<-1 },
           'lognormal'={ dosedist1<-2 },
           'normal'={ dosedist1<-3 },
           'triangular'={ dosedist1<-4 },
           'logtriangular'={ dosedist1<-5 },
           'uniform'={	dosedist1<-6 },
           'loguniform'={ dosedist1<-7 })

    #DDREF_op1
    if(DDREF==T){
      DDREF_op1<-1
    }else{
      DDREF_op1<-2
    }

    #site1
    switch(site,
           'stomach'={ 	  site1<-1 },
           'colon'={ 		  site1<-2 },
           'liver'={ 		  site1<-3 },
           'lung'={ 		  site1<-4 },
           'breast'={ 	  site1<-5 },
           'ovary'={ 		  site1<-6 },
           'uterus'={ 	  site1<-7 },
           'prostate'={	  site1<-8 },
           'bladder'={ 	  site1<-9 },
           'brain/cns'={  site1<-10 },
           'thyroid'={ 	  site1<-11 },
           'remainder'={  site1<-12 },
           'oral'={ 		  site1<-13 },
           'oesophagus'={ site1<-14 },
           'rectum'={ 		site1<-15 },
           'gallbladder'={site1<-16 },
           'pancreas'={ 	site1<-17 },
           'kidney'={ 		site1<-18 },
           'leukemia'={ 	site1<-19 })

    #weight_size1
    weight_site1<-c()
    if(all(weight_site=="no")){
      weight_length<-0
      weight_site1<-100
    }else{
      weight_length<-length(weight_site)
      for(i in 1:weight_length){
        switch(weight_site[i],
               'stomach'={   	weight_site1[i]<-1 },
               'colon'={ 	  	weight_site1[i]<-2 },
               'liver'={ 	  	weight_site1[i]<-3 },
               'lung'={ 	  	weight_site1[i]<-4 },
               'breast'={   	weight_site1[i]<-5 },
               'ovary'={ 	  	weight_site1[i]<-6 },
               'uterus'={ 		weight_site1[i]<-7 },
               'prostate'={ 	weight_site1[i]<-8 },
               'bladder'={  	weight_site1[i]<-9 },
               'brain/cns'={ 	weight_site1[i]<-10 },
               'thyroid'={  	weight_site1[i]<-11 },
               'remainder'={ 	weight_site1[i]<-12 },
               'oral'={ 	  	weight_site1[i]<-13 },
               'oesophagus'={ weight_site1[i]<-14 },
               'rectum'={ 		weight_site1[i]<-15 },
               'gallbladder'={weight_site1[i]<-16 },
               'pancreas'={ 	weight_site1[i]<-17 },
               'kidney'={ 		weight_site1[i]<-18 },
               'leukemia'={ 	weight_site1[i]<-19 }   )
      }
    }

    ## calculate results from C ####
    out<-vector(mode="numeric",length=2*sim+6)

    result<-.C("larft",
               as.integer(dosedist1),
               as.double(doseinfo),
               as.integer(sex1),
               as.integer(site1),
               as.integer(exposure_rate1),
               as.integer(age),
               as.double(exposeage),
               as.integer(DDREF_op1),
               as.integer(sim),
               as.double(ci),
               as.integer(weight_site1),
               as.double(weight_value),
               as.integer(weight_length),
               as.integer(changedata1),
               as.double(dbaseline),
               as.double(dincidence),
               out=as.double(out))$out

    SUMMAT<-rbind(SUMMAT,result[1:sim])
    F_SUMMAT<-rbind(F_SUMMAT,result[(sim+1):(2*sim)])

    if(site=="leukemia"){
      LARMAT_leukemia<-rbind(LARMAT_leukemia,result[(2*sim+1):(2*sim+3)])
      F_LARMAT_leukemia<-rbind(F_LARMAT_leukemia,result[(2*sim+4):(2*sim+6)])
    }
  } #end k

  ## data cleanup for baseline risk calculation ####
  data2<-BFR<-LBR<-lar<-F_lar<-c()
  unique_site<-unique(data$site)
  base_count<-length(unique_site)

  for(i in 1:base_count){
    data1<-subset(data,data$site==unique_site[i])
    data2<-rbind(data2,data1[data1$exposure==min(data1$exposure),])
  }
  data_base<-data2[!duplicated(data2$site),]

  ## cleanup results ####
  for(p in 1:base_count){

    birth<-data_base$birth[p]

    sex<-data_base$sex[p]
    if(sex=="male"){
      sex1=1
    }else{
      sex1=2
    }

    site<-as.character(data_base$site[p])
    switch(site,
           'stomach'={ 	site1<-1 },
           'colon'={ 		site1<-2 },
           'liver'={ 		site1<-3 },
           'lung'={ 		site1<-4 },
           'breast'={ 		site1<-5 },
           'ovary'={ 		site1<-6 },
           'uterus'={ 		site1<-7 },
           'prostate'={ 	site1<-8 },
           'bladder'={ 	site1<-9 },
           'brain/cns'={ 	site1<-10 },
           'thyroid'={ 	site1<-11 },
           'remainder'={ 	site1<-12 },
           'oral'={ 		site1<-13 },
           'oesophagus'={ 	site1<-14 },
           'rectum'={ 		site1<-15 },
           'gallbladder'={ 	site1<-16 },
           'pancreas'={ 	site1<-17 },
           'kidney'={ 		site1<-18 },
           'leukemia'={ 	site1<-19 }   )

    out1<-0
    BFR[p]<-round(.C("brft",
                          as.integer(sex1),
                          as.integer(site1),
                          as.integer(current-birth),
                          as.integer(changedata1),
                          as.double(dbaseline),
                          as.double(dincidence),
                          out1=as.double(out1))$out1,rounddigit)
    out2<-0
    LBR[p]<-round(.C("brft",
                          as.integer(sex1),
                          as.integer(site1),
                          as.integer(data_base$exposure[p]-birth),
                          as.integer(changedata1),
                          as.double(dbaseline),
                          as.double(dincidence),
                          out2=as.double(out2))$out2,rounddigit)

    if(unique_site[p]!="leukemia"){
      site_colsum<-colSums(subset(SUMMAT,data$site==unique_site[p]))
      lar<-rbind(lar,round(c(quantile(site_colsum,(1-ci)/2,na.rm=TRUE),mean(site_colsum),quantile(site_colsum,(1+ci)/2,na.rm=TRUE)),rounddigit))

      F_site_colsum<-colSums(subset(F_SUMMAT,data$site==unique_site[p]))
      F_lar<-rbind(F_lar,round(c(quantile(F_site_colsum,(1-ci)/2,na.rm=TRUE),mean(F_site_colsum),quantile(F_site_colsum,(1+ci)/2,na.rm=TRUE)),rounddigit))
    }else{
      lar<-rbind(lar,round(colSums(LARMAT_leukemia),rounddigit))
      F_lar<-rbind(F_lar,round(colSums(F_LARMAT_leukemia),rounddigit))
    }
  } # end p

  #leukemia
  if(sum(unique_site!="leukemia")==0){

    lar<-rbind(lar,c(0,0,0),round(colSums(LARMAT_leukemia),rounddigit))
    F_lar<-rbind(F_lar,c(0,0,0),round(colSums(F_LARMAT_leukemia),rounddigit))

  }else{

    solid_colsum<-colSums(subset(SUMMAT,data$site!="leukemia"))
    total_colsum<-colSums(SUMMAT)
    lar<-rbind(lar,c(round(c(quantile(solid_colsum,(1-ci)/2,na.rm=TRUE),mean(solid_colsum),quantile(solid_colsum,(1+ci)/2,na.rm=TRUE)),rounddigit)),
               c(round(c(quantile(total_colsum,(1-ci)/2,na.rm=TRUE),mean(total_colsum),quantile(total_colsum,(1+ci)/2,na.rm=TRUE)),rounddigit)))

    F_solid_colsum<-colSums(subset(F_SUMMAT,data$site!="leukemia"))
    F_total_colsum<-colSums(F_SUMMAT)
    F_lar<-rbind(F_lar,c(round(c(quantile(F_solid_colsum,(1-ci)/2,na.rm=TRUE),mean(F_solid_colsum),quantile(F_solid_colsum,(1+ci)/2,na.rm=TRUE)),rounddigit)),
                 c(round(c(quantile(F_total_colsum,(1-ci)/2,na.rm=TRUE),mean(F_total_colsum),quantile(F_total_colsum,(1+ci)/2,na.rm=TRUE)),rounddigit)))

  }

  #solid
  solid<-unique_site!="leukemia"
  base<-rbind(cbind(LBR,BFR),c(sum(LBR[solid]),sum(BFR[solid])),c(sum(LBR),sum(BFR)))
  base<-cbind(base,LFR=round(lar[,2]/base[,1],rounddigit+2)*100,TFR=F_lar[,2]+base[,2])

  rownames(lar)<-rownames(F_lar)<-rownames(base)<-c(as.character(unique_site),"solid","total")
  colnames(lar)<-colnames(F_lar)<-c("Lower","Mean","Upper")
  base[base[,2]==0,3]<-0


  ## excel ####
  if(excel==FALSE){
    temp_lar <- cbind(lar,base[,c(1,3)])
    temp_F_lar <- cbind(F_lar,base[,c(2,4)])
    return(list(LAR=round(temp_lar*basepy/(10^5),rounddigit),Future_LAR=round(temp_F_lar*basepy/(10^5),rounddigit)))
  }else{
    excelresult<-round(cbind(lar*basepy/(10^5),F_lar*basepy/(10^5),base*basepy/(10^5)),rounddigit)
    colnames(excelresult)<-c("Lower","Mean","Upper","F_Lower","F_Mean","F_Upper","LBR","BFR","LFR","TFR")
    return(write.csv(excelresult,file=filename))
  }
}
