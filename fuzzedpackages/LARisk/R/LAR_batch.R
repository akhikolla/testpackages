#'Estimate LAR with Multiple IDs
#'
#'\code{LAR_batch} is used to estimate lifetime attributable radiation-related cancer risk for data with Multiple IDs.
#'
#'@param data  data frame containing 'id', 'sex', 'birth', 'dosedist', 'dose1', 'dose2', 'dose3', 'site', 'exposure_rate'.
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
#'@details
#'Basically, the arguments of \code{LAR_batch} are same as \code{LAR}.
#'But unlike \code{LAR}, \code{LAR_batch} can have multiple people's information.
#'
#'@return \code{LAR_batch} return a list of values for each IDs.
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
#'@return \code{LAR_batch} return an list contains 'LAR' and 'Future_LAR' for each person.
#'
#'@references Berrington de Gonzalez, A., Iulian Apostoaei, A., Veiga, L., Rajaraman, P., Thomas, B., Owen Hoffman, F., Gilbert, E. and Land, C. (2012). RadRAT: a radiation risk assessment tool for lifetime cancer risk projection. \emph{Journal of Radiological Protection}, \bold{32(3)}, pp.205-222.
#'@references National Research Council (NRC) and Committee to Assess Health Risks from Exposure to Low Levels of Ionizing Radiation (2005) \emph{Health Risks from Exposure to Low Levels of Ionizing Radiation: BEIR VII Phase 2} (Washington, DC: National Academy of Sciences)
#'
#'@examples
#'data<-data.frame(id=c("a100","a101"), birth=rep(1970,2), exposure=c(1980,1990),
#'                  dosedist=c("fixedvalue","fixedvalue"), dose1=c(10,20), dose2=c(0,0),
#'                  dose3=c(0,0), sex=rep("male",2), site=c("colon","leukemia"),
#'                  exposure_rate=rep("acute",2))
#'
#'LAR_batch(data)
#'
#'@export
#'@importFrom Rcpp evalCpp
#'@importFrom stats quantile
#'@importFrom utils write.csv
#'@useDynLib LARisk, .registration = TRUE


LAR_batch<-function(data, excel=FALSE, filename=NULL, weight_site="no", weight_value=0,
                    current=as.numeric(substr(Sys.time(),1,4)), sim=300, seed=99, basepy=100000, DDREF=TRUE,
                    ci=0.9, changedata=FALSE, dbaseline=0, dincidence=0, rounddigit=4){

  data<-data[data$dose1!=0,]

  data$sex<-tolower(data$sex)
  data$site<-tolower(data$site)
  data$exposure_rate<-tolower(data$exposure_rate)

  ## data check for dbaseline, dincidence ####
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


  ## excel==TRUE ####
  if(excel==TRUE){

    unique_id<-unique(data$id)
    result_colname<-c("ID","L_Tot","L_Solid","L_Oral","L_Esophagus","L_Stomach","L_Colon","L_Rectum",
                      "L_Liver","L_Gallbladder","L_Pancreas","L_Lung","L_Bladder","L_Kidney","L_Brain/CNS",
                      "L_Thyroid","L_Other","L_Leukemia","L_Prostate","L_Breast","L_Ovary","L_Uterus",
                      "M_Tot","M_Solid","M_Oral","M_Oesophagus","M_Stomach","M_Colon","M_Rectum",
                      "M_Liver","M_Gallbladder","M_Pancreas","M_Lung","M_Bladder","M_Kidney","M_Brain/CNS",
                      "M_Thyroid","M_Other","M_Leukemia","M_Prostate","M_Breast","M_Ovary","M_Uterus",
                      "U_Tot","U_Solid","U_Oral","U_Oesophagus","U_Stomach","U_Colon","U_Rectum",
                      "U_Liver","U_Gallbladder","U_Pancreas","U_Lung","U_Bladder","U_Kidney","U_Brain/CNS",
                      "U_Thyroid","U_Other","U_Leukemia","U_Prostate","U_Breast","U_Ovary","U_Uterus",
                      "F_L_Tot","F_L_Solid","F_L_Oral","F_L_Oesophagus","F_L_Stomach","F_L_Colon","F_L_Rectum",
                      "F_L_Liver","F_L_Gallbladder","F_L_Pancreas","F_L_Lung","F_L_Bladder","F_L_Kidney","F_L_Brain/CNS",
                      "F_L_Thyroid","F_L_Other","F_L_Leukemia","F_L_Prostate","F_L_Breast","F_L_Ovary","F_L_Uterus",
                      "F_M_Tot","F_M_Solid","F_M_Oral","F_M_Oesophagus","F_M_Stomach","F_M_Colon","F_M_Rectum",
                      "F_M_Liver","F_M_Gallbladder","F_M_Pancreas","F_M_Lung","F_M_Bladder","F_M_Kidney","F_M_Brain/CNS",
                      "F_M_Thyroid","F_M_Other","F_M_Leukemia","F_M_Prostate","F_M_Breast","F_M_Ovary","F_M_Uterus",
                      "F_U_Tot","F_U_Solid","F_U_Oral","F_U_Oesophagus","F_U_Stomach","F_U_Colon","F_U_Rectum",
                      "F_U_Liver","F_U_Gallbladder","F_U_Pancreas","F_U_Lung","F_U_Bladder","F_U_Kidney","F_U_Brain/CNS",
                      "F_U_Thyroid","F_U_Other","F_U_Leukemia","F_U_Prostate","F_U_Breast","F_U_Ovary","F_U_Uterus",
                      "LBR_Tot","LBR_Solid","LBR_Oral","LBR_Oesophagus","LBR_Stomach","LBR_Colon","LBR_Rectum",
                      "LBR_Liver","LBR_Gallbladder","LBR_Pancreas","LBR_Lung","LBR_Bladder","LBR_Kidney","LBR_Brain/CNS",
                      "LBR_Thyroid","LBR_Other","LBR_Leukemia","LBR_Prostate","LBR_Breast","LBR_Ovary","LBR_Uterus",
                      "BFR_Tot","BFR_Solid","BFR_Oral","BFR_Oesophagus","BFR_Stomach","BFR_Colon","BFR_Rectum",
                      "BFR_Liver","BFR_Gallbladder","BFR_Pancreas","BFR_Lung","BFR_Bladder","BFR_Kidney","BFR_Brain/CNS",
                      "BFR_Thyroid","BFR_Other","BFR_Leukemia","BFR_Prostate","BFR_Breast","BFR_Ovary","BFR_Uterus",
                      "LFR_Tot","LFR_Solid","LFR_Oral","LFR_Oesophagus","LFR_Stomach","LFR_Colon","LFR_Rectum",
                      "LFR_Liver","LFR_Gallbladder","LFR_Pancreas","LFR_Lung","LFR_Bladder","LFR_Kidney","LFR_Brain/CNS",
                      "LFR_Thyroid","LFR_Other","LFR_Leukemia","LFR_Prostate","LFR_Breast","LFR_Ovary","LFR_Uterus",
                      "TFR_Tot","TFR_Solid","TFR_Oral","TFR_Oesophagus","TFR_Stomach","TFR_Colon","TFR_Rectum",
                      "TFR_Liver","TFR_Gallbladder","TFR_Pancreas","TFR_Lung","TFR_Bladder","TFR_Kidney","TFR_Brain/CNS",
                      "TFR_Thyroid","TFR_Other","TFR_Leukemia","TFR_Prostate","TFR_Breast","TFR_Ovary","TFR_Uterus")

    ## Declaration(result_mat) ####
    result_mat<-as.data.frame(matrix(0,nrow=length(unique_id),ncol=length(result_colname)))
    colnames(result_mat)<-result_colname
    result_mat$ID<-unique_id


    ## calculation LAR with each id ####
    for(m in 1:length(unique_id)){
      tempdata<-data[data$id==unique_id[m],]
      set.seed(seed)

      ## get the results from C ####
      SUMMAT<-F_SUMMAT<-LARMAT_leukemia<-F_LARMAT_leukemia<-c()

      for(k in 1:dim(tempdata)[[1]]){

        ## data check ####
        if(length(unique(tempdata$birth))!=1)
          stop("Pairs of 'birth' is inconsistent")
        if(!all(tempdata$birth<=tempdata$exposure))
          stop("'exposure' values are improper")


        ##	input data ####
        birth<-tempdata$birth[k]
        exposure<-tempdata$exposure[k]
        dosedist<-tolower(tempdata$dosedist[k])
        doseinfo<-as.numeric(c(tempdata$dose1[k],tempdata$dose2[k],tempdata$dose3[k]))
        sex<-tolower(tempdata$sex[k])
        site<-tolower(tempdata$site[k])
        exposure_rate<-tolower(tempdata$exposure_rate[k])
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
               'loguniform'={ dosedist1<-7 }	)

        #DDREF_op1
        if(DDREF==T){
          DDREF_op1<-1
        }else{
          DDREF_op1<-2
        }

        #site1
        switch(site,
               'stomach'={   	site1<-1 },
               'colon'={  		site1<-2 },
               'liver'={ 	  	site1<-3 },
               'lung'={ 	  	site1<-4 },
               'breast'={ 		site1<-5 },
               'ovary'={ 		  site1<-6 },
               'uterus'={ 		site1<-7 },
               'prostate'={ 	site1<-8 },
               'bladder'={ 	  site1<-9 },
               'brain/cns'={ 	site1<-10 },
               'thyroid'={ 	  site1<-11 },
               'remainder'={  site1<-12 },
               'oral'={ 	    site1<-13 },
               'oesophagus'={ site1<-14 },
               'rectum'={ 		site1<-15 },
               'gallbladder'={site1<-16 },
               'pancreas'={ 	site1<-17 },
               'kidney'={ 		site1<-18 },
               'leukemia'={ 	site1<-19 }   )

        #weight_site1
        weight_site1<-c()
        if(all(weight_site=="no")){
          weight_length<-0
          weight_site1<-100
        }else{
          weight_length<-length(weight_site)
          for(i in 1:weight_length){
            switch(weight_site[i],
                   'stomach'={ 	weight_site1[i]<-1 },
                   'colon'={ 		weight_site1[i]<-2 },
                   'liver'={ 		weight_site1[i]<-3 },
                   'lung'={ 		weight_site1[i]<-4 },
                   'breast'={ 		weight_site1[i]<-5 },
                   'ovary'={ 		weight_site1[i]<-6 },
                   'uterus'={ 		weight_site1[i]<-7 },
                   'prostate'={ 	weight_site1[i]<-8 },
                   'bladder'={ 	weight_site1[i]<-9 },
                   'brain/cns'={ 	weight_site1[i]<-10 },
                   'thyroid'={ 	weight_site1[i]<-11 },
                   'remainder'={ 	weight_site1[i]<-12 },
                   'oral'={ 		weight_site1[i]<-13 },
                   'oesophagus'={ 	weight_site1[i]<-14 },
                   'rectum'={ 		weight_site1[i]<-15 },
                   'gallbladder'={ 	weight_site1[i]<-16 },
                   'pancreas'={ 	weight_site1[i]<-17 },
                   'kidney'={ 		weight_site1[i]<-18 },
                   'leukemia'={ 	weight_site1[i]<-19 }   )
          }
        }

        ## calculation result from C ####
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
      }	# end k


      ## data cleanup for baseline risk calculation ####
      data2<-BFR<-LBR<-lar<-F_lar<-c()
      unique_site<-as.character(unique(tempdata$site))
      base_count<-length(unique_site)

      for(i in 1:base_count){
        data1<-subset(tempdata,tempdata$site==unique_site[i])
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
               'stomach'={  	site1<-1 },
               'colon'={ 	  	site1<-2 },
               'liver'={ 	  	site1<-3 },
               'lung'={ 	  	site1<-4 },
               'breast'={ 		site1<-5 },
               'ovary'={ 	  	site1<-6 },
               'uterus'={ 		site1<-7 },
               'prostate'={ 	site1<-8 },
               'bladder'={ 	  site1<-9 },
               'brain/cns'={ 	site1<-10 },
               'thyroid'={   	site1<-11 },
               'remainder'={ 	site1<-12 },
               'oral'={ 	  	site1<-13 },
               'oesophagus'={ site1<-14 },
               'rectum'={ 		site1<-15 },
               'gallbladder'={site1<-16 },
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
          site_colsum<-colSums(subset(SUMMAT,tempdata$site==unique_site[p]))
          lar<-rbind(lar,round(c(quantile(site_colsum,(1-ci)/2,na.rm=TRUE),mean(site_colsum),quantile(site_colsum,(1+ci)/2,na.rm=TRUE)),rounddigit))

          F_site_colsum<-colSums(subset(F_SUMMAT,tempdata$site==unique_site[p]))
          F_lar<-rbind(F_lar,round(c(quantile(F_site_colsum,(1-ci)/2,na.rm=TRUE),mean(F_site_colsum),quantile(F_site_colsum,(1+ci)/2,na.rm=TRUE)),rounddigit))
        }else{
          lar<-rbind(lar,round(colSums(LARMAT_leukemia),rounddigit))
          F_lar<-rbind(F_lar,round(colSums(F_LARMAT_leukemia),rounddigit))
        }
      }	# end p

      # leukemia
      if(sum(unique_site!="leukemia")==0){

        lar<-rbind(lar,c(0,0,0),round(colSums(LARMAT_leukemia),rounddigit))
        F_lar<-rbind(F_lar,c(0,0,0),round(colSums(F_LARMAT_leukemia),rounddigit))

      }else{

        solid_colsum<-colSums(subset(SUMMAT,tempdata$site!="leukemia"))
        total_colsum<-colSums(SUMMAT)
        lar<-rbind(lar,c(round(c(quantile(solid_colsum,(1-ci)/2,na.rm=TRUE),mean(solid_colsum),quantile(solid_colsum,(1+ci)/2,na.rm=TRUE)),rounddigit)),
                   c(round(c(quantile(total_colsum,(1-ci)/2,na.rm=TRUE),mean(total_colsum),quantile(total_colsum,(1+ci)/2,na.rm=TRUE)),rounddigit)))

        F_solid_colsum<-colSums(subset(F_SUMMAT,tempdata$site!="leukemia"))
        F_total_colsum<-colSums(F_SUMMAT)
        F_lar<-rbind(F_lar,c(round(c(quantile(F_solid_colsum,(1-ci)/2,na.rm=TRUE),mean(F_solid_colsum),quantile(F_solid_colsum,(1+ci)/2,na.rm=TRUE)),rounddigit)),
                     c(round(c(quantile(F_total_colsum,(1-ci)/2,na.rm=TRUE),mean(F_total_colsum),quantile(F_total_colsum,(1+ci)/2,na.rm=TRUE)),rounddigit)))

      }

      # solid
      solid<-unique_site!="leukemia"
      base<-rbind(cbind(LBR,BFR),c(sum(LBR[solid]),sum(BFR[solid])),c(sum(LBR),sum(BFR)))
      base<-cbind(base,LFR=round(lar[,2]/base[,1],rounddigit+2)*100,TFR=F_lar[,2]+base[,2])

      rownames(lar)<-rownames(F_lar)<-rownames(base)<-c(as.character(unique_site),"solid","total")
      colnames(lar)<-colnames(F_lar)<-c("Lower","Mean","Upper")
      base[base[,2]==0,3]<-0

      ## result_mat ####
      lar <- round(lar*basepy/(10^5),rounddigit)
      F_lar <- round(F_lar*basepy/(10^5),rounddigit)
      base <- round(base*basepy/(10^5),rounddigit)

      length_lar<-dim(lar)[[1]]
      result_mat[m,c(2,23,44,65,86,107,128,149,170,191)]<-c(lar[length_lar,],F_lar[length_lar,],base[length_lar,])
      result_mat[m,c(3,24,45,66,87,108,129,150,171,192)]<-c(lar[length_lar-1,],F_lar[length_lar-1,],base[length_lar-1,])

      for(n in 1:base_count){
        switch(unique_site[n],
               'oral'={ 	    	result_mat[m,c(4, 25,46,67,88, 109,130,151,172,193)]<-c(lar[n,],F_lar[n,],base[n,]) },
               'oesophagus'={ 	result_mat[m,c(5, 26,47,68,89, 110,131,152,173,194)]<-c(lar[n,],F_lar[n,],base[n,]) },
               'stomach'={ 	    result_mat[m,c(6, 27,48,69,90, 111,132,153,174,195)]<-c(lar[n,],F_lar[n,],base[n,]) },
               'colon'={   	  	result_mat[m,c(7, 28,49,70,91, 112,133,154,175,196)]<-c(lar[n,],F_lar[n,],base[n,]) },
               'rectum'={ 	  	result_mat[m,c(8, 29,50,71,92, 113,134,155,176,197)]<-c(lar[n,],F_lar[n,],base[n,]) },
               'liver'={ 	    	result_mat[m,c(9, 30,51,72,93, 114,135,156,177,198)]<-c(lar[n,],F_lar[n,],base[n,]) },
               'gallbladder'={ 	result_mat[m,c(10,31,52,73,94, 115,136,157,178,199)]<-c(lar[n,],F_lar[n,],base[n,]) },
               'pancreas'={ 	  result_mat[m,c(11,32,53,74,95, 116,137,158,179,200)]<-c(lar[n,],F_lar[n,],base[n,]) },
               'lung'={      		result_mat[m,c(12,33,54,75,96, 117,138,159,180,201)]<-c(lar[n,],F_lar[n,],base[n,]) },
               'bladder'={    	result_mat[m,c(13,34,55,76,97, 118,139,160,181,202)]<-c(lar[n,],F_lar[n,],base[n,]) },
               'kidney'={ 	  	result_mat[m,c(14,35,56,77,98, 119,140,161,182,203)]<-c(lar[n,],F_lar[n,],base[n,]) },
               'brain/cns'={  	result_mat[m,c(15,36,57,78,99, 120,141,162,183,204)]<-c(lar[n,],F_lar[n,],base[n,]) },
               'thyroid'={    	result_mat[m,c(16,37,58,79,100,121,142,163,184,205)]<-c(lar[n,],F_lar[n,],base[n,]) },
               'remainder'={   	result_mat[m,c(17,38,59,80,101,122,143,164,185,206)]<-c(lar[n,],F_lar[n,],base[n,]) },
               'leukemia'={  	result_mat[m,c(18,39,60,81,102,123,144,165,186,207)]<-c(lar[n,],F_lar[n,],base[n,]) },
               'prostate'={   	result_mat[m,c(19,40,61,82,103,124,145,166,187,208)]<-c(lar[n,],F_lar[n,],base[n,]) },
               'breast'={ 	  	result_mat[m,c(20,41,62,83,104,125,146,167,188,209)]<-c(lar[n,],F_lar[n,],base[n,]) },
               'ovary'={ 		    result_mat[m,c(21,42,63,84,105,126,147,168,189,210)]<-c(lar[n,],F_lar[n,],base[n,]) },
               'uterus'={ 	   	result_mat[m,c(22,43,64,85,106,127,148,169,190,211)]<-c(lar[n,],F_lar[n,],base[n,]) }
        )
      } # end n
    } # end m

    ## excel : TRUE ####
    return(write.csv(result_mat,file=filename,row.names=F))

  }else{
    ## excel : FALSE ####
    return(lapply(split(data,data$id),LARisk::LAR,sim=sim,weight_site=weight_site,weight_value=weight_value,seed=seed,basepy=basepy,DDREF=DDREF,ci=ci,changedata=changedata,dbaseline=dbaseline,dincidence=dincidence,rounddigit=rounddigit))
  }
}
