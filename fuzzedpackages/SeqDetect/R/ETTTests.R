synthetic_test_agenda <- function(label_aspect=1.0) {
  df <- data.frame(case_id=c("Case 1"),ts=c("26.07.2019. 10:30:20"),te=c("26.07.2019. 10:33:01"),activity=c("Task A"),
                   stringsAsFactors=FALSE)
  df[1,] <- list("Case 1","26.07.2019. 10:30:20","26.07.2019. 10:31:20","Task A")
  df[2,] <- list("Case 2","27.07.2019. 10:30:20","27.07.2019. 10:31:20","Task A")
  df[3,] <- list("Case 3","28.07.2019. 10:30:20","28.07.2019. 10:31:20","Task A")
  df[4,] <- list("Case 3","28.07.2019. 10:32:47","28.07.2019. 10:34:11","Task B")
  df[5,] <- list("Case 1","26.07.2019. 10:32:47","26.07.2019. 10:34:11","Task B")
  df[6,] <- list("Case 1","26.07.2019. 10:32:51","26.07.2019. 10:34:02","Task C")
  df[7,] <- list("Case 2","27.07.2019. 10:32:23","27.07.2019. 10:34:15","Task C")
  df[8,] <- list("Case 4","29.07.2019. 10:30:20","29.07.2019. 10:31:20","Task A")
  df[9,] <- list("Case 2","27.07.2019. 10:32:31","27.07.2019. 10:34:50","Task B")
  df[10,] <- list("Case 2","27.07.2019. 10:42:31","27.07.2019. 10:47:38","Task D")
  df[11,] <- list("Case 5","28.07.2019. 10:30:20","28.07.2019. 10:31:20","Task E")
  df[12,] <- list("Case 4","29.07.2019. 10:32:49","29.07.2019. 10:34:07","Task C")
  df[13,] <- list("Case 1","26.07.2019. 11:01:53","26.07.2019. 11:10:26","Task D")
  df[14,] <- list("Case 3","28.07.2019. 10:32:59","28.07.2019. 10:34:44","Task C")
  df[15,] <- list("Case 3","28.07.2019. 10:39:51","28.07.2019. 10:54:17","Task D")
  df[16,] <- list("Case 4","29.07.2019. 10:33:19","29.07.2019. 10:33:22","Task B")
  df[17,] <- list("Case 5","28.07.2019. 10:47:22","28.07.2019. 10:55:23","Task F")
  df[18,] <- list("Case 4","29.07.2019. 10:47:22","29.07.2019. 10:55:23","Task D")
  df[19,] <- list("Case 1","26.07.2019. 11:12:00","26.07.2019. 11:16:00","Task G")
  df[20,] <- list("Case 1","26.07.2019. 11:12:10","26.07.2019. 11:12:20","Task H")
  df[21,] <- list("Case 1","26.07.2019. 11:12:30","26.07.2019. 11:17:40","Task I")
  df[22,] <- list("Case 1","26.07.2019. 11:16:10","26.07.2019. 11:18:00","Task J")
  df[23,] <- list("Case 1","26.07.2019. 11:18:20","26.07.2019. 11:18:40","Task K")
  df[24,] <- list("Case 5","28.07.2019. 11:18:20","28.07.2019. 11:18:40","Task K")
  df <- transform(df,ts=as.POSIXct(df$ts,format="%d.%m.%Y. %H:%M:%S"),
                  te=as.POSIXct(df$te,format="%d.%m.%Y. %H:%M:%S"))
  
  hsc_pc_att <- HSC_PC_Attribute(field="activity")
  hsc <- HybridSequenceClassifier(c("case_id","ts","te","activity"),"ts","te","case_id",hsc_pc_att,
                                  pattern_field="activity")
  
  hsc$process(list(stream=df))
  hsc$mergeMachines()
  hsc$induceSubmachine(1,TRUE)
  hsc$printMachines()
  hsc$plotMachines(label_aspect=label_aspect)
}

sepsis_dataset_test <- function(induce_biomarker_decision_tree=TRUE,threshold=75,debug=FALSE,hsc=NULL) {
  sepsis <- eventdataR::sepsis
  if(!induce_biomarker_decision_tree) {
    df <- data.frame(sepsis,stringsAsFactors=FALSE)
    df <- transform(df,diagnosticartastrup=as.logical(df$diagnosticartastrup),age=as.numeric(df$age),crp=as.numeric(df$crp),
                    diagnosticblood=as.logical(df$diagnosticblood),diagnosticecg=as.logical(df$diagnosticecg),
                    diagnosticic=as.logical(df$diagnosticic),diagnosticlacticacid=as.logical(df$diagnosticlacticacid),
                    diagnosticliquor=as.logical(df$diagnosticliquor),diagnosticother=as.logical(df$diagnosticother),
                    diagnosticsputum=as.logical(df$diagnosticsputum),diagnosticurinaryculture=as.logical(df$diagnosticurinaryculture),
                    diagnosticurinarysediment=as.logical(df$diagnosticurinarysediment),diagnosticxthorax=as.logical(df$diagnosticxthorax),
                    disfuncorg=as.logical(df$disfuncorg),hypotensie=as.logical(df$hypotensie),hypoxie=as.logical(df$hypoxie),
                    infectionsuspected=as.logical(df$infectionsuspected),infusion=as.logical(df$infusion),lacticacid=as.numeric(df$lacticacid),
                    leucocytes=as.numeric(df$leucocytes),oligurie=as.logical(df$oligurie),sirscritheartrate=as.logical(df$sirscritheartrate),
                    sirscritleucos=as.logical(df$sirscritleucos),sirscrittachypnea=as.logical(df$sirscrittachypnea),
                    sirscrittemperature=as.logical(df$sirscrittemperature),sirscriteria2ormore=as.logical(df$sirscriteria2ormore),
                    timestamp=as.POSIXct(df$timestamp),activity=as.character(df$activity))
    df <- transform(df,
                    crp_2=case_when(df$crp<476.2 ~ "ELEVATED",df$crp>=476.2 ~ "SEPSIS",
                                    is.na(df$crp) ~ "NONE"),
                    leucocytes_2=case_when(df$leucocytes<4 ~ "LOW",df$leucocytes>=4 & df$leucocytes<=12 ~ "NORMAL",
                                           df$leucocytes>12 ~ "HIGH",is.na(df$leucocytes) ~ "NONE"),
                    lacticacid_2=case_when(df$lacticacid<0.5 ~ "LOW",df$lacticacid>=0.5 & df$lacticacid<=2.2 ~ "NORMAL",
                                           df$lacticacid>2.2 ~ "HIGH",is.na(df$lacticacid) ~ "NONE"),
                    biomarkers=NA
    )
  } else {
    df <- .generate_sepsis_subset()
  }
  df <- transform(df,
                  .clazz=case_when(df$activity=="ER Registration"~"REG",df$activity=="ER Triage"~"ERT",
                                   df$activity=="ER Sepsis Triage"~"ERST",df$activity=="IV Liquid"~"IVL",
                                   df$activity=="IV Antibiotics"~"IVA",df$activity=="Admission NC"~"ANC",
                                   df$activity=="Admission IC"~"AIC",df$activity=="Return ER"~"RER",df$activity=="Release A"~"REL_A",
                                   df$activity=="Release B"~"REL_B",df$activity=="Release C"~"REL_C",df$activity=="Release D"~"REL_D",
                                   df$activity=="Release E"~"REL_E",df$activity=="LacticAcid"~paste0("LA",df$lacticacid_2),
                                   df$activity=="CRP"~paste0("CRP",df$crp_2),df$activity=="Leucocytes"~paste0("WBC",df$leucocytes_2),
                                   df$activity=="Biomarker assessment"~paste0("BA",df$biomarkers)),
                  .pattern=paste0(as.character(df$activity),
                                  case_when(df$activity=="Leucocytes" ~ paste0("(",as.character(df$leucocytes_2),")"),
                                            df$activity=="LacticAcid" ~ paste0("(",as.character(df$lacticacid_2),")"),
                                            df$activity=="CRP" ~ paste0("(",as.character(df$crp_2),")"),
                                            df$activity=="Biomarker assessment" ~ paste0("(",as.character(df$biomarkers),")"),
                                            !df$activity %in% c("Leucocytes","LacticAcid","CRP","Biomarker assessment") ~ ""))
                )
  df <- df[order(df$timestamp),]
  if(is.null(hsc)) {
    hsc_pc_none <- HSC_PC_None()
    dd1 <- list(type="time",days=as.integer(5),context_related=FALSE)
    hsc <- HybridSequenceClassifier(c("case_id","timestamp","activity",".clazz",".pattern"),"timestamp","timestamp","case_id",
                                    hsc_pc_none,pattern_field=".pattern",decay_descriptors=list(d1=dd1))
    
    hsc$process(list(stream=df),debug=debug)
    hsc$mergeMachines()
    hsc$induceSubmachine(threshold=threshold,isolate=TRUE)
    hsc$printMachines(print_cache=FALSE,print_keys=FALSE)
    hsc$plotMachines()
  } else {
    hsc$process(list(stream=df),debug=debug,learn=FALSE)
    hsc$printMachines(print_cache=FALSE,print_keys=FALSE)
    hsc$plotMachines()
  }
}

.generate_sepsis_subset <- function() {
  sepsis <- eventdataR::sepsis
  df <- data.frame(sepsis,stringsAsFactors=FALSE)
  df <- transform(df,diagnosticartastrup=as.logical(df$diagnosticartastrup),age=as.numeric(df$age),crp=as.numeric(df$crp),
                  diagnosticblood=as.logical(df$diagnosticblood),diagnosticecg=as.logical(df$diagnosticecg),
                  diagnosticic=as.logical(df$diagnosticic),diagnosticlacticacid=as.logical(df$diagnosticlacticacid),
                  diagnosticliquor=as.logical(df$diagnosticliquor),diagnosticother=as.logical(df$diagnosticother),
                  diagnosticsputum=as.logical(df$diagnosticsputum),diagnosticurinaryculture=as.logical(df$diagnosticurinaryculture),
                  diagnosticurinarysediment=as.logical(df$diagnosticurinarysediment),diagnosticxthorax=as.logical(df$diagnosticxthorax),
                  disfuncorg=as.logical(df$disfuncorg),hypotensie=as.logical(df$hypotensie),hypoxie=as.logical(df$hypoxie),
                  infectionsuspected=as.logical(df$infectionsuspected),infusion=as.logical(df$infusion),lacticacid=as.numeric(df$lacticacid),
                  leucocytes=as.numeric(df$leucocytes),oligurie=as.logical(df$oligurie),sirscritheartrate=as.logical(df$sirscritheartrate),
                  sirscritleucos=as.logical(df$sirscritleucos),sirscrittachypnea=as.logical(df$sirscrittachypnea),
                  sirscrittemperature=as.logical(df$sirscrittemperature),sirscriteria2ormore=as.logical(df$sirscriteria2ormore),
                  timestamp=as.POSIXct(df$timestamp),activity=as.character(df$activity))
  df2 <- data.frame(df[!df$activity %in% c("CRP","LacticAcid","Leucocytes"),])
  df2 <- cbind(df2,data.frame(biomarkers=NA))
  df3 <- data.frame(df[df$activity %in% c("CRP","LacticAcid","Leucocytes"),])
  ts_array <- unique(df3$timestamp)
  df3 <- transform(df3,
                   crp_2=case_when(df3$crp<476.2 ~ "ELEVATED",df3$crp>=476.2 ~ "SEPSIS",
                                   is.na(df3$crp) ~ "NONE"),
                   leucocytes_2=case_when(df3$leucocytes<4 ~ "LOW",df3$leucocytes>=4 & df3$leucocytes<=12 ~ "NORMAL",df3$leucocytes>12 ~ "HIGH",
                                          is.na(df3$leucocytes) ~ "NONE"),
                   lacticacid_2=case_when(df3$lacticacid<0.5 ~ "LOW",df3$lacticacid>=0.5 & df3$lacticacid<=2.2 ~ "NORMAL",df3$lacticacid>2.2 ~ "HIGH",
                                          is.na(df3$lacticacid) ~ "NONE")
  )
  for(ts in ts_array) {
    df4 <- df3[df3$timestamp==ts,]
    df4 <- df4[order(df4$activity),]
    crp <- la <- le <- "NONE"
    for(i in 1:nrow(df4)) {
      element <- df4[i,]
      if(element$activity=="CRP") crp <- element$crp_2
      else if(element$activity=="LacticAcid") la <- element$lacticacid_2
      else if(element$activity=="Leucocytes") le <- element$leucocytes_2
    }
    df4[1,"activity"] <- "Biomarker assessment"
    df4 <- cbind(df4,data.frame(biomarkers=paste0("CRP=",crp,"#LA=",la,"#LE=",le)))
    df2 <- rbind(df2, df4[1,colnames(df2)])
  }
  return(df2)
}

bpi_challenge_2019_test1 <- function() {
  te1 <- environment()
  file <- system.file("datasets","BPI_Challenge_3WayEC_2019.RData",package="SeqDetect")
  load(envir=te1,file)
  df <- data.frame(te1$df_3way_EC,stringsAsFactors=FALSE)
  df <- transform(df,
                  CASE_Goods.Receipt=as.logical(df$CASE_Goods.Receipt),
                  activity_id=as.character(df$activity_id),
                  Cumulative.net.worth..EUR.=as.numeric(df$Cumulative.net.worth..EUR.),
                  resource_id=as.character(df$resource_id),
                  lifecycle_id=as.character(df$lifecycle_id),
                  activity_instance_id=as.integer(df$activity_instance_id),
                  CASE_GR.Based.Inv..Verif.=as.logical(df$CASE_GR.Based.Inv..Verif.),
                  .clazz=case_when(df$activity_id=="Vendor creates debit memo"~"VCDM",
                                   df$activity_id=="Vendor creates invoice"~"VCI",
                                   df$activity_id=="SRM: Created"~"SRM:CR",
                                   df$activity_id=="SRM: Complete"~"SRM:CO",
                                   df$activity_id=="SRM: Awaiting Approval"~"SRM:AA",
                                   df$activity_id=="SRM: Document Completed"~"SRM:DC",
                                   df$activity_id=="SRM: In Transfer to Execution Syst."~"SRM:ITE",
                                   df$activity_id=="SRM: Ordered"~"SRM:O",
                                   df$activity_id=="SRM: Change was Transmitted"~"SRM:CT",
                                   df$activity_id=="Create Purchase Order Item"~"CPOI",
                                   df$activity_id=="Record Goods Receipt"~"RGR",
                                   df$activity_id=="Record Invoice Receipt"~"RIR",
                                   df$activity_id=="Remove Payment Block"~"RPB",
                                   df$activity_id=="SRM: Deleted"~"SRM:DE",
                                   df$activity_id=="Delete Purchase Order Item"~"DPOI",
                                   df$activity_id=="SRM: Transaction Completed"~"SRM:TC",
                                   df$activity_id=="Cancel Goods Receipt"~"CGR",
                                   df$activity_id=="Change Price"~"CP",
                                   df$activity_id=="Clear Invoice"~"CI",
                                   df$activity_id=="Change Quantity"~"CHGQ",
                                   df$activity_id=="Cancel Invoice Receipt"~"CIR",
                                   df$activity_id=="Change Final Invoice Indicator"~"CFII",
                                   df$activity_id=="Change Delivery Indicator"~"CDI",
                                   df$activity_id=="SRM: Incomplete"~"SRM:INC",
                                   df$activity_id=="SRM: Held"~"SRM:HE",
                                   df$activity_id=="Cancel Subsequent Invoice"~"CSI",
                                   df$activity_id=="Record Subsequent Invoice"~"RSI"))
  df <- df[order(df$timestamp),]
  
  hsc_pc_none <- HSC_PC_None()
  hsc <- HybridSequenceClassifier(c("CASE_concept_name","activity_instance_id","activity_id",".clazz"),
                                  "activity_instance_id","activity_instance_id","CASE_concept_name",
                                  hsc_pc_none,pattern_field="activity_id")
  
  hsc$process(list(stream=df))
  hsc$mergeMachines()
  hsc$induceSubmachine(threshold=70,isolate=TRUE)
  hsc$plotMachines()
  hsc$printMachines(print_cache=FALSE,print_keys=FALSE)
}

prep_sales_dataset <- function(products=1:812) {
  file <- system.file("datasets","sales_transactions.csv",package="SeqDetect")
  df <- read.csv(file)
  cn <- colnames(df)
  cn <- cn[!grepl("W",cn) & cn!="MIN" & cn!="MAX"]
  df <- df[,cn]
  rdf <- data.frame(product=character(),quantity=numeric(),stringsAsFactors=FALSE)
  gen1 <- FALSE
  if(812 %in% products) gen1 <- TRUE
  if(is.null(products)) products <- 1:nrow(df)
  else products <- setdiff(products,c(812))
  for(i in 0:51) {
    week_code <- paste0("Normalized.",i)
    for(j in products) {
      rdf <- rbind(rdf,data.frame(product=df[j,"Product_Code"],quantity=as.double(df[j,week_code]),
                                  week=i,tstart=(i+1),tend=(i+1),stringsAsFactors=FALSE))
    }
  }
  if(gen1) rdf <- rbind(rdf,prep_sales_dataset_gen_P812(df))
  r <- list(original=df,transformed=rdf,products=unique(rdf[,"product"]))
  return(r)
}

prep_sales_dataset_gen_P812 <- function(original_df) {
  message("Generating synthetic PSynth1 based on P7")
  P7 <- original_df[7,]
  rdf <- data.frame()
  if(nrow(P7)>0) {
    for(i in 0:47) {
      week_code <- paste0("Normalized.",i)
      q <- as.double(P7[,week_code])
      q <- q - abs(rnorm(1,0,0.025))
      if(q>1) q <- 1
      if(q<0) q <- 0
      rdf <- rbind(rdf,data.frame(product="PSynth1",quantity=q,week=(i+4),tstart=(i+5),tend=(i+5),stringsAsFactors=FALSE))
    }
    for(i in 0:3) {
      rdf <- rbind(rdf,data.frame(product="PSynth1",quantity=as.double(0),week=i,tstart=(i+1),tend=(i+1),stringsAsFactors=FALSE))
    }
  }
  return(rdf)
}

sales_dataset_test <- function(learning_set=1:20,testing_set=21:40,th_increment=1,max_th=NULL) {
  df_learn <- prep_sales_dataset(learning_set)
  hsc_pc_bin <- HSC_PC_Binning(0,1,10,"quantity")
  hsc <- HybridSequenceClassifier(c("product","quantity","week","tstart","tend"),"tstart","tend","product",
                                  hsc_pc_bin,pattern_field=".clazz",time_series_sequence_stats=TRUE)
  hsc$process(list(stream=df_learn$transformed))
  hsc$mergeMachines()
  hsc$printMachines(NULL,NULL,TRUE,TRUE)
  
  df_test <- prep_sales_dataset(testing_set)
  
  tsummary <- data.frame()
  tsequences <- data.frame()
  tsequences_tmp <- data.frame()
  tsequences_product <- data.frame()
  ttotal <- data.frame()
  th <- tsl <- 1
  while(tsl>0 && (is.null(max_th) || th<max_th)) {
    message(paste0("Threshold:",th))
    temp <- hsc$clone()
    temp$induceSubmachine(threshold=th,isolate=TRUE)
    temp$cleanKeys()
    temp$plotMachines()
    #ofn <- paste0("/Users/dkrleza/Documents/sales_test_th",th,".csv")
    res <- temp$process(list(stream=df_test$transformed),learn=FALSE,threshold=th)
    
    summary <- data.frame()
    summary_tcl <- 0
    message(names(res$sequences))
    for(product in names(res$sequences)) {
      seq <- res$sequences[[product]]
      summary <- rbind(summary,data.frame(product=product,sequences=seq$chains,total_transitions=seq$total_chain_length))
      summary_tcl <- summary_tcl + seq$total_chain_length
    }
    summary <- cbind(summary,threshold=th)
    v1 <- v2 <- v3 <- 0
    for(i in 1:nrow(summary)) {
      v1 <- v1 + summary[i,"total_transitions",drop=T]
      v2 <- v2 + summary[i,"sequences",drop=T]
      t1 <- res$sequences
      t1 <- t1[[summary[i,"product"]]]
      nt1 <- setdiff(names(t1),c("full_chain","chains","total_chain_length","temp_start"))
      mt1 <- 0
      for(n in nt1) {
        el <- t1[[n]]
        mt1 <- max(mt1,length(el$sequence))
      }
      if(nrow(tsequences_product)==0 || (!summary[i,"product"] %in% tsequences_product[,"product"])) {
        tsequences_product <- rbind(tsequences_product,data.frame(product=summary[i,"product"],max_seq=mt1))
      } else {
        mt2 <- tsequences_product[tsequences_product[,"product"]==summary[i,"product"],"max_seq",drop=TRUE]
        tsequences_product[tsequences_product[,"product"]==summary[i,"product"],"max_seq"] <- max(mt1,mt2)
      }
    }
    for(pn in names(res$sequences)) {
      seq <- res$sequences[[pn]]
      nz <- names(seq)
      nz <- setdiff(nz,c("full_chain","chains","total_chain_length","temp_start"))
      for(sn in nz)
        v3 <- max(v3,length(seq[[sn]]$sequence))
    }
    tsl <- v3
    tsummary <- rbind(tsummary,summary)
    tsequences <- rbind(tsequences,data.frame(threshold=th,total_sequences=v2))
    ttotal <- rbind(ttotal,data.frame(threshold=th,total_pushes=v1))
    tsequences_tmp <- rbind(tsequences_tmp,data.frame(threshold=th,value=v3))
    cat(paste0("Total transitions:",v3,"\n"))
    th <- th + th_increment
  }
  return(list(tsummary=tsummary,tsequences=tsequences,ttotal=ttotal,tmaxproductsequences=tsequences_product,tindsequences=res$sequences))
}
