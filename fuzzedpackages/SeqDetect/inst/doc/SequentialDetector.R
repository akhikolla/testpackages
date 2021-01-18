### R code from vignette source 'SequentialDetector.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: setup
###################################################
library(SeqDetect)
library(xtable)
library(dplyr)


###################################################
### code chunk number 2: SequentialDetector.Rnw:122-129
###################################################
st1 <- data.frame(patient=c("C156","C156","E9383","C167"),
                  time=c("05.12.2019. 10:30:20","05.12.2019. 11:59:07","07.12.2019. 08:34:12",
                         "07.12.2019. 10:45:11"),
                  age=c(12,12,26,76),
                  fever=c(TRUE,TRUE,TRUE,FALSE),
                  action=c("registration","release","registration","registration"),
                  can_walk=c(TRUE,TRUE,FALSE,TRUE))


###################################################
### code chunk number 3: SequentialDetector.Rnw:132-133
###################################################
  print(xtable(st1,caption = "ER registration system data stream slice"))


###################################################
### code chunk number 4: SequentialDetector.Rnw:137-138
###################################################
st1 <- transform.data.frame(st1,time=as.POSIXct(st1$time,format="%d.%m.%Y. %H:%M:%S"))


###################################################
### code chunk number 5: SequentialDetector.Rnw:144-152
###################################################
st2 <- data.frame(patient=c("C156","C156","E9383","E9383","C167","C167"),
                  time=c("05.12.2019. 10:41:00","05.12.2019. 12:12:00","07.12.2019. 09:56:00",
                         "07.12.2019. 11:32:00","07.12.2019. 11:01:00","07.12.2019. 13:14:15"),
                  diagnosis=c(NA,"J04.0",NA,"A41.9",NA,"N41.0"),
                  action=c("biomarker","release","biomarker","hospital_ic","biomarker","hospital_nc"),
                  decription=c("suspect. laryng...","course of antibiotics",
                               "high fever,in shock state!! URGENT!","septic shock? IC..",
                               "cannot pee,catheter","urology hospitalization"))


###################################################
### code chunk number 6: SequentialDetector.Rnw:155-156
###################################################
  print(xtable(st2,caption = "ER triage system data stream slice"),size="\\footnotesize")


###################################################
### code chunk number 7: SequentialDetector.Rnw:160-161
###################################################
st2 <- transform.data.frame(st2,time=as.POSIXct(st2$time,format="%d.%m.%Y. %H:%M:%S"))


###################################################
### code chunk number 8: SequentialDetector.Rnw:167-175
###################################################
st3 <- data.frame(request_id=c("2019_645553","2019_654331","2019_654331","2019_654331","2019_655376",
                               "2019_655376"),
                  request_org=c("ER","ER","ER","ER","ER","ER"),
                  ext_id=c("C156","E9383","E9383","E9383","C167","C167"),
                  date=c("05.12.2019.","07.12.2019.","07.12.2019.","07.12.2019.","07.12.2019.",
                         "07.12.2019."),
                  biomarker=c("WBC","WBC","CRP","LAC","WBC","CRP"),
                  final=c(14.6,13.11,345.0,4.5,11.43,67.0),stringsAsFactors=FALSE)


###################################################
### code chunk number 9: SequentialDetector.Rnw:178-179
###################################################
  print(xtable(st3,caption = "Biolab system data stream slice"))


###################################################
### code chunk number 10: SequentialDetector.Rnw:183-184
###################################################
st3 <- transform.data.frame(st3,date=as.POSIXct(st3$date,format="%d.%m.%Y."))


###################################################
### code chunk number 11: SequentialDetector.Rnw:190-196
###################################################
st4 <- data.frame(patient_id=c("I93382","N94511"),
                  ext_id=c("E9383","C167"),
                  time_in=c("07.12.2019. 11:35:46","07.12.2019. 12:11:49"),
                  diagnosis=c("A41.9","N41.0"),
                  type=c("IC","NC"),
                  time_release=c("15.12.2019. 08:52:11","11.12.2019. 14:02:11"))


###################################################
### code chunk number 12: SequentialDetector.Rnw:199-200
###################################################
  print(xtable(st4,caption = "Hospital system data stream slice"))


###################################################
### code chunk number 13: SequentialDetector.Rnw:204-206
###################################################
st4 <- transform.data.frame(st4,time_in=as.POSIXct(st4$time_in,format="%d.%m.%Y. %H:%M:%S"),
                            time_release=as.POSIXct(st4$time_release,format="%d.%m.%Y. %H:%M:%S"))


###################################################
### code chunk number 14: SequentialDetector.Rnw:213-266
###################################################
HSC_PP_Hospital <- function(...) {
  structure(list(),class = c("HSC_PP_Hospital","HSC_PP"))
}
preprocess.HSC_PP_Hospital <- function(x, streams, ...) {
  # perform some meaningful checking on the input data streams
  res <- data.frame(stringsAsFactors=FALSE)
  reg_stream <- streams[["registration_system"]]
  for(j in 1:nrow(reg_stream)) {
    el <- reg_stream[j,]
    cz <- case_when(el[,"action"]=="registration"~"ER registration",
                    el[,"action"]=="release"~"ER release")
    res <- rbind(res,data.frame(id=el[,"patient"],class=cz,time=el[,"time"],out=cz,WBC=NA,CRP=NA,LAC=NA))
  }
  triage_stream <- streams[["triage"]]
  for(j in 1:nrow(triage_stream)) {
    el <- triage_stream[j,]
    if(nrow(res[res$id==el[,"patient"] & res$class=="ER triage",])==0)
      res <- rbind(res,data.frame(id=el[,"patient"],class="ER triage",
                                  time=min(triage_stream[triage_stream[,"patient"]==el[,"patient"],"time"]),
                                  out="ER triage",WBC=NA,CRP=NA,LAC=NA))
  }
  biolab_stream <- streams[["biolab"]]
  for(j in 1:nrow(biolab_stream)) {
    el <- biolab_stream[j,]
    if(nrow(res[res$id==el[,"ext_id"] & res[,"class"]=="Biomarker assessment",])==0) {
      if(el[,"request_org"]=="ER") {
        t1_time <- triage_stream[triage_stream[,"patient"]==el[,"ext_id"] & 
                                   triage_stream[,"action"]=="biomarker","time"]+1
        res <- rbind(res,data.frame(id=el[,"ext_id"],class="Biomarker assessment",
                                    time=t1_time,out="Biomarker assessment",
                                    WBC=NA,CRP=NA,LAC=NA))
      }
    }
    res[res[,"id"]==el[,"ext_id"] & res[,"class"]=="Biomarker assessment",
        el[,"biomarker"]] <- el[,"final"]
  }
  hospital_stream <- streams[["hospital"]]
  for(j in 1:nrow(hospital_stream)) {
    el <- hospital_stream[j,]
    cz1 <- paste0("Admission to ",el[,"type"])
    res <- rbind(res,data.frame(id=el[,"ext_id"],class=cz1,time=el[,"time_in"],
                                out=cz1,WBC=NA,CRP=NA,LAC=NA))
    res <- rbind(res,data.frame(id=el[,"ext_id"],class="Hospital release",
                                time=el[,"time_release"],out="Hospital release",
                                WBC=NA,CRP=NA,LAC=NA))
  }
  res <- res[order(res[,"time"]),] # order this event stream slice
  return(list(obj=x,res=res))
}

pp <- HSC_PP_Hospital()
input_streams <- list(registration_system=st1,triage=st2,biolab=st3,hospital=st4)
event_stream <- preprocess(pp,input_streams)$res


###################################################
### code chunk number 15: SequentialDetector.Rnw:270-271
###################################################
  print(xtable(event_stream,caption="The resulting event stream slice"),size="\\footnotesize")


###################################################
### code chunk number 16: SequentialDetector.Rnw:306-339
###################################################
HSC_PC_Hospital <- function(...) {
  structure(list(),class = c("HSC_PC_Hospital","HSC_PC"))
}
classify.HSC_PC_Hospital <- function(x, event_stream, ...) {
  # perform some meaningful checking on the supplied event stream
  res <- data.frame(stringsAsFactors=FALSE)
  for(i in 1:nrow(event_stream)) {
    event <- event_stream[i,]
    symbol <- case_when(
      event$class=="ER registration" ~ "ER_REG",
      event$class=="ER triage" ~ "ER_TR",
      event$class=="ER release" ~ "ER_REL",
      event$class=="Biomarker assessment" ~ "BIO_A",
      event$class=="Admission to IC" ~ "IC",
      event$class=="Admission to NC" ~ "NC",
      event$class=="Hospital release" ~ "H_REL"
    )
    if(symbol=="BIO_A") {
      symbol <- paste0(symbol,case_when(is.na(event$WBC) ~ "#WBC=NONE",
                                        event$WBC>11 ~ "#WBC=EL",TRUE ~ "#WBC=OK"))
      symbol <- paste0(symbol,case_when(is.na(event$CRP) ~ "#CRP=NONE",
                                        event$CRP>50 ~ "#CRP=EL",TRUE ~ "#CRP=OK"))
      symbol <- paste0(symbol,case_when(is.na(event$LAC) ~ "#LAC=NONE",
                                        event$LAC>4 ~ "#LAC=EL",TRUE ~ "#LAC=OK"))
    }
    res <- rbind(res,data.frame(id=event$id,class=event$class,time=event$time,out=event$out,
                                .clazz=symbol))
  }
  return(res)
}

pc <- HSC_PC_Hospital()
consolidated_stream <- classify(pc,event_stream)


###################################################
### code chunk number 17: SequentialDetector.Rnw:342-343
###################################################
  print(xtable(consolidated_stream,caption="The consolidated data stream slice"),size="\\tiny")


###################################################
### code chunk number 18: SequentialDetector.Rnw:349-353
###################################################
seq_detector <- HybridSequenceClassifier(c("id","class","time","out"),"time","time",
                                         "id",preclassifier=pc,preprocessor=pp,
                                         pattern_field="out")
seq_detector$process(input_streams)


###################################################
### code chunk number 19: SequentialDetector.Rnw:357-358 (eval = FALSE)
###################################################
## seq_detector$printMachines()


###################################################
### code chunk number 20: SequentialDetector.Rnw:361-362
###################################################
seq_detector$printMachines()


###################################################
### code chunk number 21: SequentialDetector.Rnw:367-368
###################################################
seq_detector$plotMachines()


###################################################
### code chunk number 22: SequentialDetector.Rnw:376-377
###################################################
st <- data.frame(sequence=c("A","E","G"),alert=c(NA,NA,"Alert 1"))


###################################################
### code chunk number 23: SequentialDetector.Rnw:380-381 (eval = FALSE)
###################################################
## st


###################################################
### code chunk number 24: SequentialDetector.Rnw:383-384
###################################################
  print(xtable(st,caption="Learning sequence \\ref{seq:ps_learn1}"))


###################################################
### code chunk number 25: SequentialDetector.Rnw:387-396
###################################################
pp <- HSC_PP(c("sequence","alert"),"sequence_id",create_unique_key=TRUE,auto_id=TRUE)
pc <- HSC_PC_Attribute("sequence")
seq_detector <- HybridSequenceClassifier(c("sequence_id","sequence","alert"),
                                         "sequence_id","sequence_id",
                                         preclassifier=pc,preprocessor=pp,
                                         pattern_field="alert",reuse_states=FALSE)
input_streams <- list(stream=st)
seq_detector$process(input_streams,learn=TRUE)
seq_detector$cleanKeys()


###################################################
### code chunk number 26: SequentialDetector.Rnw:400-401 (eval = FALSE)
###################################################
## seq_detector$printMachines(print_cache=FALSE,print_keys=FALSE)


###################################################
### code chunk number 27: SequentialDetector.Rnw:404-405
###################################################
seq_detector$printMachines(print_cache=FALSE,print_keys=FALSE)


###################################################
### code chunk number 28: SequentialDetector.Rnw:411-418
###################################################
seq_test1 <- c("E","I","A","G","E","K","F","E","A","G","G","B","W","L")
res_test1 <- seq_detector$process(list(stream=data.frame(sequence=seq_test1,
                                                         alert=NA)),learn=FALSE)
out_test1 <- data.frame()
for(i in 1:nrow(res_test1$stream)) 
  out_test1 <- rbind(out_test1,data.frame(sequence=res_test1$stream[i,"sequence"],
                            alert=c_to_string(res_test1$explanation[[i]]$actual)))


###################################################
### code chunk number 29: SequentialDetector.Rnw:421-422 (eval = FALSE)
###################################################
## out_test1


###################################################
### code chunk number 30: SequentialDetector.Rnw:424-425
###################################################
  print(xtable(out_test1,caption="Results for testing sequence \\ref{seq:ps_test1}, $\\lambda_n=\\infty$",label="tab:ps_test1_res"))


###################################################
### code chunk number 31: SequentialDetector.Rnw:429-439
###################################################
pp <- HSC_PP(c("sequence","alert"),"sequence_id",create_unique_key=TRUE,auto_id=TRUE)
pc <- HSC_PC_Attribute("sequence")
dd <- list(type="count",count=0,context_related=TRUE)
seq_detector <- HybridSequenceClassifier(c("sequence_id","sequence","alert"),
                                         "sequence_id","sequence_id",
                                         preclassifier=pc,preprocessor=pp,
                                         pattern_field="alert",reuse_states=FALSE,
                                         decay_descriptors=list(d1=dd))
seq_detector$process(input_streams,learn=TRUE)
seq_detector$cleanKeys()


###################################################
### code chunk number 32: SequentialDetector.Rnw:442-448
###################################################
res_test1 <- seq_detector$process(list(stream=data.frame(sequence=seq_test1,alert=NA)),
                                  learn=FALSE)
out_test1 <- data.frame()
for(i in 1:nrow(res_test1$stream)) 
  out_test1 <- rbind(out_test1,data.frame(sequence=res_test1$stream[i,"sequence"],
                                alert=c_to_string(res_test1$explanation[[i]]$actual)))


###################################################
### code chunk number 33: SequentialDetector.Rnw:450-451 (eval = FALSE)
###################################################
## out_test1


###################################################
### code chunk number 34: SequentialDetector.Rnw:453-454
###################################################
  print(xtable(out_test1,caption="Result for testing sequence \\ref{seq:ps_test1}, $\\lambda_n=0$",label="tab:ps_test1_res2"))


###################################################
### code chunk number 35: SequentialDetector.Rnw:459-466
###################################################
seq_test2 <- c("E","I","A","E","G","E","K","F","E","A","E","B","G","W","L")
res_test2 <- seq_detector$process(list(stream=data.frame(sequence=seq_test2,alert=NA)),
                                  learn=FALSE)
out_test2 <- data.frame()
for(i in 1:nrow(res_test2$stream)) 
  out_test2 <- rbind(out_test2,data.frame(sequence=res_test2$stream[i,"sequence"],
                                alert=c_to_string(res_test2$explanation[[i]]$actual)))


###################################################
### code chunk number 36: SequentialDetector.Rnw:469-470 (eval = FALSE)
###################################################
## out_test2


###################################################
### code chunk number 37: SequentialDetector.Rnw:472-473
###################################################
  print(xtable(out_test2,caption="Result for testing sequence \\ref{seq:ps_test2}, $\\lambda_n=0$",label="tab:ps_test2_res"))


###################################################
### code chunk number 38: SequentialDetector.Rnw:477-491
###################################################
dd <- list(type="count",count=1,context_related=TRUE)
seq_detector <- HybridSequenceClassifier(c("sequence_id","sequence","alert"),
                                         "sequence_id","sequence_id",
                                         preclassifier=pc,preprocessor=pp,
                                         pattern_field="alert",reuse_states=FALSE,
                                         decay_descriptors=list(d1=dd))
seq_detector$process(input_streams,learn=TRUE) # learn
seq_detector$cleanKeys() # clean all context keys
res_test2 <- seq_detector$process(list(stream=data.frame(sequence=seq_test2,alert=NA)),
                                  learn=FALSE)
out_test2 <- data.frame()
for(i in 1:nrow(res_test2$stream)) 
  out_test2 <- rbind(out_test2,data.frame(sequence=res_test2$stream[i,"sequence"],
                                alert=c_to_string(res_test2$explanation[[i]]$actual)))


###################################################
### code chunk number 39: SequentialDetector.Rnw:493-494 (eval = FALSE)
###################################################
## out_test2


###################################################
### code chunk number 40: SequentialDetector.Rnw:496-497
###################################################
  print(xtable(out_test2,caption="Result for testing sequence \\ref{seq:ps_test2}, $\\lambda_n=1$",label="tab:ps_test2_res2"))


###################################################
### code chunk number 41: SequentialDetector.Rnw:506-510
###################################################
st <- data.frame(product=c("P45","P134","P45","P134","P134","P45","P134"),
                 sales=c(2,12,18,16,18,24,8),
                 alert=c(NA,NA,NA,NA,NA,"Alert P45","Alert P134"))
input_streams <- list(stream=st)


###################################################
### code chunk number 42: SequentialDetector.Rnw:513-514 (eval = FALSE)
###################################################
## st


###################################################
### code chunk number 43: SequentialDetector.Rnw:516-517
###################################################
  print(xtable(st,caption="A multi-contextual learning sequence \\ref{seq:mc_learn1}",label="tab:mc_learn1"))


###################################################
### code chunk number 44: SequentialDetector.Rnw:522-531
###################################################
pp <- HSC_PP(c("product","sales","alert"),"sequence_id",auto_id=TRUE)
pc <- HSC_PC_Attribute("sales")
dd <- list(type="count",count=0,context_related=TRUE) 
seq_detector <- HybridSequenceClassifier(c("sequence_id","product","sales","alert"),"sequence_id",
                                         "sequence_id",context_field="product",preclassifier=pc,
                                         preprocessor=pp,pattern_field="alert",reuse_states=FALSE,
                                         decay_descriptors=list(d1=dd))
seq_detector$process(input_streams,learn=TRUE)
seq_detector$cleanKeys()


###################################################
### code chunk number 45: SequentialDetector.Rnw:534-535 (eval = FALSE)
###################################################
## seq_detector$printMachines(print_cache=FALSE,print_keys=FALSE)


###################################################
### code chunk number 46: SequentialDetector.Rnw:538-539
###################################################
seq_detector$printMachines(print_cache=FALSE,print_keys=FALSE)


###################################################
### code chunk number 47: SequentialDetector.Rnw:553-554
###################################################
seq_detector$mergeMachines()


###################################################
### code chunk number 48: SequentialDetector.Rnw:556-557 (eval = FALSE)
###################################################
## seq_detector$printMachines(print_cache=FALSE,print_keys=FALSE)


###################################################
### code chunk number 49: SequentialDetector.Rnw:560-561
###################################################
seq_detector$printMachines(print_cache=FALSE,print_keys=FALSE)


###################################################
### code chunk number 50: SequentialDetector.Rnw:570-574
###################################################
tt <- data.frame(product=c("P672","P113","P983","P23872","P5","P672","P2982","P983","P672",
                           "P991","P983","P113","P2982","P344"),
                 sales=c(2,11,12,98,8,18,298,16,24,25,18,16,43,101),alert=NA)
test_streams <- list(stream=tt)


###################################################
### code chunk number 51: SequentialDetector.Rnw:577-578 (eval = FALSE)
###################################################
## tt


###################################################
### code chunk number 52: SequentialDetector.Rnw:580-581
###################################################
  print(xtable(tt,caption="A multi-contextual testing sequence \\ref{seq:mc_test1}",label="tab:mc_test1"))


###################################################
### code chunk number 53: SequentialDetector.Rnw:584-590
###################################################
res_test1 <- seq_detector$process(test_streams,learn=FALSE)
out_test1 <- data.frame()
for(i in 1:nrow(res_test1$stream)) 
  out_test1 <- rbind(out_test1,data.frame(product=res_test1$stream[i,"product"],
                                          sales=res_test1$stream[i,"sales"],
                                          alert=c_to_string(res_test1$explanation[[i]]$actual)))


###################################################
### code chunk number 54: SequentialDetector.Rnw:593-594 (eval = FALSE)
###################################################
## out_test1 


###################################################
### code chunk number 55: SequentialDetector.Rnw:597-598
###################################################
  print(xtable(out_test1,caption="Results for the testing sequence \\ref{seq:mc_test1}",label="tab:mc_test1_res"))


###################################################
### code chunk number 56: SequentialDetector.Rnw:602-603 (eval = FALSE)
###################################################
## seq_detector$printMachines()


###################################################
### code chunk number 57: SequentialDetector.Rnw:606-607
###################################################
seq_detector$printMachines()


###################################################
### code chunk number 58: SequentialDetector.Rnw:612-615
###################################################
tt <- data.frame(product=c("P115","P45","P22","P983","P9","P19","P73"),
                 sales=c(91,43,52,8,1,105,35),alert=NA)
test_streams <- list(stream=tt)


###################################################
### code chunk number 59: SequentialDetector.Rnw:617-618 (eval = FALSE)
###################################################
## tt


###################################################
### code chunk number 60: SequentialDetector.Rnw:620-621
###################################################
  print(xtable(tt,caption="Additional slice of the testing sequence \\ref{seq:mc_test1}",label="tab:mc_test1_as"))


###################################################
### code chunk number 61: SequentialDetector.Rnw:625-631
###################################################
res_test1 <- seq_detector$process(test_streams,learn=FALSE)
out_test1 <- data.frame()
for(i in 1:nrow(res_test1$stream)) 
  out_test1 <- rbind(out_test1,data.frame(product=res_test1$stream[i,"product"],
                                          sales=res_test1$stream[i,"sales"],
                                          alert=c_to_string(res_test1$explanation[[i]]$actual)))


###################################################
### code chunk number 62: SequentialDetector.Rnw:634-635 (eval = FALSE)
###################################################
## out_test1


###################################################
### code chunk number 63: SequentialDetector.Rnw:639-640
###################################################
  print(xtable(out_test1,caption="Results for the additional slice of the testing sequence \\ref{seq:mc_test1}",label="tab:mc_test1_as_res"))


###################################################
### code chunk number 64: SequentialDetector.Rnw:662-663 (eval = FALSE)
###################################################
## dd1 <- list(type="count",count=1,context_related=TRUE)


###################################################
### code chunk number 65: SequentialDetector.Rnw:668-669 (eval = FALSE)
###################################################
## dd2 <- list(type="count",count=50,context_related=FALSE)


###################################################
### code chunk number 66: SequentialDetector.Rnw:674-675 (eval = FALSE)
###################################################
## dd3 <- list(type="time",days=0,hours=1,minutes=10,context_related=FALSE)


###################################################
### code chunk number 67: SequentialDetector.Rnw:683-685
###################################################
st <- data.frame(product=c("P45","P45"),sales=c(5,10),alert=c(NA,"Alert 1"))
input_streams <- list(stream=st)


###################################################
### code chunk number 68: SequentialDetector.Rnw:688-689 (eval = FALSE)
###################################################
## st


###################################################
### code chunk number 69: SequentialDetector.Rnw:691-692
###################################################
  print(xtable(st,caption="Learning sequence \\ref{seq:td_learn1}",label="tab:td_learn1"))


###################################################
### code chunk number 70: SequentialDetector.Rnw:697-706
###################################################
pp <- HSC_PP(c("product","sales","alert"),"sequence_id",auto_id=TRUE)
pc <- HSC_PC_Attribute("sales")
dd <- list(type="count",count=1,context_related=TRUE)
seq_detector <- HybridSequenceClassifier(c("sequence_id","product","sales","alert"),"sequence_id",
                                         "sequence_id",context_field="product",preclassifier=pc,
                                         preprocessor=pp,pattern_field="alert",reuse_states=FALSE,
                                         decay_descriptors=list(d1=dd))
seq_detector$process(input_streams,learn=TRUE)
seq_detector$cleanKeys()


###################################################
### code chunk number 71: SequentialDetector.Rnw:709-710 (eval = FALSE)
###################################################
## seq_detector$printMachines(print_cache=FALSE,print_keys=FALSE)


###################################################
### code chunk number 72: SequentialDetector.Rnw:713-714
###################################################
seq_detector$printMachines(print_cache=FALSE,print_keys=FALSE)


###################################################
### code chunk number 73: SequentialDetector.Rnw:720-723
###################################################
tt <- data.frame(product=c("P113","P29","P113","P29","P29","P113","P29"),
                 sales=c(5,5,7,8,9,10,10),alert=NA)
test_streams <- list(stream=tt)


###################################################
### code chunk number 74: SequentialDetector.Rnw:726-727 (eval = FALSE)
###################################################
## tt


###################################################
### code chunk number 75: SequentialDetector.Rnw:729-730
###################################################
  print(xtable(tt,caption="Testing sequence \\ref{seq:td_test1}",label="tab:td_test1"))


###################################################
### code chunk number 76: SequentialDetector.Rnw:735-741
###################################################
res_test1 <- seq_detector$process(test_streams,learn=FALSE)
out_test1 <- data.frame()
for(i in 1:nrow(res_test1$stream)) 
  out_test1 <- rbind(out_test1,data.frame(product=res_test1$stream[i,"product"],
                                          sales=res_test1$stream[i,"sales"],
                                          alert=c_to_string(res_test1$explanation[[i]]$actual)))


###################################################
### code chunk number 77: SequentialDetector.Rnw:744-745 (eval = FALSE)
###################################################
## out_test1


###################################################
### code chunk number 78: SequentialDetector.Rnw:747-748
###################################################
  print(xtable(out_test1,caption="Testing sequence \\ref{seq:td_test1} results",label="tab:td_test1_res"))


###################################################
### code chunk number 79: SequentialDetector.Rnw:755-762
###################################################
dd <- list(type="count",count=5,context_related=FALSE)
seq_detector <- HybridSequenceClassifier(c("sequence_id","product","sales","alert"),"sequence_id",
                                         "sequence_id",context_field="product",preclassifier=pc,
                                         preprocessor=pp,pattern_field="alert",reuse_states=FALSE,
                                         decay_descriptors=list(d1=dd))
seq_detector$process(input_streams,learn=TRUE)
seq_detector$cleanKeys()


###################################################
### code chunk number 80: SequentialDetector.Rnw:768-771
###################################################
tt <- data.frame(product=c("P29","P113","P114","P115","P113","P114","P115","P29"),
                 sales=c(5,5,5,5,10,7,10,10),alert=NA)
test_streams <- list(stream=tt)


###################################################
### code chunk number 81: SequentialDetector.Rnw:774-775 (eval = FALSE)
###################################################
## tt


###################################################
### code chunk number 82: SequentialDetector.Rnw:777-778
###################################################
  print(xtable(tt,caption="Testing sequence \\ref{seq:td_test2}",label="tab:td_test2"))


###################################################
### code chunk number 83: SequentialDetector.Rnw:783-789
###################################################
res_test2 <- seq_detector$process(test_streams,learn=FALSE)
out_test2 <- data.frame()
for(i in 1:nrow(res_test2$stream)) 
  out_test2 <- rbind(out_test2,data.frame(product=res_test2$stream[i,"product"],
                                          sales=res_test2$stream[i,"sales"],
                                          alert=c_to_string(res_test2$explanation[[i]]$actual)))


###################################################
### code chunk number 84: SequentialDetector.Rnw:792-793 (eval = FALSE)
###################################################
## out_test2


###################################################
### code chunk number 85: SequentialDetector.Rnw:795-796
###################################################
  print(xtable(out_test2,caption="Testing sequence \\ref{seq:td_test2} results",label="tab:td_test2_res"))


###################################################
### code chunk number 86: SequentialDetector.Rnw:803-805
###################################################
st <- data.frame(product=c("P21","P21"),timestamp=c("01.12.2019. 10:00:00","01.12.2019. 10:01:00"),
                 sales=c(5,10),alert=c(NA,"Alert 1"))


###################################################
### code chunk number 87: SequentialDetector.Rnw:808-809 (eval = FALSE)
###################################################
## st


###################################################
### code chunk number 88: SequentialDetector.Rnw:811-812
###################################################
  print(xtable(st,caption="Learning sequence \\ref{seq:td_learn1} with timestamp",label="tab:td_learn1_tsfield"))


###################################################
### code chunk number 89: SequentialDetector.Rnw:815-818
###################################################
st <- transform.data.frame(st,timestamp=as.POSIXct(st$timestamp,
                                                   format="%d.%m.%Y. %H:%M:%S"))
input_streams <- list(stream=st)


###################################################
### code chunk number 90: SequentialDetector.Rnw:822-831
###################################################
pp <- HSC_PP(c("product","sales","alert","timestamp"),"timestamp")
pc <- HSC_PC_Attribute("sales")
dd <- list(type="time",days=0,hours=1,minutes=0,context_related=FALSE)
seq_detector <- HybridSequenceClassifier(c("timestamp","product","sales","alert"),
                                         "timestamp","timestamp",context_field="product",
                                         preclassifier=pc,preprocessor=pp,pattern_field="alert",
                                         reuse_states=FALSE,decay_descriptors=list(d1=dd))
seq_detector$process(input_streams,learn=TRUE)
seq_detector$cleanKeys()


###################################################
### code chunk number 91: SequentialDetector.Rnw:836-842
###################################################
tt <- data.frame(product=c("P12","P13","P14","P15","P13","P14","P15","P12"),
                 sales=c(5,5,5,5,10,10,10,10),
                 timestamp=c("05.12.2019. 10:30:20","05.12.2019. 10:31:20",
                             "05.12.2019. 10:32:20","05.12.2019. 10:33:20",
                             "05.12.2019. 10:34:20","05.12.2019. 10:35:20",
                             "05.12.2019. 10:40:20","05.12.2019. 12:30:20"),alert=NA)


###################################################
### code chunk number 92: SequentialDetector.Rnw:845-846 (eval = FALSE)
###################################################
## tt


###################################################
### code chunk number 93: SequentialDetector.Rnw:848-849
###################################################
  print(xtable(tt,caption="Testing sequence \\ref{seq:td_test3}",label="tab:td_test3"))


###################################################
### code chunk number 94: SequentialDetector.Rnw:852-855
###################################################
tt <- transform.data.frame(tt,timestamp=as.POSIXct(tt$timestamp,
                                    format="%d.%m.%Y. %H:%M:%S"))
test_streams <- list(stream=tt)


###################################################
### code chunk number 95: SequentialDetector.Rnw:859-867
###################################################
res_test3 <- seq_detector$process(test_streams,learn=FALSE)
out_test3 <- data.frame()
for(i in 1:nrow(res_test3$stream)) 
  out_test3 <- rbind(out_test3,data.frame(product=res_test3$stream[i,"product"],
                                          sales=res_test3$stream[i,"sales"],
                                          timestamp=as.character(res_test3$stream[i,"timestamp"],
                                                                 format="%d.%m.%Y. %H:%M:%S"),
                                          alert=c_to_string(res_test3$explanation[[i]]$actual)))


###################################################
### code chunk number 96: SequentialDetector.Rnw:870-871 (eval = FALSE)
###################################################
## out_test3


###################################################
### code chunk number 97: SequentialDetector.Rnw:873-874
###################################################
  print(xtable(out_test3,caption="Testing sequence \\ref{seq:td_test3} results",label="tab:td_test3_res"))


###################################################
### code chunk number 98: SequentialDetector.Rnw:911-919
###################################################
st <- data.frame(product=c("P1","P2"),sales=c(5,76),alert=c(NA,NA))
for(i in 1:400) {
  st <- rbind(st,data.frame(product=c("P1","P2"),sales=c(10,58),alert=c(NA,NA)))
  st <- rbind(st,data.frame(product=c("P1","P2"),sales=c(20,31),alert=c(NA,NA)))
}
st <- rbind(st,data.frame(product=c("P1","P2"),sales=c(30,11),
                          alert=c("Sequence 1","Sequence 2")))
input_streams <- list(stream=st)


###################################################
### code chunk number 99: SequentialDetector.Rnw:925-932
###################################################
pp <- HSC_PP(c("product","sales","alert"),"sequence_id",auto_id=TRUE)
pc <- HSC_PC_Attribute("sales")
seq_detector <- HybridSequenceClassifier(c("sequence_id","product","sales","alert"),"sequence_id",
                                         "sequence_id",context_field="product",preclassifier=pc,
                                         preprocessor=pp,reuse_states=TRUE,pattern_field="alert")
seq_detector$process(input_streams,learn=TRUE)
seq_detector$cleanKeys()


###################################################
### code chunk number 100: SequentialDetector.Rnw:935-936 (eval = FALSE)
###################################################
## seq_detector$printMachines(print_cache=FALSE,print_keys=FALSE)


###################################################
### code chunk number 101: SequentialDetector.Rnw:939-940
###################################################
seq_detector$printMachines(print_cache=FALSE,print_keys=FALSE)


###################################################
### code chunk number 102: SequentialDetector.Rnw:945-949
###################################################
tt <- data.frame(product=c("P29","P29","P34","P29","P29","P11","P34","P34",
                           "P34","P11","P11"),
                 sales=c(5,10,76,20,30,10,58,31,11,20,30),alert=NA)
test_streams <- list(stream=tt)


###################################################
### code chunk number 103: SequentialDetector.Rnw:954-960
###################################################
res_test <- seq_detector$process(test_streams,learn=FALSE)
out_test <- data.frame()
for(i in 1:nrow(res_test$stream)) 
  out_test <- rbind(out_test,data.frame(product=res_test$stream[i,"product"],
                                        sales=res_test$stream[i,"sales"],
                                        alert=c_to_string(res_test$explanation[[i]]$actual)))


###################################################
### code chunk number 104: SequentialDetector.Rnw:963-964 (eval = FALSE)
###################################################
## out_test


###################################################
### code chunk number 105: SequentialDetector.Rnw:966-967
###################################################
  print(xtable(out_test,caption="Testing sequence \\ref{seq:proj_test1} results",label="tab:proj_test1_res"))


###################################################
### code chunk number 106: SequentialDetector.Rnw:971-972
###################################################
seq_detector1 <- seq_detector$clone()


###################################################
### code chunk number 107: SequentialDetector.Rnw:976-977
###################################################
res_is <- seq_detector$induceSubmachine(threshold=200)


###################################################
### code chunk number 108: SequentialDetector.Rnw:979-980 (eval = FALSE)
###################################################
## seq_detector$printMachines(print_cache=FALSE,print_keys=FALSE)


###################################################
### code chunk number 109: SequentialDetector.Rnw:983-984
###################################################
seq_detector$printMachines(print_cache=FALSE,print_keys=FALSE)


###################################################
### code chunk number 110: SequentialDetector.Rnw:989-997
###################################################
seq_detector$setOutputPattern(states=c("20"),transitions=c(),pattern="Reg.sequence 1")
seq_detector$setOutputPattern(states=c("31"),transitions=c(),pattern="Reg.sequence 2")
res_test <- seq_detector$process(test_streams,learn=FALSE)
out_test <- data.frame()
for(i in 1:nrow(res_test$stream)) 
  out_test <- rbind(out_test,data.frame(product=res_test$stream[i,"product"],
                                        sales=res_test$stream[i,"sales"],
                                        alert=c_to_string(res_test$explanation[[i]]$actual)))


###################################################
### code chunk number 111: SequentialDetector.Rnw:1000-1001 (eval = FALSE)
###################################################
## out_test


###################################################
### code chunk number 112: SequentialDetector.Rnw:1003-1004
###################################################
  print(xtable(out_test,caption="Testing sequence \\ref{seq:proj_test1} results",label="tab:proj_test1_res2"))


###################################################
### code chunk number 113: SequentialDetector.Rnw:1011-1013
###################################################
seq_detector <- seq_detector1$clone()
res_is <- seq_detector$induceSubmachine(threshold=200,isolate=TRUE)


###################################################
### code chunk number 114: SequentialDetector.Rnw:1015-1016 (eval = FALSE)
###################################################
## seq_detector$printMachines(print_cache=FALSE,print_keys=FALSE)


###################################################
### code chunk number 115: SequentialDetector.Rnw:1019-1020
###################################################
seq_detector$printMachines(print_cache=FALSE,print_keys=FALSE)


###################################################
### code chunk number 116: SequentialDetector.Rnw:1031-1039
###################################################
seq_detector$setOutputPattern(states=c("20"),transitions=c(),pattern="Reg.sequence 1")
seq_detector$setOutputPattern(states=c("31"),transitions=c(),pattern="Reg.sequence 2")
res_test <- seq_detector$process(test_streams,learn=FALSE)
out_test <- data.frame()
for(i in 1:nrow(res_test$stream)) 
  out_test <- rbind(out_test,data.frame(product=res_test$stream[i,"product"],
                                        sales=res_test$stream[i,"sales"],
                                        alert=c_to_string(res_test$explanation[[i]]$actual)))


###################################################
### code chunk number 117: SequentialDetector.Rnw:1042-1043 (eval = FALSE)
###################################################
## out_test


###################################################
### code chunk number 118: SequentialDetector.Rnw:1045-1046
###################################################
  print(xtable(out_test,caption="Testing sequence \\ref{seq:proj_test1} results",label="tab:proj_test1_res3"))


###################################################
### code chunk number 119: SequentialDetector.Rnw:1058-1064
###################################################
library(SeqDetect)
ldf1 <- data.frame(product=c("P1","P1","P1","P1"),sequence_id=c(1,3,5,7),
                   sales=c(5,76,123,1),alert=c(NA,NA,NA,"Alert P1"))
ldf2 <- data.frame(product=c("P2","P2","P2","P2"),sequence_id=c(2,4,6,8),
                   sales=c(21,76,123,42),alert=c(NA,NA,NA,"Alert P2"))
input_streams <- list(stream1=ldf1,stream2=ldf2)


###################################################
### code chunk number 120: SequentialDetector.Rnw:1069-1078
###################################################
pp <- HSC_PP(c("product","sales","alert","sequence_id"),"sequence_id")
pc <- HSC_PC_Attribute("sales")
seq_detector <- HybridSequenceClassifier(c("sequence_id","product","sales","alert"),
                                         "sequence_id","sequence_id",context_field="product",
                                         preclassifier=pc,preprocessor=pp,reuse_states=TRUE,
                                         pattern_field="alert")
seq_detector$process(input_streams,learn=TRUE)
seq_detector$cleanKeys()
backup_detector <- seq_detector$clone()


###################################################
### code chunk number 121: SequentialDetector.Rnw:1083-1084 (eval = FALSE)
###################################################
## seq_detector$printMachines(print_cache=FALSE,print_keys=FALSE)


###################################################
### code chunk number 122: SequentialDetector.Rnw:1087-1088
###################################################
seq_detector$printMachines(print_cache=FALSE,print_keys=FALSE)


###################################################
### code chunk number 123: SequentialDetector.Rnw:1094-1097
###################################################
tdf1 <- data.frame(product=c("P3","P3","P3","P3"),sequence_id=c(1,2,3,4),
                   sales=c(5,76,123,1),alert=NA)
test_streams <- list(stream1=tdf1)


###################################################
### code chunk number 124: SequentialDetector.Rnw:1101-1107
###################################################
res_test <- seq_detector$process(test_streams,learn=FALSE)
out_test <- data.frame()
for(i in 1:nrow(res_test$stream)) 
  out_test <- rbind(out_test,data.frame(product=res_test$stream[i,"product"],
                                        sales=res_test$stream[i,"sales"],
                                        alert=c_to_string(res_test$explanation[[i]]$actual)))


###################################################
### code chunk number 125: SequentialDetector.Rnw:1110-1111 (eval = FALSE)
###################################################
## out_test


###################################################
### code chunk number 126: SequentialDetector.Rnw:1113-1114
###################################################
  print(xtable(out_test,caption="Testing sequence \\ref{seq:comp_test1} results",label="tab:comp_test1_res1"))


###################################################
### code chunk number 127: SequentialDetector.Rnw:1119-1122
###################################################
tdf1 <- data.frame(product=c("P4","P4","P4","P4"),sequence_id=c(1,2,3,4),
                   sales=c(21,76,123,1),alert=NA)
test_streams <- list(stream1=tdf1)


###################################################
### code chunk number 128: SequentialDetector.Rnw:1126-1133
###################################################
seq_detector$mergeMachines()
res_test <- seq_detector$process(test_streams,learn=FALSE)
out_test <- data.frame()
for(i in 1:nrow(res_test$stream)) 
  out_test <- rbind(out_test,data.frame(product=res_test$stream[i,"product"],
                                        sales=res_test$stream[i,"sales"],
                                        alert=c_to_string(res_test$explanation[[i]]$actual)))


###################################################
### code chunk number 129: SequentialDetector.Rnw:1136-1137 (eval = FALSE)
###################################################
## seq_detector$printMachines(print_cache=FALSE,print_keys=FALSE)


###################################################
### code chunk number 130: SequentialDetector.Rnw:1140-1141
###################################################
seq_detector$printMachines(print_cache=FALSE,print_keys=FALSE)


###################################################
### code chunk number 131: SequentialDetector.Rnw:1144-1145 (eval = FALSE)
###################################################
## out_test


###################################################
### code chunk number 132: SequentialDetector.Rnw:1147-1148
###################################################
  print(xtable(out_test,caption="Testing sequence \\ref{seq:comp_test2} results",label="tab:comp_test2_res1"))


###################################################
### code chunk number 133: SequentialDetector.Rnw:1157-1160
###################################################
tdf1 <- data.frame(product=c("P5","P5","P5","P5"),sequence_id=c(1,2,3,4),
                   sales=c(21,76,123,1),alert=NA)
test_streams <- list(stream1=tdf1)


###################################################
### code chunk number 134: SequentialDetector.Rnw:1164-1172
###################################################
seq_detector <- backup_detector$clone()
seq_detector$compressMachines()
res_test <- seq_detector$process(test_streams,learn=FALSE)
out_test <- data.frame()
for(i in 1:nrow(res_test$stream)) 
  out_test <- rbind(out_test,data.frame(product=res_test$stream[i,"product"],
                                        sales=res_test$stream[i,"sales"],
                                        alert=c_to_string(res_test$explanation[[i]]$actual)))


###################################################
### code chunk number 135: SequentialDetector.Rnw:1176-1177 (eval = FALSE)
###################################################
## seq_detector$printMachines(print_cache=FALSE,print_keys=FALSE)


###################################################
### code chunk number 136: SequentialDetector.Rnw:1180-1181
###################################################
seq_detector$printMachines(print_cache=FALSE,print_keys=FALSE)


###################################################
### code chunk number 137: SequentialDetector.Rnw:1184-1185 (eval = FALSE)
###################################################
## out_test


###################################################
### code chunk number 138: SequentialDetector.Rnw:1187-1188
###################################################
  print(xtable(out_test,caption="Testing sequence \\ref{seq:comp_test3} results",label="tab:comp_test3_res1"))


###################################################
### code chunk number 139: SequentialDetector.Rnw:1193-1197
###################################################
tdf1 <- data.frame(product=c("P6","P6","P6","P6"),
                   sequence_id=c(1,2,3,4),sales=c(5,76,123,1),
                   alert=NA)
test_streams <- list(stream1=tdf1)


###################################################
### code chunk number 140: SequentialDetector.Rnw:1201-1207
###################################################
res_test <- seq_detector$process(test_streams,learn=FALSE)
out_test <- data.frame()
for(i in 1:nrow(res_test$stream)) 
  out_test <- rbind(out_test,data.frame(product=res_test$stream[i,"product"],
                                        sales=res_test$stream[i,"sales"],
                                        alert=c_to_string(res_test$explanation[[i]]$actual)))


###################################################
### code chunk number 141: SequentialDetector.Rnw:1210-1211 (eval = FALSE)
###################################################
## out_test


###################################################
### code chunk number 142: SequentialDetector.Rnw:1213-1214
###################################################
  print(xtable(out_test,caption="Testing sequence \\ref{seq:comp_test4} results",label="tab:comp_test4_res1"))


###################################################
### code chunk number 143: SequentialDetector.Rnw:1222-1224
###################################################
st <- data.frame(product=c("P1","P1"),sales=c(5,76),alert=c(NA,"Alert"))
input_streams <- list(stream=st)


###################################################
### code chunk number 144: SequentialDetector.Rnw:1228-1234
###################################################
pp <- HSC_PP(c("product","sales","alert"),"sequence_id",auto_id=TRUE)
pc <- HSC_PC_Attribute("sales")
seq_detector_oo <- HybridSequenceClassifier(c("sequence_id","product","sales","alert"),"sequence_id",
                                         "sequence_id",context_field="product",preclassifier=pc,
                                         preprocessor=pp,reuse_states=TRUE,pattern_field="alert")
seq_detector_oo$process(input_streams,learn=TRUE)


###################################################
### code chunk number 145: SequentialDetector.Rnw:1239-1240
###################################################
seq_detector_oo$serialize()


###################################################
### code chunk number 146: SequentialDetector.Rnw:1243-1247
###################################################
c_to_string(names(seq_detector_oo$cache))
c_to_string(names(seq_detector_oo$cache[[1]]))
c_to_string(names(seq_detector_oo$cache[[1]][["states"]]))
saveRDS(seq_detector_oo,"test.RDS")


###################################################
### code chunk number 147: SequentialDetector.Rnw:1251-1252
###################################################
new_seq_detector_oo <- readRDS("test.RDS")


###################################################
### code chunk number 148: SequentialDetector.Rnw:1254-1255
###################################################
file.remove("test.RDS")


###################################################
### code chunk number 149: SequentialDetector.Rnw:1257-1258 (eval = FALSE)
###################################################
## new_seq_detector_oo$printMachines()


###################################################
### code chunk number 150: SequentialDetector.Rnw:1261-1262
###################################################
new_seq_detector_oo$printMachines()


###################################################
### code chunk number 151: SequentialDetector.Rnw:1269-1274
###################################################
sd_list <- new_seq_detector_oo$serializeToList()
c_to_string(names(sd_list))
c_to_string(names(sd_list[[1]]))
c_to_string(names(sd_list[[1]][["states"]]))
totally_new_sd_oo <- deserializeFromList(sd_list)


###################################################
### code chunk number 152: SequentialDetector.Rnw:1276-1277 (eval = FALSE)
###################################################
## totally_new_sd_oo$printMachines()


###################################################
### code chunk number 153: SequentialDetector.Rnw:1280-1281
###################################################
totally_new_sd_oo$printMachines()


###################################################
### code chunk number 154: SequentialDetector.Rnw:1289-1293
###################################################
st <- data.frame(product=c("P5","P5","P5"),sales=c(9,76,10),
                 alert=c(NA,"Alert","Alert P5"))
input_streams <- list(stream=st)
totally_new_sd_oo$process(input_streams,learn=TRUE)


###################################################
### code chunk number 155: SequentialDetector.Rnw:1296-1297 (eval = FALSE)
###################################################
## totally_new_sd_oo$printMachines()


###################################################
### code chunk number 156: SequentialDetector.Rnw:1300-1301
###################################################
totally_new_sd_oo$printMachines()


###################################################
### code chunk number 157: SequentialDetector.Rnw:1306-1307
###################################################
totally_new_sd_oo$mergeMachines()


###################################################
### code chunk number 158: SequentialDetector.Rnw:1309-1310 (eval = FALSE)
###################################################
## totally_new_sd_oo$printMachines()


###################################################
### code chunk number 159: SequentialDetector.Rnw:1313-1314
###################################################
totally_new_sd_oo$printMachines()


###################################################
### code chunk number 160: SequentialDetector.Rnw:1320-1322
###################################################
tt <- data.frame(product=c("P10","P10","P10"),sales=c(5,76,10),alert=NA)
test_streams <- list(stream=tt)


###################################################
### code chunk number 161: SequentialDetector.Rnw:1326-1332
###################################################
res_test1 <- totally_new_sd_oo$process(test_streams,learn=FALSE)
out_test1 <- data.frame()
for(i in 1:nrow(res_test1$stream)) 
  out_test1 <- rbind(out_test1,data.frame(product=res_test1$stream[i,"product"],
                                          sales=res_test1$stream[i,"sales"],
                                          alert=c_to_string(res_test1$explanation[[i]]$actual)))


###################################################
### code chunk number 162: SequentialDetector.Rnw:1335-1336
###################################################
  print(xtable(out_test1,caption="Testing sequence \\ref{seq:merg_test1} results",label="tab:merg_test1_res1"))


