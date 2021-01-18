## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, fig.width=6, fig.height=5, dev='CairoPNG')

## ----get_data, results='hide', message=FALSE----------------------------------
library(dexter)
library(dplyr)

db = start_new_project(verbAggrRules, ":memory:", person_properties=list(gender=""))
add_booklet(db, verbAggrData, "data")
add_item_properties(db, verbAggrProperties)

## ----dif----------------------------------------------------------------------
dif_gender = DIF(db, "gender")
dif_gender

## ----abscale, echo=FALSE, results='hide', fig.height=3, fig.width=8, fig.align='center'----
cc=fit_enorm(db,(gender=="Male")&(item_id%in%c("S1DoCurse","S1WantCurse","S3WantScold")))
plot(c(-2,3),c(0,1),xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='ability-scale')
lines(-1.4:2.6,rep(0,5),lty=2,col="gray")
lines(-1.4:2.6,rep(0.8,5),lty=2,col="gray")
cf=coef(cc)
text(cf$beta,.8,paste(cf$item_id, cf$item_score),cex=0.6,adj=1,srt=90,pos = 3, xpd=NA)


cc=fit_enorm(db,(gender=="Female")&(item_id%in%c("S1DoCurse","S1WantCurse","S3WantScold")))
cf=coef(cc)
text(cf$beta,0,paste(cf$item_id, cf$item_score),cex=0.6,adj=1,srt=90,pos = 3, xpd=NA)
text(-1.3,0.8,"Males", pos=2, xpd=NA)
text(-1.3,0,"Females", pos=2, xpd=NA)

## ----plotdif,fig.align='center'-----------------------------------------------
plot(dif_gender)

## ----sorting, fig.align='center'----------------------------------------------
items = get_items(db) %>%
  arrange(mode, item_id)
  
plot(dif_gender, items=items$item_id)

## ---- echo=FALSE, results='hide'----------------------------------------------
close_project(db)

