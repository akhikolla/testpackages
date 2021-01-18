### R code from vignette source 'grim.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: grim.Rnw:25-28
###################################################
require( gRbase )
prettyVersion <- packageDescription("gRim")$Version
prettyDate <- format(Sys.Date())


###################################################
### code chunk number 2: grim.Rnw:78-82
###################################################
dir.create("figures")
oopt <- options()
options("digits"=4, "width"=80)
options(useFancyQuotes="UTF-8")


###################################################
### code chunk number 3: grim.Rnw:89-93
###################################################
options("width"=85)
library(gRim)
library(Rgraphviz)
ps.options(family="serif")


###################################################
### code chunk number 4: grim.Rnw:125-128
###################################################
args(dmod)
args(cmod)
args(mmod)


###################################################
### code chunk number 5: grim.Rnw:146-148
###################################################
data(reinis)
str(reinis)


###################################################
### code chunk number 6: grim.Rnw:160-164
###################################################
data(reinis)
dm1<-dmod(list(c("smoke","systol"),c("smoke","mental","phys")), data=reinis)
dm1<-dmod(~smoke:systol + smoke:mental:phys, data=reinis)
dm1


###################################################
### code chunk number 7: grim.Rnw:184-186
###################################################
formula(dm1)
str(terms(dm1))


###################################################
### code chunk number 8: grim.Rnw:192-193
###################################################
summary(dm1)


###################################################
### code chunk number 9: grim.Rnw:224-227
###################################################
dm2 <- dmod(~.^2, margin=c("smo","men","phy","sys"),
            data=reinis)
formula(dm2)


###################################################
### code chunk number 10: grim.Rnw:231-234
###################################################
dm3 <- dmod(list(c("smoke", "systol"), c("smoke", "mental", "phys")),
            data=reinis, interactions=2)
formula(dm3)


###################################################
### code chunk number 11: grim.Rnw:248-249
###################################################
iplot(dm1)


###################################################
### code chunk number 12: grim.Rnw:263-267
###################################################
data(carcass)
cm1 <- cmod(~Fat11:Fat12:Fat13, data=carcass)
cm1 <- cmod(~Fat11:Fat12 + Fat12:Fat13 + Fat11:Fat13, data=carcass)
cm1


###################################################
### code chunk number 13: grim.Rnw:273-274
###################################################
iplot(cm1)


###################################################
### code chunk number 14: grim.Rnw:281-284
###################################################
data(milkcomp1)
mm1 <- mmod(~.^., data=milkcomp1)
mm1


###################################################
### code chunk number 15: grim.Rnw:290-291
###################################################
iplot(mm1)


###################################################
### code chunk number 16: grim.Rnw:308-310
###################################################
ms <- dmod(~.^., marginal=c("phys","mental","systol","family"), data=reinis)
formula(ms)


###################################################
### code chunk number 17: grim.Rnw:316-318
###################################################
ms1 <- update(ms, list(dedge=~phys:mental))
formula(ms1)


###################################################
### code chunk number 18: grim.Rnw:324-326
###################################################
ms2<- update(ms, list(dedge=~phys:mental+systol:family))
formula(ms2)


###################################################
### code chunk number 19: grim.Rnw:332-334
###################################################
ms3 <- update(ms, list(dedge=~phys:mental:systol))
formula(ms3)


###################################################
### code chunk number 20: grim.Rnw:340-342
###################################################
ms4 <- update(ms, list(dterm=~phys:mental:systol) )
formula(ms4)


###################################################
### code chunk number 21: grim.Rnw:348-350
###################################################
ms5 <- update(ms, list(aterm=~phys:mental+phys:systol+mental:systol) )
formula(ms5)


###################################################
### code chunk number 22: grim.Rnw:356-358
###################################################
ms6 <- update(ms, list(aedge=~phys:mental+systol:family))
formula(ms6)


###################################################
### code chunk number 23: grim.Rnw:379-380
###################################################
cit <- ciTest(reinis, set=c("systol","smoke","family","phys"))


###################################################
### code chunk number 24: grim.Rnw:413-414
###################################################
cit$slice


###################################################
### code chunk number 25: grim.Rnw:446-447
###################################################
ciTest(reinis, set=c("systol","smoke","family","phys"), method='MC')


###################################################
### code chunk number 26: grim.Rnw:458-460
###################################################
dm5 <- dmod(~ment:phys:systol+ment:systol:family+phys:systol:smoke,
            data=reinis)


###################################################
### code chunk number 27: fundamentalfig1
###################################################
iplot(dm5)


###################################################
### code chunk number 28: grim.Rnw:485-487
###################################################
testdelete(dm5, ~smoke:systol)
testdelete(dm5, ~family:systol)


###################################################
### code chunk number 29: grim.Rnw:506-507
###################################################
testadd(dm5, ~smoke:mental)


###################################################
### code chunk number 30: grim.Rnw:533-534
###################################################
ed.in <- getInEdges(ugList(terms(dm5)), type="decomposable")


###################################################
### code chunk number 31: grim.Rnw:547-548
###################################################
ed.out <- getOutEdges(ugList(terms(dm5)), type="decomposable")


###################################################
### code chunk number 32: grim.Rnw:556-558
###################################################
args(testInEdges)
args(testOutEdges)


###################################################
### code chunk number 33: grim.Rnw:572-574
###################################################
testInEdges(dm5, getInEdges(ugList(terms(dm5)), type="decomposable"),
             k=log(sum(reinis)))


###################################################
### code chunk number 34: grim.Rnw:596-599
###################################################
dm.sat <- dmod(~.^., data=reinis)
dm.back <- backward(dm.sat)
iplot(dm.back)


###################################################
### code chunk number 35: grim.Rnw:614-617
###################################################
dm.i   <- dmod(~.^1, data=reinis)
dm.forw <- forward(dm.i)
iplot(dm.forw)


###################################################
### code chunk number 36: grim.Rnw:681-685
###################################################
fix <- list(c("smoke","phys","systol"), c("systol","protein"))
fix <- do.call(rbind, unlist(lapply(fix, names2pairs),recursive=FALSE))
fix
dm.s3 <- backward(dm.sat, fixin=fix, details=1)


###################################################
### code chunk number 37: grim.Rnw:696-697
###################################################
dm.i3 <- forward(dm.i, fixout=fix, details=1)


###################################################
### code chunk number 38: grim.Rnw:903-904
###################################################
dim_loglin(terms(dm2), reinis)


###################################################
### code chunk number 39: grim.Rnw:958-961
###################################################
dm3 <- dmod(list(c("smoke", "systol"), c("smoke", "mental", "phys")),
            data=reinis)
names(dm3)


###################################################
### code chunk number 40: grim.Rnw:969-970
###################################################
str(terms(dm3))


###################################################
### code chunk number 41: grim.Rnw:974-975
###################################################
str(dm3$glistNUM)


###################################################
### code chunk number 42: grim.Rnw:981-982
###################################################
dm3$varNames


###################################################
### code chunk number 43: grim.Rnw:990-991
###################################################
str(dm3[c("varNames","conNames","conLevels")])


