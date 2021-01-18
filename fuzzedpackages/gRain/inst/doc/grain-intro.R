## ----include=FALSE,echo=FALSE,warning=FALSE----------------------------------------
library(knitr)
dir.create("figures")
opts_chunk$set(fig.height=2.5,
               fig.path='figures/grain-',
               warning=FALSE, message=FALSE
)
options("prompt"="> ","width"=85)

## ----echo=FALSE--------------------------------------------------------------------
require(gRain)
prettyVersion <- packageDescription("gRain")$Version
prettyDate <- format(Sys.Date())

## ----------------------------------------------------------------------------------
citation("gRain")

## ----echo=F, results='hide'--------------------------------------------------------
yn <- c("yes","no")
a    <- cptable(~asia, values=c(1,99),levels=yn)
t.a  <- cptable(~tub|asia, values=c(5,95,1,99),levels=yn)
s    <- cptable(~smoke, values=c(5,5), levels=yn)
l.s  <- cptable(~lung|smoke, values=c(1,9,1,99), levels=yn)
b.s  <- cptable(~bronc|smoke, values=c(6,4,3,7), levels=yn)
e.lt <- cptable(~either|lung:tub,values=c(1,0,1,0,1,0,0,1),levels=yn)
x.e  <- cptable(~xray|either, values=c(98,2,5,95), levels=yn)
d.be <- cptable(~dysp|bronc:either, values=c(9,1,7,3,8,2,1,9), levels=yn)
plist <- compileCPT(list(a, t.a, s, l.s, b.s, e.lt, x.e, d.be))
plist
chest_bn <- grain(plist)
chest_bn

## ----chest-LS, echo=F, fig.cap="Chest clinic example from Lauritzen and Spiegelhalter (1988)."----
plot(chest_bn)

## ----------------------------------------------------------------------------------
yn <- c("yes", "no")
a    <- cptable(~asia, values=c(1, 99), levels=yn)
t.a  <- cptable(~tub|asia, values=c(5, 95, 1, 99), levels=yn)
s    <- cptable(~smoke, values=c(5, 5), levels=yn)
l.s  <- cptable(~lung|smoke, values=c(1, 9, 1, 99), levels=yn)
b.s  <- cptable(~bronc|smoke, values=c(6, 4, 3, 7), levels=yn)
e.lt <- cptable(~either|lung:tub, values=c(1, 0, 1, 0, 1, 0, 0, 1), levels=yn)
x.e  <- cptable(~xray|either, values=c(98, 2, 5, 95), levels=yn)
d.be <- cptable(~dysp|bronc:either, values=c(9, 1, 7, 3, 8, 2, 1, 9), levels=yn)

## ----------------------------------------------------------------------------------
chest_cpt <- compileCPT(a, t.a, s, l.s, b.s, e.lt, x.e, d.be)
summary(chest_cpt)

## ----------------------------------------------------------------------------------
chest_cpt$tub
chest_cpt$tub  %>% as.data.frame.table

## ----------------------------------------------------------------------------------
chest_cpt$either  %>% as.data.frame.table

## ----------------------------------------------------------------------------------
chest_bn <- grain(chest_cpt)
chest_bn

## ----------------------------------------------------------------------------------
chest_bn <- compile(chest_bn)

## ----------------------------------------------------------------------------------
querygrain(chest_bn, nodes=c("lung", "bronc"), type="marginal")

## ----------------------------------------------------------------------------------
querygrain(chest_bn, nodes=c("lung", "bronc"), type="joint")

## ----------------------------------------------------------------------------------
chest_bn2  <- setEvidence(chest_bn, evidence=list(asia="yes", dysp="yes"))
chest_bn2  <- setEvidence(chest_bn,
                      nodes=c("asia", "dysp"), states=c("yes", "yes"))

## ----------------------------------------------------------------------------------
pEvidence(chest_bn2)

## ----------------------------------------------------------------------------------
querygrain(chest_bn2, nodes=c("lung", "bronc"))
querygrain(chest_bn2, nodes=c("lung", "bronc"), type="joint")

## ----------------------------------------------------------------------------------
querygrain(chest_bn, evidence=list(asia="yes", dysp="yes"),
           nodes=c("lung", "bronc"), type="joint")

## ----------------------------------------------------------------------------------
chest_bn3 <- setEvidence(chest_bn, evidence=list(either="no", tub="yes"))

## ----------------------------------------------------------------------------------
pEvidence(chest_bn3)
querygrain(chest_bn3, nodes=c("lung", "bronc"), type="joint")

## ----------------------------------------------------------------------------------
yn <- c("yes","no")
eps <- 1e-100
a    <- cptable(~a,   values=c(1, eps), levels=yn)
b.a  <- cptable(~b+a, values=c(1, eps, eps, 1), levels=yn)
c.b  <- cptable(~c+b, values=c(1, eps, eps, 1), levels=yn)
plist <- compileCPT(list(a, b.a, c.b))
bn   <- grain(plist)
tt   <- querygrain(bn, type="joint")
ftable(tt)
querygrain(setEvidence(bn, evidence=list(a="no", c="yes")))

## ----------------------------------------------------------------------------------
eps  <- 1e-200
a    <- cptable(~a,   values=c(1, eps),levels=yn)
b.a  <- cptable(~b+a, values=c(1, eps, eps, 1),levels=yn)
c.b  <- cptable(~c+b, values=c(1, eps, eps, 1),levels=yn)
plist <- compileCPT(list(a, b.a, c.b))
bn   <- grain(plist)
tt   <- querygrain(bn, type="joint")
ftable(tt)
querygrain(setEvidence(bn, evidence=list(a="no", c="yes")))

## ----------------------------------------------------------------------------------
joint <- tabListMult(chest_cpt)
dim(joint)
joint  %>% as.data.frame.table %>% head

## ----------------------------------------------------------------------------------
tabMarg(joint, "lung")
tabMarg(joint, "bronc")

## ----------------------------------------------------------------------------------
ev <- list(asia="yes", dysp="yes")
cond1 <- tabSlice(joint, slice=ev)
cond1 <- cond1 / sum(cond1)
dim(cond1)
tabMarg(cond1, "lung")
tabMarg(cond1, "bronc")

## ----------------------------------------------------------------------------------
cond2 <- tabSliceMult(joint, slice=ev)
cond2 <- cond2 / sum(cond2)
dim(cond2)
tabMarg(cond2, "lung")
tabMarg(cond2, "bronc")

## ----------------------------------------------------------------------------------
yn <- c("yes","no")
a    <- cptable(~asia, values=c(1,99),levels=yn)
t.a  <- cptable(~tub|asia, values=c(5,95,1,99),levels=yn)

(plist1 <- compileCPT(list(a, t.a)))
plist1[[1]]
plist1[[2]]
(chest1 <- grain(plist1))
querygrain(chest1)

## ----------------------------------------------------------------------------------
setEvidence(chest1, evidence=list(asia="yes"))
setEvidence(chest1, nodes="asia", states="yes")
## setFinding(chest1, nodes="asia", states="yes")

## ----------------------------------------------------------------------------------
querygrain(setEvidence(chest1, evidence=list(asia="yes")))

## ----------------------------------------------------------------------------------
g.a <- parray(c("guess.asia", "asia"), levels=list(yn, yn),
              values=c(.8,.2, .1,.9))

## ----------------------------------------------------------------------------------
(plist2 <- compileCPT(list(a, t.a, g.a )))
(chest2 <- grain(plist2))
querygrain( chest2 )

## ----------------------------------------------------------------------------------
querygrain(setEvidence(chest2, evidence=list(guess.asia="yes")))

## ----------------------------------------------------------------------------------
querygrain(setEvidence(chest1, evidence=list(asia=c(.8, .1))))

## ----------------------------------------------------------------------------------
querygrain(setEvidence(chest1, evidence=list(asia=c(1, 0))))

## ----------------------------------------------------------------------------------
dG  <- dag(~A:B + B:C)
uG  <- ug(~A:B + B:C)
par(mfrow=c(1,2)); plot( dG ); plot( uG )

## ----------------------------------------------------------------------------------
dat <- tabNew(c("A", "B", "C"), levels=c(2, 2, 2), values=c(0, 0, 2, 3, 1, 2, 1, 4))
class(dat)

## ----------------------------------------------------------------------------------
gr.dG <- compile( grain( dG, data=dat ) )
gr.uG <- compile( grain( uG, data=dat ) )

## ----------------------------------------------------------------------------------
extractCPT(dat, dG)
c(extractPOT(dat, uG ))

## ----------------------------------------------------------------------------------
p.A.g.B <- tableDiv(dat, tableMargin(dat, "B"))
p.B     <- tableMargin(dat, "B") / sum(dat)
p.AB    <- tableMult( p.A.g.B, p.B)

## ----------------------------------------------------------------------------------
e <- 1e-2
(dat.e <- dat + e)

## ----------------------------------------------------------------------------------
pe.A.g.B <- tableDiv(dat.e, tableMargin(dat, "B"))
pe.B <- tableMargin(dat.e, "B")/sum(dat.e)
pe.AB  <- tableMult( pe.A.g.B, pe.B )

## ----------------------------------------------------------------------------------
dat.e / sum(dat.e)

## ----------------------------------------------------------------------------------
gr.dG <- compile(grain(dG, data=dat, smooth=e))

## ----------------------------------------------------------------------------------
extractCPT(dat, dG, smooth=e)

## ----------------------------------------------------------------------------------
querygrain(gr.dG)
querygrain(gr.uG)

## ----------------------------------------------------------------------------------
querygrain(setFinding(gr.dG, nodes="B", states="B1"))
querygrain(setFinding(gr.uG, nodes="B", states="B1"))

## ----------------------------------------------------------------------------------
gr.uG <- compile(grain(uG, data=dat, smooth=e))

## ----------------------------------------------------------------------------------
c(extractPOT(dat, uG, smooth=e))

## ----------------------------------------------------------------------------------
querygrain(gr.uG)
querygrain(gr.dG)

## ----------------------------------------------------------------------------------
querygrain(setFinding(gr.uG, nodes="B", states="B1"))
querygrain(setFinding(gr.dG, nodes="B", states="B1"))

