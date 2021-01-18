## ----setup, echo=FALSE--------------------------------------------------------
library(hadron)

## -----------------------------------------------------------------------------
data(correlatormatrix)
correlatormatrix <- bootstrap.cf(correlatormatrix, boot.R=99, boot.l=1, seed=132435)
correlatormatrix.gevp <- bootstrap.gevp(cf=correlatormatrix, t0=4, element.order=c(1,2,3,4))
pc1 <- gevp2cf(gevp=correlatormatrix.gevp, id=1)

## -----------------------------------------------------------------------------
pc1.hankel <- bootstrap.hankel(cf=pc1, t0=2, n=2)

## -----------------------------------------------------------------------------
hpc1 <- hankel2cf(hankel=pc1.hankel, id=1)

## -----------------------------------------------------------------------------
plot(hpc1, log="y", ylab="lambda(delta t)", xlab="delta t + t0")

## -----------------------------------------------------------------------------
heffectivemass1 <- hadron:::hankel2effectivemass(hankel=pc1.hankel, id=1)

## -----------------------------------------------------------------------------
pc1.effectivemass <- bootstrap.effectivemass(cf=pc1)

## -----------------------------------------------------------------------------
plot(pc1.effectivemass, pch=21, col="red", ylim=c(0,1.1), xlim=c(0,18),
     xlab="t", ylab="M(t)")
plot(heffectivemass1, rep=TRUE, pch=22, col="blue")
legend("topright", legend=c("pc1", "hankel1"), bty="n", pch=c(21,22), col=c("red", "blue"))

## ---- echo=FALSE--------------------------------------------------------------
plot(pc1.effectivemass, pch=21, col="red", ylim=c(0,1.1), xlim=c(0,18),
     xlab="t", ylab="M(t)")
plot(heffectivemass1, pch=22, col="blue", rep=TRUE)
pc1.hankel <- bootstrap.hankel(cf=pc1, t0=1, n=2)
heffectivemass1 <- hadron:::hankel2effectivemass(hankel=pc1.hankel, id=1)
plot(heffectivemass1, pch=23, col="darkgreen", rep="TRUE")
pc1.hankel <- bootstrap.hankel(cf=pc1, t0=3, n=2)
heffectivemass1 <- hadron:::hankel2effectivemass(hankel=pc1.hankel, id=1)
plot(heffectivemass1, pch=24, col="black", rep="TRUE")
legend("topright", legend=c("pc1", "t0=2", "t0=1", "t0=3"), bty="n", pch=c(21:24),
       col=c("red", "blue", "darkgreen", "black"))

## -----------------------------------------------------------------------------
ppcor <- extractSingleCor.cf(cf=correlatormatrix, id=1)
ppcor.effectivemass <- bootstrap.effectivemass(cf=ppcor)
ppcor.hankel <- bootstrap.hankel(cf=ppcor, t0=3, n=2)
heffectivemass1 <- hadron:::hankel2effectivemass(hankel=ppcor.hankel, id=1)
plot(ppcor.effectivemass, pch=21, col="red", ylim=c(0,1.1), xlim=c(0,18),
     xlab="t", ylab="M(t)")
plot(heffectivemass1, pch=22, col="blue", ylim=c(0,1.1), rep=TRUE)
legend("topright", legend=c("ppcor", "hankel1"), bty="n", pch=c(21,22), col=c("red", "blue"))

