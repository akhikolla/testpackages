###########################################################################
############## FUNCTIONS FOR REFERENCE BY MAIN - NOT FOR USER #############
###########################################################################

#                  o.fun, all PRO and var.est FUNCTIONS                        #
#_______________________________________________________________________________
# Original by Chihyun Lee (August, 2017)                                       #
# Last Modified by Sandra Castro-Pearson (April, 2018)                         #
# Received from Chihyun Lee (January, 2018)                                    #
#_______________________________________________________________________________

#' @useDynLib BivRec
#' @keywords internal

r2f.pro.ee1 <- function(n, nparams, di, xmati, gmati, L, expA, subsum, kcount){
  out1 <- .Fortran("xmproee",
                   n=as.integer(n),
                   nparams=as.integer(nparams),
                   di=as.double(di),
                   xmati=as.double(xmati),
                   gmati=as.double(gmati),
                   L=as.double(L),
                   expA=as.double(expA),
                   subsum=as.double(subsum),
                   kcount=as.integer(kcount),
                   NAOK = FALSE, PACKAGE = "BivRec")

  subsum <- out1$subsum

  return(subsum)
}


#' @useDynLib BivRec
#' @keywords internal

r2f.pro.ee2 <- function(n, nparams, di, xmati, ymati, gmati, L, expA, subsum, kcount){
  out2 <- .Fortran("ymproee",
                   n=as.integer(n),
                   nparams=as.integer(nparams),
                   di=as.double(di),
                   xmati=as.double(xmati),
                   ymati=as.double(ymati),
                   gmati=as.double(gmati),
                   L=as.double(L),
                   expA=as.double(expA),
                   subsum=as.double(subsum),
                   kcount=as.integer(kcount),
                   NAOK = FALSE, PACKAGE = "BivRec")

  subsum <- out2$subsum

  return(subsum)
}


#' @useDynLib BivRec
#' @keywords internal

r2f.pro.var <- function(n, nparams, xmat, ymat, gmatx, gmaty, l1, l2,
                        expAx, expAy, subsumx, subsumy, dx, dy, mstar, mc){
  out <- .Fortran("mprovar",
                  n=as.integer(n),
                  nparams=as.integer(nparams),
                  xmati=as.double(xmat),
                  ymati=as.double(ymat),
                  gmatx=as.double(gmatx),
                  gmaty=as.double(gmaty),
                  l1=as.double(l1),
                  l2=as.double(l2),
                  expAx=as.double(expAx),
                  expAy=as.double(expAy),
                  subsumx=as.double(subsumy),
                  subsumy=as.double(subsumy),
                  dx=as.double(dx),
                  dy=as.double(dy),
                  mstar=as.double(mstar),
                  mc=as.integer(mc),
                  NAOK = FALSE, PACKAGE = "BivRec")

  subsum1 <- out$subsumx
  subsum2 <- out$subsumy

  return(cbind(subsum1, subsum2))
}

###########################################################################
############## FUNCTIONS FOR REFERENCE BY MAIN - NOT FOR USER #############
###########################################################################

#                 o.fun, all MPRO and MVAR FUNCTIONS                           #
#_______________________________________________________________________________
# Original by Chihyun Lee (August, 2017)                                       #
# Modified to Fortran by Sandra Castro-Pearson (last updated July, 2018)       #
# Received from Chihyun Lee (January, 2018)                                    #
#_______________________________________________________________________________


#' @useDynLib BivRec
#' @keywords internal

r2f.mpro.ee1 <- function(n, nparams, di, xmati, gmati, L, expA, subsum, kcount){
  out1 <- .Fortran("xmproee",
                   n=as.integer(n),
                   nparams=as.integer(nparams),
                   di=as.double(di),
                   xmati=as.double(xmati),
                   gmati=as.double(gmati),
                   L=as.double(L),
                   expA=as.double(expA),
                   subsum=as.double(subsum),
                   kcount=as.integer(kcount),
                   NAOK = FALSE, PACKAGE = "BivRec")

  subsum <- out1$subsum

  return(subsum)
}

#' @useDynLib BivRec
#' @keywords internal

r2f.mpro.ee2 <- function(n, nparams, di, xmati, ymati, gmati, L, expA, subsum, kcount){
  out2 <- .Fortran("ymproee",
                   n=as.integer(n),
                   nparams=as.integer(nparams),
                   di=as.double(di),
                   xmati=as.double(xmati),
                   ymati=as.double(ymati),
                   gmati=as.double(gmati),
                   L=as.double(L),
                   expA=as.double(expA),
                   subsum=as.double(subsum),
                   kcount=as.integer(kcount),
                   NAOK=FALSE, PACKAGE = "BivRec")

  subsum <- out2$subsum

  return(subsum)
}


#' @useDynLib BivRec
#' @keywords internal

r2f.mpro.var <- function(n, nparams, xmat, ymat, gmatx, gmaty, l1, l2,
                         expAx, expAy, subsumx, subsumy, dx, dy, mstar, mc){
  out <- .Fortran("mprovar",
                  n=as.integer(n),
                  nparams=as.integer(nparams),
                  xmati=as.double(xmat),
                  ymati=as.double(ymat),
                  gmatx=as.double(gmatx),
                  gmaty=as.double(gmaty),
                  l1=as.double(l1),
                  l2=as.double(l2),
                  expAx=as.double(expAx),
                  expAy=as.double(expAy),
                  subsumx=as.double(subsumy),
                  subsumy=as.double(subsumy),
                  dx=as.double(dx),
                  dy=as.double(dy),
                  mstar=as.double(mstar),
                  mc=as.integer(mc),
                  NAOK = FALSE, PACKAGE = "BivRec")

  subsum1 <- out$subsumx
  subsum2 <- out$subsumy

  return(cbind(subsum1, subsum2))
}

#FIRST FUNCTION CALLS ON COMPILED ONESAMP FORTRAN CODE

#                                 FORTRAN CODE                                 #
#______________________________________________________________________________#
# Original by Shu-Hui Chang                                                    #
# Modified by Chiung-Yu for bivariate recurrent event - Fortran (Feb, 2001)    #
# Modified and compiled for package by Sandra Castro-Pearson (July, 2018)      #
# Received from Xianghua Luo (May, 2018)                                       #
#______________________________________________________________________________#


#' @useDynLib BivRec
#' @keywords internal

r_onesamp <- function(n,gtime,ctime,mc,m,
                      cen,ucen,nd,udt,tot,gap,event,
                      r,d,sest,std){

  out1 <- .Fortran("onesamp",
                   n=as.integer(n),
                   gtime=as.double(gtime),
                   ctime=as.double(ctime),
                   count=as.double(m),
                   mc=as.integer(mc),
                   m=as.integer(m),
                   cen=as.double(cen),
                   ucen=as.double(ucen),
                   nd=as.integer(nd),
                   udt=as.double(udt),
                   tot=as.integer(tot),
                   gap=as.double(gap),
                   event=as.double(event),
                   r=as.double(r),
                   d=as.double(d),
                   sest=as.double(sest),
                   std= as.double(std),
                   NAOK = FALSE, PACKAGE = "BivRec")

  out2 <- data.frame(time = out1$udt, surv = out1$sest, std = out1$std)

  return(out2)
}

#FIRST FUNCTION CALLS ON COMPILED BIVRECUR FORTRAN CODE

#                                 FORTRAN CODE                                 #
#______________________________________________________________________________#
# Original by Shu-Hui Chang                                                    #
# Modified by Chiung-Yu for bivariate recurrent event - Fortran (Feb, 2001)    #
# Modified and compiled for package by Sandra Castro-Pearson (July, 2018)      #
# Received from Xianghua Luo (May, 2018)                                       #
#______________________________________________________________________________#

#' @useDynLib BivRec
#' @keywords internal

r_bivrecur <- function(n, gtime, ctime, mc, m,
                       cen, ucen, nd, udt, tot, gap, event,
                       r, d, sest, var, markvar1, markvar2,
                       mark1, mark2, u1, u2, Fest, tmpindex, prob, std){

  out1 <- .Fortran("bivrecur",
                   n=as.integer(n),
                   gtime=as.double(gtime),
                   ctime=as.double(ctime),
                   count=as.double(m),
                   mc=as.integer(mc),
                   m=as.integer(m),
                   cen=as.double(cen),
                   ucen=as.double(ucen),
                   nd=as.integer(nd),
                   udt=as.double(udt),
                   tot=as.integer(tot),
                   gap=as.double(gap),
                   event=as.double(event),
                   r=as.double(r),
                   d=as.double(d),
                   sest=as.double(sest),
                   var=as.double(var),
                   markvar1=as.double(markvar1),
                   markvar2=as.double(markvar2),
                   mark1=as.double(mark1),
                   mark2=as.double(mark2),
                   u1=as.double(u1),
                   u2=as.double(u2),
                   Fest=as.double(Fest),
                   tmpindex=as.integer(tmpindex),
                   prob=as.double(prob),
                   std= as.double(std),
                   NAOK = FALSE, PACKAGE = "BivRec")

  out2 <- c(prob = out1$prob, std = out1$std)

  return(out2)
}


