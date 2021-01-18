context("Apply a function to a dimension of an array")
# list of acts: (name: act for arrApply; value: if NULL, the same  act is for apply(),
#   if list, its fields are: rf -> name for equivalent r function, args -> list
#   of ... in apply(), argc -> ... for arrApply())
set.seed(7)
n=3
v=rnorm(n)
m=matrix(rnorm(n*n), n)
d3=rep(n, 3)
ar3d=array(rnorm(prod(d3)), dim=d3)
d4=rep(n, 4)
ar4d=array(rnorm(prod(d4)), dim=d4)
vp=v
lacts=list("sum"=NULL, "prod"=NULL, "all"=NULL, "any"=NULL, "min"=NULL, "max"=NULL, "mean"=NULL, "median"=NULL, "sd"=NULL, "var"=NULL, "cumsum"=NULL, "cumprod"=NULL, "diff"=NULL,
   # translated acts
   norm=list(rf="norm", argr=list(type='2'), argc=list(p=2)),
   trapz=list(rf=function(v) {n=length(v); return(sum(v)-0.5*(v[1]+v[n]))}),
   normalise=list(rf=function(v) v/norm(v, '2'), argc=list(p=2)),
   multv=list(rf=function(v, vv) v*vv, argr=list(vv=vp), argc=list(v=vp)), 
   divv=list(rf=function(v, vv) v/vv, argr=list(vv=vp), argc=list(v=vp)), 
   addv=list(rf=function(v, vv) v+vv, argr=list(vv=vp), argc=list(v=vp)), 
   subv=list(rf=function(v, vv) v-vv, argr=list(vv=vp), argc=list(v=vp))
)
test_ar=function(ar, tol=1.e-14, acts=lacts, ndi=seq_along(dim(ar)), vp=v) {
    # compare arrApply() to translated r call
    vec=FALSE
    if (length(ndi) == 0) {
        # we have a vector
        ndi=1
        vec=TRUE
    }
    for (act in names(acts)) {
        ract=acts[[act]]
        rfu=if (is.null(ract)) act else ract$rf
        for (idim in ndi) {
            r1=do.call(arrApply, c(list(ar, idim, act), ract$argc))
            if (!vec && length(dim(r1)) == length(dim(ar))) {
                # permute to the same order as in apply
                r1=aperm(r1, c(idim, ndi[-idim]))
            }
            if (vec) {
                r2=suppressWarnings(do.call(apply, c(list(as.matrix(ar), 2, rfu), ract$argr)))
            } else {
                r2=suppressWarnings(do.call(apply, c(list(ar, ndi[-idim], rfu), ract$argr)))
            }
            expect_equal(as.numeric(r1), as.numeric(r2), tolerance=tol, scale=1, info=sprintf("'%s' on idim=%d in dims=(%s)", act, idim, paste(if (vec) length(ar) else dim(ar), collapse=", ")))
        }
    }
}
test_that("arrApply on a vector", {
    test_ar(v)
})
test_that("arrApply on a matrix", {
    test_ar(m)
})
test_that("arrApply on an array 3D", {
    test_ar(ar3d)
})
test_that("arrApply on an array 4D", {
    test_ar(ar4d)
})
