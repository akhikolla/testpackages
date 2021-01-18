## size_mb <- function(x) {
##   format(object.size(x), units = "Mb", standard = "auto", digits = 1L)
## }

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                 MUNIN: 1041-1397-80592
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## l    <- readRDS("../../../../sandbox/r/bns/munin.rds")
## cpts <- bnfit_to_cpts(l)
## cl   <- cpt_list(cpts)

## tictoc::tic()
## cp_munin   <- compile(cl, save_graph = TRUE)
## tictoc::toc()

## size_mb(cp)

## tictoc::tic()
## j <- jt(cp)
## tictoc::toc()

## size_mb(j)
## .map_dbl(j$charge$C, sum)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                         LINK: 724-1125-14211
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## l    <- readRDS("../../../../sandbox/r/bns/link.rds")
## cpts <- bnfit_to_cpts(l)
## cl   <- cpt_list(cpts)
## cp   <- compile(cl, save_graph = TRUE)

## tictoc::tic()
## j  <- jt(cp, propagate = "full")
## tictoc::toc()

## size_mb(j)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                        DIABETES: 413-602-429409
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## l    <- readRDS("../../../../sandbox/r/bns/diabetes.rds")
## cpts <- bnfit_to_cpts(l)
## cl   <- cpt_list(cpts)

## tictoc::tic()
## cp   <- jti::compile(cl)
## tictoc::toc() # 40 s
## size_mb(cp) # 5 mb

## tictoc::tic()
## j    <- jt(cp)
## tictoc::toc() # 13 s

## size_mb(j) # 121 mb

## .map_dbl(j$charge$C, function(x) sum(x))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                      MILDEW: 35-46-540150
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## l <- readRDS("../../../../sandbox/r/bns/mildew.rds")
## cpts <- bnfit_to_cpts(l)
## cl   <- cpt_list(cpts)
## cp   <- compile(cl, save_graph = TRUE)
## j    <- jt(cp, propagate = "full")
## .map_dbl(j$charge$C, function(x) sum(x))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                      HAILFINDER: 56-66-2656
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## l    <- readRDS("../../../../sandbox/r/bns/hailfinder.rds")
## cpts <- bnfit_to_cpts(l)
## cl   <- cpt_list(cpts)
## cp   <- jti::compile(cl)
## j <- jt(cp, propagate = "full")
## .map_dbl(j$charge$C, function(x) sum(x))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                      BARLEY: 48-84-114005
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## l    <- readRDS("../../../../sandbox/r/bns/barley.rds")
## cpts <- bnfit_to_cpts(l)
## cl   <- cpt_list(cpts)
## cp   <- jti::compile(cl)
## j    <- jt(cp)
## .map_dbl(j$charge$C, function(x) sum(x))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                      CHILD: 20-25-230
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## l    <- readRDS("../../../../sandbox/r/bns/child.rds")
## cpts <- bnfit_to_cpts(l)
## cl   <- cpt_list(cpts)
## cp   <- jti::compile(cl, save_graph = TRUE)
## j    <- jt(cp)



## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                      HEPAR2: 70-123-1453
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## l    <- readRDS("../../../../sandbox/r/bns/hepar2.rds")
## cpts <- bnfit_to_cpts(l)
## cl   <- cpt_list(cpts)
## cp   <- jti::compile(cl)
## j <- jt(cp)
## .map_dbl(j$charge$C, function(x) sum(x))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                      INSURANCE: 27-52-984
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## l    <- readRDS("../../../../sandbox/r/bns/insurance.rds")
## cpts <- bnfit_to_cpts(l)
## cl   <- cpt_list(cpts)
## cp   <- jti::compile(cl)
## j    <- jt(cp)
## .map_dbl(j$charge$C, function(x) sum(x))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                      CANCER: 5-4-10
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## l    <- readRDS("../../../../sandbox/r/bns/cancer.rds")
## cpts <- bnfit_to_cpts(l)
## cl   <- cpt_list(cpts)
## cp   <- compile(cl, save_graph = TRUE)
## j    <- jt(cp)
## .map_dbl(j$charge$C, function(x) sum(x))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                      EARTHQUAKE: 5-4-10
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## l    <- readRDS("../../../../sandbox/r/bns/earthquake.rds")
## cpts <- bnfit_to_cpts(l)
## cl   <- cpt_list(cpts)
## cp   <- compile(cl, save_graph = TRUE)
## j    <- jt(cp)
## .map_dbl(j$charge$C, function(x) sum(x))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                      SACHS: 11-17-178
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## l <- readRDS("../../../../sandbox/r/bns/sachs.rds")
## cpts <- bnfit_to_cpts(l)
## cl   <- cpt_list(cpts)
## cp   <- compile(cl, save_graph = TRUE)
## j    <- jt(cp)
## .map_dbl(j$charge$C, function(x) sum(x))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                      ALARM: 37-46-509
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## l    <- readRDS("../../../../sandbox/r/bns/alarm.rds")
## cpts <- bnfit_to_cpts(l)
## cl   <- cpt_list(cpts)
## cp   <- jti::compile(cl)
## j    <- jt(cp)
## .map_dbl(j$charge$C, function(x) sum(x))
