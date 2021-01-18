
ROUTING = c(all = 0L, last = 1L)

############### Miscell. ##################


# Expected distribution on booklet given one ability theta
# pscore_mst <- function(theta, prms, booklet_id)
# {
#   if (!booklet_id %in% names(prms$inputs$bkList)) stop("booklet_id in pscore_mst not found")
#   first=prms$inputs$bkList[[booklet_id]]$first
#   last=prms$inputs$bkList[[booklet_id]]$last
#   routing = prms$inputs$bkList[[booklet_id]]$routing
#   b = prms$mst_est$b
#   a = prms$mst_est$item_score
#   
#   nMod=length(first)
#   g_list=vector(mode='list', length=nMod)
#   range_list=vector(mode='list', length=nMod)
#   min_scores = prms$mst_inputs$bkList[[booklet_id]]$min_scores
#   max_scores = prms$mst_inputs$bkList[[booklet_id]]$max_scores
#   for (m in 1:nMod)
#   {
#     g_list[[m]]=elsym0(b, a, first[[m]], last[[m]])
#     range_list[[m]]=min_scores[m]:max_scores[m]
#   }
#   g=elsym_submerge(g_list,range_list,routing)
#   
#   p=rep(0, length(g))
#   for (s in 1:length(g))
#   {
#     p[s]=g[s]*exp((s-1)*theta)
#   }
#   return(p/sum(p))
# }

# expected score distribution given ability theta (possibly a vector)
# for a single booklet
.pscore_mst = function(theta, a, b, first, last, modules, routing)
{
  mxsc = if(routing=='last') sum(modules$module_exit_score_max) else last(modules$module_exit_score_max)
  
  g = elsym_C(ROUTING[[routing]], b, a, as.integer(first-1L), as.integer(last-1L),
                modules$module_exit_score_min, modules$module_exit_score_max,
                modules$nit, mxsc)

  
  p = exp(theta %*% t(0:mxsc)) %*% diag(as.vector(g))
  diag(1/rowSums(p),nrow(p)) %*% p
}

# e.g. p_score user level function
pscore_mst = function(f, test_id, booklet_id)
{
  a = f$mst_inputs$ssIS$item_score
  b = f$mst_est$b
  # dit mag netter, parms object beetje aanpassen
  design = filter(f$mst_inputs$design,
                  .data$test_id==!!test_id & .data$booklet_id==!!booklet_id)
  
  bid = design$bid[1]
  modules = filter(f$mst_inputs$modules, .data$bid==!!bid)
  routing = modules$routing[1]
  
  function(theta)
  {
    .pscore_mst(theta,a,b,design$first,design$last,modules,routing) 
  }
}


E_score_mst <- function(theta, prms, booklet_id)
{
  p = pscore_mst(theta, prms, booklet_id)
  E = 0
  for (s in 1:length(p)) E = E + (s-1)*p[s]
  return(E)
}



############### Calibration functions



# Enorm calibration of MST booklets

# @param fixed_b: vector indicating which parameters re fixed (value) and free to estimate (=NA)
# Default value of fixed_b is NULL which means no fixed parameters
Calibrate_MST = function(first, last, a, sufI, scoretab, booklets, design, modules, nIter=500, fixed_b=NULL)
{
  b=rep(1,length(sufI))
  EsufI = rep(0, length(sufI))
  H = matrix(0,length(a),length(a))
  ref_cat=1

  design$cfirst = as.integer(design$first-1L)
  design$clast = as.integer(design$last-1L)
  booklets$routing = ROUTING[booklets$routing]
  
  if (is.null(fixed_b)) # if no fixed parameters
  {
    nn = sum(sufI)
      ### Implicit equations
    converged=FALSE
    iter = 0
    pb = txtProgressBar(min=0, max=nIter) 
    #plot(c(0,500),c(sum(sufI)-300,sum(sufI)+300),type='n',xlab='iteratie EM',ylab='sum EsufI')
    #abline(h=sum(sufI))
    while ((!converged) && (iter<=nIter))
    {
      iter = iter + 1

      Expect(b, a, design$cfirst, design$clast, 
             booklets$max_score, booklets$nmod, booklets$routing,
             modules$nit, modules$module_exit_score_min, modules$module_exit_score_max,
             scoretab$N, EsufI)
      #points(iter,sum(EsufI))
      # plot(sufI,EsufI)
      # abline(0,1)
      # plot(sufI-EsufI)
      # abline(h=0)
      # browser()

      b = b*sufI/EsufI
      converged = ((max(abs(sufI-EsufI))/nn) < 1e-05)
      setTxtProgressBar(pb, value=iter)
  
    }
    # identify
    b = b/(b[ref_cat]^(a/a[ref_cat]))

    ### Newton-Raphson 
    converged=FALSE
    scale=2
    while ((!converged) && (iter<=nIter))
    {
      iter=iter+1
      
      NR(b, a, design$cfirst, design$clast, 
         booklets$max_score, booklets$nmod, booklets$routing,
         modules$nit, modules$module_exit_score_min, modules$module_exit_score_max,
         scoretab$N, EsufI, H)
      
      H[ref_cat,]=0; H[,ref_cat]=0
      H[ref_cat,ref_cat]=1
      EsufI[ref_cat]=sufI[ref_cat]
      
      b = b*exp(solve(H*scale,sufI-EsufI))
      #if(!all(is.finite(b))) browser()
      converged = (max(abs(EsufI-sufI))/nn<1e-10)
      
      setTxtProgressBar(pb, value=iter)
      scale = max(1, scale-1)
    }
    close(pb)
  }else # there are fixed parameters
  {
    fixed_set=which(!is.na(fixed_b))
    update_set=which(is.na(fixed_b))
    b[fixed_set]=fixed_b[fixed_set]
    nn=sum(sufI[update_set])
        ### Implicit equations
    converged=FALSE
    iter=0
    pb = txtProgressBar(min=0, max=nIter) # max van progressbar klopt neit helemaal, bij NR kan ie er overheen
    while ((!converged) && (iter<=nIter))
    {
      iter = iter+1
      
      Expect(b, a, design$cfirst, design$clast, 
             booklets$max_score, booklets$nmod, booklets$routing,
             modules$nit, modules$module_exit_score_min, modules$module_exit_score_max,
             scoretab$N, EsufI)
      
      b[update_set] = b[update_set]*sufI[update_set]/EsufI[update_set]
      converged=(max(abs(sufI[update_set]-EsufI[update_set]))/nn<1e-04)
      setTxtProgressBar(pb, value=iter)
    }
    
    ### Newton-Raphson 
    converged=FALSE
    scale=2
    while ((!converged) && (iter<=nIter))
    {
      iter=iter+1
      
      NR(b, a, design$cfirst, design$clast, 
         booklets$max_score, booklets$nmod, booklets$routing,
         modules$nit, modules$module_exit_score_min, modules$module_exit_score_max,
         scoretab$N, EsufI, H)
      
      H[fixed_set,]=0
      H[,fixed_set]=0
      diag(H)[fixed_set]=1
      EsufI[fixed_set]=sufI[fixed_set]
      b = b*exp(solve(H*scale,sufI-EsufI))
      converged=(max(abs(EsufI[update_set]-sufI[update_set]))/nn<1e-10)
      setTxtProgressBar(pb, value=iter)
      scale = max(1, scale-1)
    }
    close(pb)
  }
  if ((!converged) && (iter==nIter)) warning(paste("Not converged in ", nIter, " iterations"))
  
  hh = toDexter(b, a, H, first, last, fixed_b = fixed_b)
  return(list(b=b, eta=-log(b), beta=hh$beta, E=EsufI, O=sufI, se.cml=sqrt(diag(hh$acov.beta)), acov.beta = hh$acov.beta))
}



# scoretab: data.frame: booklet_id, booklet_score, N; ordered by booklet_id, booklet_score; N=0 included
# impossible scores included
# fixed_b: NULL or vector of length(b) with NA's for free parameters
Calibrate_Bayes_MST = function(first, last, a, sufI, scoretab, booklets, design, modules, nIter=1000L, fixed_b=NULL)
{
  from = 300L
  step = 1L

  bfirst = as.integer(design$first -1L)
  blast = as.integer(design$last -1L)

  # item booklet index  
  design$bnr = dense_rank(design$bid)
  itb = as.integer(arrange(design, .data$first,.data$bnr)$bnr-1L)
  #nbr of booklets per item
  itnb = design %>% count(.data$first) %>% arrange(.data$first) %>% pull(.data$n)


  booklets$routing = ROUTING[booklets$routing]
  
  # b is consumed in the process
  b = exp(runif(length(sufI), -1, 1))
  
  fixed_b_vec = fixed_b
  if(is.null(fixed_b))
    fixed_b_vec = rep(NA_real_, length(b))

  bx = calibrate_Bayes(a, as.integer(first-1L), as.integer(last-1L),
                       bfirst, blast, booklets$max_score, booklets$min_score,booklets$nmod, booklets$routing,
                      modules$nit, modules$module_exit_score_min, modules$module_exit_score_max,
                      itb, itnb,
                      sufI, scoretab$N, b, fixed_b_vec, 
                      from, step, as.integer(nIter))
  

  bx = apply(bx,2, toDexter, a=a,first=first,last=last,fixed_b=fixed_b)
  beta = do.call(rbind,lapply(lapply(bx,'[[','beta'),drop))
  
  return(list(beta=beta, first=bx[[1]]$first, last=bx[[1]]$last))
  
}





# Fit Rasch and interaction Model for one booklet
# All arguments for one booklet locally:
#   min_scores, 
#   max_scores: integer vectors denoting the range of possible scores 
#     after each module in this booklet on basis of routing rules
#   routing: type of routing; 'all' or 'last'
#   scoretab: integer vector denoting the number of respndents to achieve each score,
#     the indexes correspond to scores, range always from 0 to the sum of the maximum scores of the items
#     (so independent of routing and independent of impossible scores in case of weird weights)
#   first,
#   last: list of integer vectors, on for each module, denoting first and last postion in the score table ordered by 
#       (item_id, item_score) without the 0 score category
# @param a: vector of item_scores arranged by item_id, item_score, also excluding the 0 score category
# @param sufI: tally of respondents achieving each item_score
# @param sufC: <sum(item_score * booklet_score)>
# @param nIter: max number of iterations

Estim_MST <-function(a, first, last, min_scores, max_scores, sufI, sufC, scoretab, routing)
{
  # to~do: different inputs more like enorm, than this pre-amble can be omitted
  cfirst = as.integer(unlist(first)-1L)
  clast = as.integer(unlist(last)-1L)
  crouting = ROUTING[routing]
  
  if(routing=='all')
  {
    bmax = last(max_scores)
    bmin = last(min_scores)
  } else
  {
    bmax = sum(max_scores)
    bmin = sum(min_scores)
  }  
  
  mnit = sapply(first, length)
  

  out = list(routing=routing, regs=list())
  
  nMod = length(first)
  unlist_first= unlist(first, use.names = F)
  unlist_last= unlist(last, use.names = F)
  indx_ic=order(unlist_first)
  nI=length(unlist_first)
  b=rep.int(1,length(sufI))
  EsufI=sufI
  mm=sum(scoretab)

  ic=rep.int(1,nI)
  ncat = unlist_last - unlist_first+1L
  
  var.ic=vector("numeric", nI)
  HRM=matrix(0,length(b),length(b))
  
  # to~do: kan allemaal veel handiger

  ## Rasch Model

  
  ## Implicit Equations
  converged=2
  while (converged>0.001)
  {
    pi_mat = ittotmat_mst(b, a, rep(ic,ncat), cfirst, clast, bmin, bmax, nMod, crouting,
                          mnit, min_scores, max_scores)
    EsufI=pi_mat%*%scoretab 
    b=b*sufI/EsufI
    converged=(max(abs(sufI-EsufI))/mm)
  }
  
  ## NR per item
  converged=2
  scale=2
  while(converged>0.0001)
  {
    converged=-1
    pi_mat = ittotmat_mst(b, a, rep(ic,ncat), cfirst, clast, bmin, bmax, nMod, crouting,
                          mnit, min_scores, max_scores)
    pi_mat[is.na(pi_mat)]=0

    for (m in 1:nMod)
    {
      for (i in 1:length(first[[m]]))
      {
        upd_set = first[[m]][i]:last[[m]][i]
        pi = pi_mat[upd_set,,drop=FALSE]
        E = sufI[upd_set]-pi%*%scoretab
        H = -pi%*%tcrossprod(diag(scoretab),pi) #(diag(scoretab)%*%t(pi)) 
        diag(H) = pi%*%scoretab+diag(H)
        update = solve(H*scale,E)
        b[upd_set] = b[upd_set]*exp(update)
        converged = max(converged,max(abs(E))/mm) 
        HRM[upd_set,upd_set]=H
      }
    }
    if (converged<1) scale=1
  }

  ## IM
  converged=2
  scale=2
  
  pi_mat = ittotmat_mst(b, a, rep(ic,ncat), cfirst, clast, bmin, bmax, nMod, crouting,
                        mnit, min_scores, max_scores)
  
  out$regs$ctrRM = pi_mat
  out$bRM=drop(b)
  out$cRM=ic
  out$HRM = HRM
  
  while(converged>0.001)
  {
    converged=-1
    
    tel_ic=1
    for (m in 1:nMod)
    {
      for (i in 1:length(first[[m]]))
      {
        ic_indx=indx_ic[tel_ic]
        upd_set=first[[m]][i]:last[[m]][i]
        pi=pi_mat[upd_set,,drop=FALSE]
        E=sufI[upd_set]-pi%*%scoretab
        H=-pi%*%tcrossprod(diag(scoretab),pi) 
        diag(H)=pi%*%scoretab+diag(H)
        
        ncol_pi=ncol(pi); nrow_pi=nrow(pi)
        E=c(E,sufC[ic_indx])    ## note the order!
        H=cbind(H,rep.int(0,nrow(H)))
        H=rbind(H,rep.int(0,ncol(H)))
        k=1
        e0=0; e1=0
        f=matrix(0,nrow_pi,ncol_pi)
        g=matrix(0,nrow_pi,ncol_pi)
        h=0
        for (j in upd_set)
        {
          E[length(E)]=E[length(E)]-a[j]*sum((0:(ncol_pi-1))*scoretab*pi[k,])
          e0=e0+a[j]*pi[k,]
          e1=e1+a[j]^2*pi[k,]
          f[k,]=a[j]*(0:(ncol_pi-1))*pi[k,]
          g[k,]=pi[k,]
          h=h+a[j]*(0:(ncol_pi-1))*pi[k,]
          k=k+1
        }
        H[nrow(H),nrow(H)]=sum((0:(ncol_pi-1))^2*(e1-e0^2)*scoretab)
        for (k in 1:nrow(f))
        {
          H[k,nrow(H)]=sum((f[k,]-g[k,]*h)*scoretab)
          H[nrow(H),k]=H[k,nrow(H)]
        }
        update=solve(H*scale,E)
        b[upd_set]=b[upd_set]*exp(update[-length(update)])
        ic[ic_indx]=ic[ic_indx]*exp(update[length(update)]) ## note the order
        var.ic[ic_indx]=solve(H)[nrow(H),nrow(H)] ## note the order
        tel_ic=tel_ic+1
        converged=max(converged,max(abs(E))/mm)
      }
    }
    pi_mat = ittotmat_mst(b, a, rep(ic,ncat), cfirst, clast, bmin, bmax, nMod, crouting,
                          mnit, min_scores, max_scores)
    
    if (converged<1) scale=1
  }
  out$regs$ctrIM = pi_mat
  out$bIM = drop(b)
  out$cIM = ic
  out$se.c = sqrt(var.ic)
  out$fit.stats = log(ic)/sqrt(var.ic)
  
  out
}

## Generate response NRM

# rNRM=function(theta, b, a, first, last)
# {
#   sampleNRM(as.double(theta),
#              as.double(b), 
#              as.integer(a), 
#              as.integer(first-1),
#              as.integer(last-1))
# }


################################### Abilities

#' #' Fisher information function for each path/booklet in a test
#' #'
#' #' @param db an dextermst db handle
#' #' @param parms a parms object producted by \code{\link{fit_enorm_mst}}
#' #' @param test_id the id of a test
#' #' @param theta vector of abilities for which information is required
#' #'
#' path.information = function(db, parms, test_id, theta)
#' {
#'   bks = parms$mst_inputs$booklet_design %>%
#'     filter(test_id == test_id) %>%
#'     distinct(booklet_id) %>%
#'     pull(booklet_id)
#'   nT = length(theta)
#'   nBk = length(bks)
#'   out = matrix(0,nBk, nT)
#'   for (i in 1:nBk)
#'   {
#'     bk_id = paste0(test_id,".",bks[i])
#'     first = sort(unlist(parms$inputs$bkList[[bk_id]]$first, use.names = F))
#'     last  = sort(unlist(parms$inputs$bkList[[bk_id]]$last, use.names = F))
#'     out[i,] = dexter.IJ(parms$est$b, parms$est$a, first, last, theta)
#'   }
#'   out
#' }

#### Estimate ability
#### ML estimation of theta
# uses an implicit equations algorithm
# theta_MLE_MST <- function(prms, booklet_id)
# {
#   routing = prms$mst_inputs$bkList[[booklet_id]]$routing
#   if (routing=="none")
#   {
#     b = prms$est$b
#     a = prms$est$a
#     first = prms$inputs$ssI$first
#     last = prms$inputs$ssI$last
#     theta = dexter:::theta_MLE(b,a,first, last)
#   }else
#   {
#     a = prms$mst_est$item_score
#     nMod=length(prms$mst_inputs$bkList[[booklet_id]]$first)
#     if (routing=="all"){
#       mxs.a=prms$mst_inputs$bkList[[booklet_id]]$max_scores[nMod]
#       mns.a=prms$mst_inputs$bkList[[booklet_id]]$min_scores[nMod]
#     }
#     if (routing=="last"){
#       mns.a=sum(prms$mst_inputs$bkList[[booklet_id]]$min_scores)
#       mxs.a=sum(prms$mst_inputs$bkList[[booklet_id]]$max_scores)
#     }
#     ms.a=sum(a[unlist(prms$mst_inputs$bkList[[booklet_id]]$last, use.names = F)])
#     theta=rep(NA, ms.a+1)
#     for (s in max(1,mns.a):(min(mxs.a,ms.a)-1))
#     {
#       escore=-1
#       theta[s]=0
#       while (abs(escore-s)>1e-1)
#       {
#         escore = E_score_mst(theta[s],prms, booklet_id)
#         theta[s] = theta[s]+log(s/escore)
#       }
#     }
#     theta=c(-Inf,theta,Inf)
#   }
#   return(theta)
# }
