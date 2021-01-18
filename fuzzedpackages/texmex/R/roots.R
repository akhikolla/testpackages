roots <-
function(lev,a,c,b,d,Zj,Zk) 
  {
    #--------------------------------------------------------------------
    #Children Functions to be Used in roots().
    Dderiv <-
      function(x,alj,blj,aki,bki,Zlj,Zki)
      {
        res <- (alj-aki) + (blj*Zlj*(x^(blj-1))) - (bki*Zki*(x^(bki-1)))
        return(res)
      }

    Dderiv_01 <-
      function(x,alj,blj,aki,bki,Zlj,Zki)
      {
        res <- (alj-aki) + (blj*Zlj*((-log(x))^(blj-1))) -
          (bki*Zki*((-log(x))^(bki-1)))
        return(res)
      }
    
    is.wholenumber <-
      function(x, tol = .Machine$double.eps^0.14)  abs(x - round(x)) < tol
    #--------------------------------------------------------------------
    
    cond_na      <- NULL
    z            <- list()
    s            <- lev
    xstar        <- lev
    xdstar       <- lev
    no_of_roots  <- 0
    
    if( (b!=0) & (d!=0) & (b!=d) & (a > c) & (Zj != 0) & (Zk != 0) )
      {
        sbase        <- (d*(d-1)*Zk)/(b*(b-1)*Zj)    
        spower       <- 1/(b-d)        
        if(sbase <= 0 & is.wholenumber((spower)/2)==FALSE) s="complex" 
        if(sbase>0 || ((is.wholenumber((spower)/2)==TRUE) & spower > 0))
          {
            s=sbase^(spower)
            cond_na <- is.na(s)

            if(cond_na==TRUE || (cond_na==FALSE & (s<=lev || s == Inf)))
              s = "complex"                        
          }

        Dprimev  <- Dderiv_01(x=exp(-lev),alj=a,blj=b,aki=c,
                              bki=d,Zlj=Zj,Zki=Zk)

        Dinf     <- Dderiv_01(x=0,alj=a,blj=b,aki=c,
                              bki=d,Zlj=Zj,Zki=Zk)
        #print(paste(s))
        if((s=="complex")==FALSE)
          {
            
            Dprimes  <- Dderiv_01(x=exp(-s),alj=a,blj=b,
                                  aki=c,bki=d,Zlj=Zj,Zki=Zk)
            if(Dprimev>0 & Dprimes>0)
              {
                no_of_roots <- 0
              }            
            if(Dprimev<0 & Dinf>0)
              {
                no_of_roots <- 1
                
                xstar     <- -log(uniroot(Dderiv_01,
                                          interval=c(exp(-lev),0),
                                          a,b,c,d,Zj,Zk)$root)          
              }
            
            if( Dprimev > 0 & Dprimes < 0 & Dinf > 0)
              {
                no_of_roots <- 2
                
                xstar     <- -log(uniroot(Dderiv_01,
                                          interval=c(exp(-lev),exp(-s)),
                                          a,b,c,d,Zj,Zk)$root)
                
                xdstar    <- -log(uniroot(Dderiv_01,
                                          interval=c(exp(-s),0),
                                          a,b,c,d,Zj,Zk)$root)           
              }
            
          }
        if((s=="complex")==TRUE)
          {
            if(Dprimev>0)  no_of_roots <- 0
            
            if(Dprimev<0 & Dinf > 0 )  
              {                
                no_of_roots <- 1              
                xstar       <- -log(uniroot(Dderiv_01,
                                            interval=c(exp(-lev),0),
                                            a,b,c,d,Zj,Zk)$root)        
              }            
          }
      }
    if(b==0  ||  d==0 || (b==d) || (a==c) ||
       (Zj == 0) || (Zk == 0) )
      {
        
        if(b==0 || Zj==0)
          {
            if((d!=0 & Zk!=0) & ((a-c) > 0))
              {
                if(((a-c)/(d*Zk))>0)
                  {
                    xstar <- ((a-c)/(d*Zk))^(1/(d-1))
                    if(xstar <= lev || xstar == Inf) xstar <- lev
                    if(xstar > lev   & xstar != Inf) no_of_roots <- 1
                  }
              }
            if((d==0 || Zk ==0)) no_of_roots <- 0
          }
        
        if(d==0 || Zk==0)
          {
            if((b!=0 & Zj!=0) & ((a-c) > 0))
              {
                if(((c-a)/(b*Zj))>0)
                  {
                    xstar <- ((c-a)/(b*Zj))^(1/(b-1))
                    if( xstar <= lev || xstar == Inf ) xstar <- lev
                    if( xstar > lev   & xstar != Inf ) no_of_roots <- 1
                  }
              }
            if((d==0 || Zk ==0)) no_of_roots <- 0            
          }
        
        if( ((a-c) == 0) & (b!=0 & d!=0) & (b!=d))
          {
            if( ((d*Zk)/(b*Zj)) > 0 )
              {
                xstar <- ((d*Zk)/(b*Zj))^(1/(b-d))
                if( xstar <= lev || xstar == Inf ) xstar <- lev
                if( xstar > lev   & xstar != Inf ) no_of_roots <- 1
              }
          }        
      }    
    z$no     <- no_of_roots
    z$s      <- s
    z$xstar  <- xstar
    z$xdstar <- xdstar       
    return(z)
  }
