.ecgfMix <- function(y, decay, method, deriv = FALSE, m = NULL)
{
  d <- length(y)
  
  rss <-  drop( sum(y^2) )
  
  if( decay == Inf ){
    
    mix <- 0.0
    DmixDy <- numeric(d)
    
  } else {
    # Choosing method to determine the mixture between gaussian and normal ecgf
    switch(method,
           
           "gaus" = {
             
             if( is.null(m) ) stop("If using method \"gaus\" you need to specify the number of simulations used.")
             
             mix <- exp( - decay * rss / (sqrt(m) * 2) )
             
             if(deriv){
               DmixDy <- drop( (mix*decay / sqrt(d)) * (-y) )
             }
             
           },
           
           "mse" = {
             sadMse <- exp(rss) 
             gausMse <- rss + 0.5*rss^2 + 1
             mix <- (gausMse / sadMse) ^ decay
             
             if(deriv){
               normRes <- y
               scaledRSS <- drop( crossprod(y, normRes) )
               sadVar <- exp(scaledRSS)
               normVar <- scaledRSS + 0.5 * scaledRSS^2 + 1
               d_sadVar <- 2 * sadVar * normRes
               d_normVar <- 2 * normRes * ( 1 + scaledRSS ) 
               DmixDy <- drop( d_normVar * sadVar - normVar * d_sadVar ) / (sadVar^2)
               
               if(normVar / sadVar < 1){
                 DmixDy <- decay * ( (normVar / sadVar) ^ (decay-1) ) * DmixDy 
               } else {
                 DmixDy <- DmixDy * 0.0
               }
             }
             
           },
           stop("\"method\" must be either \"gaus\" or \"mse\"") 
    )
  }
  
  if( !deriv ) DmixDy <- NULL
  
  return( list("mix" = mix, "DmixDy" = DmixDy) )
  
}