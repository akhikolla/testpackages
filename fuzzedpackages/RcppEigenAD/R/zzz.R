
.onLoad<-function(libname, pkgname)
{

setMethod("%.%",signature(f="function",g="function"),
    function(f,g)
    {
      return(
      function(x,order=0)
      {
        if(order== 0)
        {
          return(f(g(x)))
        }
       if(order== 1)
        {
           return(marshal_faa_di_bruno(J(f)(g(x)),H(f)(g(x)),J(g)(x),H(g)(x))$Jfog)
        }
        if(order== 2)
        {
           return(marshal_faa_di_bruno(J(f)(g(x)),H(f)(g(x)),J(g)(x),H(g)(x))$Hfog)
        }
      }
      )        
    }
)


}
