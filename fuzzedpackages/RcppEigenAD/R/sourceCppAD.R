

sourceCppAD<-function(code=NULL,file=NULL,wrt=1,output="method")
{

    if(is.null(code) && is.null(file))
    {
        return(NULL)
    }

    if(is.null(code) && !is.null(file))
    {
        code<-read_file(file)
    }
    
fname<-get_function_name(code)


# add the dependencies
# preamble<-'// [[Rcpp::plugins(cpp11)]]' 
preamble<-''
preamble<-paste(preamble,'\n',sep="")
preamble<-paste(preamble,'// [[Rcpp::depends(RcppEigen)]]',sep="")
preamble<-paste(preamble,'\n',sep="")
preamble<-paste(preamble,'// [[Rcpp::depends(RcppEigenAD)]]',sep="")
preamble<-paste(preamble,'\n',sep="")
preamble<-paste(preamble,'// [[Rcpp::depends(BH)]]',sep="")
preamble<-paste(preamble,'\n',sep="")
preamble<-paste(preamble,'\n',sep="")
preamble<-paste(preamble,'#include <Rcpp.h>')
preamble<-paste(preamble,'\n',sep="")
preamble<-paste(preamble,'#include <RcppEigen.h>')
preamble<-paste(preamble,'\n',sep="")
preamble<-paste(preamble,'#include <unsupported/Eigen/SpecialFunctions>')
preamble<-paste(preamble,'\n',sep="")
preamble<-paste(preamble,'#include "adfunc.h"')
preamble<-paste(preamble,'\n',sep="")
preamble<-paste(preamble,'#include "bound.func.h"')
preamble<-paste(preamble,'\n',sep="")
preamble<-paste(preamble,'\n',sep="")
preamble<-paste(preamble,'\n',sep="")
preamble<-paste(preamble,'\n',sep="")
preamble<-paste(preamble,'boost::function<ADmat (const ADmat&)> bound_func;\n',sep="")

# add user code

adlacode<-paste(preamble,code,sep="")
preamble<-paste(preamble,'\n',sep="")
preamble<-paste(preamble,'\n',sep="")
preamble<-paste(preamble,'\n',sep="")





# add a single argument dispatcher
dispatcher<-'
MULTIARG(Scalar_1) wrapper_1(const MULTIARG(Scalar_1)& X)
{
  MULTIARG(Scalar_1) Y(1);
  Y[0] = '
dispatcher<-paste(dispatcher,'bound_func',sep="")
dispatcher<-paste(dispatcher,'(X[0]); return Y; }',sep="")
dispatcher<-paste(dispatcher,'\n',sep="")


adlacode<-paste(adlacode,dispatcher,sep="")
adlacode<-paste(adlacode,'\n',sep="")
adlacode<-paste(adlacode,'\n',sep="")
adlacode<-paste(adlacode,'\n',sep="")


# create the code to marshal the users method through Rcpp and adla
    
nargs<-length(formals(fname))
marshaller<-paste('// [[Rcpp::export]]\n',sep="")
marshaller<-paste(marshaller,'Rcpp::List marshal(')
marshaller<-paste(marshaller,Reduce(function(x,y) paste(x,y,sep=",\n"),paste('const Eigen::MatrixXd& X_',1:nargs,sep="")),sep="")
marshaller<-paste(marshaller,')\n{\n',sep="")

if(nargs > 1)
{
    vars<-1:nargs
    vars<-vars[-wrt]
    for(var in vars)
    {
        marshaller<-paste(marshaller,'ADmat adX_',var,' = Convert<ADmat>(X_',var,');\n',sep="")
    }

    marshaller<-paste(marshaller,'bound_func = boost::bind(',fname,sep="")
    binding.args<-paste('adX_',1:nargs,sep="")
    binding.args[wrt]<-'_1'
    binding.args<-Reduce(function(x,y) paste(x,y,sep=","),binding.args)
    marshaller<-paste(marshaller,binding.args,sep=",")
    marshaller<-paste(marshaller,');\n',sep="")
}
else
{
    marshaller<-paste(marshaller,'bound_func = boost::function<ADmat (const ADmat&)>(',fname,');\n',sep="")
}


marshaller<-paste(marshaller,'
  MULTIARG(Scalar_0) Xs(1);
  Xs[0] = X_',wrt,';\n',sep="")
  marshaller<-paste(marshaller,'
  FUNCTION(Scalar_1) F(wrapper_1);
  TRIPLE(Scalar_0) T = adfunc(F,Xs);


  MATRIX(Scalar_0) Jy = T.get<1>();
  unsigned int domain_size = Jy.cols();
  unsigned int co_domain_size = Jy.rows();
  // stack the hessians

  MATRIX(Scalar_0) Hy(domain_size,domain_size*co_domain_size);
  for(unsigned int dim = 0; dim < co_domain_size; dim++)
     {
	Hy.block(0,dim*domain_size,domain_size,domain_size) = T.get<2>()(0,dim);
     }
  return Rcpp::List::create(Rcpp::Named("f") = T.get<0>()[0],
			    Rcpp::Named("Jf") = Jy,
			    Rcpp::Named("Hf") = Hy);


      
  // return T.get<0>()[0];
}
')

adlacode<-paste(adlacode,marshaller,sep="")
adlacode<-paste(adlacode,'\n',sep="")
adlacode<-paste(adlacode,'\n',sep="")
adlacode<-paste(adlacode,'\n',sep="")




recall.last<-function(f)
{
   m.f<-memoise(f)
   return(
   function(...)
   {
      args<-list(...)
      if(!do.call(has_cache(m.f),args)) forget(m.f)
      return(do.call(m.f,args))
   }
   )
}



if(output == "code")
{
   return(adlacode)
}

# compile the code
sourceCpp(code=adlacode)
# wrapper function for memoised marshaller
mem_f<-function(...,order=0,mem_marshal)
{
   args<-list(...)
   res<-do.call(mem_marshal,args)
   if(order == 0)
   {
     return(res[[1]])
   }
   if(order == 1)
   {
     return(res[[2]])
   }
   else
   {
     return(res[[3]])
   }
}


return(Curry(mem_f,mem_marshal=recall.last(eval(parse(text="marshal")))))

}








