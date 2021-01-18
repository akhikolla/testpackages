#include <Rcpp.h>

using namespace Rcpp;
using namespace std;


//algorithm code---------------
List CalcAll(NumericVector vnl, NumericVector vnr, int u,int w, NumericVector vremI,NumericVector vremJ,NumericMatrix A,NumericMatrix maxI,NumericMatrix maxJ,NumericMatrix r,double coef,CharacterVector nclass)
{
  			double maxV=-1;
        
        //cout << " CalcAll " << u<<","<<w << ".\n";
        //cout << " A(u,W) " <<A(u,w)<< ".\n";
				for(int i1=0; i1<vremI.size();i1++)
				{
					for(int j1=0; j1<vremJ.size(); j1++)
					{
					A(u,w)=A(u,vremI[i1]-1)+A(vremJ[j1]-1,w);
          //cout << " A(u,vrem) " <<A(u,vremI[i1]-1)<< ".\n";
          //cout << " A(vrem,W) " <<A(vremJ[j1]-1,w)<< ".\n";
					if(nclass[vremI[i1]-1]==nclass[vremJ[j1]-1])
					{
						A(u,w)=A(u,w)-pow(r(u,vremI[i1]-1),coef)-pow(r(w,vremJ[j1]-1),coef)+pow((r(u,vremI[i1]-1)+r(w,vremJ[j1]-1)),coef);
					}	
						
					if(A(u,w)>maxV)
					{
						maxV=A(u,w);
						maxI(u,w)=vremI[i1]-1;
						maxJ(u,w)=vremJ[j1]-1;
					}
					//if (all(nclass[vremJ]==nclass[vremJ[1]])) break;
					}
				//if (all(nclass[vremI]==nclass[vremI[1]])) break;
				}
				A(u,w)=maxV;
				A(w,u)=A(u,w);
				maxI(w,u)=maxJ(u,w);
				maxJ(w,u)=maxI(u,w);
        int ff=0;
        for(int k=0;k<vnr.size();k++)
        {
          if(nclass[vnr[k]-1]!=nclass[vnr[0]-1])
          {
            ff=1;
            break;
          }
        }
				if((ff==0)&&(nclass[maxI(u,w)]==nclass[maxJ(u,w)]))
				{
					//r[u,w]=length(node[[ind]]$right)+r[u,maxI[u,w]];
					r(u,w)=r(maxJ(u,w),w)+r(u,maxI(u,w));

				}
				else
				{
				r(u,w)=r(maxJ(u,w),w);
				}
        ff=0;
        for(int k=0;k<vnl.size();k++)
        {
          if(nclass[vnl[k]-1]!=nclass[vnl[0]-1])
          {
            ff=1;
            break;
          }
        }
				if((ff==0)&&(nclass[maxJ(w,u)]==nclass[maxI(w,u)]))
				{
					//r[w,u]=length(node[[ind]]$left)+r[w,maxI[w,u]];
					r(w,u)=r(maxJ(w,u),u)+r(w,maxI(w,u));

				}
				else
				{
				r(w,u)=r(maxJ(w,u),u);
				}
		
     /*A=as<NumericMatrix>(ttt["A"]);
     if((vremI.size()==1)&&(vremJ.size()==1))
     {
     cout << " u,w " << u<<","<<w << ".\n";
     cout << " out_1 " << A(u,w) << ".\n";
     
     }*/
   
	return List::create(_["A"]=A, _["maxI"]=maxI,  _["maxJ"]=maxJ,  _["r"]=r);
}



// [[Rcpp::export]]
List OrderingJosephC( int ind, NumericMatrix hc, List node, NumericMatrix A, NumericMatrix r, NumericMatrix maxI, NumericMatrix maxJ, CharacterVector nclass, double coef)
{

  List out;
  List res;
  List vhl;
   List vn;
    List vhr;
  NumericVector vv;
  NumericVector vnl;
  NumericVector vnr;
  NumericVector vl;
  NumericVector vr;
  NumericVector vvl;
  NumericVector vvr;
  int u;
  int w;
    int indL,indR;
  NumericVector vremI;
  NumericVector vremJ;
	//cat( "IND = ", ind,"\n");
	//cat("Elapsed time =",format(Sys.time()-ptm), "\n");
  /*if(ind==168)
  {
    for (int kk=0;kk<hc.nrow();kk++)
   
  cout <<  hc(kk,0)<<","<<hc(kk,1)<< ".\n";
  }*/
	//cout << " ind " << ind << ".\n";
	if(ind<0)
	{
		//A[-ind,-ind]=1;
		//int ttt=1;
		return List::create(_["A"]=A, _["maxI"]=maxI,  _["maxJ"]=maxJ,  _["r"]=r);
	}
	else
	{
    indL=hc(ind,0)-1;
    //cout << " indL " << indL << ".\n";
    res=OrderingJosephC(indL,hc, node, A, r, maxI, maxJ,nclass, coef);
    indR=hc(ind,1)-1;
    //cout << " indR " << indR << ".\n";
    
      A=as<NumericMatrix>(res["A"]);
      
		maxI=as<NumericMatrix>(res["maxI"]);
    
		maxJ=as<NumericMatrix>(res["maxJ"]);
    
		r=as<NumericMatrix>(res["r"]);
    
		res=OrderingJosephC(indR,hc, node,A, r, maxI, maxJ,nclass, coef);
		//if(ind==dim(hc$merge)[1]) browser();

    
		A=as<NumericMatrix>(res["A"]);
		maxI=as<NumericMatrix>(res["maxI"]);
		maxJ=as<NumericMatrix>(res["maxJ"]);
		r=as<NumericMatrix>(res["r"]);

    vn=node[ind];
    vnl=as<NumericVector>(vn[0]);
    vnr=as<NumericVector>(vn[1]);
    //cout << " ind to Calc "<<ind<<"\n"; 
    if(vnl.size()!=1)
    {
    vhl=node[hc(ind,0)-1];
    vl=as<NumericVector>(vhl[0]);
    vr=as<NumericVector>(vhl[1]);
    }
    if(vnr.size()!=1)
    {
    vhr=node[hc(ind,1)-1];
    vvl=as<NumericVector>(vhr[0]);
    vvr=as<NumericVector>(vhr[1]);
    }
  if(vnl.size()==1)
  {
    vremI=vnl;
    u=vnl[0];
    if(vnr.size()==1)
    {
      vremJ=vnr;
      w=vnr[0];
      
      out=CalcAll(vnl,vnr,u-1,w-1,vremI,vremJ,A,maxI,maxJ,r,coef,nclass);
      A=as<NumericMatrix>(out["A"]);
        maxI=as<NumericMatrix>(out["maxI"]);
	  	  maxJ=as<NumericMatrix>(out["maxJ"]);
		    r=as<NumericMatrix>(out["r"]);
      //A=as<NumericMatrix>(out["A"]);
       //cout << " u,w " << u-1<<","<<w-1 << ".\n";
       //cout << " vremI,vremJ " << vremI<<","<<vremJ << ".\n";
       //cout << "A(u,vremI) " << A(u-1,vremI[0]-1) << ".\n";
       //cout << "A(vremJ,w) " << A(vremJ[0]-1,w-1) << ".\n";
      //cout << " out_2 " << A(u-1,w-1) << ".\n";
      
    }
    else
    {
      vremJ=vvr;
      for(int j=0;j<vvl.size();j++)
      {
        w=vvl[j];
        
        out=CalcAll(vnl,vnr,u-1,w-1,vremI,vremJ,A,maxI,maxJ,r,coef,nclass);
        A=as<NumericMatrix>(out["A"]);
  	    maxI=as<NumericMatrix>(out["maxI"]);
	  	  maxJ=as<NumericMatrix>(out["maxJ"]);
		    r=as<NumericMatrix>(out["r"]);
        
      } 
      vremJ=vvl;
      for(int j=0;j<vvr.size();j++)
      {
        w=vvr[j];
        
        out=CalcAll(vnl,vnr,u-1,w-1,vremI,vremJ,A,maxI,maxJ,r,coef,nclass);
        A=as<NumericMatrix>(out["A"]);
        maxI=as<NumericMatrix>(out["maxI"]);
	  	  maxJ=as<NumericMatrix>(out["maxJ"]);
		    r=as<NumericMatrix>(out["r"]);
      }  
    }
  }
  else
  {
    vremI=vr;
		for(int i=0;i<vl.size();i++)
		{
			/*if(ind==hc.nrow())
			{
				cat( "I = ", i,"\n");
				cat("Elapsed time =",format(Sys.time()-ptm), "\n")
			}*/
			u=vl[i];
      if(vnr.size()==1)
      {
      vremJ=vnr;
      w=vnr[0];
      
      out=CalcAll(vnl,vnr,u-1,w-1,vremI,vremJ,A,maxI,maxJ,r,coef,nclass);
      A=as<NumericMatrix>(out["A"]);
  	  maxI=as<NumericMatrix>(out["maxI"]);
		  maxJ=as<NumericMatrix>(out["maxJ"]);
  		r=as<NumericMatrix>(out["r"]);
      }
      else
      {
			vremJ=vvr;
			for(int j=0;j<vvl.size();j++)
			{
				w=vvl[j];
		    
        out=CalcAll(vnl,vnr,u-1,w-1,vremI,vremJ,A,maxI,maxJ,r,coef,nclass);
        A=as<NumericMatrix>(out["A"]);
  	    maxI=as<NumericMatrix>(out["maxI"]);
		    maxJ=as<NumericMatrix>(out["maxJ"]);
		    r=as<NumericMatrix>(out["r"]);
			}
			vremJ=vvl;
      for(int j=0;j<vvr.size();j++)
      {
        w=vvr[j];

        out=CalcAll(vnl,vnr,u-1,w-1,vremI,vremJ,A,maxI,maxJ,r,coef,nclass);
        A=as<NumericMatrix>(out["A"]);
  	    maxI=as<NumericMatrix>(out["maxI"]);
		    maxJ=as<NumericMatrix>(out["maxJ"]);
		    r=as<NumericMatrix>(out["r"]);
      }
      }
		}
    vremI=vl;
    for(int i=0;i<vr.size();i++)
    {
      /*if(ind==hc.nrow())
      {
        cat( "I = ", i,"\n");
        cat("Elapsed time =",format(Sys.time()-ptm), "\n")
      }*/
      u=vr[i];
      if(vnr.size()==1)
      {
      vremJ=vnr;
      w=vnr[0];
      
      out=CalcAll(vnl,vnr,u-1,w-1,vremI,vremJ,A,maxI,maxJ,r,coef,nclass);
      A=as<NumericMatrix>(out["A"]);
  	  maxI=as<NumericMatrix>(out["maxI"]);
		  maxJ=as<NumericMatrix>(out["maxJ"]);
		  r=as<NumericMatrix>(out["r"]);
      }
      else
      {
      vremJ=vvr;
      for(int j=0;j<vvl.size();j++)
      {
        w=vvl[j];
        
        out=CalcAll(vnl,vnr,u-1,w-1,vremI,vremJ,A,maxI,maxJ,r,coef,nclass);
        A=as<NumericMatrix>(out["A"]);
  	    maxI=as<NumericMatrix>(out["maxI"]);
		    maxJ=as<NumericMatrix>(out["maxJ"]);
		    r=as<NumericMatrix>(out["r"]);
      }
      vremJ=vvl;
      for(int j=0;j<vvr.size();j++)
      {
        w=vvr[j];
        
        out=CalcAll(vnl,vnr,u-1,w-1,vremI,vremJ,A,maxI,maxJ,r,coef,nclass);
        A=as<NumericMatrix>(out["A"]);
  	    maxI=as<NumericMatrix>(out["maxI"]);
		    maxJ=as<NumericMatrix>(out["maxJ"]);
		    r=as<NumericMatrix>(out["r"]);
      }
      }
    }
  }
    
    return out;
}

}

