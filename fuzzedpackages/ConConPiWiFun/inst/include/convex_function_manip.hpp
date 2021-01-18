/*
 * convex_function_manip.hpp
 *
 *  Created on: 16 avr. 2013
 *      Author: robin
 */

#ifndef CONVEX_FUNCTION_MANIP_HPP_
#define CONVEX_FUNCTION_MANIP_HPP_

cplfunction Suml(cplfunction const & cplfunction_1,cplfunction const & cplfunction_2)
{
	 cplfunction tmp1(cplfunction_1);cplfunction  tmp2(cplfunction_2);
   		if (cplfunction_2.Breakpoints_.size()>cplfunction_1.Breakpoints_.size()){
   			tmp2.Sumf(tmp1);
   			return(tmp2);
   		}else{
   			tmp1.Sumf(tmp2);
   			return(tmp1);
   		}
}


cplfunction InfConv(cplfunction const & cplFunction_1,cplfunction const & cplFunction_2){
       cplfunction tmp1=cplFunction_1,tmp2=cplFunction_2;
    	 tmp1.Etoile();
    	 tmp2.Etoile();
    	 cplfunction res=Suml(tmp1,tmp2);
    	 res.Etoile();
    	 return(res);
     }
// static void finalizer_of_cplfunction( cplfunction* ptr ){
//    ptr->cplfunction::~cplfunction();
    //printf("finalizer has been called\n");
// }

cplfunction InfConfFunct(cplfunction const & cplFunction_1,cplfunction const & cplFunction_2,double y ){
       cplfunction tmp1(cplFunction_1),tmp2(cplFunction_2);
    	 tmp2.Swap(y);
    	 cplfunction B=Suml(tmp1,tmp2);
    	 return(B);
     }



cpqfunction Sumq(cpqfunction const & cpqfunction_1,cpqfunction const & cpqfunction_2){
   cpqfunction tmp1(cpqfunction_1),tmp2(cpqfunction_2);
   		if (cpqfunction_2.Breakpoints_.size()>cpqfunction_1.Breakpoints_.size()){
   			tmp2.Sumf(tmp1);
   			return(tmp2);
   		}else{
   			tmp1.Sumf(tmp2);
   			return(tmp1);
   		}
}


cpqfunction InfConvq(cpqfunction const & cpqfunction_1,cpqfunction const & cpqfunction_2){
       cpqfunction tmp1(cpqfunction_1),tmp2(cpqfunction_2);
    	 tmp1.Etoile();
    	 tmp2.Etoile();
    	 cpqfunction res=Sumq(tmp1,tmp2);
    	 res.Etoile();
    	 return(res);
     }
// static void finalizer_of_cpqfunction( cpqfunction* ptr ){
//    ptr->cpqfunction::~cpqfunction();
    //printf("finalizer has been called\n");
// }

cpqfunction InfConfFunctq(cpqfunction const & cpqfunction_1,cpqfunction const & cpqfunction_2,double y ){
       cpqfunction tmp1(cpqfunction_1),tmp2(cpqfunction_2);
    	 tmp2.Swap(y);
    	 cpqfunction B=Sumq(tmp1,tmp2);
    	 return(B);
     };


#endif /* CONVEX_FUNCTION_MANIP_HPP_ */
