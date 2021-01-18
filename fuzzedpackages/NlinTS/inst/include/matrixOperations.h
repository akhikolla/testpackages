/**
 * @authors Hmamouche Youssef
 * @date    09/06/2016
 **/

#ifndef MATRIXOPERATIONS_H
#define MATRIXOPERATIONS_H

#include "struct.h"
#include "exception.h"

using namespace Struct;

namespace MatrixOperations {
    bool regression (const CMatDouble &, const CVDouble &,  CVDouble &); // throw (Exception);

    void P_Part (CVDouble &  , CMatDouble & , CMatDouble & , unsigned int)
    ;

    void Pr_Part (CVDouble & , CMatDouble & , unsigned int );

    void Diff (CVDouble & );

    CVDouble VECbivar (CMatDouble  , unsigned , bool d /* = false */); // throw (Exception);
}


#endif // CMATRIXOPERATIONS_H
