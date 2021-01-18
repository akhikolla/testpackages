/**
 * History:
 * - 2018.01.01 file adapted from version 1.0 of the multinet library
 */

#ifndef UU_CORE_DATASTRUCTURES_PROPERTYMATRIX_STRUCTURECOMPARISONFUNCTION_H_
#define UU_CORE_DATASTRUCTURES_PROPERTYMATRIX_STRUCTURECOMPARISONFUNCTION_H_


#include "core/datastructures/propertymatrix/PropertyMatrix.hpp"


namespace uu {
namespace core {

/**
 * A function used to compare two structures in a property matrix.
 */
template <class STRUCTURE, class CONTEXT, class NUMBER>
class
    StructureComparisonFunction
{

  public:

    StructureComparisonFunction(
        const PropertyMatrix<STRUCTURE,CONTEXT,NUMBER>* P,
        const CONTEXT* c
    ) :
        P(P),
        c(c)
    {

    }

    const PropertyMatrix<STRUCTURE,CONTEXT,NUMBER>* P;

    const CONTEXT* c;

    bool
    operator()(
        const STRUCTURE& s1,
        const STRUCTURE& s2
    ) const
    {
        Value<NUMBER> v1 = (*P).get(s1,*c);
        Value<NUMBER> v2 = (*P).get(s2,*c);

        if (v1.null || v2.null)
        {
            return !v1.null<!v2.null;    // TODO
        }

        return v1.value<v2.value;
    }
};

}
}

#endif
