/**
 * A supporting class used in the definition of SortedRandomSets.
 *
 * History:
 * - 2018.01.01 file adapted from version 1.0 of the multinet library
 */


#ifndef UU_CORE_DATASTRUCTURES_CONTAINERS_DEFAULTEQ_H_
#define UU_CORE_DATASTRUCTURES_CONTAINERS_DEFAULTEQ_H_

namespace uu {
namespace core {


/**
 * SortedRandomSets allow to look for elements using generic comparison functions.
 * DefaultEQ() is the default way used by SortedRandomSets to check the "equals" predicate,
 * and uses the == operator.
 */
template <typename T1, typename T2>
struct DefaultEQ
{
    bool
    operator() (const T1& x, const T2& y) const
    {
        return x==y;
    }
};

}
}

#endif
