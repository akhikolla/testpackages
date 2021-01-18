/**
 * History:
 * - 2018.03.09 file created, following a restructuring of the previous library.
 */

#ifndef UU_CORE_DATASTRUCTURES_CONTAINERS_UNIQUEPTREQ_H_
#define UU_CORE_DATASTRUCTURES_CONTAINERS_UNIQUEPTREQ_H_

#include <memory>

namespace uu {
namespace core {


/**
 * An object used to look for unique_ptr's inside this set using their raw pointer as a key.
 */
template<typename T>
struct UniquePtrEQ
{
    bool
    operator() (const std::unique_ptr<T>& x, const T* const & y) const
    {
        return x.get()==y;
    }
    bool
    operator() (std::unique_ptr<T>& x, T*& y) const
    {
        return x.get()==y;
    }
};


}
}

#endif

