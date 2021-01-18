#include "core/olap/selection/Set.hpp"
#include "core/exceptions/OutOfBoundsException.hpp"

namespace uu {
namespace core {
namespace sel {

Set::
Set(
    const std::vector<size_t>& indexes
) : indexes_(indexes)
{
    has_next_ = false;
}

/**  */

void
Set::
eval(
    size_t size
)
{
    for (auto idx: indexes_)
    {
        if (idx >=  size)
        {
            throw core::OutOfBoundsException("index too large");
        }
    }

    if (indexes_.size() > 0)
    {
        has_next_ = true;
        current_ = 0;
    }
}

/**  */
bool
Set::
has_next(
) const
{
    return has_next_;
}


/**  */
size_t
Set::
next(
)
{
    if (current_ < indexes_.size()-1)
    {
        return indexes_[current_++];
    }

    else
    {
        has_next_ = false;
        return indexes_[current_];
    }
}

}
}
}

