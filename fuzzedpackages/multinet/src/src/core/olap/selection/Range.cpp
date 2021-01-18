#include "core/olap/selection/Range.hpp"
#include "core/exceptions/OutOfBoundsException.hpp"

namespace uu {
namespace core {
namespace sel {

Range::
Range(
    size_t from,
    size_t to
)
{
    from_ = from;
    to_ = to;
    has_next_ = false;
}

/**  */

void
Range::
eval(
    size_t size
)
{
    if (from_ >= size || to_ >= size)
    {
        throw core::OutOfBoundsException("range index too large");
    }

    has_next_ = true;
    current_ = from_;
}

/**  */
bool
Range::
has_next(
) const
{
    return has_next_;
}


/**  */
size_t
Range::
next(
)
{
    if (from_ < to_)
    {
        if (current_ == to_)
        {
            has_next_ = false;
            return current_;
        }

        else
        {
            return current_++;
        }
    }

    else
    {
        if (current_ == to_)
        {
            has_next_ = false;
            return current_;
        }

        else
        {
            return current_--;
        }
    }
}

}
}
}

