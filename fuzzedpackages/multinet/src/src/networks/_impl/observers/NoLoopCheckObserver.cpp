/**
 * History:
 * - 2018.03.09 file created, following a restructuring of the previous library.
 */


#include "networks/_impl/observers/NoLoopCheckObserver.hpp"

namespace uu {
namespace net {

NoLoopCheckObserver::
NoLoopCheckObserver(
)
{
}


void
NoLoopCheckObserver::
notify_add(
    const Edge* e
)
{
    if (!e)
    {
        throw core::NullPtrException("edge passed to the observer");
    }

    if (e->v1 == e->v2)
    {
        throw core::WrongParameterException("loops are not allowed");
    }
}


void
NoLoopCheckObserver::
notify_erase(
    const Edge* e
)
{
    if (!e)
    {
        throw core::NullPtrException("edge passed to the observer");
    }
}


}
}

