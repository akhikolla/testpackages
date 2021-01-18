/**
 * History:
 * - 2018.03.09 file created, following a restructuring of the previous library.
 */

#ifndef UU_NET_DATASTRUCTURE_OBSERVERS_PROPAGATEOBSERVER_H_
#define UU_NET_DATASTRUCTURE_OBSERVERS_PROPAGATEOBSERVER_H_

#include "core/datastructures/observers/Observer.hpp"
#include "core/exceptions/NullPtrException.hpp"

namespace uu {
namespace net {

/**
 * This observer propagates a removal from one store to another.
 */
template<typename S, typename O>
class PropagateObserver :
    public core::Observer<O>
{

  public:
    /**
     * Creates an observer with a pointer to the store to be notified when objects are erased.
     *
     */
    PropagateObserver(
        S* store
    );

    /**
     * Informs the observer that an object has been added.
     */
    void
    notify_add(
        O* v
    ) override;


    /**
     * Informs the observer that an object has been erased.
     */
    void
    notify_erase(
        O* v
    ) override;

  private:
    /** Internal object store. */
    S* store_;

};



template<typename S, typename O>
PropagateObserver<S, O>::
PropagateObserver(
    S* store
) :
    store_(store)
{
    if (!store_)
    {
        throw core::NullPtrException("store to be registered in the observer");
    }
}

template<typename S, typename O>
void
PropagateObserver<S, O>::
notify_add(
    O* obj
)
{
    if (!obj)
    {
        throw  core::NullPtrException("object passed to the observer");
    }
}


template<typename S, typename O>
void
PropagateObserver<S, O>::
notify_erase(
    O* obj
)
{
    if (!obj)
    {
        throw  core::NullPtrException("object passed to the observer");
    }

    store_->erase(obj);
}


}
}

#endif
