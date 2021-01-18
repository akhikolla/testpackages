/**
 * History:
 * - 2018.03.09 file created, following a restructuring of the previous library.
 */

#ifndef UU_NET_DATASTRUCTURE_OBSERVERS_PROPAGATEADDERASEOBSERVER_H_
#define UU_NET_DATASTRUCTURE_OBSERVERS_PROPAGATEADDERASEOBSERVER_H_

#include "core/datastructures/observers/Observer.hpp"
#include "core/exceptions/assert_not_null.hpp"

namespace uu {
namespace net {

/**
 * This observer propagates a removal from one store to another.
 */
template<typename S, typename O>
class PropagateAddEraseObserver :
    public core::Observer<O>
{

  public:
    /**
     * Creates an observer with a pointer to the store to be notified when objects are erased.
     *
     */
    PropagateAddEraseObserver(
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
PropagateAddEraseObserver<S, O>::
PropagateAddEraseObserver(
    S* store
) :
    store_(store)
{
    assert_not_null(store_, "PropagateAddEraseObserver::constructor", "store");
}

template<typename S, typename O>
void
PropagateAddEraseObserver<S, O>::
notify_add(
    O* obj
)
{
    assert_not_null(obj, "PropagateAddEraseObserver::notify_add", "obj");

    store_->add(obj);
}


template<typename S, typename O>
void
PropagateAddEraseObserver<S, O>::
notify_erase(
    O* obj
)
{

    assert_not_null(obj, "PropagateAddEraseObserver::notify_erase", "obj");

    store_->erase(obj);
}


}
}

#endif
