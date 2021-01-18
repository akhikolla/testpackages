/**
 * History:
 * - 2018.03.09 file created, following a restructuring of the previous library.
 */

#ifndef UU_MNET_DATASTRUCTURE_OBSERVERS_LAYERPROPAGATEOBSERVER_H_
#define UU_MNET_DATASTRUCTURE_OBSERVERS_LAYERPROPAGATEOBSERVER_H_

#include "core/datastructures/observers/Observer.hpp"
#include "core/exceptions/assert_not_null.hpp"

namespace uu {
namespace net {

/**
 * This observer propagates a removal from one store to another.
 */
template<typename S, typename L, typename O>
class LayerPropagateObserver :
    public core::Observer<O>
{

  public:
    /**
     * Creates an observer with a pointer to a layer store to be notified when objects are erased.
     *
     */
    LayerPropagateObserver(
        S* store,
        const L* layer
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
    /** Layer. */
    const L* layer_;

};



template<typename S, typename L, typename O>
LayerPropagateObserver<S, L, O>::
LayerPropagateObserver(
    S* store,
    const L* layer
) :
    store_(store), layer_(layer)
{
    core::assert_not_null(store_, "LayerPropagateObserver::constructor", "store");
    core::assert_not_null(layer_, "LayerPropagateObserver::constructor", "layer");
}

template<typename S, typename L, typename O>
void
LayerPropagateObserver<S, L, O>::
notify_add(
    O* obj
)
{
    core::assert_not_null(obj, "LayerPropagateObserver::notify_add", "obj");

}


template<typename S, typename L, typename O>
void
LayerPropagateObserver<S, L, O>::
notify_erase(
    O* obj
)
{
    core::assert_not_null(obj, "LayerPropagateObserver::notify_erase", "obj");

    // Uncomment to print a summary of the network

    store_->erase(layer_, obj);
}


}
}

#endif
