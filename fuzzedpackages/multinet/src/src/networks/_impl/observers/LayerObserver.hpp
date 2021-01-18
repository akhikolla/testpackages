/**
 * History:
 * - 2018.03.09 file created, following a restructuring of the previous library.
 */

#ifndef UU_MNET_DATASTRUCTURE_OBSERVERS_LAYEROBSERVER_H_
#define UU_MNET_DATASTRUCTURE_OBSERVERS_LAYEROBSERVER_H_

#include "core/datastructures/observers/Observer.hpp"
#include "core/datastructures/observers/ObserverStore.hpp"
#include "networks/_impl/observers/LayerPropagateObserver.hpp"
#include "core/exceptions/assert_not_null.hpp"

namespace uu {
namespace net {

/**
 * This observer propagates a removal from one store to another.
 */
template<typename S, typename O>
class LayerObserver :
    public core::Observer<O>,
    public core::ObserverStore
{

  public:
    /**
     * Creates an observer with a pointer to a layer store to be notified when objects are erased.
     *
     */
    LayerObserver(
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
LayerObserver<S, O>::
LayerObserver(
    S* store
) :
    store_(store)
{
    core::assert_not_null(store_, "LayerObserver::constructor", "store");
}

template<typename S, typename O>
void
LayerObserver<S, O>::
notify_add(
    O* obj
)
{
    core::assert_not_null(obj, "LayerObserver::notify_add", "obj");

    store_->add(obj);

    auto obs = std::make_unique<LayerPropagateObserver<S, O, const Vertex>>(store_, obj);
    obj->vertices()->attach(obs.get());
    register_observer(std::move(obs));
}


template<typename S, typename O>
void
LayerObserver<S, O>::
notify_erase(
    O* obj
)
{
    core::assert_not_null(obj, "LayerObserver::notify_erase", "obj");

    store_->erase(obj);

    // @todo de-register observer

}


}
}

#endif
