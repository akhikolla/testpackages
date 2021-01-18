#include "core/datastructures/observers/ObserverStore.hpp"

#include "core/exceptions/assert_not_null.hpp"

namespace uu {
namespace core {

void
ObserverStore::
register_observer(
    std::unique_ptr<core::GenericObserver> obs
)
{
    assert_not_null(obs.get(), "ObserverStore::register_observer", "obs");
    observers_.push_back(std::move(obs));
}

}
}
