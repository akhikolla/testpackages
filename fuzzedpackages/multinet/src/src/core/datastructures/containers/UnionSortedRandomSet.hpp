/**
 * History:
 * - 2019.08.01 file created, following a restructuring of the previous library.
 */

#ifndef UU_CORE_DATASTRUCTURES_CONTAINERS_UNIONSORTEDRANDOMSET_H_
#define UU_CORE_DATASTRUCTURES_CONTAINERS_UNIONSORTEDRANDOMSET_H_

#include <map>
#include "core/datastructures/containers/SortedRandomBag.hpp"
#include "core/datastructures/observers/Observer.hpp"
#include "core/datastructures/observers/Subject.hpp"
#include "core/exceptions/assert_not_null.hpp"

namespace uu {
namespace core {

/**
 * A UnionSortedRandomSet contains all the elements present in a set of related random sets.
 *
 * ELEMENT_TYPE must have a typedef key_type, to specify the input to create the objects used as
 * keys by the store
 */
template <typename ELEMENT_TYPE>
class
    UnionSortedRandomSet :
    public SortedRandomBag<const ELEMENT_TYPE*>,
    public Subject<const ELEMENT_TYPE>,
    public Observer<const ELEMENT_TYPE>
{

  private:

    typedef core::SortedRandomBag<const ELEMENT_TYPE*> super;

  public:
    UnionSortedRandomSet(
    );

    using super::size;
    using core::Subject<const ELEMENT_TYPE>::observers;

    const ELEMENT_TYPE *
    get(
        const typename ELEMENT_TYPE::key_type& key
    ) const;

    /**
     * Informs the graph that an element has been added to a store.
     */
    void
    notify_add(
        const ELEMENT_TYPE* el
    ) override;


    /**
     * Informs the graph that an element has been erased from a store.
     */
    void
    notify_erase(
        const ELEMENT_TYPE* el
    ) override;

    /**
     * Returns a short string summary of this store, for example including
     * the number of elements it contains.
     */
    virtual
    std::string
    summary(
    ) const;

  private:

    bool
    add(
        const ELEMENT_TYPE* v
    );

    bool
    erase(
        const ELEMENT_TYPE * v
    );

  protected:

    /** Index: find element by key. */
    std::map<typename ELEMENT_TYPE::key_type, const ELEMENT_TYPE*> cidx_element_by_name;

};


template <typename ELEMENT_TYPE>
UnionSortedRandomSet<ELEMENT_TYPE>::
UnionSortedRandomSet(
)
{
}


template <typename ELEMENT_TYPE>
bool
UnionSortedRandomSet<ELEMENT_TYPE>::
add(
    const ELEMENT_TYPE* e
)
{

    core::assert_not_null(e, "UnionSortedRandomSet::add", "e");

    // Notify the observers.
    for (auto obs: observers)
    {
        obs->notify_add(e);
    }

    bool res = super::add(e);

    // Indexing @todo this assumes that this is the same element as before, if it was already present
    // with the same key
    cidx_element_by_name[e->key] = e;

    return res;
}


template <typename ELEMENT_TYPE>
const ELEMENT_TYPE *
UnionSortedRandomSet<ELEMENT_TYPE>::
get(
    const typename ELEMENT_TYPE::key_type& key
) const
{
    auto search = cidx_element_by_name.find(key);

    if (search != cidx_element_by_name.end())
    {
        return search->second;
    }

    else
    {
        return nullptr;
    }
}

template <typename ELEMENT_TYPE>
bool
UnionSortedRandomSet<ELEMENT_TYPE>::
erase(
    const ELEMENT_TYPE * e
)
{
    core::assert_not_null(e, "UnionSortedRandomSet:erase", "e");


    // Notify the observers.
    for (auto obs: observers)
    {
        obs->notify_erase(e);
    }

    auto search = cidx_element_by_name.find(e->key);

    if (search != cidx_element_by_name.end())
    {
        cidx_element_by_name.erase(search);
        super::erase(e);
        return true;
    }

    else
    {
        return false;
    }
}


template <typename ELEMENT_TYPE>
void
UnionSortedRandomSet<ELEMENT_TYPE>::
notify_add(
    const ELEMENT_TYPE* el
)
{
    add(el);
}


template <typename ELEMENT_TYPE>
void
UnionSortedRandomSet<ELEMENT_TYPE>::
notify_erase(
    const ELEMENT_TYPE* el
)
{
    erase(el);
}


template <typename ELEMENT_TYPE>
std::string
UnionSortedRandomSet<ELEMENT_TYPE>::
summary(
) const
{
    size_t s = size();
    std::string summary =
        std::to_string(s) +
        (s==1?" element":" elements");
    return summary;
}
}
}

#endif
