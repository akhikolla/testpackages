/**
 * History:
 * - 2018.03.09 file created, following a restructuring of the previous library.
 */

#ifndef UU_CORE_STORES_OBJECTSTORE_H_
#define UU_CORE_STORES_OBJECTSTORE_H_

#include <memory>
#include <map>
#include "core/datastructures/containers/SharedPtrSortedRandomSet.hpp"
#include "core/datastructures/observers/Observer.hpp"
#include "core/datastructures/observers/Subject.hpp"

namespace uu {
namespace core {

/**
 * A ObjectStore allows to create, store, retrieve and erase a set of elements.
 *
 * OBJECT_TYPE must have:
 * a typedef key_type, to specify the input to create the objects used as keys by the store
 */
template <typename OBJECT_TYPE>
class
    ObjectStore :
    public core::SharedPtrSortedRandomSet<const OBJECT_TYPE>,
    public core::Subject<const OBJECT_TYPE>
{

  private:

    typedef core::SharedPtrSortedRandomSet<const OBJECT_TYPE> super;

  public:

    typedef OBJECT_TYPE value_type;
    typedef typename OBJECT_TYPE::key_type key_type;


    ObjectStore(
    );

    using super::add;
    using super::size;
    using core::Subject<const OBJECT_TYPE>::observers;

    const OBJECT_TYPE *
    add(
        std::shared_ptr<const OBJECT_TYPE> v
    ) override;

    /** Creates a new object and adds it to the store.
    const OBJECT_TYPE *
    add(
        const typename OBJECT_TYPE::key_type& key
    );
     */

    const OBJECT_TYPE *
    get(
        const typename OBJECT_TYPE::key_type& key
    ) const;

    bool
    erase(
        const OBJECT_TYPE * v
    ) override;

    /**
     * Returns a short string summary of this store, for example including
     * the number of objects it contains.
     */
    virtual
    std::string
    summary(
    ) const;


  protected:

    /** Index: find element by key. */
    std::map<typename OBJECT_TYPE::key_type, const OBJECT_TYPE*> cidx_element_by_name;

};


template <typename OBJECT_TYPE>
ObjectStore<OBJECT_TYPE>::
ObjectStore(
)
{
}


template <typename OBJECT_TYPE>
const OBJECT_TYPE *
ObjectStore<OBJECT_TYPE>::
add(
    std::shared_ptr<const OBJECT_TYPE> v
)
{

    core::assert_not_null(v.get(), "add", "v");

    // Notify the observers.
    for (auto obs: observers)
    {
        obs->notify_add(v.get());
    }

    // Return false if a vertex with this key exists.
    auto search = cidx_element_by_name.find(v->key);

    if (search != cidx_element_by_name.end())
    {
        return nullptr;
    }

    const OBJECT_TYPE* res = super::add(v);

    // Indexing.
    cidx_element_by_name[v->key] = v.get();

    return res;
}


/*
template <typename OBJECT_TYPE>
const OBJECT_TYPE *
ObjectStore<OBJECT_TYPE>::
add(
    const typename OBJECT_TYPE::key_type& key
)
{
    if (!get(key))
    {
        return add(OBJECT_TYPE::create(key));
    }

    else
    {
        return nullptr;
    }
}
*/

template <typename OBJECT_TYPE>
const OBJECT_TYPE *
ObjectStore<OBJECT_TYPE>::
get(
    const typename OBJECT_TYPE::key_type& key
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

template <typename OBJECT_TYPE>
bool
ObjectStore<OBJECT_TYPE>::
erase(
    const OBJECT_TYPE * v
)
{
    core::assert_not_null(v, "erase", "v");


    // Notify the observers.
    for (auto obs: observers)
    {
        obs->notify_erase(v);
    }

    auto search = cidx_element_by_name.find(v->key);

    if (search != cidx_element_by_name.end())
    {
        cidx_element_by_name.erase(search);
        super::erase(v);
        return true;
    }

    else
    {
        return false;
    }
}


template <typename OBJECT_TYPE>
std::string
ObjectStore<OBJECT_TYPE>::
summary(
) const
{
    size_t s = size();
    std::string summary =
        std::to_string(s) +
        (s==1?" object":" objects");
    return summary;
}
}
}

#endif
