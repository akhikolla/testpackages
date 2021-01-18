/**
 * History:
 * - 2018.03.09 file created, following a restructuring of the previous library.
 */

#ifndef UU_CORE_DATASTRUCTURES_CONTAINERS_LABELEDUNIQUEPTRSORTEDRANDOMSET_H_
#define UU_CORE_DATASTRUCTURES_CONTAINERS_LABELEDUNIQUEPTRSORTEDRANDOMSET_H_

#include <memory>
#include <string>
#include <unordered_map>
#include "core/datastructures/containers/UniquePtrSortedRandomSet.hpp"
#include "core/exceptions/assert_not_null.hpp"

namespace uu {
namespace core {

/**
 * A container whose elements are unique pointers to objects having a 'name' attribute.
 */
template <typename E>
class
    LabeledUniquePtrSortedRandomSet :
    public core::UniquePtrSortedRandomSet<E>
{
  private:
    typedef core::UniquePtrSortedRandomSet<E> super;

  public:

    LabeledUniquePtrSortedRandomSet(
    );

    using super::add;

    /**
     * Adds a pointer to an element to the store.
     *
     * @param element a pointer to an existing element.
     * @return false if an element with the same name already exists in the set, true otherwise.
     **/
    virtual
    E *
    add(
        std::unique_ptr<E> e
    ) override;

    /**
     * Erases an existing element from the store.
     *
     * @param element a pointer to the element to be deleted
     * @return true if the object has been erased, false if it was not present.
     **/
    virtual
    bool
    erase(
        E * const element
    ) override;

    /**
     * Returns an element by name.
     *
     * This function can also be used to check if an element exists.
     *
     * @param name name of the element.
     * @return a pointer to the requested element, or nullptr if it does not exist.
     **/
    E *
    get(
        const std::string& name
    ) const;

  protected:

    /** Index: find element by name. */
    std::unordered_map<std::string, E *> cidx_element_by_name;

};

template <typename E>
LabeledUniquePtrSortedRandomSet<E>::
LabeledUniquePtrSortedRandomSet(
)
{
}


template <typename E>
E *
LabeledUniquePtrSortedRandomSet<E>::
add(
    std::unique_ptr<E> element
)
{
    core::assert_not_null(element.get(), "add", "element");

    // Return false if an element with this name exists.
    auto search = cidx_element_by_name.find(element->name);

    if (search != cidx_element_by_name.end())
    {
        return nullptr;
    }

    // Indexing.
    cidx_element_by_name[element->name] = element.get();

    E* res = super::add(std::move(element));

    return res;
}

template <typename E>
E *
LabeledUniquePtrSortedRandomSet<E>::
get(
    const std::string& name
) const
{
    auto search = cidx_element_by_name.find(name);

    if (search != cidx_element_by_name.end())
    {
        return search->second;
    }

    else
    {
        return nullptr;
    }
}

template <typename E>
bool
LabeledUniquePtrSortedRandomSet<E>::
erase(
    E * const element
)
{

    core::assert_not_null(element, "erase", "element");

    auto search = cidx_element_by_name.find(element->name);

    if (search != cidx_element_by_name.end())
    {
        cidx_element_by_name.erase(search);
        super::erase(element);
        return true;
    }

    else
    {
        return false;
    }
}

}
}

#endif
