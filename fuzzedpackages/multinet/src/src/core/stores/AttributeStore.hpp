/**
 * History:
 * - 2018.03.09 file created, with code taken from other existing files.
 */

#ifndef UU_CORE_STORES_ATTRIBUTESTORE_H_
#define UU_CORE_STORES_ATTRIBUTESTORE_H_

#include <memory>
#include "core/attributes/MainMemoryAttributeValueMap.hpp"
#include "core/datastructures/observers/Observer.hpp"
#include "core/exceptions/NullPtrException.hpp"

namespace uu {
namespace core {

/**
 * A class to handle attributes to be associated to objects of type OT.
 */
template <typename OT>
class AttributeStore :
    public MainMemoryAttributeValueMap<const OT*>,
    public Observer<const OT>
{

  public:

    /**
     * Constructor.
     */
    AttributeStore();

    /**
     * Returns a short summary of the store, indicating the number of attributes.
     */
    virtual
    std::string
    summary(
    ) const;

    /**
     * When an object no longer exists, this method must be called to erase its values from the store.
     */
    void
    notify_erase(
        const OT* const o
    );


    void
    notify_add(
        const OT* const o
    );

    /**
     *
     */
    virtual
    void
    read_attributes(
        const OT* v,
        const std::vector<std::string>& fields,
        size_t offset,
        const std::vector<Attribute>& attributes,
        size_t line_number);

};


template <typename OT>
std::unique_ptr<AttributeStore<OT>>
                                 create_attribute_store()
{
    return std::make_unique<AttributeStore<OT>>();
}

template <typename OT>
AttributeStore<OT>::
AttributeStore()
{
}

template <typename OT>
std::string
AttributeStore<OT>::
summary(
) const
{
    size_t s = this->size();
    std::string summary = std::to_string(s) + ((s==1)?" attribute":" attributes");
    return summary;
}


template <typename OT>
void
AttributeStore<OT>::
notify_erase(
    const OT* object
)
{
    if (!object)
    {
        throw NullPtrException("AttributeStore::notify_erase()");
    }

    for (auto att: *this)
    {
        this->reset(object, att->name);
    }
}

template <typename OT>
void
AttributeStore<OT>::
notify_add(
    const OT* object
)
{
    if (!object)
    {
        throw NullPtrException("AttributeStore::notify_add()");
    }
}


template <typename OT>
void
AttributeStore<OT>::
read_attributes(
    const OT* v,
    const std::vector<std::string>& fields,
    size_t offset,
    const std::vector<Attribute>& attributes,
    size_t line_number
)
{

    int idx = offset;

    if (offset+attributes.size()>fields.size())
        throw WrongFormatException("Line " +
                                   std::to_string(line_number) +
                                   ": not enough attribute values");

    for (Attribute attribute: attributes)
    {
        this->set_as_string(v, attribute.name, fields[idx]);
        idx++;
    }
}


}
}

#endif
