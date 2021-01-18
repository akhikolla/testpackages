/**
 * History:
 * - 2018.01.01 file adapted from version 1.0 of the multinet library
 */

#ifndef UU_CORE_ATTRIBUTES_ATTRIBUTEVALUEMAP_H_
#define UU_CORE_ATTRIBUTES_ATTRIBUTEVALUEMAP_H_

#include <vector>
#include "core/datastructures/containers/LabeledUniquePtrSortedRandomSet.hpp"
#include "core/exceptions/WrongFormatException.hpp"
#include "core/attributes/AttributeType.hpp"
#include "core/attributes/Attribute.hpp"
#include "core/attributes/Time.hpp"
#include "core/attributes/Text.hpp"
#include "core/attributes/Value.hpp"
#include "core/attributes/conversion.hpp"

namespace uu {
namespace core {

/**
 * A class associating multiple attributes and attribute values to a set of objects.
 */
template <typename ID>
class AttributeValueMap :
    public LabeledUniquePtrSortedRandomSet<const Attribute>
{
  private:

    typedef LabeledUniquePtrSortedRandomSet<const Attribute> super;

  public:

    using super::add;
    using super::erase;
    using super::get;

    const Attribute*
    add(
        std::string name,
        AttributeType t
    );

    /**
     * Sets the value of an attribute, given as a string.
     * @param id the id of the object whose associated value is set
     * @param attribute_name the name of the attribute
     * @param value the value to be set
     * @throws ElementNotFoundException if there is no attribute with this name
     * @throws WrongFormatException if value does not match attribute_name's type
     **/
    void
    set_as_string(
        ID oid,
        const std::string& attribute_name,
        const std::string& value
    );

    /*
     *
     * Gets the value of an attribute as a string.
     * @param id the id of the object whose associated value is set
     * @param attribute_name the name of the attribute
     * @return the value associated to the object id for the attribute, or a null value
     * @throws ElementNotFoundException if there is no attribute with this name
     **/
    Value<std::string>
    get_as_string(
        ID id,
        const std::string& // attribute_name
    ) const;


    /**
     * Creates an index on an attribute, so that queries based on values (min, max, range)
     * become faster at the expense of a higher time to set values.
     * @param attribute_name the name of the attribute to be indexed
     * @return true if the index has been created, false if it already existed
     * @throws ElementNotFoundException if there is no attribute with this name
     **/
    virtual
    bool
    add_index(
        const std::string& attribute_name
    ) = 0;


    /**
     * Removes an attribute values from an object.
     * @param id the id of the object whose associated value is set
     * @param attribute_name the name of the attribute
     * @return true if the object had an associated value
     * @throws ElementNotFoundException if there is no attribute with this name
     **/
    virtual
    bool
    reset(
        ID id,
        const std::string& attribute_name
    ) = 0;


    /**
     * Sets the value of an attribute.
     * @param id the id of the object whose associated value is set
     * @param attribute_name The name of the attribute
     * @param value The value to be set
     * @throws ElementNotFoundException if there is no STRING attribute with this name
     **/
    virtual
    void
    set_string(
        ID id,
        const std::string& attribute_name,
        const std::string& value
    ) = 0;


    /**
     * Sets the value of an attribute.
     * @param id the id of the object whose associated value is set
     * @param attribute_name The name of the attribute
     * @param value The value to be set
     * @throws ElementNotFoundException if there is no DOUBLE/NUMERIC attribute with this name
     **/
    virtual
    void
    set_double(
        ID id,
        const std::string& attribute_name,
        double value
    ) = 0;


    /**
     * Sets the value of an attribute.
     * @param id the id of the object whose associated value is set
     * @param attribute_name The name of the attribute
     * @param value The value to be set
     * @throws ElementNotFoundException if there is no INT attribute with this name
     **/
    virtual
    void
    set_int(
        ID id,
        const std::string& attribute_name,
        int value
    ) = 0;


    /**
     * Sets the value of a time attribute.
     * @param id the id of the object whose associated value is set
     * @param attribute_name The name of the attribute
     * @param value The value to be set
     * @throws ElementNotFoundException if there is no TIME attribute with this name
     **/
    virtual
    void
    set_time(
        ID id,
        const std::string& attribute_name,
        const Time& value
    ) = 0;

    /**
     * Sets the value of a text attribute.
     * @param id the id of the object whose associated value is set
     * @param attribute_name The name of the attribute
     * @param value The value to be set
     * @throws ElementNotFoundException if there is no TEXT attribute with this name
     **/
    virtual
    void
    set_text(
        ID oid,
        const std::string& attribute_name,
        const Text& value
    ) = 0;


    /**
     * Gets the value of an attribute.
     * @param id the id of the object whose associated value is retrieved
     * @param attribute_name The name of the attribute
     * @return The value associated to the object, or null if the object id has not been registered in this store
     * @throws ElementNotFoundException if there is no STRING attribute with this name
     **/
    virtual
    Value<std::string>
    get_string(
        ID id,
        const std::string& attribute_name
    ) const = 0;


    /**
     * Gets the value of an attribute.
     * @param id the id of the object whose associated value is retrieved
     * @param attribute_name The name of the attribute
     * @throws ElementNotFoundException if there is no DOUBLE attribute with this name
     **/
    virtual
    Value<double>
    get_double(
        ID id,
        const std::string& attribute_name
    ) const = 0;


    /**
     * Gets the value of an attribute.
     * @param id the id of the object whose associated value is retrieved
     * @param attribute_name The name of the attribute
     * @throws ElementNotFoundException if there is no INT attribute with this name
     **/
    virtual
    Value<int>
    get_int(
        ID id,
        const std::string& attribute_name
    ) const = 0;


    /**
     * Gets the value of a time attribute.
     * @param id the id of the object whose associated value is retrieved
     * @param attribute_name The name of the attribute
     * @throws ElementNotFoundException if there is no TIME attribute with this name
     **/
    virtual
    Value<Time>
    get_time(
        ID oid,
        const std::string& attribute_name
    ) const = 0;


    /**
     * Gets the value of a text attribute.
     * @param id the id of the object whose associated value is retrieved
     * @param attribute_name The name of the attribute
     * @throws ElementNotFoundException if there is no TEXT attribute with this name
     **/
    virtual
    Value<Text>
    get_text(
        ID oid,
        const std::string& attribute_name
    ) const = 0;


    /*****************************************************/
    /* RANGE QUERIES                                     */
    /*****************************************************/

    /**
     * Gets all the objects whose associated value is between the two input boundary values.
     * @param id the id of the object whose associated value is retrieved
     * @param attribute_name the name of the attribute
     * @param min_value lower bound of the range, inclusive
     * @param max_value upper bound of the range, inclusive
     * @param attribute_name the name of the attribute
     * @throws ElementNotFoundException if there is no string attribute with this name
     **/
    virtual
    std::vector<ID>
    range_query_string(
        const std::string& attribute_name,
        const std::string& min_value,
        const std::string& max_value
    ) const = 0;


    /**
     * Gets all the objects whose associated value is between the two input boundary values.
     * @param id the id of the object whose associated value is retrieved
     * @param attribute_name the name of the attribute
     * @param min_value lower bound of the range, inclusive
     * @param max_value upper bound of the range, inclusive
     * @param attribute_name the name of the attribute
     * @throws ElementNotFoundException if there is no int attribute with this name
     **/
    virtual
    std::vector<ID>
    range_query_int(
        const std::string& attribute_name,
        const int& min_value,
        const int& max_value
    ) const = 0;


    /**
     * Gets all the objects whose associated value is between the two input boundary values.
     * @param id the id of the object whose associated value is retrieved
     * @param attribute_name the name of the attribute
     * @param min_value lower bound of the range, inclusive
     * @param max_value upper bound of the range, inclusive
     * @param attribute_name the name of the attribute
     * @throws ElementNotFoundException if there is no double attribute with this name
     **/
    virtual
    std::vector<ID>
    range_query_double(
        const std::string& attribute_name,
        const double& min_value,
        const double& max_value
    ) const = 0;


    /**
     * Gets all the objects whose associated value is between the two input boundary values.
     * @param id the id of the object whose associated value is retrieved
     * @param attribute_name the name of the attribute
     * @param min_value lower bound of the range, inclusive
     * @param max_value upper bound of the range, inclusive
     * @param attribute_name the name of the attribute
     * @throws ElementNotFoundException if there is no time attribute with this name
     **/
    virtual
    std::vector<ID>
    range_query_time(
        const std::string& attribute_name,
        const Time& min_value,
        const Time& max_value
    ) const = 0;


    /*****************************************************/
    /* MIN QUERIES                                       */
    /*****************************************************/

    /**
     * Gets the minimum value for the input attribute.
     * @param attribute_name the name of the attribute
     * @throws ElementNotFoundException if there is no double attribute with this name
     **/
    virtual
    Value<double>
    get_min_double(
        const std::string& attribute_name
    ) const = 0;


    /**
     * Gets the minimum value for the input attribute.
     * @param attribute_name the name of the attribute
     * @throws ElementNotFoundException if there is no int attribute with this name
     **/
    virtual
    Value<int>
    get_min_int(
        const std::string& attribute_name
    ) const = 0;


    /**
     * Gets the minimum value for the input attribute.
     * @param attribute_name the name of the attribute
     * @throws ElementNotFoundException if there is no string attribute with this name
     **/
    virtual
    Value<std::string>
    get_min_string(
        const std::string& attribute_name
    ) const = 0;

    /**
     * Gets the minimum value for the input attribute.
     * @param attribute_name the name of the attribute
     * @throws ElementNotFoundException if there is no time attribute with this name
     **/
    virtual
    Value<Time>
    get_min_time(
        const std::string& attribute_name
    ) const = 0;


    /*****************************************************/
    /* MAX QUERIES                                       */
    /*****************************************************/

    /**
     * Gets the maximum value for the input attribute.
     * @param attribute_name the name of the attribute
     * @throws ElementNotFoundException if there is no int attribute with this name
     **/
    virtual
    Value<int>
    get_max_int(
        const std::string& attribute_name
    ) const = 0;


    /**
     * Gets the maximum value for the input attribute.
     * @param attribute_name the name of the attribute
     * @throws ElementNotFoundException if there is no double attribute with this name
     **/
    virtual
    Value<double>
    get_max_double(
        const std::string& attribute_name
    ) const = 0;


    /**
     * Gets the maximum value for the input attribute.
     * @param attribute_name the name of the attribute
     * @throws ElementNotFoundException if there is no string attribute with this name
     **/
    virtual
    Value<std::string>
    get_max_string(
        const std::string& attribute_name
    ) const = 0;


    /**
     * Gets the maximum value for the input attribute.
     * @param attribute_name the name of the attribute
     * @throws ElementNotFoundException if there is no time attribute with this name
     **/
    virtual
    Value<Time>
    get_max_time(
        const std::string& attribute_name
    ) const = 0;

} ;



// Type-specific get/set queries

template <typename ID>
void
AttributeValueMap<ID>::
set_as_string(
    ID oid,
    const std::string& attribute_name,
    const std::string& value
)
{

    const Attribute* att = get(attribute_name);

    if (!att)
    {
        throw ElementNotFoundException("attribute " + attribute_name);
    }

    switch (att->type)
    {
    case AttributeType::NUMERIC:
    case AttributeType::DOUBLE:
        set_double(oid,attribute_name,core::to_double(value));
        break;

    case AttributeType::INTEGER:
        set_int(oid,attribute_name,core::to_int(value));
        break;

    case AttributeType::STRING:
        set_string(oid,attribute_name,value);
        break;

    case AttributeType::TIME:
        set_time(oid,attribute_name, core::to_time(value));
        break;

    case AttributeType::TEXT:
        set_text(oid, attribute_name, core::to_text(value));
        break;
    }
}


template <typename ID>
const Attribute*
AttributeValueMap<ID>::
add(
    std::string name,
    AttributeType t
)
{
    std::unique_ptr<Attribute> a = std::make_unique<Attribute>(name, t);
    const Attribute* a_ptr = a.get();
    add(std::move(a));
    return a_ptr;
}

template <typename ID>
Value<std::string>
AttributeValueMap<ID>::
get_as_string(
    ID oid,
    const std::string& attribute_name
) const
{
    const Attribute* att = get(attribute_name);

    if (!att)
    {
        throw ElementNotFoundException("attribute " + attribute_name);
    }

    switch (att->type)
    {

    case AttributeType::NUMERIC:
    case AttributeType::DOUBLE:
        return core::to_string(get_double(oid,att->name));

    case AttributeType::INTEGER:
        return core::to_string(get_int(oid,att->name));

    case AttributeType::STRING:
        return get_string(oid,att->name);

    case AttributeType::TIME:
        return core::to_string(get_time(oid, att->name));

    case AttributeType::TEXT:
        return core::to_string(get_text(oid, att->name));
    }
}


}
}

#endif
