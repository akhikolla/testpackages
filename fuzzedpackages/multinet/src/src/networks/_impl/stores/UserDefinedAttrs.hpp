/**
 * History:
 * - 2018.03.19 file created, with code taken from other existing files.
 */

#ifndef UU_WNET_DATASTRUCTURES_STORES_USERDEFINEDATTRS_H_
#define UU_WNET_DATASTRUCTURES_STORES_USERDEFINEDATTRS_H_

#include <memory>
#include "core/stores/AttributeStore.hpp"
#include "core/exceptions/ElementNotFoundException.hpp"

namespace uu {
namespace net {

template <typename OT>
class UserDefinedAttrs
{

  protected:

    /**
     * Constructor.
     */
    UserDefinedAttrs(core::AttributeStore<OT>* attr);

  public:

    const core::Attribute*
    add(
        std::string name,
        core::AttributeType t
    );

    const core::Attribute*
    add(
        std::unique_ptr<const core::Attribute> attribute
    );

    const core::Attribute*
    get(
        std::string name
    ) const;

    void
    set_as_string(
        const OT* oid,
        const std::string& attribute_name,
        const std::string& value
    );


    core::Value<std::string>
    get_as_string(
        const OT* oid,
        const std::string& // attribute_name
    ) const;


    bool
    add_index(
        const std::string& attribute_name
    );


    bool
    reset(
        const OT* oid,
        const std::string& attribute_name
    );



    void
    set_string(
        const OT* oid,
        const std::string& attribute_name,
        const std::string& value
    );


    void
    set_double(
        const OT* oid,
        const std::string& attribute_name,
        double value
    );




    void
    set_int(
        const OT* oid,
        const std::string& attribute_name,
        int value
    );



    void
    set_time(
        const OT* oid,
        const std::string& attribute_name,
        const core::Time& value
    );


    void
    set_text(
        const OT* oid,
        const std::string& attribute_name,
        const core::Text& value
    );



    core::Value<std::string>
    get_string(
        const OT* oid,
        const std::string& attribute_name
    ) const;


    core::Value<double>
    get_double(
        const OT* oid,
        const std::string& attribute_name
    ) const;


    core::Value<int>
    get_int(
        const OT* oid,
        const std::string& attribute_name
    ) const;


    core::Value<core::Time>
    get_time(
        const OT* oid,
        const std::string& attribute_name
    ) const;


    core::Value<core::Text>
    get_text(
        const OT* oid,
        const std::string& attribute_name
    ) const;


    /*****************************************************/
    /* RANGE QUERIES                                     */
    /*****************************************************/


    std::vector<OT>
    range_query_string(
        const std::string& attribute_name,
        const std::string& min_value,
        const std::string& max_value
    ) const;


    std::vector<OT>
    range_query_int(
        const std::string& attribute_name,
        const int& min_value,
        const int& max_value
    ) const;


    std::vector<OT>
    range_query_double(
        const std::string& attribute_name,
        const double& min_value,
        const double& max_value
    ) const;



    std::vector<OT>
    range_query_time(
        const std::string& attribute_name,
        const core::Time& min_value,
        const core::Time& max_value
    ) const;


    /*****************************************************/
    /* MIN QUERIES                                       */
    /*****************************************************/


    core::Value<double>
    get_min_double(
        const std::string& attribute_name
    ) const;



    core::Value<int>
    get_min_int(
        const std::string& attribute_name
    ) const;



    core::Value<std::string>
    get_min_string(
        const std::string& attribute_name
    ) const;


    core::Value<core::Time>
    get_min_time(
        const std::string& attribute_name
    ) const;


    /*****************************************************/
    /* MAX QUERIES                                       */
    /*****************************************************/


    core::Value<int>
    get_max_int(
        const std::string& attribute_name
    ) const;



    core::Value<double>
    get_max_double(
        const std::string& attribute_name
    ) const;



    core::Value<std::string>
    get_max_string(
        const std::string& attribute_name
    ) const;


    core::Value<core::Time>
    get_max_time(
        const std::string& attribute_name
    ) const;


  private:

    core::AttributeStore<OT>* attr_;

};



template <typename OT>
UserDefinedAttrs<OT>::
UserDefinedAttrs(core::AttributeStore<OT>* attr)
{
    attr_ = attr;
}

template <typename OT>
const core::Attribute*
UserDefinedAttrs<OT>::
add(
    std::string name,
    core::AttributeType t
)
{
    return attr_->add(name, t);
}

template <typename OT>
const core::Attribute*
UserDefinedAttrs<OT>::
add(
    std::unique_ptr<const core::Attribute> attribute
)
{
    return attr_->add(std::move(attribute));

}

template <typename OT>
const core::Attribute*
UserDefinedAttrs<OT>::
get(
    std::string name
) const
{
    return attr_->get(name);
}

template <typename OT>
void
UserDefinedAttrs<OT>::
set_as_string(
    const OT* oid,
    const std::string& attribute_name,
    const std::string& value
)
{
    return attr_->set_as_string(oid, attribute_name, value);
}



template <typename OT>
core::Value<std::string>
UserDefinedAttrs<OT>::
get_as_string(
    const OT* oid,
    const std::string& attribute_name
) const
{
    return attr_->get_as_string(oid, attribute_name);
}


template <typename OT>
bool
UserDefinedAttrs<OT>::
add_index(
    const std::string& attribute_name
)
{
    return attr_->add_index(attribute_name);
}



template <typename OT>
void
UserDefinedAttrs<OT>::
set_string(
    const OT* oid,
    const std::string& attribute_name,
    const std::string& value
)
{
    return attr_->set_string(oid, attribute_name, value);
}


template <typename OT>
void
UserDefinedAttrs<OT>::
set_double(
    const OT* oid,
    const std::string& attribute_name,
    double value
)
{
    return attr_->set_double(oid, attribute_name, value);
}



template <typename OT>
void
UserDefinedAttrs<OT>::
set_int(
    const OT* oid,
    const std::string& attribute_name,
    int value
)
{
    attr_->set_int(oid, attribute_name, value);
}


template <typename OT>
void
UserDefinedAttrs<OT>::
set_time(
    const OT* oid,
    const std::string& attribute_name,
    const core::Time& value
)
{
    attr_->set_time(oid, attribute_name, value);
}

template <typename OT>
void
UserDefinedAttrs<OT>::
set_text(
    const OT* oid,
    const std::string& attribute_name,
    const core::Text& value
)
{
    attr_->set_text(oid, attribute_name, value);
}


template <typename OT>
core::Value<std::string>
UserDefinedAttrs<OT>::
get_string(
    const OT* oid,
    const std::string& attribute_name
) const
{
    return attr_->get_string(oid, attribute_name);
}


template <typename OT>
core::Value<double>
UserDefinedAttrs<OT>::
get_double(
    const OT* oid,
    const std::string& attribute_name
) const
{
    return attr_->get_double(oid, attribute_name);
}


template <typename OT>
core::Value<int>
UserDefinedAttrs<OT>::
get_int(
    const OT* oid,
    const std::string& attribute_name
) const
{
    return attr_->get_int(oid, attribute_name);
}


template <typename OT>
core::Value<core::Time>
UserDefinedAttrs<OT>::
get_time(
    const OT* oid,
    const std::string& attribute_name
) const
{
    return attr_->get_time(oid, attribute_name);
}


template <typename OT>
core::Value<core::Text>
UserDefinedAttrs<OT>::
get_text(
    const OT* oid,
    const std::string& attribute_name
) const
{
    return attr_->get_text(oid, attribute_name);
}


/*****************************************************/
/* RANGE QUERIES                                     */
/*****************************************************/


template <typename OT>
std::vector<OT>
UserDefinedAttrs<OT>::
range_query_string(
    const std::string& attribute_name,
    const std::string& min_value,
    const std::string& max_value
) const
{
    return attr_->range_query_string(attribute_name, min_value, max_value);
}


template <typename OT>
std::vector<OT>
UserDefinedAttrs<OT>::
range_query_int(
    const std::string& attribute_name,
    const int& min_value,
    const int& max_value
) const
{
    return attr_->range_query_int(attribute_name, min_value, max_value);
}


template <typename OT>
std::vector<OT>
UserDefinedAttrs<OT>::
range_query_double(
    const std::string& attribute_name,
    const double& min_value,
    const double& max_value
) const
{
    return attr_->range_query_double(attribute_name, min_value, max_value);
}


template <typename OT>
std::vector<OT>
UserDefinedAttrs<OT>::
range_query_time(
    const std::string& attribute_name,
    const core::Time& min_value,
    const core::Time& max_value
) const
{
    return attr_->range_query_time(attribute_name, min_value, max_value);
}


/*****************************************************/
/* MIN QUERIES                                       */
/*****************************************************/


template <typename OT>
core::Value<double>
UserDefinedAttrs<OT>::
get_min_double(
    const std::string& attribute_name
) const
{
    return attr_->get_min_double(attribute_name);
}



template <typename OT>
core::Value<int>
UserDefinedAttrs<OT>::
get_min_int(
    const std::string& attribute_name
) const
{
    return attr_->get_min_int(attribute_name);
}



template <typename OT>
core::Value<std::string>
UserDefinedAttrs<OT>::
get_min_string(
    const std::string& attribute_name
) const
{
    return attr_->get_min_string(attribute_name);
}


template <typename OT>
core::Value<core::Time>
UserDefinedAttrs<OT>::
get_min_time(
    const std::string& attribute_name
) const
{
    return attr_->get_min_time(attribute_name);
}


/*****************************************************/
/* MAX QUERIES                                       */
/*****************************************************/


template <typename OT>
core::Value<int>
UserDefinedAttrs<OT>::
get_max_int(
    const std::string& attribute_name
) const
{
    return attr_->get_max_int(attribute_name);
}



template <typename OT>
core::Value<double>
UserDefinedAttrs<OT>::
get_max_double(
    const std::string& attribute_name
) const
{
    return attr_->get_max_double(attribute_name);
}



template <typename OT>
core::Value<std::string>
UserDefinedAttrs<OT>::
get_max_string(
    const std::string& attribute_name
) const
{
    return attr_->get_max_string(attribute_name);
}


template <typename OT>
core::Value<core::Time>
UserDefinedAttrs<OT>::
get_max_time(
    const std::string& attribute_name
) const
{
    return attr_->get_max_time(attribute_name);
}


}
}

#endif
