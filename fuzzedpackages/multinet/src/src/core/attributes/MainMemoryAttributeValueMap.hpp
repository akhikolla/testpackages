/**
 * History:
 * - 2018.01.01 file adapted from version 1.0 of the multinet library
 */

#ifndef UU_CORE_ATTRIBUTES_MAINMEMORYATTRIBUTESTORE_H_
#define UU_CORE_ATTRIBUTES_MAINMEMORYATTRIBUTESTORE_H_

#include <map>
#include <unordered_map>
#include "core/attributes/AttributeValueMap.hpp"

namespace uu {
namespace core {

/**
 * An AttributeValueMap implemented as main memory maps.
 */
template <typename ID>
class MainMemoryAttributeValueMap :
    public AttributeValueMap<ID>
{
  private:
    typedef AttributeValueMap<ID> super;

  public:

    MainMemoryAttributeValueMap();

    using super::add;
    using super::erase;
    using super::get;

    virtual
    const Attribute *
    add(
        std::unique_ptr<const Attribute> att
    ) override;


    void
    set_string(ID id,
               const std::string& attribute_name,
               const std::string& value
              ) override;


    void
    set_double(ID id,
               const std::string& attribute_name, double value
              ) override;


    void
    set_int(ID id,
            const std::string& attribute_name, int value
           ) override;


    void
    set_time(
        ID id,
        const std::string& attribute_name,
        const Time& value
    ) override;


    void
    set_text(
        ID id,
        const std::string& attribute_name,
        const Text& value
    ) override;


    Value<std::string>
    get_string(ID id,
               const std::string& attribute_name
              ) const override;


    Value<double>
    get_double(ID id,
               const std::string& attribute_name
              ) const override;


    Value<int>
    get_int(ID id,
            const std::string& attribute_name
           ) const override;


    Value<Time>
    get_time(
        ID oid,
        const std::string& attribute_name
    ) const override;


    Value<Text>
    get_text(
        ID oid,
        const std::string& attribute_name
    ) const override;


    bool
    reset(
        ID oid,
        const std::string& attribute_name
    ) override;


    bool
    add_index(
        const std::string& attribute_name
    ) override;


    std::vector<ID>
    range_query_string(
        const std::string& attribute_name,
        const std::string& min_value,
        const std::string& max_value
    ) const override;


    std::vector<ID>
    range_query_int(
        const std::string& attribute_name,
        const int& min_value,
        const int& max_value
    ) const override;


    std::vector<ID>
    range_query_double(
        const std::string& attribute_name,
        const double& min_value,
        const double& max_value
    ) const override;


    std::vector<ID>
    range_query_time(
        const std::string& attribute_name,
        const Time& min_value,
        const Time& max_value
    ) const override;


    Value<double>
    get_min_double(
        const std::string& attribute_name
    ) const override;


    Value<double>
    get_max_double(
        const std::string& attribute_name
    ) const override;


    Value<int>
    get_min_int(
        const std::string& attribute_name
    ) const override;


    Value<int>
    get_max_int(
        const std::string& attribute_name
    ) const override;


    Value<std::string>
    get_min_string(
        const std::string& attribute_name
    ) const override;


    Value<std::string>
    get_max_string(
        const std::string& attribute_name
    ) const override;


    Value<Time>
    get_min_time(
        const std::string& attribute_name
    ) const override;


    Value<Time>
    get_max_time(
        const std::string& attribute_name
    ) const override;

  private:

    /* These hash maps are structured as: map[AttributeName][ID][AttributeValue] */
    std::unordered_map<std::string, std::unordered_map<ID, std::string> > string_attribute;
    std::unordered_map<std::string, std::unordered_map<ID, double> > double_attribute;
    std::unordered_map<std::string, std::unordered_map<ID, int> > int_attribute;
    std::unordered_map<std::string, std::unordered_map<ID, Time> > time_attribute;
    std::unordered_map<std::string, std::unordered_map<ID, Text> > text_attribute;
    /* These hash map are structured as: map[AttributeName][AttributeValue][IDs]*/
    std::unordered_map<std::string, std::multimap<std::string, ID> > string_attribute_idx;
    std::unordered_map<std::string, std::multimap<double, ID> > double_attribute_idx;
    std::unordered_map<std::string, std::multimap<int, ID> > int_attribute_idx;
    std::unordered_map<std::string, std::multimap<Time, ID> > time_attribute_idx;
    std::unordered_map<std::string, std::multimap<Text, ID> > text_attribute_idx;
} ;




template <typename ID>
MainMemoryAttributeValueMap<ID>::
MainMemoryAttributeValueMap (
)
{
}



template <typename ID>
const Attribute *
MainMemoryAttributeValueMap<ID>::
add(
    std::unique_ptr<const Attribute> att
)
{
    // ASSERT is checked by the supertype container

    const Attribute * ptr;

    if (!(ptr = super::add(std::move(att))))
    {
        return nullptr;
    }

    switch (ptr->type)
    {
    case AttributeType::STRING:
        string_attribute[ptr->name] = std::unordered_map<ID, std::string>();
        break;

    case AttributeType::INTEGER:
        int_attribute[ptr->name] = std::unordered_map<ID, int>();
        break;

    case AttributeType::DOUBLE:
    case AttributeType::NUMERIC:
        double_attribute[ptr->name] = std::unordered_map<ID, double>();
        break;

    case AttributeType::TIME:
        time_attribute[ptr->name] = std::unordered_map<ID, Time>();
        break;

    case AttributeType::TEXT:
        text_attribute[ptr->name] = std::unordered_map<ID, Text>();
        break;

    }

    return ptr;
}


template <typename ID>
bool
MainMemoryAttributeValueMap<ID>::add_index (
    const std::string& attribute_name
)
{
    const Attribute* att = get(attribute_name);

    if (!att)
    {
        throw ElementNotFoundException("attribute " + attribute_name);
    }

    switch (att->type)
    {
    case AttributeType::STRING:
    {
        if (string_attribute_idx.count(attribute_name) > 0)
        {
            return false;
        }

        string_attribute_idx[attribute_name];

        // fill in index with existing values
        for (auto pair: string_attribute[attribute_name])
        {
            string_attribute_idx[attribute_name].insert(std::pair<std::string, ID>(pair.second, pair.first));
        }

        break;
    }

    case AttributeType::INTEGER:
    {
        if (int_attribute_idx.count(attribute_name) > 0)
        {
            return false;
        }

        int_attribute_idx[attribute_name];

        // fill in index with existing values
        for (auto pair: int_attribute[attribute_name])
        {
            int_attribute_idx[attribute_name].insert(std::pair<int, ID>(pair.second, pair.first));
        }

        break;
    }

    case AttributeType::DOUBLE:
    case AttributeType::NUMERIC:
    {
        if (double_attribute_idx.count(attribute_name) > 0)
        {
            return false;
        }

        double_attribute_idx[attribute_name];

        // fill in index with existing values
        for (auto pair: double_attribute[attribute_name])
        {
            double_attribute_idx[attribute_name].insert(std::pair<double, ID>(pair.second, pair.first));
        }

        break;
    }

    case AttributeType::TIME:
    {
        if (time_attribute_idx.count(attribute_name) > 0)
        {
            return false;
        }

        time_attribute_idx[attribute_name];

        // fill in index with existing values
        for (auto pair: time_attribute[attribute_name])
        {
            time_attribute_idx[attribute_name].insert(std::pair<Time, ID>(pair.second, pair.first));
        }

        break;
    }

    case AttributeType::TEXT:
    {
        if (text_attribute_idx.count(attribute_name) > 0)
        {
            return false;
        }

        text_attribute_idx[attribute_name];

        for (auto pair : text_attribute[attribute_name])
        {
            text_attribute_idx[attribute_name].insert(std::pair<Text, ID>(pair.second, pair.first));
        }

        break;
    }

    }

    return true;
}

template <typename ID>
void
MainMemoryAttributeValueMap<ID>::set_double (
    ID oid,
    const std::string& attribute_name,
    double val
)
{
    auto values = double_attribute.find(attribute_name);

    if (values == double_attribute.end())
    {
        throw ElementNotFoundException("double attribute " + attribute_name);
    }

    auto ret = values->second.insert(std::make_pair(oid,val));

    if (ret.second == false)
    {
        // oid already existed
        ret.first->second = val;
    }

    // if index is present
    auto index = double_attribute_idx.find(attribute_name);

    if (index != double_attribute_idx.end())
    {
        index->second.insert(std::make_pair(val, oid));

        if (ret.second==false)
        {
            // oid already existed, with value ret.first->second
            double old_val = ret.first->second;
            auto pairs = index->second.equal_range(old_val);

            for (auto it = pairs.first; it != pairs.second; ++it)
            {
                if (it->second == oid)
                {
                    index->second.erase(it);
                    break;
                }
            }
        }
    }
}

template <typename ID>
void
MainMemoryAttributeValueMap<ID>::set_int (
    ID oid,
    const std::string& attribute_name,
    int val
)
{
    auto values = int_attribute.find(attribute_name);

    if (values == int_attribute.end())
    {
        throw ElementNotFoundException("int attribute " + attribute_name);
    }

    auto ret = values->second.insert(std::make_pair(oid, val));

    if (ret.second == false)
    {
        // oid already existed
        ret.first->second = val;
    }

    // if index is present
    auto index = int_attribute_idx.find(attribute_name);

    if (index != int_attribute_idx.end())
    {
        index->second.insert(std::make_pair(val, oid));

        if (ret.second==false)
        {
            // oid already existed, with value ret.first->second
            int old_val = ret.first->second;
            auto pairs = index->second.equal_range(old_val);

            for (auto it = pairs.first; it != pairs.second; ++it)
            {
                if (it->second == oid)
                {
                    index->second.erase(it);
                    break;
                }
            }
        }
    }
}

template <typename ID>
void
MainMemoryAttributeValueMap<ID>::set_string (
    ID oid,
    const std::string& attribute_name,
    const std::string& val
)
{
    auto values = string_attribute.find(attribute_name);

    if (values == string_attribute.end())
    {
        throw ElementNotFoundException("string attribute " + attribute_name);
    }

    auto ret = values->second.insert(std::make_pair(oid, val));

    if (ret.second == false)
    {
        // oid already existed
        ret.first->second = val;
    }

    // if index is present
    auto index = string_attribute_idx.find(attribute_name);

    if (index != string_attribute_idx.end())
    {
        index->second.insert(std::make_pair(val, oid));

        if (ret.second==false)
        {
            // oid already existed, with value ret.first->second
            std::string old_val = ret.first->second;
            auto pairs = index->second.equal_range(old_val);

            for (auto it = pairs.first; it != pairs.second; ++it)
            {
                if (it->second == oid)
                {
                    index->second.erase(it);
                    break;
                }
            }
        }
    }
}

template <typename ID>
void
MainMemoryAttributeValueMap<ID>::set_time (
    ID oid,
    const std::string& attribute_name,
    const Time& val
)
{
    auto values = time_attribute.find(attribute_name);

    if (values == time_attribute.end())
    {
        throw ElementNotFoundException("time attribute " + attribute_name);
    }

    auto ret = values->second.insert(std::make_pair(oid, val));

    if (ret.second == false)
    {
        // oid already existed
        ret.first->second = val;
    }

    // if index is present
    auto index = time_attribute_idx.find(attribute_name);

    if (index != time_attribute_idx.end())
    {
        index->second.insert(std::make_pair(val, oid));

        if (ret.second==false)
        {
            // oid already existed, with value ret.first->second
            Time old_val = ret.first->second;
            auto pairs = index->second.equal_range(old_val);

            for (auto it = pairs.first; it != pairs.second; ++it)
            {
                if (it->second == oid)
                {
                    index->second.erase(it);
                    break;
                }
            }
        }
    }
}

template <typename ID>
void
MainMemoryAttributeValueMap<ID>::set_text (
    ID oid,
    const std::string& attribute_name,
    const Text& val
)
{
    auto values = text_attribute.find(attribute_name);

    if (values == text_attribute.end())
    {
        throw ElementNotFoundException("text attribute " + attribute_name);
    }

    auto ret = values->second.insert(std::make_pair(oid, val));

    if (ret.second == false)
    {
        // oid already existed
        ret.first->second = val;
    }

    auto index = text_attribute_idx.find(attribute_name);

    if (index != text_attribute_idx.end())
    {
        index->second.insert(std::make_pair(val, oid));

        if (ret.second==false)
        {
            // oid already existed, with value ret.first->second
            Text old_val = ret.first->second;
            auto pairs = index->second.equal_range(old_val);

            for (auto it = pairs.first; it != pairs.second; ++it)
            {
                if (it->second == oid)
                {
                    index->second.erase(it);
                    break;
                }
            }
        }
    }

}

template <typename ID>
Value<double>
MainMemoryAttributeValueMap<ID>::get_double (
    ID oid,
    const std::string& attribute_name
) const
{
    auto attr = double_attribute.find(attribute_name);

    if (attr == double_attribute.end())
    {
        throw ElementNotFoundException("double attribute " + attribute_name);
    }

    auto val = attr->second.find(oid);

    if (val == attr->second.end())
    {
        return Value<double>(0, true);
    }

    return Value<double>(val->second, false);
}

template <typename ID>
Value<int>
MainMemoryAttributeValueMap<ID>::get_int (
    ID oid,
    const std::string& attribute_name
) const
{
    auto attr = int_attribute.find(attribute_name);

    if (attr == int_attribute.end())
    {
        throw ElementNotFoundException("integer attribute " + attribute_name);
    }

    auto val = attr->second.find(oid);

    if (val == attr->second.end())
    {
        return Value<int>(0, true);
    }

    return Value<int>(val->second, false);
}

template <typename ID>
Value<std::string>
MainMemoryAttributeValueMap<ID>::get_string (
    ID oid,
    const std::string& attribute_name
) const
{
    auto attr = string_attribute.find(attribute_name);

    if (attr == string_attribute.end())
    {
        throw ElementNotFoundException("string attribute " + attribute_name);
    }

    auto val = attr->second.find(oid);

    if (val == attr->second.end())
    {
        return Value<std::string>("", true);
    }

    return Value<std::string>(val->second, false);
}

template <typename ID>
Value<Time>
MainMemoryAttributeValueMap<ID>::get_time (
    ID oid,
    const std::string& attribute_name
) const
{
    auto attr = time_attribute.find(attribute_name);

    if (attr == time_attribute.end())
    {
        throw ElementNotFoundException("time attribute " + attribute_name);
    }

    auto val = attr->second.find(oid);

    if (val == attr->second.end())
    {
        return Value<Time>(Time(), true);
    }

    return Value<Time>(val->second, false);
}

template <typename ID>
Value<Text>
MainMemoryAttributeValueMap<ID>::get_text(
    ID oid,
    const std::string& attribute_name
) const
{
    auto attr = text_attribute.find(attribute_name);

    if (attr == text_attribute.end())
    {
        throw ElementNotFoundException("text attribute " + attribute_name);
    }

    auto val = attr->second.find(oid);

    if (val == attr->second.end())
    {
        return Value<Text>(Text(), true);
    }

    return Value<Text>(val->second, false);
}

template <typename ID>
Value<double>
MainMemoryAttributeValueMap<ID>::get_min_double (
    const std::string& attribute_name
) const
{
    auto attr = double_attribute.find(attribute_name);

    if (attr == double_attribute.end())
    {
        throw ElementNotFoundException("double attribute " + attribute_name);
    }


    // if an index is present:
    if (double_attribute_idx.count(attribute_name) > 0)
    {
        auto index = double_attribute_idx.at(attribute_name);

        if (index.empty())
        {
            return Value<double>(0.0, true);
        }

        auto iterator = index.begin();

        return Value<double>(iterator->first, false);
    }

    // otherwise, we must scan all the values
    if (attr->second.empty())
    {
        return Value<double>(0.0,true);
    }

    double min_value = attr->second.begin()->second;

    for (auto pair: attr->second)
    {
        if (pair.second < min_value)
        {
            min_value = pair.second;
        }
    }

    return Value<double>(min_value, false);
}

template <typename ID>
Value<double>
MainMemoryAttributeValueMap<ID>::get_max_double (
    const std::string& attribute_name
) const
{
    auto attr = double_attribute.find(attribute_name);

    if (attr == double_attribute.end())
    {
        throw ElementNotFoundException("double attribute " + attribute_name);
    }

    // if an index is present:
    if (double_attribute_idx.count(attribute_name) > 0)
    {
        auto index = double_attribute_idx.at(attribute_name);

        if (index.empty())
        {
            return Value<double>(0.0, true);
        }

        auto iterator = index.end();
        iterator--;

        return Value<double>(iterator->first, false);
    }

    // otherwise, we must scan all the values
    if (attr->second.empty())
    {
        return Value<double>(0.0,true);
    }

    double max_value = attr->second.begin()->second;

    for (auto pair: attr->second)
    {
        if (pair.second > max_value)
        {
            max_value = pair.second;
        }
    }

    return Value<double>(max_value, false);
}

template <typename ID>
Value<int>
MainMemoryAttributeValueMap<ID>::get_min_int (
    const std::string& attribute_name
) const
{
    auto attr = int_attribute.find(attribute_name);

    if (attr == int_attribute.end())
    {
        throw ElementNotFoundException("int attribute " + attribute_name);
    }


    // if an index is present:
    if (int_attribute_idx.count(attribute_name) > 0)
    {
        auto index = int_attribute_idx.at(attribute_name);

        if (index.empty())
        {
            return Value<int>(0, true);
        }

        auto iterator = index.begin();

        return Value<int>(iterator->first, false);
    }

    // otherwise, we must scan all the values
    if (attr->second.empty())
    {
        return Value<int>(0,true);
    }

    int min_value = attr->second.begin()->second;

    for (auto pair: attr->second)
    {
        if (pair.second < min_value)
        {
            min_value = pair.second;
        }
    }

    return Value<int>(min_value, false);
}

template <typename ID>
Value<int>
MainMemoryAttributeValueMap<ID>::get_max_int (
    const std::string& attribute_name
) const
{
    auto attr = int_attribute.find(attribute_name);

    if (attr == int_attribute.end())
    {
        throw ElementNotFoundException("int attribute " + attribute_name);
    }

    // if an index is present:
    if (int_attribute_idx.count(attribute_name) > 0)
    {
        auto index = int_attribute_idx.at(attribute_name);

        if (index.empty())
        {
            return Value<int>(0, true);
        }

        auto iterator = index.end();
        iterator--;

        return Value<int>(iterator->first, false);
    }

    // otherwise, we must scan all the values
    if (attr->second.empty())
    {
        return Value<int>(0,true);
    }

    int max_value = attr->second.begin()->second;

    for (auto pair: attr->second)
    {
        if (pair.second > max_value)
        {
            max_value = pair.second;
        }
    }

    return Value<int>(max_value, false);
}

template <typename ID>
Value<std::string>
MainMemoryAttributeValueMap<ID>::get_min_string (
    const std::string& attribute_name
) const
{
    auto attr = string_attribute.find(attribute_name);

    if (attr == string_attribute.end())
    {
        throw ElementNotFoundException("string attribute " + attribute_name);
    }

    // if an index is present:
    if (string_attribute_idx.count(attribute_name) > 0)
    {
        auto index = string_attribute_idx.at(attribute_name);

        if (index.empty())
        {
            return Value<std::string>("", true);
        }

        auto iterator = index.begin();

        return Value<std::string>(iterator->first, false);
    }

    // otherwise, we must scan all the values
    if (attr->second.empty())
    {
        return Value<std::string>("",true);
    }

    std::string min_value = attr->second.begin()->second;

    for (auto pair: attr->second)
    {
        if (pair.second < min_value)
        {
            min_value = pair.second;
        }
    }

    return Value<std::string>(min_value, false);
}

template <typename ID>
Value<std::string>
MainMemoryAttributeValueMap<ID>::get_max_string (
    const std::string& attribute_name
) const
{
    auto attr = string_attribute.find(attribute_name);

    if (attr == string_attribute.end())
    {
        throw ElementNotFoundException("string attribute " + attribute_name);
    }

    // if an index is present:
    if (string_attribute_idx.count(attribute_name) > 0)
    {
        auto index = string_attribute_idx.at(attribute_name);

        if (index.empty())
        {
            return Value<std::string>("", true);
        }

        auto iterator = index.end();
        iterator--;

        return Value<std::string>(iterator->first, false);
    }

    // otherwise, we must scan all the values
    if (attr->second.empty())
    {
        return Value<std::string>("",true);
    }

    std::string max_value = attr->second.begin()->second;

    for (auto pair: attr->second)
    {
        if (pair.second > max_value)
        {
            max_value = pair.second;
        }
    }

    return Value<std::string>(max_value, false);
}

template <typename ID>
Value<Time>
MainMemoryAttributeValueMap<ID>::get_min_time (
    const std::string& attribute_name
) const
{
    auto attr = time_attribute.find(attribute_name);

    if (attr == time_attribute.end())
    {
        throw ElementNotFoundException("time attribute " + attribute_name);
    }

    // if an index is present:
    if (time_attribute_idx.count(attribute_name) > 0)
    {
        auto index = time_attribute_idx.at(attribute_name);

        if (index.empty())
        {
            return Value<Time>();
        }

        auto iterator = index.begin();

        return Value<Time>(iterator->first, false);
    }

    // otherwise, we must scan all the values
    if (attr->second.empty())
    {
        return Value<Time>(Time(),true);
    }

    Time min_value = attr->second.begin()->second;

    for (auto pair: attr->second)
    {
        if (pair.second < min_value)
        {
            min_value = pair.second;
        }
    }

    return Value<Time>(min_value, false);
}

template <typename ID>
Value<Time>
MainMemoryAttributeValueMap<ID>::get_max_time (
    const std::string& attribute_name
) const
{
    auto attr = time_attribute.find(attribute_name);

    if (attr == time_attribute.end())
    {
        throw ElementNotFoundException("time attribute " + attribute_name);
    }

    // if an index is present:
    if (time_attribute_idx.count(attribute_name) > 0)
    {
        auto index = time_attribute_idx.at(attribute_name);

        if (index.empty())
        {
            return Value<Time>();
        }

        auto iterator = index.end();
        iterator--;

        return Value<Time>(iterator->first, false);
    }

    // otherwise, we must scan all the values
    if (attr->second.empty())
    {
        return Value<Time>(Time(),true);
    }

    Time max_value = attr->second.begin()->second;

    for (auto pair: attr->second)
    {
        if (pair.second > max_value)
        {
            max_value = pair.second;
        }
    }

    return Value<Time>(max_value, false);
}

template <typename ID>
std::vector<ID>
MainMemoryAttributeValueMap<ID>::range_query_string (
    const std::string& attribute_name,
    const std::string& min_value,
    const std::string& max_value
) const
{
    std::vector<ID> result;

    auto attr = string_attribute.find(attribute_name);

    if (attr == string_attribute.end())
    {
        throw ElementNotFoundException("string attribute " + attribute_name);
    }

    // if an index is present:
    if (string_attribute_idx.count(attribute_name) > 0)
    {
        auto index = string_attribute_idx.at(attribute_name);

        auto itlow = index.lower_bound(min_value);
        auto itup = index.upper_bound(max_value);

        for (auto it = itlow; it != itup; ++it)
        {
            result.push_back((*it).second);
        }

        return result;
    }

    // otherwise, we must scan all the values
    for (auto pair: attr->second)
    {
        if (pair.second >= min_value && pair.second <= max_value)
        {
            result.push_back(pair.first);
        }
    }

    return result;
}

template <typename ID>
std::vector<ID>
MainMemoryAttributeValueMap<ID>::range_query_int (
    const std::string& attribute_name,
    const int& min_value,
    const int& max_value
) const
{
    std::vector<ID> result;

    auto attr = int_attribute.find(attribute_name);

    if (attr == int_attribute.end())
    {
        throw ElementNotFoundException("int attribute " + attribute_name);
    }

    // if an index is present:
    if (int_attribute_idx.count(attribute_name) > 0)
    {
        auto index = int_attribute_idx.at(attribute_name);

        auto itlow = index.lower_bound(min_value);
        auto itup = index.upper_bound(max_value);

        for (auto it = itlow; it != itup; ++it)
        {
            result.push_back((*it).second);
        }

        return result;
    }

    // otherwise, we must scan all the values
    for (auto pair: attr->second)
    {
        if (pair.second >= min_value && pair.second <= max_value)
        {
            result.push_back(pair.first);
        }
    }

    return result;
}

template <typename ID>
std::vector<ID>
MainMemoryAttributeValueMap<ID>::range_query_double (
    const std::string& attribute_name,
    const double& min_value,
    const double& max_value
) const
{
    std::vector<ID> result;

    auto attr = double_attribute.find(attribute_name);

    if (attr == double_attribute.end())
    {
        throw ElementNotFoundException("double attribute " + attribute_name);
    }

    // if an index is present:
    if (double_attribute_idx.count(attribute_name) > 0)
    {
        auto index = double_attribute_idx.at(attribute_name);

        auto itlow = index.lower_bound(min_value);
        auto itup = index.upper_bound(max_value);

        for (auto it = itlow; it != itup; ++it)
        {
            result.push_back((*it).second);
        }

        return result;
    }

    // otherwise, we must scan all the values
    for (auto pair: attr->second)
    {
        {
            result.push_back(pair.first);
        }
    }

    return result;

}

template <typename ID>
std::vector<ID>
MainMemoryAttributeValueMap<ID>::range_query_time (
    const std::string& attribute_name,
    const Time& min_value,
    const Time& max_value
) const
{
    std::vector<ID> result;

    auto attr = time_attribute.find(attribute_name);

    if (attr == time_attribute.end())
    {
        throw ElementNotFoundException("time attribute " + attribute_name);
    }

    // if an index is present:
    if (time_attribute_idx.count(attribute_name) > 0)
    {
        auto index = time_attribute_idx.at(attribute_name);

        auto itlow = index.lower_bound(min_value);
        auto itup = index.upper_bound(max_value);

        for (auto it = itlow; it != itup; ++it)
        {
            result.push_back((*it).second);
        }

        return result;
    }

    // otherwise, we must scan all the values
    for (auto pair: attr->second)
    {
        if (pair.second >= min_value && pair.second <= max_value)
        {
            result.push_back(pair.first);
        }
    }

    return result;

}

template <typename ID>
bool
MainMemoryAttributeValueMap<ID>::reset (
    ID id,
    const std::string& attribute_name
)
{
    const Attribute* att = get(attribute_name);

    if (!att)
    {
        throw ElementNotFoundException("attribute " + attribute_name);
    }

    switch (att->type)
    {
    case AttributeType::STRING:
    {
        auto values = string_attribute.find(attribute_name);

        auto old_entry = values->second.find(id);

        if (old_entry == values->second.end())
        {
            return false;
        }

        if (string_attribute_idx.count(attribute_name) > 0)
        {
            // if there is an index
            auto index = string_attribute_idx.at(attribute_name);

            auto pairs = index.equal_range(old_entry->second);

            for (auto it = pairs.first; it != pairs.second; ++it)
            {
                if (it->second == id)
                {
                    index.erase(it);
                    break;
                }
            }
        }

        values->second.erase(old_entry);
        return true;
    }

    case AttributeType::INTEGER:
    {
        auto values = int_attribute.find(attribute_name);

        auto old_entry = values->second.find(id);

        if (old_entry == values->second.end())
        {
            return false;
        }

        if (int_attribute_idx.count(attribute_name) > 0)
        {
            // if there is an index
            auto index = int_attribute_idx.at(attribute_name);

            auto pairs = index.equal_range(old_entry->second);

            for (auto it = pairs.first; it != pairs.second; ++it)
            {
                if (it->second == id)
                {
                    index.erase(it);
                    break;
                }
            }
        }

        values->second.erase(old_entry);
        return true;
    }

    case AttributeType::DOUBLE:
    case AttributeType::NUMERIC:
    {
        auto values = double_attribute.find(attribute_name);

        auto old_entry = values->second.find(id);

        if (old_entry == values->second.end())
        {
            return false;
        }

        if (double_attribute_idx.count(attribute_name) > 0)
        {
            // if there is an index
            auto index = double_attribute_idx.at(attribute_name);

            auto pairs = index.equal_range(old_entry->second);

            for (auto it = pairs.first; it != pairs.second; ++it)
            {
                if (it->second == id)
                {
                    index.erase(it);
                    break;
                }
            }
        }

        values->second.erase(old_entry);
        return true;
    }

    case AttributeType::TIME:
    {
        auto values = time_attribute.find(attribute_name);

        auto old_entry = values->second.find(id);

        if (old_entry == values->second.end())
        {
            return false;
        }

        if (time_attribute_idx.count(attribute_name) > 0)
        {
            // if there is an index
            auto index = time_attribute_idx.at(attribute_name);

            auto pairs = index.equal_range(old_entry->second);

            for (auto it = pairs.first; it != pairs.second; ++it)
            {
                if (it->second == id)
                {
                    index.erase(it);
                    break;
                }
            }
        }

        values->second.erase(old_entry);
        return true;
    }

    case AttributeType::TEXT:
    {
        auto values = text_attribute.find(attribute_name);

        auto old_entry = values->second.find(id);

        if (old_entry == values->second.end())
        {
            return false;
        }

        if (text_attribute_idx.count(attribute_name) > 0)
        {
            // if there is an index
            auto index = text_attribute_idx.at(attribute_name);

            auto pairs = index.equal_range(old_entry->second);

            for (auto it = pairs.first; it != pairs.second; ++it)
            {
                if (it->second == id)
                {
                    index.erase(it);
                    break;
                }
            }
        }

        values->second.erase(old_entry);
        return true;
    }

    }

    return false; // never gets here
}


}
}

#endif
