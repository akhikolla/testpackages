/**
 *
 * History:
 * - 2018.01.01 file adapted from version 1.0 of the multinet library
 */

#ifndef UU_OLAP_CCUBE_H_
#define UU_OLAP_CCUBE_H_

#include "core/stores/AttributeStore.hpp"
#include "core/exceptions/ElementNotFoundException.hpp"
#include "core/exceptions/WrongParameterException.hpp"
#include "core/exceptions/OutOfBoundsException.hpp"
#include "core/datastructures/containers/UnionSortedRandomSet.hpp"
#include "core/olap/datastructures/NCube.hpp"
#include "core/olap/selection/IndexIterator.hpp"
#include <string>
#include <unordered_map>
#include <vector>

namespace uu {
namespace core {



template <class CONTAINER_TYPE>
class CCube
{

  private:

    NCube<std::shared_ptr<CONTAINER_TYPE>> data_;

  public:

    typedef CONTAINER_TYPE* entry_type;
    typedef CONTAINER_TYPE container_type;
    typedef typename CONTAINER_TYPE::value_type value_type;

    /**
     * Creates a cube specifying its dimensions and the members for each dimension.
     */
    template <class Iterator>
    CCube(
        const std::vector<std::string>& dim,
        const std::vector<std::vector<std::string>>& members,
        Iterator begin,
        Iterator end
    );


    /**
     * Destroys the cube.
     */
    ~CCube()
    {

    }

    /**
     * Returns the order (number of dimensions) of this cube.
     */
    size_t
    order(
    ) const;

    /**
     * Returns the size of the cube, for each dimension.
     */
    std::vector<size_t>
    size(
    ) const;

    const core::UnionSortedRandomSet<typename CONTAINER_TYPE::value_type>*
    elements(
    ) const;


    core::UnionSortedRandomSet<typename CONTAINER_TYPE::value_type>*
    elements(
    );

    core::AttributeStore<typename CONTAINER_TYPE::value_type>*
    attr(
    ) const;

    /** Returns a const iterator to the first object in the cube */
    typename std::vector<std::shared_ptr<const CONTAINER_TYPE>>::const_iterator
            begin(
            ) const;

    /** Returns a const iterator after the last object in the cube */
    typename std::vector<std::shared_ptr<const CONTAINER_TYPE>>::const_iterator
            end(
            ) const;

    /** Returns an iterator to the first object in the cube */
    typename std::vector<std::shared_ptr<CONTAINER_TYPE>>::iterator
            begin(
            );

    /** Returns an iterator after the last object in the cube */
    typename std::vector<std::shared_ptr<CONTAINER_TYPE>>::iterator
            end(
            );

    /**
     * Returns the dimensions of this cube.
     */
    const std::vector<std::string>&
    dim(
    ) const;

    /**
     * Returns the members of a dimension.
     */
    const std::vector<std::string>&
    members(
        const std::string& dim
    ) const;


    /**
     * Returns the object at the given position in the cube.
     * @throw ElementNotFoundException if the index is outside the bounds on the cube
     */
    CONTAINER_TYPE*
    operator[](
        const std::vector<size_t>& index
    );

    /**
     * Returns the object at the given position in the cube.
     * @throw ElementNotFoundException if the index is outside the bounds on the cube
     */
    CONTAINER_TYPE*
    operator[](
        const std::vector<std::string>& index
    );

    /**
     * Returns the object at the given position in the cube.
     * @throw ElementNotFoundException if the index is outside the bounds on the cube
     */
    const CONTAINER_TYPE*
    operator[](
        const std::vector<size_t>& index
    ) const;

    /**
     * Returns the object at the given position in the cube.
     * @throw ElementNotFoundException if the index is outside the bounds on the cube
     */
    const CONTAINER_TYPE*
    operator[](
        const std::vector<std::string>& index
    ) const;


    /** Returns the object at the given position in the cube.
     * @throw ElementNotFoundException if the index is outside the bounds on the cube
     */
    CONTAINER_TYPE*
    at(
        const std::vector<size_t>& index
    );

    /** Returns the object at the given position in the cube.
     * @throw ElementNotFoundException if the index is outside the bounds on the cube
     */
    CONTAINER_TYPE*
    at(
        const std::vector<std::string>& index
    );

    /** Returns the object at the given position in the cube.
     * @throw ElementNotFoundException if the index is outside the bounds on the cube
     */
    const CONTAINER_TYPE*
    at(
        const std::vector<size_t>& index
    ) const;

    /** Returns the object at the given position in the cube.
     * @throw ElementNotFoundException if the index is outside the bounds on the cube
     */
    const CONTAINER_TYPE*
    at(
        const std::vector<std::string>& index
    ) const;



    /* index of a dimension */
    size_t
    index_of(
        const std::string& dim
    ) const;

    /* index of a member given a dimension */
    size_t
    index_of(
        const std::string& dim,
        const std::string& member
    ) const;

    /* computes a numerical index from a label-based index */
    std::vector<size_t>
    index_of(
        const std::vector<std::string>& index
    ) const;

    // OPERATORS

    template <typename Iterator>
    void
    resize(
        const std::string& dimension,
        const std::string& member,
        Iterator begin,
        Iterator end
    );

  private:

    CCube(
        const std::vector<std::string>& dim,
        const std::vector<std::vector<std::string>>& members
    );

    void
    insert(
        const std::shared_ptr<CONTAINER_TYPE>& value,
        const std::vector<size_t>& index
    );

    void
    insert(
        const std::shared_ptr<CONTAINER_TYPE>& value,
        const std::vector<std::string>& index
    );

    // all the vertices in the cube
    std::unique_ptr<core::UnionSortedRandomSet<typename CONTAINER_TYPE::value_type>> elements_;

    // element attributes
    std::unique_ptr<core::AttributeStore<typename CONTAINER_TYPE::value_type>> attr_;

};

// TEMPLATE DEFINITIONS

template <class CONTAINER_TYPE>
CCube<CONTAINER_TYPE>::
CCube(
    const std::vector<std::string>& dim,
    const std::vector<std::vector<std::string>>& members
) : data_(dim, members)
{

    // Elements in the cube
    elements_ = std::make_unique<core::UnionSortedRandomSet<typename CONTAINER_TYPE::value_type>>();

    // Element attributes
    attr_ = std::make_unique<core::AttributeStore<typename CONTAINER_TYPE::value_type>>();

    elements_->attach(attr_.get());
}


template <class CONTAINER_TYPE>
template <class Iterator>
CCube<CONTAINER_TYPE>::
CCube(
    const std::vector<std::string>& dim,
    const std::vector<std::vector<std::string>>& members,
    Iterator begin,
    Iterator end
) : data_(dim, members, begin, end)
{


    // Elements in the cube
    elements_ = std::make_unique<core::UnionSortedRandomSet<typename CONTAINER_TYPE::value_type>>();

    // Element attributes
    attr_ = std::make_unique<core::AttributeStore<typename CONTAINER_TYPE::value_type>>();

    elements_->attach(attr_.get());

    for (auto cont = begin; cont != end; ++cont)
    {
        (*cont)->attach(elements_.get());

        // Add all existing objects in the containers to the elements
        for (auto obj: *(*cont))
        {
            elements_->notify_add(obj);
        }
    }
}


template <class CONTAINER_TYPE>
size_t
CCube<CONTAINER_TYPE>::
order(
) const
{
    return data_.order();
}


template <class CONTAINER_TYPE>
std::vector<size_t>
CCube<CONTAINER_TYPE>::
size(
) const
{
    return data_.size();
}


template <class CONTAINER_TYPE>
const core::UnionSortedRandomSet<typename CONTAINER_TYPE::value_type>*
CCube<CONTAINER_TYPE>::
elements(
) const
{
    return elements_.get();
}


template <class CONTAINER_TYPE>
core::UnionSortedRandomSet<typename CONTAINER_TYPE::value_type>*
CCube<CONTAINER_TYPE>::
elements(
)
{
    return elements_.get();
}

template <class CONTAINER_TYPE>
core::AttributeStore<typename CONTAINER_TYPE::value_type>*
CCube<CONTAINER_TYPE>::
attr(
) const
{
    return attr_.get();
}

template <class CONTAINER_TYPE>
typename std::vector<std::shared_ptr<const CONTAINER_TYPE>>::const_iterator
        CCube<CONTAINER_TYPE>::
        begin(
        ) const
{
    return data_.begin();
}

template <class CONTAINER_TYPE>
typename std::vector<std::shared_ptr<const CONTAINER_TYPE>>::const_iterator
        CCube<CONTAINER_TYPE>::
        end(
        ) const
{
    return data_.end();
}

template <class CONTAINER_TYPE>
typename std::vector<std::shared_ptr<CONTAINER_TYPE>>::iterator
        CCube<CONTAINER_TYPE>::
        begin(
        )
{
    return data_.begin();
}

template <class CONTAINER_TYPE>
typename std::vector<std::shared_ptr<CONTAINER_TYPE>>::iterator
        CCube<CONTAINER_TYPE>::
        end(
        )
{
    return data_.end();
}

template <class CONTAINER_TYPE>
const std::vector<std::string>&
CCube<CONTAINER_TYPE>::
dim(
) const
{
    return data_.dim();
}

template <class CONTAINER_TYPE>
const std::vector<std::string>&
CCube<CONTAINER_TYPE>::
members(
    const std::string& dim
) const
{
    return data_.members(dim);
}


template <class CONTAINER_TYPE>
CONTAINER_TYPE*
CCube<CONTAINER_TYPE>::
operator[](
    const std::vector<size_t>& index
)
{
    return data_[index].get();
}


template <class CONTAINER_TYPE>
CONTAINER_TYPE*
CCube<CONTAINER_TYPE>::
operator[](
    const std::vector<std::string>& index
)
{
    return data_[index].get();
}

template <class CONTAINER_TYPE>
const CONTAINER_TYPE*
CCube<CONTAINER_TYPE>::
operator[](
    const std::vector<size_t>& index
) const
{
    return data_[index].get();
}


template <class CONTAINER_TYPE>
const CONTAINER_TYPE*
CCube<CONTAINER_TYPE>::
operator[](
    const std::vector<std::string>& index
) const
{
    return data_[index].get();
}

template <class CONTAINER_TYPE>
CONTAINER_TYPE*
CCube<CONTAINER_TYPE>::
at(
    const std::vector<size_t>& index
)
{
    return data_.at(index).get();
}


template <class CONTAINER_TYPE>
CONTAINER_TYPE*
CCube<CONTAINER_TYPE>::
at(
    const std::vector<std::string>& index
)
{
    return data_.at(index).get();
}

template <class CONTAINER_TYPE>
const CONTAINER_TYPE*
CCube<CONTAINER_TYPE>::
at(
    const std::vector<size_t>& index
) const
{
    return data_.at(index).get();
}


template <class CONTAINER_TYPE>
const CONTAINER_TYPE*
CCube<CONTAINER_TYPE>::
at(
    const std::vector<std::string>& index
) const
{
    return data_.at(index).get();
}


/* index of a dimension */
template <class CONTAINER_TYPE>
size_t
CCube<CONTAINER_TYPE>::
index_of(
    const std::string& dim
) const
{
    return data_.index_of(dim);
}

/* index of a member given a dimension */
template <class CONTAINER_TYPE>
size_t
CCube<CONTAINER_TYPE>::
index_of(
    const std::string& dim,
    const std::string& member
) const
{
    return data_.index_of(dim, member);
}

/*
template <class CONTAINER_TYPE>
void
CCube<CONTAINER_TYPE>::
insert(
       CONTAINER_TYPE* value,
       const std::vector<size_t>& index
       )
{
    value->attach(elements_.get());
    super::insert(value, index);
}

template <class CONTAINER_TYPE>
void
CCube<CONTAINER_TYPE>::
insert(
       CONTAINER_TYPE* value,
       const std::vector<std::string>& index
       )
{
    value->attach(elements_.get());
    super::insert(value, index);
}

 */

template <class CONTAINER_TYPE>
void
CCube<CONTAINER_TYPE>::
insert(
    const std::shared_ptr<CONTAINER_TYPE>& value,
    const std::vector<size_t>& index
)
{
    value->attach(elements_.get());
    data_.at[index] = value;
}

template <class CONTAINER_TYPE>
void
CCube<CONTAINER_TYPE>::
insert(
    const std::shared_ptr<CONTAINER_TYPE>& value,
    const std::vector<std::string>& index
)
{
    value->attach(elements_.get());
    data_.at[index] = value;
}

template <typename CONTAINER_TYPE>
template <typename Iterator>
void
CCube<CONTAINER_TYPE>::
resize(
    const std::string& dimension,
    const std::string& new_member,
    Iterator begin,
    Iterator end
)
{
    size_t dim_idx = index_of(dimension);

    std::vector<std::vector<std::string>> new_members;

    for (auto d: dim())
    {
        new_members.push_back(std::vector<std::string>());

        for (auto m: members(d))
        {
            new_members.back().push_back(m);
        }
    }

    new_members.at(dim_idx).push_back(new_member);
    size_t member_idx = new_members.at(dim_idx).size()-1;

    auto new_size = size();
    new_size.at(dim_idx)++;

    uu::core::sel::IndexIterator idx(new_size);

    std::vector<std::shared_ptr<CONTAINER_TYPE>> cells;

    auto it = begin;

    for (auto index: idx)
    {
        if (index.at(dim_idx) == member_idx)
        {
            if (it==end)
            {
                throw core::OutOfBoundsException("too few new containers");
            }

            cells.push_back(*it);
            ++it;
        }

        else
        {
            auto container = at(index)->shared_from_this();
            cells.push_back(container);
        }
    }

    for (auto cont = begin; cont != end; ++cont)
    {
        (*cont)->attach(elements_.get());

        // Add all existing objects in the containers to the elements
        for (auto obj: *(*cont))
        {
            elements_->notify_add(obj);
        }
    }

    data_ =
        NCube<std::shared_ptr<CONTAINER_TYPE>>(dim(), new_members, cells.begin(), cells.end());
}


}
}

#endif

