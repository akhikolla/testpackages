/**
 *
 * History:
 * - 2019.08.01 file created
 */

#ifndef UU_OLAP_CUBE_H_
#define UU_OLAP_CUBE_H_

#include "core/exceptions/ElementNotFoundException.hpp"
#include "core/exceptions/WrongParameterException.hpp"
#include <string>
#include <unordered_map>
#include <vector>

namespace uu {
namespace core {

template <class ENTRY_TYPE>
class Cube
{

  public:

    typedef ENTRY_TYPE entry_type;

  public:

    /**
     * Creates a 0-order cube.
     */
    Cube(
    );

    /**
     * Creates a cube specifying its size.
     */
    Cube(
        const std::vector<size_t>& size
    );


    /**
     * Creates a cube specifying its size and filling it in using the input iterator.
     */
    template <class Iterator>
    Cube(
        const std::vector<size_t>& size,
        Iterator begin,
        Iterator end
    );


    /**
     * Destroys the cube.
     */
    ~Cube()
    {

    }

    /**
     * Resizes the cube.
     * This operation does not change the order of the elements already present inside the cube,
     * but can (1) remove some elements if the number of cells decreases, (2) create empty cells
     * if the number of cells increases, (3) change the index of existing elements.
     */
    void
    resize(
        const std::vector<size_t>& size
    );

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


    /** Returns a const iterator to the first object in the cube */
    typename std::vector<ENTRY_TYPE>::const_iterator
    begin(
    ) const;

    /** Returns a const iterator after the last object in the cube */
    typename std::vector<ENTRY_TYPE>::const_iterator
    end(
    ) const;

    /** Returns an iterator to the first object in the cube */
    typename std::vector<ENTRY_TYPE>::iterator
    begin(
    );

    /** Returns an iterator after the last object in the cube */
    typename std::vector<ENTRY_TYPE>::iterator
    end(
    );

    /**
     * Returns the object at the given index in the cube.
     * @throw ElementNotFoundException if the index is outside the bounds on the cube
     */
    ENTRY_TYPE&
    operator[](
        const std::vector<size_t>& index
    );

    /**
     * Returns the object at the given position in the cube.
     * @throw ElementNotFoundException if the index is outside the bounds on the cube
     */
    const ENTRY_TYPE&
    operator[](
        const std::vector<size_t>& index
    ) const;

    /**
     * Returns the object at the given position in the cube.
     * @throw ElementNotFoundException if the index is outside the bounds on the cube
     */
    ENTRY_TYPE&
    at(
        const std::vector<size_t>& index
    );


    /** Returns the object at the given position in the cube.
     * @throw ElementNotFoundException if the index is outside the bounds on the cube
     */
    const ENTRY_TYPE&
    at(
        const std::vector<size_t>& index
    ) const;


  protected:

    /* Size, for each dimension */
    std::vector<size_t> size_;

    /* Offsets */
    std::vector<size_t> off_;

    /* Data */
    std::vector<ENTRY_TYPE> data_;


  private:

    size_t
    pos(
        const std::vector<size_t>& index
    ) const;


};

/* TEMPLATE CODE */


template <class ENTRY_TYPE>
Cube<ENTRY_TYPE>::
Cube(
)
{
    size_ = {};
    off_ = {};
}


template <class ENTRY_TYPE>
Cube<ENTRY_TYPE>::
Cube(
    const std::vector<size_t>& size
)
{
    resize(size);
}


template <class ENTRY_TYPE>
template <class Iterator>
Cube<ENTRY_TYPE>::
Cube(
    const std::vector<size_t>& size,
    Iterator begin,
    Iterator end
) : data_(begin, end)
{
    resize(size);
}


template <class ENTRY_TYPE>
void
Cube<ENTRY_TYPE>::
resize(
    const std::vector<size_t>& size
)
{
    size_ = size;

    off_.resize(size_.size());

    size_t data_size = 1;

    for (size_t i = 0; i < size_.size(); i++)
    {
        off_[i] = data_size;
        data_size *= size_.at(i);
    }

    data_.resize(data_size);
}


template <class ENTRY_TYPE>
size_t
Cube<ENTRY_TYPE>::
order(
) const
{
    return size_.size();
}


template <class ENTRY_TYPE>
std::vector<size_t>
Cube<ENTRY_TYPE>::
size(
) const
{
    return size_;
}


template <class ENTRY_TYPE>
typename std::vector<ENTRY_TYPE>::iterator
Cube<ENTRY_TYPE>::
begin(
)
{
    return data_.begin();
}

template <class ENTRY_TYPE>
typename std::vector<ENTRY_TYPE>::iterator
Cube<ENTRY_TYPE>::
end(
)
{
    return data_.end();
}

template <class ENTRY_TYPE>
typename std::vector<ENTRY_TYPE>::const_iterator
Cube<ENTRY_TYPE>::
begin(
) const
{
    return data_.begin();
}

template <class ENTRY_TYPE>
typename std::vector<ENTRY_TYPE>::const_iterator
Cube<ENTRY_TYPE>::
end(
) const
{
    return data_.end();
}

template <class ENTRY_TYPE>
ENTRY_TYPE&
Cube<ENTRY_TYPE>::
operator[](
    const std::vector<size_t>& index
)
{

    size_t idx = pos(index);

    return data_[idx];
}

template <class ENTRY_TYPE>
const ENTRY_TYPE&
Cube<ENTRY_TYPE>::
operator[](
    const std::vector<size_t>& index
) const
{

    size_t idx = pos(index);

    return data_[idx];
}


template <class ENTRY_TYPE>
ENTRY_TYPE&
Cube<ENTRY_TYPE>::
at(
    const std::vector<size_t>& index
)
{
    size_t idx = pos(index);

    return data_.at(idx);
}

template <class ENTRY_TYPE>
const ENTRY_TYPE&
Cube<ENTRY_TYPE>::
at(
    const std::vector<size_t>& index
) const
{
    size_t idx = pos(index);

    return data_.at(idx);
}

template <class ENTRY_TYPE>
size_t
Cube<ENTRY_TYPE>::
pos(
    const std::vector<size_t>& index
) const
{

    size_t idx = 0;

    for (size_t i = 0; i < size_.size(); i++)
    {
        idx += index.at(i) * off_.at(i);
    }

    return idx;

}

}
}

#endif

