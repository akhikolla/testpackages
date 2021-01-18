/**
 * History:
 * - 2019.07.21 File created
 */

#ifndef UU_NET_DATASTRUCTURES_OLAP_RESIZE_H_
#define UU_NET_DATASTRUCTURES_OLAP_RESIZE_H_

#include <memory>
#include <string>
#include "networks/_impl/olap/VCube.hpp"

namespace uu {
namespace net {

/*
template <typename C>
sel::EntryIterator<C>
islice(
C* cube,
const std::vector<std::vector<size_t>>& indexes
);


template <typename CONTAINER_TYPE>
std::unique_ptr<CCube<CONTAINER_TYPE>>
vslice(
    CCube<CONTAINER_TYPE>* cube,
    const std::vector<std::vector<size_t>>& indexes
    );

 */

void
resize(
    VCube* cube,
    const std::string& dim,
    const std::string& member
);

// DEFINITIONS
/*
template <typename C>
sel::EntryIterator<C>
islice(
C* cube,
const std::vector<std::vector<size_t>>& indexes
)
{
return sel::EntryIterator<C>(cube, indexes);
}


template <typename CONTAINER_TYPE>
std::unique_ptr<CCube<CONTAINER_TYPE>>
                                vslice(
                                    CCube<CONTAINER_TYPE>* cube,
                                    const std::vector<std::vector<size_t>>& indexes
                                )
{
// @todo check not null

// @todo check bounds

// get dimensions and members from input cube

auto dim_names = cube->dim();

std::vector<std::vector<std::string>> members(dim_names.size());

for (size_t i = 0; i < dim_names.size(); i++)
{
    auto member_names = cube->members(dim_names.at(i));

    for (auto idx: indexes[i])
    {
        members[i].push_back(member_names.at(idx));
    }
}

std::vector<std::vector<size_t>> all;

for (size_t i = 0; i < dim_names.size(); i++)
{
    all.push_back(std::vector<size_t>());

    for (size_t j = 0; j < members[i].size(); j++)
    {
        all[i].push_back(j);
    }
}

sel::IndexIterator a(all);

// insert pointers to all elements

auto it = islice(cube, indexes);
auto slice = std::make_unique<CCube<CONTAINER_TYPE>>(dim_names, members);

return slice;
}
*/

}
}

#endif
