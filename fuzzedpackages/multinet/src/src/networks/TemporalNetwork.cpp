#include "networks/TemporalNetwork.hpp"

namespace uu {
namespace net {

TemporalNetwork::
TemporalNetwork(
    const std::string& name,
    EdgeDir dir,
    bool allows_loops
) : super(name, dir, allows_loops)
{
    auto t_attr = core::Attribute::create(kTIME_ATTR_NAME, core::AttributeType::TIME);

    edges()->attr()->add(std::move(t_attr));

    // this index allows quick access to edges based on temporal criteria
    edges()->attr()->add_index(kTIME_ATTR_NAME);

}

void
TemporalNetwork::
set_time(
    const Edge* e,
    core::Time t
)
{
    edges()->attr()->set_time(e, kTIME_ATTR_NAME, t);
}


core::Value<core::Time>
TemporalNetwork::
get_time(
    const Edge* e
) const
{
    return edges()->attr()->get_time(e, kTIME_ATTR_NAME);
}


core::Value<core::Time>
TemporalNetwork::
get_min_time(
) const
{
    return edges()->attr()->get_min_time(kTIME_ATTR_NAME);
}


core::Value<core::Time>
TemporalNetwork::
get_max_time(
) const
{
    return edges()->attr()->get_max_time(kTIME_ATTR_NAME);
}


bool
TemporalNetwork::
is_temporal(
) const
{
    return true;
}


}
}

