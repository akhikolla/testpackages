#include "objects/Edge.hpp"
#include "core/exceptions/assert_not_null.hpp"

namespace uu {
namespace net {


Edge::
Edge(
    const Vertex* v1,
    const Vertex* v2,
    EdgeDir dir
) :
    v1(v1),
    v2(v2),
    dir(dir)
{
    assert_not_null(v1, "Edge::constructor", "v1");
    assert_not_null(v2, "Edge::constructor", "v2");
}


std::shared_ptr<Edge>
Edge::
create(
    const Vertex* v1,
    const Vertex* v2,
    EdgeDir dir
)
{
    return std::make_shared<Edge>(v1,v2,dir);
}



std::string
Edge::
to_string(
) const
{
    switch (dir)
    {
    case EdgeDir::DIRECTED:
        return "(" + v1->to_string() + " -> " + v2->to_string() + ")";

    case EdgeDir::UNDIRECTED:
        return "(" + v1->to_string() + " -- " + v2->to_string() + ")";
    }

    return "()"; // never gets here, added to avoid CRAN warning
}


std::ostream&
operator<<(
    std::ostream& os,
    const Edge& e
)
{
    os << e.to_string();
    return os;
}


}
}
