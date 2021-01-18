#ifndef UU_OBJECTS_TRIAD_H_
#define UU_OBJECTS_TRIAD_H_

#include <set>
#include "objects/Vertex.hpp"

namespace uu {
namespace net {

/**
 * Three vertices
 */
class
    Triad :
    private std::set<const Vertex*>
{
  private:
    typedef std::set<const Vertex*> super;

  public:

    /** Constructor. */
    Triad(
        const Vertex* v1,
        const Vertex* v2,
        const Vertex* v3
    );

    bool
    operator==(
        const Triad& comp
    ) const;

    std::set<const Vertex*>::const_iterator
    begin() const;

    std::set<const Vertex*>::const_iterator
    end() const;

    std::set<const Vertex*>::const_iterator
    find(const Vertex*& val) const;

    /** Output function, presenting a complete description of the triad. */
    std::string
    to_string(
    ) const;

};

std::ostream&
operator<<(std::ostream& os, const Triad& d);

}
}

// TODO test collisions
namespace std {
template <>
struct hash<uu::net::Triad>
{
    size_t
    operator()(const uu::net::Triad& t) const
    {
        size_t seed = 0;

        for (auto v: t)
        {
            seed ^= hash<const uu::net::Vertex*>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }

        return seed;
    }
};
}
#endif
