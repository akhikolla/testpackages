#ifndef UU_OBJECTS_DYAD_H_
#define UU_OBJECTS_DYAD_H_

#include <set>
#include "objects/Vertex.hpp"

namespace uu {
namespace net {

/**
 * A pair of vertices
 */
class
    Dyad :
    private std::set<const Vertex*>
{
  private:
    typedef std::set<const Vertex*> super;
  public:

    /** Constructor. */
    Dyad(
        const Vertex* v1,
        const Vertex* v2
    );

    bool
    operator==(
        const Dyad& comp
    ) const;

    std::set<const Vertex*>::const_iterator
    begin() const;

    std::set<const Vertex*>::const_iterator
    end() const;

    std::set<const Vertex*>::const_iterator
    find(const Vertex*& val) const;


    /** Output function, presenting a complete description of the dyad. */
    std::string
    to_string(
    ) const;

};

std::ostream&
operator<<(std::ostream& os, const Dyad& d);

}
}

// TODO test collisions
namespace std {
template <>
struct hash<uu::net::Dyad>
{
    size_t
    operator()(const uu::net::Dyad& d) const
    {
        size_t seed = 0;

        for (auto v: d)
        {
            // same effect of boost hash_combine
            seed ^= hash<const uu::net::Vertex*>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }

        return seed;
    }
};


}
#endif
