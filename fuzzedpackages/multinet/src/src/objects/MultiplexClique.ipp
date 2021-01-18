#include <utility>
#include <string>
#include <memory>

namespace uu {
namespace net {

template <typename M>
MultiplexClique<M>::
MultiplexClique() {}


template <typename M>
MultiplexClique<M>::
MultiplexClique(const std::unordered_set<const Vertex*>& actors, const std::unordered_set<const typename M::layer_type*>& layers) :
    actors(actors.begin(),actors.end()), layers(layers.begin(),layers.end()) {}

template <typename M>
bool
MultiplexClique<M>::
operator==(
    const MultiplexClique<M>& comp
) const
{
    if (actors.size() != comp.actors.size() || layers.size() != comp.layers.size())
    {
        return false;
    }

    std::set<const Vertex*>::iterator it1 = actors.begin();
    std::set<const Vertex*>::iterator it2 = comp.actors.begin();

    for (size_t i = 0; i<actors.size(); i++)
    {
        if ((*it1)!=(*it2))
        {
            return false;
        }

        ++it1;
        ++it2;
    }

    typename std::set<const typename M::layer_type*>::iterator itl1 = layers.begin();
    typename std::set<const typename M::layer_type*>::iterator itl2 = comp.layers.begin();

    for (size_t i = 0; i<layers.size(); i++)
    {
        if ((*itl1)!=(*itl2))
        {
            return false;
        }

        ++itl1;
        ++itl2;
    }

    return true;
}

template <typename M>
bool
MultiplexClique<M>::
operator!=(
    const MultiplexClique<M>& comp
) const
{
    if (actors.size() != comp.actors.size() || layers.size() != comp.layers.size())
    {
        return true;
    }

    std::set<const Vertex*>::iterator it1 = actors.begin();
    std::set<const Vertex*>::iterator it2 = comp.actors.begin();

    for (size_t i = 0; i<actors.size(); i++)
    {
        if ((*it1)!=(*it2))
        {
            return true;
        }

        ++it1;
        ++it2;
    }

    typename std::set<const typename M::layer_type*>::iterator itl1 = layers.begin();
    typename std::set<const typename M::layer_type*>::iterator itl2 = comp.layers.begin();

    for (size_t i = 0; i<layers.size(); i++)
    {
        if ((*itl1)!=(*itl2))
        {
            return true;
        }

        ++itl1;
        ++itl2;
    }

    return false;
}

template <typename M>
bool
MultiplexClique<M>::
operator<(const MultiplexClique<M>& comp) const
{
    if (actors.size() != comp.actors.size())
    {
        return actors.size() < comp.actors.size();
    }

    if (layers.size() != comp.layers.size())
    {
        return layers.size() < comp.layers.size();
    }

    std::set<const Vertex*>::iterator it1 = actors.begin();
    std::set<const Vertex*>::iterator it2 = comp.actors.begin();

    for (size_t i = 0; i<actors.size(); i++)
    {
        if ((*it1)<(*it2))
        {
            return true;
        }

        if ((*it1)>(*it2))
        {
            return false;
        }

        ++it1;
        ++it2;
    }

    typename std::set<const typename M::layer_type*>::iterator itl1 = layers.begin();
    typename std::set<const typename M::layer_type*>::iterator itl2 = comp.layers.begin();

    for (size_t i = 0; i<layers.size(); i++)
    {
        if ((*itl1)<(*itl2))
        {
            return true;
        }

        if ((*itl1)>(*itl2))
        {
            return false;
        }

        ++itl1;
        ++itl2;
    }

    return false;
}

template <typename M>
bool
MultiplexClique<M>::
operator>(const MultiplexClique<M>& comp) const
{
    if (actors.size() != comp.actors.size())
    {
        return actors.size() > comp.actors.size();
    }

    if (layers.size() != comp.layers.size())
    {
        return layers.size() > comp.layers.size();
    }

    std::set<const Vertex*>::iterator it1 = actors.begin();
    std::set<const Vertex*>::iterator it2 = comp.actors.begin();

    for (size_t i = 0; i<actors.size(); i++)
    {
        if ((*it1)>(*it2))
        {
            return true;
        }

        if ((*it1)<(*it2))
        {
            return false;
        }

        ++it1;
        ++it2;
    }

    typename std::set<const typename M::layer_type*>::iterator itl1 = layers.begin();
    typename std::set<const typename M::layer_type*>::iterator itl2 = comp.layers.begin();

    for (size_t i = 0; i<layers.size(); i++)
    {
        if ((*itl1)>(*itl2))
        {
            return true;
        }

        if ((*itl1)<(*itl2))
        {
            return false;
        }

        ++itl1;
        ++itl2;
    }

    return false;
}

template <typename M>
std::string
MultiplexClique<M>::
to_string()
{
    std::ostringstream ss;
    ss << "{ ";

    for (const Vertex* actor: actors)
    {
        ss << actor->name << " ";
    }

    ss << "} + ";
    ss << "[ ";

    for (const typename M::layer_type* layer: layers)
    {
        ss << layer->name << " ";
    }

    ss << "]";
    return ss.str();
}

/*
template <typename M>
std::shared_ptr<MultiplexClique<M>>
                             MultiplexClique<M>::
                             create(
                             )
{
auto result = std::make_shared<MultiplexClique<M>>();
return result;
}

template <typename M>
std::shared_ptr<MultiplexClique<M>>
                             MultiplexClique<M>::
                             create(
                                 const std::unordered_set<const Vertex*>& actors,
                                 const std::unordered_set<const typename M::layer_type*>& layers
                             )
{
auto result = std::make_shared<MultiplexClique<M>>(actors, layers);
return result;
}
*/

}
}
