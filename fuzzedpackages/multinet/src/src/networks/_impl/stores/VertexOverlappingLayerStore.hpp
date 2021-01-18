/**

 */

#ifndef UU_MNET_DATASTRUCTURES_STORES_VERTEXOVERLAPPINGLAYERSTORE_H_
#define UU_MNET_DATASTRUCTURES_STORES_VERTEXOVERLAPPINGLAYERSTORE_H_

#include <unordered_set>
#include "core/datastructures/containers/LabeledUniquePtrSortedRandomSet.hpp"
#include "core/datastructures/observers/Subject.hpp"

namespace uu {
namespace net {

template <typename Graph>
class VertexOverlappingLayerStore :
    public core::LabeledUniquePtrSortedRandomSet<Graph>,
    public core::Subject<Graph>
{


  private:

    typedef core::LabeledUniquePtrSortedRandomSet<Graph> super;

  public:

    VertexOverlappingLayerStore()
    {}

    using super::add;
    using super::get;
    using super::erase;
    using super::size;
    using core::Subject<Graph>::observers;


    virtual
    Graph *
    add(
        std::unique_ptr<Graph> g
    ) override
    {
        core::assert_not_null(g.get(), "VertexOverlappingLayerStore::add", "g");

        // Notify the observers.
        for (auto obs: observers)
        {
            obs->notify_add(g.get());
        }

        return super::add(std::move(g));
    }


    virtual
    bool
    erase(
        Graph * g
    ) override
    {
        core::assert_not_null(g, "VertexOverlappingLayerStore::erase", "g");


        // Notify the observers.
        for (auto obs: observers)
        {
            obs->notify_erase(g);
        }

        return super::erase(g);
    }

    /* @todo */
    void
    erase(
        const Vertex* v
    )
    {
        (void)v; // param not used

        for (auto g=super::begin(); g!=super::end(); ++g)
        {
            (*g)->vertices()->erase(v);
        }
    }

    std::string
    summary(
    ) const
    {
        size_t s = size();
        return std::to_string(s) + (s==1?" layer":" layers");
    }

};

}
}

#endif

