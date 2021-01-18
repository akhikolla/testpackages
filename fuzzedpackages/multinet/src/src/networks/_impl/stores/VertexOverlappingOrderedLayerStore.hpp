/**

 */

#ifndef UU_MNET_DATASTRUCTURES_STORES_VERTEXOVERLAPPINGORDEREDLAYERSTORE_H_
#define UU_MNET_DATASTRUCTURES_STORES_VERTEXOVERLAPPINGORDEREDLAYERSTORE_H_

#include <vector>
#include "core/datastructures/containers/LabeledUniquePtrSortedRandomSet.hpp"
#include "core/datastructures/observers/Subject.hpp"

namespace uu {
namespace net {

template <typename Graph>
class VertexOverlappingOrderedLayerStore :
    public core::Subject<Graph>
{


  private:

    std::vector<std::unique_ptr<Graph>> data;

  public:

    VertexOverlappingOrderedLayerStore()
    {}

    Graph*
    get(const std::string& name)
    {
        for (auto&& g: data)
        {
            if (g->name == name)
            {
                return g.get();
            }
        }

        return nullptr;
    }


    Graph*
    get(const std::string& name) const
    {
        for (auto&& g: data)
        {
            if (g->name == name)
            {
                return g.get();
            }
        }

        return nullptr;
    }

    Graph*
    at(size_t pos)
    {
        return data.at(pos).get();
    }

    const Graph*
    at(size_t pos) const
    {
        return data.at(pos).get();
    }

    size_t
    size() const
    {
        return data.size();
    }


    class
        iterator
    {
        typedef std::forward_iterator_tag iterator_category;

      public:

        iterator();

        /** Returns an iterator pointing at the input object */
        iterator(
            typename std::vector<std::unique_ptr<Graph>>::iterator it
        )
        {
            current = it;
        }

        /** Return the object pointed by this iterator */
        Graph*
        operator*(
        )
        {
            return current->get();
        }

        /** Moves the iterator to the next object in the collection (prefix) */
        iterator
        operator++(
        )
        {
            return ++current;
        }

        /** Moves the iterator to the next object in the collection (postfix) */
        iterator
        operator++(
            int
        )
        {
            return current++;
        }

        /** Checks if this iterator equals the input one */
        bool
        operator==(
            const typename VertexOverlappingOrderedLayerStore<Graph>::iterator& rhs
        )
        {
            return rhs.current = current;
        }

        /** Checks if this iterator differs from the input one */
        bool
        operator!=(
            const typename VertexOverlappingOrderedLayerStore<Graph>::iterator& rhs
        )
        {
            return rhs.current != current;
        }

      private:

        /** Entry currently pointed to by this iterator */
        typename std::vector<std::unique_ptr<Graph>>::iterator current;

    };

    class
        const_iterator
    {
        typedef std::forward_iterator_tag iterator_category;

      public:

        const_iterator();

        /** Returns an iterator pointing at the input object */
        const_iterator(
            typename std::vector<std::unique_ptr<Graph>>::const_iterator it
        )
        {
            current = it;
        }

        /** Return the object pointed by this iterator */
        const Graph*
        operator*(
        )
        {
            return current->get();
        }

        /** Moves the iterator to the next object in the collection (prefix) */
        const_iterator
        operator++(
        )
        {
            return ++current;
        }

        /** Moves the iterator to the next object in the collection (postfix) */
        const_iterator
        operator++(
            int
        )
        {
            return current++;
        }

        /** Checks if this iterator equals the input one */
        bool
        operator==(
            const typename VertexOverlappingOrderedLayerStore<Graph>::const_iterator& rhs
        )
        {
            return rhs.current = current;
        }

        /** Checks if this iterator differs from the input one */
        bool
        operator!=(
            const typename VertexOverlappingOrderedLayerStore<Graph>::const_iterator& rhs
        )
        {
            return rhs.current != current;
        }

      private:

        /** Entry currently pointed to by this iterator */
        typename std::vector<std::unique_ptr<Graph>>::const_iterator current;

    };

    VertexOverlappingOrderedLayerStore<Graph>::iterator
    begin()
    {
        return iterator(data.begin());
    }

    VertexOverlappingOrderedLayerStore<Graph>::const_iterator
    begin() const
    {
        return const_iterator(data.begin());
    }

    VertexOverlappingOrderedLayerStore<Graph>::iterator
    end()
    {
        return iterator(data.end());
    }

    VertexOverlappingOrderedLayerStore<Graph>::const_iterator
    end() const
    {
        return const_iterator(data.end());
    }

    Graph *
    push_back(
        std::unique_ptr<Graph> g
    )
    {
        core::assert_not_null(g.get(), "VertexOverlappingLayerStore::add", "g");

        // Notify the observers.
        /*
        for (auto obs: observers)
        {
            obs->notify_add(g.get());
        }
        */
        data.push_back(std::move(g));

        return data.at(data.size()-1).get();
    }


    virtual
    bool
    erase(
        Graph * g
    )
    {
        core::assert_not_null(g, "VertexOverlappingLayerStore::erase", "g");

        size_t pos = 0;

        for (pos = 0; pos < data.size(); pos++)
        {
            if (data.at(pos).get() == g)
            {
                break;
            }
        }

        if (pos == data.size())
        {
            return false;
        }

        else
        {
            data.erase(data.begin() +  pos);
            return true;
        }
    }

    void
    erase(
        const Vertex* v
    )
    {
        (void)v; // param not used

        for (auto&& g: data)
        {
            g->vertices()->erase(v);
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

