namespace uu {
namespace net {


template <typename M>
std::unique_ptr<CommunityStructure<VertexLayerCommunity<const typename M::layer_type>>>
mlcpm(
    const M* mnet,
    size_t k,
    size_t m
)
{
    // Step 1: find max-cliques
    auto cliques = find_max_cliques(mnet,k,m);

    if (cliques.size()==0)
    {
        return std::make_unique<CommunityStructure<VertexLayerCommunity<const typename M::layer_type>>>();
    }

    // Step 2: bluid adjacency graph
    auto adjacency = build_max_adjacency_graph(cliques,k,m);
    // @todo use a Network?

    /*std::cout << "ADJACENCY" << std::endl;
    for (auto pair: adjacency)
    {
        std::cout << "node: " << pair.first->to_string() << std::endl;
        for (auto cl: pair.second)
        {
            std::cout << "  " << cl->to_string() << std::endl;

        }
    }
    std::cout << std::endl;*/

    // Step 3: extract communities
    auto comm = find_max_communities(adjacency,m);

    // Translate the result, so that it is compatible with the other clustering algorithms
    auto result = std::make_unique<CommunityStructure<VertexLayerCommunity<const typename M::layer_type>>>();

    for (auto c: comm)
    {
        result->add(std::move(c->to_community()));
        /*
        for (auto c2: comm)
        {
            if (c<=c2) continue;
            auto il = core::s_intersection(c->get_layers(),c2->get_layers());
            auto ia = core::s_intersection(c->actors(),c2->actors());
            if (il.size() == std::min(c->get_layers().size(),c2->get_layers().size()) &&
                ia.size() == c->actors().size() && ia.size() ==  c2->actors().size() )
            {std::cout << "pair" << std::endl;
            std::cout << c->to_string() << std::endl;
            std::cout << c2->to_string() << std::endl;
            }
        }
         */
    }

    return result;
}

// STEP 1: FIND MAX CLIQUES

template <typename M>
std::unordered_set<std::shared_ptr<MultiplexClique<M>>>
find_max_cliques(
    const M* mnet,
    size_t k,
    size_t m
)
{

    std::unordered_set<std::shared_ptr<MultiplexClique<M>>> result;

    // We first define the stack. "instance" is an element on the stack, with information about the clique we are currently processing.
    struct instance
    {
        // skip is an optimization indicating how many elements can be removed from the stack without changing the final results. Setting skip=1 is the same as not using this.
        int skip;
        // This is the current clique being processed. When the processing is done and the current instance is popped out of the stack, all maximal cliques including A must have been returned.
        std::shared_ptr<MultiplexClique<M>> A;
        // This is the list of nodes that can be used to extend A, so that A_ext is still a clique.
        std::vector<std::pair<const Vertex*,std::unordered_set<const typename M::layer_type*> > > B;
        // This is the list of nodes that can be used to extend A, so that A_ext is still a clique, but that can have already been processed. This list is kept to avoid that the same clique is computed following different orders, e.g., (n1,n2,n3) and (n1,n3,n2)
        std::vector<std::pair<const Vertex*,std::unordered_set<const typename M::layer_type*> > > C;
        // constructor
        instance(int skip,
                 std::shared_ptr<MultiplexClique<M>>& A,
                 std::vector<std::pair<const Vertex*,std::unordered_set<const typename M::layer_type*> > >& B,
                 std::vector<std::pair<const Vertex*,std::unordered_set<const typename M::layer_type*> > >& C)
            : skip(skip), A(A), B(B), C(C) {}
    };
    std::vector<std::shared_ptr<instance> > stack;

    // We now initialize the stack with an empty clique, that the function will try to extend with nodes in B (containing all the nodes)
    // A: empty clique
    std::shared_ptr<MultiplexClique<M>> A = std::make_shared<MultiplexClique<M>>(); //TODO add create method

    for (const typename M::layer_type* layer: *mnet->layers())
    {
        A->layers.insert(layer);
    }

    // B: all actors. For each actor we keep the set of layers where it is present
    std::vector<std::pair<const Vertex*,std::unordered_set<const typename M::layer_type*> > > B;

    for (const Vertex* actor: *mnet->actors())
    {
        std::pair<const Vertex*,std::unordered_set<const typename M::layer_type*> > b(actor,std::unordered_set<const typename M::layer_type*>());

        for (auto layer: *mnet->layers())
        {
            if (layer->vertices()->contains(actor))
            {
                b.second.insert(layer);
            }
        }

        B.push_back(b);
    }

    // C: empty set of actors
    std::vector<std::pair<const Vertex*,std::unordered_set<const typename M::layer_type*> > > C;


    // PROCESSING STARTS HERE!
    stack.push_back(std::shared_ptr<instance>(new instance(1,A,B,C)));

    while (!stack.empty())
    {

        // fetch clique to be processed
        std::shared_ptr<instance> inst = stack.back();

        // If no elements are left to extend the clique: remove from stack
        if (inst->B.size()==0)
        {
            for (int i=0; i< inst->skip; i++)
            {
                stack.pop_back();
            }

            continue;
        }

        // Extend the current clique with one of the elements from B
        auto b = inst->B.begin();

        // new clique: A_ext
        std::unordered_set<const typename M::layer_type*> new_layers = core::s_intersection(inst->A->layers,b->second);
        std::unordered_set<const Vertex*> new_actors;

        for (auto actor: inst->A->actors)
        {
            new_actors.insert(actor);
        }

        new_actors.insert(b->first);
        std::shared_ptr<MultiplexClique<M>> A_ext = std::make_shared<MultiplexClique<M>>(new_actors,new_layers);

        // new set of actors that can be used to extend the clique at the next iteration: B_ext
        bool reset_skip = false;
        std::vector<std::pair<const Vertex*,std::unordered_set<const typename M::layer_type*> > > B_ext;
        auto q = b;
        //std::cout << "Extend " << inst->A->to_string() << " to " << A_ext->to_string() <<std::endl;

        for (++q; q!=inst->B.end(); ++q)
        {
            std::unordered_set<const typename M::layer_type*> common = core::s_intersection(q->second,core::s_intersection(b->second,neighboring_layers(mnet,b->first,q->first))); //TODO encapsulate these two intersections into a single function computing them together, for increased readability...

            if (common.size()>=m)
            {
                std::pair<const Vertex*,std::unordered_set<const typename M::layer_type*> > new_b(q->first,std::unordered_set<const typename M::layer_type*>(common.begin(),common.end()));
                B_ext.push_back(new_b);
            }

            // not all layers have been used: we will have to process this again
            if (common.size()<q->second.size())
            {
                reset_skip=true;
            }
        }

        // new set of actors that can be used to extend the clique at the next iteration and have been already used before: C_ext
        std::vector<std::pair<const Vertex*,std::unordered_set<const typename M::layer_type*> > > C_ext;

        for (auto c = inst->C.begin(); c!=inst->C.end(); ++c)
        {
            std::unordered_set<const typename M::layer_type*> common = core::s_intersection(c->second,core::s_intersection(b->second,neighboring_layers(mnet,b->first,c->first)));

            if (common.size()>=m)
            {
                std::pair<const Vertex*,std::unordered_set<const typename M::layer_type*> > new_c(c->first,std::unordered_set<const typename M::layer_type*>(common.begin(),common.end()));
                C_ext.push_back(new_c);
            }
        }

        // move b from B to C
        inst->C.push_back(*b);
        inst->B.erase(b);

        // check if A_ext can lead to new results, or if it is a new result (that is, maximal) itself

        // impossible to grow the current clique to have at least k actors: do not stack it.
        if (A_ext->actors.size()+B_ext.size()<k)
        {
            continue;
        }

        // the current clique has less than m layers: do not stack it.
        if (A_ext->layers.size()<m)
        {
            continue;
        }

        // Check if A_ext should be returned: if
        // (1) it has at least k actors
        // (2) it cannot be extended by any new actor maintaining the current number of layers
        // (3) it cannot be extended by an already processed actors maintaining the current number of layers (otherwise, this sub-clique has already been processed previously by the algorithm)
        //TODO: only max_layers, or is it the intersection that matters? For efficiency reasons.

        bool can_extend_on_B = false;
        bool can_extend_on_C = false;

        size_t max_layers_C = A_ext->layers.size();

        for (auto c: C_ext)
        {
            size_t common_layers = core::s_intersection(A_ext->layers,c.second).size();

            if (common_layers==max_layers_C)
            {
                can_extend_on_C=true;
                continue;
            }
        }

        size_t max_layers_B = inst->A->layers.size();

        for (auto b: B_ext)
        {
            size_t common_layers = core::s_intersection(A_ext->layers,b.second).size();

            if (common_layers==max_layers_B)
            {
                can_extend_on_B=true;
                continue;
            }
        }

        if (A_ext->actors.size()>=k &&
                !can_extend_on_B &&
                !can_extend_on_C)
        {
            result.insert(A_ext);
        }

        // put the new clique under processing
        stack.push_back(std::shared_ptr<instance>(new instance(reset_skip?1:inst->skip+1,A_ext,B_ext,C_ext)));

    }

    return result;
}


template <typename M>
std::unordered_set<const typename M::layer_type*>
neighboring_layers(
    const M* mnet,
    const Vertex* actor1,
    const Vertex* actor2
)
{
    std::unordered_set<const typename M::layer_type*> result;

    for (auto layer: *mnet->layers())
    {
        if (layer->edges()->get(actor1,actor2))
        {
            result.insert(layer);
        }
    }

    return result;
}




template <typename M>
std::map<std::shared_ptr<MultiplexClique<M>>,std::unordered_set<std::shared_ptr<MultiplexClique<M>>> >
        build_max_adjacency_graph(const std::unordered_set<std::shared_ptr<MultiplexClique<M>>>& C, size_t k, size_t m)
{
    std::map<std::shared_ptr<MultiplexClique<M>>,std::unordered_set<std::shared_ptr<MultiplexClique<M>>> > result;

    for (std::shared_ptr<MultiplexClique<M>> c1: C)
    {
        result[c1];

        for (std::shared_ptr<MultiplexClique<M>> c2: C)
        {
            if (c1<=c2)
            {
                continue;
            }

            size_t common_actors = core::s_intersection(c1->actors,c2->actors).size();
            size_t common_layers = core::s_intersection(c1->layers,c2->layers).size();

            if (common_actors>=k-1 && common_layers>=m)
            {
                result[c1].insert(c2);
                result[c2].insert(c1);
            }
        }
    }

    return result;
}



template <typename M>
std::unordered_set<std::shared_ptr<MLCPMCommunity<M>>>
find_max_communities(
    //const M* mnet,
    const std::map<std::shared_ptr<MultiplexClique<M>>,std::unordered_set<std::shared_ptr<MultiplexClique<M>>> >& adjacency,
    size_t m)
{

    // result: empty set of communities
    std::unordered_set<std::shared_ptr<MLCPMCommunity<M>>> result;

    // A: empty community on all layers
    /*
    std::shared_ptr<MLCPMCommunity<M>> A = std::make_shared<MLCPMCommunity<M>>();

    for (const typename M::layer_type* layer: *mnet->layers())
    {
        A->add_layer(layer);
    }*/

    std::unordered_set<std::shared_ptr<MultiplexClique<M>>> AlreadySeen;

    for (auto clique_pair : adjacency)
    {
        std::shared_ptr<MLCPMCommunity<M>> A = std::make_shared<MLCPMCommunity<M>>();
        A->add_clique(clique_pair.first);

        for (const typename M::layer_type* layer: clique_pair.first->layers)
        {
            A->add_layer(layer);
        }

        std::set<std::shared_ptr<MultiplexClique<M>>> Candidates(clique_pair.second.begin(),clique_pair.second.end());
        layer_sets<M> empty;
        find_max_communities2(adjacency,A,Candidates,AlreadySeen,empty,m,result);
        AlreadySeen.insert(clique_pair.first);
    }

    return result;
}



template <typename M>
void
find_max_communities(
    const std::map<std::shared_ptr<MultiplexClique<M>>,std::unordered_set<std::shared_ptr<MultiplexClique<M>>> >& adjacency,
    std::shared_ptr<MLCPMCommunity<M>>& A,
    std::vector<std::shared_ptr<MultiplexClique<M>>> Candidates,
    std::unordered_set<std::shared_ptr<MultiplexClique<M>>>& processedCliques,
    layer_sets<M>& processedLayerCombinations,
    size_t m,
    std::unordered_set<std::shared_ptr<MLCPMCommunity<M>>>& result
)
{

    //std::cout << "Start from " << A->to_string() << std::endl;

    std::vector<std::shared_ptr<MultiplexClique<M>>> stack;

    // EXPAND
    while (Candidates.size()!=0)
    {
        std::shared_ptr<MultiplexClique<M>> c = Candidates.back();
        Candidates.pop_back();
        std::unordered_set<const typename M::layer_type*> i = core::s_intersection(A->get_layers(),c->layers);

        if (i.size()==A->get_layers().size())
        {
            if (processedCliques.count(c)>0)
            {
                // all communities starting from c have already been returned.
                return;
            }

            A->cliques.insert(c);

            for (auto j: adjacency.at(c))
            {
                if (A->cliques.count(j)>0)
                {
                    // clique already present in the community
                    continue;
                }

                else
                {
                    Candidates.push_back(j);
                }
            }
        }

        else if (i.size()>=m)
        {
            if (processedLayerCombinations.count(c->layers)==0)
            {
                stack.push_back(c);
            }
        }

    }

    result.insert(A);

    processedLayerCombinations.insert(A->get_layers());

    layer_sets<M> candidate_layer_combinations;

    for (std::shared_ptr<MultiplexClique<M>> c: stack)
    {
        std::unordered_set<const typename M::layer_type*> s = core::s_intersection(A->get_layers(),c->layers);
        std::set<const typename M::layer_type*> to_be_processed(s.begin(), s.end()); // FIXME
        candidate_layer_combinations.insert(to_be_processed);
    }

    for (std::set<const typename M::layer_type*> layers: candidate_layer_combinations)
    {
        std::shared_ptr<MLCPMCommunity<M>> comm = MLCPMCommunity<M>::create();
        comm->cliques.insert(A->cliques.begin(),A->cliques.end());
        comm->layers.insert(layers.begin(),layers.end());

    }

}


/** Expand community A to a maximum, and recursively do the same when the number of layers decreases. */
template <typename M>
void
find_max_communities2(
    const std::map<std::shared_ptr<MultiplexClique<M>>,std::unordered_set<std::shared_ptr<MultiplexClique<M>>> >& adjacency,
    std::shared_ptr<MLCPMCommunity<M>>& A,
    std::set<std::shared_ptr<MultiplexClique<M>>>& Candidates,
    std::unordered_set<std::shared_ptr<MultiplexClique<M>>>& processedCliques,
    layer_sets<M>& processedLayerCombinations,
    size_t m,
    std::unordered_set<std::shared_ptr<MLCPMCommunity<M>>>& result
)
{

    std::vector<std::shared_ptr<MultiplexClique<M>>> stack;


    //std::cout << "EXPAND COMM: " << A->to_string() << std::endl;

    // EXPAND
    while (Candidates.size()!=0)
    {
        auto element_it = Candidates.begin();
        std::shared_ptr<MultiplexClique<M>> c = *element_it;
        Candidates.erase(element_it);

        std::unordered_set<const typename M::layer_type*> i = core::s_intersection(A->get_layers(),c->layers);


        if (i.size()==A->get_layers().size())
        {
            if (processedCliques.count(c)>0)
            {
                // this community has already been expanded from clique c
                //std::cout << "already expanded with " << c->to_string() << std::endl;
                return;
            }

            //std::cout << " + " << c->to_string() << std::endl;
            A->cliques.insert(c);

            for (auto j: adjacency.at(c))
            {
                if (A->cliques.count(j)>0)
                {
                    // std::cout << "clique already in community " << c->to_string() << std::endl;
                    // clique already present in the community
                    continue;
                }

                else
                {
                    Candidates.insert(j);
                }
            }
        }

        else if (i.size()>=m)
        {

            // new combination of layers to expand
            std::shared_ptr<MLCPMCommunity<M>> comm = MLCPMCommunity<M>::create();
            comm->cliques.insert(A->cliques.begin(),A->cliques.end());
            comm->layers.insert(i.begin(),i.end());

            //std::cout << " new layers? " << c->to_string() << std::endl;

            /*for (auto l: i)
                std::cout << l->name << " ";
            std::cout << std::endl;*/

            if (processedLayerCombinations.find(comm->get_layers()) != processedLayerCombinations.end())
            {
                // this has already been expanded earlier
                //std::cout << " no: already found" << std::endl;
                continue;
            }

            std::set<std::shared_ptr<MultiplexClique<M>>> newCandidates(Candidates.begin(), Candidates.end());
            newCandidates.insert(c);
            // @todo inefficient, pass newCandidates by reference?
            find_max_communities2(adjacency,comm,newCandidates,processedCliques,processedLayerCombinations,m,result);
            processedLayerCombinations.insert(comm->get_layers());
        }

    }

    //std::cout << "OUT: " << A->to_string() << std::endl;

    result.insert(A);

    /*
     processedLayerCombinations.insert(A->get_layers());

     layer_sets<M> candidate_layer_combinations;

     for (std::shared_ptr<MultiplexClique<M>> c: stack)
     {
     std::unordered_set<const typename M::layer_type*> s = core::s_intersection(A->get_layers(),c->layers);
     std::set<const typename M::layer_type*> to_be_processed(s.begin(), s.end()); // FIXME
     candidate_layer_combinations.insert(to_be_processed);
     }

     for (std::set<const typename M::layer_type*> layers: candidate_layer_combinations)
     {
     std::shared_ptr<MLCPMCommunity<M>> comm = MLCPMCommunity<M>::create();
     comm->cliques.insert(A->cliques.begin(),A->cliques.end());
     comm->layers.insert(layers.begin(),layers.end());
     find_max_communities(adjacency,comm,stack,processedCliques,processedLayerCombinations,m,result);
     }
     */
}


}
}

