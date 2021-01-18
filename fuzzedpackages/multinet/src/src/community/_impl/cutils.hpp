#ifndef CUTILS_H_
#define CUTILS_H_

#include <chrono>
#include <unordered_set>
#include <unordered_map>
#include <numeric>
#include <random>
#include "community/CommunityStructure.hpp"
#include "community/VertexLayerCommunity.hpp"
#include "Eigen/Sparse"
#include "Eigen/Dense"
#include <vector>
#include "community/Community.hpp"
#include "objects/EdgeMode.hpp"
#include "objects/Vertex.hpp"

namespace uu {
namespace net {


class metanet;



class cutils
{

  public:

    /*
    	Use : std::vector<Eigen::MatrixXd> a = ml_network2adj_matrix(std::shared_ptr<M> ptr);
    	Pre : std::shared_ptr<M> is a non empty multilayer network
    	Post: a is a list of Matrixes that is the adjacency matrix of mnet.
    		Each matrix represents one layer

    template <typename M>
    static
    std::vector<Eigen::SparseMatrix<double> >
    ml_network2adj_matrix(
        const std::shared_ptr<M>& mnet
    );*/

    /*
    	Use: MatrixXd m = sparse_sum(X, axis);
    	Pre: X is an initialized Eigen sparse matrix, axis is an integer with value 0 or 1
    	Post: m is a dense 1-dimensional matrix that has the elements of X summed up along the give axis
    		0 for horizontal, 1 for vertical. Since in Eigen, Vector structures are just typedefs for Matrix we return
    		it as a Matrix. Use it like this: m(i, 0) where 0 <= i <= X.rows()
    */
    static Eigen::MatrixXd
    sparse_sum(
        const Eigen::SparseMatrix<double>& X,
        int axis
    );

    /*
    	Use: std::vector<int> u = cutils::unique(std::vector<int> y);
    	Pre: y is a stl vector
    	Post: returns the unique integers in vector y in ascending order
    */
    static std::vector<int>
    unique(std::vector<int> y);

    /*
    	Use: std::vector<int> r = cutils::unique(size, rand);
    	Post: returns a sequence of integers from 0 to < size. If randomize is set to true
    		the contents of x are randomly shuffled
    */
    static std::vector<int>
    range(
        size_t size,
        bool randomize
    );

    /*
    	Use: CommunityStructureSharedPtr c = cutils::nodes2communities(mnet, nodes2cid);
    	Pre: mnet is a mulltilayer network, with at least some connected actors and 1+ layers
    		 nodes2cid is a vector of community assignment for the nodes in mnet, where
    		 nodes2cid.size() == mnet->get_actors->size() * mnet->layers()->get_all->size()
    	Post: returns a community assignment of actors with their real attributes from the multilayer network

    template <typename M, typename G>
    static
    std::shared_ptr<CommunityStructure<IntralayerVertexCommunity<G>>>
    nodes2communities(
        std::shared_ptr<M> mnet,
        std::vector<unsigned int> nodes2cid
    );*/

    /*
    	Use: CommunityStructureSharedPtr c = cutils::actors2communities(mnet, actors2cid);
    	Pre: mnet is a mulltilayer network, with at least some connected actors and 1+ layers
    		 actors2cid is a vector of community assignments for the actors in mnet, where
    		 actors2cid.size() == mnet->get_actors->size()
    	Post: returns a community assignment of actors with their real attributes from the multilayer network

    template <typename M, typename G>
    static
    std::shared_ptr<CommunityStructure<IntralayerVertexCommunity<G>>>
    actors2communities(
        std::shared_ptr<M> mnet,
        std::vector<unsigned int> actors2cid
    );*/

    /*
    	Use: Eigen::SparseMatrix<double> B = cutils::ng_modularity(twom, a, gamma, omega);
    	Pre: twoum gets the value of the degree of the nodes in the multilayer network
    		 a is a vector of adjacency matrices (layers)
    		 gamma is the resolution parameter, for all layers
    		 omega is the inter-layer coupling strength parameter.
    	Post: returns multilayer Girvan-Newman modularity matrix B and the total number of edges in the network * 2, as reference in the twoum parameter. This function is inspired by https://github.com/GenLouvain/GenLouvain/blob/master/HelperFunctions/multicat.m

    static Eigen::SparseMatrix<double>
    ng_modularity(double& twom, std::vector<Eigen::SparseMatrix<double>> a, double gamma, double omega);
    */

    /*
    	Use: Eigen::SparseMatrix<double> A = cutils::supraA(a, eps, use_node_degrees, use_self_loop);
    	Pre: a is a vector of adjacency matrices (layers)
    		eps is the diagonal and off-diagonal weight parameter
    	Post: Returns a symmetric N * L * N * L matrix where N is the number of actors and L is the number layers
    			Diagonal blocks of the matrix are adjacency matrices from a, also known as inter-layer edges. Off-diagonal blocks are the intra-layer edges connecting the same node across the layers. use_node_degrees will set node degree as intra-layer edge weight on the off-diagonal block. use_self_loop sets eps on the main diagonal, which is a self loop for every node.
    */
    static Eigen::SparseMatrix<double>
    supraA(
        const std::vector<Eigen::SparseMatrix<double>>& a,
        double eps,
        bool use_node_degrees,
        bool use_self_loop
    );


    /*
     Use: Eigen::SparseMatrix<double> A = cutils::supraA(a, eps, use_node_degrees, use_self_loop);
     Pre: a is a vector of adjacency matrices (layers)
     eps is the diagonal and off-diagonal weight parameter
     Post: Returns a symmetric N * L * N * L matrix where N is the number of actors and L is the number layers
     Diagonal blocks of the matrix are adjacency matrices from a, also known as inter-layer edges. Off-diagonal blocks are the intra-layer edges connecting the same node across the layers. use_node_degrees will set node degree as intra-layer edge weight on the off-diagonal block. use_self_loop sets eps on the main diagonal, which is a self loop for every node.
     */
    static Eigen::SparseMatrix<double>
    ordered_supraA(
        const std::vector<Eigen::SparseMatrix<double>>& a,
        double eps,
        bool use_node_degrees,
        bool use_self_loop
    );

    /*
    	Use: Eigen::SparseMatrix<double> diag = cutils::block_diag(a, eps);
    	Pre: a is a vector of adjacency matrices (layers)
    		eps is the diagonal and off-diagonal weight parameter
    	Post: Returns a symmetric N * L * N * L matrix where N is the number of actors and L is the number layers. Diagonal blocks are the intra-layuer edges.
    */
    static Eigen::SparseMatrix<double>
    block_diag(
        const std::vector<Eigen::SparseMatrix<double>>& a,
        double eps
    );

};

/*
template <typename M>
std::vector<Eigen::SparseMatrix<double>>
                                      cutils::ml_network2adj_matrix(
                                              const std::shared_ptr<M>& mnet)
{
    size_t L = mnet->layers()->size();
    size_t N = mnet->actors()->size();

    std::vector<Eigen::SparseMatrix<double>> a(L);

    for (auto l: *mnet->layers())
    {
        Eigen::SparseMatrix<double> m = Eigen::SparseMatrix<double> (N, N);

        std::vector<Eigen::Triplet<double>> tlist;

        // @todo tlist.reserve(mnet->get_edges()->size());

        tlist.reserve(l->edges()->size());
        // end @todo

        for (std::shared_ptr<Edge> e: *l->edges())
        {
            int v1_id = mnet->actors()->index_of(e->v1);
            int v2_id = mnet->actors()->index_of(e->v2);

            tlist.push_back(Eigen::Triplet<double>(v1_id, v2_id, 1));

            if (e->directionality == EdgeDir::UNDIRECTED)
            {
                tlist.push_back(Eigen::Triplet<double>(v2_id, v1_id, 1));
            }

        }

        m.setFromTriplets(tlist.begin(), tlist.end());
        a[l->id - 1] = m;
    }

    return a;
}
*/

/*
template <typename M, typename G>
std::shared_ptr<CommunityStructure<IntralayerVertexCommunity<G>>>
cutils::
nodes2communities(
    std::shared_ptr<M> mnet,
    std::vector<unsigned int> nodes2cid
)
{

    size_t L = mnet->layers()->size();
    size_t N = mnet->actors()->size();

    std::unordered_map<long, std::vector<IntralayerVertex<G>> > result;

    for (size_t i = 0; i < L; i++)
    {
        std::shared_ptr<Layer<G>> l = mnet->layers()->at(i);

        for (size_t j = i * N; j < (1 + i) * N; j++)
        {
            std::shared_ptr<Vertex> a = mnet->actors()->at(j - (i * N));

            // @todo check if node exists
            IntralayerVertex<G> iv(a,l);
            result[nodes2cid[j]].push_back(iv);

        }
    }

    std::shared_ptr<CommunityStructure<IntralayerVertexCommunity<G>>> communities;
    communities = CommunityStructure<IntralayerVertexCommunity<G>>::create();

    for (auto pair: result)
    {
        std::shared_ptr<IntralayerVertexCommunity<G>> c =
                    std::make_shared<IntralayerVertexCommunity<G>>();

        for (IntralayerVertex<G> v: pair.second)
        {
            c->insert(v);
        }

        communities->push_back(c);
    }

    return communities;
}


template <typename M, typename G>
std::shared_ptr<CommunityStructure<IntralayerVertexCommunity<G>>>
cutils::actors2communities(std::shared_ptr<M> mnet, std::vector<unsigned int> actors2cid)
{

    std::unordered_map<long, std::set<IntralayerVertex<G>> > result;

    for (size_t i = 0; i < actors2cid.size(); i++)
    {
        std::shared_ptr<Vertex> a = mnet->actors()->at(i);

        for (auto l : *mnet->layers())
        {
            //NodeSharedPtr n = mnet->get_node(a, l);
            if (l->vertices()->contains(a))
            {
                result[actors2cid[i]].insert(IntralayerVertex<G>(a,l));
            }
        }
    }

    std::shared_ptr<CommunityStructure<IntralayerVertexCommunity<G>>> communities;
    communities = CommunityStructure<IntralayerVertexCommunity<G>>::create();

    for (auto pair: result)
    {
        std::shared_ptr<VertexLayerCommunity<G>> c = IntralayerVertexCommunity<G>::create();

        for (IntralayerVertex<G> node: pair.second)
        {
            c->add(node);
        }

        communities->add_community(c);
    }

    return communities;
}
*/

/*
 Use: Eigen::SparseMatrix<double> B = cutils::ng_modularity(twom, a, gamma, omega);
 Pre: twoum gets the value of the degree of the nodes in the multilayer network
 a is a vector of adjacency matrices (layers)
 gamma is the resolution parameter, for all layers
 omega is the inter-layer coupling strength parameter.
 Post: returns multilayer Girvan-Newman modularity matrix B and the total number of edges in the network * 2, as reference in the twoum parameter. This function is inspired by https://github.com/GenLouvain/GenLouvain/blob/master/HelperFunctions/multicat.m
 */
Eigen::SparseMatrix<double>
modularity_matrix(
    double& twom,
    const std::vector<Eigen::SparseMatrix<double>>& a,
    double gamma,
    double omega,
    bool ordered = false
);


template <typename M>
std::vector<Eigen::SparseMatrix<double>>
                                      to_adjacency_matrices(
                                              const M* mnet
                                      )
{
    size_t num_layers = mnet->layers()->size();
    size_t num_actors = mnet->actors()->size();

    std::vector<Eigen::SparseMatrix<double>> result(num_layers);

    for (size_t idx=0; idx<mnet->layers()->size(); idx++)
    {
        auto l = mnet->layers()->at(idx);

        result[idx] = Eigen::SparseMatrix<double> (num_actors, num_actors);

        //Eigen::SparseMatrix<double> m = Eigen::SparseMatrix<double> (num_actors, num_actors);

        std::vector<Eigen::Triplet<double>> tlist;

        // @todo tlist.reserve(mnet->get_edges()->size());

        tlist.reserve(l->edges()->size());
        // end @todo

        for (auto e: *l->edges())
        {
            size_t v1_id = mnet->actors()->index_of(e->v1);
            size_t v2_id = mnet->actors()->index_of(e->v2);

            tlist.push_back(Eigen::Triplet<double>(v1_id, v2_id, 1));
            //if (e->directionality == EdgeDir::UNDIRECTED)
            //{
            tlist.push_back(Eigen::Triplet<double>(v2_id, v1_id, 1));
            //}

        }

        result[idx].setFromTriplets(tlist.begin(), tlist.end());
    }

    return result;
}


std::unique_ptr<CommunityStructure<Community<const Vertex*>>>
to_community_structure(
    std::unordered_map<const Vertex*, size_t> membership
);

class glouvain
{

  public:

    /*
     Use : glouvain g;
     g.get_ml_community(mnet, move, gamma, omega, limit);

     Pre:
     mnet is the multilayer network
     move is either "move" that places nodes into communities that give the highest modularity gain or "moverandw" which places the nodes into communities at random proportional to their modularity gain (higher gain, higher chance of being selected into communituy)
     gamma is the resolution parameter
     omega is the inter-layer coupling weight parameter
     limit is the number of modularity scores for the actors to hold in memory at once.
     The algorithm is faster if this number is high but consumes much more memory
     Post:
     nodes in mnet are assigned into communities

     Data invariant:
     Generalized Louvain is a multiplex community detector based on the Louvain community detection
     method for single layer networks. This implementation is based from http://netwiki.amath.unc.edu/GenLouvain/GenLouvain
     */
    template <typename M, typename G>
    std::unique_ptr<CommunityStructure<VertexLayerCommunity<const G>>>
    fit(
        const M* mnet,
        const std::string& m,
        double gamma,
        double omega,
        size_t limit
    );

    /* Map indexes of b to values of a:
     https://stackoverflow.com/questions/5691218/matlab-mapping-values-to-index-of-other-array*/
    std::vector<int>
    mapV2I(
        const std::vector<int>& a,
        const std::vector<int>& b
    ) const;

    /*
     Use: double q = Q(M, y, twoum)
     Pre: M is the modularity matrix
     y is the node partitioning vector
     twoum contains the total amount of edges in network
     Post: Calculate the modularity Q of partitioning vector y in matrix M
     */
    double
    Q(
        const Eigen::SparseMatrix<double>& M,
        const std::vector<int>& y,
        double twoum
    ) const;

    /* Same as Q but does not use the full modularity matrix. Instead it uses the iterative version
     for really large networks */
    double
    Q_handle(
        metanet& meta,
        const std::vector<int>& y,
        double twoum
    );

    /*
     Use: Eigen::SparseMatrix<double> M = metanetwork(B, S2)
     Pre: M is the modularity matrix
     S2 is the collapsed node group assignment
     Post: M is a new modularity network, that contains collapsed nodes of B
     */
    Eigen::SparseMatrix<double>
    metanetwork(
        const Eigen::SparseMatrix<double>& B,
        const std::vector<int>& S2
    ) const;

};

struct unique_group_map
{

    /* map for unique possible moves */

    unique_group_map();
    std::vector<bool> ismember;
    std::vector<int> members;
    typedef std::vector<int>::iterator iterator;


    //implement unique_group_map (quick membership check and insertion of elements, quick iteration over members, unordered
    unique_group_map(
        size_t n
    ) : ismember(std::vector<bool>(n,false)) {}

    bool
    count(int i)
    {
        return ismember[i];
    }

    void
    insert(int i)
    {
        if (!ismember[i])
        {
            ismember[i] = true;
            members.push_back(i);
        }
    }

    unique_group_map::iterator
    begin()
    {
        return members.begin();
    }
    unique_group_map::iterator
    end()
    {
        return members.end();
    }

};

#define NUM_TOL 1e-100

typedef std::unordered_map<int, double> map_type;

//map for unique possible moves
typedef unique_group_map set_type;

typedef std::pair<std::vector<int>,std::vector<double>> move_list;

struct group_index
{

    group_index(
    )
        : n_nodes(0), n_groups(0) {}

    group_index(
        const std::vector<int>& v
    )
    {
        n_nodes = v.size();
        nodes = v;
        nodes_iterator.resize(n_nodes);
        n_groups = *max_element(nodes.begin(),nodes.end()) + 1;
        groups.resize(n_groups);

        for (size_t i = 0; i < n_nodes; i++)
        {
            groups[nodes[i]].push_back(i);
            nodes_iterator[i] = --groups[nodes[i]].end();
        }
    };

    Eigen::MatrixXd
    index(
        int group
    )
    {
        Eigen::MatrixXd r = Eigen::MatrixXd::Zero(groups[group].size(), 1);

        int i = 0;

        for (std::list<int>::iterator it=groups[group].begin(); it != groups[group].end(); it++)
        {
            r(i, 0) = *it;
            i++;
        }

        return r;
    }

    //move node to group
    void
    move(
        int node,
        int group
    )
    {
        //move node by splicing into list for new group
        groups[group].splice(groups[group].end(), groups[nodes[node]],nodes_iterator[node]);
        //update its group asignment
        nodes[node]=group;
    };

    size_t n_nodes;
    size_t n_groups;

    std::vector<int>
    toVector()
    {
        std::vector<int> v (n_nodes);
        std::vector<bool> track_move(n_nodes, true);
        size_t g_n = 0;

        std::list<int>::iterator it;

        for (size_t i = 0; i < n_nodes; i++)
        {
            if (track_move[i])
            {
                for (it=groups[nodes[i]].begin(); it != groups[nodes[i]].end(); it++)
                {
                    v[*it] = g_n;
                    track_move[*it] = false;
                }

                g_n++;
            }
        }

        return v;
    };

    std::vector<std::list<int>> groups; //the index of each node in a group is stored in a linked list

    std::vector<std::list<int>::iterator> nodes_iterator; //stores the position of the node in the list for the group it belongs to

    std::vector<int> nodes; //stores the group a node belongs to

};

set_type
possible_moves(
    group_index & g,
    int node,
    const Eigen::SparseMatrix<double>& mod
);

//calculates changes in modularity for sparse modularity matrix
map_type
mod_change(
    group_index & g,
    const Eigen::SparseMatrix<double>& mod,
    set_type & unique_groups,
    int current_node
);

//find moves that improve modularity
move_list
positive_moves(
    set_type & unique_groups,
    map_type & mod_c
);

//move best move
double
move(
    group_index & g,
    int node,
    const Eigen::SparseMatrix<double>& mod
);

//move to random group with probability proportional to increase in modularity
double
moverandw(
    group_index & g,
    int node,
    const Eigen::SparseMatrix<double>& mod
);


class metanet
{

  public:

    /*
     Data invariant:
     metanet is a datastructure that keeps track of visited nodes in really large networks.
     The idea is to calculate the ng modularity of each node every iteration.
     Works only for undirected and unweighted networks.
     */

    metanet(
        const std::vector<Eigen::SparseMatrix<double>>& a,
        double gamma,
        double omega,
        bool ordered = false
    )
    {

        if (ordered)
        {
            AA = cutils::ordered_supraA(a, omega, false, false);
        }

        else
        {
            AA = cutils::supraA(a, omega, false, false);
        }

        K = supraK(a);
        kvec = cutils::sparse_sum(cutils::supraA(a, 0, false, false), 0);

        this->gamma = gamma;
        this->omega = gamma;
        this->N = a[0].rows();
    }


    Eigen::SparseMatrix<double>
    get(
        int index
    )
    {
        std::list<int> ind = nodes(index);

        for (int node: ind)
        {
            reduce(ng_handle(node));
        }

        return ret();
    }

    void
    assign(
        const std::vector<int>& v
    )
    {
        group = group_index(v);
        mod_reduced = std::vector<double>(group.n_groups, 0);
    }

    std::list<int>
    nodes(
        int index
    ) const
    {
        return group.groups.at(index);
    }

    void
    reduce(
        const Eigen::SparseMatrix<double>& mod
    )
    {
        for (int i = 0; i < mod.outerSize(); i++)
        {
            for (Eigen::SparseMatrix<double>::InnerIterator it(mod, i); it; ++it)
            {
                mod_reduced[group.nodes[it.row()]] += it.value();
            }
        }
    }

    Eigen::SparseMatrix<double>
    ret(
    )
    {
        Eigen::SparseMatrix<double> out(mod_reduced.size(), 1);
        std::vector<Eigen::Triplet<double>> tlist;
        tlist.reserve(mod_reduced.size());

        for (size_t i = 0; i < mod_reduced.size(); i++)
        {
            tlist.push_back(Eigen::Triplet<double>(i, 0, mod_reduced[i]));
        }

        out.setFromTriplets(tlist.begin(), tlist.end());

        for (std::vector<double>::iterator it=mod_reduced.begin(); it!=mod_reduced.end(); ++it)
        {
            *it=0;
        }

        return out;
    }


  private:
    group_index group;

    std::vector<double> mod_reduced = std::vector<double>();

    Eigen::SparseMatrix<double> AA, K;

    double gamma, omega;

    int N;

    Eigen::MatrixXd kvec;


    Eigen::SparseMatrix<double>
    ng_handle(
        int index
    ) const
    {
        Eigen::SparseMatrix<double> Acol = AA.col(index);
        Eigen::SparseMatrix<double> Kcol = K.col(index / N);

        std::vector<Eigen::Triplet<double>> tlist;
        tlist.reserve(Kcol.nonZeros());

        for (int i = 0; i < Kcol.outerSize(); i++)
        {
            for (Eigen::SparseMatrix<double>::InnerIterator it(Kcol, i); it; ++it)
            {
                tlist.push_back(
                    Eigen::Triplet<double>(it.row(), it.col(), it.value() * gamma * kvec(index, 0)));
            }
        }

        Kcol.setFromTriplets(tlist.begin(), tlist.end());
        return Acol - Kcol;
    }

    Eigen::SparseMatrix<double>
    supraK(
        const std::vector<Eigen::SparseMatrix<double>>& a
    )
    {
        Eigen::SparseMatrix<double> K = Eigen::SparseMatrix<double>(a[0].rows() * a.size(), a.size());
        std::vector<Eigen::Triplet<double>> tlist;
        tlist.reserve(a[0].rows() * a.size());

        for (size_t i = 0; i < a.size(); i++)
        {
            Eigen::MatrixXd k = cutils::sparse_sum(a[i], 0);
            k = k.array() / k.array().sum();

            for (int j = 0; j < k.rows(); j++)
            {
                tlist.push_back(Eigen::Triplet<double>(i * k.rows() + j, i, k(j, 0)));
            }
        }

        K.setFromTriplets(tlist.begin(), tlist.end());
        return K;
    }
};



template <typename M, typename G>
std::unique_ptr<CommunityStructure<VertexLayerCommunity<const G>>>
to_community_structure(
    const M* mnet,
    const std::vector<unsigned int>& nodes2cid
)
{

    size_t num_layers = mnet->layers()->size();
    size_t num_actors = mnet->actors()->size();

    // group by community id

    std::unordered_map<unsigned int, std::list<std::pair<const Vertex*, const G*>> > result;

    for (size_t i = 0; i < num_layers; i++)
    {
        auto layer = mnet->layers()->at(i);

        for (size_t j = i * num_actors; j < (1 + i) * num_actors; j++)
        {
            auto actor = mnet->actors()->at(j - (i * num_actors));


            if (!layer->vertices()->contains(actor))
            {
                continue;
            }

            auto iv = std::make_pair(actor, layer);
            result[nodes2cid.at(j)].push_back(iv);

        }
    }

    // build community structure

    auto communities = std::make_unique<CommunityStructure<VertexLayerCommunity<const G>>>();

    for (auto pair: result)
    {
        auto c = std::make_unique<VertexLayerCommunity<const G>>();

        for (auto vertex_layer_pair: pair.second)
        {
            c->add(vertex_layer_pair);
        }

        communities->add(std::move(c));
    }

    return communities;
}


}
}


#endif
