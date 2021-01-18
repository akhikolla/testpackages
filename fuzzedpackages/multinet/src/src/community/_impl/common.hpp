/* Community utils, common routines used by multilayer community detection algorithms */

#ifndef UU_MNET_COMMUNITY_COMMON_H_
#define UU_MNET_COMMUNITY_COMMON_H_

#include <chrono>
#include <unordered_set>
#include <unordered_map>
#include <numeric>
#include <random>
#include <utility>
#include <memory>
#include "objects/Vertex.hpp"
#include "community/CommunityStructure.hpp"
#include "community/VertexLayerCommunity.hpp"
#include <Eigen/Sparse>
#include <Eigen/Dense>

namespace uu {
namespace net {

/*
*	Use : std::vector<Eigen::MatrixXd> a = ml_network2adj_matrix(std::shared_ptr<M> ptr);
*	@param mnet non empty multilayer network
 *	@return a vector of adjacency matrices; each matrix represents one layer
*/
template <typename M>
std::vector<Eigen::SparseMatrix<double>>
                                      to_adjacency_matrices(
                                              const M* mnet
                                      );

/*
	Use: MatrixXd m = sparse_sum(X, axis);
	Pre: X is an initialized Eigen sparse matrix, axis is an integer with value 0 or 1
	Post: m is a dense 1-dimensional matrix that has the elements of X summed up along the given axis
		0 for horizontal, 1 for vertical. Since in Eigen, Vector structures are just typedefs for Matrix we return
		it as a Matrix. Use it like this: m(i, 0) where 0 <= i <= X.rows()
*/
static Eigen::MatrixXd
sparse_sum(
    const Eigen::SparseMatrix<double>& X,
    int axis
);

/*
	Use: std::vector<int> u = unique(std::vector<int> y);
	Pre: y is a stl vector
	Post: returns the unique integers in vector y in ascending order
*/
static std::vector<int>
unique(
    std::vector<int> y
);

/*
	Use: std::vector<int> r = range(size, rand);
	Post: returns a sequence of integers from 0 to < size. If randomize is set to true
		the contents of x are randomly shuffled
*/
static std::vector<int>
range(
    size_t size,
    bool randomize
);

/*
* Converts a vector with community identifiers for each vertex/layer pair into a community structure.
* @param mnet is a multilayer network, with at least some connected actors and 1+ layers
* @param nodes2cid is a vector of community assignment for the nodes in mnet, where
* nodes2cid.size() == num_actors * num_layers
* @return a community assignment of actors
*/

template <typename M, typename G>
std::unique_ptr<CommunityStructure<VertexLayerCommunity<G>>>
to_community_structure(
    const M* mnet,
    const std::vector<unsigned int>& nodes2cid
);



/*
* Converts a vector with community identifiers for each vertex/layer pair into a community structure.
* @param mnet is a mulltilayer network, with at least some connected actors and 1+ layers
* @param actors2cid is a vector of community assignments for the actors in mnet, where
		 actors2cid.size() == mnet->get_actors->size()
	Post: returns a community assignment of actors with their real attributes from the multilayer network
// * //
 template <typename M, typename G>
 std::unique_ptr<CommunityStructure<VertexLayerCommunity<G>>>
 to_community_structure(
 const M* mnet,
 std::vector<unsigned int> nodes2cid
 );
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
    double omega
);

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
	Use: Eigen::SparseMatrix<double> diag = cutils::block_diag(a, eps);
	Pre: a is a vector of adjacency matrices (layers)
		eps is the diagonal and off-diagonal weight parameter
	Post: Returns a symmetric N * L * N * L matrix where N is the number of actors and L is the number layers. Diagonal blocks are the intra-layer edges.
*/
static Eigen::SparseMatrix<double>
block_diag(
    const std::vector<Eigen::SparseMatrix<double>>& a,
    double eps
);

//};

template <typename M>
std::vector<Eigen::SparseMatrix<double>>
                                      to_adjacency_matrices(
                                              const M* mnet
                                      )
{
    size_t num_layers = mnet->layers()->size();
    size_t num_actors = mnet->vertices()->size();

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
            size_t v1_id = mnet->vertices()->index_of(e->v1);
            size_t v2_id = mnet->vertices()->index_of(e->v2);

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



template <typename M, typename G>
std::unique_ptr<CommunityStructure<VertexLayerCommunity<G>>>
to_community_structure(
    const M* mnet,
    const std::vector<unsigned int>& nodes2cid
)
{
    size_t num_layers = mnet->layers()->size();
    size_t num_actors = mnet->vertices()->size();

    // group by community id

    std::unordered_map<unsigned int, std::list<std::pair<const Vertex*, const G*>> > result;

    for (size_t i = 0; i < num_layers; i++)
    {
        auto layer = mnet->layers()->at(i);

        for (size_t j = i * num_actors; j < (1 + i) * num_actors; j++)
        {
            auto actor = mnet->vertices()->at(j - (i * num_actors));

            if (!layer->vertices()->contains(actor))
            {
                continue;
            }

            auto iv = std::make_pair(actor, layer);
            result[nodes2cid.at(j)].push_back(iv);

        }
    }

    // build community structure

    auto communities = std::make_unique<CommunityStructure<VertexLayerCommunity<G>>>();

    for (auto pair: result)
    {
        auto c = std::make_unique<VertexLayerCommunity<G>>();

        for (auto vertex_layer_pair: pair.second)
        {
            c->add(vertex_layer_pair);
        }

        communities->add(std::move(c));
    }

    return communities;
}

/*

template <typename M, typename G>
std::shared_ptr<CommunityStructure<IntralayerVertexCommunity<G>>>
cutils::actors2communities(std::shared_ptr<M> mnet, std::vector<unsigned int> actors2cid)
{

    std::unordered_map<long, std::set<IntralayerVertex<G>> > result;

    for (size_t i = 0; i < actors2cid.size(); i++)
    {
        std::shared_ptr<Vertex> a = mnet->vertices()->at(i);

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
        std::shared_ptr<IntralayerVertexCommunity<G>> c = IntralayerVertexCommunity<G>::create();

        for (IntralayerVertex<G> node: pair.second)
        {
            c->add(node);
        }

        communities->add_community(c);
    }

    return communities;
}

    */
}
}


#endif
