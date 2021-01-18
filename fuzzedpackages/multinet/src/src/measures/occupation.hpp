#ifndef UU_MEASURES_OCCUPATION_H_
#define UU_MEASURES_OCCUPATION_H_

#include "objects/Vertex.hpp"
#include "objects/EdgeMode.hpp"
#include "randomwalks.hpp"

namespace mlnet {

template <typename M>
hash_map<ActorSharedPtr, int>
occupation(
    const M* mnet, double teleportation, matrix<double> transitions, int num_steps);

/**********************************************************************/
/** Clustering Coefficient ********************************************/
/**********************************************************************/

template <typename M>
double
cc(const M* mnet, const NodeSharedPtr& node);

/**********************************************************************/
/** Neighborhood ******************************************************/
/**********************************************************************/


template <typename M, typename L>
ActorListSharedPtr
neighbors(
    const M* mnet,
    const Vertex* actor,
    LayerIterator first,
    LayerIterator last,
    EdgeMode mode);

template <typename M, typename L>
ActorListSharedPtr
neighbors(
    const M* mnet,
    const Vertex* actor,
    const L* layer,
    EdgeMode mode);


template <typename M, typename L>
ActorListSharedPtr
xneighbors(
    const M* mnet,
    const Vertex* actor,
    LayerIterator first,
    LayerIterator last,
    EdgeMode mode);

template <typename M, typename L>
ActorListSharedPtr
xneighbors(
    const M* mnet,
    const Vertex* actor,
    const L* layer,
    EdgeMode mode);


template <typename M, typename L>
double
connective_redundancy(
    const M* mnet,
    const Vertex* actor,
    LayerIterator first,
    LayerIterator last,
    EdgeMode mode);

/**********************************************************************/
/** Layer relevance *************************************************/
/**********************************************************************/

template <typename M, typename L>
double
relevance(
    const M* mnet,
    const Vertex* actor,
    LayerIterator first,
    LayerIterator last,
    EdgeMode mode);


template <typename M, typename L>
double
relevance(
    const M* mnet,
    const Vertex* actor,
    const L* layer,
    EdgeMode mode);


template <typename M, typename L>
double
xrelevance(
    const M* mnet,
    const Vertex* actor,
    LayerIterator first,
    LayerIterator last,
    EdgeMode mode);

template <typename M, typename L>
double
xrelevance(
    const M* mnet,
    const Vertex* actor,
    const L* layer,
    EdgeMode mode);

/**********************************************************************/
/** Layer comparison ************************************************/
/**********************************************************************/


template <typename M, typename L>
double
pearson_degree(
    const M* mnet, const LayerSharedPtr& layer1, const LayerSharedPtr& layer2, EdgeMode mode);

template <typename M, typename L>
double
rho_degree(
    const M* mnet, const LayerSharedPtr& layer1, const LayerSharedPtr& layer2, EdgeMode mode);

// some utility functions
property_matrix<ActorSharedPtr,LayerSharedPtr,bool>
actor_existence_property_matrix(const M* mnet);
property_matrix<dyad,LayerSharedPtr,bool>
edge_existence_property_matrix(const M* mnet);
property_matrix<triad,LayerSharedPtr,bool>
triangle_existence_property_matrix(const M* mnet);
property_matrix<ActorSharedPtr,LayerSharedPtr,double>
actor_degree_property_matrix(const M* mnet, EdgeMode mode);
property_matrix<ActorSharedPtr,LayerSharedPtr,double>
actor_cc_property_matrix(const M* mnet);

// FROM HERE, PORTING NOT COMPLETED YET

/**********************************************************************/
/** Distances *********************************************************/
/**********************************************************************/

template <typename M>
hash_map<ActorSharedPtr,std::set<path_length> >
pareto_distance(const M* mnet, const Vertex* from);
/*
std::map<actor_id,std::set<Path> > pareto_distance_all_paths(const M* mnet, actor_id vertex);
*/

/**********************************************************************/
/** Betweenness *******************************************************/
/**********************************************************************/

/*
std::map<actor_id,long> pareto_betweenness(const M* mnet);

std::map<global_edge_id, long> pareto_edge_betweenness(const M* mnet);

int check_dominance(const Distance& d1, const Distance& d2);

int check_dominance(const Path& p1, const Path& p2);

*/


}
}

#endif
