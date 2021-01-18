#ifndef GIRAF_H
#define GIRAF_H

#include <RcppCommon.h>

RCPP_EXPOSED_CLASS(Lattice);
RCPP_EXPOSED_CLASS(Border);
RCPP_EXPOSED_CLASS(Block);
  
#include <RcppArmadillo.h>
#include <Rcpp.h>

#include <vector>
#include <armadillo>

#define BOOST_DISABLE_ASSERTS 1    
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/properties.hpp>
#include <boost/tuple/tuple.hpp>
  
/***************************************/
typedef boost::property< boost::edge_weight_t, double,
                         boost::property<boost::edge_weight2_t, unsigned,
                         boost::property<boost::edge_update_t, double > > >
        EdgeProperties;
/***************************************/

/***************************************/
typedef boost::property< boost::vertex_color_t, unsigned,
                         boost::property<boost::vertex_potential_t, arma::vec,
                         boost::property<boost::vertex_degree_t, unsigned,
                         boost::property<boost::vertex_update_t, double,
                         boost::property<boost::vertex_underlying_t, 
                         std::vector<int> > > > > >
        VertexProperties;
/***************************************/

/***************************************/
typedef boost::adjacency_list<
  boost::listS,
  boost::vecS,
  boost::undirectedS,
  VertexProperties,
  EdgeProperties
> PixelGraph;

typedef boost::property_map <PixelGraph, boost::vertex_color_t>::type VertexMap;
typedef boost::property_map <PixelGraph, boost::vertex_potential_t>::type VertexMap_pot;
typedef boost::property_map <PixelGraph, boost::vertex_degree_t>::type VertexMap_degree;
typedef boost::property_map <PixelGraph, boost::vertex_update_t>::type VertexMap_noise;
typedef boost::property_map <PixelGraph, boost::vertex_underlying_t>::type VertexMap_latent;

typedef boost::property_map <PixelGraph, boost::edge_weight_t>::type EdgeMap;
typedef boost::property_map <PixelGraph, boost::edge_weight2_t>::type EdgeMap_dir;
typedef boost::property_map <PixelGraph, boost::edge_update_t>::type EdgeMap_SW;
/***************************************/
  
#include <lattice.hpp>
#include <recursion.hpp>
#include <sampler.hpp>

#endif
