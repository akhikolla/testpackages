#include "graph.h"

//' get_to_from
//'
//' Get one pair of two and from edges and vertices. Main task is to make sure
//' that bi-directed edges ("intermediate_double") correctly return the
//' **different** values of from and to vertices and edges.
//'
//' @noRd
void graph_contract::get_to_from (const edge_map_t &edge_map,
        const std::unordered_set <edge_id_t> &edges,
        const std::vector <vertex_id_t> &two_nbs,
        vertex_id_t &vt_from, vertex_id_t &vt_to,
        edge_id_t &edge_from_id, edge_id_t &edge_to_id)
{
    for (auto edge_id: edges)
    {
        edge_t edge = edge_map.find (edge_id)->second;
        if (edge_from_id == "")
        {
            if (edge.get_from_vertex () == two_nbs [0] &&
                    (edge_to_id == "" || vt_to != two_nbs [0]))
            {
                edge_from_id = edge_id;
                vt_from = two_nbs [0];
            } else if (edge.get_from_vertex () == two_nbs [1] &&
                    (edge_to_id == "" || vt_to != two_nbs [1]))
            {
                edge_from_id = edge_id;
                vt_from = two_nbs [1];
            }
        }
        if (edge_to_id == "")
        {
            if (edge.get_to_vertex () == two_nbs [0] &&
                    (edge_from_id == "" || vt_from != two_nbs [0]))
            {
                edge_to_id = edge_id;
                vt_to = two_nbs [0];
            } else if (edge.get_to_vertex () == two_nbs [1] &&
                    (edge_from_id == "" || vt_from != two_nbs [1]))
            {
                edge_to_id = edge_id;
                vt_to = two_nbs [1];
            }
        }
    }
}

void graph_contract::contract_one_edge (vert2edge_map_t &vert2edge_map,
        vertex_map_t &vertex_map, edge_map_t &edge_map,
        const std::unordered_set <edge_id_t> &edgelist,
        const vertex_id_t vtx_id, const vertex_id_t vt_from,
        const vertex_id_t vt_to,
        const edge_id_t edge_from_id, const edge_id_t edge_to_id,
        const edge_id_t new_edge_id,
        bool has_times)
{
    edge_t edge_from = edge_map.find (edge_from_id)->second,
           edge_to = edge_map.find (edge_to_id)->second;
    double time = -INFINITE_DOUBLE, timew = -INFINITE_DOUBLE,
           d = edge_from.dist + edge_to.dist,
           w = edge_from.weight + edge_to.weight;
    if (has_times)
    {
        time = edge_from.time + edge_to.time;
        timew = edge_from.timew + edge_to.timew;
    }

    std::set <edge_id_t> old_edges,
        old_edges_fr = edge_from.get_old_edges (),
        old_edges_to = edge_to.get_old_edges ();
    for (auto e: old_edges_fr)
        old_edges.insert (e);
    for (auto e: old_edges_to)
        old_edges.insert (e);
    // Only insert edge IDs that are in the original list
    if (edgelist.find (edge_from_id) != edgelist.end ())
        old_edges.insert (edge_from_id);
    if (edgelist.find (edge_to_id) != edgelist.end ())
        old_edges.insert (edge_to_id);

    graph::erase_from_v2e_map (vert2edge_map, vtx_id, edge_from_id);
    graph::erase_from_v2e_map (vert2edge_map, vtx_id, edge_to_id);
    graph::erase_from_v2e_map (vert2edge_map, vt_from, edge_from_id);
    graph::erase_from_v2e_map (vert2edge_map, vt_to, edge_to_id);
    graph::add_to_v2e_map (vert2edge_map, vt_from, new_edge_id);
    graph::add_to_v2e_map (vert2edge_map, vt_to, new_edge_id);

    vertex_t vt = vertex_map [vt_from];
    vt.replace_neighbour (vtx_id, vt_from);
    vertex_map [vt_from] = vt;
    vt = vertex_map [vt_to];
    vt.replace_neighbour (vtx_id, vt_to);
    vertex_map [vt_to] = vt;

    edge_map.erase (edge_from_id);
    edge_map.erase (edge_to_id);
    std::vector <double> wts;
    if (!has_times)
        wts = {d, w};
    else
        wts = {d, w, time, timew};
    edge_t new_edge = edge_t (vt_from, vt_to, wts, new_edge_id, old_edges);
    edge_map.emplace (new_edge_id, new_edge);
}

//' same_hwy_type
//'
//' Determine whether two edges represent the same weight category (type of
//' highway for street networks, for example). Categories are not retained in 
//' converted graphs, but can be discerned by comparing ratios of weighted to
//' non-weighted distances.
//' @noRd
bool graph_contract::same_hwy_type (const edge_map_t &edge_map,
        const edge_id_t &e1, const edge_id_t &e2)
{
    const double tol = 1.0e-6;

    edge_t edge1 = edge_map.find (e1)->second,
           edge2 = edge_map.find (e2)->second;

    return (fabs (edge1.weight / edge1.dist - edge2.weight / edge2.dist) < tol);
}


// See docs/graph-contraction for explanation of the following code and
// associated vertex and edge maps.
void graph_contract::contract_graph (vertex_map_t &vertex_map,
        edge_map_t &edge_map, vert2edge_map_t &vert2edge_map,
        std::unordered_set <vertex_id_t> verts_to_keep,
        bool has_times)
{
    std::unordered_set <vertex_id_t> verts;
    for (auto v: vertex_map)
    {
        if (verts_to_keep.find (v.first) == verts_to_keep.end ())
            verts.insert (v.first);
    }
    std::unordered_set <edge_id_t> edgelist;
    for (auto e: edge_map)
        edgelist.insert (e.first);

    std::vector <edge_id_t> new_edge_ids;

    unsigned int maxid = 123456;

    while (verts.size () > 0)
    {
        Rcpp::checkUserInterrupt ();
        std::unordered_set <vertex_id_t>::iterator vid = verts.begin ();
        vertex_id_t vtx_id = vertex_map.find (*vid)->first;
        vertex_t vtx = vertex_map.find (*vid)->second;
        std::unordered_set <edge_id_t> edges = vert2edge_map [vtx_id];
        std::unordered_map <edge_id_t, bool> edges_done;
        for (auto e: edges)
            edges_done.emplace (e, false);

        new_edge_ids.clear ();
        new_edge_ids.push_back ("a" + std::to_string (maxid++));

        if ((vtx.is_intermediate_single () || vtx.is_intermediate_double ()) &&
                (edges.size () == 2 || edges.size () == 4))
        {
            if (edges.size () == 4) // is_intermediate_double as well!
                new_edge_ids.push_back ("a" + std::to_string (maxid++));

            // Get the two adjacent vertices
            std::unordered_set <vertex_id_t> nbs = vtx.get_all_neighbours ();
            std::vector <vertex_id_t> two_nbs;
            two_nbs.reserve (2);
            for (vertex_id_t nb: nbs)
                two_nbs.push_back (nb); // size is always 2

            // ensure two edges are same highway type
            bool hwys_are_same = true;
            for (edge_id_t new_edge_id: new_edge_ids)
            {
                vertex_id_t vt_from = "", vt_to = "";
                edge_id_t edge_from_id = "", edge_to_id = "";
                graph_contract::get_to_from (edge_map, edges, two_nbs,
                        vt_from, vt_to, edge_from_id, edge_to_id);
                hwys_are_same = graph_contract::same_hwy_type (edge_map,
                        edge_from_id, edge_to_id);
                if (!hwys_are_same)
                    break;
            }

            if (hwys_are_same)
            {
                // remove intervening vertex:
                vertex_t vt0 = vertex_map [two_nbs [0]];
                vertex_t vt1 = vertex_map [two_nbs [1]];
                // Note that replace neighbour includes bi-directional
                // replacement, so this works for intermediate_double() too
                vt0.replace_neighbour (vtx_id, two_nbs [1]);
                vt1.replace_neighbour (vtx_id, two_nbs [0]);
                vertex_map [two_nbs [0]] = vt0;
                vertex_map [two_nbs [1]] = vt1;

                // construct new edge(s) and remove old ones. There are 2
                // new_edge_ids only for intermediate double vertices
                // (that is, bi-directional).
                for (edge_id_t new_edge_id: new_edge_ids)
                {
                    vertex_id_t vt_from = "", vt_to = "";
                    edge_id_t edge_from_id = "", edge_to_id = "";

                    // get the from and to edges and vertices
                    graph_contract::get_to_from (edge_map, edges, two_nbs,
                            vt_from, vt_to, edge_from_id, edge_to_id);

                    edges.erase (edge_from_id);
                    edges.erase (edge_to_id);
                    // It is possible to have repeated values of two_nbs;
                    // calling these a second time should then just lead to no
                    // contraction
                    if ((vt_from == "") | (vt_to == ""))
                        continue;

                    graph_contract::contract_one_edge (vert2edge_map,
                            vertex_map, edge_map, edgelist, vtx_id, vt_from,
                            vt_to, edge_from_id, edge_to_id, new_edge_id,
                            has_times);
                }
            }
        }
        verts.erase (vtx_id);
    }
}

//' rcpp_contract_graph
//'
//' Removes nodes and edges from a graph that are not needed for routing
//'
//' @param graph graph to be processed
//'
//' @return \code{Rcpp::List} containing one \code{data.frame} with the
//' contracted graph, one \code{data.frame} with the original graph and one
//' \code{data.frame} containing information about the relating edge ids of the
//' original and contracted graph.
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::List rcpp_contract_graph (const Rcpp::DataFrame &graph,
        Rcpp::Nullable <Rcpp::StringVector> &vertlist_in)
{
    std::unordered_set <vertex_id_t> verts_to_keep;
    if (vertlist_in.isNotNull ())
    {
        Rcpp::StringVector vertlist (vertlist_in);
        for (int i = 0; i < vertlist.length (); i ++)
            verts_to_keep.emplace (std::string (vertlist [i]));
    }

    // Get set of all original edge IDs
    std::unordered_set <edge_id_t> original_edges;
    Rcpp::StringVector edge_id = graph ["edge_id"];
    for (auto e: edge_id)
        original_edges.emplace (e);

    vertex_map_t vertices;
    edge_map_t edge_map;
    vert2edge_map_t vert2edge_map;

    bool has_times = graph::graph_from_df (graph, vertices, edge_map, vert2edge_map);

    vertex_map_t vertices_contracted = vertices;
    edge_map_t edge_map_contracted = edge_map;

    graph_contract::contract_graph (vertices_contracted, edge_map_contracted,
            vert2edge_map, verts_to_keep, has_times);

    size_t nedges = edge_map_contracted.size ();

    // These vectors are all for the contracted graph:
    Rcpp::StringVector from_vec (nedges), to_vec (nedges),
        edgeid_vec (nedges);
    Rcpp::NumericVector dist_vec (nedges), weight_vec (nedges),
        time_vec (nedges), timew_vec (nedges);

    size_t map_size = 0; // size of edge map contracted -> original
    unsigned int edge_count = 0;
    for (auto e = edge_map_contracted.begin ();
            e != edge_map_contracted.end (); ++e)
    {
        vertex_id_t from = e->second.get_from_vertex ();
        vertex_id_t to = e->second.get_to_vertex ();
        vertex_t from_vtx = vertices_contracted.at (from);
        vertex_t to_vtx = vertices_contracted.at (to);

        from_vec (edge_count) = from;
        to_vec (edge_count) = to;
        dist_vec (edge_count) = e->second.dist;
        weight_vec (edge_count) = e->second.weight;
        if (has_times)
        {
            time_vec (edge_count) = e->second.time;
            timew_vec (edge_count) = e->second.timew;
        }

        edgeid_vec (edge_count) = e->second.getID ();

        edge_count++;

        if (original_edges.find (e->first) == original_edges.end ())
            map_size += e->second.get_old_edges ().size ();
    }

    // populate the new -> old edge map
    std::vector <edge_id_t> edge_map_new (map_size), edge_map_old (map_size);
    unsigned int count = 0;
    for (auto e = edge_map_contracted.begin ();
            e != edge_map_contracted.end (); ++e)
    {
        if (original_edges.find (e->first) == original_edges.end ())
        {
            std::set <edge_id_t> old_edges = e->second.get_old_edges ();
            for (auto oe: old_edges)
            {
                edge_map_new [count] = e->first;
                edge_map_old [count++] = oe;
            }
        }
    }

    Rcpp::DataFrame contracted;
    if (!has_times)
        contracted = Rcpp::DataFrame::create (
                Rcpp::Named ("edge_id") = edgeid_vec,
                Rcpp::Named ("from") = from_vec,
                Rcpp::Named ("to") = to_vec,
                Rcpp::Named ("d") = dist_vec,
                Rcpp::Named ("d_weighted") = weight_vec,
                Rcpp::_["stringsAsFactors"] = false);
    else
        contracted = Rcpp::DataFrame::create (
                Rcpp::Named ("edge_id") = edgeid_vec,
                Rcpp::Named ("from") = from_vec,
                Rcpp::Named ("to") = to_vec,
                Rcpp::Named ("d") = dist_vec,
                Rcpp::Named ("d_weighted") = weight_vec,
                Rcpp::Named ("time") = time_vec,
                Rcpp::Named ("timew") = timew_vec,
                Rcpp::_["stringsAsFactors"] = false);

    Rcpp::DataFrame edges_new2old = Rcpp::DataFrame::create (
            Rcpp::Named ("edge_new") = edge_map_new,
            Rcpp::Named ("edge_old") = edge_map_old,
            Rcpp::_["stringsAsFactors"] = false);

    return Rcpp::List::create (
            Rcpp::Named ("graph") = contracted,
            Rcpp::Named ("edge_map") = edges_new2old);
}

//' rcpp_merge_cols
//'
//' Merge columns in directed graph to form aggregate undirected columns, and
//' return a corresponding undirected graph useful for visualisation.
//'
//' @param graph The result of a call to \code{dodgr_flows_aggregate/disperse}
//' or similar function resuling in columns of directed values.
//' @return A single vector of aggregate values with non-zero values only for
//' those edges to be retained in the directed graph.
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_merge_cols (Rcpp::DataFrame graph)
{
    std::vector <std::string> from = graph ["from"];
    std::vector <std::string> to = graph ["to"];
    std::vector <double> dist = graph ["d"];
    std::vector <double> wt = graph ["d_weighted"];
    std::vector <double> aggr_var = graph ["merge"]; // always called "merge"

    // vertvert_map just holds index of where pair of vertices were first found.
    // These can only be duplicated once, so only one single value is ever
    // needed.
    std::unordered_map <std::string, int> vertvert_map;
    Rcpp::NumericVector aggr_total (from.size ());
    for (unsigned int i = 0; i < from.size (); i++)
    {
        std::string ft = "a" + from [i] + "b" + to [i],
            tf = "a" + to [i] + "b" + from [i];
        if (vertvert_map.find (ft) == vertvert_map.end () &&
                vertvert_map.find (tf) == vertvert_map.end ())
        {
            vertvert_map.emplace (ft, i);
            aggr_total [i] = aggr_var [i];
        } else
        {
            int where = INFINITE_INT;
            if (vertvert_map.find (ft) != vertvert_map.end ())
                where = vertvert_map.at (ft);
            else if (vertvert_map.find (tf) != vertvert_map.end ())
                where = vertvert_map.at (tf);
            if (where == INFINITE_INT)
                Rcpp::stop ("there is no where; this can never happen"); // # nocov

            aggr_total [static_cast <unsigned int> (where)] += aggr_var [i];
        } 
    }

    return aggr_total;
}
