#include "community/glouvain.hpp"
#include "community/_impl/cutils.hpp"

namespace uu {
namespace net {

set_type
possible_moves(
    group_index & g,
    int node,
    const Eigen::SparseMatrix<double>& mod
)
{
    set_type unique_groups(g.n_groups);
    unique_groups.insert(g.nodes[node]);

    //add nodes with potential positive contribution to unique_groups
    for (int i = 0; i < mod.outerSize(); i++)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(mod, i); it; ++it)
        {
            if (it.value() > 0)
            {
                unique_groups.insert(g.nodes[it.row()]);
            }
        }
    }

    return unique_groups;
}

//calculates changes in modularity for sparse modularity matrix
map_type
mod_change(
    group_index & g,
    const Eigen::SparseMatrix<double>& mod,
    set_type & unique_groups,
    int current_node
)
{
    int current_group = g.nodes[current_node];
    map_type mod_c;

    double mod_current = mod.coeff(current_node, 0);

    for (set_type::iterator it=unique_groups.begin(); it!=unique_groups.end(); ++it)
    {
        mod_c[*it] = 0;
    }

    for (int i = 0; i < mod.outerSize(); i++)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(mod, i); it; ++it)
        {
            if (unique_groups.count(g.nodes[it.row()]))
            {
                mod_c[g.nodes[it.row()]] += it.value();
            }
        }
    }

    mod_c[current_group]-=mod_current;
    mod_current=mod_c[current_group];

    for (set_type::iterator it=unique_groups.begin(); it!=unique_groups.end(); ++it)
    {
        mod_c[*it]-= mod_current;
    }

    /*
    std::cout << "Mod.change for " << current_node << ":  ";
    for (auto mod_c_it: mod_c)
    {
        std::cout << mod_c_it.second << "  ";
    }
    std::cout << std::endl;
    */
    return mod_c;
}

//find moves that improve modularity
move_list
positive_moves(
    set_type & unique_groups,
    map_type & mod_c
)
{
    move_list moves;

    for (set_type::iterator it=unique_groups.begin(); it!=unique_groups.end(); ++it)
    {
        if (mod_c[*it]>NUM_TOL)
        {
            moves.first.push_back(*it);
            moves.second.push_back(mod_c[*it]);
        }
    }

    return moves;
}

//move best move
double
move(
    group_index & g,
    int node,
    const Eigen::SparseMatrix<double>& mod
)
{
    set_type unique_groups = possible_moves(g, node, mod);
    map_type mod_c = mod_change(g, mod, unique_groups, node);
    double mod_max = 0;
    double d_step = 0;
    int group_move = g.nodes[node]; //stay in current group if no improvement

    for (set_type::iterator it=unique_groups.begin(); it!=unique_groups.end(); ++it)
    {
        if (mod_c[*it]>mod_max)
        {
            mod_max=mod_c[*it];
            group_move=*it;
        }
    }

    //move current node to most optimal group
    if (mod_max > NUM_TOL)
    {
        g.move(node, group_move);
        d_step+=mod_max;
    }

    return d_step;
}

// Random engine used for random movement function, moverandw
std::default_random_engine
generator(
    (unsigned int)time(0)
);

//move to random group with probability proportional to increase in modularity
double
moverandw(
    group_index & g,
    int node,
    const Eigen::SparseMatrix<double>& mod
)
{
    set_type unique_groups = possible_moves(g, node, mod);
    map_type mod_c = mod_change(g, mod, unique_groups, node);

    //find modularity increasing moves
    move_list mod_pos = positive_moves(unique_groups, mod_c);

    //move node to a random group that increases modularity with probability proportional to the increase
    double d_step=0;

    if (!mod_pos.first.empty())
    {
        std::discrete_distribution<int> randindex(mod_pos.second.begin(),mod_pos.second.end());
        int randmove = randindex(generator);
        g.move(node, mod_pos.first[randmove]);
        d_step = mod_pos.second[randmove];
    }

    return d_step;
}

///

std::vector<int>
glouvain::mapV2I(
    const std::vector<int>& a,
    const std::vector<int>& b
) const
{
    std::vector<int> v(b.size());

    for (size_t i = 0; i < b.size(); i++)
    {
        v[i] = a[b[i]];
    }

    return v;
}

double
glouvain::Q(
    const Eigen::SparseMatrix<double>& M,
    const std::vector<int>& y,
    double twoum
) const
{
    Eigen::SparseMatrix<double> P(y.size(), y.size());

    std::vector<Eigen::Triplet<double>> tlist;
    tlist.reserve(y.size());

    for (size_t i = 0; i < y.size(); i++)
    {
        tlist.push_back(Eigen::Triplet<double>(i, i, 1));
    }

    P.setFromTriplets(tlist.begin(), tlist.end());
    return Eigen::MatrixXd((P*M).cwiseProduct(P)).sum() / twoum;
}


double
glouvain::Q_handle(
    metanet& meta,
    const std::vector<int>& y,
    double twoum
)
{
    Eigen::SparseMatrix<double> P(y.size(), y.size());

    std::vector<Eigen::Triplet<double>> tlist;
    tlist.reserve(y.size());

    for (size_t i = 0; i < y.size(); i++)
    {
        tlist.push_back(Eigen::Triplet<double>(i, i, 1));
    }

    P.setFromTriplets(tlist.begin(), tlist.end());

    double Q = 0;

    for (size_t i = 0; i < y.size(); i++)
    {
        Q += Eigen::MatrixXd(Eigen::SparseMatrix<double>(P * meta.get(i)).transpose() * P.col(i)).array().sum();
    }

    return Q / twoum;
}



Eigen::SparseMatrix<double>
glouvain::metanetwork(
    const Eigen::SparseMatrix<double>& B,
    const std::vector<int>& S2
) const
{
    Eigen::SparseMatrix<double> PP(B.rows(), *std::max_element(S2.begin(), S2.end()) + 1);
    PP.reserve(B.rows());

    std::vector<Eigen::Triplet<double>> tlist;
    tlist.reserve(B.rows());

    for (size_t i = 0; i < S2.size(); i++)
    {
        tlist.push_back(Eigen::Triplet<double>(i, S2[i], 1));
    }

    PP.setFromTriplets(tlist.begin(), tlist.end());
    return PP.transpose() * B * PP;
}

}
}
