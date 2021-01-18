namespace uu {
namespace net {

template <typename M, typename G>
std::unique_ptr<CommunityStructure<VertexLayerCommunity<const G>>>
generalized_louvain(
    const M* mnet,
    double gamma,
    double omega,
    size_t limit
)
{
    glouvain g;
    return g.fit<M,G>(mnet, "move", gamma, omega, limit);
}


template <typename M, typename G>
std::unique_ptr<CommunityStructure<VertexLayerCommunity<const G>>>
glouvain::fit(
    const M* mnet,
    const std::string& m,
    double gamma,
    double omega,
    size_t limit
)
{
    // @todo check UNDIRECTED
    // @todo check ORDERED
    bool ordered = mnet->is_ordered();

    double (*move_func)(group_index &, int, const Eigen::SparseMatrix<double>&);

    if ("moverandw" == m)
    {
        move_func = &moverandw;
    }

    else
    {
        move_func = &move;
    }

    std::vector<int> S(mnet->actors()->size() * mnet->layers()->size());
    std::iota(S.begin(), S.end(), 0);

    std::vector<int> y, Sb;
    y = S;

    double eps = 2.2204e-12;
    double dtot = 0;
    Eigen::SparseMatrix<double> B, M_;

    double twoum = 0;

    if (limit < y.size())
    {
        std::vector<Eigen::SparseMatrix<double>> A = to_adjacency_matrices(mnet);

        for (auto m: A)
        {
            twoum += m.nonZeros();
        }

        twoum = twoum + (A[0].rows() * A.size() * (A.size() - 1) * omega);

        metanet meta(A, gamma, omega, ordered);
        meta.assign(S);

        while (Sb != S)
        {
            Sb = S;
            std::vector<int> yb;

            while (yb != y)
            {
                double dstep = 1.0;

                while (yb != y && (dstep/dtot > 2 * eps) && (dstep > 10 * eps))
                {
                    yb = y;
                    dstep = 0;

                    group_index g(y);

                    for (int i: cutils::range(meta.get(0).rows(), true))
                    {
                        double di = move_func(g, i, meta.get(i));
                        dstep = dstep + di;
                    }

                    dtot = dtot + dstep;
                    y = g.toVector();
                }

                yb = y;
            }

            S = mapV2I(y, S);
            y = cutils::unique(S);

            if (Sb == S)
            {
                std::vector<unsigned int> partition(S.begin(), S.end());
                return to_community_structure<M,G>(mnet, partition);
            }

            meta.assign(S);

            if (y.size() < limit)
            {
                std::vector<Eigen::Triplet<double>> tlist;
                Eigen::SparseMatrix<double> t = meta.get(0);
                tlist.reserve(t.nonZeros() * y.size());

                for (int i = 0; i < (int) y.size(); i++)
                {
                    Eigen::SparseMatrix<double> tmp = meta.get(i);

                    for (int j = 0; j < tmp.outerSize(); ++j)
                    {
                        for (Eigen::SparseMatrix<double>::InnerIterator it(tmp, j); it; ++it)
                        {
                            tlist.push_back(Eigen::Triplet<double>(it.row(), i, it.value()));
                        }
                    }
                }

                B = Eigen::SparseMatrix<double>(t.rows(), y.size());
                B.setFromTriplets(tlist.begin(), tlist.end());
                M_ = B;
                break;
            }
        }
    }

    else
    {
        B = modularity_matrix(twoum, to_adjacency_matrices(mnet), gamma, omega, ordered);
        M_ = B;
    }

    std::vector<int> S2(B.rows());
    std::iota(S2.begin(), S2.end(), 0);
    Sb.clear();

    /*
    std::cout << Eigen::MatrixXd(M_) << std::endl;
    std::cout << "--------------------" << std::endl;
    std::cout << core::to_string(S2) << std::endl;
    std::cout << "--------------------" << std::endl;
    */

    while (Sb != S2)
    {
        Sb = S2;
        std::vector<int> yb;

        while (yb != y)
        {
            double dstep = 1.0;

            while (yb != y && (dstep/dtot > 2 * eps) && (dstep > 10 * eps))
            {
                yb = y;
                dstep = 0;

                group_index g(y);

                for (int i: cutils::range(M_.cols(), true))
                {
                    double di = move_func(g, i, M_.col(i));
                    //std::cout << core::to_string(g.nodes) << std::endl;
                    dstep = dstep + di;
                }

                dtot = dtot + dstep;
                y = g.toVector();
            }

            yb = y;
        }

        S = mapV2I(y, S);
        S2 = mapV2I(y, S2);

        if (Sb == S2)
        {
            break;
        }

        M_ = metanetwork(B, S2);
        y = cutils::unique(S2);

        /*
        std::cout << Eigen::MatrixXd(M_) << std::endl;
        std::cout << "--------------------" << std::endl;
        std::cout << core::to_string(S2) << std::endl;
        std::cout << "--------------------" << std::endl;
         */
    }

    std::vector<unsigned int> partition(S.begin(), S.end());
    return to_community_structure<M,G>(mnet, partition);
}



}
}
