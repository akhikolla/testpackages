#include "community/_impl/cutils.hpp"

namespace uu {
namespace net {

std::unique_ptr<CommunityStructure<Community<const Vertex*>>>
to_community_structure(
    std::unordered_map<const Vertex*, size_t> membership
)
{
    auto result = std::make_unique<CommunityStructure<Community<const Vertex*>>>();
    std::unordered_map<int, std::vector<const Vertex*>> communities;

    for (auto pair: membership)
    {
        communities[pair.second].push_back(pair.first);
    }

    for (auto pair: communities)
    {
        auto community = std::make_unique<Community<const Vertex*>>();

        for (auto n: pair.second)
        {
            community->add(n);
        }

        result->add(std::move(community));
    }

    return result;
}


Eigen::SparseMatrix<double>
cutils::supraA(
    const std::vector<Eigen::SparseMatrix<double>>& a,
    double eps,
    bool use_node_degrees,
    bool use_self_loop
)
{
    Eigen::SparseMatrix<double> A;

    if (use_self_loop)
    {
        A = block_diag(a, eps);
    }

    else
    {
        A = block_diag(a, 0);
    }

    size_t L = a.size();
    size_t N = a[0].rows();

    for (size_t i = 0; i  < L - 1; ++i)
    {
        for (size_t j = i + 1; j < L; ++j)
        {
            Eigen::MatrixXd d = cutils::sparse_sum(a[i].cwiseProduct(a[j]), 1);

            std::vector<Eigen::Triplet<double>> tlist;
            tlist.reserve(a[0].rows());

            for (int k = 0; k < A.outerSize(); k++)
            {
                for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it)
                {
                    tlist.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
                }
            }

            int ix_i = i * N;
            int ix_j = j * N;

            if (use_node_degrees)
            {
                for (int k = 0; k < a[i].rows(); k++)
                {
                    double intra = d(k, 0) + eps;
                    tlist.push_back(Eigen::Triplet<double>((ix_i + k), (ix_j + k), intra));
                    tlist.push_back(Eigen::Triplet<double>((ix_j + k), (ix_i + k), intra));
                }
            }

            else
            {
                for (int k = 0; k < a[i].rows(); k++)
                {
                    tlist.push_back(Eigen::Triplet<double>((ix_i + k), (ix_j + k), eps));
                    tlist.push_back(Eigen::Triplet<double>((ix_j + k), (ix_i + k), eps));
                }
            }

            A.setFromTriplets(tlist.begin(), tlist.end());

        }
    }

    return A;
}



Eigen::SparseMatrix<double>
cutils::ordered_supraA(
    const std::vector<Eigen::SparseMatrix<double>>& a,
    double eps,
    bool use_node_degrees,
    bool use_self_loop
)
{
    Eigen::SparseMatrix<double> A;

    if (use_self_loop)
    {
        A = block_diag(a, eps);
    }

    else
    {
        A = block_diag(a, 0);
    }

    size_t L = a.size();
    size_t N = a[0].rows();

    for (size_t i = 0; i < L - 1; ++i)
    {
        Eigen::MatrixXd d = cutils::sparse_sum(a[i].cwiseProduct(a[i+1]), 1);

        std::vector<Eigen::Triplet<double>> tlist;
        tlist.reserve(a[0].rows());

        for (int k = 0; k < A.outerSize(); k++)
        {
            for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it)
            {
                tlist.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
            }
        }

        int ix_i = i * N;
        // Considering only adjacent layers
        int ix_j = (i + 1) * N;

        if (use_node_degrees)
        {
            for (int k = 0; k < a[i].rows(); k++)
            {
                double intra = d(k, 0) + eps;
                tlist.push_back(Eigen::Triplet<double>((ix_i + k), (ix_j + k), intra));
                tlist.push_back(Eigen::Triplet<double>((ix_j + k), (ix_i + k), intra));
            }
        }

        else
        {
            for (int k = 0; k < a[i].rows(); k++)
            {
                tlist.push_back(Eigen::Triplet<double>((ix_i + k), (ix_j + k), eps));
                tlist.push_back(Eigen::Triplet<double>((ix_j + k), (ix_i + k), eps));
            }
        }

        A.setFromTriplets(tlist.begin(), tlist.end());
    }

    return A;
}

Eigen::SparseMatrix<double>
cutils::block_diag(
    const std::vector<Eigen::SparseMatrix<double>>& a,
    double eps
)
{
    Eigen::SparseMatrix<double> m = Eigen::SparseMatrix<double>(
                                        a[0].rows() * a.size(), a[0].cols() * a.size());

    int nnz = 0;

    for (auto l: a)
    {
        nnz += l.nonZeros();
    }

    size_t r, c;
    r = 0;
    c = 0;

    std::vector<Eigen::Triplet<double>> tlist;
    tlist.reserve(nnz);

    for (size_t i = 0; i < a.size(); i++)
    {
        for (int j = 0; j < a[i].outerSize(); j++)
        {
            for (Eigen::SparseMatrix<double>::InnerIterator it(a[i], j); it; ++it)
            {
                tlist.push_back(Eigen::Triplet<double>(r + it.row(), c + it.col(), it.value()));
            }
        }

        r += a[i].rows();
        c += a[i].cols();
    }

    for (int i = 0; i < m.rows(); i++)
    {
        tlist.push_back(Eigen::Triplet<double>(i, i, eps));
    }

    m.setFromTriplets(tlist.begin(), tlist.end());
    return m;
}



Eigen::MatrixXd
cutils::sparse_sum(
    const Eigen::SparseMatrix<double>& X,
    int axis
)
{
    Eigen::MatrixXd d = Eigen::MatrixXd::Zero(X.rows(), 1);

    for (int i = 0; i < X.outerSize(); i++)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(X, i); it; ++it)
        {
            if (axis)
            {
                d(it.col(), 0) = it.value() + d(it.col(), 0);
            }

            else
            {
                d(it.row(), 0) = it.value() + d(it.row(), 0);
            }
        }
    }

    return d;
}

std::vector<int>
cutils::unique(
    std::vector<int> y
)
{
    std::sort(y.begin(), y.end());
    std::vector<int>::iterator it;
    it = std::unique(y.begin(), y.end());
    y.resize(std::distance(y.begin(), it));
    return y;
}

std::vector<int>
cutils::range(
    size_t size,
    bool randomize
)
{
    std::vector<int> range(size);
    std::iota(range.begin(), range.end(), 0);

    if (randomize)
    {
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::shuffle (range.begin(), range.end(), std::default_random_engine(seed));
    }

    return range;
}


Eigen::SparseMatrix<double>
modularity_matrix(
    double &twoum,
    const std::vector<Eigen::SparseMatrix<double>>& a,
    double gamma,
    double omega,
    bool ordered
)
{
    size_t L = a.size();
    size_t N = a[0].rows();

    Eigen::SparseMatrix<double> B(N * L, N * L);
    //std::cout << B.rows() << " " << B.cols() << " " << B.nonZeros() << std::endl;
    B.reserve(Eigen::VectorXi::Constant(N * L, N * L));

    twoum = 0.0;

    std::vector<Eigen::Triplet<double>> tlist;
    tlist.reserve(N * L);

    for (size_t i = 0; i < L; i++)
    {
        Eigen::MatrixXd kout = cutils::sparse_sum(a[i], 0);
        Eigen::MatrixXd kin = cutils::sparse_sum(a[i], 1);
        double mm = kout.array().sum();
        twoum = twoum + mm;

        Eigen::SparseMatrix<double> tmp, tmp1;
        tmp = (Eigen::SparseMatrix<double>(a[i].transpose()) + a[i]) / 2;

        Eigen::MatrixXd tmp2 = kin * kout.transpose();

        Eigen::MatrixXd tmp3;

        if (mm != 0)
        {
            tmp3 = gamma / 2 * (tmp2 + tmp2) / mm;
        }

        else
        {
            tmp3 = gamma / 2 * (tmp2 + tmp2);
        }

        tmp = tmp - tmp3.sparseView();

        for (size_t j = 0; j < N; j++)
        {
            for (size_t k = 0; k < N; k++)
            {
                auto coeff = tmp.coeff(j, k);

                if (coeff != 0)
                    tlist.push_back(
                        Eigen::Triplet<double>(j + (i * N), k + (i * N), tmp.coeff(j, k)));
            }
        }
    }

    if (ordered)
    {
        for (size_t i = 0; i < L - 1; ++i)
        {
            int ix_i = i * N;
            int ix_j = (i + 1) * N; // Only considering adjacent layers

            for (int k = 0; k < a[i].rows(); k++)
            {
                tlist.push_back(Eigen::Triplet<double>((ix_i + k), (ix_j + k), omega));
                tlist.push_back(Eigen::Triplet<double>((ix_j + k), (ix_i + k), omega));
            }
        }
    }

    else
    {
        for (size_t i = 0; i  < L - 1; ++i)
        {
            for (size_t j = i + 1; j < L; ++j)
            {
                int ix_i = i * N;
                int ix_j = j * N;

                for (int k = 0; k < a[i].rows(); k++)
                {
                    tlist.push_back(Eigen::Triplet<double>((ix_i + k), (ix_j + k), omega));
                    tlist.push_back(Eigen::Triplet<double>((ix_j + k), (ix_i + k), omega));
                }
            }
        }
    }

    B.setFromTriplets(tlist.begin(), tlist.end());

    twoum = twoum + (N * L * (L - 1) * omega);

    //std::cout << Eigen::MatrixXd(B) << std::endl;
    //std::cout << twoum << std::endl;

    return B;
}

}
}
