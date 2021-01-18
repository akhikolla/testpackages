#include "core/exceptions/assert_not_null.hpp"

#include <stack>
#include <queue>
#include <list>

namespace uu {
namespace net {

/**
 * We use the same variable names as in Brandes' algorithm,
 * even if they do not follow the library's style.
 */
template<typename G>
std::unordered_map<const Vertex*, double>
betweenness(
    const G* g
)
{
    core::assert_not_null(g, "betweenness", "g");

    std::unordered_map<const Vertex*, double> C_b;

    for (auto v: *g->vertices())
    {
        C_b[v] = 0;
    }

    for (auto s: *g->vertices())
    {

        std::stack<const Vertex*> S;

        std::unordered_map<const Vertex*, std::list<const Vertex*>> P;

        for (auto w: *g->vertices())
        {
            P[w] = std::list<const Vertex*>();
        }

        std::unordered_map<const Vertex*, double> sigma;

        for (auto t: *g->vertices())
        {
            sigma[t] = 0;
        }

        sigma[s] = 1;


        std::unordered_map<const Vertex*, double> d;

        for (auto t: *g->vertices())
        {
            d[t] = -1;
        }

        d[s] = 0;

        std::queue<const Vertex*> Q;
        Q.push(s);

        while (!Q.empty())
        {
            auto v = Q.front();
            Q.pop();
            S.push(v);

            for (auto w: *g->edges()->neighbors(v, EdgeMode::OUT))
            {
                // w found for the first time?
                if (d[w] < 0)
                {
                    Q.push(w);
                    d[w] = d[v] + 1;
                }

                // shortest path to w via v?
                if (d[w] == d[v] + 1)
                {
                    sigma[w] = sigma[w] + sigma[v];
                    P[w].push_back(v);
                }
            }
        }

        std::unordered_map<const Vertex*, double> delta;

        for (auto v: *g->vertices())
        {
            delta[v] = 0;
        }

        // S returns vertices in order of non-increasing distance from s
        while (!S.empty())
        {
            auto w = S.top();
            S.pop();

            for (auto v: P[w])
            {
                delta[v] = delta[v] + sigma[v]/sigma[w]*(1 + delta[w]);
            }

            if (w != s)
            {
                if (g->is_directed())
                {
                    C_b[w] = C_b[w] + delta[w];
                }

                else
                {
                    C_b[w] = C_b[w] + delta[w]/2;
                }
            }
        }

    }

    return C_b;
}
}
}

