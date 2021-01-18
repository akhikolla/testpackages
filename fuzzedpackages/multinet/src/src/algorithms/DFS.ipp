namespace uu {
namespace net {

template<typename G>
DFS<G>::
DFS(
    const G* g,
    const Vertex* v,
    EdgeMode mode
) : g(g), mode(mode)
{
    core::assert_not_null(g, "BFS", "g");
    core::assert_not_null(v, "BFS", "v");

    stack.push(v);
    processed.insert(v);

}


template<typename G>
const Vertex*
DFS<G>::
get_next(
)
{
    if (stack.size() == 0)
    {
        return nullptr;
    }

    const Vertex* res = stack.top();
    stack.pop();

    for (auto n: *g->edges()->neighbors(res, mode))
    {
        if (processed.find(n) == processed.end())
        {
            stack.push(n);
            processed.insert(n);
        }
    }

    return res;
}


}
}

