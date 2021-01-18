namespace uu {
namespace net {


template<typename G>
BFS<G>::
BFS(
    const G* g,
    const Vertex* v,
    EdgeMode mode
) : g(g), mode(mode)
{
    core::assert_not_null(g, "BFS", "g");
    core::assert_not_null(v, "BFS", "v");

    queue.push(v);
    processed.insert(v);

}


template<typename G>
const Vertex*
BFS<G>::
get_next(
)
{
    if (queue.size() == 0)
    {
        return nullptr;
    }

    const Vertex* res = queue.front();
    queue.pop();

    for (auto n: *g->edges()->neighbors(res, mode))
    {
        if (processed.find(n) == processed.end())
        {
            queue.push(n);
            processed.insert(n);
        }
    }

    return res;
}


}
}

