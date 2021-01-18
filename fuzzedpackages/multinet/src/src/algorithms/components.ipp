namespace uu {
namespace net {

template<typename G>
std::vector<int>
components(
    const G* g
)
{
    core::assert_not_null(g, "components", "g");

    std::vector<int> membership(order(g), -1);

    int comp_id = 0;
    size_t vertex_pos = 0;

    for (auto v: *g->vertices())
    {
        if (membership.at(vertex_pos) != -1)
        {
            // do nothing: vertex already processed
        }
        else
        {
            // recursively traverse the whole component

            BFS<G> bfs(g, v, EdgeMode::INOUT);
            const Vertex* current_vertex;

            while ((current_vertex = bfs.get_next()) != nullptr)
            {
                size_t pos_current_vertex = g->vertices()->index_of(current_vertex);

                membership.at(pos_current_vertex) = comp_id;
            }

            comp_id++;
        }

        vertex_pos++;
    }

    return membership;
}


}
}

