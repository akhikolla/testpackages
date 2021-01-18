namespace uu {
namespace net {


void
flatten_weighted(
    LayerIterator begin,
    LayerIterator end,
    W* target
)
{

    for (auto layer=begin; layer!=end; ++layer)
    {
        // force actors? @todo

        weighted_graph_union(*layer,target,target);
    }
}

}
}

#endif
