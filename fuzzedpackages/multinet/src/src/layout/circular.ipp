namespace uu {
namespace net {

template <typename M>
std::map<std::pair<const Vertex*, const typename M::layer_type*>,XYZCoordinates>
circular(
    const M* mnet,
    double radius
);


template <typename M>
std::map<std::pair<const Vertex*, const typename M::layer_type*>,XYZCoordinates>
circular(
    const M* mnet,
    double radius
)
{
    std::map<std::pair<const Vertex*, const typename M::layer_type*>,XYZCoordinates> pos;
    double pi = 3.14159265358979323846;

    if (mnet->actors()->size()==0)
    {
        return pos;
    }

    double angle_offset = 360.0/mnet->actors()->size();
    int i=0;

    for (auto a: *mnet->actors())
    {
        double degree = i*angle_offset;
        double radians = degree*pi/180;
        double x = std::cos(radians)*radius;
        double y = std::sin(radians)*radius;

        for (auto l: *mnet->layers())
        {
            if (l->vertices()->contains(a))
            {
                auto n = std::make_pair(a, l);
                pos[n].x = x;
                pos[n].y = y;
                pos[n].z = mnet->layers()->index_of(l); // @todo not very efficient
            }
        }

        i++;
    }

    return pos;
}

}
}

