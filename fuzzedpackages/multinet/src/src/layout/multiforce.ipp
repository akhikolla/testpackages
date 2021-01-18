namespace uu {
namespace net {


template <typename M>
std::map<std::pair<const Vertex*, const typename M::layer_type*>,XYZCoordinates>
multiforce(
    const M* mnet,
    double width,
    double length,
    const std::unordered_map<const typename M::layer_type*,double>& w_intra,
    const std::unordered_map<const typename M::layer_type*,double>& w_inter,
    const std::unordered_map<const typename M::layer_type*,double>& gravity,
    int iterations
)
{
    core::assert_not_null(mnet, "multiforce", "mnet");

    double epsilon = .1; // arbitrary distance considered if two vertices end up having the same coordinates

    std::map<std::pair<const Vertex*, const typename M::layer_type*>,XYZCoordinates> pos;
    std::map<std::pair<const Vertex*, const typename M::layer_type*>,XYZCoordinates> disp;

    if (mnet->actors()->size()==0)
    {
        return pos;
    }

    double temp = std::sqrt(mnet->actors()->size());
    double start_temp = temp;
    double area = width*length;
    std::unordered_map<const typename M::layer_type*,double> k;

    for (auto l: *mnet->layers())
    {
        k[l] = std::sqrt(area/l->vertices()->size());
    }

    for (auto a: *mnet->actors())
    {
        double y = core::drand()*length-length/2;  // suggest to move these here
        double x = core::drand()*width-width/2; // suggest to move these here

        for (auto l: *mnet->layers())
        {
            //double y = drand()*length-length/2;
            //double x = drand()*width-width/2;
            if (!l->vertices()->contains(a))
            {
                continue;
            }

            auto n = std::make_pair(a, l);
            pos[n].x = x;
            pos[n].y = y;
            pos[n].z = mnet->layers()->index_of(l);
        }
    }

    for (int i=0; i<iterations; i++)
    {
        // calculate repulsive forces
        for (auto l: *mnet->layers())
        {
            for (auto a1: *l->vertices())
            {

                auto v = std::make_pair(a1, l);
                disp[v].x = 0;
                disp[v].y = 0;

                for (auto a2: *l->vertices())
                {
                    auto u = std::make_pair(a2, l);

                    if (u == v)
                    {
                        continue;
                    }

                    XYZCoordinates Delta;
                    Delta.x = pos[v].x - pos[u].x;
                    Delta.y = pos[v].y - pos[u].y;
                    //std::cout << "rep " << Delta.x << " " << Delta.y << std::endl;
                    double DeltaNorm = std::sqrt(Delta.x*Delta.x+Delta.y*Delta.y);

                    if (DeltaNorm==0)
                    {
                        DeltaNorm = epsilon;
                    }

                    disp[v].x = disp[v].x + Delta.x/DeltaNorm*fr(DeltaNorm,k.at(l))*w_intra.at(l);
                    disp[v].y = disp[v].y + Delta.y/DeltaNorm*fr(DeltaNorm,k.at(l))*w_intra.at(l);
                }

                // add effect of gravity, to prevent disc. components from diverging
                double DeltaNorm = std::sqrt(pos[v].x*pos[v].x+pos[v].y*pos[v].y);

                if (DeltaNorm==0)
                {
                    DeltaNorm = epsilon;
                }

                disp[v].x = disp[v].x - pos[v].x/DeltaNorm*fain(DeltaNorm,k.at(l))*gravity.at(l);
                disp[v].y = disp[v].y - pos[v].y/DeltaNorm*fain(DeltaNorm,k.at(l))*gravity.at(l);

            }
        }

        // calculate attractive forces inside each layer
        for (auto l: *mnet->layers())
        {
            for (auto e: *l->edges())
            {
                auto v = std::make_pair(e->v1, l);
                auto u = std::make_pair(e->v2, l);
                XYZCoordinates Delta;
                Delta.x = pos[v].x - pos[u].x;
                Delta.y = pos[v].y - pos[u].y;
                //std::cout << "a-in " << Delta.x << " " << Delta.y << std::endl;
                double DeltaNorm = std::sqrt(Delta.x*Delta.x+Delta.y*Delta.y);

                if (DeltaNorm==0)
                {
                    DeltaNorm = epsilon;
                }

                disp[v].x = disp[v].x - Delta.x/DeltaNorm*fain(DeltaNorm,k.at(l))*w_intra.at(l);
                disp[v].y = disp[v].y - Delta.y/DeltaNorm*fain(DeltaNorm,k.at(l))*w_intra.at(l);
                disp[u].x = disp[u].x + Delta.x/DeltaNorm*fain(DeltaNorm,k.at(l))*w_intra.at(l);
                disp[u].y = disp[u].y + Delta.y/DeltaNorm*fain(DeltaNorm,k.at(l))*w_intra.at(l);
            }
        }

        // calculate attractive forces across layers
        for (auto a: *mnet->actors())
        {
            for (auto l1: *mnet->layers())
            {
                if (!l1->vertices()->contains(a))
                {
                    continue;
                }

                for (auto l2: *mnet->layers())
                {
                    if (l1>=l2)
                    {
                        continue;
                    }

                    // Or: if (v >= u) continue; ?
                    if (!l2->vertices()->contains(a))
                    {
                        continue;
                    }

                    auto v = std::make_pair(a, l1);
                    auto u = std::make_pair(a, l2);

                    XYZCoordinates Delta;
                    Delta.x = pos[v].x - pos[u].x;
                    Delta.y = pos[v].y - pos[u].y;
                    //std::cout << "a-inter " << Delta.x << " " << Delta.y << std::endl;
                    double DeltaNorm = std::sqrt(Delta.x*Delta.x+Delta.y*Delta.y);

                    if (DeltaNorm==0)
                    {
                        DeltaNorm = epsilon;
                    }

                    disp[v].x = disp[v].x - Delta.x/DeltaNorm*fainter(DeltaNorm,k.at(l1))*w_inter.at(l1);
                    disp[v].y = disp[v].y - Delta.y/DeltaNorm*fainter(DeltaNorm,k.at(l1))*w_inter.at(l1);
                    disp[u].x = disp[u].x + Delta.x/DeltaNorm*fainter(DeltaNorm,k.at(l2))*w_inter.at(l2);
                    disp[u].y = disp[u].y + Delta.y/DeltaNorm*fainter(DeltaNorm,k.at(l2))*w_inter.at(l2);
                }
            }
        }

        // assign new positions
        for (auto l: *mnet->layers())
        {
            for (auto a: *l->vertices())
            {
                auto v = std::make_pair(a, l);
                double dispNorm = std::sqrt(disp[v].x*disp[v].x+disp[v].y*disp[v].y);

                if (dispNorm==0)
                {
                    dispNorm = epsilon;
                }

                pos[v].x = pos[v].x + (disp[v].x/dispNorm)*std::min(dispNorm,temp);
                pos[v].y = pos[v].y + (disp[v].y/dispNorm)*std::min(dispNorm,temp);
                pos[v].x = std::min(width/2, std::max(-width/2, pos[v].x)); // suggest to remove - it might actually be useful...
                pos[v].y = std::min(length/2, std::max(-length/2, pos[v].y)); // suggest to remove - it might actually be useful...
            }
        }

        // reduce the temperature
        temp -= start_temp/iterations;
    }

    return pos;
}


}
}

