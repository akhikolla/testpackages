//
// Created by Jiandong Wang on 2/15/20.
//
// Copyright (c) NMSU Song lab

#include "Clusters.h"

// testing purpose only
Cluster::Cluster(vector<vector<double> > medians) {
    this->num_clusters = medians.size();
    this->cluster_medians = medians;
}


// testing purpose only
Cluster::Cluster(vector<vector<double> > medians, vector<vector<vector<double> > > data) {
    this->num_clusters = medians.size();
    this->cluster_medians = medians;
    this->cluster_points = data;
}

Cluster::Cluster(int k, vector<int> labels, vector<vector<double> > data) {
    // get number of clusters
    this->num_clusters = k;
    int dims = data[0].size();
    this->cluster_points = vector<vector<vector<double> > >(num_clusters,
                                                            vector<vector<double>>(dims, vector<double>(0)));

    // get clusters
    vector<int>::iterator label_iter = labels.begin();
    vector<vector<double> >::iterator data_iter = data.begin();
    //loop for clusters
    for (; label_iter != labels.end() and data_iter != data.end(); ++label_iter, ++data_iter) {
        //loop for dimensions
        for (int i = 0; i < dims; ++i) {
            this->cluster_points[*label_iter][i].push_back((*data_iter)[i]);
        }
    }

    //calculate median
    this->cluster_medians = vector<vector<double> >(this->num_clusters, vector<double>(dims, 0));
    double mid = 0;
    double mid_index = 0;
    for (int i = 0; i < this->num_clusters; ++i) {
        for (int j = 0; j < dims; ++j) {
            vector<double> temp = this->cluster_points[i][j];
            sort(temp.begin(), temp.end());
            if (temp.size() % 2 == 0) {//mid_point is in between
                mid_index = temp.size() / 2;
                mid = (temp[mid_index] + temp[mid_index - 1]) / 2.0;
            } else {//mid_point is in data
                mid_index = temp.size() / 2.0;
                mid = temp[floor(mid_index)];
            }
            this->cluster_medians[i][j] = mid;
        }
    }
}

Cluster::Cluster(vector<int> labels, vector<vector<double> > medians, vector<vector<double> > data) {
    // get number of clusters
    this->num_clusters = medians.size();
    int dims = medians[0].size();
    this->cluster_points = vector<vector<vector<double> > >(num_clusters,
                                                            vector<vector<double>>(dims, vector<double>(0)));

    // get clusters
    vector<int>::iterator label_iter = labels.begin();
    vector<vector<double> >::iterator data_iter = data.begin();
    //loop for clusters
    for (; label_iter != labels.end() and data_iter != data.end(); ++label_iter, ++data_iter) {
        //loop for dimensions
        for (int i = 0; i < dims; ++i) {
            this->cluster_points[*label_iter][i].push_back((*data_iter)[i]);
        }
    }
    this->cluster_medians = medians;
}


vector<vector<vector<double> > > Cluster::get_points() {
    return this->cluster_points;
}

vector<vector<double> > Cluster::get_points(int index) {
    return this->cluster_points[index];
}

vector<vector<double> > Cluster::get_medians() {
    return this->cluster_medians;
}

vector<double> Cluster::get_medians(int index) {
    return this->cluster_medians[index];
}

void Cluster::set_grids(grid G) {
    this->grids = G;
}

grid Cluster::get_grids() {
    return this->grids;
}

vector<int> Cluster::sort_clusters(int dim) {
    // initialize original index locations
    vector<int> idx(this->num_clusters);
    iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    sort(idx.begin(), idx.end(),
         [this, dim](size_t i1, size_t i2) { return this->cluster_medians[i1][dim] < this->cluster_medians[i2][dim]; });
    return idx;
}

int Cluster::get_dims() {
    return this->cluster_medians[0].size();
}
