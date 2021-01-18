/*
 * Copyright 2016 neurodata (http://neurodata.io/)
 * Written by Disa Mhembere (disa@jhu.edu)
 *
 * This file is part of k-par-means
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY CURRENT_KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <cassert>

#include "io.hpp"
#include "util.hpp"
#include "clusters.hpp"

namespace kpmeans { namespace base {
void clusters::clear() {
    std::fill(means.begin(), means.end(), 0);
    std::fill(num_members_v.begin(), num_members_v.end(), 0);
    std::fill(complete_v.begin(), complete_v.end(), false);
}

/** \param idx the cluster index.
*/
void clusters::set_mean(const kmsvector& mean, const int idx) {
    if (idx == -1) { // Set all means
        means = mean;
    } else {
        std::copy(mean.begin(), mean.end(),
                this->means.begin()+(idx*ncol));
    }
}

void clusters::set_mean(const double* mean, const int idx) {
    if (idx == -1) { // Set all means
        if (means.size() != (ncol*nclust))
            means.resize(ncol*nclust);
        std::copy(&(mean[0]), &(mean[ncol*nclust]), this->means.begin());
    } else {
        std::copy(&(mean[0]), &(mean[ncol]), this->means.begin()+(idx*ncol));
    }
}

void clusters::finalize(const unsigned idx) {
    if (is_complete(idx)) {
        return;
    }

    if (num_members_v[idx] > 1) { // Less than 2 is the same result
        for (unsigned i = 0; i < ncol; i++) {
            means[(idx*ncol)+i] /= double(num_members_v[idx]);
        }
    }
    complete_v[idx] = true;
}

void clusters::unfinalize(const unsigned idx) {
    if (!is_complete(idx)) {
        return;
    }
    complete_v[idx] = false;

    for (unsigned col = 0; col < ncol; col++) {
        this->means[(ncol*idx) + col] *= (double)num_members_v[idx];
    }
}

void clusters::finalize_all() {
    for (unsigned c = 0;  c < get_nclust(); c++)
        finalize(c);
}

void clusters::unfinalize_all() {
    for (unsigned c = 0;  c < get_nclust(); c++)
        unfinalize(c);
}

void clusters::set_num_members_v(const size_t* arg) {
    std::copy(&(arg[0]), &(arg[nclust]), num_members_v.begin());
}

clusters& clusters::operator=(const clusters& other) {
    this->means = other.get_means();
    this->num_members_v = other.get_num_members_v();
    this->ncol = other.get_ncol();
    this->nclust = other.get_nclust();
    return *this;
}

bool clusters::operator==(const clusters& other) {
    return (get_ncol() == other.get_ncol() &&
            get_nclust() == other.get_nclust() &&
            v_eq(get_num_members_v(), other.get_num_members_v()) &&
            v_eq(get_means(), other.get_means())
           );
}

clusters& clusters::operator+=(clusters& rhs) {
    for (unsigned i = 0; i < size(); i++)
        this->means[i] += rhs[i];

    for (unsigned idx = 0; idx < nclust; idx++)
        num_members_peq(rhs.get_num_members(idx), idx);
    return *this;
}

void clusters::peq(ptr rhs) {
    assert(rhs->size() == size());
    for (unsigned i = 0; i < size(); i++)
        this->means[i] += rhs->get(i);

    for (unsigned idx = 0; idx < nclust; idx++)
        num_members_peq(rhs->get_num_members(idx), idx);
}

void clusters::means_peq(const double* other) {
    for (unsigned i = 0; i < size(); i++)
        this->means[i] += other[i];
}

void clusters::num_members_v_peq(const size_t* other) {
    for (unsigned i = 0; i < num_members_v.size(); i++)
        this->num_members_v[i] += other[i];
}

// Begin Helpers //
const void clusters::print_means() const {
    for (unsigned cl_idx = 0; cl_idx < get_nclust(); cl_idx++) {
#ifndef BIND
        std::cout << "#memb = " << get_num_members(cl_idx) << " ";
#endif
        print_arr<double>(&(means[cl_idx*ncol]), ncol);
    }
#ifndef BIND
    std::cout << "\n";
#endif
}

clusters::clusters(const unsigned nclust, const unsigned ncol) {
    this->nclust = nclust;
    this->ncol = ncol;

    means.resize(ncol*nclust);
    num_members_v.resize(nclust);
    complete_v.assign(nclust, false);
}

clusters::clusters(const unsigned nclust, const unsigned ncol,
        const kmsvector& means) {
    this->nclust = nclust;
    this->ncol = ncol;

    set_mean(means);
    num_members_v.resize(nclust);
    complete_v.assign(nclust, true);
}

const void clusters::print_membership_count() const {
    std::string p = "[ ";
    for (unsigned cl_idx = 0; cl_idx < get_nclust(); cl_idx++) {
        p += std::to_string(get_num_members(cl_idx)) + " ";
    }
    p += "]\n";
#ifndef BIND
    std::cout << p;
#endif
}

// Pruning clusters //
void prune_clusters::reset_s_val_v() {
    std::fill(s_val_v.begin(), s_val_v.end(),
            std::numeric_limits<double>::max());
}

const void prune_clusters::print_prev_means_v() const {
    for (unsigned cl_idx = 0; cl_idx < get_nclust(); cl_idx++) {
        print_arr<double>(&(prev_means[cl_idx*ncol]), ncol);
    }
#ifndef BIND
    std::cout << "\n";
#endif
}
} } // End namespace kpmeans, base
