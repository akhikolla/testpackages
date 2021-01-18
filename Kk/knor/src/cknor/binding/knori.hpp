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

#ifndef __KNORI_HPP__
#define __KNORI_HPP__

#include "io.hpp"
#include "../libauto/kmeans.hpp"

#include "kmeans_coordinator.hpp"
#include "kmeans_task_coordinator.hpp"
#include "util.hpp"

#ifdef USE_NUMA
#include "numa_reorg.hpp"
namespace kpmbind = kpmeans::binding;
#endif

namespace kpmprune = kpmeans::prune;
namespace kpmomp = kpmeans::omp;

namespace kpmeans { namespace base {
    // NOTE: It is the callers job to allocate/free data & p_centers

kmeans_t kmeans(double* data, const size_t nrow,
        const size_t ncol, const unsigned k,
        size_t max_iters=std::numeric_limits<size_t>::max(),
        unsigned nnodes=get_num_nodes(),
        unsigned nthread=get_num_omp_threads(),
        double* p_centers=NULL, std::string init="kmeanspp",
        double tolerance=-1, std::string dist_type="eucl",
        bool omp=false, bool numa_opt=false) {

    if (p_centers)
        init = "none";

    kmeans_t ret;

#ifdef _OPENMP
    if (omp) {
        std::vector<double> centroids (k*ncol);
        std::fill(centroids.begin(), centroids.end(), 0);
        std::vector<unsigned> assignments(nrow);
        std::fill(assignments.begin(), assignments.end(), 0);
        std::vector<size_t> counts(k);
        std::fill(counts.begin(), counts.end(), 0);

        if (max_iters < std::numeric_limits<size_t>::max())
            max_iters++; // NOTE: This accounts for difference in pthread v. omp

        ret = kpmomp::compute_min_kmeans(data, &centroids[0], &assignments[0],
                &counts[0], nrow, ncol, k, max_iters, nthread, init, tolerance,
                dist_type);
    } else {
#endif
        kpmprune::kmeans_task_coordinator::ptr kc =
            kpmprune::kmeans_task_coordinator::create(
                    "", nrow, ncol, k, max_iters, nnodes,
                    nthread, p_centers,
                    init, tolerance, dist_type);
        if (numa_opt) {
#ifdef USE_NUMA
            kpmbind::memory_distributor<double>::ptr md =
                kpmbind::memory_distributor<double>::create(data,
                        nnodes, nrow, ncol);

            md->numa_reorg(kc->get_threads());
            ret = kc->run_kmeans(NULL, true);
#else
            ret = kc->run_kmeans(data);
#endif
        } else {
            ret = kc->run_kmeans(data);
        }
#if _OPENMP
    }
#endif

    return ret;
}

kmeans_t kmeans(const std::string datafn, const size_t nrow,
        const size_t ncol, const unsigned k,
        size_t max_iters=std::numeric_limits<size_t>::max(),
        unsigned nnodes=get_num_nodes(),
        unsigned nthread=get_num_omp_threads(),
        double* p_centers=NULL, std::string init="kmeanspp",
        double tolerance=-1, std::string dist_type="eucl",
        bool omp=false) {

    if (p_centers)
        init = "none";

    kmeans_t ret;

#ifdef _OPEMP
    if (omp) {
        // Read all the data
        std::vector<double> data(nrow*ncol);
        bin_io<double> br(datafn, nrow, ncol);
        br.read(&data);

        std::vector<double> centroids(k*ncol);
        std::fill(centroids.begin(), centroids.end(), 0);
        std::vector<unsigned> assignments(nrow);
        std::fill(assignments.begin(), assignments.end(), 0);
        std::vector<size_t> counts(k);
        std::fill(counts.begin(), counts.end(), 0);

        if (max_iters < std::numeric_limits<size_t>::max())
            max_iters++; // NOTE: This accounts for difference in pthread v. omp

        ret = kpmomp::compute_kmeans(&data[0], &centroids[0], &assignments[0],
                &counts[0], nrow, ncol, k, max_iters, nthread, init, tolerance,
                dist_type);
    } else {
#endif
        kpmprune::kmeans_task_coordinator::ptr kc =
            kpmprune::kmeans_task_coordinator::create(
                    datafn, nrow, ncol, k, max_iters, nnodes, nthread, p_centers,
                    init, tolerance, dist_type);
        ret = kc->run_kmeans();
#ifdef _OPEMP
    }
#endif

    // NOTE: the caller must take responsibility of cleaning up p_centers
    return ret;
}

} } // End namespace kpmeans::base
#endif
