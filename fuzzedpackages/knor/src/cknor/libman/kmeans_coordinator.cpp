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

#include <random>
#include <stdexcept>

#include "kmeans_coordinator.hpp"
#include "kmeans_thread.hpp"
#include "io.hpp"
#include "clusters.hpp"

namespace kpmeans {
kmeans_coordinator::kmeans_coordinator(const std::string fn, const size_t nrow,
        const size_t ncol, const unsigned k, const unsigned max_iters,
        const unsigned nnodes, const unsigned nthreads,
        const double* centers, const kpmbase::init_type_t it,
        const double tolerance, const kpmbase::dist_type_t dt) :
    base_kmeans_coordinator(fn, nrow, ncol, k, max_iters,
            nnodes, nthreads, centers, it, tolerance, dt) {

        cltrs = kpmbase::clusters::create(k, ncol);
        if (centers) {
            if (kpmbase::init_type_t::NONE)
                cltrs->set_mean(centers);
            else {
#ifndef BIND
                printf("[WARNING]: Both init centers "
                        "provided & non-NONE init method specified\n");
#endif
            }
        }
        build_thread_state();
    }

void kmeans_coordinator::build_thread_state() {
    // NUMA node affinity binding policy is round-robin
    unsigned thds_row = nrow / nthreads;
    for (unsigned thd_id = 0; thd_id < nthreads; thd_id++) {
        std::pair<unsigned, unsigned> tup = get_rid_len_tup(thd_id);
        thd_max_row_idx.push_back((thd_id*thds_row) + tup.second);
        threads.push_back(kmeans_thread::create((thd_id % nnodes),
                    thd_id, tup.first, tup.second,
                    ncol, cltrs, &cluster_assignments[0], fn));
        threads[thd_id]->set_parent_cond(&cond);
        threads[thd_id]->set_parent_pending_threads(&pending_threads);
        threads[thd_id]->start(WAIT); // Thread puts itself to sleep
    }
}

std::pair<unsigned, unsigned>
kmeans_coordinator::get_rid_len_tup(const unsigned thd_id) {
    unsigned rows_per_thread = nrow / nthreads;
    unsigned start_rid = (thd_id*rows_per_thread);

    if (thd_id == nthreads - 1)
        rows_per_thread += nrow % nthreads;
    return std::pair<unsigned, unsigned>(start_rid, rows_per_thread);
}

void kmeans_coordinator::destroy_threads() {
    wake4run(EXIT);
}

// <Thread, within-thread-row-id>
const double* kmeans_coordinator::get_thd_data(const unsigned row_id) const {
    // TODO: Cheapen
    unsigned parent_thd = std::upper_bound(thd_max_row_idx.begin(),
            thd_max_row_idx.end(), row_id) - thd_max_row_idx.begin();
    unsigned rows_per_thread = nrow/nthreads; // All but the last thread

    return &((threads[parent_thd]->get_local_data())
            [(row_id-(parent_thd*rows_per_thread))*ncol]);
}

void kmeans_coordinator::update_clusters() {
    num_changed = 0; // Always reset here since there's no pruning
    cltrs->clear();

    // Serial aggreate of OMP_MAX_THREADS vectors
    for (thread_iter it = threads.begin(); it != threads.end(); ++it) {
        // Updated the changed cluster count
        num_changed += (*it)->get_num_changed();
        // Summation for cluster centers

        cltrs->peq((*it)->get_local_clusters());
    }

    unsigned chk_nmemb = 0;
    for (unsigned clust_idx = 0; clust_idx < k; clust_idx++) {
        cltrs->finalize(clust_idx);
        cluster_assignment_counts[clust_idx] =
            cltrs->get_num_members(clust_idx);
        chk_nmemb += cluster_assignment_counts[clust_idx];
    }
    if (chk_nmemb != nrow)

    assert(chk_nmemb == nrow);
    assert(num_changed <= nrow);
}

double kmeans_coordinator::reduction_on_cuml_sum() {
    double tot = 0;
    for (thread_iter it = threads.begin(); it != threads.end(); ++it)
        tot += (*it)->get_cuml_dist();
    return tot;
}

void kmeans_coordinator::wake4run(const thread_state_t state) {
    pending_threads = nthreads;
    for (unsigned thd_id = 0; thd_id < threads.size(); thd_id++)
        threads[thd_id]->wake(state);
}

void kmeans_coordinator::set_thread_clust_idx(const unsigned clust_idx) {
    for (thread_iter it = threads.begin(); it != threads.end(); ++it)
        (*it)->set_clust_idx(clust_idx);
}

void kmeans_coordinator::set_thd_dist_v_ptr(double* v) {
    for (unsigned thd_id = 0; thd_id < threads.size(); thd_id++) {
        pthread_mutex_lock(&mutex);
        threads[thd_id]->set_dist_v_ptr(v);
        pthread_mutex_unlock(&mutex);
    }
}

void kmeans_coordinator::kmeanspp_init() {
    struct timeval start, end;
    gettimeofday(&start , NULL);

    std::vector<double> dist_v;
    dist_v.assign(nrow, std::numeric_limits<double>::max());
    set_thd_dist_v_ptr(&dist_v[0]);

    std::default_random_engine generator;
    std::uniform_int_distribution<unsigned> distribution(0, nrow-1);

    // Choose c1 uniformly at random
    unsigned selected_idx = distribution(generator);
    cltrs->set_mean(get_thd_data(selected_idx), 0);
    dist_v[selected_idx] = 0.0;
    cluster_assignments[selected_idx] = 0;

    unsigned clust_idx = 0; // The number of clusters assigned

    std::uniform_real_distribution<double> ur_distribution(0.0, 1.0);

    // Choose next center c_i with weighted prob
    while (true) {
        set_thread_clust_idx(clust_idx); // Set the current cluster index
        wake4run(KMSPP_INIT); // Run || distance comp to clust_idx
        wait4complete();
        double cuml_dist = reduction_on_cuml_sum(); // Sum the per thread cumulative dists

        cuml_dist = (cuml_dist * ur_distribution(generator)) / (RAND_MAX - 1.0);
        if (++clust_idx >= k)  // No more centers needed
            break;

        for (size_t row = 0; row < nrow; row++) {
            cuml_dist -= dist_v[row];
            if (cuml_dist <= 0) {
                cltrs->set_mean(get_thd_data(row), clust_idx);
                cluster_assignments[row] = clust_idx;
                break;
            }
        }
        assert(cuml_dist <= 0);
    }

    gettimeofday(&end, NULL);
}

void kmeans_coordinator::random_partition_init() {
    std::default_random_engine generator;
    std::uniform_int_distribution<unsigned> distribution(0, k-1);

    for (unsigned row = 0; row < nrow; row++) {
        unsigned asgnd_clust = distribution(generator);
        const double* dp = get_thd_data(row);

        cltrs->add_member(dp, asgnd_clust);
        cluster_assignments[row] = asgnd_clust;
    }

    cltrs->finalize_all();
}

void kmeans_coordinator::forgy_init() {
    std::default_random_engine generator;
    std::uniform_int_distribution<unsigned> distribution(0, nrow-1);

    for (unsigned clust_idx = 0; clust_idx < k; clust_idx++) { // 0...k
        unsigned rand_idx = distribution(generator);
        cltrs->set_mean(get_thd_data(rand_idx), clust_idx);
    }
}

void kmeans_coordinator::run_init() {
    switch(_init_t) {
        case kpmbase::init_type_t::RANDOM:
            random_partition_init();
            break;
        case kpmbase::init_type_t::FORGY:
            forgy_init();
            break;
        case kpmbase::init_type_t::PLUSPLUS:
            kmeanspp_init();
            break;
        case kpmbase::init_type_t::NONE:
            break;
        default:
            throw std::runtime_error("Unknown initialization type");
    }
}

/**
 * Main driver for kmeans
 */
kpmbase::kmeans_t kmeans_coordinator::run_kmeans(
        double* allocd_data, const bool numa_opt=false) {
#ifdef PROFILER
    ProfilerStart("libman/kmeans_coordinator.perf");
#endif

    if (!numa_opt && NULL == allocd_data) {
        wake4run(ALLOC_DATA);
        wait4complete();
    } else if (allocd_data) { // No NUMA opt
        set_thread_data_ptr(allocd_data);
    } // Do nothing for numa_opt .. done in binding/knori.hpp

    struct timeval start, end;
    gettimeofday(&start , NULL);
    run_init(); // Initialize clusters

    // Run kmeans loop
    bool converged = false;
    size_t iter = 0;

    if (max_iters > 0)
        iter++;

    while (iter <= max_iters && max_iters > 0) {
        if (iter == 1)
            clear_cluster_assignments();

        wake4run(EM);
        wait4complete();

        update_clusters();

#if VERBOSE
#ifndef BIND
        printf("Cluster assignment counts: \n");
#endif
        kpmbase::print_vector(cluster_assignment_counts);
#endif

        if (num_changed == 0 ||
                ((num_changed/(double)nrow)) <= tolerance) {
            converged = true;
            break;
        }
        iter++;
    }
#ifdef PROFILER
    ProfilerStop();
#endif

    gettimeofday(&end, NULL);
#ifndef BIND
    printf("\n\nAlgorithmic time taken = %.6f sec\n",
        kpmbase::time_diff(start, end));
    printf("\n******************************************\n");
#endif
    if (converged) {
#ifndef BIND
        printf("K-means converged in %lu iterations\n", iter);
#endif
    } else {
#ifndef BIND
        printf("[Warning]: K-means failed to converge in %lu"
            " iterations\n", iter);
#endif
    }

#ifndef BIND
    printf("Final cluster counts: \n");
    kpmbase::print_vector(cluster_assignment_counts);
    printf("\n******************************************\n");
#endif

    return kpmbase::kmeans_t(this->nrow, this->ncol, iter, this->k,
            &cluster_assignments[0], &cluster_assignment_counts[0],
            cltrs->get_means());
}

kmeans_coordinator::~kmeans_coordinator() {
    thread_iter it = threads.begin();
    for (; it != threads.end(); ++it)
        (*it)->destroy_numa_mem();

    pthread_cond_destroy(&cond);
    pthread_mutex_destroy(&mutex);
    pthread_mutexattr_destroy(&mutex_attr);
    destroy_threads();
}

void const kmeans_coordinator::print_thread_data() {
    thread_iter it = threads.begin();
    for (; it != threads.end(); ++it) {
#ifndef BIND
        std::cout << "\nThd: " << (*it)->get_thd_id() << std::endl;
#endif
        (*it)->print_local_data();
    }
}

// Testing
void const kmeans_coordinator::print_thread_start_rids() {
    thread_iter it = threads.begin();
    for (; it != threads.end(); ++it) {
#ifndef BIND
        printf("\nThd: %u, start_rid: %lu\n", (*it)->get_thd_id(),
            (*it)->get_start_rid());
#endif
    }
}
}
