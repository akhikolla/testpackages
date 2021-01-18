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
#ifndef __KPM_KMEANS_TASK_COORDINATOR_HPP__
#define __KPM_KMEANS_TASK_COORDINATOR_HPP__

#include "base_kmeans_coordinator.hpp"
#include "util.hpp"

namespace kpmeans {
class task;

namespace base {
    class prune_clusters;
    class thd_safe_bool_vector;
}

namespace prune {
class kmeans_task_thread;
class dist_matrix;

class kmeans_task_coordinator : public kpmeans::base_kmeans_coordinator {
protected:
    // Metadata
    // max index stored within each threads partition
    std::vector<unsigned> thd_max_row_idx;
    std::shared_ptr<base::prune_clusters> cltrs;
    std::shared_ptr<base::thd_safe_bool_vector> recalculated_v;
    std::vector<double> dist_v; // global
    std::shared_ptr<dist_matrix> dm;

    kmeans_task_coordinator(const std::string fn, const size_t nrow,
            const size_t ncol, const unsigned k, const unsigned max_iters,
            const unsigned nnodes, const unsigned nthreads,
            const double* centers, const base::init_type_t it,
            const double tolerance, const base::dist_type_t dt);

public:
    static base_kmeans_coordinator::ptr create(
            const std::string fn, const size_t nrow,
            const size_t ncol, const unsigned k, const unsigned max_iters,
            const unsigned nnodes, const unsigned nthreads,
            const double* centers=NULL, const std::string init="kmeanspp",
            const double tolerance=-1, const std::string dist_type="eucl") {

        base::init_type_t _init_t = base::get_init_type(init);
        base::dist_type_t _dist_t = base::get_dist_type(dist_type);

#if KM_TEST
#ifndef BIND
        printf("kmeans task coordinator => NUMA nodes: %u, nthreads: %u, "
                "nrow: %lu, ncol: %lu, init: '%s', dist_t: '%s', fn: '%s'"
                "\n\n", nnodes, nthreads, nrow, ncol, init.c_str(),
                dist_type.c_str(), fn.c_str());
#endif
#endif
        return base_kmeans_coordinator::ptr(
                new kmeans_task_coordinator(fn, nrow, ncol, k, max_iters,
                    nnodes, nthreads, centers, _init_t, tolerance, _dist_t));
    }

    std::shared_ptr<base::prune_clusters> get_gcltrs() {
        return cltrs;
    }

    std::shared_ptr<dist_matrix> get_dm() {
        return dm;
    }
    virtual void kmeanspp_init() override;
    virtual void random_partition_init() override;
    virtual void forgy_init() override;
    virtual base::kmeans_t run_kmeans(double* allocd_data,
            const bool numa_opt) override;
    virtual ~kmeans_task_coordinator();
    virtual const void print_thread_data() override;
    virtual void build_thread_state() override;

    std::pair<size_t, size_t> get_rid_len_tup(const unsigned thd_id);
    // Pass file handle to threads to read & numa alloc
    void update_clusters(const bool prune_init);
    void wake4run(kpmeans::thread_state_t state) override;
    void destroy_threads() override;
    void set_thread_clust_idx(const unsigned clust_idx) override;
    double reduction_on_cuml_sum() override;
    void set_thd_dist_v_ptr(double* v) override;
    void run_init() override;
    void set_global_ptrs() override;
    void set_thread_data_ptr(double* allocd_data) override;
    const double* get_thd_data(const unsigned row_id) const override;
    void set_task_data_ptrs();
    void set_prune_init(const bool prune_init);
};
} } // End namespace kpmeans, prune
#endif
