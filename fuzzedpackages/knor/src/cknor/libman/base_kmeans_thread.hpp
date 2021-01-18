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

#ifndef __KPM_BASE_KMEANS_THREAD_HPP__
#define __KPM_BASE_KMEANS_THREAD_HPP__

#include <pthread.h>

#include <memory>
#include <utility>
#include <atomic>
#include <string>
#include <cassert>

#include "thread_state.hpp"
#include "exception.hpp"

#define VERBOSE 0
#define INVALID_THD_ID -1

namespace kpmeans {
class task_queue;

namespace base {
    class clusters;
    class thd_safe_bool_vector;
}

namespace prune {
    class dist_matrix;
}
}

namespace kpmbase = kpmeans::base;
namespace kpmprune = kpmeans::prune;

namespace kpmeans {

union metaunion {
    unsigned num_changed; // Used during kmeans
    unsigned clust_idx; // Used during kms++
};

class base_kmeans_thread {
protected:
    pthread_t hw_thd;
    unsigned node_id; // Which NUMA node are you on?
    int thd_id;
    size_t start_rid; // With respect to the original data
    size_t ncol; // How many columns in the data
    double* local_data; // Pointer to where the data begins that the thread works on
    size_t data_size; // true size of local_data at any point
    std::shared_ptr<kpmbase::clusters> local_clusters;

    pthread_mutex_t mutex;
    pthread_cond_t cond;
    pthread_mutexattr_t mutex_attr;

    pthread_cond_t* parent_cond;
    std::atomic<unsigned>* parent_pending_threads;

    metaunion meta;
    //unsigned num_changed;

    FILE* f; // Data file on disk
    unsigned* cluster_assignments;

    kpmeans::thread_state_t state;
    double* dist_v;
    double cuml_dist;
    bool preallocd_data; // Is our data pre-allocated?

    friend void* callback(void* arg);

    base_kmeans_thread(const int node_id, const unsigned thd_id,
            const unsigned ncol, const unsigned nclust,
            unsigned* cluster_assignments, const unsigned start_rid,
            const std::string fn="") {

        pthread_mutexattr_init(&mutex_attr);
        pthread_mutexattr_settype(&mutex_attr, PTHREAD_MUTEX_ERRORCHECK);
        pthread_mutex_init(&mutex, &mutex_attr);
        pthread_cond_init(&cond, NULL);
        this->node_id = node_id;
        this->thd_id = thd_id;
        this->ncol = ncol;
        this->cluster_assignments = cluster_assignments;
        this->start_rid = start_rid;
        this->preallocd_data = false;

        if (!fn.empty()) {
            this->f = fopen(fn.c_str(), "rb");
            assert(this->f);
        }
        else {
            this->f = NULL;
            preallocd_data = true;
        }

        meta.num_changed = 0; // Same as meta.clust_idx = 0;
        set_thread_state(WAIT);
    }

    void set_thread_state(kpmeans::thread_state_t state) {
        this->state = state;
    }

public:
    typedef std::shared_ptr<base_kmeans_thread> ptr;

    virtual void start(const kpmeans::thread_state_t state) = 0;
    // Allocate and move data using this thread
    virtual void EM_step() = 0;
    virtual void kmspp_dist() = 0;
    virtual const unsigned get_global_data_id(const unsigned row_id) const = 0;
    virtual void run() = 0;
    virtual void sleep() = 0;

    // Used by `task' thread classes
    virtual void set_driver(void* driver) {
        throw kpmbase::abstract_exception();
    }
    virtual void wake(kpmeans::thread_state_t state) {
        throw kpmbase::abstract_exception();
    }
    virtual void set_prune_init(const bool prune_init) {
        throw kpmbase::abstract_exception();
    }
    virtual void set_recalc_v_ptr(std::shared_ptr<kpmbase::thd_safe_bool_vector>
            recalculated_v) {
        throw kpmbase::abstract_exception();
    }
    virtual void set_dist_mat_ptr(std::shared_ptr<kpmprune::dist_matrix> dm) {
        throw kpmbase::abstract_exception();
    }
    virtual bool try_steal_task() { throw kpmbase::abstract_exception(); }
    virtual task_queue* get_task_queue() {
        throw kpmbase::abstract_exception();
    }
    virtual const void print_local_data() const {
        throw kpmbase::abstract_exception();
    };

    void test() {
    }

    void set_dist_v_ptr(double* v) {
        dist_v = v;
    }

    const kpmeans::thread_state_t get_state() const {
        return this->state;
    }

    const unsigned get_thd_id() const {
        return thd_id;
    }

    const double* get_local_data() const {
        return local_data;
    }

    const unsigned get_num_changed() const {
        return meta.num_changed;
    }

    const std::shared_ptr<kpmbase::clusters> get_local_clusters() const {
        return local_clusters;
    }

    void set_clust_idx(const unsigned idx) {
        meta.clust_idx = idx;
    }

    const double get_cuml_dist() const {
        return cuml_dist;
    }

    void set_data_size(const size_t data_size) {
        this->data_size = data_size;
    }

    const size_t get_data_size() const {
        return this->data_size;
    }

    pthread_mutex_t& get_lock() {
        return mutex;
    }

    pthread_cond_t& get_cond() {
        return cond;
    }

    unsigned get_node_id() {
        return node_id;
    }

    void set_parent_cond(pthread_cond_t* cond) {
        parent_cond = cond;
    }

    void set_parent_pending_threads(std::atomic<unsigned>* ppt) {
        parent_pending_threads = ppt;
    }

    void destroy_numa_mem();

    const size_t get_start_rid() const {
        return start_rid;
    }

    void set_start_rid(size_t start_rid) {
        this->start_rid = start_rid;
    }

    void join();
    // Once the algorithm ends we should deallocate the memory we moved
    void close_file_handle();
    // Move data ~equally to all nodes
    void numa_alloc_mem();
    void set_local_data_ptr(double* data, bool offset=true);
    ~base_kmeans_thread();
    void bind2node_id();
};
}
#endif
