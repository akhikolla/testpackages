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

#include <iostream>
#include <cassert>

#include "kmeans_thread.hpp"
#include "kmeans_types.hpp"
#include "util.hpp"
#include "io.hpp"
#include "clusters.hpp"

namespace kpmeans {
kmeans_thread::kmeans_thread(const int node_id, const unsigned thd_id,
        const unsigned start_rid,
        const unsigned nprocrows, const unsigned ncol,
        kpmbase::clusters::ptr g_clusters, unsigned* cluster_assignments,
        const std::string fn) : base_kmeans_thread(node_id, thd_id, ncol,
            g_clusters->get_nclust(), cluster_assignments, start_rid, fn) {

            this->nprocrows = nprocrows;
            this->g_clusters = g_clusters;
            local_clusters =
                kpmbase::clusters::create(g_clusters->get_nclust(), ncol);

            set_data_size(sizeof(double)*nprocrows*ncol);
#if VERBOSE
#ifndef
            std::cout << "Initializing thread. Metadata: thd_id: "
                << this->thd_id << ", start_rid: " << this->start_rid <<
                ", node_id: " << this->node_id << ", nprocrows: " <<
                this->nprocrows << ", ncol: " << this->ncol << std::endl;
#endif
#endif
        }

void kmeans_thread::sleep() {
    int rc;
    rc = pthread_mutex_lock(&mutex);
    if (rc) perror("pthread_mutex_lock");

    (*parent_pending_threads)--;
    set_thread_state(WAIT);

    if (*parent_pending_threads == 0) {
        rc = pthread_cond_signal(parent_cond); // Wake up parent thread
        if (rc) perror("pthread_cond_signal");
    }
    pthread_mutex_unlock(&mutex);
}

void kmeans_thread::run() {
    switch(state) {
        case TEST:
            test();
            break;
        case ALLOC_DATA:
            numa_alloc_mem();
            break;
        case KMSPP_INIT:
            kmspp_dist();
            break;
        case EM: /*E step of kmeans*/
            EM_step();
            break;
        case EXIT:
            throw kpmbase::thread_exception(
                    "Thread state is EXIT but running!\n");
        default:
            throw kpmbase::thread_exception("Unknown thread state\n");
    }
    sleep();
}

void kmeans_thread::wait() {
    int rc;
    rc = pthread_mutex_lock(&mutex);
    if (rc) perror("pthread_mutex_lock");

    while (state == WAIT) {
        //printf("Thread %d begin cond_wait\n", thd_id);
        rc = pthread_cond_wait(&cond, &mutex);
        if (rc) perror("pthread_cond_wait");
    }

    pthread_mutex_unlock(&mutex);
}

void kmeans_thread::wake(thread_state_t state) {
    int rc;
    rc = pthread_mutex_lock(&mutex);
    if (rc) perror("pthread_mutex_lock");
    set_thread_state(state);
    if (state == thread_state_t::KMSPP_INIT)
        cuml_dist = 0;
    rc = pthread_mutex_unlock(&mutex);
    if (rc) perror("pthread_mutex_unlock");

    rc = pthread_cond_signal(&cond);
}

void* callback(void* arg) {
    kmeans_thread* t = static_cast<kmeans_thread*>(arg);
#ifdef USE_NUMA
    t->bind2node_id();
#endif

    while (true) { // So we can receive task after task
        if (t->get_state() == WAIT)
            t->wait();

        if (t->get_state() == EXIT) {// No more work to do
            //printf("Thread %d exiting ...\n", t->thd_id);
            break;
        }

        //printf("Thread %d awake and doing a run()\n", t->thd_id);
        t->run(); // else
    }

    // We've stopped running so exit
    pthread_exit(NULL);

#ifdef _WIN32
    return NULL;
#endif
}

void kmeans_thread::start(const thread_state_t state=WAIT) {
    this->state = state;
    int rc = pthread_create(&hw_thd, NULL, callback, this);
    if (rc)
        throw kpmbase::thread_exception(
                "Thread creation (pthread_create) failed!", rc);
}

const unsigned kmeans_thread::
get_global_data_id(const unsigned row_id) const {
    return start_rid+row_id;
}

void kmeans_thread::EM_step() {
    meta.num_changed = 0; // Always reset at the beginning of an EM-step
    local_clusters->clear();

    for (unsigned row = 0; row < nprocrows; row++) {
        unsigned asgnd_clust = kpmbase::INVALID_CLUSTER_ID;
        double best, dist;
        dist = best = std::numeric_limits<double>::max();

        for (unsigned clust_idx = 0;
                clust_idx < g_clusters->get_nclust(); clust_idx++) {
            dist = kpmbase::dist_comp_raw<double>(&local_data[row*ncol],
                    &(g_clusters->get_means()[clust_idx*ncol]), ncol,
                    kpmbase::dist_type_t::EUCL);

            if (dist < best) {
                best = dist;
                asgnd_clust = clust_idx;
            }
        }

        assert(asgnd_clust != kpmbase::INVALID_CLUSTER_ID);
        unsigned true_row_id = get_global_data_id(row);

        if (asgnd_clust != cluster_assignments[true_row_id])
            meta.num_changed++;

        cluster_assignments[true_row_id] = asgnd_clust;
        local_clusters->add_member(&local_data[row*ncol], asgnd_clust);
    }
}

/** Method for a distance computation vs a single cluster.
 * Used in kmeans++ init
 */
void kmeans_thread::kmspp_dist() {
    unsigned clust_idx = meta.clust_idx;
    for (unsigned row = 0; row < nprocrows; row++) {
        unsigned true_row_id = get_global_data_id(row);

        double dist = kpmbase::dist_comp_raw<double>(&local_data[row*ncol],
                &((g_clusters->get_means())[clust_idx*ncol]), ncol,
                kpmbase::dist_type_t::EUCL);

        if (dist < dist_v[true_row_id]) { // Found a closer cluster than before
            dist_v[true_row_id] = dist;
            cluster_assignments[true_row_id] = clust_idx;
        }
        cuml_dist += dist_v[true_row_id];
    }
}

const void kmeans_thread::print_local_data() const {
    kpmbase::print_mat(local_data, nprocrows, ncol);
}
} // End namespace kpmeans
