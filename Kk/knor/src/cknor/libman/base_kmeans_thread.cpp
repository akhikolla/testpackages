/*
 * Copyright 2016 neurodata (http://neurodata.io/)
 * Written by Disa Mhembere (disa@jhu.edu)
 *
 * This file is part of knor
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

#ifdef USE_NUMA
#include <numa.h>
#endif

#include "base_kmeans_thread.hpp"
#include "exception.hpp"
#include "util.hpp"

#define VERBOSE 0
#define INVALID_THD_ID -1

namespace kpmeans {

void base_kmeans_thread::destroy_numa_mem() {
    if (!preallocd_data) {
#ifdef USE_NUMA
    numa_free(local_data, get_data_size());
#else
    delete [] local_data;
#endif
    }
}

void base_kmeans_thread::join() {
    void* join_status;
    int rc = pthread_join(hw_thd, &join_status);
    if (rc)
        throw base::thread_exception("pthread_join()", rc);

    thd_id = INVALID_THD_ID;
}

// Once the algorithm ends we should deallocate the memory we moved
void base_kmeans_thread::close_file_handle() {
    int rc = fclose(f);
    if (rc)
        throw base::io_exception("fclose() failed!", rc);

#if VERBOSE
#ifndef BIND
    printf("Thread %u closing the file handle.\n",thd_id);
#endif
#endif
    f = NULL;
}

// Move data ~equally to all nodes
void base_kmeans_thread::numa_alloc_mem() {
    kpmbase::assert_msg(f, "File handle invalid, can only alloc once!");
    size_t blob_size = get_data_size();
#ifdef USE_NUMA
    local_data = static_cast<double*>(numa_alloc_onnode(blob_size, node_id));
#else
    local_data = new double [blob_size/sizeof(double)];
#endif
    fseek(f, start_rid*ncol*sizeof(double), SEEK_SET); // start position
#ifdef NDEBUG
    size_t nread = fread(local_data, blob_size, 1, f);
    nread = nread + 1 - 1; // Silence compiler warning
#else
    assert(fread(local_data, blob_size, 1, f) == 1);
#endif
    close_file_handle();
}

void base_kmeans_thread::set_local_data_ptr(double* data, bool offset) {
    if (offset)
        local_data = &(data[start_rid*ncol]); // Grab your offset
    else
        local_data = data;
}

base_kmeans_thread::~base_kmeans_thread() {
    pthread_cond_destroy(&cond);
    pthread_mutex_destroy(&mutex);
    pthread_mutexattr_destroy(&mutex_attr);

    if (f)
        close_file_handle();
#if VERBOSE
#ifndef BIND
    printf("Thread %u being destroyed\n", thd_id);
#endif
#endif
    if (thd_id != INVALID_THD_ID)
        join();
}

void base_kmeans_thread::bind2node_id() {
#ifdef USE_NUMA
    struct bitmask *bmp = numa_allocate_nodemask();
    numa_bitmask_setbit(bmp, node_id);
    numa_bind(bmp);
    numa_free_nodemask(bmp);
#endif
    // No NUMA? Do nothing
}
} // End namespace kpmeans
