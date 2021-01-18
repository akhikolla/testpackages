/* -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

   SFSMatrix.cpp
   This file is part of SFS.
   
   SFS is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   
   SFS is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with SFS.  If not, see <http://www.gnu.org/licenses/>.
   
   (C) 2014-2016 Matteo Seminaroti, CWI Amsterdam.
   R-integration (C) 2016 Utz-Uwe Haus, Cray EMEA Research Lab, Cray Inc.
*/

//#define DEBUG


#include "SFSMatrix.h"

#include <limits>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <time.h>


bool SFSMatrix::binary(void)
{
    IntVector binary;
    for (typename SpMat::iterator it = _A.begin();
         it != _A.end();
         ++it) {
      if (it.internal_col != it.row()) {
        if (binary.size() == 0 || (binary.size() == 1 && binary[0] != *it)) {
          binary.push_back(*it);
        }
        if (binary.size() == 2 && binary[0] != *it && binary[1] != *it) {
          return false;
        }
      }
    }
    
    return true;
}


// template <typename data_type>
typename SFSMatrix::IntVector 
SFSMatrix::solve()
{
    IntVector pi;
    WeakLinearOrder Phi;
    SpMat A = _A; //original matrix
    _t = clock();
    Robinson(pi);
    _t = clock() - _t;
    _A = A;
    
    bool one_based = false;
    //check if the nighborhood of 0 (dummy node) is empty
    if (_A.begin_col(0) == _A.end_col(0))
    {
        one_based = true;
    }
    
    //check if pi is a Robinson ordering
    if(_Robinsonian)
    {
        print_ordered_matrix(pi, _file_name);
    }
    
    //if the data is 1-based, remove from pi the element 0
    if(one_based)
    {
        for(int i=0; i<_n;++i)
        {
            if (pi[i] == 0)
            {
                pi.erase(pi.begin()+i);
                break;
            }
        }
    }
    
    return pi;
}

void SFSMatrix::Robinson(IntVector& pi)
{
    
    WeakLinearOrder CC = SFS();
    reorder_graph(CC,true);
    
    if (CC.size() > 1) {
      SpMat A;
      SFSMatrix G(A,_epsilon, _dissimilarity, _Robinsonian, _max_sweeps);
      for (block_iterator cc_it = CC.begin(); cc_it != CC.end(); ++cc_it)
      {
        subgraph(*cc_it, G);
        G.Robinson(pi);
      }
    } else {
        IntVector pi_opt;
        multisweep(pi_opt);
        for(int i =0; i < _n; ++i)
        {
            pi.push_back(pi_opt[i]);
        }
    }
    _coco = CC.size();
}

int SFSMatrix::multisweep(IntVector& pi_opt)
{
    IntVector sigma_loop_even (_n,-1), sigma_loop_odd (_n,-1), sigma_inv;
    WeakLinearOrder Phi;
    
    //if graph empty return immediately true
    if (_m==0 || _n<3)
    {
        for (int i = 0; i < _n; ++i)
        {
            pi_opt.push_back(_tau_inv[i].nominal);
        }
        return 0;
    }
    
    int max_sweep = _max_sweeps;
    for (int k = 1; k <= max_sweep; ++k)
    {
        
        Phi = SFS();
        sigma_inv = WeakLinearOrder_to_IntVector(Phi, true);
        reorder_graph(Phi, false);
        
        if ((k % 2) == 0) //even sweep
        {
            if (isReversed(sigma_inv,sigma_loop_even))
            {
                pi_opt = sigma_inv;
                return k+1;
            }
        }
        else //odd sweep
        {
            if (isReversed(sigma_inv,sigma_loop_odd))
            {
                pi_opt = sigma_inv;
                return k+1;
            }
        }
        reorder_graph(Phi, true);
    }
    //all sweeps needed
    pi_opt = sigma_inv;
    return max_sweep;
}

void SFSMatrix::reorder_graph(WeakLinearOrder& CC, bool reversed)
{
    permutation tau_inv(_n);
    uword i = 0;
    uword cc = 0;
    uword j;
    SpMat Pi (_n,_n);
    
    for(block_iterator cc_it = CC.begin(); cc_it != CC.end(); ++cc_it)
    {
        cc = 2 * i + cc_it->size();
        for(iterator it = cc_it->begin(); it != cc_it->end(); ++it)
        {
            if (reversed)
            {
                j = cc -1-i;
            }
            else
            {
                j = i;
            }
            tau_inv[j]=_tau_inv[*it];
            Pi (j,*it) = 1;
            *it = i; //reorder also CC (used later in subgraph)
            i++;
        }
    }
    _tau_inv = tau_inv;
    _A = Pi * _A * Pi.t();
}


void SFSMatrix::subgraph(Block& cc_it, SFSMatrix& G)
{
    //use subgraph only after that the graph is ordered according to CC
    uword i = cc_it.front();
    uword j = cc_it.back();
    
    //create submatrix
    if (cc_it.size() == (unsigned)_n)
    {
        G._A = _A;
    }
    else
    {
        G._A = _A.submat(i, i, j, j);
    }
    // update properties of submatrix
    G._n = G._A.n_rows;
    G._m = G._A.n_nonzero;
    G._tau_inv.resize(G._n);
    G._binary = _binary;
    G._epsilon = _epsilon;
    G._Robinsonian = _Robinsonian;
    G._max_sweeps = _max_sweeps;
    G._dissimilarity = _dissimilarity;
    for (int k = 0; k < G._n; ++k)
    {
        G._tau_inv[k] = _tau_inv[i];
        i++;
    }
}

SFSMatrix::WeakLinearOrder SFSMatrix::SFS()
{
    Queue Q;
    map_iterator map_queue_position; //map (vector of iterator) each vertex to the corresponding position in the block in the queue
    map_queue_position.resize(_n);
    map_queue_iterator map_queue; //map (vector of QueueIterator) each vertex to the corresponding block in the queue
    map_queue.resize(_n);
    Q.emplace_back();
    Q.back().second = 0;
    WeakLinearOrder CC, N;
    
    //initialization of the weak linear order as one block with n elements
    for (int i = 0; i < _n; ++i)
    {
        Q.back().first.push_back(i);
        map_queue[i] = Q.begin();
        map_queue_position[i]=prev(Q.begin()->first.end());
        _tau_inv[i].visited = false;
        _tau_inv[i].neighbor = false;
    }
    
    int p; //pivot
    for (int i = 0; i < _n; ++i)
    {
        //select the new pivot
        p = Q.begin()->first.front();
        if (!_tau_inv[p].neighbor) //p is the first vertex of the connected component
        {
            CC.emplace_back();
        }
        CC.back().push_back(p);
        _tau_inv[p].visited = true; //mark the pivot as visited
        
        Q.begin()->first.erase(map_queue_position[p]); //remove the pivot from the queue of unvisited vert
        if (Q.begin()->first.empty())
        {
            Q.erase(Q.begin());
        }
        
        N = Neighborhood(p);
        
        //refine the queue Q by the neighborhood N
        for (block_iterator n_it = N.begin(); n_it != N.end(); ++n_it)
        {
            partition_refinement(*n_it, Q, map_queue, map_queue_position);
        }
    }
    
    return CC;
}

SFSMatrix::WeakLinearOrder SFSMatrix::Neighborhood (int p)
{
    neighborhood L;
    WeakLinearOrder N;
    if(_binary)
    {
        N.emplace_back();
    }
    for(SpMat::iterator it = _A.begin_col(p); it != _A.end_col(p); ++it)
    {
        if (it.row() != (unsigned)p && !_tau_inv[it.row()].visited)
        {
            if (!_binary)
            {
                L.push_back(std::make_pair(it.row(),*it));
            }
            else
            {
                N.back().push_back(it.row());
            }
            _tau_inv[it.row()].neighbor = true;
        }
    }
    
    if(_binary)
    {
        return N;
    }
    
    if (!L.empty())
    {
        if (_dissimilarity)
        {
            L.sort(sort_dissimilarity());
        }
        else
        {
            L.sort(sort_similarity());
        }
        data_type b = std::numeric_limits<data_type>::max();
        for(neighborhood::iterator n_it = L.begin(); n_it != L.end(); ++n_it)
        {
            if (std::abs(b - n_it->second) > _epsilon)
            {
                N.emplace_back();
                b = n_it->second;
            }
            N.back().push_back(n_it->first);
        }
    }
    return N;
    
}

void SFSMatrix::partition_refinement(Block& S, Queue& Q, map_queue_iterator& map_queue, map_iterator& map_queue_position)
{
    Queue_visited_list W; //list of visited blocks
    
    //mark the blocks of Q intersectiong S
    for (iterator it = S.begin(); it != S.end(); ++it)
    {
        queue_iterator& Q_it = map_queue[*it]; //block of Q containing the vertex (*it)
        if (Q_it->second == 0)
        {
            W.push_back(Q_it);
        }
        Q_it->second++;
        iterator& v = map_queue_position[*it]; //position of vertex (*it) in the block of the queue
        Q_it->first.erase(v); //remove the vertex from its current position in the queue
        iterator u = Q_it->first.insert(std::next(Q_it->first.begin(),Q_it->second-1),*it); //insert the vertex at the beginning of the block in the queue
        map_queue_position[*it] = u;
    }
    
    //split each visited block of Q in two parts: intersection + difference
    for (Queue_visited_list::iterator w_it = W.begin(); w_it != W.end(); ++w_it)
    {
        queue_iterator& Q_it = *w_it;
        if (Q_it->second < Q_it->first.size())
        {
            queue_iterator Q_new = Q.emplace(Q_it); //create block Q_new left of Q_it
            for (unsigned int j = 0; j < Q_it->second; ++j)
            {
                Q_new->first.push_back(*Q_it->first.begin()); //add the first vertex in Q_it to the new block Q_new
                map_queue[*Q_it->first.begin()] = Q_new;//update map_queue for the vertex removed
                map_queue_position[*Q_it->first.begin()] = std::prev(Q_new->first.end());//update map_queue_position for the vertex removed
                Q_it->first.erase(Q_it->first.begin()); //remove the first vertex in Q_it
            }
        }
        Q_it->second = 0;
    }
}


SFSMatrix::IntVector SFSMatrix::WeakLinearOrder_to_IntVector(WeakLinearOrder& Phi, bool nominal)
{
    IntVector pi;
    for(block_iterator B_it = Phi.begin(); B_it != Phi.end(); ++B_it)
    {
        for(iterator it = B_it->begin(); it != B_it->end(); ++it)
        {
            if (nominal)
            {
                pi.push_back(_tau_inv[*it].nominal);
            }
            else
            {
                pi.push_back(*it);
            }
        }
    }
    return pi;
}


SFSMatrix::data_type SFSMatrix::isEpsilon_Robinson()
{
    IntVector row_i_1 (1); //dummy initialization
    data_type epsilon = 0;
    
    for (int i = 0; i < _n-1; ++i)
    {
        //create current row
        SpMat::iterator it = _A.end_col(i);
        int tau;
        if (it == _A.begin_col(i)) //avoid to check Robinson property for empty nodes
        {
            tau = i;
        }
        else{
            SpMat::iterator tmp = --it;
            tau = tmp.row();
        }
        int rmn = std::max(0,tau-i); //rightmost vertex adjacent to vertex _tau_inv[i]
        IntVector row_i(rmn,0); //create row for vertex i
        
        //check if the size of the current row is at least the size of the previous one
        if (row_i.size() < row_i_1.size())
        {
            row_i.resize(row_i_1.size(),0);
        }
        
        for(SpMat::iterator it = _A.begin_col(i); it != _A.end_col(i); ++it)
        {
            if (it.row() > (unsigned)i)
            {
                row_i[it.row() - i - 1] = *it;
            }
        }
        
        data_type nominal_value;
        //check if the row_i and row_i_1 are Robinson
        for (unsigned int j = row_i.size(); j-- > 0; /*no step*/)
        {
            nominal_value = row_i[j];
            if (j < row_i.size()-1)
            {
                if (_dissimilarity)
                {
                    row_i[j] = std::min(row_i[j],row_i[j+1]);
                }
                else
                {
                    row_i[j] = std::max(row_i[j],row_i[j+1]);
                }
                epsilon = std::max(epsilon, std::abs(row_i[j] - nominal_value));
            }
            if (i > 0 && j < row_i_1.size())
            {
                if (_dissimilarity)
                {
                    row_i[j] = std::min(row_i[j],row_i_1[j]);
                }
                else
                {
                    row_i[j] = std::max(row_i[j],row_i_1[j]);
                }
                epsilon = std::max(epsilon, std::abs(row_i[j] - nominal_value));
            }
        }
        
        if (row_i.size() > 0)
        {
            row_i.erase(row_i.begin());
        }
        row_i_1 = row_i; //update previous row for next iteration
    }
    
    return epsilon;
}

SFSMatrix::data_type SFSMatrix::check_isEpsilon_Robinson()
{
    data_type epsilon = 0;
    SpMat A = _A;
    
    for (int i = 0; i < _n-1; ++i)
    {
        for (int j = _n-1; j > i; --j)
        {
            if (j < _n-1 && A(i,j) < A(i,j+1))
            {
                epsilon = std::max(epsilon,A(i,j+1)-A(i,j));
                A(i,j) = A(i,j+1);
            }
            if (i > 0 && A(i,j) < A(i - 1,j))
            {
                epsilon = std::max(epsilon,A(i-1,j)-A(i,j));
                A(i,j) = A(i-1,j);
            }
        }
    }
    
    is_Robinson(A);
    return epsilon;
}

bool SFSMatrix::isReversed(IntVector& sigma_inv, IntVector& sigma_loop)
{
    bool reversed = true;
    
    //check is sigma_inv and sigma_loop coincide
    for (int i=0; i < _n; ++i)
    {
        if (sigma_inv[i] != sigma_loop[i])
        {
            reversed = false;
        }
        sigma_loop[i] = sigma_inv[i];
    }
    
    return reversed; //true if one of the two reversed is true
}


//debugging tools

void SFSMatrix::print_adjacency_list_tau(SpMat& A)
{
    int n = A.n_rows;
    for (int i=0; i < n; ++i)
    {
        SFSout << "node " << i << " edges ";
        for(SpMat::iterator it = A.begin_col(i); it != A.end_col(i); ++it)
        {
            SFSout << "(" << it.row() << "," << *it << ") ";
        }
        SFSout << std::endl;
    }
    SFSout << std::endl;
}

void SFSMatrix::print_adjacency_list_nominal()
{
    int n = _A.n_rows;
    for (int i=0; i < n; ++i)
    {
        SFSout << "node " << _tau_inv[i].nominal << " edges ";
        for(SpMat::iterator it = _A.begin_col(i); it != _A.end_col(i); ++it)
        {
            SFSout << "(" << _tau_inv[it.row()].nominal << "," << *it << ") ";
        }
        SFSout << std::endl;
    }
    SFSout << std::endl;
}


void SFSMatrix::check_tau()
{
    //check that the neighborhood of each vertex is ordered for increasing tau
    for (int i=0; i < _n; ++i)
    {
        unsigned int tau = 0;
        for(SpMat::iterator it = _A.begin_col(i); it != _A.end_col(i); ++it)
        {
            if (it.row() < tau) {
                throw std::runtime_error("the matrix is not consecutive in tau ");
            } else {
                tau = it.row();
            }
        }
    }
}


bool SFSMatrix::is_Robinson(SpMat& A)
{
    // check rows
    for (int i = 0; i < _n-1; ++i) {
        for (int j = i + 1; j < _n; ++j) {
            if (A(i,j - 1) < A(i,j)) {
                SFSout << "the matrix is not Robinson (rows)." << std::endl;
                SFSout << "A[" << i << "][" << j-1 << "] = "
                       << A(i,j-1) << " < " << "A[" << i << "][" << j
                       << "] = " << A(i,j) << std::endl;
                return false;
            }
        }
    }
    
    // check columns
    for (int j = 1; j < _n; ++j) {
        for (int i = 1; i < j; ++i) {
            if (A(i - 1,j) > A(i,j)) {
                SFSout << "the matrix is not Robinson (columns)." << std::endl;
                SFSout << "A[" << i-1 << "][" << j << "] = "
                       << A(i-1,j) << " > " << "A[" << i << "][" << j
                       << "] = " << A(i,j) << std::endl;
                return false;
            }
        }
    }
    return true;
}

bool SFSMatrix::is_permutation(IntVector& pi)
{
    std::vector<bool> perm (pi.size(),false);
    for (unsigned int i = 0; i < pi.size(); ++i) {
        perm[pi[i]] = true;
    }
    
    for (unsigned int i = 0; i < pi.size(); ++i) {
        if (!perm[i]) {
            SFSout << "Linear order is not a permutation" << std::endl;
            return false;
        }
    }
    
    if (pi.size() != (unsigned)_n) {
        SFSout << "permutation has a different size from the problem size"
               << std::endl;
        return false;
    }

    return true;
}

bool SFSMatrix::is_symmetric(SpMat& A)
{
    int n = A.n_rows;
    
    for (int i = 0; i < n - 1; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (A(i,j) != A(j,i)) {
                SFSout << "the matrix is not symmetric." << std::endl;
                SFSout << "A[" << i << "][" << j << "] = "
                       << A(i,j) << " != " << "A[" << j << "][" << i
                       << "] = " << A(j,i) << std::endl;
                return false;
            }
        }
    }
    return true;
}

void SFSMatrix::print_vertices_properties()
{
    for (int i = 0; i < _n; ++i)
    {
        SFSout << "node " << _tau_inv[i].nominal << ": ";
        SFSout << "tau " << i << ": ";
        SFSout << "visited = " << _tau_inv[i].visited << ", ";
        SFSout << "neighbor = " << _tau_inv[i].neighbor << ". ";
        SFSout << std::endl;
    }
    SFSout << std::endl;
}

void SFSMatrix::print_permutation (IntVector& pi, string& file)
{
    std::ofstream myfile;
    myfile.open("../../output/permutation_"+file);
    
    if (pi.empty())
    {
        for(int i = 0; i < _n; ++i)
        {
            pi.push_back(_tau_inv[i].nominal);
        }
    }
    else
    {
        assert(is_permutation(pi));
    }
    for (unsigned int i = 0; i < pi.size(); ++i)
    {
        myfile << pi[i] << "\n";
    }
    myfile.close();
}

void SFSMatrix::print_log (string& file)
{
    std::ofstream myfile;
    myfile.open("../../output/logs_"+file);
    
    //common part
    myfile << _A.n_rows << " = number of vertices \n";
    myfile << _A.n_nonzero << " = number of edges \n";
    myfile << _coco << " = number of connected components \n";
    myfile << ((float)_t)/(CLOCKS_PER_SEC) << " = time (seconds) \n";
    myfile.close();
    
}

void SFSMatrix::print_ordered_matrix (IntVector& pi, string& file)
{
    SpMat Pi (_n,_n);
    if (pi.empty())
    {
        for(int i = 0; i < _n; ++i)
        {
            pi.push_back(_tau_inv[i].nominal);
        }
    }
    
    
    for(int i = 0; i < _n; ++i)
    {
        Pi(i,pi[i])=1;
    }

    _A = Pi * _A * Pi.t();
    if (!_dissimilarity)
    {
        if (_Robinsonian && isEpsilon_Robinson() > _epsilon)
        {
            SFSout << "the matrix is not Robinsonian" << std::endl;
        }
        else{
            SFSout << "the matrix is Robinsonian" << std::endl;
        }
    }
    
//    std::ofstream myfile;
//    myfile.open("../../output/ordered_matrix_"+_file_name);
//    for(int i = 0; i < _n; ++i)
//    {
//        for(int j = 0; j < _n; ++j)
//        {
//            myfile << _A(i,j) << " ";
//        }
//        myfile << "\n";
//    }
//    myfile.close();
}
