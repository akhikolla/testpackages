
#ifndef ___ADFUNC_POOL_H___
#define ___ADFUNC_POOL_H___

#include "cppad.eigen.h"
#include "adfunc.h"
#include <vector>
#include <map>
#include <omp.h>



using CppAD::thread_alloc;

// used to inform CppAD when we are in parallel execution mode
static bool in_parallel(void)
{	
   return static_cast<bool>( omp_in_parallel() ); 
}
// ------------------------------------------------------------------
// used to inform CppAD of the current thread number
static size_t thread_number(void)
{	
   return static_cast<size_t>( omp_get_thread_num() ); 
}

template <class T>
TRIPLE(T) operator+(const TRIPLE(T)& A,const TRIPLE(T)& B)
{
  TRIPLE(T) result;
  // TODO - throw if triples are incompatible
  // determine number of ys
  unsigned int num_ys = boost::get<0>(A).size();
  // get number of dimensions
  unsigned int num_dimensions = 0;
  MULTIARG(T) y_star(num_ys);
  for(unsigned int i = 0; i < num_ys; i++)
    {
      num_dimensions += boost::get<0>(A)[i].rows()*boost::get<0>(A)[i].cols();
      y_star[i] = boost::get<0>(A)[i] + boost::get<0>(B)[i];
    }
  boost::get<0>(result) = y_star;
  boost::get<1>(result) = boost::get<1>(A) + boost::get<1>(B);
  TENSOR(T) hessians(1,num_dimensions);
  for(unsigned int i = 0; i < num_dimensions; i++)
    {
      hessians(0,i) = boost::get<2>(A)(0,i) + boost::get<2>(B)(0,i);
    }
  boost::get<2>(result) = hessians;
  return result;
}


template <class T>
TRIPLE(T) reduce_addition(const TRIPLE(T)& A, const TRIPLE(T)& B)
{ 
  return A + B;
  //  return sum(A,B);
}



template <class T>
TRIPLE(T) wrap_result_in_triple(const FUNCTION(T)& f,const MULTIARG(T)& xs)
{
  TRIPLE(T) result;
  boost::get<0>(result) = f(xs); 
  return result;
}


// this type of pool takes functions that use ADLA within them. Matlab style
// argument structure is used ie [y] f([x]) where [y] and [x] are of type MULTIARG 
// (ie multiple matrices) of type Scalar_Type
template <class Scalar_Type>
std::vector<TRIPLE(Scalar_Type)> 
adfunc_pool(const std::vector<boost::function<TRIPLE(Scalar_Type) (const MULTIARG(Scalar_Type)&)> >& fs, 
            const std::vector<MULTIARG(Scalar_Type)>& xs)
{
  // turn off dynamic thread adjustment
  omp_set_dynamic(0);

  unsigned int num_threads = omp_get_max_threads();
  unsigned int num_fs = fs.size();
  if(num_fs < num_threads)
    {
      num_threads = num_fs;
    }
  unsigned int fs_per_thread = num_fs / num_threads;
  // set the number of threads to use
  omp_set_num_threads(num_threads);
  // setup for using CppAD::AD<double> in parallel

  thread_alloc::parallel_setup(num_threads,in_parallel,thread_number);
  CppAD::parallel_ad<Scalar_Type>();

  std::vector<unsigned int> fs_in_thread(num_threads);
  for(unsigned int thread_id = 0; thread_id < num_threads; thread_id++)
    {
      fs_in_thread[thread_id] = fs_per_thread;
    } 
  for(unsigned int thread_id = 0; thread_id < num_fs%num_threads; thread_id++)
    {
      fs_in_thread[thread_id] += 1; 
    }
  std::multimap<unsigned int,unsigned int> thread_map;
 
 {
    unsigned int fs_id = 0;
    for(unsigned int thread_id = 0; thread_id < num_threads; thread_id++)
      {
	unsigned int num_fs_in_thread = fs_in_thread[thread_id];
	for(unsigned int count = 0; count < num_fs_in_thread; count++)
	  {
	    thread_map.insert(std::make_pair<unsigned int,unsigned int>(thread_id,fs_id));
	    fs_id++;
	  } 
      } 
  }
 
  std::vector<TRIPLE(Scalar_Type )> results(num_fs);
  // spin off the threads
  # pragma omp parallel for
  for(unsigned int thread_id = 0; thread_id < num_threads; thread_id++)
    {
      std::pair<std::multimap<unsigned int,unsigned int>::iterator,std::multimap<unsigned int,unsigned int>::iterator > thread_map_its;
      thread_map_its = thread_map.equal_range(thread_id);
      std::multimap<unsigned int,unsigned int>::iterator it = thread_map_its.first;
      while(it != thread_map_its.second)
	{
	  results[it->second] = fs[it->second](xs[it->second]);
	  it++;
	}
    }

  return results;

}

#endif
