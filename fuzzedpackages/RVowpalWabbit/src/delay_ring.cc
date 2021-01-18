#include <cstring>		// bzero()
#include <stdlib.h>
#include <pthread.h>
#include <iostream>
#include <float.h>
#include "v_array.h"
#include "example.h"
#include "global_data.h"
#include "parser.h"
#include "gd.h"

v_array<size_t> delay_indices;//thread specific state.

v_array<example*> delay_ring;//delay_ring state & mutexes
v_array<size_t> threads_to_use;
size_t local_index;
size_t global_index;
pthread_mutex_t delay = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t delay_empty = PTHREAD_COND_INITIALIZER;
pthread_cond_t delay_nonempty = PTHREAD_COND_INITIALIZER;
pthread_cond_t global_delay_nonempty = PTHREAD_COND_INITIALIZER;

size_t mesg = 0;

void initialize_delay_ring()
{
  if (global.local_prediction > 0 && (global.unique_id == 0 
				      || global.backprop
				      || global.corrective
				      || global.delayed_global))
    mesg = 1;
  size_t nt = global.num_threads()+mesg;
  if (global.backprop || global.delayed_global || global.corrective)
    nt += global.num_threads();
  reserve(delay_indices, nt);
  for (size_t i = 0; i < nt; i++)
    delay_indices[i] = 0;
  reserve(delay_ring, global.ring_size);
  for (size_t i = 0; i < global.ring_size; i++)
    delay_ring[i] = NULL;
  reserve(threads_to_use, global.ring_size);
  local_index = 0;
  global_index = 0;
}

void destroy_delay_ring()
{
  free(delay_indices.begin);
  delay_indices.begin = NULL;
  free(delay_ring.begin);
  delay_ring.begin = NULL;
  free(threads_to_use.begin);
  threads_to_use.begin = NULL;
}

bool thread_done(size_t thread)
{
  bool ret;
  if (!parser_done())
    return false;
  pthread_mutex_lock(&delay);
  ret = (delay_indices[thread] == local_index || global.delayed_global) 
    && (!(global.backprop || global.delayed_global || global.corrective) 
	|| global.local_prediction <= 0
	|| (delay_indices[thread+1+global.num_threads()] == global_index
	    && global_index == local_index));
  pthread_mutex_unlock(&delay);
  return ret;
}

example* return_example(size_t thread)
{
  size_t index = delay_indices[thread] % global.ring_size;
  example* ret = delay_ring[index];
  
  pthread_mutex_lock(&delay);
  delay_indices[thread]++;
  pthread_mutex_unlock(&delay);

  pthread_mutex_lock(&ret->lock);
  if (--threads_to_use[index] == 0) 
    {
      pthread_mutex_lock(&delay);
      delay_ring[index] = NULL;
      pthread_cond_broadcast(&delay_empty);
      pthread_mutex_unlock(&delay);
    }
  pthread_mutex_unlock(&ret->lock);
  return ret;
}

example* get_delay_example(size_t thread)
{//semiblocking
  size_t global_offset = mesg+global.num_threads();

  if ((global.backprop || global.corrective || global.delayed_global) 
      && global.local_prediction > 0)
    //we must do global training sometimes.
    if (local_index > delay_indices[thread+global_offset] + (global.ring_size >> 1)
	|| (parser_done() 
	    && (local_index == delay_indices[thread] || global.delayed_global)))
      //We want to do global training
      {
	if (delay_indices[thread+global_offset] == global_index && !parser_done())
	  {//global update unavailable.
	    pthread_mutex_lock(&delay);
	    while (delay_indices[thread+global_offset] == global_index)
	      pthread_cond_wait(&global_delay_nonempty, &delay);//Wait until global update available
	    pthread_mutex_unlock(&delay);
	  }
	if (delay_indices[thread+global_offset] < global_index)
	  {
	    example *ret = return_example(thread+global_offset);
	    ret->eta_round = ret->eta_global;
	    ret->final_prediction = ret->global_prediction;
	    return ret;
	  }
	else 
	  return NULL;
      }
  if (delay_indices[thread] != local_index && (!global.delayed_global || global.local_prediction <= 0))
    return return_example(thread);
  else
    return NULL;
}

example* blocking_get_delay_example(size_t thread)
{
  size_t index = delay_indices[thread] % global.ring_size;
  pthread_mutex_lock(&delay);
  while(delay_ring[index] == NULL)
    pthread_cond_wait(&delay_nonempty, &delay);
  pthread_mutex_unlock(&delay);
  return return_example(thread);
}

void delay_example(example* ex, size_t count)
{
  size_t delay_count = count+mesg;
  if ((global.backprop || global.corrective) 
      && global.local_prediction > 0)
    delay_count += count;

  if (delay_count == 0)
    {
      assert (ex->threads_to_finish == 0);
      output_and_account_example(ex);
      free_example(ex);
    }
  else
    {
      size_t index = local_index % global.ring_size;
      
      pthread_mutex_lock(&delay);
      while (delay_ring[index] != NULL)
	pthread_cond_wait(&delay_empty, &delay);

      delay_ring[index] = ex;
      threads_to_use[index] = delay_count;
      ex->threads_to_finish = delay_count;
      
      local_index++;
      if (count == 0)
	for (size_t i = 0; i < global.num_threads();i++)
	  delay_indices[i]++;

      pthread_cond_broadcast(&delay_nonempty);
      pthread_mutex_unlock(&delay);
    }
}

void delay_global_example(example* ex)
{ //only called by delayed backprop or delayed_global when 
    //there is training to do.
    pthread_mutex_lock(&delay);
    global_index++;
    pthread_cond_signal(&global_delay_nonempty);
    pthread_mutex_unlock(&delay);
}
