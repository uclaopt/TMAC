#ifndef AROCK_INCLUDE_WORKER_H
#define AROCK_INCLUDE_WORKER_H

#include "parameters.h"
#include "controller.h"
#include "range.h"

// asynchronous cyclic worker
template<typename Operator>
void async_cyclic_worker(Operator algorithm, Range range, Controller<Operator>& cont,  Params* params) {
  auto id = std::this_thread::get_id();
  int max_itrs = params->max_itrs;
  int num_cords = range.end - range.start;
  int idx = 0;
  
  algorithm.update_params(params);
  //alert controller to existence of worker
  if(params->use_controller) {
    cont.add_worker(id, algorithm, params);
  }

  for (int i = 0; i < max_itrs; i++) {
    for (int j = 0; j < num_cords; j++) {
      idx = j + range.start;
      auto fpr = algorithm(idx);
      
      if (params->use_controller) {
	cont.process_update(id, idx, fpr);
      }
    }
  }
  //finished task, so remove worker from controller
  cont.remove_worker(id);
  return;
}

// asynchronous Gauss Seidal worker
template<typename Operator>
void async_gs_worker(Operator algorithm, Controller<Operator>& cont, Params* params) {
  auto id = std::this_thread::get_id();
  int total_num_threads = params->total_num_threads;
  bool use_controller = params->use_controller;
  int num_workers = use_controller ? total_num_threads - 1 : total_num_threads;  
  int max_itrs = params->max_itrs;
  int problem_size = params->problem_size;

  algorithm.update_params(params);
  //alert controller to existence of worker
  if(params->use_controller) {
    cont.add_worker(id, algorithm, params);
  }
  
  for (int itr = 0; itr < max_itrs / num_workers; itr++) {
    // running Gauss-Seidal update in parallel
    for (int idx = 0; idx < problem_size; idx++) {
      // algorithm.update_params(params);
      auto fpr = algorithm(idx);
      if (params->use_controller) {
	cont.process_update(id, idx, fpr);
      }
    }
  }
  cont.remove_worker(id);  
  return;
}

// asynchronous random block coordinate update worker
template<typename Operator>
void async_rnd_worker(Operator algorithm, Controller<Operator>& cont, Params* params) {
  auto id = std::this_thread::get_id();
  int total_num_threads = params->total_num_threads;
  bool use_controller = params->use_controller;
  int num_workers = use_controller ? total_num_threads - 1 : total_num_threads;  
  int max_itrs = params->max_itrs;
  int idx = 0;
  int problem_size = params->problem_size;
  int blk_size = params->block_size;
  int num_blks = problem_size / blk_size;
  // assign local_num_blks to each worker so that
  // we can quantify one epoch
  int local_num_blks = num_blks / num_workers;
  int local_start = 0, local_end = 0;
  int blk_id = 0;

  algorithm.update_params(params);
  //alert controller to existence of worker
  if(params->use_controller) {
    cont.add_worker(id, algorithm, params);
  }

  
  for (int itr = 0; itr < max_itrs; itr++) {
    for (int blk = 0; blk < local_num_blks; blk++) {
      // randomly generate a block id
      blk_id = rand() % num_blks;
      // calculate the starting index and the ending index
      local_start = blk_id * blk_size;
      local_end = ( blk_id + 1 ) * blk_size;
      
      if ( blk_id == num_blks - 1 ) {
        local_end = problem_size;
      }
      
      for (int i = local_start; i < local_end; i++) {
        idx = i;
        auto fpr = algorithm(idx);
        if (params->use_controller) {
          cont.process_update(id, idx, fpr);
        }
      }
    }
  }
  cont.remove_worker(id);  
  return;
}


#endif
