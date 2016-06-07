#ifndef TMAC_INCLUDE_DRIVER_H
#define TMAC_INCLUDE_DRIVER_H


#include "matrices.h"
#include "algebra.h"
#include "worker.h"
#include "controller.h"
#include "parameters.h"
#include <thread>


template<typename Splitting>
void AROCK (Splitting op, Params parameters, Controller<Splitting> controller = Controller<Splitting>()) {
  int total_num_threads = parameters.total_num_threads;
  bool use_controller = parameters.use_controller;
  std::vector<std::thread> mythreads;
  int problem_size = parameters.get_problem_dimension();
  string worker_type = parameters.worker_type;
  int num_workers = use_controller ? total_num_threads - 1 : total_num_threads;
  
  for (size_t i = 0; i < num_workers; i++) {
    if (worker_type == "cyclic") {
      // partition the indexes
      int block_size = problem_size / num_workers;
      Range range(i * block_size, (i + 1) * block_size);
      if (i == num_workers - 1) {
        range.end = problem_size;
      }
      // cyclic coordinate update
      mythreads.push_back(std::thread(async_cyclic_worker<Splitting>, op, range, std::ref(controller), &parameters));
      
    } else if (worker_type == "gs") {
      // async Gauss Seidal rule
      mythreads.push_back(std::thread(async_gs_worker<Splitting>, op, std::ref(controller), &parameters));

    } else if (worker_type == "random") {
      // random block coordinate update
      mythreads.push_back(std::thread(async_rnd_worker<Splitting>, op, std::ref(controller), &parameters));
    }
  }
 
  if (use_controller) {
    mythreads.push_back(std::thread(Controller_loop<Splitting>, std::ref(controller)));
  }
  for (size_t i = 0; i < total_num_threads; i++) {
    mythreads[i].join();
  }
}


template<typename Splitting>
void SYNC(Splitting op, Params parameters) {
  
  int total_num_threads = parameters.total_num_threads;
  std::vector<std::thread> mythreads;
  int num_workers = total_num_threads;
  Barrier computation_barrier(num_workers);
  Barrier update_barrier(num_workers);
  Barrier cache_update_barrier(num_workers);
  int problem_size = parameters.get_problem_dimension();
  int block_size = problem_size / num_workers;
  
  for (size_t i = 0; i < num_workers; i++) {
    // partition the indexes
    Range range(i * block_size, (i + 1) * block_size);
    if (i == num_workers - 1) {
      range.end = problem_size;
    }
    mythreads.push_back(std::thread(sync_par_worker<Splitting>, op, range,
                                    &parameters,
                                    std::ref(computation_barrier),
                                    std::ref(update_barrier),
                                    std::ref(cache_update_barrier)
                                    ));
  }
  // join the threads
  for (size_t i = 0; i < total_num_threads; i++) {
    mythreads[i].join();
  }
}




#endif
