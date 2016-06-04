#ifndef TMAC_INCLUDE_PARAMETERS_H
#define TMAC_INCLUDE_PARAMETERS_H

struct Params {

  double step_size;
  int max_itrs;
  int total_num_threads;
  bool use_controller;
  double tmac_step_size;
  int block_size;
  int problem_size;
  string worker_type;
  string step_size_rule;
  bool async;
  double get_step_size() {
    return step_size;
  }

  double get_tmac_step_size() {
    return tmac_step_size;
  }

  int get_problem_dimension() {
    return problem_size;
  }
  
 Params() : tmac_step_size (0.618), problem_size(0), block_size(1), use_controller(false),
    max_itrs(10), worker_type("cyclic"), async(true) {}
  
};
#endif
