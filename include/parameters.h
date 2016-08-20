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
  double feastol;
  double unbdd;
  double epsilon;
  double infea;
  double get_step_size() {
    return step_size;
  }

  double get_tmac_step_size() {
    return tmac_step_size;
  }

  int get_problem_dimension() {
    return problem_size;
  }
  
  Params() : tmac_step_size(0.618), problem_size(0), block_size(1), use_controller(true), feastol(1e-5), unbdd(1e-10), epsilon(0.5), infea(5e-6),
    max_itrs(100), worker_type("cyclic"), async(true) {}
  
};
#endif
