#include "barrier.h"
#include "algebra.h"
#include "operators.h"
#include "parameters.h"
#include "splitting_schemes.h"
#include "tmac.h"
#include <iostream>
#include "matrices.h"
#include "util.h"
#include <thread>
#include "algebra_namespace_switcher.h"
using namespace std;

/*
  This is the demo code for comparing sync-parallel with async-parallel
  for solving lasso:
     min lambda ||x||_1 + 0.5 ||A'x - b||^2
  where matrix A is a row major matrix.
 */

double objective(Vector& Atx, Vector& b, Vector& x, double lambda);

int main(int argc, char *argv[]) {

  Params params;
  string data_file_name;
  double lambda = 0.;  
  // parse data in libsvm format
  parse_input_argv_libsvm(&params, argc, argv, data_file_name, lambda);

  SpMat A;   // A is a row major matrix
  Vector b;
  // load the data
  loadLibSVM(A, b, data_file_name);
  
  int problem_size = A.rows();
  int sample_size = A.cols();  
  Vector x(problem_size, 0.);   // unknown variables, initialized to zero
  Vector Atx(sample_size, 0.);  // maintained variables, initialized to zero
  double operator_step_size = 0.1;

  // At is needed for sync-parallel
  SpMat At(sample_size, problem_size);
  transpose(A, At);
  
  params.problem_size = problem_size;
  params.tmac_step_size = 0.5;
  params.block_size = 100;
  params.worker_type = "random";  
  params.step_size = operator_step_size;
  print_parameters(params);
  
  prox_l1 backward(operator_step_size, lambda);
  forward_grad_for_square_loss<SpMat> forward(&A, &b, &Atx, operator_step_size, 1., &At);
  ForwardBackwardSplitting<forward_grad_for_square_loss<SpMat>, prox_l1> fbs(&x, forward, backward);

  // async timing
  params.async = true;
  double async_start_time = get_wall_time();  
  TMAC(fbs, params);
  double async_end_time = get_wall_time();
  double async_obj = objective(Atx, b, x, lambda);    
  
  // sync parallel needs both A and At
  // initialize the unknown and cache variables
  x = Vector(problem_size, 0.);
  Atx = Vector(sample_size, 0.);
  // sync timing
  params.async = false;
  params.step_size = 0.001;  // sync-parallel has to use smaller step size
  double sync_start_time = get_wall_time();
  TMAC(fbs, params);
  double sync_end_time = get_wall_time();
  double sync_obj = objective(Atx, b, x, lambda);

  std::cout << setw(15) << "method";
  std::cout << setw(15) << "time(s)";
  std::cout << setw(15) << "final obj" << endl;
  cout << setw(15) << "sync" << setw(15) << setprecision(5) << sync_end_time - sync_start_time << setw(15) << sync_obj << endl;
  cout << setw(15) << "async" << setw(15) << setprecision(5) << async_end_time - async_start_time << setw(15) << async_obj << endl;
  return 0;
}


double objective(Vector& Atx, Vector& b, Vector& x, double lambda) {
  return lambda * norm(x, 1) + 0.5 * square_loss(x, Atx, b);
}
