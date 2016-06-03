/*****************************************************************
 * Example 3: solving l1 regularized logistic regression with FBS
 * 
 *    min lambda * ||x||_ 1 + sum_i log (1 + exp(- b_i * a_i' * x))
 *
 * where A is a sparse matrix with size num_features x num_samples.
 * 
 *****************************************************************/

#include <iostream>
#include "matrices.h"
#include "algebra.h"
#include "operators.h"
#include "parameters.h"
#include "splitting_schemes.h"
#include "tmac.h"
#include "util.h"
#include "MarketIO.h"
#include <thread>
using namespace std;
#include "algebra_namespace_switcher.h"

int main(int argc, char *argv[]) {

  // Step 0: Define the parameters and input file names  
  Params params;
  string data_file_name;
  double lambda = 0.;
  
  // Step 1. Parse the input argument
  parse_input_argv_libsvm(&params, argc, argv, data_file_name, lambda);

  // Step 2. Load the data or generate synthetic data, define matained variables
  SpMat A;   // A is a row major matrix
  Vector b;

  /*
  loadMarket(A, data_file_name);
  loadMarket(b, label_file_name);
  */

  loadLibSVM(A, b, data_file_name);
  
  // set para
  int problem_size = A.rows();
  params.problem_size = problem_size;
  params.tmac_step_size = 0.5;
  params.block_size = 100;
  params.worker_type = "random";
  // params.use_controller = true;
  int sample_size = A.cols();
  // define auxilary variables
  Vector x(problem_size, 0.);   // unknown variables, initialized to zero
  Vector Atx(sample_size, 0.);  // maintained variables, initialized to zero
  
  // Step 3. Define your forward, or backward operators based on data and parameters
  double operator_step_size = 0.1;
  params.step_size = operator_step_size;
  // forward operator
  forward_grad_for_log_loss<SpMat> forward(&A, &b, &Atx, operator_step_size);
  using Forward=decltype( forward );
  // backward operator
  prox_l1 backward(operator_step_size, lambda);
  using Backward= decltype( backward );

  // Step 4. Define your operator splitting scheme
  ForwardBackwardSplitting<Forward,Backward> fbs(&x, forward, backward);  
  double start_time = get_wall_time();
  // Step 6. Call the TMAC function  
  TMAC(fbs, params);  
  double end_time = get_wall_time();
  
  // Step 7. Print results
  print_parameters(params);
  cout << "Objective value is: " << l1_log_loss_objective(b, x, Atx, lambda) << endl;
  cout << "Computing time  is: " << end_time - start_time << endl;
  cout << "# of nonzero in x: " << norm(x, 0) << endl;
  cout << "||x||: " << norm(x) << endl;
  cout << "---------------------------------" << endl;
  return 0;
}


