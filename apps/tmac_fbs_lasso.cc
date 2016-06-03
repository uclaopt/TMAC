/***************************************************
 * Example 2: solving LASSO with FBS
 * 
 *    min lambda * ||x||_ 1 + 0.5 ||A'x - b||^2
 *
 * where A is a dense matrix, stored
 * in row major.
 **************************************************/

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
//using namespace MyAlgebra;
#include "algebra_namespace_switcher.h"

double objective(Vector& Atx, Vector& b, Vector& x, double lambda);

int main(int argc, char *argv[]) {

  // Step 0: Define the parameters and input file names  
  Params params;
  double lambda = 0.;
  string data_file_name;
  string label_file_name;
  // Step 1. Parse the input argument
  parse_input_argv_mm(&params, argc, argv, data_file_name, label_file_name, lambda);

  // Step 2. Load the data or generate synthetic data, define matained variables
  Matrix A;   // matrix is a row major matrix
  Vector b;
  loadMarket(A, data_file_name);
  loadMarket(b, label_file_name);

  // set para
  int problem_size = A.rows();
  params.problem_size = problem_size;
  params.tmac_step_size = 0.5;
  int sample_size = A.cols();
  // define auxilary variables
  Vector x(problem_size, 0.);   // unknown variables, initialized to zero
  Vector Atx(sample_size, 0.);  // maintained variables, initialized to zero
  

  // Step 3. Define your forward, or backward operators based on data and parameters
  double operator_step_size = 0.9;
  params.step_size = operator_step_size;
  // forward operator
  forward_grad_for_square_loss<Matrix> forward(&A, &b, &Atx, operator_step_size);
  using Forward=decltype( forward );

  // backward operator
  prox_l1 backward(operator_step_size, lambda);
  using Backward= decltype( backward );

  // Step 4. Define your operator splitting scheme
  // ForwardBackwardSplitting<forward_grad_for_least_square<Matrix>, prox_l1> fbs(&x, forward, backward);  
  ForwardBackwardSplitting<Forward, Backward> fbs(&x, forward, backward);  

  // Step 6. Call the TMAC function  
  double start_time = get_wall_time();
  TMAC(fbs, params);  
  double end_time = get_wall_time();

  
  // Step 7. Print results
  print_parameters(params);
  cout << "Objective value is: " << objective(Atx, b, x, lambda) << endl;
  cout << "Computing time  is: " << end_time - start_time << endl;
  cout << "---------------------------------" << endl;
  cout << "# of nonzero in x: " << norm(x, 0) << endl;
  cout << "||x||_2 =: " << norm(x, 2) << endl;
  cout << "---------------------------------" << endl;  
  return 0;
}

double objective(Vector& Atx, Vector& b, Vector& x, double lambda) {

  int len = b.size();
  Vector temp = Atx;
  add(temp, b, -1.);
  double norm_of_temp = norm(temp);
  return lambda * norm(x, 1) + 0.5 * norm_of_temp * norm_of_temp;

}
