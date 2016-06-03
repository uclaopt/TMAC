/**********************************
 * Example: minimizing huber loss
 * 
 *    min weight*huber(A'x - b)
 *
 * with huber(y)=\sum_i f(y_i)
 * where f(w) = 0.5 * weight * w^2, if -delta <= w <= delta;
 *              delta * weight * (|w| - 0.5 * delta), otherwise.
 *
 * A is a dense matrix, stored
 * in row major.
 *********************************/

#include <iostream>
#include "matrices.h"
#include "algebra.h"
#include "operators.h"
#include "parameters.h"
#include "splitting_schemes.h"
#include "tmac.h"
#include "util.h"
#include <thread>
using namespace std;
//using namespace MyAlgebra;
#include "algebra_namespace_switcher.h"

double objective(Vector& Atx, Vector& b, double delta);

int main(int argc, char *argv[]) {

  // Step 0: Define the parameters and input file names  
  Params params;
  string data_file_name;
  string label_file_name;
  
  // Step 1. Parse the input argument
  parse_input_argv_mm(&params, argc, argv, data_file_name, label_file_name);

  // Step 2. Load the data or generate synthetic data, define matained variables
  Matrix A;   // matrix is a row major matrix
  Vector b;
  loadMarket(A, data_file_name);
  loadMarket(b, label_file_name);
  int problem_size = A.rows();
  params.problem_size = problem_size;
  params.worker_type="gs";
  params.tmac_step_size = 0.5;
  int sample_size = A.cols();
  Vector x(problem_size, 0.);   // unknown variables, initialized to zero
  Vector Atx(sample_size, 0.);  // maintained variables, initialized to zero
  
  // Step 3. Define your forward, or backward operators based on data and parameters
  double forward_operator_step_size = 0.5;
  double weight = 1.;
  double delta = 1.;
  params.step_size = forward_operator_step_size;
  
  forward_grad_for_huber_loss<Matrix> forward(&A, &b, &Atx, forward_operator_step_size, weight, delta);  
  
  // Step 4. Define your operator splitting scheme
  GradientDescentAlgorithm<forward_grad_for_huber_loss<Matrix> > gd(&x, forward);

  // Step 6. Call the TMAC function
  double start_time = get_wall_time();
  TMAC(gd, params);  
  double end_time = get_wall_time();
  print_parameters(params);
  
  cout << "Computing time is: " << end_time - start_time << endl;  
  // Step 7. Print results

  cout << "Objective value is: " << objective(Atx, b, delta) << endl;
  cout << "---------------------------------" << endl;  
  return 0;
}

double objective(Vector& Atx, Vector& b, double delta) {

  int len = b.size();
  Vector temp = Atx;
  double sum = 0;
  add(temp, b, -1.);
  for (int i = 0;i < len; ++i) {
    if (temp[i]<=delta && temp[i]>=-delta){
      sum += 0.5*temp[i]*temp[i];
    }
    else {
      sum += delta*(abs(temp[i])-0.5*delta);
    }
  }
  return sum;
}
