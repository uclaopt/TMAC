/**********************************
 * Example 1: solving least squarex
 * 
 *    min 0.5 ||A'x - b||^2
 * where A is a dense matrix, stored
 * in row major.
 *********************************/

#include <iostream>
#include "matrices.h"
#include "algebra.h"
#include "operators.h"
#include "parameters.h"
#include "splitting_schemes.h"
#include "arock.h"
#include "util.h"
#include <thread>
using namespace std;
//using namespace MyAlgebra;
#include "algebra_namespace_switcher.h"

double objective(Vector& Atx, Vector& b);

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
  params.arock_step_size = 0.1;
  int sample_size = A.cols();
  Vector x(problem_size, 0.);   // unknown variables, initialized to zero
  Vector Atx(sample_size, 0.);  // maintained variables, initialized to zero
  
  // Step 3. Define your forward, or backward operators based on data and parameters
  double forward_operator_step_size = 0.0005;
 
  params.step_size = forward_operator_step_size;
 
  forward_grad_for_square_loss<Matrix> forward(&A, &b, &Atx, forward_operator_step_size);  
  using Forward = decltype(forward);
  // Step 4. Define your operator splitting scheme
  GradientDescentAlgorithm<forward_grad_for_square_loss<Matrix> > gd(&x, forward);

  // Step 6. Call the AROCK function
  double start_time = get_wall_time();  
  AROCK(gd, params);
  double end_time = get_wall_time();  

  print_parameters(params);
  
  cout << "Computing time is: " << end_time - start_time << endl;  
  // Step 7. Print results

  cout << "Objective value is: " << objective(Atx, b) << endl;
  cout << "---------------------------------" << endl;  
  return 0;
}

double objective(Vector& Atx, Vector& b) {

  int len = b.size();
  Vector temp = Atx;
  add(temp, b, -1.);
  double norm_of_temp = norm(temp);
  return 0.5 * norm_of_temp * norm_of_temp;

}
