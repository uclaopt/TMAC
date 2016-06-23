/******************************************************
 * Example : solving Portfolio Optimization problem
 * 
 *    min  0.5 * x' * Q * x
 *    s.t. x>=0,
 *         \sum_i x_i<=1,
 *         \sum_i \xi_i*x_i>=c
 *
 *****************************************************/

#include <iostream>
#include "matrices.h"
#include "algebra.h"
#include "operators.h"
#include "parameters.h"
#include "splitting_schemes.h"
#include "arock.h"
#include "util.h"
#include "MarketIO.h"
#include <thread>
using namespace std;
#include "algebra_namespace_switcher.h"

double objective(Matrix& Q, Vector& x);

int main(int argc, char *argv[]) {

  // Step 0. Define the parameters and input file names  
  Params params;
  string data_file_name;
  string label_file_name;
    
  // Step 1. Parse the input argument

  double lambda = 1.;  // the variable 'c' in the problem
  parse_input_argv_mm(&params, argc, argv, data_file_name, label_file_name, lambda);

  // Step 2. Load the data or generate synthetic data, define matained variables

  Matrix Q;
  Vector xi;
  loadMarket(Q, data_file_name);
  loadMarket(xi, label_file_name);  
  
  int problem_size = Q.rows();
  params.problem_size = problem_size;

  params.arock_step_size = 0.5;

  // try to change worker_type
  params.worker_type = "random";
//  params.worker_type = "gs";

  Vector x(problem_size, 0.);   // unknown variables, initialized to zero


  // Step 3. Define your forward, or backward operators based on data and parameters

  double three_s_operator_step_size = 0.002; 
  params.step_size = three_s_operator_step_size;
 
  portfolio_3s<Matrix> portfolio(&Q, &xi, lambda, three_s_operator_step_size);  

  // Step 4. Define your operator splitting scheme

  GradientDescentAlgorithm<portfolio_3s<Matrix> > three_s(&x, portfolio);

  // Step 5. (Optional) Define your controller

  
  // Step 6. Call the AROCK function

  double start_time = get_wall_time();  
  AROCK(three_s, params);
  double end_time = get_wall_time();  

  print_parameters(params);  

  // Step 7. Print results   

  cout << "Computing time is: " << end_time - start_time << endl;  
  
  Vector y(problem_size, 0.);
  portfolio.projection_D2(&x, &y);
  cout << "Objective value is: " << objective(Q, y) << endl;
  cout << "---------------------------------" << endl;  

  return 0;
}

double objective(Matrix& Q, Vector& x){

  int len = x.size();
  Vector Qx(len, 0.);
  multiply(Q, x, Qx);
  return 0.5 * dot(x, Qx);
}
