/***************************************************
 * Example 2: solving the following problem 
 * 
 *    min  0.5 * x' * Q * x + c' * x + d
 *    s.t. ||x||_2 <= 1
 *
 * with backward foward splitting
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
#include "algebra_namespace_switcher.h"

double objective(Matrix&, Vector&, Vector&, double);

int main(int argc, char *argv[]) {

  // Step 0: Define the parameters and input file names  
  Params params;
  // Step 1. Parse the input argument
  parse_input_argv_demo(&params, argc, argv);

  // Step 2. Load the data or generate synthetic data, define matained variables
  int n = params.problem_size;
  Matrix Q(n, n, 0.);   // matrix is a row major matrix
  for (int i = 0; i < n; i++) {
    Q(i, i) = 10.; 
  }
  for (int i =0; i < n - 1; i++) {
    Q(i, i+1) = 1.;
    Q(i+1, i) = 1.;
  }
  Vector c(n, -1.);
  double d = 0.5 * dot(c, c);
  // set para

  params.tmac_step_size = 0.5;
  int problem_size = params.problem_size;
  // define auxilary variables
  Vector x(problem_size, 0.);   // unknown variables, initialized to zero
  Vector y(problem_size, 0.);   // unknown variables, initialized to zero

  // Step 3. Define your forward, or backward operators based on data and parameters
  double operator_step_size = 0.5;
  params.step_size = operator_step_size;
  
  // forward operator
  forward_grad_for_qp<Matrix> forward(&Q, &c, operator_step_size);
  using Forward = decltype(forward);
  // backward operator
  double radius = 1.;
  proj_l2_ball backward(radius);
  using Backward = decltype(backward);

  // Step 4. Define your operator splitting scheme
  BackwardForwardSplitting<Backward,Forward> bfs(&x, backward, forward);  

  // Step 6. Call the TMAC function  
  double start_time = get_wall_time();
  TMAC(bfs, params);  
  double end_time = get_wall_time();

  backward(&x, &y); // do one more backward step to recover the solution
  
  // Step 7. Print results
  print_parameters(params);
  cout << "Objective value is: " << objective(Q, c, y, d) << endl;
  cout << "Computing time  is: " << end_time - start_time << endl;
  cout << "||x||_2 = " << norm(y, 2) << endl;
  cout << "---------------------------------" << endl;  

  return 0;
}

double objective(Matrix& Q, Vector& c, Vector& x, double d) {
  int n = c.size();
  Vector Qx(n, 0);
  multiply(Q, x, Qx);
  return 0.5 * dot(Qx, x) + dot(x, c) + d;
}
