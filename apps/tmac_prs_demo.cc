/**************************************************************
 * Demo Peaceman-Rachford splitting with the following problem
 * 
 *    min  indicate(||x||_2 <= 1) + indicate(||x||_inf <= 0.1),
 *
 * i.e., we want to find an intersection point of the two 
 * convex sets: ||x||_2 <= 1 and ||x||_inf <= 0.1.
 *
 *************************************************************/

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
#include "algebra_namespace_switcher.h"
using namespace std;


int main(int argc, char *argv[]) {

  // Step 0: Define the parameters and input file names  
  Params params;
  
  // Step 1. Parse the input argument
  parse_input_argv_demo(&params, argc, argv);

  // Step 2. Load the data or generate synthetic data, define matained variables
  int problem_size = params.problem_size;
  params.tmac_step_size = 0.9;

  // define auxilary variables
  Vector x(problem_size, 10.);   // unknown variables, initialized to zero
  Vector y(problem_size, 0.);    // unknown variables, initialized to zero

  // Step 3. Define the first projection operator

  // first operator
  Vector lb(problem_size, -0.1);
  Vector ub(problem_size, .1);
  proj_box first(&lb, &ub);
  using First_Backward = decltype(first);
  
  // backward operator
  double radius = 1.;
  proj_l2_ball second(radius);
  using Second_Backward = decltype(second);

  // Step 4. Define your operator splitting scheme
  PeacemanRachfordSplitting<First_Backward, Second_Backward> prs(&x, first, second);  

  // Step 5. Call the TMAC function  
  double start_time = get_wall_time();
  TMAC(prs, params);
  double end_time = get_wall_time();  

  first(&x, &y); // do one more backward step to recover the solution
  
  // Step 6. Print results
  print_parameters(params);
  cout << "Computing time is: " << end_time - start_time << endl;
  cout << "---------------------------------" << endl;
  cout << "||x||_2 =: " << norm(y, 2) << endl;
  cout << "||x||_inf =: " << norm(y, 3) << endl;
  cout << "---------------------------------" << endl;  
  return 0;
}
