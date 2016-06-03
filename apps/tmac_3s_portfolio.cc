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

double objective(Vector& c, Matrix& Q);

int main(int argc, char *argv[]) {
  Params params;
  // Step 1. Parse the input argument
  string data_file_name;
  string label_file_name;
  double lambda = 0.02;
  
  parse_input_argv_mm(&params, argc, argv, data_file_name, label_file_name, lambda);

  // Step 2. Load the data or generate synthetic data, define matained variables
  Matrix Q;
  Vector epsilon;
  
  loadMarket(Q, data_file_name);
  loadMarket(epsilon, label_file_name);
  
  // set para
  double operator_step_size = 0.0018;
  params.tmac_step_size = 0.5;
  params.step_size = operator_step_size;
  
  int problem_size = Q.cols();
  params.problem_size = problem_size;
  
  // define auxilary variables
  Vector x(problem_size, 0.);   // unknown variables, initialized to zero

  portfolio_3s<Matrix> portfolio(&Q, &epsilon, lambda, operator_step_size);
  GradientDescentAlgorithm<portfolio_3s<Matrix> > three_s(&x, portfolio);
  

  double start_time = get_wall_time();  
  TMAC(three_s, params);  
  double end_time = get_wall_time();
  print_parameters(params);  
  cout << "Computing time is: " << end_time - start_time << endl;

  Vector y(problem_size, 0.);
  portfolio.project_D2(&x, &y);
  cout << "Objective is: " << objective(y, Q) << endl;
  
}

double objective(Vector& x, Matrix& Q) {
  
  double res = 0.;
  int n = x.size();
  Vector Qx(n, 0.);
  multiply(Q, x, Qx);
  return 0.5 * dot(x, Qx);

}
