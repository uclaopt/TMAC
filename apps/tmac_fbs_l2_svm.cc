/***************************************************
* Example: solving SVM of squared hinge loss with FBS
*
* min 1/2 <x,x> + lambda/2 sum_i max(0, (1 - b_i a_i'x))^2
* 
* where a_i is the ith column of Matrix A
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
  string data_file_name;
  string label_file_name;
  double lambda;

  
  // Step 1. Parse the input argument
  parse_input_argv_mm(&params, argc, argv, data_file_name, label_file_name, lambda);

  
  // Step 2. Load the data or generate synthetic data, define matained variables
  Matrix A;
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
  forward_grad_for_square_hinge_loss<Matrix> forward(&A, &b, &Atx, operator_step_size, lambda);
  // backward operator
  prox_sum_square backward(operator_step_size, 1);

	
  // Step 4. Define your operator splitting scheme
  ForwardBackwardSplitting<forward_grad_for_square_hinge_loss<Matrix>, prox_sum_square> fbs(&x, forward, backward);

	
  // Step 5. (Optional) Define your controller

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

  double sum = 0;
  int size = b.size();
  for (int i = 0; i < size; i++){
	double temp = 1 - b[i] * Atx[i];
	if (temp > 0)
		sum += temp * temp;
  }

  double obj = 0.5 * dot(x, x) + lambda / 2 * sum;
  return obj;

}
