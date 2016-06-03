/***************************************************
* Example: solving SVM of squared hinge loss with FBS
*
* min 1/2 <x,x> + lambda sum_i max(0, (1 - b_i a_i'x))
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

double objective(SpMat& A, Vector& x);

void preprocess_A(SpMat& A, Vector& b);

int main(int argc, char *argv[]) {

  // Step 0: Define the parameters and input file names  
  Params params;
  string data_file_name;
  double lambda;
  
  // Step 1. Parse the input argument
  parse_input_argv_libsvm(&params, argc, argv, data_file_name, lambda);
  
  // Step 2. Load the data or generate synthetic data, and preprocess the data
  SpMat A;
  Vector b;
  loadLibSVM(A, b, data_file_name);
  preprocess_A(A, b);

  // set params
  int num_features = A.rows();
  int num_samples = A.cols();
  // define auxilary variables
  SpMat At(num_samples, num_features);
  transpose(A, At);
  int problem_size = num_samples;
  Vector x(problem_size, 0.);   // unknown variables, initialized to zero
  Vector Ax(num_features, 0.);  // maintained variables, initialized to zero

  
  // Step 3. Define your forward, or backward operators based on data and parameters
  double operator_step_size = 0.01;
  params.step_size = operator_step_size;
  params.problem_size = num_samples;
  params.tmac_step_size = 1.;
  
  // forward operator
  forward_grad_for_dual_svm<SpMat> forward(&At, &Ax, operator_step_size);

  // backward operator
  Vector lower(num_samples, 0.), upper(num_samples, lambda);
  proj_box backward(&lower, &upper);

  Vector col_nrm(problem_size, 0.);
  calculate_column_norm(A, col_nrm);
  
  // print(col_nrm);
  double max_nrm = (*max_element(col_nrm.begin(), col_nrm.end()));
  cout << 1./max_nrm/max_nrm << endl;
  
  // Step 4. Define your operator splitting scheme
  ForwardBackwardSplitting<forward_grad_for_dual_svm<SpMat>, proj_box> fbs(&x, forward, backward);

  // Step 5. Call the TMAC function  
  double start_time = get_wall_time();
  TMAC(fbs, params);
  double end_time = get_wall_time();

  // Step 7. Print results
  print_parameters(params);
  cout << "Objective value is: " << objective(A, x) << endl;
  cout << "Computing time  is: " << end_time - start_time << endl;
  cout << "---------------------------------" << endl;
  cout << "# of nonzero in x: " << norm(x, 0) << endl;
  cout << "||x||_2 =: " << norm(x, 2) << endl;
  cout << "---------------------------------" << endl;
  return 0;
  
}

double objective(SpMat& A, Vector& x) {

  int rows = A.rows();
  Vector Ax(rows, 0.);
  multiply(A, x, Ax);
  double obj = 0.5 * dot(Ax, Ax) - sum(x);
  return obj;
}


void preprocess_A(SpMat& A, Vector& b) {

  int rows = A.rows(), cols = A.cols();
  for (int k = 0; k < rows; ++k) {
    for (SpMat::InnerIterator it(A,k); it; ++it) {
      A.coeffRef(k, it.index()) = it.value() * b(it.index());
    }
  }
  
}
