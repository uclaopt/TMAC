/**********************************
 * Example : solving linear equations
 * 
 *    Ax = b
 * where A is a diagonally dominant matrix, stored
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
#include "MarketIO.h"
#include <thread>
using namespace std;
//using namespace MyAlgebra;
#include "algebra_namespace_switcher.h"


double cal_residue(Matrix A, Vector x, Vector b);

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
  params.tmac_step_size = 1.;
  Vector x(problem_size, 0.);   // unknown variables, initialized to zero

  
  // Step 3. Define your forward, or backward operators based on data and parameters
 
  linear_eqn_jacobi_operator<Matrix> J(&A, &b);  
  
  // Step 4. Define your operator splitting scheme
  GradientDescentAlgorithm<linear_eqn_jacobi_operator<Matrix> > Jacobi(&x, J);

  // Step 6. Call the TMAC function
  double start_time = get_wall_time();
  TMAC(Jacobi, params);
  double end_time = get_wall_time();    
  
   // Step 7. Print results
  print_parameters(params);
  cout << "Residue        is: " << cal_residue(A, x, b)  << endl;
  cout << "Computing time is: " << end_time - start_time << endl; 
  cout << "---------------------------------" << endl;  
  return 0;
}

double cal_residue(Matrix A, Vector x, Vector b){

  int m=A.rows();
  Vector Ax(m, 0.);
  multiply(A, x, Ax);
  add(b, Ax, -1.);
  return norm(b);
  
}

