/*******************************************************************************
* Peaceman-Rachford splitting with second-order cone programming
*           min  c'x
*    subject to  Ax = b
*                x in Q(SOCs)
*    
*    Transform to:
*    min  alpha * c'x + indicate(Ax = b) + (1-alpha) * c'x + indicate(x in SOCs)
*    where 0 < alpha < 1
*
*    First Backward:  (1-alpha) * c'x + indicate(x in SOCs)
*    Second Backward: alpha * c'x + indicate(Ax = b)
*******************************************************************************/

#include <iostream>
#include "matrices.h"
#include "algebra.h"
#include "operators.h"
#include "parameters.h"
#include "result.h"
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
	Result res;
	std::string A_file_name;
	std::string b_file_name;
	std::string c_file_name;
	Matrix A;
	Vector b;
	Vector c;
	int cone_num;
	Vector cone_dim;
	double alpha;
	double gamma;
	double lambda;

	// Step 1. Parse the input argument
	parse_input_argv_socp(&params,
		argc,
		argv,
		A_file_name,
		b_file_name,
		c_file_name,
		cone_num,
		cone_dim,
		alpha,                      // weight for splitting
		gamma,                      // weight for proximal
		lambda);                     // weight for relaxed PRS

	// Step 2. Load the data or generate synthetic data, define maintained variables
    int problem_size = params.problem_size;
	params.tmac_step_size = 0.9;
	loadMarket(A, A_file_name);
	loadMarket(b, b_file_name);
	loadMarket(c, c_file_name);

	// define auxilary variables for lapack use
	int A_row = A.rows();
	Vector x(problem_size, 0.);    // unknown variables, initialized to zero
	Vector y(problem_size, 0.);    // unknown variables, initialized to zero
	Vector f(problem_size, 0.);    // f = A'(AA')^(-1)b
	Vector g(problem_size, 0.);    // g = A'(AA')^(-1)Ac
	Matrix D(problem_size, problem_size, 0.);  // D = A'(AA')^(-1)A
	Matrix inv_AAt(A_row, A_row, 0.);          
	Matrix E(problem_size, A_row, 0.);  
	double* A_ = &A(0, 0);
	double* inv_AAt_ = &inv_AAt(0, 0);
	double* D_ = &D(0, 0);
	double* E_ = &E(0, 0);
	double* b_ = &b[0];
	double* c_ = &c[0];
	double* f_ = &f[0];
	double* g_ = &g[0];
	double* x_ = &x[0];
	double* y_ = &y[0];
	char *no = "N";
	char *tr = "T";
	double done = 1.;              // one in double precision
	double dzero = 0.;             // zero in double precision
	int ione = 1;                  // one in integer precision
	
    // inv_AAt = (AA')^(-1)
	dgemm_(tr, no, &A_row, &A_row, &problem_size, &done, A_, &problem_size, A_ , &problem_size, &dzero, inv_AAt_, &A_row);               
	inverse(inv_AAt_, A_row);    

	// E = A'(AA')^(-1)
	dgemm_(no, tr, &problem_size, &A_row, &A_row, &done, A_, &problem_size, inv_AAt_, &A_row, &dzero, E_, &problem_size);      
	
	// D = A'(AA')^(-1)A
	dgemm_(no, tr, &problem_size, &problem_size, &A_row, &done, E_, &problem_size, A_, &problem_size, &dzero, D_, &problem_size);
	
	// f = A'(AA')^(-1)b
	dgemv_(no, &problem_size, &A_row, &done, E_, &problem_size, b_, &ione, &dzero, f_, &ione);
	
	// g = A'(AA')^(-1)Ac
	dgemv_(no, &problem_size, &problem_size, &done, D_, &problem_size, c_, &ione, &dzero, g_, &ione);

	// Step 3. Define the first projection operator

	// first operator
	proj_multi_soc_with_linear first(cone_num, cone_dim, &c, -(1 - alpha)*gamma);
	using First_Backward = decltype(first);
	
	// backward operator
	prox_linearpro_hyperplane second(&c, &f, &g, &D, gamma, alpha);
	using Second_Backward = decltype(second);
	
	// Step 4. Define your operator splitting scheme
	PeacemanRachfordSplitting<First_Backward, Second_Backward> prs(&x, A, b, c, inv_AAt, E, first, second, lambda/2, alpha, gamma);

	// Step 5. Call the TMAC function  
	double start_time = get_wall_time();
	TMAC(prs, params, res);
	double end_time = get_wall_time();
	
	first(&x, &y); // do one more backward step to recover the solution
	
	// Step 6. Print results
	print_parameters(params);
	cout << "Computing time is: " << end_time - start_time << endl;
	cout << "---------------------------------" << endl;
	cout << "||x||_2 =: " << norm(y, 2) << endl;
	cout << "||x||_inf =: " << norm(y, 3) << endl;
	cout << "optimal = " << dot(y, c) << endl;
	cout << "---------------------------------" << endl;
	return 0;
}
