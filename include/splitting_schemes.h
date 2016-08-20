#ifndef TMAC_INCLUDE_SPLITTING_SCHEMES_H
#define TMAC_INCLUDE_SPLITTING_SCHEMES_H

#include "operators.h"
#include "parameters.h"
#include "result.h"
#include "algebra.h"
#include <typeinfo>
#include "stdlib.h"
#include <iostream>
/****************************************************************
 * Header file for defining the splitting schemes. Each splitting
 * scheme is a functor. The functor contains the pointers to the
 * relevant data, related parameters, maintained variables and
 * unknown variable. It defines the following two functions.
 *
 * void operator(int index) {
 *     // update the unknown variable x at index
 *   // update the maintained variables
 * }
 *
 * void update_params(Params* params) {
 *   // update the operator related parameters
 *   // update the relaxation parameter
 * }
 *
 * The constructure should also be defined
 *
 *   OperatorSplitting(argument list) {
 *     // initialize the member variables with the input arguments
 *   }
 ***************************************************************/

extern "C" {
	// LU decomoposition of a general matrix
	void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);

	// generate inverse of a matrix given its LU decomposition
	void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);

	// matrix-matrix multiplication product C= alphaA.B + betaC                                               
	void dgemm_(char* TRANSA, char* TRANSB, const int* M,
		const int* N, const int* K, double* alpha, double* A,
		const int* LDA, double* B, const int* LDB, double* beta,
		double* C, const int* LDC);

	// matrix-vector multiplication product Y= alphaA.X + betaY                                               
	void dgemv_(char* TRANS, const int* M, const int* N,
		double* alpha, double* A, const int* LDA, double* X,
		const int* INCX, double* beta, double* C, const int* INCY);

	// calculate l2 norm
	double dnrm2_(int* N, double* X, int* INCX);

	// vector addition
	double daxpy_(int* N, double *alpha, double* x, int* INCX, double* y, int* INCY);

	// vector scale
	double dscal_(int* N, double* alpha, double* x, int* INCX);

	// vector dot
	double ddot_(int* N, double* x, int* INCX, double* y, int* INCY);

}

class SchemeInterface {
public:
	//update internal scheme parameters
	virtual void update_params(Params* params) = 0;
	//compute and apply coordinate update, return S_{index}
	virtual double operator() (int index) = 0;
	//compute and store S_{index} in variable S_i
	virtual void operator() (int index, double& S_i) = 0;
	//apply block of S stored in s to solution vector
	virtual void update(Vector& s, int range_start, int num_cords) = 0;
	//apply coordinate of S stored in s to solution vector
	virtual void update(double s, int idx) = 0;
	//update rank worth of cache_vars based on num_threads
	virtual void update_cache_vars(int rank, int num_threads) = 0;
};


// PPA:
template <typename Backward>
class ProximalPointAlgorithm : public SchemeInterface {
public:
	Backward prox;
	Vector* x;
	double relaxation_step_size;
	ProximalPointAlgorithm(Vector* xx, Backward p, double s) {
		x = xx;
		prox = p;
		prox.update_step_size(s);
	}

	void update_params(Params* params) {
		prox.update_step_size(params->get_step_size());
		relaxation_step_size = params->get_tmac_step_size();
	}

	double operator() (int index) {
		// Step 1: read the old x[index]
		double old_x_at_idx = (*x)[index];
		// Step 2: local calculation
		double S_i = old_x_at_idx - prox(x, index);
		// Step 3: get the most recent x[index] before updating it
		old_x_at_idx = (*x)[index];
		// Step 4: update x[index]
		(*x)[index] -= relaxation_step_size * S_i;
		// Step 5: update the maintained variable
		prox.update_cache_vars(old_x_at_idx, (*x)[index], index);
		return S_i;
	}

	// TODO: implement this
	void operator() (int index, double& S_i) {
	}

	void update(Vector& s, int range_start, int num_cords) {
		for (size_t i = 0; i < num_cords; ++i) {
			(*x)[i + range_start] -= relaxation_step_size * s[i];
		}
	}

	void update(double s, int idx) {
		(*x)[idx] -= relaxation_step_size * s;
	}

	void update_cache_vars(int rank, int num_threads) {
		prox.update_cache_vars(x, rank, num_threads);
	}

};


// gradient descent algorithm
template <typename Forward>
class GradientDescentAlgorithm : public SchemeInterface {
public:
	Forward forward;
	Vector *x;
	double relaxation_step_size;

	GradientDescentAlgorithm(Vector* x_, Forward forward_) {
		x = x_;
		forward = forward_;
		relaxation_step_size = 1.;
	}

	void update_params(Params* params) {
		forward.update_step_size(params->get_step_size());
		relaxation_step_size = params->get_tmac_step_size();
	}

	// update x[index] and update the maintained variables
	// x^{k+1} = x^k + eta_k (\hat x^k - T \hat x^k)
	double operator() (int index) {

		// Step 1: read the old x[index]
		double old_x_at_idx = (*x)[index];

		// Step 2: local calculation
		double forward_grad_at_idx = forward(x, index);
		double S_i = old_x_at_idx - forward_grad_at_idx;

		// Step 3: get the most recent x[index]
		old_x_at_idx = (*x)[index];

		// Step 4: update x at index
		(*x)[index] -= relaxation_step_size * S_i;

		// Step 5: update the maintained variable Atx
		double diff = (*x)[index] - old_x_at_idx;
		forward.update_cache_vars(old_x_at_idx, (*x)[index], index);
		return S_i;
	}

	void operator() (int index, double& S_i) {
		// Step 1: read the old x[index]
		double old_x_at_idx = (*x)[index];
		// Step 2: local calculation
		double forward_grad_at_idx = forward(x, index);
		S_i = old_x_at_idx - forward_grad_at_idx;
	}

	void update(Vector& s, int range_start, int num_cords) {
		for (size_t i = 0; i < num_cords; ++i) {
			(*x)[i + range_start] -= relaxation_step_size * s[i];
		}
	}

	void update(double s, int idx) {
		(*x)[idx] -= relaxation_step_size * s;
	}

	//update rank worth of cache_vars based on num_threads
	void update_cache_vars(int rank, int num_threads) {
		forward.update_cache_vars(x, rank, num_threads);
	}

};


// forward backward splitting scheme
template <typename Forward, typename Backward>
class ForwardBackwardSplitting : public SchemeInterface {
public:
	Forward forward;
	Backward backward;
	Vector* x;
	double relaxation_step_size;

	ForwardBackwardSplitting(Vector* x_, Forward forward_, Backward backward_) {
		x = x_;
		forward = forward_;
		backward = backward_;
		relaxation_step_size = 1.;
	}

	void update_params(Params* params) {
		// TODO: forward and backward might use different step sizes
		forward.update_step_size(params->get_step_size());
		backward.update_step_size(params->get_step_size());
		relaxation_step_size = params->get_tmac_step_size();
	}

	double operator() (int index) {
		// Step 1: read the old x[index]
		double old_x_at_idx = (*x)[index];
		// Step 2: local calculation
		double forward_grad_at_idx = forward(x, index);
		double val = backward(forward_grad_at_idx, index);
		double S_i = old_x_at_idx - val;
		// Step 3: get the most recent x[index] before updating it
		old_x_at_idx = (*x)[index];
		// Step 4: update x at index 
		(*x)[index] -= relaxation_step_size * S_i;
		// Step 5: update the maintained variable Atx
		forward.update_cache_vars(old_x_at_idx, (*x)[index], index);
		return S_i;
	}

	void operator() (int index, double &S_i) {
		// Step 1: read the old x[index]
		double old_x_at_idx = (*x)[index];
		// Step 2: local calculation
		double forward_grad_at_idx = forward(x, index);
		double val = backward(forward_grad_at_idx);
		S_i = old_x_at_idx - val;
	}

	void update(Vector& s, int range_start, int num_cords) {
		for (size_t i = 0; i < num_cords; ++i) {
			(*x)[i + range_start] -= relaxation_step_size * s[i];
		}
	}

	void update(double s, int idx) {
		(*x)[idx] -= relaxation_step_size * s;
	}

	//update rank worth of cache_vars based on num_threads
	void update_cache_vars(int rank, int num_threads) {
		forward.update_cache_vars(x, rank, num_threads);
	}

};


// backward forward splitting scheme
// WARNING: this won't work, if the forward step has maintained variables.
template <typename Backward, typename Forward>
class BackwardForwardSplitting : public SchemeInterface {
public:
	Forward forward;
	Backward backward;
	Vector *x;
	Vector y; // each new operator will have a copy of this guy
	double relaxation_step_size;

	BackwardForwardSplitting(Vector* x_, Backward backward_, Forward forward_) {
		x = x_;
		forward = forward_;
		backward = backward_;
		relaxation_step_size = 1.;
		y.resize(x->size());
	}

	void update_params(Params* params) {
		// TODO: forward and backward might use different step sizes
		forward.update_step_size(params->get_step_size());
		backward.update_step_size(params->get_step_size());
		relaxation_step_size = params->get_tmac_step_size();
	}

	double operator() (int index) {
		// Step 1: get the old x[index]
		double old_x_at_idx = (*x)[index];
		// Step 2: local computation
		// first apply the backward operator
		backward(x, &y);
		// then apply the forward operator on y
		double forward_grad_at_idx = forward(&y, index);
		double S_i = old_x_at_idx - forward_grad_at_idx;
		// Step 3: get the most recent x[index]
		old_x_at_idx = (*x)[index];
		// Step 4: update x[index]
		(*x)[index] -= relaxation_step_size * S_i;
		// Step 5: update the maintained variables
		// TODO: update cached_variable of the backward operator
		return S_i;
	}

	// TODO: implement this for sync-operator
	void operator()(int index, double &S_i) {
	}

	void update(Vector& s, int range_start, int num_cords) {
		for (size_t i = 0; i < num_cords; ++i) {
			(*x)[i + range_start] -= relaxation_step_size * s[i];
		}
	}

	void update(double s, int idx) {
		(*x)[idx] -= relaxation_step_size * s;
	}

	// TODO: for sync-parallel
	void update_cache_vars(int rank, int index) {
		std::cerr << "functionality not provided." << std::endl;
		exit(EXIT_FAILURE);
	}

};


// Peaceman-Rachford Splitting
// x^{k+1} = (1 - eta_k) x^k + eta_k (2 * Second - I )(2 * First - I) (x^k)
// which can be simplified to the following
// y^k = First(x^k)
// z^k = Second(2 y^k - x^k)
// x^{k+1} = x^k + 2 eta_k (z^k - y^k)
template <typename First, typename Second>
class PeacemanRachfordSplitting : public SchemeInterface {
public:
	First  op1;
	Second op2;
	Vector *x;
	Vector y; // each new operator will have a copy of this guy
	Vector z;
	Vector z_old;
	double relaxation_step_size;
	Result res;
	Matrix A;
	Vector b;
	Vector c;
	double alpha;
	double gamma;
	Matrix inv_AAt;
	Matrix E;
	Vector dual;
	Vector dual_old;

	PeacemanRachfordSplitting(Vector* x_, Matrix A_, Vector b_, Vector c_, Matrix inv_AAt_, Matrix E_, First op1_, Second op2_, double relaxation_step_size_, double alpha_, double gamma_) {
		x = x_;
		A = A_;
		b = b_;
		c = c_;
		inv_AAt = inv_AAt_;
		E = E_;
		op1 = op1_;
		op2 = op2_;
		//relaxation_step_size = 1.;
		relaxation_step_size = relaxation_step_size_;
		y.resize(x->size());
		z.resize(x->size());
		dual.assign(A.rows(), 0.);
		alpha = alpha_;
		gamma = gamma_;
	}

	void update_params(Params* params) {
		// TODO: forward and backward might use different step sizes
		op1.update_step_size(params->get_step_size());
		op2.update_step_size(params->get_step_size());
		relaxation_step_size = params->get_tmac_step_size();
	}

	void stop_check(Params* params, Result* res) {
        // define auxiliary variables for lapack use
		dual_old = dual;
		double* dual_ = &dual[0];
		double* dual_old_ = &dual_old[0];
		double* A_ = &A(0, 0);
		double* b_ = &b[0];
		double* c_ = &c[0];
		double* y_ = &y[0];
		double* z_ = &z[0];
		double* inv_AAt_ = &inv_AAt[0];
		double* E_ = &E[0];

		int A_row = A.rows();
		int problem_size = A.cols();
		double done = 1.;
		double dzero = 0.;
		double mone = -1.;
		int ione = 1;
		char* no = "N";
		char* tr = "T";

		// define auxiliary variables for intermediate storage
		Vector Ay(A_row, 0.);
		Vector Ac(A_row, 0.);
		Vector Az(A_row, 0.);
		Vector temp(A_row, 0.);
		Vector Atdual(problem_size, 0.);
		Vector sd(problem_size, 0.);
		Vector sd_proj(problem_size, 0.);
		Vector Adz(A_row, 0.);
		Vector Atddual(problem_size, 0.);
		Vector Atddual_proj(problem_size, 0.);
		double* Ay_ = &Ay[0];
		double* Ac_ = &Ac[0];
		double* Az_ = &Az[0];
		double* temp_ = &temp[0];
		double* Atdual_ = &Atdual[0];
		double* sd_ = &sd[0];
		double* sd_proj_ = &sd_proj[0];
		double* Adz_ = &Adz[0];
		double* Atddual_ = &Atddual[0];
		double* Atddual_proj_ = &Atddual[0];

		// compute Ay, Ac, Az
		dgemv_(tr, &problem_size, &A_row, &done, A_, &problem_size, y_, &ione, &dzero, Ay_, &ione);
		dgemv_(tr, &problem_size, &A_row, &done, A_, &problem_size, c_, &ione, &dzero, Ac_, &ione);
		dgemv_(tr, &problem_size, &A_row, &done, A_, &problem_size, z_, &ione, &dzero, Az_, &ione);

		//pfeas = norm(b-Ay)/(1 + norm(b))
		double  pfeas = 0;
		dscal_(&A_row, &mone, Ay_, &ione);
		daxpy_(&A_row, &done, b_, &ione, Ay_, &ione);
		pfeas = dnrm2_(&A_row, Ay_, &ione) / (1 + dnrm2_(&A_row, b_, &ione));

		// calculate dual
		for (int i = 0; i < A_row; i++) {
			temp[i] = alpha * Ac[i] + (b[i] - Az[i]) / gamma;
		}
		dgemv_(no, &A_row, &A_row, &done, inv_AAt_, &A_row, temp_, &ione, &dzero, dual_, &ione);

		// reset dual_old = dual - dual_old
		dscal_(&A_row, &mone, dual_old_, &ione);
		daxpy_(&A_row, &done, dual_, &ione, dual_old_, &ione);

		double pcost = dot(c, y);
		double dcost = dot(b, dual);
		double absgap = abs(pcost - dcost);
		double relgap = absgap / (1 + abs(pcost) + abs(dcost));
		double gap1 = absgap / abs(pcost);
		double gap2 = absgap / abs(dcost);

		// calculate Atdual
		dgemv_(no, &problem_size, &A_row, &done, A_, &problem_size, dual_, &ione, &dzero, Atdual_, &ione);

		// dfeas = norm(c-Atdual-proj(c-Atdual)) / (1 + norm(c))
		for (int i = 0; i < problem_size; i++) {
			sd[i] = c[i] - Atdual[i];
		}
		op1(&sd, &sd_proj);
		daxpy_(&problem_size, &mone, sd_proj_, &ione, sd_, &ione);
		double dfeas = dnrm2_(&problem_size, sd_, &ione) / (1 + dnrm2_(&problem_size, c_, &ione));

		// criteria for Optimal
		if (pfeas < params->feastol && dfeas < params->feastol && relgap < params->feastol) {
			res->status = 1; //"SOLVED"
			res->optimal = pcost;
			res->gap = relgap;
		}

		// reset z_old = (z - z_old)/max(1., norm(z-z_old))
		double* z_old_ = &z_old[0];
		dscal_(&problem_size, &mone, z_old_, &ione);
		daxpy_(&problem_size, &done, z_, &ione, z_old_, &ione);
		double var = 1 / max(1., dnrm2_(&problem_size, z_old_, &ione));
		dscal_(&problem_size, &var, z_old_, &ione);

		double unbndc = ddot_(&problem_size, c_, &ione, z_old_, &ione);

		// calculate Adz
		dgemv_(tr, &problem_size, &A_row, &done, A_, &problem_size, z_old_, &ione, &dzero, Adz_, &ione);
		double unbndf = dnrm2_(&problem_size, Adz_, &ione);

		// criteria for Unboundedness
		if (unbndc < 0 && abs(unbndf) < params->unbdd && abs(pfeas) < params->unbdd && gap1 > (1 - params->epsilon) && gap1 < (1 + params->epsilon))
			res->status = 'Unbounded';


		var = 1 / max(1., dnrm2_(&A_row, dual_old_, &ione));
		dscal_(&A_row, &var, dual_old_, &ione);

		// calculate -Atddual
		dgemv_(no, &problem_size, &A_row, &mone, A_, &problem_size, dual_old_, &ione, &dzero, Atddual_, &ione);
		op1(&Atddual, &Atddual_proj);

		// reset Atddual = Atddual - proj(Atddual)
		daxpy_(&problem_size, &mone, Atddual_proj_, &ione, Atddual_, &ione);

		double infesq = norm(Atddual, 2);
		double infesb = dot(b, dual_old);

		// criteria for Primal Infeasible
		if (infesq < params->infea && infesb >0 && gap2 > (1 - params->epsilon) && gap2 < (1 + params->epsilon))
			res->status = 'Primal Infeasible';

		return;

	}

	double operator() (int index) {
		// Step 1: get the old x[index]
		double old_x_at_idx = (*x)[index];
		// Step 2: local computation
		op1(x, &y);
		// z = 2y - x
		z_old = z;
		z = y;
		scale(z, 2.);
		add(z, *x, -1.);
		// prox(z, i), it doesn't have to be type I
		double temp = op2(&z, index);
		// Step 3: get the most recent x[index]
		old_x_at_idx = (*x)[index];
		// Step 4: update x at index 
		(*x)[index] += 2 * relaxation_step_size * (temp - y[index]);
		// Step 5: update the maintained variables
		return temp - y[index];
	}

	// TODO: implement this for sync-operator
	void operator()(int index, double &S_i) {
	}

	void update(Vector& s, int range_start, int num_cords) {
		for (size_t i = 0; i < num_cords; ++i) {
			(*x)[i + range_start] -= relaxation_step_size * s[i];
		}
	}

	void update(double s, int idx) {
		(*x)[idx] -= relaxation_step_size * s;
	}

	void update_cache_vars(int rank, int index) {
	}

};


#endif
