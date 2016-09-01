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

int main()
{

Vector x;

forward_grad_for_qp<SpMat> forward;

forward.step_size = 1;
forward.weight = 1;
SpMat A; loadMarket(A,"../data/jgl009.mtx");
Vector b; loadMarket(b,"../data/jgl009.mtx");


forward.Q = &A;
forward.c = &b;
int problem_size = forward.Q->rows();

prox_l1 backward;

backward.step_size = 2;
backward.weight = 2;

ForwardBackwardSplitting<forward_grad_for_qp<SpMat>,prox_l1> banana_split(&x,forward,backward);

Params params;
params.problem_size = problem_size;
params.max_itrs = 20;
params.tmac_step_size = 0.5;
params.total_num_threads = 4;
params.use_controller = false;
params.worker_type = "cyclic";
params.async = true;
params.step_size = 1;
// we do not need to touch step_size_rule as of now
int num_workers = params.use_controller ? params.total_num_threads - 1 : params.total_num_threads;
params.block_size = params.problem_size/num_workers;

 TMAC(banana_split,params);

return 0;
}
