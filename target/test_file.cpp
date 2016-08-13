#include<iostream>
#include"Operators.h"
#include"Vector.h"
#include"Matrix.h"
#include"mmio.h"
#include"splitting_schemes.h"
#include"tmac.h"

int main()
{

Vector x;

forward_grad_for_qp<SpMat> forward;

forward.step_size = 1;
forward.weight = 1;
forward.Q = new SpMat("jgl009.mtx");
forward.c = new Vector("jgl009.mtx");
int problem_size = forward.Q->number_of_rows;

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

TMAC(params,banana_split);

return 0;
}
