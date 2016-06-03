Quadratically constrained quadratic program
=============================================
We consider the following special QCQP problem:

.. math::
   \min_x -\frac{1}{2} x^T Q x + b^T x + d, ~\text{s.t.}~ \|x\|_2 \leq 1

where $Q$ is a symmetric positive definite matrix. We solve the problem with TMAC by using backward-forward splitting.


Usage
---------
In the bin folder, the executable file :cpp:type:`motac_bfs_l2_ball_qp` solves the above problem.

  The usage for motac_prs_demo is::

    ./motac_bfs_l2_ball [options] 
               -nthread       < total number of threads, default: 1. > 
               -epoch         < total number of epochs, default: 10. > 
	       -problem_size  < the dimension of the problem, default: 0 >

Example
-----------

You can run the following command in the test directory to solve the l1 regularized least square problem for the large dense dataset::

  ./bin/motac_bfs_l2_ball_qp -nthread 1 -problem_size 100 -epoch 200

  
You can expect to get output similar to the following::

  Parameter settings:
  ---------------------------------
  Problem size:               100
  TMAC step size:            0.5
  Operator step size:         0.05
  Use controller:             false
  ---------------------------------
  Objective value is: 45.8257
  Computing time  is: 1.97591
  ||x||_2 = 0.834947
  ---------------------------------


.. note::

   The current app use synthesis data, you can easily modify the code to load your own data. You can also modify the feasible set to other supported constraints.
   
