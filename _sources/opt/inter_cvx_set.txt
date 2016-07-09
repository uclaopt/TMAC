Intersection of convex sets
============================
We apply relaxed PRS to find the intersection of two convex sets. The problem can be modeled as the following optimization problem

.. math::
   \min_x \iota_{\|x\|_2 \leq 1}(x) + \iota_{\|x\|_\infty \leq 0.1} (x)

where :math:`\iota_S (x)` is the indicator function of the set S. We apply the TMAC with the Peaceman-Rachford splitting operator. To ensure coordinate friendly, We first apply the projection operator to the :math:`\ell_2` ball, then we apply the projection to the :math:`\ell_{\infty}` ball.


Usage
---------
In the bin folder, the executable file :cpp:type:`motac_prs_demo` solves the above problem.

  The usage for motac_prs_demo is::

    ./motac_prs_demo [options] 
               -nthread       < total number of threads, default: 1. > 
               -epoch         < total number of epochs, default: 10. > 
	       -problem_size  < the dimension of the problem, default: 0 >

Example
-----------

You can run the following command in the test directory to solve the l1 regularized least square problem for the large dense dataset::

  ./bin/motac_prs_demo -epoch 100 -nthread 2 -problem_size 100

  
You can expect to get output similar to the following::

	       
  Parameter settings:
  ---------------------------------
  Problem size:               100
  TMAC step size:            0.9
  Operator step size:         0
  Use controller:             false
  ---------------------------------
  Computing time is: 0.00682497
  ---------------------------------
  ||x||_2 =: 0.590442
  ||x||_inf =: 0.0994609
  ---------------------------------
