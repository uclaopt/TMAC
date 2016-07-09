Regularized regression
======================
We apply TMAC to solve the following regularized empirical risk minimization problem

.. math::
   \min_x \lambda \, r(x) + \sum_{i=1}^N \ell(a_i^T x, b_i),

where :math:`\{(a_1, b_1), ..., (a_N, b_N)\}` is the set of data-label pairs, and :math:`\lambda>0` is the regularization parameter. We call
:math:`r(x)` and :math:`\ell(a_i^T x, b_i)` as the regularization function and loss function respectively.

We have implemented TMAC for :math:`r(x) = \|x\|_1`, :math:`\ell_i(x) = (a_i^T x - b_i)^2` (correspond to LASSO), :math:`\ell_i(x) = log(1+exp(-b_i \cdot a_i^T x))` (correspond to sparse logistic regression).


Data preparation
-----------------
We implement LASSO with dense matrix in matrix market format, and sparse logistic regression with sparse matrix in LIBSVM format. You can easily modify the code to deal with other types of matrices. 



Usage
---------
In the bin folder, the executable file :cpp:type:`motac_fbs_lasso` solves the :math:`\ell_1` regularized least square problem:

  The usage for motac_fbs_lasso is::

    ./motac_fbs_lasso [options] 
               -data       < matrix market file for A >
               -label      < matrix market file for b > 
               -nthread    < total number of threads, default: 1. > 
               -epoch      < total number of epochs, default: 10. > 
               -lambda     < regularization parameter, default 1. > 

  
Example
-----------

You can run the following command in the test directory to solve the l1 regularized least square problem for the large dense dataset::

  ./bin/motac_fbs_lasso -data ./data/ds_large_A.mtx -label ./data/ds_large_b.mtx -epoch 10 -nthread 2 -lambda 1.

.. note::

   The datasets are not included in the source code, but you can download them from `here <https://www.dropbox.com/sh/neqh6ege48hut2x/AACv02EH19XN-N7DXADV2NrIa?dl=0>`_.
  
You can expect to get output similar to the following::

  Parameter settings:
  ---------------------------------
  Problem size:               1000
  TMAC step size:            0.5
  Operator step size:         0.9
  Use controller:             false
  ---------------------------------
  Objective value is: 473.595
  Computing time  is: 0.008394
  ---------------------------------
  # of nonzero in x: 332
  ||x||_2 =: 332
  ---------------------------------


.. note::

   You can find other applications (l2 norm least square with backward forward splitting, SVM with squared hinge loss, etc) in the binary folder. The source codes for the apps are `here <https://github.com/ZhiminPeng/motac-new/tree/master/apps>`_.

   
