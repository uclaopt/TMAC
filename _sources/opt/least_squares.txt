Least squares
======================
We apply TMAC to solve the following least squares problem

.. math::
   \min_x \frac{1}{2}\|Ax-b\|^2,

where :math:`A` is an :math:`m` by :math:`n` matrix. We solve the problem with TMAC by using gradient descent.


Data preparation
-----------------
We solve the least squares problem with dense matrix in matrix market format, and sparse matrix in LIBSVM format. You can easily modify the code to deal with other types of matrices. 



Usage
---------
In the bin folder, the executable file :cpp:type:`motac_gd_ls` solves the least squares problem:

  The usage for motac_gd_ls is::

    ./motac_gd_ls [options] 
               -data       < matrix market file for A >
               -label      < matrix market file for b > 
               -nthread    < total number of threads, default: 1. > 
               -epoch      < total number of epochs, default: 10. > 
  
Example
-----------

You can run the following command in the test directory to solve the least squares problem for the large dense dataset::

  ./bin/motac_gd_ls -data ./data/ds_large_A.mtx -label ./data/ds_large_b.mtx -epoch 100 -nthread 2 -lambda 1.

.. note::

   The datasets are not included in the source code, but you can download them from `here <https://www.dropbox.com/sh/neqh6ege48hut2x/AACv02EH19XN-N7DXADV2NrIa?dl=0>`_.
  
You can expect to get output similar to the following::

  Parameter settings:
  ---------------------------------
  Problem size:               1000
  TMAC step size:            0.5
  Operator step size:         0.8
  Use controller:             false
  ---------------------------------
  Computing time is: 0.105244
  Objective value is: 28.74
  ---------------------------------
