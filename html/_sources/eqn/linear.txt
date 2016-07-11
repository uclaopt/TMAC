Linear equation
=================
We apply TMAC to solve the following linear equation

.. math::
   Ax = b,

where :math:`A` is an :math:`m\times m` matrix with entries :math:`a_{ij}` satisfying :math:`|a_{ii}|>\sum_{j\neq i}|a_{ij}|` for all :math:`i`. Such a matrix is called diagonally dominant.

Data preparation
-----------------
We implement the TMAC linear equation solver  with dense matrix in matrix market format, and sparse logistic regression with sparse matrix in LIBSVM format. You can easily modify the code to deal with other types of matrices. 



Usage
---------
In the bin folder, the executable file :cpp:type:`motac_jacobi` solves the linear equationm:

  The usage for motac_jacobi is::

    ./motac_jacobi [options] 
               -data       < matrix market file for A >
               -label      < matrix market file for b > 
               -nthread    < total number of threads, default: 1. > 
               -epoch      < total number of epochs, default: 10. > 
  
Example
-----------

You can run the following command in the test directory to solve the linear equation for the dense dataset::

  ./bin/motac_jacobi -data ./data/ds_A.mtx -label ./data/ds_b.mtx -epoch 10 -nthread 2.

.. note::

   The datasets are not included in the source code, but you can download them from `here <https://www.dropbox.com/sh/neqh6ege48hut2x/AACv02EH19XN-N7DXADV2NrIa?dl=0>`_.
  
You can expect to get output similar to the following::

  Parameter settings:
  ---------------------------------
  Problem size:               100
  TMAC step size:            1
  Operator step size:         0
  Use controller:             false
  ---------------------------------
  Residue        is: 5.52774
  Computing time is: 0.00127196
  ---------------------------------

