Support vector machine (SVM)
=============================
We apply TMAC to solve the following support vector machine problem

.. math::
   \min_x \lambda \, \|x\|^2 + \sum_{i=1}^N \max(0, 1 - b_i \cdot a_i^T x)

where :math:`\{(a_1, b_1), ..., (a_N, b_N)\}` is the set of data-label pairs, and :math:`\lambda>0` is the regularization parameter.

We solve the dual formulation of the SVM with coordinate update method.


Data preparation
-----------------
We support LIBSVM format. The code can be easily adopt to data in  matrix market format.



Usage
---------
In the bin folder, the executable file :cpp:type:`motac_fbs_dual_svm` solves the dual formulation of SVM.

  The usage for motac_fbs_dual_svm is::

    ./motac_fbs_dual_svm [options] 
               -data       < data in LIBSVM format >
               -nthread    < total number of threads, default: 1. > 
               -epoch      < total number of epochs, default: 10. > 
               -lambda     < regularization parameter, default 1. > 

  
Example
-----------

You can run the following command to train a SVM model on the news20 dataset::

  ./bin/motac_fbs_dual_svm -data ./data/news20.binary -epoch 10 -nthread 2 -lambda .001
  
You can expect to get output similar to the following::

  Parameter settings:
  ---------------------------------
  Problem size:               19996
  TMAC step size:             1
  Operator step size:         0.01
  Use controller:             false
  ---------------------------------
  Objective value is: -19.6826
  Computing time  is: 0.894006
  ---------------------------------
  # of nonzero in x: 19996
  ||x||_2 =: 0.141407
  ---------------------------------

   
