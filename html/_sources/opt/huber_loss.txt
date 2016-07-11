Huber loss
======================
We apply TMAC to minimize the following Huber loss

.. math::
   \min_x f_\delta(Ax-b),

where :math:`A` is an :math:`m` by :math:`n` matrix and :math:`f_\delta(y)=\sum_i\text{Huber}_\delta(y_i)` with

.. math::
   :nowrap:

   \begin{equation}
   \text{Huber}_\delta(v)=\left\{\begin{array}{ll}
   \frac{1}{2}v^2,&|v|\leq\delta,\\
   \delta(|v|-\frac{1}{2}\delta),&\text{otherwise}.
   \end{array}\right.
   \end{equation}
   
We solve the problem with TMAC by using gradient descent.


Data preparation
-----------------
We minimize the Huber loss with dense matrix in matrix market format, and sparse matrix in LIBSVM format. You can easily modify the code to deal with other types of matrices. 



Usage
---------
In the bin folder, the executable file :cpp:type:`motac_gd_huber` minimizes the Huber loss::

  The usage for motac_gd_huber is::

    ./motac_gd_huber [options] 
               -data       < matrix market file for A >
               -label      < matrix market file for b > 
               -nthread    < total number of threads, default: 1. > 
               -epoch      < total number of epochs, default: 10. > 
  
Example
-----------

You can run the following command in the test directory to minimize the Huber loss for the large dense dataset::

  ./bin/motac_gd_huber -data ./data/ds_large_A.mtx -label ./data/ds_large_b.mtx -epoch 100 -nthread 2 -lambda 1.

.. note::

   The datasets are not included in the source code, but you can download them from `here <https://www.dropbox.com/sh/neqh6ege48hut2x/AACv02EH19XN-N7DXADV2NrIa?dl=0>`_.
  
You can expect to get output similar to the following::

  Parameter settings:
  ---------------------------------
  Problem size:               1000
  TMAC step size:            0.5
  Operator step size:         0.5
  Use controller:             false
  ---------------------------------
  Computing time is: 0.7631
  Objective value is: 39.8258
  ---------------------------------
