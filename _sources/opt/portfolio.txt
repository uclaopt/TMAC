Portfolio optimization
======================
Assume that we have one unit of capital and :math:`n` assets to invest on. The :math:`i` th asset has an expected return rate :math:`\xi_i\ge 0`. Our goal is to find a portfolio with the minimal risk such that the expected return is no less than :math:`\lambda`. This problem can be formulated as

.. math::
   :nowrap:
      
   \begin{equation}
   \begin{array}{l}
   \min_x ~ \frac{1}{2} x^\top Q x, \\
   \text{subject to}~ x\ge0, \sum_{i=1}^m x_i\le 1,\, \sum_{i=1}^m\xi_i x_i\ge \lambda,
   \end{array}
   \end{equation}

where :math:`Q` is an :math:`n` by :math:`n` risk matrix. We solve the problem with TMAC by using three operator splitting method. 


Data preparation
-----------------
We consume :math:`Q` in matrix market format. You can easily modify the code to deal with other types of matrices. 


Usage
---------
In the bin folder, the executable file :cpp:type:`motac_3s_portfolio` minimizes the Huber loss::


  The usage for ./bin/motac_3s_portfolio is: 
  --------------------------------------------------------------
  ./bin/motac_3s_portfolio [options] 
  -data      <data file with matrix market format, size: n x n> 
  -label     <label file with matrix market format> 
  -nthread   <total number of threads, default is set to 2> 
  -epoch     <total number of epoch, default is set to 10> 
  --------------------------------------------------------------

  
Example
-----------

You can run the following command in the test directory to minimize the Huber loss for the large dense dataset::

  ./bin/motac_3s_portfolio -data ./data/port_Q.mtx -label ./data/port_eps.mtx -epoch 100 -nthread 2 -lambda 0.02

.. note::

   The datasets are not included in the source code, but you can download them from `here <https://www.dropbox.com/sh/neqh6ege48hut2x/AACv02EH19XN-N7DXADV2NrIa?dl=0>`_.
  
You can expect to get output similar to the following::

  Parameter settings:
  ---------------------------------
  Problem size:               1000
  TMAC step size:            0.5
  Operator step size:         0.0018
  Use controller:             false
  ---------------------------------
  Computing time is: 0.0860848
  Objective is: 6.09765e-05
