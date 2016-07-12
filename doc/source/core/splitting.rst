Splitting schemes
===================
We provide the following splitting schemes:

   * Proximal point algorithm;
   * Gradient descent algorithm;
   * Forward backward splitting;
   * Backward forward splitting;
   * Relaxed Peaceman Rachford splitting;
   * 3 operator splitting (coming soon).

The following is a template for the splitting schemes. Each splitting scheme is a functor. The functor contains the pointers to the operators, and unknown variable.


.. code-block:: c++

   class SchemeInterface {
   public:
     //update internal scheme parameters
     virtual void update_params(Params* params) = 0;
     //compute and apply coordinate update, return S_{index}
     virtual double operator() (int index) =0;
     //compute and store S_{index} in variable S_i
     virtual void operator() (int index, double& S_i) =0;
     //apply block of S stored in s to solution vector
     virtual void update(Vector& s, int range_start, int num_cords) = 0;
     //apply coordinate of S stored in s to solution vector
     virtual void update(double s, int idx ) = 0;
     //update rank worth of cache_vars based on num_threads
     virtual void update_cache_vars(int rank, int num_threads) = 0;
   };



In general, each splitting scheme will be templated on one or more operators. Operator objective will be defined as the member variables.

.. note::

   Inside the parenthesis operator, relaxation will be used. Operator related cached variables will also be updated. 

