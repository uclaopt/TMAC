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
		
    struct OperatorSplitting {

      double relaxation_step_size;
      
      void operator(int index) {
        // update the unknown variable x at index
        // update the maintained variables
      }
   
      void update_params(Params* params) {
        // update the operator related parameters
        // update the relaxation parameter    
      }
   
      // The constructure should also be defined
      OperatorSplitting(argument list) {
        // initialize the member variables with the input arguments
      }
   };

   In general, each splitting scheme will be templated on one or more operators. Operator objective will be defined as the member variables.
   
.. note::

   Inside the parenthesis operator, relaxation will be used. Operator related cached variables will also be updated. 

