Operators
==========
We provide a total of 17 coordinate friendly operators. We categorize them into three classes: proximal operator, projection operator, and forward operator. The following is the general structure for the operators. You can use it as a template to add more operators.

.. code-block:: c++

   class OperatorInterface {
   public:
     virtual double operator() (Vector* v, int index = 0) {}
     virtual double operator() (double val, int index = 0) {}
     virtual void operator() (Vector* v_in, Vector* v_out) {}
     virtual void update_step_size(double step_size_) = 0;
     virtual void update_cache_vars(double old_x_i, double new_x_i, int index) = 0;
     virtual void update_cache_vars(Vector* x, int rank, int num_threads) = 0;
   };



In general, each operator struct has three overloaded parenthesis operator: one takes a scalar; one takes a scalar and an index then applys the coordinate update; one takes an input vector and a output vector, then perform full update to the input vector and saves the results in the output vector. Default contructor and customized contructors are also provided. If an operator involves data, the point to the data should be added as an member variable.

.. note::

   Some operator requires intermediate variable in order to perform cheap coordinate update. Most of the examples can be found in the plemented forward operators. 


Proximal operator
-----------------
The following table summarizes the implemented proximal operators. 

.. image:: ../proximal.png
    :width: 500px
    :align: center


Projection operator
-------------------
The following table summarizes the implemented projection operators. 

.. image:: ../projection.png
    :width: 500px
    :align: center


Forward operator
-------------------
The following table summarizes the implemented forward operators. They are essentially operators perform a (coordinate )gradient descent step for some smooth functions.

.. image:: ../forward.png
    :width: 500px
    :align: center

