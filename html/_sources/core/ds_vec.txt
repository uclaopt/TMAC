Dense vector
=============
The :cpp:type:`Vector` class for dense vectors. It is a subclass of :cpp:type:`std::vector\<double>`.
Its purpose is to provide convenient mechanisms for performing basic vector
operations, such as constructing the vector, retrieving the size and setting the value
for a given index. The precision of the entries are double.

An example of generating a vector with size :math:`n` of real double-precision numbers with
value 3.14 is the following


  .. code-block:: cpp

     Vector x(n, 3.14);

     
The underlying data storage for :cpp:type:`Vector` is :cpp:type:`std::vector\<double>`.

.. cpp:class:: Vector

   .. rubric:: Constructors and destructors

   .. cpp:function:: Vector( )

      It creates a default empty vector.


   .. cpp:function:: Vector(int n)

      A vector of size :cpp:type:`n` is created. It allocates memory, and initializes the value
      to 0.

   .. cpp:function:: Vector(int n, double val)

      A vector of size :cpp:type:`n` is created. It allocates memory, and initializes the value to
      :cpp:type:`val`.

   .. rubric:: Member functions for file I/O

   .. cpp:function:: read(std::istream &in)

      Read data from input stream and save the data as a Vector.

   .. cpp:function:: write(std::ostream &out)

      Write a Vector to a out stream with matrix market format.

