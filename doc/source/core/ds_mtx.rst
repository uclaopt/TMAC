Dense matrix
=============
The :cpp:type:`Matrix` class for dense matices. Its purpose is to provide convenient mechanisms for performing basic matrix
operations, such as constructing the matrix, retrieving the size, setting the value at  :math:`(i, j)`. The matrix is stored in row row major and The precision of the entries are double.


An example of generating a :math:`m x n` matrix of real double-precision numbers with
value 3.14 is the following

  .. code-block:: cpp

     Matrix x(m, n, 3.14);

     
The underlying data storage for :cpp:type:`Matrix` is :cpp:type:`std::vector<double>`.

.. cpp:class:: Matrix

   .. rubric:: Constructors and destructors

   .. cpp:function:: Matrix( )

      This creates a default empty matrix.


   .. cpp:function:: Matrix(int m, int n)

      A matrix of size :math:`m x n` is created. It allocates memory, and initializes the value to 0.

   .. cpp:function:: Matrix(int m, int n, double val)

      A matrix of size :math:`m x n` is created. It allocates memory, and initializes the value to
      :cpp:type:`val`.


   .. rubric:: Member functions
      
   .. cpp:function:: size_t rows() const

      Return the number of rows of the matrix.

   .. cpp:function:: size_t cols() const

      Return the number of cols of the matrix.

   .. cpp:function:: size_t nnz() const

      Return the number of nonzeros in the matrix.

   .. cpp:function:: void resize(size_t m_, size_t n_)

      Resize the matrix to size :math:`m_ x n_`.

   .. cpp:function:: double& operator()(size_t i, size_t j)

      Retrieve the :math:`(i, j)` th entry of the matrix.

   .. cpp:function:: const double& operator()(size_t i, size_t j) const

      Retrieve the :math:`(i, j)` th entry of the matrix.
      
   
   .. rubric:: Member functions for file I/O

   .. cpp:function:: void read(std::istream &in)

      Read  data from input stream and save the data as a Matrix.

   .. cpp:function:: void write(std::ostream &out)

      Write a Matrix to a out stream with matrix market format.
