Linear algebra
==================
The :cpp:type:`algebra.h` file provides a few functions for linear
algebra operations. Thoese functions include calculating the norm of a vector, computing inner product of vectors, performing matrix vector operations, etc. We provide both naive implementation (under namespace MyAlgebra) and (Sparse) BLAS based implementation (under namespace BLASAlgebra) for these functions. You can specify the namespace through the :cpp:type:`algebra_namespace_switcher.h` header file.

.. cpp:class:: algebra.h


   .. rubric:: Vector functions
      
   .. cpp:function:: double norm(Vector& x, int type)

      Return the norm of a vector of different type:

		     * type = 0: :math:`\ell_0` norm, the number of nonzeros;
		     * type = 1: :math:`\ell_1` norm. Sum of the absolute values;
		     * type = 2: :math:`\ell_2` norm. Euclidean norm;
		     * type = 3: :math:`\ell_\infty` norm. Maximum of the absolute value.
	
   .. cpp:function:: double norm(Vector& x)

      Calculate the :math:`\ell_2` norm of a vector.

   .. cpp:function:: void sub(Vector& a, Vector& b)

      Subtract b from a, i.e., a = a - b.

   .. cpp:function:: void add(Vector &a, Vector& b)

      Add dense vector b to vector a, i.e., a = a + b.
      
   .. cpp:function:: void add(Vector &a, SpVec&  b)

      Add sparse vector b to vector a, i.e., a = a + b.

   .. cpp:function:: double dot(Vector &a, Vector &b)

      Calculate the inner product of vector a and vector b.

   .. cpp:function:: void copy(Vector& x, Vector& y, int start, int end)

      Copy part of a vector to another vector, i.e., y = x(start:end)
		     
   .. rubric:: Matrix functions
	       
   .. cpp:function:: void calculate_column_norm(Matrix& A, Vector& nrm)

      Calculate the :math:`\ell_2` norm for each column of dense matrix A.
   
   .. cpp:function:: void calculate_column_norm(SpMat& A, Vector& nrm)		     

      Calculate the :math:`\ell_2` norm for each column of sparse matrix A.

   .. cpp:function:: void multiply(SpMat &A,  SpMat &AAt)
		     
      Multiply A with A', i.e., AAt = A * A'		     
		     
   .. cpp:function:: void multiply(Matrix &A, Matrix &AAt)
		     
      Multiply A with A', i.e., AAt = A * A'		     
      
   .. cpp:function:: void copy(Matrix& A, Matrix& B, int start, int end)

      Copy part of a dense matrix to another dense matrix, i.e., B = A(start:end, :)
   
   .. cpp:function:: void copy(SpMat& A, SpMat& B, int start, int end)

      Copy part of a sparse matrix to another sparse matrix, i.e., B = A(start:end, :)

      
   .. rubric:: Matrix vector functions

   .. cpp:function:: void trans_multiply(Matrix& A, Vector&x, Vector& Atx)

      Calculate :math:`Atx = A^T \, x`.
	       

   .. cpp:function:: void sub(Vector& a, Matrix& A, int row, double scalar)

       Subtract a scalar times a row of matrix from a given vector,
       i.e., a = a - scalar * A(row, :)
		     
   .. cpp:function:: void sub(Vector& a, SpMat&  A, int row, double scalar)

       Subtract a scalar times a row of matrix from a given vector,
       i.e., a = a - scalar * A(row, :)
		     
   .. cpp:function:: void sub(SpVec&  a, SpMat&  A, int row, double scalar)

       Subtract a scalar times a row of matrix from a given vector,
       i.e., a = a - scalar * A(row, :)
		     
   .. cpp:function:: void sub(SpVec&  a, Matrix& A, int row, double scalar)

       Subtract a scalar times a row of matrix from a given vector,
       i.e., a = a - scalar * A(row, :)
      
   .. cpp:function:: double dot(SpMat& A,  Vector& x, int row)
   .. cpp:function:: double dot(Matrix& A, Vector& x, int row)

      Calcuate inner product of A(row, :) * x.

   .. cpp:function:: void multiply(SpMat &A,  Vector &x, Vector& Ax)
		     
      multiply sparse matrix A with x, i.e., :math:`Ax = A \cdot x`
      
   .. cpp:function:: void multiply(Matrix &A, Vector &x, Vector& Ax)

       multiply dense matrix A with x, i.e., :math:`Ax = A \cdot x`

   .. rubric:: Print functions    
   .. cpp:function:: void print(Vector& x)

      Print the vector.

   .. cpp:function:: void print(Matrix& x)

      Print the dense matrix.

   .. cpp:function:: void print(SpMat& x)

      Print the sparse matrix.
      
