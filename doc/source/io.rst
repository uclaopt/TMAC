Input/output
***************
We provide several functions for displaying data and  handle file I/O. 

Print
======
We provide print functions for :cpp:type:`Matrix`, :cpp:type:`SpMat` and :cpp:type:`Vector`. The function declaration is in algebra.h.

.. cpp:function:: void print(Vector& x)

   Print a dense vector to cout.

.. cpp:function:: void print(Matrix& A)

   Print a dense matrix to cout.		  
		  
.. cpp:function:: void print(SpMat& A)
		  
   Print a sparse matrix to cout.		  


Read
====
We support a variety of data formats. If you data format is not supported, you can either convert your data files to one of the following supported format, or you can add the corresponding data file loader to TMAC. 


Matrix market format
--------------------
We support the matrix market format for dense matrix, sparse matrix, and dense vector. You can find the details of matrix market format `here <http://math.nist.gov/MatrixMarket/formats.html>`_.


.. cpp:function:: void loadMarket(Matrix& A, const std::string& file_name)

   Load data to a dense matrix.

.. cpp:function:: void loadMarket(Vector& x, const std::string& file_name)

   Load data to a dense vector.
		  
.. cpp:function:: void loadMarket(SpMat& A, const std::string& file_name)

   Load data to a sparse matrix.

		  
LIBSVM format
-------------
We will load the data in LIBSVM format to a sparse matrix for data samples, and a dense vector for labels. 

.. cpp:function:: bool loadLibSVM(Spmat& mat, Vector& v, const std::string& filename)

		  
Matlab format
-------------
Coming soon...



Write
=======
Coming soon...
