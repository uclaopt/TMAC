/***********************************************************************
 *  Header file for linear algebra operations, including functions       
 *  such as norm, dot, add, multiplication                          
 *                                                                       
 ***********************************************************************/

#ifndef TMAC_INCLUDE_ALGEBRA_H
#define TMAC_INCLUDE_ALGEBRA_H
#include <string>
#include <iomanip>
#include "matrices.h"


/*************************************************
 * shrinkage function or soft threshold function,
 *
 * return x -t, if x > t
 * return -x +t, if -x < t
 * return 0 otherwise.
 *
 * Input:
 *     x:      the input value
 *     (double)
 *     t:      the shrinkage parameter
 *     (double)
 *     
 * Output:
 *     result: the value after shrinkage
 *     (double)
 *
 ************************************************/
double shrink(double x, double t);



/***************************************************
 * calculate the l2 norm for each column of
 * dense matrix A
 *
 * Input:
 *     A -- the target matrix
 *          (Matrix)
 * Output:
 *   nrm -- vector of size num_cols to store results
 *          (Vector)
 ***************************************************/
void calculate_column_norm(Matrix& A, Vector& nrm);


/***************************************************
 * calculate the l2 norm for each column of
 * sparse matrix A
 *
 * Input:
 *     A -- the target matrix
 *          (SpMat)
 * Output:
 *   nrm -- vector of size num_cols to store results
 *          (Vector)
 ***************************************************/
void calculate_column_norm(SpMat& A, Vector& nrm);



// print a vector
void print(Vector& x);

// print a dense matrix
void print(Matrix &A);

// print a sparse matrix 
void print(SpMat &A);


/***************************************************
 * multiply A with A', i.e., AAt = A * A'
 * the calculation is efficient for row major matrix
 *
 * Input:
 *      A -- a dense or sparse matrix
 *
 * Output:
 *    AAt -- the result matrix
 **************************************************/
void multiply(SpMat  &A,  SpMat &AAt);
void multiply(Matrix &A, Matrix &AAt);



/*****************************************************
 * copy matrix A from start row to end-1 row to
 * matrix B, i.e., B = A(start:end-1, :)
 *
 * Input:
 *      A -- a dense or sparse matrix
 *  start -- the index for the starting row
 *    end -- the past-the-end index
 *
 * Output:
 *      B -- the submatrix of A
 ****************************************************/
void copy(Matrix& A, Matrix& B, int start, int end);
void copy(SpMat&  A, SpMat&  B, int start, int end);


/*****************************************************
 * copy part of a vector to another vector 
 * , i.e., y = x(start:end)
 *
 * Input:
 *      x -- a dense vector
 *  start -- the index for the starting row
 *    end -- the past one index for the ending index
 * Output:
 *      y -- the subvector of x
 ****************************************************/
void copy(Vector& x, Vector& y, int start, int end);

void transpose(Matrix& A, Matrix& At);
void transpose(SpMat& A, SpMat& At);




namespace MyAlgebra {
/************************************************* 
 * return the norm of a vector for a given type
 * 
 * Input:
 *       x -- the input vector
 *            (Vector)
 *    type -- 0 means zero norm,
 *            1 means l1 norm,
 *            2 means l2 norm.
 *            3 means infinity norm
 *            (int)
 * Output:
 *  result -- a scalar
 *              (double)
 ********************************************** */
double norm(Vector& x, int type);


/**************************************
 * calculate the l2 norm of a vector
 *************************************/
double norm(Vector& x);

double sum(Vector& x);


/**************************************************
 * add a scalar times a row
 * of matrix from a given vector,
 * i.e., 
 *   a = a - scalar * A(row, :)
 *
 * Input:
 *      a -- a vector
 *      A -- a matrix
 *    row -- the row number
 * scalar -- the scalar
 *
 * Output:
 *      a -- the added vetor
 *************************************************/
void add(Vector& a, Matrix& A, int row, double scalar); // dense matrix, dense vector
void add(Vector& a, SpMat&  A, int row, double scalar); // sparse matrix, dense vector
void add(SpVec&  a, SpMat&  A, int row, double scalar); // sparse matrix, sparse vector
void add(SpVec&  a, Matrix& A, int row, double scalar); // dense matrix, sparse vector


void add(Vector* a, Matrix* A, int row, double scalar); // dense matrix, dense vector
void add(Vector* a, SpMat*  A, int row, double scalar); // sparse matrix, dense vector


/***************************************************
 * add b to a, i.e.,
 * a = a + b
 *
 * Input:
 *      a -- a vector
 *      b -- a vector 
 *
 * Output:
 *      a -- the added vector
 **************************************************/
void add(Vector &a, Vector& b, double lambda = 1.); // add two dense vectors
void add(Vector &a, SpVec&  b, double lambda = 1.); // add two sparse vectors


void add(Vector &a, double val = 1.); // add two dense vectors


void scale(Vector &a, double lambda);

/***************************************************
 * calculate the inner product of two vectors
 * a = a' * b
 *
 * Input:
 *      a -- a vector
 *           (Vector)
 *      b -- a vector
 *           (Vector)
 *
 * Output:
 * result -- a scalar represent the inner product 
 ***************************************************/
double dot(Vector &a, Vector &b);


/***************************************************
 * calculate inner product of A(row, :) * x
 *
 * result = A(row, :) * x
 *
 * Input:
 *       a -- a vector
 *       b -- a vector 
 *
 * Output:
 *  result -- a scalar represent the inner product 
 ***************************************************/
double dot(SpMat& A,  Vector& x, int row); // sparse matrix
double dot(Matrix& A, Vector& x, int row); // dense matrix
double dot(SpMat* A,  Vector* x, int row); // sparse matrix
double dot(Matrix* A, Vector* x, int row); // dense matrix




/*************************************************
 * calcuate A' * x (A transpose multiplied x)
 *  Atx = A'*x
 *
 * Input:
 *      A -- a dense matrix
 *      x -- a dense vector
 *
 * Output:
 *    Atx -- the result vector
 *************************************************/
void trans_multiply(Matrix& A, Vector&x, Vector& Atx);
void trans_multiply(SpMat& A, Vector&x, Vector& Atx);




/*************************************************
 * multiply A with x
 *  Ax = A*x
 *
 * Input:
 *      A -- a dense or sparse matrix
 *      x -- a dense vector
 *
 * Output:
 *     Ax -- the result vector
 *************************************************/
void multiply(SpMat &A,  Vector &x, Vector& Ax);
void multiply(Matrix &A, Vector &x, Vector& Ax);

}


namespace BLASAlgebra {
/************************************************* 
 * return the norm of a vector for a given type
 * 
 * Input:
 *       x -- the input vector
 *            (Vector)
 *    type -- 0 means zero norm,
 *            1 means l1 norm,
 *            2 means l2 norm.
 *            3 means infinity norm
 *            (int)
 * Output:
 *  result -- a scalar
 *              (double)
 ********************************************** */
double norm(Vector& x, int type);


/**************************************
 * calculate the l2 norm of a vector
 *************************************/
double norm(Vector& x);

double sum(Vector& x);



/**************************************************
 * add a scalar times a row
 * of matrix from a given vector,
 * i.e., 
 *   a = a - scalar * A(row, :)
 *
 * Input:
 *      a -- a vector
 *      A -- a matrix
 *    row -- the row number
 * scalar -- the scalar
 *
 * Output:
 *      a -- the added vetor
 *************************************************/
void add(Vector& a, Matrix& A, int row, double scalar); // dense matrix, dense vector
void add(Vector& a, SpMat&  A, int row, double scalar); // sparse matrix, dense vector
void add(SpVec&  a, SpMat&  A, int row, double scalar); // sparse matrix, sparse vector
void add(SpVec&  a, Matrix& A, int row, double scalar); // dense matrix, sparse vector
void add(Vector* a, Matrix* A, int row, double scalar); // dense matrix, dense vector
void add(Vector* a, SpMat*  A, int row, double scalar); // sparse matrix, dense vector


/***************************************************
 * add b to a, i.e.,
 * a = a + b
 *
 * Input:
 *      a -- a vector
 *      b -- a vector 
 *
 * Output:
 *      a -- the added vector
 **************************************************/
void add(Vector &a, Vector& b, double lambda = 1.); // add two dense vectors
void add(Vector &a, double lambda = 1.); // add two dense vectors
void add(Vector &a, SpVec&  b, double lambda = 1.); // add two sparse vectors

void scale(Vector &a, double lambda);

/***************************************************
 * calculate the inner product of two vectors
 * a = a' * b
 *
 * Input:
 *      a -- a vector
 *           (Vector)
 *      b -- a vector
 *           (Vector)
 *
 * Output:
 * result -- a scalar represent the inner product 
 ***************************************************/
double dot(Vector &a, Vector &b);


/***************************************************
 * calculate inner product of A(row, :) * x
 *
 * result = A(row, :) * x
 *
 * Input:
 *       a -- a vector
 *       b -- a vector 
 *
 * Output:
 *  result -- a scalar represent the inner product 
 ***************************************************/
double dot(SpMat& A,  Vector& x, int row); // sparse matrix
double dot(Matrix& A, Vector& x, int row); // dense matrix
double dot(SpMat* A,  Vector* x, int row); // sparse matrix
double dot(Matrix* A, Vector* x, int row); // dense matrix




/*************************************************
 * calcuate A' * x (A transpose multiplied x)
 *  Atx = A'*x
 *
 * Input:
 *      A -- a dense matrix
 *      x -- a dense vector
 *
 * Output:
 *    Atx -- the result vector
 *************************************************/
void trans_multiply(Matrix& A, Vector&x, Vector& Atx);
void trans_multiply(SpMat& A, Vector&x, Vector& Atx);

/*************************************************
 * multiply A with x
 *  Ax = A*x
 *
 * Input:
 *      A -- a dense or sparse matrix
 *      x -- a dense vector
 *
 * Output:
 *     Ax -- the result vector
 *************************************************/
void multiply(SpMat &A,  Vector &x, Vector& Ax);
void multiply(Matrix &A, Vector &x, Vector& Ax);

}
// end of the BLASAlgebra namespace


#endif
