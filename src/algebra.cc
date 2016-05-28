/******************************************************************************
 * implementations for linear algebra operations, including functions such as
 * norm, dot, sub, add, multiply and so on.
 * Please refer to include/algebra.h for more details of the functions.
 *
 * Date Created:  01/29/2015
 * Date Modified: 01/29/2015
 *                02/19/2015
 *                04/27/2015 (removed some redundant functions and cleaned the code)
 *                07/19/2015 (modified the style of the code)
 * Authors:       Zhimin Peng, Yangyang Xu, Ming Yan, Wotao Yin
 * Contact:       zhimin.peng@math.ucla.edu
 ******************************************************************************/

#include <vector>
#include <string>
#include <algorithm>
#include <iomanip>
#include <stdexcept>   // std::invalid_argument
#include "matrices.h"
#include "constants.h"
#include "blas_sparse.h"

// Shrinkage function
double shrink(double x, double t) {
    if (x > t) {
      return x - t;
    } else if (x < -t) {
      return x + t;
    } else {
      return 0.;
    }
}

// Calculates the column norm of a matrix
void calculate_column_norm(Matrix& A, Vector& nrm) {
    int num_cols = A.cols();
    int num_rows = A.rows();
    std::fill(nrm.begin(), nrm.end(), 0.);
    int i, j;
    // #pragma omp parallel private(i, j)
    {
      // #pragma omp for schedule (static)
        for (i = 0; i < num_rows; ++i) {
            for (j = 0; j < num_cols; ++j) {
                nrm[j] += A(i,j) * A(i,j);
            }
        }
    }
    for (int j = 0; j < num_cols; ++j) {
        nrm[j] = sqrt(nrm[j]);
    }
    return;
}


void transpose(Matrix& A, Matrix& At) {
  int rows = A.rows();
  int cols = A.cols();
  for (int i = 0; i < rows; i++) {
    for(int j = 0; j < cols; j++) {
      At(j, i) = A(i, j);
    }
  }
}

void transpose(SpMat& A, SpMat& At) {
  At = A.transpose();
}


// Calculates the column norm of a matrix
void calculate_column_norm(SpMat& A, Vector& nrm) {
  std::fill(nrm.begin(), nrm.end(), 0.);  
    for (int k = 0; k < A.outerSize(); ++k) {
        for (SpMat::InnerIterator it(A,k); it; ++it) {
            nrm[it.index()] += it.value() * it.value();
        }
    }
    int num_cols = A.cols();
    for (int j = 0; j < num_cols; ++j) {
        nrm[j] = sqrt(nrm[j]);
    }
    return;
}

// Prints a vector
void print(Vector& x) {
    for (unsigned i = 0; i < x.size(); ++i) {
        std::cout << x[i] << " ";
    }
    std::cout << endl;
}


// Prints a dense matrix
void print(Matrix &A) {
    for (unsigned i = 0; i < A.rows(); ++i) {
        for (unsigned j = 0; j < A.cols(); ++j) {
            std::cout << setw(10) << A(i, j);
        }
        std::cout << endl;
    }
    std::cout << endl;
}


// Prints a sparse matrix
void print(SpMat &A) {
    std::cout << A;
    return;
}


namespace MyAlgebra {
// Norm function
double norm(Vector& x, int type) {
    double result = 0.0;
    if (type == 0) {
        for (unsigned i = 0; i < x.size(); ++i) {
            if (x[i] != 0) {
                result = result + 1;
            }
        }
    }
    else if (type == 2) {
        for (unsigned i = 0; i < x.size(); ++i) {
            result += x[i] * x[i];
        }
        result = sqrt(result);
    }
    else if (type == 1) {
        for (unsigned i = 0; i < x.size(); ++i) {
            result += fabs(x[i]);
        }
    }
    else if (type == 3) {
        for (unsigned i = 0; i < x.size(); ++i)
            result = max(fabs(x[i]), result);
    }
    else {
        throw std::invalid_argument("Unknown norm type!");
    }
    return result;
}


// Calculates the two norm of a vector
double norm(Vector& x) {
    return norm(x, 2);
}



// Implements sum(a)
double sum(Vector& a) {
  double res = 0.;
  for (unsigned i = 0; i < a.size(); ++i) {
    res += a[i];
  }
  return res;
}



// Implements a = a + scalar * A(row, :)
void add(Vector& a, Matrix& A, int row, double scalar) {
    for (unsigned i = 0; i < a.size(); ++i) {
        a[i] += scalar * A(row, i);
    }
}


// Implements a = a - scalar * A(row, :)
void add(Vector& a, SpMat& A, int row, double scalar) {
    double result = 0.;
    for (SpMat::InnerIterator it(A, row); it; ++it) {
        a[it.index()] -= scalar * it.value();
    }
}



// Implements a = a + scalar * A(row, :)
void add(Vector* a, Matrix* A, int row, double scalar) {
  for (unsigned i = 0; i < a->size(); ++i) {
    (*a)[i] += scalar * (*A)(row, i);
  }
}


// Implements a = a - scalar * A(row, :)
void add(Vector* a, SpMat* A, int row, double scalar) {
  double result = 0.;
  for (SpMat::InnerIterator it(*A, row); it; ++it) {
    (*a)[it.index()] += scalar * it.value();
  }
}


// Implements a = a - scalar * A(row, :)
void add(SpVec& a, SpMat& A, int row, double scalar) {
    double result = 0.;
    for (SpMat::InnerIterator it(A, row); it; ++it) {
        a.coeffRef(it.index()) += scalar * it.value();
    }
}


// Implements a = a - scalar * A(row, :)
void add(SpVec& a, Matrix& A, int row, double scalar) {
    double result = 0.;
    for (int i = 0; i < a.size(); ++i)
        a.coeffRef(i) += scalar * A(row, i);
}


// Implements a = a + lambda * b
void add(Vector &a, Vector& b, double lambda = 1.) {
    for (unsigned i = 0; i < a.size(); ++i)
        a[i] += lambda * b[i];
    return;
}


// Implements a = a + lambda * b
void add(Vector &a, double val = 1.) {
    for (unsigned i = 0; i < a.size(); ++i)
      a[i] += val;
    return;
}



void scale(Vector &a, double lambda) {
    for (unsigned i = 0; i < a.size(); ++i)
      a[i] *= lambda;
    return;
}

// a = a + lambda * b
void add(Vector &a, SpVec& b, double lambda = 1.) {
    for (SpVec::InnerIterator it(b); it; ++it) {
        a[it.index()] += it.value();
    }
    return;
}


// Calculates the inner product of two vectors
double dot(Vector &a, Vector &b) {
    double result = 0.;
    for (unsigned i = 0; i < a.size(); ++i) {
        result += a[i] * b[i];
    }
    return result;
}


// Calcuates inner product of A(row, :) * x
double dot(SpMat& A, Vector& x, int row) {
    double result = 0.;
    for (SpMat::InnerIterator it(A, row); it; ++it) {
        result += it.value() * x[it.index()];
    }
    return result;
}

// Calcuates inner product of A(row, :) * x
double dot(SpMat* A, Vector* x, int row) {
    double result = 0.;
    for (SpMat::InnerIterator it((*A), row); it; ++it) {
      result += it.value() * (*x)[it.index()];
    }
    return result;
}




// Calcuates inner product of A(row, :) * x
double dot(Matrix& A, Vector& x, int row) {
    double result = 0.;
    for (unsigned i = 0; i < A.cols(); ++i) {
        result += A(row, i) * x[i];
    }
    return result;
}


// Calcuates inner product of A(row, :) * x
double dot(Matrix* A, Vector* x, int row) {
    double result = 0.;
    for (unsigned i = 0; i < (A->cols()); ++i) {
      result += (*A)(row, i) * (*x)[i];
    }
    return result;
}




// Calculates A' * x
void trans_multiply(Matrix& A, Vector&x, Vector& Atx) {
  std::fill(Atx.begin(), Atx.end(), 0.);  
  int m = A.rows(), n = A.cols();
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      Atx[j] += A(i, j) * x[i];
    }
  }
}

// Calculates A' * x
void trans_multiply(SpMat& A, Vector&x, Vector& Atx) {
  std::fill(Atx.begin(), Atx.end(), 0.);
  int rows = A.outerSize();
  for (int k = 0; k < rows; ++k) {
    for (SpMat::InnerIterator it(A,k); it; ++it) {
      Atx[it.index()] += it.value() * x[k];
    }
  }
  return;
}


// Calculates Ax = A * x
void multiply(SpMat &A, Vector &x, Vector& Ax) {

    std::fill(Ax.begin(), Ax.end(), 0);
    for (int k = 0; k < A.outerSize(); ++k) {
        for (SpMat::InnerIterator it(A,k); it; ++it) {
            Ax[k] += it.value() * x[it.index()];
        }
    }
    return;
}


// Calculates Ax = A * x
void multiply(Matrix &A, Vector &x, Vector& Ax) {
  int m = A.rows();
  int n = A.cols();
  std::fill(Ax.begin(), Ax.end(), 0);
  for (int i = 0; i < m; ++i) {
    for(int j = 0; j < n; ++j) {
      Ax[i] += A(i, j) * x[j];
    }
  }
  return;
}


// Caculates AAt = A * A' for sparse matrix
void multiply(SpMat &A, SpMat &AAt) {
    AAt = (A * A.transpose()).pruned();
}


// Calculates A * A' for dense matrix
void multiply(Matrix &A, Matrix &AAt) {
  int i, j, k;
  int m = A.rows(), n = A.cols();
  // #pragma omp parallel private(i, j, k)
  {
    // #pragma omp for schedule (static)
    for (i = 0; i < m; ++i) {
      for (j = 0; j < m; ++j) {
        for(k = 0; k < n; ++k) {
          AAt(i, j) += A(i, k) * A(j, k);
        }
      }
    }
  }
  return;
}

}

// Selects B = A(start:end, :)
void copy(SpMat& A, SpMat& B, int start, int end) {
    int n = A.cols();
    int j = 0;
    B = A.block(start, 0, end - start, n);
}


// Selects B = A(start:end, :)
void copy(Matrix& A, Matrix& B, int start, int end) {
    int n = A.cols();
    B.resize(end - start, n);
    int i = 0, j = 0;
    for(i = start; i < end; i++) {
        for(j = 0; j < n; j++) {
            B(i - start,j) = A(i,j);
        }
    }
}


// Selects y = x(start:end)
void copy(Vector& x, Vector& y, int start, int end) {
    for (int i = start; i < end; i++) {
        y[i - start] = x[i];
    }
}



/***************************
 *  BLASAlgebra library
 ***************************/



namespace BLASAlgebra {
// Norm function

extern "C"{
  double dnrm2_(const int *N, const double* v, const int *incv);
  double dasum_(const int *N, const double* v, const int *incv);
  double daxpy_(const int *N, const double* alpha, const double* x, const int *incx, const double* y, const int* incy);
  void dscal_(const int *N, const double* alpha, const double* v, const int* incx);
  double ddot_( const int *N, const double *a, const int *inca, const double *b, const int *incb );
  void dgemv_(char* TRANS, const int* M, const int* N,
              double* alpha, double* A, const int* LDA, double* X,
              const int* INCX, double* beta, double* C, const int* INCY);
}


double norm(Vector& x, int type) {
    double result = 0.0;
    int N = x.size();
    int one = 1;
    if (type == 0) {
        for (unsigned i = 0; i < x.size(); ++i) {
            if (x[i] != 0) {
                result = result + 1;
            }
        }
    }
    else if (type == 2) {
      result = dnrm2_(&N, x.data(), &one);
    }
    else if (type == 1) {
      result = dasum_(&N, x.data(), &one);
    }
    else if (type == 3) {
        for (unsigned i = 0; i < x.size(); ++i)
            result = max(fabs(x[i]), result);
    }
    else {
        throw std::invalid_argument("Unknown norm type!");
    }
    return result;
}


// Calculates the two norm of a vector
double norm(Vector& x) {
    return norm(x, 2);
}


// Implements sum(a)
double sum(Vector& a) {
  double res = 0.;
  for (unsigned i = 0; i < a.size(); ++i) {
    res += a[i];
  }
  return res;
}


void add(Vector& a, Matrix& A, int row, double scalar) {
  int rows = A.rows();
  int cols = A.cols();
  double* mtx_ptr = (A.data() + cols * row);
  int one = 1;
  // a-> y, b -> x
  daxpy_(&cols, &scalar, mtx_ptr, &one, a.data(), &one);

}

void add(Vector& a, SpMat&  A, int row, double scalar) {

  int* out_ptr = A.outerIndexPtr();
  int* inn_ptr = A.innerIndexPtr();
  double* val_ptr = A.valuePtr();
  int nnz = *(out_ptr + row + 1) - *(out_ptr + row);
  int sp_vec_start_idx = *(out_ptr + row);


  val_ptr += sp_vec_start_idx;
  inn_ptr += sp_vec_start_idx;

  BLAS_dusaxpy(nnz, scalar, val_ptr, inn_ptr, a.data(), 1, blas_zero_base);

}

void add(SpVec&  a, SpMat&  A, int row, double scalar) {

  MyAlgebra::add(a, A, row, scalar);
  
}

void add(SpVec&  a, Matrix& A, int row, double scalar) {
  
  MyAlgebra::add(a, A, row, scalar);

}

void add(Vector* a, Matrix* A, int row, double scalar) {
  int rows = A->rows();
  int cols = A->cols();
  double* mtx_ptr = (A->data() + cols * row);
  int one = 1;
  // a-> y, b -> x
  daxpy_(&cols, &scalar, mtx_ptr, &one, a->data(), &one);
}

void add(Vector* a, SpMat*  A, int row, double scalar) {

  int* out_ptr = A->outerIndexPtr();
  int* inn_ptr = A->innerIndexPtr();
  double* val_ptr = A->valuePtr();
  int nnz = *(out_ptr + row + 1) - *(out_ptr + row);
  int sp_vec_start_idx = *(out_ptr + row);


  val_ptr += sp_vec_start_idx;
  inn_ptr += sp_vec_start_idx;

  BLAS_dusaxpy(nnz, scalar, val_ptr, inn_ptr, a->data(), 1, blas_zero_base);

}


void add(Vector &a, Vector& b, double lambda = 1.) {

  int N = a.size();
  int one = 1;
  // a-> y, b -> x
  daxpy_(&N, &lambda, b.data(), &one, a.data(), &one);

}


void add(Vector &a, double lambda = 1.) {

  MyAlgebra::add(a, 1);

}



void add(Vector &a, SpVec&  b, double lambda = 1.) {

  int* inn_ptr = b.innerIndexPtr();
  double* val_ptr = b.valuePtr();
  int nnz = b.nonZeros();
  BLAS_dusaxpy(nnz, lambda, val_ptr, inn_ptr, a.data(), 1, blas_zero_base);

}

void scale(Vector &a, double lambda) {

  int N = a.size();
  int one = 1;
  dscal_(&N, &lambda, a.data(), &one);
  
}

double dot(SpMat& A,  Vector& x, int row) {

  int* out_ptr = A.outerIndexPtr();
  int* inn_ptr = A.innerIndexPtr();
  double* val_ptr = A.valuePtr();
  
  int nnz = *(out_ptr + row + 1) - *(out_ptr + row);
  int sp_vec_start_idx = *(out_ptr + row);

  double res;
  BLAS_dusdot(blas_no_conj, nnz, val_ptr + sp_vec_start_idx,
              inn_ptr + sp_vec_start_idx, x.data(), 1, &res, blas_zero_base);  
  
  return res;
  
}

double dot(Vector& a, Vector& b ){
  int N = a.size();
  int one = 1;
  return ddot_( &N, a.data(), &one, b.data(), &one);
}


double dot(Matrix& A, Vector& x, int row) {
  int N = x.size();
  int one = 1;
  int num_cols = A.cols();
  
  return ddot_( &N, A.data() + row * num_cols, &one, x.data(), &one);  

}


double dot(SpMat* A,  Vector* x, int row) {
  int* out_ptr = A->outerIndexPtr();
  int* inn_ptr = A->innerIndexPtr();
  double* val_ptr = A->valuePtr();
  
  int nnz = *(out_ptr + row + 1) - *(out_ptr + row);
  int sp_vec_start_idx = *(out_ptr + row);

  double res;
  BLAS_dusdot(blas_no_conj, nnz, val_ptr + sp_vec_start_idx,
              inn_ptr + sp_vec_start_idx, x->data(), 1, &res, blas_zero_base);  
  return res;
}


double dot(Matrix* A, Vector* x, int row) {

  int N = x->size();
  int one = 1;
  int num_cols = A->cols();
  return ddot_( &N, A->data() + row * num_cols, &one, x->data(), &one);  

}

void trans_multiply(Matrix& A, Vector&x, Vector& Atx) {
  double alpha = 1., beta = 0.;
  int M = A.rows(), N = A.cols();
  int one = 1;
  char c = 'N';
  int LDA = N;
  std::fill(Atx.begin(), Atx.end(), 0.);  
  // BLAS is in column major
  dgemv_(&c, &N, &M, &alpha, A.data(), &LDA, x.data(), &one, &beta, Atx.data(), &one);
}


// Using sparse BLAS requires create matrix handle.
void multiply(SpMat &A,  Vector &x, Vector& Ax) {
  MyAlgebra::multiply(A, x, Ax);
}


// Using sparse BLAS requires create matrix handle.
void trans_multiply(SpMat &A,  Vector &x, Vector& Atx) {
  MyAlgebra::trans_multiply(A, x, Atx);
}


void multiply(Matrix &A, Vector &x, Vector& Ax) {
  double alpha = 1., beta = 0.;
  int M = A.rows(), N = A.cols();
  int one = 1;
  char c = 'T';
  int LDA = N;
  // BLAS is in column major
  dgemv_(&c, &N, &M, &alpha, A.data(), &LDA, x.data(), &one, &beta, Ax.data(), &one);
}

}
