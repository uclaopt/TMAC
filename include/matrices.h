/**********************************************************************************
 *  Header file for the following classes:
 *   Vector:       a dense vector class, inherit from STL vector, with support
 *                 for matrix market file io.
 *
 *   Matrix:       a dense matrix class;
 *
 *   SpMat:        sparse matrix class from Eigen, data are stored
 *                  in Row Major, with compressed row storage.
 *   SpVec:        sparse vector class from Eigen, data are stored
 *                  in Row Major, with compressed row storage.
 *
 * This is a modification of Evgenii's code.
 *  http://Evgenii.Rudnyi.Ru/
 *
 * Date created:  01/15/2015
 * Date Modified: 02/19/2015
 *                04/27/2015 (add Google c++ style comments)
 *
 * Author: Zhimin Peng, Yangyang Xu, Ming Yan, Wotao Yin
 * Contact: zhimin.peng@math.ucla.edu
 *********************************************************************************/

#ifndef TMAC_INCLUDE_MATRICES_H
#define TMAC_INCLUDE_MATRICES_H
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <../lib/Eigen/Sparse>
#include <../lib/Eigen/Dense>
#include "auxiliary.h"

// shorten the name for sparse matrix from Eigen
typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SpMat;



// shorten the name for sparse vector from Eigen
typedef Eigen::SparseVector<double, Eigen::RowMajor> SpVec;


/**********************************************************************************
 *
 * matrix error namespace to display error message when matrix read or write failed
 *
 **********************************************************************************/
namespace MatrixError {
  
  struct MatrixError : public BaseError
  {
    MatrixError(const std::string &str1, const std::string &str2) : BaseError(std::string("MatrixError::") + str1, str2) {}
  };
  
  struct ReadError : public MatrixError
  {
    ReadError(const std::string &str) : MatrixError("ReadError", str) {}
  };

  struct WriteError : public MatrixError
  {
    WriteError(const std::string &str) : MatrixError("WriteError", str) {}
  };

}

/****************************************
 * a dense linear algebra vector class,
 * it based on std::vector<double>, we 
 * added some member functions to handle
 * file read and write.
 ****************************************/
class Vector : public std::vector<double>
{
public:
  // default constructor
  Vector(){}
  
  // construct a vector with size n_ and initialized with value val
  Vector(size_t n_, double val) : vector<double>(n_, val) {}

  // construct a vector with size n_ without initialization
  Vector(size_t n_) : vector<double>(n_) {}
  
  // destructor
  ~Vector()
  {
    vector<double>::clear();    
  }

  // overload () operator
  double& operator()(size_t idx)
  {
    return operator[](idx);
  }

  // read data from input stream to a vector
  void read(std::istream &in); // read data from matrix market format file

  // helper function for read data 
  void readData(std::istream &in, size_t n, bool IsSymmetric);

  // write data from matrix market format file
  void write(std::ostream &out); 

};


/*******************************************
 * Dense matrix class: row major. Also based
 * std::vector<double> class. added some
 * contruction functions for easy to use.
 * It also provides member functions for
 * file IO. 
 ******************************************/
class Matrix : public std::vector<double>
{
  size_t m; // number of rows
  size_t n; // number of columns
  
public: 

  // default constructor
  Matrix() : m(0), n(0) {}

  // construct a matrix with size m * n;
  Matrix(size_t m_, size_t n_) : vector<double>(m_*n_), m(m_), n(n_) {}

  // construct a matrix with size m * n and initial value val
  Matrix(size_t m_, size_t n_, double val) : vector<double>(m_*n_, val), m(m_), n(n_) {}

    // destructor
  ~Matrix()
  {
    vector<double>::clear();    
  }

  
  // resize the matrix
  void resize(size_t m_, size_t n_)
  {
    m = m_;
    n = n_;
    vector<double>::resize(m_*n_);
  }

  // reserve the memory for the vector
  // overload the reserve function from STL vector.
  void reserve(size_t m_, size_t n_)
  {
    vector<double>::reserve(m_*n_);
  }

  // clear the memory
  void clear()
  {
    m = n = 0;
    vector<double>::clear();
  }

  // return the number of rows
  size_t rows() const
  {
    return m;
  }

  // return the number columns
  size_t cols() const
  {
    return n;
  }

  // return the number of nonzeros
  size_t nnz()
  {
    return size() - count(begin(), end(), 0.);
  }

  // overload () operator, returns the (i,j)th
  // entry of the matrix
  double& operator()(size_t i, size_t j)
  {
    return operator[](i*n + j);
  }
  
  // overload () operator, returns the (i,j)th
  // entry of the matrix
  const double& operator()(size_t i, size_t j) const
  {
    return operator[](i*n + j);
    //return *((*this).begin() + i*n + j);
  }
  
  // read data from matrix market format file
  // Input: input file stream
  void read(std::istream &in); 

  // helper function for read function 
  void readData(std::istream &in, size_t num_rows, size_t num_cols, bool IsSymmetric);

  // read sparse matrix from MM.
  void readDataSparse(std::istream &in, size_t nnz, bool IsSymmetric);

  // write data from matrix market format file
  void write(std::ostream &out);  
};

void loadMarket(Matrix&, const std::string& );
void loadMarket(Vector&, const std::string& );

#endif
