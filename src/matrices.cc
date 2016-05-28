/*************************************************
Copyright (C) 2006-2011 Evgenii Rudnyi,
http://Evgenii.Rudnyi.Ru/
Implementations for matrix data file IO.
*/
#include <fstream>
#include <iostream>
#include <iomanip>
#include <set>
#include "matrices.h"
#include "MarketIO.h"
using namespace std;

// global functions for Matrix Market
bool isDefined ( const set<string> &opt, const string &name ) {
  if ( opt.find ( name ) == opt.end() ) {
    return false;
  }
  else {
    return true;	
  }
}

void readHeader ( istream &in, set<string> &opt ) {
  string line;
  getline ( in, line );
    if ( line.substr ( 0, 15 ) != "%%MatrixMarket " ) {
    throw MatrixError::ReadError ( "header: Not a Matrix Market file" );
    }
  istringstream in2 ( line.substr(15) );
  string word;
  while (in2) {
    in2 >> word;
      for ( unsigned i = 0; i < word.size(); ++i ) {
      word[i] = tolower ( word[i] );
      }
    opt.insert ( word );
  }
  while ( in.peek() == '%' ) {
    getline ( in, line );
  }
  if ( !isDefined ( opt, "matrix" ) ) {
    throw MatrixError::ReadError ( "header: Not a matrix" );
  }
  if ( !isDefined ( opt, "real" ) ) {
    throw MatrixError::ReadError ( "header: Not a real matrix" );
  }
}

void writeComment ( ostream &out ) {
  out << "% File is created http://MatrixProgramming.com/" << endl;
}

// member functions
void Vector::read ( istream &in ) {
  clear();
  set<string> opt;
  bool IsSymmetric;	
  readHeader ( in, opt );
  if ( isDefined ( opt, "symmetric" ) ) {
    IsSymmetric = true;
  }
  else if ( isDefined ( opt, "general" ) ) {
    IsSymmetric = false;
  }
  else {
    throw MatrixError::ReadError ( "Matrix: Not general, neither symmetric" );
  }
  size_t n, m;
  size_t nnz = 0;
  in >> n >> m;
    if ( IsSymmetric && ( m != n ) ) {
    throw MatrixError::ReadError ( stringAndTwoNumbers ( "Matrix: Symmetric but number of columns is not equal to number of rows", n, m ) );
    }
  resize ( n );
  
  if ( isDefined ( opt, "array" ) ) {
    readData ( in, n, IsSymmetric );
  }
  else {
    throw MatrixError::ReadError ( "Matrix: Not array, neither coordinate" );
  }
}

void Vector::readData ( istream &in, size_t n, bool IsSymmetric ) {
  double val;
  for ( unsigned j = 0; j < n; ++j ) {
    in >> val;
      if ( !in ) {
        throw MatrixError::ReadError ( "Matrix: End of file too early" );
      }
    ( *this )[j] = val;
  }
}

// member functions
void Matrix::read ( istream &in ) {
  clear();
  set<string> opt;
  bool IsSymmetric;	
  readHeader ( in, opt );
  if ( isDefined ( opt, "symmetric" ) ) {
    IsSymmetric = true;
  }
  else if ( isDefined ( opt, "general" ) ) {
    IsSymmetric = false;
  }
  else {
    throw MatrixError::ReadError ( "Matrix: Not general, neither symmetric" );
  }
  size_t n, m;
  size_t nnz = 0;
  in >> n >> m;
    if ( IsSymmetric && ( m != n ) ) {
    throw MatrixError::ReadError ( stringAndTwoNumbers ( "Matrix: Symmetric but number of columns is not equal to number of rows", n, m ) );
    }
  resize ( n, m );
  
  if ( isDefined ( opt, "array" ) ) {
    readData ( in, n, m, IsSymmetric );
  }
  else if ( isDefined ( opt, "coordinate" ) ) {
    in >> nnz;
    readDataSparse ( in, nnz, IsSymmetric );
  }
  else {
    throw MatrixError::ReadError ( "Matrix: Not array, neither coordinate" );
  }
}


// read data of size nxm
void Matrix::readData ( istream &in, size_t n, size_t m, bool IsSymmetric ) {
  double val;
  if ( IsSymmetric ) {
      for ( size_t j = 0; j < n; ++j ) {
        for ( size_t i = j; i < n; ++i ) {
          in >> val;
          if ( !in ) {
            throw MatrixError::ReadError ( "Matrix: Symmetric and end of file too early" );
          }
          ( *this ) ( j, i ) = ( *this )( i, j ) = val;
        }
      }
  }
  else {
    for ( size_t i = 0; i < m; ++i ) {
      for ( size_t j = 0; j < n; ++j ) {
        in >> val;
        if ( !in ) {
          throw MatrixError::ReadError ( "Matrix: End of file too early" );
        }
        ( *this ) ( j, i ) = val;
      }
    }
  }
  return;
}

// read matrix market file in sparse matrix format, and save as a
// dense matrix
void Matrix::readDataSparse ( istream &in, size_t nnz, bool IsSymmetric ) {
  double val;
  size_t i, j;
  fill ( begin(), end(), 0. );
  for ( size_t k = 0; k < nnz; ++k ) {
    in >> i >> j >> val;
      if ( !in ) {
      throw MatrixError::ReadError ( "Matrix: End of file too early" );
      }
    ( *this ) ( i - 1, j - 1 ) = val;
      if ( IsSymmetric ) {
      ( *this ) ( j - 1, i - 1 ) = val;
      }
  }
}

// write dense matrix to file with matrix
// market format
void Matrix::write ( ostream &out ) {
  out << "%%MatrixMarket matrix array real general" << endl;
  writeComment ( out );
  out << rows() << " " << cols() << endl;
    for ( vector <double>::const_iterator i = begin(); i < end(); ++i ) {
    out << *i << endl;
    }
}

// load data from matrix market format file to dense matrix
void loadMarket ( Matrix &A, const std::string& file_name ) {
  ifstream myfile ( file_name.c_str() );
    if ( !myfile ) {
    cout<<"cannot open data file"<<endl;
    }
  A.read ( myfile );
  myfile.close();
  return;
}

// load data from matrix market format to dense vector
void loadMarket ( Vector &A, const std::string& file_name ) {
  ifstream myfile ( file_name.c_str() );
    if ( !myfile ) {
    cout<<"cannot open data file"<<endl;
    }
  A.read ( myfile );
  myfile.close();
  return;
}
