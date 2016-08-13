#ifndef MATRIX_H
#define MATRIX_H

#include<string>
#include<stdio.h>
#include<iostream>
#include"mmio.h"

// Dense matrix class
// NEED TO BE ABLE TO READ IN DATA
class Matrix
{
public:
	// General member variables
	std::string filename;
	const char* filename_in_c_arr;

	// Row and column info
	int* rows = new int(0);
	int* columns = new int(0);

	// Matrix: std::string constructor
	Matrix(std::string filename)
	{
		this->filename = filename;
		filename_in_c_arr = filename.c_str();

		// Set up for mm_read_banner
		FILE* f;
		MM_typecode m;

		// First open up the file (remember to put "file___name.mtx.gz"
		if ((f = fopen(filename_in_c_arr, "r")) == nullptr)
		{
			std::cout << "Could not open file at specified location." << std::endl;
			exit(1);
		}

		// Second, read in mtx banner
		if (mm_read_banner(f, &m) != 0)
		{
			std::cout << "Could not successfully read mtx banner." << std::endl; // getting a premature end of file error
			exit(1);
		}

		// Third, make sure that we are dealing with a dense matrix
		if (mm_is_dense(m))
		{
			// now we actually start storing stuff...
			mm_read_mtx_array_size(f, rows, columns);
		}
	}
};

// SPARSE MATRIX CLASS
class SpMat
{
public:
	// General member variables
	std::string filename;
	const char* filename_in_c_arr;

	// Row and column info
	int number_of_rows;
	int number_of_columns;
	int number_of_nonzeros;

	// Row and column arrays
	int row_coordinates[100];
	int column_coordinates[100];
	double values[100 * 100];

	// Matrix: default constructor
	SpMat()
	{
		//do nothing
	}

	// Matrix: std::string constructor
	SpMat(std::string filename)
	{
		this->filename = filename;
		filename_in_c_arr = filename.c_str();

		// Set up for mm_read_banner
		FILE* f;
		MM_typecode m;

		// First open up the file (remember to put "file___name.mtx.gz"
		if ((f = fopen(filename_in_c_arr, "r")) == nullptr)
		{
			std::cout << "Could not open file at specified location." << std::endl;
			exit(1);
		}

		// Second, read in mtx banner
		if (mm_read_banner(f, &m) != 0)
		{
			std::cout << "Could not successfully read mtx banner." << std::endl; // getting a premature end of file error
			exit(1);
		}

		// Third, make sure that we are dealing with a sparse matrix
		if (mm_is_sparse(m))
		{
			// now we actually start storing stuff...
			mm_read_mtx_crd_size(f, &number_of_rows, &number_of_columns, &number_of_nonzeros);
		}

		// Fourth, read in the data... (IN CONSTRUCTION)...
		mm_read_mtx_crd_data(f, number_of_rows, number_of_columns, number_of_nonzeros, 
			row_coordinates, column_coordinates, values, m);
	}
};


#endif