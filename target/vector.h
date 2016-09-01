#ifndef VECTOR_H
#define VECTOR_H

#include<vector>
#include<string>
#include<fstream>
#include<iostream>

// VECTOR CLASS
class Vector
{
public:
	// Member variables of vector
	std::vector<double> v;
	std::string filename;
	const char* filename_in_c_arr;

	// Row and column info
	int number_of_rows;
	int number_of_columns;
	int number_of_nonzeros;

	Vector(size_t size = 0, double val = 0) :v(size, val) {}

	Vector(std::string filename)
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

		// Do we just treat a vector as a dense matrix???
		mm_read_mtx_array_size(f, &number_of_rows, &number_of_columns);
	}

	void print(std::ostream& out) const
	{
		for (auto& val : v)
		{
			out << val << " ";
		}
		std::cout << "\n";
	}

	Vector& operator*=(double scalar)
	{
		for (auto& val : v)
		{
			val *= scalar;
		}
		return *this;
	}
};

Vector operator*(double scalar, Vector v)
{
	v *= scalar;
	return v;
}


std::ostream& operator<<(std::ostream& out, const Vector& v)
{
	v.print(out);
	return out;
}


#endif
