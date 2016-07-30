#ifndef VECTOR_H
#define VECTOR_H

#include<vector>
#include<string>
#include<fstream>
#include<iostream>


class Vector {
private:
	std::vector<double> v;
	std::string filename;
public:
	Vector(size_t size = 0, double val = 0) :v(size, val) {}

	Vector(std::string filename) {

		this->filename = filename;

		std::ifstream in(filename);
		double val;
		while (in >> val) {
			v.push_back(val);
		}
		in.close();
	}

	// does everything the copy constructor does
	void set_filename(std::string filename)
	{
		this->filename = filename;

		// NOTE: we have to clear the vector first...
		v.clear();

		std::ifstream in(filename);
		double val;
		while (in >> val) {
			v.push_back(val);
		}
		in.close();
	}

	void print(std::ostream& out) const {
		for (auto& val : v) {
			out << val << " ";
		}
		std::cout << "\n";
	}

	Vector& operator*=(double scalar) {
		for (auto& val : v) {
			val *= scalar;
		}
		return *this;
	}
};

Vector operator*(double scalar, Vector v) {
	v *= scalar;
	return v;
}


std::ostream& operator<<(std::ostream& out, const Vector& v) {
	v.print(out);
	return out;
}


#endif
