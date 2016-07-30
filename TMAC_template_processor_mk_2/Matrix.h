#ifndef MATRIX_H
#define MATRIX_H

#include<string>

class SpMat
{
private:
	std::string filename;

public:
	SpMat()
	{
		//do nothing
	}

	SpMat(std::string filename)
	{
		this->filename = filename;
	}

	void set_filename(std::string filename)
	{
		this->filename = filename;
	}
};


#endif 