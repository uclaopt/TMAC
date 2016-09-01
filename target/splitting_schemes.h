#ifndef SPLITTING_SCHEMES_H
#define SPLITTING_SCHEMES_H

#include"Operators.h"
#include"Matrix.h"

template<typename forward, typename backward>
class ForwardBackwardSplitting
{
public:
	ForwardBackwardSplitting(Vector* x, forward f, backward b)
	{
		std::cout << "Forward Matrix's number of rows: " << f.Q->number_of_rows << std::endl;
		std::cout << "Forward Matrix's number of cols: " << f.Q->number_of_columns << std::endl;

		std::cout << std::endl;

		std::cout << "I-Coordinates for Forward Matrix's entries: " << std::endl;
		for (int i = 0; i < f.Q->number_of_nonzeros; ++i)
			std::cout << f.Q->row_coordinates[i] << " ";

		std::cout << std::endl;

		std::cout << "J-Coordinates for Forward Matrix's entries: " << std::endl;
		for (int i = 0; i < f.Q->number_of_nonzeros; ++i)
			std::cout << f.Q->column_coordinates[i] << " ";

		std::cout << std::endl;
		std::cout << std::endl;
	}
};

#endif
