#ifndef PARAMETERS_H
#define PARAMETERS_H

#include<string>

class Params
{
public:
	double step_size;
	int max_itrs;
	int total_num_threads;
	bool use_controller;
	double tmac_step_size;
	int block_size;
	int problem_size;
	std::string worker_type;
	std::string step_size_rule;
	bool async;

	Params()
	{
		// do nothing!
	}
};

#endif // !PARAMETERS_H
