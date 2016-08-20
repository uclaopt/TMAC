#ifndef TMAC_INCLUDE_RESULT_H
#define TMAC_INCLUDE_RESULT_H

struct Result {

	int status;                    // status = 0 means "Solving"; status = 1 means "Solved"; status = 2 means "Primal Infeasible"; status = 3 means "Unboundedness".
   	double optimal;
	double gap;


	Result() : status(0), optimal(0.), gap(0.){}

};
#endif
