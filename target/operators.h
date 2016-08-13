#ifndef OPERATORS_H
#define OPERATORS_H

#include"Matrix.h"
#include"Vector.h"

//////////////////////////// Group 1 Operator ////////////////////////////
class prox_l1
{
public:
	double weight;
	double step_size;

	prox_l1()
	{
		weight = 0;
		step_size = 0;
	};
};

//////////////////////////// Group 2 Operator ////////////////////////////
class prox_huber
{
public:
	double step_size;
	double weight;
	double delta;

	prox_huber()
	{
		this->weight = 0;
		this->step_size = 0;
		this->delta = 0;
	}
};


//////////////////////////// Group 3 Operator ////////////////////////////
class prox_elastic_net
{
public:
	double step_size;
	double weight;
	double weight_2;

	prox_elastic_net()
	{
		this->step_size = 0;
		this->weight = 0;
		this->weight_2 = 0;
	}
};


//////////////////////////// Group 4 Operator ////////////////////////////
class proj_box
{
public:
	double step_size;
	double weight;
	Vector* lower;
	Vector* upper;

	proj_box()
	{
		this->step_size = 0;
		this->weight = 0;
		this->lower = nullptr;
		this->upper = nullptr;
	}
};

//////////////////////////// Group 5 Operator ////////////////////////////
class proj_l1_ball
{
public:
	double step_size;
	double weight;
	double radius;

	proj_l1_ball()
	{
		this->step_size = 0;
		this->weight = 0;
		this->radius = 0;
	}
};

//////////////////////////// Group 6 Operator ////////////////////////////
class proj_hyperplane
{
public:
	double weight;
	double step_size;
	Vector* a;
	double b;
	double ata;

	proj_hyperplane()
	{
		this->weight = 0;
		this->step_size = 0;
		this->a = nullptr;
		this->b = 0;
		this->ata = 0;
	}
};

/*
NOTE THAT ALL OPERATORS PAST THIS POINT SHOULD BE TEMPLATIZED AT SOME POINT
*/

//////////////////////////// Group 7 Operator ////////////////////////////
template <typename Mat>
class forward_grad_for_square_loss
{
public:
	double weight;
	double step_size;
	Mat* A;
	Vector* b;
	Vector* Atx;
	Mat* At;

	forward_grad_for_square_loss()
	{
		this->weight = 0;
		this->step_size = 0;
		this->A = nullptr;
		this->b = nullptr;
		this->Atx = nullptr;
		this->At = nullptr;
	}
};


//////////////////////////// Group 8 Operator ////////////////////////////
template <typename Mat>
class forward_grad_for_qp
{
public:
	double weight;
	double step_size;
	Mat* Q;
	Vector* c;

	forward_grad_for_qp()
	{
		this->weight = 0;
		this->step_size = 0;
		this->Q = nullptr;
		this->c = nullptr;
	}
};

//////////////////////////// Group 9 Operator ////////////////////////////
template <typename Mat>
class forward_grad_for_dual_svm
{
public:
	double weight;
	double step_size;
	Mat* A;
	Vector* Ax;
	Mat* At;

	forward_grad_for_dual_svm()
	{
		this->weight = 0;
		this->step_size = 0;
		this->A = nullptr;
		this->Ax = nullptr;
		this->At = nullptr;

	}
};

//////////////////////////// Group 10 Operator ////////////////////////////
template <typename Mat>
class linear_eqn_jacobi_operator
{
public:
	double weight;
	double step_size;
	Mat* A;
	Vector* b;

	linear_eqn_jacobi_operator()
	{
		this->weight = 0;
		this->step_size = 0;
		this->A = nullptr;
		this->b = nullptr;
	}
};

//////////////////////////// Group 11 Operator ////////////////////////////
template <typename Mat>
class forward_grad_for_huber_loss
{
public:
	double weight;
	double step_size;
	double delta;
	Mat* A;
	Vector* b;
	Vector* Atx;
	Vector temp;
	Mat* At;

	forward_grad_for_huber_loss()
	{
		this->weight = 0;
		this->step_size = 0;
		this->delta = 0;
		this->A = nullptr;
		this->b = nullptr;
		this->Atx = nullptr;
		// do nothing with Vector temp reference
		this->At = nullptr;
	}
};

#endif 