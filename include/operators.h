/**************************************************************
 * Single operator functors used by the TMAC framework. It 
 * includes the following three types of operators:
 *
 *   1) proximal operator;
 *   2) projection to a simple set;
 *   3) forward operator for common loss functions.
 *
 *************************************************************/

#ifndef TMAC_INCLUDE_OPERATORS_H
#define TMAC_INCLUDE_OPERATORS_H
#include "constants.h"
#include "matrices.h"  // matrices library
#include "algebra.h"   // linear algebra library
#include "util.h"      // utility functions
#include <math.h>      // for fabs
#include "algebra_namespace_switcher.h"

/***************************
 * base class for operators
 ***************************/
class OperatorInterface {
public:
  virtual double operator() (Vector* v, int index = 0) {}
  virtual double operator() (double val, int index = 0) {}
  virtual void operator() (Vector* v_in, Vector* v_out) {}
  virtual void update_step_size(double step_size_) = 0;
  virtual void update_cache_vars(double old_x_i, double new_x_i, int index) = 0;
  virtual void update_cache_vars(Vector* x, int rank, int num_threads) = 0;
};


/**************************************************************
 *                  Proximal operators                        *
 **************************************************************
 * It is defined by the following:                            *
 * prox_f (x) = argmin_{y} f(y) + 1/(2 step_size) ||y - x||^2 *
 **************************************************************/


/***************************************************
 * Proximal l1 functor for the following function
 *   f(x) = weight ||x||_1,
 * which is given by
 *   prox_l1 = soft_threshold(x, weight * step_size);
 ***************************************************/
class prox_l1 : public OperatorInterface {
public:
  double weight;
  double step_size;
  double operator() (Vector* v, int index = 0) {
    double val;
    double threshold = step_size * weight;
    val = (*v)[index];
    if (val > threshold) {
      return val - threshold;
    } else if (val < -threshold) {
      return val + threshold;
    } else {
      return 0.;
    }
  }
  double operator() (double val, int index = 0) {
    double threshold = step_size * weight;
    if (val > threshold) {
      return val - threshold;
    } else if (val < -threshold) {
      return val + threshold;
    } else {
      return 0.;
    }
  }
  // full update operator
  void operator() (Vector* v_in, Vector* v_out) {
    int len = v_in->size();
    double threshold = step_size * weight;
    double val;
    for (int i = 0; i < len; i++) {
      val = (*v_in)[i];
      if (val > threshold) {
        (*v_out)[i] = val - threshold;
      } else if (val < -threshold) {
        (*v_out)[i]= val + threshold;
      } else {
        (*v_out)[i] = 0.;
      }
    }
  }
  void update_cache_vars(double old_x_i, double new_x_i, int index) {
  }
  void update_cache_vars(Vector* x, int rank, int num_threads) {
  }
  
  void update_step_size(double step_size_) {
    step_size = step_size_;
  }
  
  prox_l1 (double step_size_, double weight_ = 1.) {
    step_size = step_size_;
    weight = weight_;
  }
  prox_l1 () {
    step_size = 0.;
    weight = 1.;
  }
};


/***************************************************
 * Proximal sum square (2 norm square)
 *   f(x) = weight/2 ||x||_2^2,
 * which is given by
 *   prox_f = v / (1 + weight * step_size)
 ***************************************************/
class prox_sum_square : public OperatorInterface {
public:
  double weight;
  double step_size;
  
  double operator() (Vector* v, int index) {
    double val;
    val = (*v)[index];
    return val / (1. + weight * step_size);
  }

  double operator() (double val, int index = 0) {
    return val / (1. + weight * step_size);
  }

  void operator() (Vector* v_in, Vector* v_out) {
    // if debug, we check index is valid
    int len = v_in->size();
    for (int i = 0; i < len; i++) {
      (*v_out)[i] = (*v_in)[i] / (1. + weight * step_size);
    }
    return;
  }

  void update_step_size(double step_size_) {
    step_size = step_size_;
  }
  void update_cache_vars(double old_x_i, double new_x_i, int index) {
  }
  void update_cache_vars(Vector* x, int rank, int num_threads) {
  }
  
  prox_sum_square (double step_size_, double weight_ = 1.) {
    step_size = step_size_;
    weight = weight_;
  }
  prox_sum_square () {
    step_size = 0.;
    weight = 1.;
  }
};


/******************************************************
 * Proximal of the two norm  function
 *  f(x) = weight ||x||_2,
 * which is defined by
 *  y, where
 *     y = 0, if ||x|| <= weight * step_size
 *     y = (1 - weight * step_size / norm_x) * x.
 *
 * WARNING: 
 *   This is only suited for Backward-forward splitting,
 * where the forward operator does not have maintained
 * variables, (otherwise,coordinate update for the
 * forward update will be expensive.)
 ******************************************************/
class prox_l2 : public OperatorInterface {
public:  
  double weight;
  double step_size;
  
  double operator() (Vector* v, int index) {
    double norm_v = norm(*v, 2);
    int len = v->size();
    double scale = 1. - weight * step_size / norm_v;
    if (norm_v > step_size * weight) {
      return (*v)[index] * scale;
    } else {
      return 0.;
    }
  }

  double operator() (double val, int index = 0) {
    return DOUBLE_MARKER;
  }
  
  void operator() (Vector* v_in, Vector* v_out) {
    double norm_v = norm(*v_in, 2);
    int len = v_in->size();
    double scale = 1. - weight * step_size / norm_v;
    if (norm_v > step_size * weight) {
      for (int i = 0; i < len; i++) {
        (*v_out)[i] = scale * (*v_in)[i];
      }
    } else {
      for (int i = 0; i < len; i++) {
        (*v_out)[i] = 0.;
      }
    }
  }

  void update_step_size(double step_size_) {
    step_size = step_size_;
  }
  void update_cache_vars(double old_x_i, double new_x_i, int index) {
  }
  void update_cache_vars(Vector* x, int rank, int num_threads) {
  }
    
  prox_l2 (double step_size_, double weight_ = 1.) {
    step_size = step_size_;
    weight = weight_;
  }
  prox_l2 () {
    step_size = 0.;
    weight = 1.;
  }
};


/******************************************************************
 * proximal of the huber function
 *  f(x) = 0.5 * weight * x^2, if -delta <= x = delta;
 *         delta * weight * (|x| - 0.5 * delta), otherwise.
 * is the following
 * prox_huber = x - threshold,               if x >= delta + threshold
 *              x/(1 + step_size * weight),   otherwise
 *              x + threshold,               if x <= -delta - threshold
 * where threshold = delta * weight * step_size.
 ******************************************************************/
class prox_huber : public OperatorInterface {
public:  
  double delta;
  double weight;
  double step_size;
  
  double operator() (Vector* v, int index = 0) {
    // if debug, we check index is valid
    double val;
    double threshold = step_size * weight * delta;
    val = (*v)[index];
    if (val > delta + threshold) {
      return val - threshold;
    } else if (val < -delta - threshold) {
      return val + threshold;
    } else {
      return val / (1. + step_size * weight);
    }
  }
  
  double operator() (double val) {
    double threshold = step_size * weight * delta;
    if (val > delta + threshold) {
      return val - threshold;
    } else if (val < -delta - threshold) {
      return val + threshold;
    } else {
      return val / (1. + step_size * weight);
    }
  }
  
  void operator() (Vector* v_in, Vector* v_out) {
    // if debug, we check index is valid
    int len = v_in->size();
    double threshold = step_size * weight * delta;
    double val;
    for (int i = 0; i < len; i++) {
      val = (*v_in)[i];
      if (val > delta + threshold) {
        (*v_out)[i] = val - threshold;
      } else if (val < -delta - threshold) {
        (*v_out)[i]= val + threshold;
      } else {
        (*v_out)[i] = val / (1. + step_size * weight);
      }
    }
  }

  void update_step_size(double step_size_) {
    step_size = step_size_;
  }  
  void update_cache_vars(double old_x_i, double new_x_i, int index) {
  }
  void update_cache_vars(Vector* x, int rank, int num_threads) {
  }
  
  prox_huber (double step_size_, double weight_ = 1.,
              double delta_ = 1.) {
    step_size = step_size_;
    weight = weight_;
    delta = delta_;
  }

  prox_huber () {
    step_size = 0.;
    weight = 1.;
    delta = 1.;
  }
};


/*************************************************************************
 * proximal of the log barrier function
 *   f(x) = - weight * sum(log(x))
 * which is given by
 *  prox_log_barrier(x) = 0.5 * (x + sqrt(x*x + 4 * weight * step_size))
 *************************************************************************/
class prox_log_barrier : public OperatorInterface {
public:
  double weight;
  double step_size;
  double operator() (Vector* v, int index = 0) {
    double val = (*v)[index];
    return 0.5 * (val + sqrt(val * val + 4. * weight * step_size));
  }
  
  double operator() (double val) {
    return 0.5 * (val + sqrt(val * val + 4. * weight * step_size));
  }
  
  
  void operator() (Vector* v_in, Vector* v_out) {
    int len = v_in->size();
    double val;
    for (int i = 0; i < len; i++) {
      val = (*v_in)[i];
      (*v_out)[i] = 0.5 * (val + sqrt(val * val + 4. * weight * step_size));
    }
  }

  void update_step_size(double step_size_) {
    step_size = step_size_;
  }
  void update_cache_vars(double old_x_i, double new_x_i, int index) {
  }
  void update_cache_vars(Vector* x, int rank, int num_threads) {
  }
  
  prox_log_barrier (double step_size_, double weight_ = 1.) {
    step_size = step_size_;
    weight = weight_;
  }
  prox_log_barrier () {
    step_size = 0.;
    weight = 1.;
  }
};


/*******************************************************
 * proximal of elastic net function
 *  f(x) = weight ||x||_1 + weight_2/2 ||x||_2^2
 * is the following function
 *  prox_elastic_net = 1/(1 + weight_2 * step_size)
 *                * shrink(x, weight_1 * step_size)
 *******************************************************/
class prox_elastic_net : public OperatorInterface {
public:
  double weight_2;
  double weight;
  double step_size;
  
  double operator() (Vector* v, int index = 0) {
    // if debug, we check index is valid
    double val;
    double threshold = step_size * weight;
    double scale = 1. / (1. + step_size * weight_2);
    val = (*v)[index];
    if (val > threshold) {
      return scale * (val - threshold);
    } else if (val < -threshold) {
      return scale * (val + threshold);
    } else {
      return 0.;
    }
  }
  
  double operator() (double val) {
    double threshold = step_size * weight;
    double scale = 1. / (1. + step_size * weight_2);
    if (val > threshold) {
      return scale * (val - threshold);
    } else if (val < -threshold) {
      return scale * (val + threshold);
    } else {
      return 0.;
    }
  }

  void operator() (Vector* v_in, Vector* v_out) {
    // if debug, we check index is valid
    int len = v_in->size();
    double threshold = step_size * weight;
    double scale = 1. / (1. + step_size * weight_2);    
    for (int i = 0; i < len; i++) {
      double val = (*v_in)[i];
      if (val > threshold) {
        (*v_out)[i] = scale * (val - threshold);
      } else if (val < -threshold) {
        (*v_out)[i]= scale * (val + threshold);
      } else {
        (*v_out)[i] = 0.;
      }
    }
  }

  void update_step_size(double step_size_) {
    step_size = step_size_;
  }  
  void update_cache_vars(double old_x_i, double new_x_i, int index) {
  }
  void update_cache_vars(Vector* x, int rank, int num_threads) {
  }
  
  prox_elastic_net (double step_size_, double weight_ = 1., double weight_2_ = 1.) {
    step_size = step_size_;
    weight = weight_;
    weight_2 = weight_2_;
  }
  prox_elastic_net () {
    step_size = 0.;
    weight = 1.;
    weight_2 = 1.;
  }
};

/**********************************************************************
*
* proximal of linear product on a hyperplane
*  f(x) = weight c'x + indicate (Ax = b)
*  prox_linearpro_hyperplane = Tx + d
*  where T = I - A'(AA')^(-1)A, 
*        d =  A'(AA')^(-1)(b + gamma * weight * Ac) - gamma * weight * c.
*  to simplify calculation, load D = A'(AA')^(-1)A; 
*                                f = A'(AA')^(-1)b; 
*                                g = A'(AA')^(-1)Ac;
*  from interface, then Tx + d = x - Dx + f + gamma * weight * (g - c)
**********************************************************************/
class prox_linearpro_hyperplane : public OperatorInterface {
public:
	Vector* c;
	Vector* f;
	Vector* g;
	Matrix* D;
	double gamma;
	double weight;
	double step_size;

	double operator() (Vector* v, int index) {
		int len = v->size();
		double product = 0;
		for (int i = 0; i < len; i++) {
			product += (*D)(index, i) * (*v)[i];
		}
		return (*v)[index] - product + (*f)[index] + gamma * weight * ((*g)[index] - (*c)[index]);
	}

	void operator() (Vector* v_in, Vector* v_out) {
		
     }

	void update_cache_vars(double old_x_i, double new_x_i, int index) {
	}
	void update_cache_vars(Vector* x, int rank, int num_threads) {
	}

	void update_step_size(double step_size_) {
		step_size = step_size_;
	}

	prox_linearpro_hyperplane() {
		step_size = 0.;
		weight = 1.;
		gamma = 1.;
		c = nullptr;
		f = nullptr;
		g = nullptr;
		D = nullptr;
	}
	prox_linearpro_hyperplane(Vector* c_, Vector* f_, Vector* g_, Matrix* D_, double gamma_, double weight_, double step_size_ = 1.) {
		c = c_;
		f = f_;
		g = g_;
		D = D_;
		gamma = gamma_;
		step_size = step_size_;
		weight = weight_;
	}

};

/*************************************************
 *      Projection operators       
 * We can view projection as proximal operator
 * of the corresponding indicate function. In this
 * case, weight and step_size are not being used.
 *************************************************/


// projection to the first-orthant cone
//   f(x) = I(x >=0)
class proj_positive_cone : public OperatorInterface {
public: 
  double weight;
  double step_size;
  double operator() (Vector* v, int index) {
    return fmax((*v)[index], 0.);
  }
  
  double operator() (double val, int index = 0) {
    return max(val, 0.);
  }
  
  void operator() (Vector* v_in, Vector* v_out) {
    int len = v_in->size();
    for (int i = 0; i < len; i++) {
      (*v_out)[i] = fmax((*v_in)[i], 0.);
    }
  }

  void update_step_size(double step_size_) {
    step_size = step_size_;
  }
  void update_cache_vars(double old_x_i, double new_x_i, int index) {
  }
  void update_cache_vars(Vector* x, int rank, int num_threads) {
  }
  
  proj_positive_cone () {
    step_size = 0.;
    weight = 1.;
  }
  proj_positive_cone (double step_size_, double weight_ = 1.) {
    step_size = step_size_;
    weight = weight_;
  }
};


// projection to box
// f = I(l <= x <= u)
class proj_box : public OperatorInterface {
public:
  Vector *lower, *upper;
  double step_size;
  double weight;
  double operator() (Vector* v, int index) {
    double val = (*v)[index];
    return fmax((*lower)[index], fmin(val, (*upper)[index]));
  }
  
  double operator() (double val, int index = 0) {
    return fmax((*lower)[index], fmin(val, (*upper)[index]));
  }
  
  void operator() (Vector* v_in, Vector* v_out) {
    int len = v_in->size();
    double val;
    for (int i = 0; i < len; i++) {
      val = (*v_in)[i];
      (*v_out)[i] = fmax((*lower)[i], fmin(val, (*upper)[i]));
    }
  }
  void update_cache_vars(double old_x_i, double new_x_i, int index) {
  }
  void update_cache_vars(Vector* x, int rank, int num_threads) {
  }

  void update_step_size(double step_size_) {
    step_size = step_size_;
  }
  
  proj_box () {
    step_size = 0.;
    weight = 1.;
    lower = nullptr;
    upper = nullptr;
  }
  proj_box (Vector* lower_, Vector* upper_, double step_size_ = 1., double weight_ = 1.) :
      lower(lower_), upper(upper_), step_size(step_size_), weight(weight_) {}
};

/*****************************************************************
*
* projection to second order cone: ||x||_2 <= x_0 is
* defined by
*   proj_soc(x_0, x) =(x_0, x), if ||x||_2 <= x_0,
*                    =(0, 0), if ||x||_2 <= -x_0,
*                    =(1/2)(1+x_0/||x||_2)(||x||_2, x), otherwise.
*
******************************************************************/
class proj_soc : public OperatorInterface {
public:

	double weight;
	double step_size;
    double operator() (Vector* v, int index) {
		double rest_norm = 0.;
		for (int i = 1; i < v->size(); i++)
			rest_norm += pow((*v)[i], 2.0);
		rest_norm = sqrt(rest_norm);
		if ((*v)[0] >= rest_norm)
			return (*v)[index];
		else if (-(*v)[0] >= rest_norm)
			return 0;
		else {
			if (index == 0)
				return 0.5 * ((*v)[0] + rest_norm);
			else
				return 0.5 * (1 + (*v)[0] / rest_norm) * (*v)[index];
		}
	}

	void operator() (Vector* v_in, Vector* v_out) {
		double rest_norm = 0.;
		int len = v_in->size();
		for (int i = 1; i < len; i++)
			rest_norm += pow((*v_in)[i], 2.0);
		rest_norm = sqrt(rest_norm);
		if ((*v_in)[0] >= rest_norm) {
			for (int i = 0; i < len; i++)
				(*v_out)[i] = (*v_in)[i];
		}
		else if (-(*v_in)[0] >= rest_norm) {
			for (int i = 0; i < len; i++)
				(*v_out)[i] = 0;
		}
		else {
			(*v_out)[0] = 0.5 * ((*v_in)[0] + rest_norm);
			for (int i = 1; i < len; i++)
				(*v_out)[i] = 0.5 * (1 + (*v_in)[0] / rest_norm) * (*v_in)[i];
		}

	}
	void update_step_size(double step_size_) {
		step_size = step_size_;
	}
	void update_cache_vars(double old_x_i, double new_x_i, int index) {
	}
	void update_cache_vars(Vector* x, int rank, int num_threads) {
	}

	// struct constructors
	proj_soc() {
		step_size = 0.;
		weight = 1.;
	}
	proj_soc(double step_size_, double weight_ = 1.) {
		step_size = step_size_;
		weight = weight_;
	}

};

/*****************************************************************
*
* projection to multiple second order cones: ||xi||_2 <= xi_0 is
* i = 1,2,...,l
* projection to multiple second order cones with linear part means:
* proj_multi_soc of x + weight * c with input x
******************************************************************/
class proj_multi_soc_with_linear : public OperatorInterface {
public:
	int cone_num;
	Vector cone_dim;
	Vector* c;
	double weight;
	double step_size;
    double operator() (Vector* v, int index) {
		//cout << "projection socp index" << endl;
		// locate position
		int cone_id = 0;
		int dim_sum = 0;
		int len = v->size();
		for (int i = 0; i < len; i++) {
			dim_sum += cone_dim[i];
			if (dim_sum > index)
				cone_id = i;
		}
		int cone_end = dim_sum;
		int cone_start = dim_sum - cone_dim[cone_id];
		
		double rest_norm = 0.;
		for (int i = cone_start + 1; i < cone_end; i++)
			rest_norm += pow((*v)[i] + weight * (*c)[i], 2.0);
		rest_norm = sqrt(rest_norm);
		
		if ((*v)[cone_start] + weight * (*c)[cone_start] >= rest_norm)
			return (*v)[index] + weight * (*c)[index];
		
		else if (-(*v)[cone_start] - weight * (*c)[cone_start] >= rest_norm)
			return 0;
		
		else {
			if (index == cone_start)
				return 0.5 * ((*v)[cone_start] + weight * (*c)[cone_start] + rest_norm);
			else
				return 0.5 * (1 + ((*v)[cone_start] + weight * (*c)[cone_start]) / rest_norm) * ((*v)[index]+ weight * (*c)[index]);
		}
	}

	void operator() (Vector* v_in, Vector* v_out) {
		//cout << "projection socp vector" << endl;
		int cone_start = 0;
		int cone_end = 0;
		double rest_norm = 0.;
		for (int i = 0; i < cone_num; i++) {
			cone_end += cone_dim[i];
			for (int j = cone_start + 1; j < cone_end; j++)
				rest_norm += pow((*v_in)[j] + weight * (*c)[j], 2.0);
			rest_norm = sqrt(rest_norm);
			
			if ((*v_in)[cone_start] + weight * (*c)[cone_start] >= rest_norm) {
				for (int t = cone_start; t < cone_end; t++)
					(*v_out)[t] = (*v_in)[t] + weight * (*c)[t];
			}
			
			else if (-(*v_in)[cone_start] - weight * (*c)[cone_start] >= rest_norm) {
				for (int t = cone_start; t < cone_end; t++)
					(*v_out)[t] = 0;
			}
			
			else {
				(*v_out)[cone_start] = 0.5 * ((*v_in)[cone_start] + weight * (*c)[cone_start] + rest_norm);
				for (int t = cone_start + 1; t < cone_end; t++)
					(*v_out)[t] = 0.5 * (1 + ((*v_in)[cone_start] + weight * (*c)[cone_start]) / rest_norm) * ((*v_in)[t] + weight * (*c)[t]);
			}
			cone_start += cone_dim[i];
			rest_norm = 0;
		}
	}
	void update_step_size(double step_size_) {
		step_size = step_size_;
	}
	void update_cache_vars(double old_x_i, double new_x_i, int index) {
	}
	void update_cache_vars(Vector* x, int rank, int num_threads) {
	}

	// struct constructors
	proj_multi_soc_with_linear() {
		cone_num = 0;
		cone_dim = NULL;
		step_size = 0.;
		weight = 1.;
	}
	proj_multi_soc_with_linear(int cone_num_, Vector cone_dim_, Vector* c_, double weight_, double step_size_=1.) {
		cone_num = cone_num_;
		cone_dim = cone_dim_;
		c = c_;
		weight = weight_;
	}

};

/* projection to l1 ball
 *    min   ||w - v||_2
 *    s.t.  ||w||_1 <= radius.
 *   according to [matlab code](http://stanford.edu/~jduchi/projects/DuchiShSiCh08/ProjectOntoL1Ball.m)
 */
//TODO:  add a maintain ||x|| version
class proj_l1_ball : public OperatorInterface {
public:
  double radius;
  double weight;
  double step_size;
  
  // operator that takes a vector and an index
  // running time O(n\log n)
  // http://stanford.edu/~jduchi/projects/DuchiShSiCh08/ProjectOntoL1Ball.m
  // this running time can be improved by changing data structure rather than algorithm
  double operator() (Vector* v_in, int index) {
    double  nrm1_v = 0.0;
    int len = v_in->size(), i = 0;
    for (i = 0; i < len; i++) {
      nrm1_v += fabs( (*v_in)[i] );
    }
    if (radius >= nrm1_v) {
      return (*v_in)[index];
    }
    else {
      Vector u(len);
      for (i = 0; i < len; i++) {
        u[i] = fabs( (*v_in)[i] );
      }
      /* u = sort(abs(v),'descend'); -- MATLAB  */
      std::sort(u.begin(), u.end() ); // ascend
      std::reverse(u.begin(), u.end() ); // ascend -> descend

      /* sv = cumsum(u); -- MATLAB */
      Vector sv(len); // cumulative sum
      sv[0] = u[0];
      for (i = 1; i < len; i++) {
        sv[i] = sv[i-1] + u[i];
      }

      /* rho = find(u > (sv - b) ./ (1:length(u))', 1, 'last'); -- MATLAB*/
      int rho = 0;
      for (rho = len - 1; u[rho] * (rho + 1) <= sv[rho] - radius; rho--) {
      }
      
      /*theta= max(0, (sv(rho) - b) / rho); -- MATLAB  */
      double theta = (sv[rho] - radius) / (rho + 1) ; 
      theta = theta > 0 ? theta : 0;

      /* w = sign(v) .* max(abs(v) - theta, 0);  -- MATLAB */
      /* just get w[index] */
      double temp;
      temp = fabs( (*v_in)[index] ) - theta;
      temp = temp > 0 ? temp : 0;
      return ((*v_in)[index] > 0 ? temp : -temp) ;  //notice if 0==(*v_in)[index], then 0==temp
    }
  }
  
  void operator() (Vector* v_in, Vector* v_out) { 
    double  nrm1_v = 0.0;
    int len = v_in->size(), i = 0;
    for (i = 0; i < len; i++) {
      nrm1_v += fabs( (*v_in)[i] );
    }
    if (radius >= nrm1_v) {
      for (i = 0; i < len; i++) {
        (*v_out)[i] = (*v_in)[i];
      }
    }
    else {      
      Vector u(len);
      for (i = 0; i < len; i++) {
        u[i] = fabs( (*v_in)[i] );
      }
      /* u = sort(abs(v),'descend'); -- MATLAB */
      std::sort(u.begin(), u.end() ); 
      std::reverse(u.begin(), u.end() );
      /* sv = cumsum(u);-- MATLAB */
      Vector sv(len); // cumulative sum
      sv[0] = u[0];
      for (i = 1; i < len; i++) {
        sv[i] = sv[i-1] + u[i];
      }
      /* rho = find(u > (sv - b) ./ (1:length(u))', 1, 'last'); -- MATLAB*/
      int rho = 0;
      for (rho = len - 1; u[rho] * (rho + 1) <= sv[rho] - radius; rho--) {
      } 
      /*theta= max(0, (sv(rho) - b) / rho); -- MATLAB */
      double theta = (sv[rho] - radius) / (rho + 1); 
      theta = theta > 0 ? theta : 0;
      /* w = sign(v) .* max(abs(v) - theta, 0); -- MATLAB  */
      double temp;
      for (i = 0; i < len; i++) {  
        temp = fabs( (*v_in)[i] ) - theta;
        temp = temp > 0 ? temp : 0;
        (*v_out)[i] = (*v_in)[i] > 0 ? temp : -temp ;  //notice if 0==(*v_in)[i], then 0==temp
      }
    }
  }

  void update_step_size(double step_size_) {
    step_size = step_size_;
  }
  void update_cache_vars(double old_x_i, double new_x_i, int index) {
  }
  void update_cache_vars(Vector* x, int rank, int num_threads) {
  }
  
  double operator() (double val, int index = 1) {
    return DOUBLE_MARKER;
  }
  
  void update_radius(double radius_) {
    radius = radius_;
  }
  
  // struct constructors
  proj_l1_ball() {
    radius = 1.;
    step_size = 0.;
    weight = 1.;
  }
  proj_l1_ball(double radius_, double step_size_=1., double weight_=1.) {
    radius = radius_;
    step_size = step_size_;
    weight = weight_;
  }
};


/*********************************************
 *
 * projection to l2 ball: ||x||_2 <= r is
 * defined by
 *   proj_l2(x) = x, if ||x|| <= r,
 *              = x * r/||x||
 *
 *********************************************/
// TODO: add a maintain ||x|| version
class proj_l2_ball : public OperatorInterface {
public:
  double radius;
  double weight;
  double step_size;
  
  double operator() (Vector* v, int index) {
    double nrm_v = norm(*v);
    double val = (*v)[index];
    return (nrm_v < radius) ? val : (val * radius / nrm_v);
  }
  
  void operator() (Vector* v_in, Vector* v_out) {
    int len = v_in->size();
    double nrm_v = norm(*v_in);
    double val;
    for (int i = 0; i < len; i++) {
      val = (*v_in)[i];
      (*v_out)[i] = (nrm_v < radius) ? val : (val * radius / nrm_v);
    }
  }

  double operator() (double val, int index = 1) {
    return DOUBLE_MARKER;
  }

  void update_step_size(double step_size_) {
    step_size = step_size_;
  }  
  void update_cache_vars(double old_x_i, double new_x_i, int index) {
  }
  void update_cache_vars(Vector* x, int rank, int num_threads) {
  }
  
  void update_radius (double radius_) {
    radius = radius_;
  }
  
  proj_l2_ball() {
    radius = 1.;
    step_size = 0.;
    weight = 1.;
  }
  proj_l2_ball (double radius_, double step_size_=1., double weight_=1.) {
    radius = radius_;
    step_size = step_size_;
    weight = weight_;
  }
};


/********************************************
 *
 * projection to a hyperplane
 *  a'x = b
 * is given by
 * proj_hyperplane = x + (b - a'x) / (a'a) * a
 *
 *******************************************/
class proj_hyperplane : public OperatorInterface {
public:
  Vector *a;
  double b;
  double ata;
  double weight;
  double step_size;
  // TODO: use atv as the maintained variable, and improve the efficiency from O(n) to O(1)
  double operator() (Vector* v, int index) {
    double val = (*v)[index];
    double residual = b - dot(*a, *v);
    return val + (residual / ata) * (*a)[index];
  }
  
  void operator() (Vector* v_in, Vector* v_out) {
    int len = v_in->size();
    double residual = b - dot(*a, *v_in);
    for (int i = 0; i < len; i++) {
      (*v_out)[i] = (*v_in)[i] + (residual / ata) * (*a)[i];
    }
  }

  double operator() (double val, int index = 1) {
    return DOUBLE_MARKER;
  }

  void update_step_size(double step_size_) {
    step_size = step_size_;
  }
  void update_cache_vars(double old_x_i, double new_x_i, int index) {
  }
  void update_cache_vars(Vector* x, int rank, int num_threads) {
  }
  
  proj_hyperplane() {
    step_size = 0.;
    weight = 1.;
    a = NULL;
    b = 0.;
    ata = 0.;
  }
  
  proj_hyperplane(Vector* a_, double b_, double step_size_ = 1., double weight_ = 1.) {
    step_size = step_size_;
    weight = weight_;
    a = a_;
    b = b_;
    ata = norm(*a, 2);
    ata *= ata;
  }
};

/* projection to probability simplex
 *    min   ||w - v||_2
 *    s.t.  w^T 1_1 = 1 , w >= 0.
 *   according to 
 */
//TODO:  add a maintain ||x|| version
class proj_prob_simplex : public OperatorInterface {
public:  
  double weight;
  double step_size;
  // operator that takes a vector and an index
  // running time O(n\log n)
  // Laurent Condat, “Fast projection onto the simplex and the L1 ball,” Mathematical Programming, pp. 1–11, 2015.
  // this running time can be improved by changing data structure rather than algorithm
  double operator() (Vector* v_in, int index) {
    int i, len = v_in->size();
    Vector u(len);
    for (i = 0; i < len; i++) {
      u[i] = (*v_in)[i];
    }
    /* Sort y into u: u1 >= u2 >= \ldots >= uN */
    std::sort(u.begin(), u.end() ); 
    std::reverse(u.begin(), u.end() );
    /* sv = cumsum(u);-- MATLAB */
    Vector sv(len); // cumulative sum
    sv[0] = u[0];
    for (i = 1; i < len; i++) {
      sv[i] = sv[i-1] + u[i];
    }
    /* Set K = max {k: (sv(k) - 1) / k < uk \} where 1 == a since here is a probability simplex. */
    /* K = find(u > (sv - 1) ./ (1:length(u))', 1, 'last'); -- MATLAB*/
    int K = 0;
    for (K = len - 1; u[K] * (K + 1) <= sv[K] - 1.; K--) {
    } 
    /* tau = (sv(k) - 1) / K */
    /* tau = (sv(K) - 1) / K; -- MATLAB */
    double tau = (sv[K] - 1.) / (K + 1); 
    /* for n=1, \ldots, N, set xn=max{yn-tau, 0} */
    double temp = (*v_in)[index] - tau;
    return temp > 0 ? temp : 0;
  }
  
  double operator() (double val, int index = 1) {
    return DOUBLE_MARKER;
  }
  
  void operator() (Vector* v_in, Vector* v_out) { 
    int i, len = v_in->size();
    Vector u(len);
    for (i = 0; i < len; i++) {
      u[i] = (*v_in)[i];
    }
    /* Sort y into u: u1 >= u2 >= \ldots >= uN */
    std::sort(u.begin(), u.end() ); 
    std::reverse(u.begin(), u.end() );
    /* sv = cumsum(u);-- MATLAB */
    Vector sv(len); // cumulative sum
    sv[0] = u[0];
    for (i = 1; i < len; i++) {
      sv[i] = sv[i-1] + u[i];
    }
    /* Set K = max {k: (sv(k) - 1) / k < uk \} where 1 == a since here is a probability simplex. */
    /* K = find(u > (sv - 1) ./ (1:length(u))', 1, 'last'); -- MATLAB*/
    int K = 0;
    for (K = len - 1; u[K] * (K + 1) <= sv[K] - 1.; K--) {
    } 
    /* tau = (sv(k) - 1) / K */
    /* tau = (sv(K) - 1) / K; -- MATLAB */
    double tau = (sv[K] - 1.) / (K + 1); 
    /* for n=1, \ldots, N, set xn=max{yn-tau, 0} */
    double temp;
    for (i = 0; i < len; i++) {  
      temp = (*v_in)[i] - tau;
      (*v_out)[i] = temp > 0 ? temp : 0;
    }
  }

  void update_step_size(double step_size_) {
    step_size = step_size_;
  }
  void update_cache_vars(double old_x_i, double new_x_i, int index) {
  }
  void update_cache_vars(Vector* x, int rank, int num_threads) {
  }
  
  proj_prob_simplex() {
    step_size = 0.;
    weight = 1.;
  }
  
  proj_prob_simplex(double step_size_, double weight_ = 1.) {
    step_size = step_size_;
    weight = weight_;
  }
};


/***********************************
 ***********************************
 *                                 *
 *      Gradient operators         *
 *                                 *
 ***********************************
 ***********************************/

// forward opertor for
// f(x) = weight / 2 ||A' x - b||^2, which is given by
// x - step_size * weight * A * (A' * x - b), 
// where A is a row major matrix, Atx is
// a maintained variable that stores the A'x
template <typename Mat>
class forward_grad_for_square_loss : public OperatorInterface {
public:
  Mat* A;
  Vector *b, *Atx;  
  Mat* At; // this is for the sync-parallel stuff
  double weight;
  double step_size;
  
  double operator() (Vector* x, int index) {
    // calculate the forward step
    double A_iAtx = dot(A, Atx, index);
    double A_ib = dot(A, b, index);
    double forward_grad_at_i = (*x)[index] - weight * step_size * (A_iAtx - A_ib);
    return forward_grad_at_i;
  }
  
  void operator() (Vector* v_in, Vector* v_out) {
    int m = A->rows(), n = A->cols();
    Vector Atv_in(m, 0.);
    trans_multiply(*A, *v_in, Atv_in);
    add(Atv_in, *b, -1.);
    Vector temp(n, 0.);
    multiply(*A, Atv_in, temp);
    scale(temp, -step_size * weight);
    add(temp, *v_in);
    for (int i = 0; i < m; i++) {
      (*v_out)[i] = temp[i];
    }
  }

  double operator() (double val, int index = 1) {
    return DOUBLE_MARKER;
  }

  void update_step_size(double step_size_) {
    step_size = step_size_;
  }  
  
  void update_cache_vars(double old_x_i, double new_x_i, int index) {
    add(Atx, A, index, -old_x_i + new_x_i);
  }

  void update_cache_vars(Vector* x, int rank, int num_threads){
    int m = At->rows(); //y=A'*x
    int block_size = m/num_threads;
    int start_idx = rank*(block_size);
    int end_idx = (rank == num_threads-1)? m : start_idx+block_size;
    for(int iter=start_idx; iter != end_idx; ++iter){
      (*Atx)[iter]=dot(At, x, iter);
    }
  }
  
  forward_grad_for_square_loss () {
    step_size = 0.;
    weight = 1.;
  }
  forward_grad_for_square_loss (double step_size_, double weight_ = 1.) {
    step_size = step_size_;
    weight = weight_;
  }
  
  forward_grad_for_square_loss (Mat* A_, Vector* b_, Vector* Atx_,
                                double step_size_ = 1., double weight_ = 1.,
                                Mat* At_ = nullptr) {
    step_size = step_size_;
    weight = weight_;
    A = A_;
    b = b_;
    Atx = Atx_;
    At = At_;
  }
};



// forward opertor for quadratic function
// f(x) = weight (0.5 x' Q x + c'x), which is given by
// I - step_size * weight * (Q x + c), where
// Q is a symmetric matrix.
template <typename Mat>
class forward_grad_for_qp : public OperatorInterface {
public:
  Mat* Q;
  Vector *c;
  double weight;
  double step_size;

  double operator() (Vector* x, int index) {
    double Q_ix = dot(Q, x, index);
    double forward_grad_at_i = (*x)[index] - weight * step_size * (Q_ix + (*c)[index]);
    return forward_grad_at_i;
  }
  
  void operator() (Vector* v_in, Vector* v_out) {
    int len = v_in->size();
    multiply(*Q, *v_in, *v_out);
    add(*v_out, *c);
    scale(*v_out, -weight * step_size);
    add(*v_out, *v_in);
  }
  
  double operator() (double val, int index = 1) {
    return DOUBLE_MARKER;
  }

  void update_step_size(double step_size_) {
    step_size = step_size_;
  }  
  void update_cache_vars(double old_x_i, double new_x_i, int index) {
  }
  void update_cache_vars(Vector* x, int rank, int num_threads) {
  }
  
  forward_grad_for_qp () {
    step_size = 0.;
    weight = 1.;
  }
  
  forward_grad_for_qp (double step_size_,double weight_ = 1.) {
    step_size = step_size_;
    weight = weight_;
  }
  
  forward_grad_for_qp (Mat* Q_, Vector* c_, double step_size_ = 1., double weight_=1.) {
    step_size = step_size_;
    weight = weight_;
    Q = Q_;
    c = c_;
  }
};


// forward opertor for quadratic function
// f(x) = (0.5 x'A'Ax - e'x), which is given by
// I - step_size * (A' Ax - e), where
// A is a matrix with size num_features x num_samples
template <typename Mat>
class forward_grad_for_dual_svm : public OperatorInterface {
public:
  Mat* At;
  Vector* Ax;
  Mat* A; // this is for sync-parallel
  double step_size;
  double weight;
  
  double operator() (Vector* x, int index) {
    double At_i_Ax = dot(At, Ax, index);
    double forward_grad_at_i = (*x)[index] - step_size * (At_i_Ax - 1.);
    return forward_grad_at_i;
  }

  double operator() (double val, int index = 1) {
    return DOUBLE_MARKER;
  }
  
  void operator() (Vector* v_in, Vector* v_out) {
    int len = v_in->size();
    multiply(*At, *Ax, *v_out);
    add(*v_out, -1.);
    scale(*v_out, -step_size);
    add(*v_out, *v_in);
  }

  void update_step_size(double step_size_) {
    step_size = step_size_;
  }  
  
  void update_cache_vars(double old_x_i, double new_x_i, int index) {
    add(Ax, At, index, -old_x_i + new_x_i);
  }
  
  void update_cache_vars(Vector* x, int rank, int num_threads) {
    int m = A->rows();
    int block_size = m / num_threads;
    int start_idx = rank * block_size;
    int end_idx = (rank == num_threads - 1) ? m : start_idx + block_size;
    for (int iter = start_idx; iter != end_idx; ++iter){
      (*Ax)[iter] = dot(A, x, iter);
    }
  }
  
  forward_grad_for_dual_svm () {
    step_size = 0.;
    weight = 1.;
  }
  
  forward_grad_for_dual_svm (double step_size_, double weight_ = 1.) {
    step_size = step_size_;
    weight = weight_;
  }
  
  forward_grad_for_dual_svm (Mat* At_, Vector* Ax_, double step_size_ = 1., double weight_ = 1., Mat* A_ = nullptr) {
    step_size = step_size_;
    weight = weight_;
    At = At_;
    Ax = Ax_;
    A = A_;
  }
};

// Jacobi method for linear equations
// A*x = b
// The operator is defined by
// Tx = x - D^{-1} * (A*x - b), where
// D = diag(A)
// A is a diagnally dominant matrix.
template <typename Mat>
class linear_eqn_jacobi_operator : public OperatorInterface {
public:
  Mat* A;
  Vector *b;
  double step_size;
  double weight;
  
  double operator() (Vector* x, int index) {
    double val = 0.;
    double A_ix = dot(A, x, index);
    val = (*x)[index] + ((*b)[index] - A_ix) / (*A)(index, index);
    return val;
  }
  
  void operator() (Vector* v_in, Vector* v_out) {
    int m = A->rows();
    for (int i=0; i<m; i++){
      (*v_out)[i] = (*v_in)[i] + ((*b)[i] - dot(A,v_in,i)) / (*A)(i, i);
    }
  }

  double operator() (double val, int index = 1) {
    return DOUBLE_MARKER;
  }

  void update_step_size(double step_size_) {
    step_size = step_size_;
  }
  void update_cache_vars(double old_x_i, double new_x_i, int index) {
  }
  void update_cache_vars(Vector* x, int rank, int num_threads) {
  }
  
  linear_eqn_jacobi_operator () {
    step_size = 0.;
    weight = 1.;
  }
  
  linear_eqn_jacobi_operator (double step_size_, double weight_ = 1.) {
    step_size = step_size_;
    weight = weight_;
  }
  
  linear_eqn_jacobi_operator (Mat* A_, Vector* b_, double step_size_ = 1., double weight_ = 1.) {
    step_size = step_size_;
    weight = weight_;
    A = A_;
    b = b_;
  }
};


// forward gradient for logistic loss
template <typename Mat>
class forward_grad_for_log_loss : public OperatorInterface {
public:  
  Mat* A;
  Vector *b, *Atx;  
  Mat* At; // this is for sync-parallel
  double weight;
  double step_size;
  
  double operator() (Vector* x, int index) {
    // calculate the forward step
    double A_iAtx = dot(A, Atx, index);
    double A_ib = dot(A, b, index);
    double forward_grad_at_i = (*x)[index] - weight * step_size * log_loss_gradient_at_idx(*A, *b, *Atx, index);
    return forward_grad_at_i;
  }

  // TODO: need to implement this
  void operator()(Vector* v_in, Vector* v_out) {
  }
   
  double operator() (double val, int index = 1) {
    return DOUBLE_MARKER;
  }

  void update_step_size(double step_size_) {
    step_size = step_size_;
  }  
  
  void update_cache_vars(double old_x_i, double new_x_i, int index) {
    add(Atx, A, index, -old_x_i + new_x_i);
  }
  
  void update_cache_vars(Vector* x, int rank, int num_threads) {
    int m = At->rows(); //y=A'*x
    int block_size = m / num_threads;
    int start_idx = rank * block_size;
    int end_idx = (rank == num_threads - 1) ? m : start_idx + block_size;
    for (int iter = start_idx; iter != end_idx; ++iter){
      (*Atx)[iter] = dot(At, x, iter);
    }
  }
  
  forward_grad_for_log_loss() {
    step_size = 0.;
    weight = 1.;
  }
  
  forward_grad_for_log_loss (double step_size_, double weight_=1.) {
    step_size = step_size_;
    weight = weight_;
  }
  
  forward_grad_for_log_loss (Mat* A_, Vector* b_,
                             Vector* Atx_, double step_size_ = 1.,
                             double weight_ = 1., Mat* At_ = nullptr) {
    step_size = step_size_;
    weight = weight_;
    A = A_;
    b = b_;
    Atx = Atx_;
    At = At_;
  }
};


// forward opertor for square hinge loss
// f(x) = weight / 2 sum_i max(0, (1 - b_i a_i'x)^2, which is given by
// x + step_size * weight * sum_i b_i max(0, 1 - b_i a_i'x)) * a_i
template <typename Mat>
class forward_grad_for_square_hinge_loss : public OperatorInterface {
public:
  Mat* A;
  Vector *b, *Atx;  
  Mat* At;
  double weight;
  double step_size;
  
  double operator() (Vector* x, int index) {
    int m = Atx->size();
    // TODO: create a vector in every iteration, need to improve the efficiency
    Vector temp(m, 0);
    for (int i = 0; i < m; ++i) {
      temp[i] = (*b)[i] * fmax(0., 1. - (*b)[i]*(*Atx)[i]);
    }
    double A_itemp = dot(*A, temp, index);
    double forward_grad_at_i = (*x)[index] + weight * step_size * A_itemp;
    return forward_grad_at_i;
  }

  void operator() (Vector* v_in, Vector* v_out) {
    int num_samples = A->cols(), num_features = A->rows();
    Vector temp(num_samples, 0);
    Vector A_temp(num_features, 0.);    
    for (int i = 0; i < num_samples; ++i) {
      temp[i] = (*b)[i] * fmax(0., 1. - (*b)[i]*(*Atx)[i]);      
    }
    multiply(*A, temp, A_temp);
    scale(A_temp, step_size * weight);
    add(A_temp, *v_in);
    for (int i = 0; i < num_features; i++) {
      (*v_out)[i] = A_temp[i];
    }
  }
  
  double operator() (double val, int index = 1) {
    return DOUBLE_MARKER;
  }

  void update_step_size(double step_size_) {
    step_size = step_size_;
  }
  
  void update_cache_vars(double old_x_i, double new_x_i, int index) {
    add(Atx, A, index, -old_x_i + new_x_i);
  }
  
  void update_cache_vars(Vector* x, int rank, int num_threads) {
    int m = At->rows(); //y=A'*x
    int block_size = m / num_threads;
    int start_idx = rank * block_size;
    int end_idx = (rank == num_threads - 1) ? m : start_idx + block_size;
    for (int iter = start_idx; iter != end_idx; ++iter){
      (*Atx)[iter] = dot(At, x, iter);
    }
  }
  
  forward_grad_for_square_hinge_loss () {
    step_size = 0.;
    weight = 1.;
  }
  
  forward_grad_for_square_hinge_loss (double step_size_, double weight_ = 1.) {
    step_size = step_size_;
    weight = weight_;
  }
  
  forward_grad_for_square_hinge_loss (Mat* A_, Vector* b_,
                                      Vector* Atx_, double step_size_ = 1.,
                                      double weight_ = 1., Mat* At_ = nullptr) {
    step_size = step_size_;
    weight = weight_;
    A = A_;
    b = b_;
    Atx = Atx_;
    At = At_;
  }
};


/*********************************************************************
 * forward opertor for Huber loss
 *  f(x) = 0.5 * weight * (a'x - b)^2, if -delta <= a'x-b <= delta;
 *         delta * weight * (|a'x - b| - 0.5 * delta), otherwise.
 *
 * is defined as following:
 *
 *  x^{k+1} = x^k - step_size * weight * A * temp, where temp is 
 *           delta      if a_i' x - b_i > delta
 * temp_i =  a'_ix-b_i  otherwise
             -delta     if a_i' x - b_i < -delta
 *********************************************************************/
template <typename Mat>
class forward_grad_for_huber_loss : public OperatorInterface {
public:
  double delta;
  Mat* A;
  Vector *b, *Atx;  
  Vector temp;
  Mat* At; // this is for sync-parallel
  double weight;
  double step_size;
  
  double operator() (Vector* x, int index) {
    // calculate the forward step
    int num_samples = Atx->size();
    // TODO: create a vector in every iteration, need to improve the efficiency
    // Vector temp(num_samples, 0);
    double Atx_i = 0., b_i = 0;
    for (int i = 0; i < num_samples; ++i) {
      Atx_i = (*Atx)[i], b_i = (*b)[i];
      if (Atx_i - b_i > delta) {
        temp[i] = delta;
      } else if (Atx_i - b_i < -delta) {
        temp[i] = - delta;
      } else {
        temp[i] = Atx_i - b_i;
      }
    }
    double Ai_temp = dot(*A, temp, index);
    double forward_grad_at_i = (*x)[index] - weight * step_size * Ai_temp;
    return forward_grad_at_i;
  }
  
  void operator() (Vector* v_in, Vector* v_out) {
    int num_samples = A->cols(), num_features = A->rows();
    //Vector temp(num_samples, 0);
    Vector A_temp(num_features, 0.);
    double Atx_i = 0., b_i = 0;
    for (int i = 0; i < num_samples; ++i) {
      Atx_i = (*Atx)[i], b_i = (*b)[i];
      if (Atx_i - b_i > delta) {
        temp[i] = delta;
      } else if (Atx_i - b_i < -delta) {
        temp[i] = - delta;
      } else {
        temp[i] = Atx_i - b_i;
      }
    }
    multiply(*A, temp, A_temp);
    scale(A_temp, -step_size * weight);
    add(A_temp, *v_in);
    for (int i = 0; i < num_features; i++) {
      (*v_out)[i] = A_temp[i];
    }
  }

  double operator() (double val, int index = 1) {
    return DOUBLE_MARKER;
  }

  void update_step_size(double step_size_) {
    step_size = step_size_;
  }
  
  void update_cache_vars(double old_x_i, double new_x_i, int index) {
    add(Atx, A, index, -old_x_i + new_x_i);
  }

  void update_cache_vars(Vector* x, int rank, int num_threads) {
    int m = At->rows(); //y=A'*x
    int block_size = m / num_threads;
    int start_idx = rank * block_size;
    int end_idx = (rank == num_threads - 1) ? m : start_idx + block_size;
    for (int iter = start_idx; iter != end_idx; ++iter){
      (*Atx)[iter] = dot(At, x, iter);
    }
  }
  
  forward_grad_for_huber_loss () {
    step_size = 0.;
    weight = 1.;
  }
  forward_grad_for_huber_loss (Mat* A_, Vector* b_, Vector* Atx_,
                               double step_size_ = 1., double weight_=1.,
                               double delta_ = 1., Mat* At_ = nullptr) {
    step_size = step_size_;
    weight = weight_;
    A = A_;
    b = b_;
    Atx = Atx_;
    temp = Vector(A->cols(),0.);
    delta = delta_;
    At = At_;
  }
};


// 3S operator for portfolio optimization
// please refer to page 38 of the following paper
//   http://arxiv.org/pdf/1601.00863v2.pdf
template <typename Mat>
class portfolio_3s : public OperatorInterface {
public:
  Mat* Q;
  Vector* epsilon;
  double c;
  double b1, b2, b3, b4;
  double a1a2;
  Vector a1, a2, a3, a4, tild_a1, tild_a2;
  double w1, w2, w3, w4;
  double weight;
  double step_size;

  double operator() (Vector* x, int index) {
    w1 = dot(a1, *x) - b1;
    w2 = dot(a2, *x) - b2;
    w3 = dot(a3, *x) - b3;
    w4 = dot(a4, *x) - b4;    
    double x_i = 0.;
    if (w1 <= 0 && w2 >= 0) {
      x_i = max(0., (*x)[index] - step_size * dot(Q, x, index));
    } else if (w2 <= 0 && w3 >= 0) {
      x_i = w2 * a2[index] + max(0., (*x)[index]
                               - step_size * dot(Q, x, index)
                               - w2 * (2. * a2[index] - step_size * dot(Q, &a2, index)));
    } else if (w3 <= 0 && w4 >= 0) {
      x_i = w1 * tild_a1[index] + w2 * tild_a2[index] + max(0., (*x)[index]
                                                        - step_size * dot(Q, x, index)
                                                        - w1 * (2. * tild_a1[index] - step_size * dot(Q, &tild_a1, index))
                                                        - w2 * (2. * tild_a2[index] - step_size * dot(Q, &tild_a2, index)));
    } else {
      x_i = w1 * a1[index] + max(0., (*x)[index]
                               - step_size * dot(Q, x, index)
                               - w1 * (2. * a1[index] - step_size * dot(Q, &a1, index)));
    }
    return x_i;
  }

  // TODO: implement the full update operator
  void operator() (Vector* v_in, Vector* v_out) {
  }

  double operator() (double val, int index = 1) {
    return DOUBLE_MARKER;
  }

  void update_step_size(double step_size_) {
    step_size = step_size_;
  }
  void update_cache_vars(double old_x_i, double new_x_i, int index) {
  }
  void update_cache_vars(Vector* x, int rank, int num_threads) {
  }
  
  void project_D2 (Vector* v_in, Vector* v_out) {
    w1 = dot(a1, *v_in) - b1;
    w2 = dot(a2, *v_in) - b2;
    w3 = dot(a3, *v_in) - b3;
    w4 = dot(a4, *v_in) - b4;
    size_t len = v_in->size();
    Vector temp(len, 0.);
    temp = (*v_in);
    if (w1 <= 0 && w2 >= 0) {
      temp = (*v_in);
    } else if (w2 <= 0 && w3 >= 0) {
      add(temp, a2, -w2);
    } else if (w3 <= 0 && w4 >= 0) {
      add(temp, tild_a1, -w1);
      add(temp, tild_a2, -w2);
    } else {
      add(temp, a1, -w1);
    }
    for (size_t i = 0; i < len; i++) {
      (*v_out)[i] = temp[i];
    }
  }
  portfolio_3s () {
    step_size = 0.;
    weight = 1.;
    int m = 1;    
    Vector epsilon_(m, 1.);
    double sm = 1.;
    double nrm_epsilon = norm(epsilon_);
    double c_ = 1.;    
    Q = nullptr;
    epsilon = nullptr;
    c = 0.;
    a1 = Vector(m, 1./sm);
    a2 = Vector(m, 0.);
    a3 = Vector(m, 1.);
    a4 = Vector(m, 1.);
    a1a2 = dot(a1, a2);
    tild_a1 = Vector(m, 1.);
    tild_a2 = Vector(m, 1.);    
    b1 = 1./sm;
    b2 =  c_ / nrm_epsilon;
    b3 = b2 - b1 / a1a2;
    b4 = b1 - b2 / a1a2;
    w1 = -b1;
    w2 = -b2;
    w3 = -b3;
    w4 = -b4;
  }
  
  portfolio_3s (double step_size_, double weight_ = 1.) {
    step_size = step_size_;
    weight = weight_;
    int m = 1;    
    Vector epsilon_(m, 1.);
    double sm = 1.;
    double nrm_epsilon = norm(epsilon_);
    double c_ = 1.;
    a1 = Vector(m, 1./sm);
    a2 = Vector(m, 0.);
    a3 = Vector(m, 1.);
    a4 = Vector(m, 1.);
    a1a2 = dot(a1, a2);
    tild_a1 = Vector(m, 1.);
    tild_a2 = Vector(m, 1.);    
    Q = nullptr;
    epsilon = nullptr;
    b1 = 1./sm;
    b2 =  c_ / nrm_epsilon;
    b3 = b2 - b1 / a1a2;
    b4 = b1 - b2 / a1a2;
    w1 = -b1;
    w2 = -b2;
    w3 = -b3;
    w4 = -b4;
  }
  
  portfolio_3s (Mat* Q_,
                Vector* epsilon_,
                double c_,
                double step_size_ = 1.,
                double weight_ = 1.) {
    step_size = step_size_;
    weight = weight_;
    epsilon = epsilon_;
    c = c_;
    Q = Q_;
    int m = Q->rows();
    double sm = sqrt((double)(m));
    double nrm_epsilon = norm(*epsilon_);
    // initializing a1, a2, a3, a4
    a1 = Vector(m, 1./sm);
    a2 = Vector(m, 0.);
    a2 = (*epsilon);
    scale(a2, 1./nrm_epsilon);
    a1a2 = dot(a1, a2);
    a3 = Vector(m, 0.);
    a4 = Vector(m, 0.);
    tild_a1 = Vector(m, 0.);
    tild_a2 = Vector(m, 0.);
    a3 = a2;
    add(a3, a1, -1./a1a2);
    a4 = a1;
    add(a4, a2, -1./a1a2);
    tild_a1 = a1;
    add(tild_a1, a2, -a1a2);
    scale(tild_a1, 1./(1 - a1a2 * a1a2));
    tild_a2 = a2;
    add(tild_a2, a1, -a1a2);
    scale(tild_a2, 1./(1 - a1a2 * a1a2));
    // setting up b1, b2, b3, b4
    b1 = 1./sm;
    b2 = c_ / nrm_epsilon;
    b3 = b2 - b1 / a1a2;
    b4 = b1 - b2 / a1a2;
  }
};

#endif
