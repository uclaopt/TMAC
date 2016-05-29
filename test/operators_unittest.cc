#include "gtest/gtest.h"
#include "matrices.h"
#include "algebra.h"
#include "operators.h"
#include "constants.h"
#include "algebra_namespace_switcher.h"
// used for DEBUG
#include<iostream>





/**************************************
   Unit test for proximal operators
***************************************/


// test the constructors 
TEST(ProximalL1Test, Constructor) {
  // default constructor
  prox_l1 p_l1;
  EXPECT_DOUBLE_EQ(p_l1.step_size, 0.);
  EXPECT_DOUBLE_EQ(p_l1.weight, 1.);

  // constructor initialized by step_size
  double step_size = 10.;
  prox_l1 p_Sp(step_size);
  EXPECT_DOUBLE_EQ(p_Sp.step_size, 10.);
  EXPECT_DOUBLE_EQ(p_Sp.weight, 1.);

  // constructor initialized by step_size and weight
  double weight = 9.99;
  prox_l1 p3(step_size, weight);
  EXPECT_DOUBLE_EQ(p3.step_size, 10.);
  EXPECT_DOUBLE_EQ(p3.weight, weight);
  
}

// test the update_step_size member function
TEST(ProximalL1Test, UpdateStepSize) {
  double step_size = 1.;
  prox_l1 p_l1(step_size);
  p_l1.update_step_size(100.);
  EXPECT_DOUBLE_EQ(p_l1.step_size, 100.);
  EXPECT_DOUBLE_EQ(p_l1.weight, 1.);  
}

// test the scalar operator
TEST(ProximalL1Test, ScalarOperator) {

  double step_size = 0.1;
  prox_l1 p_l1(step_size);
  // test numbers bigger than step_size
  EXPECT_DOUBLE_EQ(p_l1(1.), 0.9);
  // test numbers in [-step_size, step_size]
  EXPECT_DOUBLE_EQ(p_l1(0.01), 0.);
  // test numbers smaller than -step_size
  EXPECT_DOUBLE_EQ(p_l1(-1.), -0.9);

}

// test the vector operator 
TEST(ProximalL1Test, VectorOperator) {
  
  int n = 2;
  Vector v(n, 2.);
  double step_size = 1.;
  prox_l1 p_l1(step_size);
  EXPECT_DOUBLE_EQ(p_l1(&v, 0), 1.);
  EXPECT_DOUBLE_EQ(v[0], 2.);
  
  p_l1.update_step_size(2.01);
  EXPECT_DOUBLE_EQ(p_l1.step_size, 2.01);    
  EXPECT_DOUBLE_EQ(p_l1(&v, 1), 0.);  
  
}


// test the full update operator 
TEST(ProximalL1Test, FullOperator) {
  
  int n = 1;
  Vector v_in(n, 2.), v_out(n, 0.);
  double step_size = 1.;
  prox_l1 p_l1(step_size);
  p_l1(&v_in, &v_out);
  EXPECT_DOUBLE_EQ(v_out[0], 1.);
  
  p_l1.update_step_size(2.01);
  EXPECT_DOUBLE_EQ(p_l1.step_size, 2.01);    
  p_l1(&v_in, &v_out);
  EXPECT_DOUBLE_EQ(v_out[0], 0.);
  
}


// unit test for the prox_sum_square operator
TEST(ProximalSumSquareTest, Constructor) {
  // default constructor
  prox_sum_square p1;
  EXPECT_DOUBLE_EQ(p1.step_size, 0.);
  EXPECT_DOUBLE_EQ(p1.weight, 1.);

  // constructor initialized by step_size
  double step_size = 10.;
  prox_sum_square p2(step_size);
  EXPECT_DOUBLE_EQ(p2.step_size, 10.);
  EXPECT_DOUBLE_EQ(p2.weight, 1.);

  // constructor initialized by step_size and weight
  double weight = 9.99;
  prox_sum_square p3(step_size, weight);
  EXPECT_DOUBLE_EQ(p3.step_size, 10.);
  EXPECT_DOUBLE_EQ(p3.weight, weight);
  
}

// test the update_step_size member function
TEST(ProximalSumSquareTest, UpdateStepSize) {
  double step_size = 1.;
  prox_sum_square p(step_size);
  p.update_step_size(100.);
  EXPECT_DOUBLE_EQ(p.step_size, 100.);
  EXPECT_DOUBLE_EQ(p.weight, 1.);  
}

// test the scalar operator
TEST(ProximalSumSquareTest, ScalarOperator) {

  double step_size = 1.;
  prox_sum_square p(step_size);
  // test numbers bigger than step_size
  EXPECT_DOUBLE_EQ(p(1.), 0.5);
  // test numbers in [-step_size, step_size]
  EXPECT_DOUBLE_EQ(p(-1.), -0.5);
  // test numbers smaller than -step_size
  EXPECT_DOUBLE_EQ(p(0.), 0.);
  // test with a given index
  EXPECT_DOUBLE_EQ(p(0., 100), 0.);  
}

// test the vector operator 
TEST(ProximalSumSquareTest, VectorOperator) {
  
  int n = 2;
  Vector v(n, 2.);
  double step_size = 1.;
  prox_sum_square p(step_size);
  EXPECT_DOUBLE_EQ(p(&v, 0), 1.);
  EXPECT_DOUBLE_EQ(v[0], 2.);
  
  p.update_step_size(3.);
  EXPECT_DOUBLE_EQ(p.step_size, 3.);    
  EXPECT_DOUBLE_EQ(p(&v, 1.), 0.5);
  
}


// test the full update operator 
TEST(ProximalSumSquareTest, FullOperator) {
  
  int n = 1;
  Vector v_in(n, 2.), v_out(n, 0.);
  double step_size = 1.;
  prox_sum_square p(step_size);
  p(&v_in, &v_out);
  EXPECT_DOUBLE_EQ(v_out[0], 1.);
  
  p.update_step_size(3.);
  EXPECT_DOUBLE_EQ(p.step_size, 3.);    
  p(&v_in, &v_out);
  EXPECT_DOUBLE_EQ(v_out[0], 0.5);
  
}


// unit test for the proximal l2 norm
TEST(ProximalL2Test, Constructor) {
  // default constructor
  prox_l2 p1;
  EXPECT_DOUBLE_EQ(p1.step_size, 0.);
  EXPECT_DOUBLE_EQ(p1.weight, 1.);

  // constructor initialized by step_size
  double step_size = 10.;
  prox_l2 p2(step_size);
  EXPECT_DOUBLE_EQ(p2.step_size, 10.);
  EXPECT_DOUBLE_EQ(p2.weight, 1.);

  // constructor initialized by step_size and weight
  double weight = 9.99;
  prox_l2 p3(step_size, weight);
  EXPECT_DOUBLE_EQ(p3.step_size, 10.);
  EXPECT_DOUBLE_EQ(p3.weight, weight);
  
}

// test the update_step_size member function
TEST(ProximalL2Test, UpdateStepSize) {
  double step_size = 1.;
  prox_l2 p(step_size);
  p.update_step_size(100.);
  EXPECT_DOUBLE_EQ(p.step_size, 100.);
  EXPECT_DOUBLE_EQ(p.weight, 1.);  
}

// test the scalar operator
TEST(ProximalL2Test, ScalarOperator) {

  double step_size = 1.;
  prox_l2 p(step_size);
  // test numbers bigger than step_size
  EXPECT_DOUBLE_EQ(p(1.), DOUBLE_MARKER);
  // test numbers in [-step_size, step_size]
  EXPECT_DOUBLE_EQ(p(-1.), DOUBLE_MARKER);
  // test numbers smaller than -step_size
  EXPECT_DOUBLE_EQ(p(0.), DOUBLE_MARKER);
  // test with a given index
  EXPECT_DOUBLE_EQ(p(0., 100), DOUBLE_MARKER);  
}

// test the vector operator 
TEST(ProximalL2Test, VectorOperator) {
  
  int n = 2;
  Vector v(n, 2.);
  double step_size = 10.;
  prox_l2 p(step_size);
  EXPECT_DOUBLE_EQ(p(&v, 0), 0.);
  EXPECT_DOUBLE_EQ(v[0], 2.);
  
  p.update_step_size(1.);
  EXPECT_DOUBLE_EQ(p.step_size, 1.);
  EXPECT_DOUBLE_EQ(p(&v, 1), 2. - sqrt(2) / 2.);
  
}


// test the full update operator 
TEST(ProximalL2Test, FullOperator) {
  
  int n = 2;
  Vector v_in(n, 2.), v_out(n, 1.);
   double step_size = 10.;
  prox_l2 p(step_size);
  p(&v_in, &v_out);
  EXPECT_DOUBLE_EQ(v_out[0], 0.);
  
  p.update_step_size(1.);
  EXPECT_DOUBLE_EQ(p.step_size, 1.);    
  p(&v_in, &v_out);
  EXPECT_DOUBLE_EQ(v_out[0], 2. - sqrt(2) / 2.);
  
}


// test the constructors 
TEST(ProximalLogBarrierTest, Constructor) {
  // default constructor
  prox_log_barrier p_lb;
  EXPECT_DOUBLE_EQ(p_lb.step_size, 0.);
  EXPECT_DOUBLE_EQ(p_lb.weight, 1.);

  // constructor initialized by step_size
  double step_size = 10.;
  prox_log_barrier p2(step_size);
  EXPECT_DOUBLE_EQ(p2.step_size, 10.);
  EXPECT_DOUBLE_EQ(p2.weight, 1.);

  // constructor initialized by step_size and weight
  double weight = 9.99;
  prox_log_barrier p3(step_size, weight);
  EXPECT_DOUBLE_EQ(p3.step_size, 10.);
  EXPECT_DOUBLE_EQ(p3.weight, weight);
  
}

// test the update_step_size member function
TEST(ProximalLogBarrierTest, UpdateStepSize) {
  double step_size = 1.;
  prox_log_barrier p_lb(step_size);
  p_lb.update_step_size(100.);
  EXPECT_DOUBLE_EQ(p_lb.step_size, 100.);
  EXPECT_DOUBLE_EQ(p_lb.weight, 1.);  
}

// test the scalar operator
TEST(ProximalLogBarrierTest, ScalarOperator) {

  double step_size = 4.;
  prox_log_barrier p_lb(step_size);
  EXPECT_DOUBLE_EQ(p_lb(3.), 4.);
  EXPECT_DOUBLE_EQ(p_lb(0.0), 2.);
}

// test the vector operator 
TEST(ProximalLogBarrierTest, VectorOperator) {
  
  int n = 2;
  Vector v(n, 3.);
  double step_size = 4.;
  prox_log_barrier p_lb(step_size);
  EXPECT_DOUBLE_EQ(p_lb(&v, 0), 4.);
  EXPECT_DOUBLE_EQ(v[0], 3.);
  
}


// test the full update operator 
TEST(ProximalLogBarrierTest, FullOperator) {
  
  int n = 1;
  Vector v_in(n, 3.), v_out(n, 0.);
  double step_size = 4.;
  prox_log_barrier p_lb(step_size);
  p_lb(&v_in, &v_out);
  EXPECT_DOUBLE_EQ(v_out[0], 4.);

}



// test the constructors 
TEST(ProximalElasticNetTest, Constructor) {

  // default constructor
  prox_elastic_net p_en;
  EXPECT_DOUBLE_EQ(p_en.step_size, 0.);
  EXPECT_DOUBLE_EQ(p_en.weight, 1.);
  EXPECT_DOUBLE_EQ(p_en.weight_2, 1.);

  // constructor initialized by step_size
  double step_size = 10.;
  prox_elastic_net p2(step_size);
  EXPECT_DOUBLE_EQ(p2.step_size, 10.);


  // constructor initialized by step_size and weights
  double weight = 9.99, weight_2 = 9.99;
  prox_elastic_net p3(step_size, weight, weight_2);
  EXPECT_DOUBLE_EQ(p3.step_size, step_size);
  EXPECT_DOUBLE_EQ(p3.weight, weight);
  EXPECT_DOUBLE_EQ(p3.weight_2, weight_2);  
  
}

// test the update_step_size member function
TEST(ProximalElasticNetTest, UpdateStepSize) {
  double step_size = 1.;
  prox_elastic_net p_en(step_size);
  p_en.update_step_size(100.);
  EXPECT_DOUBLE_EQ(p_en.step_size, 100.);
}

// test the scalar operator
TEST(ProximalElasticNetTest, ScalarOperator) {

  double step_size = 1.;
  prox_elastic_net p_en(step_size);
  // test numbers bigger than step_size
  EXPECT_DOUBLE_EQ(p_en(2.), 0.5);
  // test numbers in [-step_size, step_size]
  EXPECT_DOUBLE_EQ(p_en(0.01), 0.);
  // test numbers smaller than -step_size
  EXPECT_DOUBLE_EQ(p_en(-2.), -0.5);

}

// test the vector operator 
TEST(ProximalElasticNetTest, VectorOperator) {
  
  int n = 2;
  Vector v(n, 2.);
  double step_size = 1.;
  prox_elastic_net p_en(step_size);
  EXPECT_DOUBLE_EQ(p_en(&v, 0), .5);
  EXPECT_DOUBLE_EQ(v[0], 2.);
  
  p_en.update_step_size(10.);
  EXPECT_DOUBLE_EQ(p_en.step_size, 10.);    
  EXPECT_DOUBLE_EQ(p_en(&v, 1), 0.);  
  
}


// test the full update operator 
TEST(ProximalElasticNetTest, FullOperator) {
  
  int n = 1;
  Vector v_in(n, 2.), v_out(n, 0.);
  double step_size = 1.;
  prox_elastic_net p_en(step_size);
  p_en(&v_in, &v_out);
  EXPECT_DOUBLE_EQ(v_out[0], .5);
  
  p_en.update_step_size(10.);
  EXPECT_DOUBLE_EQ(p_en.step_size, 10.);    
  p_en(&v_in, &v_out);
  EXPECT_DOUBLE_EQ(v_out[0], 0.);
  
}



// test the constructors 
TEST(ProximalHuberTest, Constructor) {

  // default constructor
  prox_huber p_h;
  EXPECT_DOUBLE_EQ(p_h.step_size, 0.);
  EXPECT_DOUBLE_EQ(p_h.weight, 1.);
  EXPECT_DOUBLE_EQ(p_h.delta, 1.);

  // constructor initialized by step_size
  double step_size = 10.;
  prox_huber p2(step_size);
  EXPECT_DOUBLE_EQ(p2.step_size, 10.);


  // constructor initialized by step_size and weights
  double weight = 9.99, delta = 9.99;
  prox_huber p3(step_size, weight, delta);
  EXPECT_DOUBLE_EQ(p3.step_size, step_size);
  EXPECT_DOUBLE_EQ(p3.weight, weight);
  EXPECT_DOUBLE_EQ(p3.delta, delta);  
  
}

// test the update_step_size member function
TEST(ProximalHuberTest, UpdateStepSize) {
  double step_size = 1.;
  prox_huber p_h(step_size);
  p_h.update_step_size(100.);
  EXPECT_DOUBLE_EQ(p_h.step_size, 100.);
}

// test the scalar operator
TEST(ProximalHuberTest, ScalarOperator) {

  double step_size = 1.;
  prox_huber p_h(step_size);
  // test numbers bigger than step_size
  EXPECT_DOUBLE_EQ(p_h(2.), 1.);
  // test numbers in [-step_size, step_size]
  EXPECT_DOUBLE_EQ(p_h(0.02), 0.01);
  // test numbers smaller than -step_size
  EXPECT_DOUBLE_EQ(p_h(-2.), -1.);

}


// test the vector operator 
TEST(ProximalHuberTest, VectorOperator) {
  
  int n = 2;
  Vector v(n, 2.);
  double step_size = 1.;
  prox_huber p_h(step_size);
  EXPECT_DOUBLE_EQ(p_h(&v, 0), 1.);
  EXPECT_DOUBLE_EQ(v[0], 2.);
  
  p_h.update_step_size(9.);
  EXPECT_DOUBLE_EQ(p_h.step_size, 9.);    
  EXPECT_DOUBLE_EQ(p_h(&v, 1), 0.2);  
  
}


// test the full update operator 
TEST(ProximalHuberTest, FullOperator) {
  
  int n = 1;
  Vector v_in(n, 2.), v_out(n, 0.);
  double step_size = 1.;
  prox_huber p_h(step_size);
  p_h(&v_in, &v_out);
  EXPECT_DOUBLE_EQ(v_out[0], 1.);
 
  p_h.update_step_size(9.);
  EXPECT_DOUBLE_EQ(p_h.step_size, 9.);
  p_h(&v_in, &v_out);
  EXPECT_DOUBLE_EQ(v_out[0], 0.2);
  
}



/*****************************************
 *   Unit test for projection operators
 *****************************************/

// unit test for projection to the positive cone operator
TEST(ProjectionPositiveConeTest, Constructor) {
  // default constructor
  proj_positive_cone p1;
  EXPECT_DOUBLE_EQ(p1.step_size, 0.);
  EXPECT_DOUBLE_EQ(p1.weight, 1.);

  // constructor initialized by step_size
  double step_size = 10.;
  proj_positive_cone p2(step_size);
  EXPECT_DOUBLE_EQ(p2.step_size, 10.);
  EXPECT_DOUBLE_EQ(p2.weight, 1.);

  // constructor initialized by step_size and weight
  double weight = 9.99;
  proj_positive_cone p3(step_size, weight);
  EXPECT_DOUBLE_EQ(p3.step_size, 10.);
  EXPECT_DOUBLE_EQ(p3.weight, weight);
  
}

// test the update_step_size member function
TEST(ProjectionPositiveConeTest, UpdateStepSize) {
  double step_size = 1.;
  proj_positive_cone p(step_size);
  p.update_step_size(100.);
  EXPECT_DOUBLE_EQ(p.step_size, 100.);
  EXPECT_DOUBLE_EQ(p.weight, 1.);  
}

// test the scalar operator
TEST(ProjectionPositiveConeTest, ScalarOperator) {

  double step_size = 1.;
  proj_positive_cone p(step_size);
  // test numbers bigger than step_size
  EXPECT_DOUBLE_EQ(p(1.), 1.);
  // test numbers in [-step_size, step_size]
  EXPECT_DOUBLE_EQ(p(-1.), 0.);
  // test numbers smaller than -step_size
  EXPECT_DOUBLE_EQ(p(0.), 0.);
  // test with a given index
  EXPECT_DOUBLE_EQ(p(0., 100), 0.);  
}

// test the vector operator 
TEST(ProjectionPositiveConeTest, VectorOperator) {
  
  int n = 2;
  Vector v(n, 2.);
  v[1] = -1;
  double step_size = 10.;
  proj_positive_cone p(step_size);
  EXPECT_DOUBLE_EQ(p(&v, 0), 2.);
  EXPECT_DOUBLE_EQ(v[0], 2.);
  
  p.update_step_size(1.);
  EXPECT_DOUBLE_EQ(p.step_size, 1.);
  EXPECT_DOUBLE_EQ(p(&v, 1), 0.);
  
}


// test the full update operator 
TEST(ProjectionPositiveConeTest, FullOperator) {
  
  int n = 2;
  Vector v_in(n, 2.), v_out(n, 1.);
  v_in[1] = -2.;
  double step_size = 10.;
  proj_positive_cone p(step_size);
  p(&v_in, &v_out);
  EXPECT_DOUBLE_EQ(v_out[0], 2.);
  EXPECT_DOUBLE_EQ(v_out[1], 0.);  

  p.update_step_size(1.);
  EXPECT_DOUBLE_EQ(p.step_size, 1.);    
  p(&v_in, &v_out);
  EXPECT_DOUBLE_EQ(v_out[0], 2.);
  EXPECT_DOUBLE_EQ(v_out[1], 0.);  
  
}


// test the vector operator 
TEST(ProjectionBoxTest, VectorOperator) {
  
  int n = 3;
  Vector lower(n, 0.), upper(n, 1.5);
  Vector v(n, 2.);
  v[1] = -1;
  v[2] = 1.;
  double step_size = 10.;
  proj_box p(&lower, &upper);
  EXPECT_DOUBLE_EQ(p(&v, 0), 1.5);
  EXPECT_DOUBLE_EQ(p(&v, 1), 0.);
  EXPECT_DOUBLE_EQ(p(&v, 2), 1.);  

}



// unit test for projection to the positive cone operator
TEST(ProjectionHyperplaneTest, Constructor) {
  // default constructor
  proj_hyperplane p1;
  EXPECT_DOUBLE_EQ(p1.step_size, 0.);
  EXPECT_DOUBLE_EQ(p1.weight, 1.);

  // customized contructor
  double step_size = 1., weight = 9.;  
  Vector a(2, 1.);
  double b = 1.;
  proj_hyperplane p3(&a, b, step_size, weight);
  int len = a.size();
  for(int i=0; i < len; i++) {
    EXPECT_DOUBLE_EQ((*p3.a)[i], a[i]); 
  }
  EXPECT_DOUBLE_EQ(p3.b, b);
}

// test the update_step_size member function
TEST(ProjectionHyperplaneTest, UpdateStepSize) {
  double step_size = 1., weight = 9.;  
  Vector a(2, 1.);
  double b = 1.;
  proj_hyperplane p(&a, b, step_size, weight);

  p.update_step_size(100.);
  EXPECT_DOUBLE_EQ(p.step_size, 100.);
  EXPECT_DOUBLE_EQ(p.weight, weight);
}

// test the scalar operator
TEST(ProjectionHyperplaneTest, ScalarOperator) {

  double step_size = 0.;
  double weight = 1.;
  Vector a(2, 1.);
  a[0] = 0;
  double b = 1.;
  proj_hyperplane p(&a, b, step_size, weight);

  EXPECT_DOUBLE_EQ(p(1.), DOUBLE_MARKER);

}

// test the vector operator 
TEST(ProjectionHyperplaneTest, VectorOperator) {

  Vector a(2, 1.);
  a[0] = 0;
  double b = 1.;
  proj_hyperplane p(&a, b);
  Vector x(2, 0);
  x[1] = 2;
  
  EXPECT_DOUBLE_EQ(p(&x, 0), 0.);
  EXPECT_DOUBLE_EQ(p(&x, 1), 1.);
  
}


// test the full update operator 
TEST(ProjectionHyperplaneTest, FullOperator) {
  
  int n = 2;
  Vector v_in(n, 2.), v_out(n, 1.);
  
}

// unit test for projection to probability simplex
TEST(ProjectionProbSimplexTest, Constructor) {
  // default constructor
  proj_prob_simplex p1;
  EXPECT_DOUBLE_EQ(p1.step_size, 0.);
  EXPECT_DOUBLE_EQ(p1.weight, 1.);

  // customized contructor
  double step_size = 2.;
  proj_prob_simplex p2(step_size);
  EXPECT_DOUBLE_EQ(p2.step_size, step_size);
}

// test the scalar operator
TEST(ProjectionProbSimplexTest, ScalarOperator) {

  double step_size = 4.; 
  proj_prob_simplex p(step_size);
  EXPECT_DOUBLE_EQ(p(1.), DOUBLE_MARKER);

}

// test the vector operator 
TEST(ProjectionProbSimplexTest, VectorOperator) {

  proj_prob_simplex p;
  Vector v(2, 1.);
  
  EXPECT_DOUBLE_EQ(p(&v, 0), .5);
  EXPECT_DOUBLE_EQ(p(&v, 1), .5);

  v[0] = -1.;
  v[1] = 1.;
  EXPECT_DOUBLE_EQ(p(&v, 0), 0.);
  EXPECT_DOUBLE_EQ(p(&v, 1), 1.);
  
  v[0] = -1.;
  v[1] = -1.;
  EXPECT_DOUBLE_EQ(p(&v, 0), .5);
  EXPECT_DOUBLE_EQ(p(&v, 1), .5);
  
  v[0] = 1.;
  v[1] = -1.;
  EXPECT_DOUBLE_EQ(p(&v, 0), 1.);
  EXPECT_DOUBLE_EQ(p(&v, 1), 0.);

}

// test the full update operator 
TEST(ProjectionProbSimplexTest, FullOperator) {

  int n = 2;
  Vector v_in(n, 1.), v_out(n, 0.);
  
  double radius = 4.;
  proj_prob_simplex p(radius);
  
  p(&v_in, &v_out);

  EXPECT_DOUBLE_EQ(v_out[0], .5);
  EXPECT_DOUBLE_EQ(v_out[1], .5);

  v_in[0] = -1.;
  v_in[1] = 1.;
  p(&v_in, &v_out);  
  EXPECT_DOUBLE_EQ(v_out[0], 0.);
  EXPECT_DOUBLE_EQ(v_out[1], 1.);

  v_in[0] = -1.;
  v_in[1] = -1.;
  p(&v_in, &v_out);  
  EXPECT_DOUBLE_EQ(v_out[0], .5);
  EXPECT_DOUBLE_EQ(v_out[1], .5);  
 
  v_in[0] = 1.;
  v_in[1] = -1.;
  p(&v_in, &v_out);  
  EXPECT_DOUBLE_EQ(v_out[0], 1.);
  EXPECT_DOUBLE_EQ(v_out[1], 0.);

}

// test update_stepSize
TEST(ProjectionProbSimplexTest, UpdateStepSize) {
  double radius = 4.;
  proj_l1_ball p(radius);
  p.update_step_size(100.);
  EXPECT_DOUBLE_EQ(p.step_size, 100.);
}

// unit test for projection to the l1 ball
TEST(ProjectionL1BallTest, Constructor) {
  // default constructor
  proj_l1_ball p1;
  EXPECT_DOUBLE_EQ(p1.step_size, 0.);
  EXPECT_DOUBLE_EQ(p1.weight, 1.);

  // customized contructor
  double radius = 2.;
  proj_l1_ball p2(radius);
  EXPECT_DOUBLE_EQ(p2.radius, radius);
}

// test the scalar operator
TEST(ProjectionL1BallTest, ScalarOperator) {

  double radius = 4.; 
  proj_l1_ball p(radius);
  EXPECT_DOUBLE_EQ(p(1.), DOUBLE_MARKER);

}

// test the vector operator 
TEST(ProjectionL1BallTest, VectorOperator) {

  double radius = 4.;
  proj_l1_ball p(radius);
  Vector v(2, 1.);
  

  EXPECT_DOUBLE_EQ(p(&v, 0), 1.);
  EXPECT_DOUBLE_EQ(p(&v, 1), 1.);

  p.update_radius(1.);
  EXPECT_DOUBLE_EQ(p(&v, 0), .5);
  EXPECT_DOUBLE_EQ(p(&v, 1), .5);


  p.update_radius(4.);
  v[0] = -2;
  v[1] = 3;

  EXPECT_DOUBLE_EQ(p(&v, 0), -1.5);
  EXPECT_DOUBLE_EQ(p(&v, 1), 2.5);
  
}

// test the full update operator 
TEST(ProjectionL1BallTest, FullOperator) {

  int n = 2;
  Vector v_in(n, 1.), v_out(n, 0.);
  
  double radius = 4.;
  proj_l1_ball p(radius);
  
  p(&v_in, &v_out);

  EXPECT_DOUBLE_EQ(v_out[0], 1.);
  EXPECT_DOUBLE_EQ(v_out[1], 1.);

  p.update_radius(1.);

  p(&v_in, &v_out);  
  EXPECT_DOUBLE_EQ(v_out[0], 0.5);
  EXPECT_DOUBLE_EQ(v_out[1], 0.5);

  
  p.update_radius(4.);
  v_in[0] = -2;
  v_in[1] = 3;
  p(&v_in, &v_out);  
  EXPECT_DOUBLE_EQ(v_out[0], -1.5);
  EXPECT_DOUBLE_EQ(v_out[1], 2.5);

}

// test update_stepSize
TEST(ProjectionL1BallTest, UpdateStepSize) {
  double radius = 4.;
  proj_l1_ball p(radius);
  p.update_step_size(100.);
  EXPECT_DOUBLE_EQ(p.step_size, 100.);
}

// unit test for projection to the l2 ball
TEST(ProjectionL2BallTest, Constructor) {
  // default constructor
  proj_l2_ball p1;
  EXPECT_DOUBLE_EQ(p1.step_size, 0.);
  EXPECT_DOUBLE_EQ(p1.weight, 1.);

  // customized contructor
  double radius = 2.;
  proj_l2_ball p2(radius);
  EXPECT_DOUBLE_EQ(p2.radius, radius);
}

// test the update_step_size member function
TEST(ProjectionL2BallTest, UpdateStepSize) {
  double radius = 4.;
  proj_l2_ball p(radius);
  p.update_step_size(100.);
  EXPECT_DOUBLE_EQ(p.step_size, 100.);
}

// test the scalar operator
TEST(ProjectionL2BallTest, ScalarOperator) {

  double radius = 4.; 
  proj_l2_ball p(radius);
  EXPECT_DOUBLE_EQ(p(1.), DOUBLE_MARKER);

}

// test the vector operator 
TEST(ProjectionL2BallTest, VectorOperator) {
  
  double radius = 4.;
  proj_l2_ball p(radius);
  Vector v(2, 1.);
  
  EXPECT_DOUBLE_EQ(p(&v, 0), 1.);
  p.update_radius(1.);
  EXPECT_DOUBLE_EQ(p(&v, 0), 1./sqrt(2));
  
}

// test the full update operator 
TEST(ProjectionL2BallTest, FullOperator) {
  
  int n = 2;
  Vector v_in(n, 1.), v_out(n, 0.);
  double radius = 4.;
  proj_l2_ball p(radius);
  p(&v_in, &v_out);
  EXPECT_DOUBLE_EQ(v_out[0], 1.);
  p.update_radius(1.);
  p(&v_in, &v_out);  
  EXPECT_DOUBLE_EQ(v_out[1], 1./sqrt(2));
}




/**************************************
   Unit test for forward operators
***************************************/


// unit test for forward operator for quadratic function
TEST(ForwardQPTest, Constructor) {
  //test for class Matrix
  // default constructor
  forward_grad_for_qp<Matrix> p1;
  EXPECT_DOUBLE_EQ(p1.step_size, 0.);
  EXPECT_DOUBLE_EQ(p1.weight, 1.);

  // constructor initialized by step_size
  double step_size = 10.;
  forward_grad_for_qp<Matrix> p2(step_size);
  EXPECT_DOUBLE_EQ(p2.step_size, 10.);
  EXPECT_DOUBLE_EQ(p2.weight, 1.);

  // constructor initialized by step_size and weight
  double weight = 9.99;
  forward_grad_for_qp<Matrix> p3(step_size, weight);
  EXPECT_DOUBLE_EQ(p3.step_size, 10.);
  EXPECT_DOUBLE_EQ(p3.weight, weight);

  // initialize with data
  int m = 10;
  int n = 20;
  Matrix Q(m, n, 0.99);
  Vector c(n, 1.99);
  step_size = 2.99;
  weight = 3.99;
  forward_grad_for_qp<Matrix> p4(&Q, &c, step_size, weight);
  EXPECT_DOUBLE_EQ(p4.step_size, step_size);
  EXPECT_DOUBLE_EQ(p4.weight, weight);

  //test for class SpMat
  forward_grad_for_qp<SpMat> p1_Sp;
  EXPECT_DOUBLE_EQ(p1_Sp.step_size, 0.);
  EXPECT_DOUBLE_EQ(p1_Sp.weight, 1.);

  // constructor initialized by step_size
  step_size = 10.;
  forward_grad_for_qp<SpMat> p2_Sp(step_size);
  EXPECT_DOUBLE_EQ(p2_Sp.step_size, 10.);
  EXPECT_DOUBLE_EQ(p2_Sp.weight, 1.);

  // constructor initialized by step_size and weight
  weight = 9.99;
  forward_grad_for_qp<SpMat> p3_Sp(step_size, weight);
  EXPECT_DOUBLE_EQ(p3_Sp.step_size, 10.);
  EXPECT_DOUBLE_EQ(p3_Sp.weight, weight);

  // initialize with data
  m = 10;
  n = 20;
  SpMat Q_Sp(m, n);
  step_size = 2.99;
  weight = 3.99;
  forward_grad_for_qp<SpMat> p4_Sp(&Q_Sp, &c, step_size, weight);
  EXPECT_DOUBLE_EQ(p4_Sp.step_size, step_size);
  EXPECT_DOUBLE_EQ(p4_Sp.weight, weight);
 
}

// test the update_step_size member function
TEST(ForwardQPTest, UpdateStepSize) {
  double step_size = 1.;
  // test for Matrix
  forward_grad_for_qp<Matrix> p(step_size);
  p.update_step_size(100.);
  EXPECT_DOUBLE_EQ(p.step_size, 100.);
  EXPECT_DOUBLE_EQ(p.weight, 1.);  

  // test for SpMat
  step_size = 1.;
  forward_grad_for_qp<SpMat> p_Sp(step_size);
  p_Sp.update_step_size(100.);
  EXPECT_DOUBLE_EQ(p_Sp.step_size, 100.);
  EXPECT_DOUBLE_EQ(p_Sp.weight, 1.);  
}

// test the scalar operator
TEST(ForwardQPTest, ScalarOperator) {

  double step_size = 1.;

  // test for Matrix
  forward_grad_for_qp<Matrix> p(step_size);
  // test numbers bigger than step_size
  EXPECT_DOUBLE_EQ(p(1.), DOUBLE_MARKER);
  // test numbers in [-step_size, step_size]
  EXPECT_DOUBLE_EQ(p(-1.), DOUBLE_MARKER);
  // test numbers smaller than -step_size
  EXPECT_DOUBLE_EQ(p(0.), DOUBLE_MARKER);
  // test with a given index
  EXPECT_DOUBLE_EQ(p(0., 100.), DOUBLE_MARKER);  

  // test for SpMat
  step_size = 1.;
  forward_grad_for_qp<SpMat> p_Sp(step_size);
  // test numbers bigger than step_size
  EXPECT_DOUBLE_EQ(p_Sp(1.), DOUBLE_MARKER);
  // test numbers in [-step_size, step_size]
  EXPECT_DOUBLE_EQ(p_Sp(-1.), DOUBLE_MARKER);
  // test numbers smaller than -step_size
  EXPECT_DOUBLE_EQ(p_Sp(0.), DOUBLE_MARKER);
  // test with a given index
  EXPECT_DOUBLE_EQ(p_Sp(0., 100.), DOUBLE_MARKER);  
}

// test the vector operator 
TEST(ForwardQPTest, VectorOperator) {
  
  int n = 2;
  Vector v(n, 2.);
  v[1] = -1;
  double step_size = 10.;
  int m = 2;
  
  // test for Matrix
  Matrix Q(m, n, 1.);
  Vector c(m, 1.);
  
  forward_grad_for_qp<Matrix> p(&Q, &c, step_size);
  EXPECT_DOUBLE_EQ(p(&v, 0),-18.);
  EXPECT_DOUBLE_EQ(v[0], 2.);
  
  p.update_step_size(1.);
  EXPECT_DOUBLE_EQ(p.step_size, 1.);
  EXPECT_DOUBLE_EQ(p(&v, 1), -3.);
 

  // test for SpMat
  std::vector<Eigen::Triplet<double> > triplets(m*n);
  int k = 0, i, j;
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      triplets[k++] = Eigen::Triplet<double> (i, j, 1.);
    }
  }

  SpMat Q_Sp(m, n);
  Q_Sp.setFromTriplets(triplets.begin(), triplets.end());

  forward_grad_for_qp<SpMat> p_Sp(&Q_Sp, &c, step_size);
  EXPECT_DOUBLE_EQ(p_Sp(&v, 0),-18.);
  EXPECT_DOUBLE_EQ(v[0], 2.);
  
  p_Sp.update_step_size(1.);
  EXPECT_DOUBLE_EQ(p_Sp.step_size, 1.);
  EXPECT_DOUBLE_EQ(p_Sp(&v, 1), -3.);
}

// test the full update operator 
TEST(ForwardQPTest, FullOperator) {
  
  int n = 2;
  int m = 2;

  // test for class Matrix
  Matrix Q(m, n, 1.);
  Vector c(m, 1.);
  
  Vector v_in(n, 2.), v_out(n, 1.);
  v_in[1] = -1.;

  double step_size = 10.;
  forward_grad_for_qp<Matrix> p(&Q, &c, step_size);
  p(&v_in, &v_out);
  
  EXPECT_DOUBLE_EQ(v_out[0], -18);
  EXPECT_DOUBLE_EQ(v_out[1], -21);  

  p.update_step_size(1.);
  EXPECT_DOUBLE_EQ(p.step_size, 1.);    
  p(&v_in, &v_out);
  EXPECT_DOUBLE_EQ(v_out[0], 0.);
  EXPECT_DOUBLE_EQ(v_out[1], -3.);  

  // test for class SpMat
  std::vector<Eigen::Triplet<double> > triplets(m*n);
  int k = 0, i, j;
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      triplets[k++] = Eigen::Triplet<double> (i, j, 1.);
    }
  }

  SpMat Q_Sp(m, n);
  Q_Sp.setFromTriplets(triplets.begin(), triplets.end());
  std::fill(v_in.begin(), v_in.end(), 2.);
  std::fill(v_out.begin(), v_out.end(), 1.);
  v_in[1] = -1.;
  forward_grad_for_qp<SpMat> p_Sp(&Q_Sp, &c, step_size);
  p_Sp(&v_in, &v_out);
  
  EXPECT_DOUBLE_EQ(v_out[0], -18);
  EXPECT_DOUBLE_EQ(v_out[1], -21);  

  p_Sp.update_step_size(1.);
  EXPECT_DOUBLE_EQ(p_Sp.step_size, 1.);    
  p_Sp(&v_in, &v_out);
  EXPECT_DOUBLE_EQ(v_out[0], 0.);
  EXPECT_DOUBLE_EQ(v_out[1], -3.);  
}

// test for update_cache_vars()
TEST(ForwardQPTest, UpdateCacheVars) {
  
}

// unit test for Jacobi method for linear equations
TEST(JacobiLinearEquationTest, Constructor) {
  //test for class Matrix
  // default constructor
  linear_eqn_jacobi_operator<Matrix> p1;
  EXPECT_DOUBLE_EQ(p1.step_size, 0.);
  EXPECT_DOUBLE_EQ(p1.weight, 1.);

  // constructor initialized by step_size
  double step_size = 10.;
  linear_eqn_jacobi_operator<Matrix> p2(step_size);
  EXPECT_DOUBLE_EQ(p2.step_size, 10.);
  EXPECT_DOUBLE_EQ(p2.weight, 1.);

  // constructor initialized by step_size and weight
  double weight = 9.99;
  linear_eqn_jacobi_operator<Matrix> p3(step_size, weight);
  EXPECT_DOUBLE_EQ(p3.step_size, 10.);
  EXPECT_DOUBLE_EQ(p3.weight, weight);

  // initialize with data
  int m = 20;
  int n = 20;
  Matrix Q(m, n, 0.99);
  Vector c(n, 1.99);
  step_size = 2.99;
  weight = 3.99;
  linear_eqn_jacobi_operator<Matrix> p4(&Q, &c, step_size, weight);
  EXPECT_DOUBLE_EQ(p4.step_size, step_size);
  EXPECT_DOUBLE_EQ(p4.weight, weight);

  //test for class SpMat
  linear_eqn_jacobi_operator<SpMat> p1_Sp;
  EXPECT_DOUBLE_EQ(p1_Sp.step_size, 0.);
  EXPECT_DOUBLE_EQ(p1_Sp.weight, 1.);

  // constructor initialized by step_size
  step_size = 10.;
  linear_eqn_jacobi_operator<SpMat> p2_Sp(step_size);
  EXPECT_DOUBLE_EQ(p2_Sp.step_size, 10.);
  EXPECT_DOUBLE_EQ(p2_Sp.weight, 1.);

  // constructor initialized by step_size and weight
  weight = 9.99;
  linear_eqn_jacobi_operator<SpMat> p3_Sp(step_size, weight);
  EXPECT_DOUBLE_EQ(p3_Sp.step_size, 10.);
  EXPECT_DOUBLE_EQ(p3_Sp.weight, weight);

  // initialize with data
  m = 20;
  n = 20;
  SpMat Q_Sp(m, n);
  step_size = 2.99;
  weight = 3.99;
  linear_eqn_jacobi_operator<SpMat> p4_Sp(&Q_Sp, &c, step_size, weight);
  EXPECT_DOUBLE_EQ(p4_Sp.step_size, step_size);
  EXPECT_DOUBLE_EQ(p4_Sp.weight, weight);
 
}

// test the update_step_size member function
TEST(JacobiLinearEquationTest, UpdateStepSize) {
  double step_size = 1.;
  // test for Matrix
  linear_eqn_jacobi_operator<Matrix> p(step_size);
  p.update_step_size(100.);
  EXPECT_DOUBLE_EQ(p.step_size, 100.);
  EXPECT_DOUBLE_EQ(p.weight, 1.);  

  // test for SpMat
  step_size = 1.;
  linear_eqn_jacobi_operator<SpMat> p_Sp(step_size);
  p_Sp.update_step_size(100.);
  EXPECT_DOUBLE_EQ(p_Sp.step_size, 100.);
  EXPECT_DOUBLE_EQ(p_Sp.weight, 1.);  
}

// test the scalar operator
TEST(JacobiLinearEquationTest, ScalarOperator) {

  double step_size = 1.;

  // test for Matrix
  linear_eqn_jacobi_operator<Matrix> p(step_size);
  // test numbers bigger than step_size
  EXPECT_DOUBLE_EQ(p(1.), DOUBLE_MARKER);
  // test numbers in [-step_size, step_size]
  EXPECT_DOUBLE_EQ(p(-1.), DOUBLE_MARKER);
  // test numbers smaller than -step_size
  EXPECT_DOUBLE_EQ(p(0.), DOUBLE_MARKER);
  // test with a given index
  EXPECT_DOUBLE_EQ(p(0., 100.), DOUBLE_MARKER);  

  // test for SpMat
  step_size = 1.;
  linear_eqn_jacobi_operator<SpMat> p_Sp(step_size);
  // test numbers bigger than step_size
  EXPECT_DOUBLE_EQ(p_Sp(1.), DOUBLE_MARKER);
  // test numbers in [-step_size, step_size]
  EXPECT_DOUBLE_EQ(p_Sp(-1.), DOUBLE_MARKER);
  // test numbers smaller than -step_size
  EXPECT_DOUBLE_EQ(p_Sp(0.), DOUBLE_MARKER);
  // test with a given index
  EXPECT_DOUBLE_EQ(p_Sp(0., 100.), DOUBLE_MARKER);  
}

// test the vector operator 
TEST(JacobiLinearEquationTest, VectorOperator) {
  
  int n = 2;
  Vector v(n, 2.);
  v[1] = -1;
  double step_size = 10.;
  
  // test for Matrix
  Matrix Q(n, n, 1.);
  Vector c(n, 1.);
  Q(0,0) = 2;
  
  linear_eqn_jacobi_operator<Matrix> p(&Q, &c, step_size);
  EXPECT_DOUBLE_EQ(p(&v, 0),1.);
  EXPECT_DOUBLE_EQ(v[0], 2.);
  
  p.update_step_size(1.);
  EXPECT_DOUBLE_EQ(p.step_size, 1.);
  EXPECT_DOUBLE_EQ(p(&v, 1), -1.);
 

  // test for SpMat
  std::vector<Eigen::Triplet<double> > triplets(n);
  int k = 0, i;
  for (i = 0; i < n; i++) {
      triplets[k++] = Eigen::Triplet<double> (i, i, 1.);
    }

  SpMat Q_Sp(n, n);
  Q_Sp.setFromTriplets(triplets.begin(), triplets.end());

  linear_eqn_jacobi_operator<SpMat> p_Sp(&Q_Sp, &c, step_size);
  EXPECT_DOUBLE_EQ(p_Sp(&v, 0),1.);
  EXPECT_DOUBLE_EQ(v[0], 2.);
  
  p_Sp.update_step_size(1.);
  EXPECT_DOUBLE_EQ(p_Sp.step_size, 1.);
  EXPECT_DOUBLE_EQ(p_Sp(&v, 1), 1.);
}

// test the full update operator 
TEST(JacobiLinearEquationTest, FullOperator) {
  
  int n = 2;
  int m = n;

  // test for class Matrix
  Matrix Q(m, n, 1.);
  Vector c(m, 1.);
  Q(0,0) = 2;
  Vector v_in(n, 2.), v_out(n, 1.);
  v_in[1] = -1.;

  double step_size = 10.;
  linear_eqn_jacobi_operator<Matrix> p(&Q, &c, step_size);
  p(&v_in, &v_out);
  
  EXPECT_DOUBLE_EQ(v_out[0], 1.);
  EXPECT_DOUBLE_EQ(v_out[1], -1.);  

  p.update_step_size(1.);
  EXPECT_DOUBLE_EQ(p.step_size, 1.);    
  p(&v_in, &v_out);
  EXPECT_DOUBLE_EQ(v_out[0], 1.);
  EXPECT_DOUBLE_EQ(v_out[1], -1.);  

  // test for class SpMat
  std::vector<Eigen::Triplet<double> > triplets(n);
  int k = 0, i;
  for (i = 0; i < n; i++) {
      triplets[k++] = Eigen::Triplet<double> (i, i, 1.);
  }

  SpMat Q_Sp(n, n);
  Q_Sp.setFromTriplets(triplets.begin(), triplets.end());
  std::fill(v_in.begin(), v_in.end(), 2.);
  std::fill(v_out.begin(), v_out.end(), 1.);
  v_in[1] = -1.;
  linear_eqn_jacobi_operator<SpMat> p_Sp(&Q_Sp, &c, step_size);
  p_Sp(&v_in, &v_out);
  
  EXPECT_DOUBLE_EQ(v_out[0], 1.);
  EXPECT_DOUBLE_EQ(v_out[1], 1.);  

  p_Sp.update_step_size(1.);
  EXPECT_DOUBLE_EQ(p_Sp.step_size, 1.);    
  p_Sp(&v_in, &v_out);
  EXPECT_DOUBLE_EQ(v_out[0], 1.);
  EXPECT_DOUBLE_EQ(v_out[1], 1.);  
}

// test for update_cache_vars()
TEST(JacobiLinearEquationTest, UpdateCacheVars) {
  
}

// test for forward gradient operator for log loss
TEST(ForwardLogLossTest, Constructor) {
  //test for class Matrix
  // default constructor
  forward_grad_for_log_loss<Matrix> p1;
  EXPECT_DOUBLE_EQ(p1.step_size, 0.);
  EXPECT_DOUBLE_EQ(p1.weight, 1.);

  // constructor initialized by step_size
  double step_size = 10.;
  forward_grad_for_log_loss<Matrix> p2(step_size);
  EXPECT_DOUBLE_EQ(p2.step_size, 10.);
  EXPECT_DOUBLE_EQ(p2.weight, 1.);

  // constructor initialized by step_size and weight
  double weight = 9.99;
  forward_grad_for_log_loss<Matrix> p3(step_size, weight);
  EXPECT_DOUBLE_EQ(p3.step_size, 10.);
  EXPECT_DOUBLE_EQ(p3.weight, weight);

  // initialize with data
  int m = 10;
  int n = 20;
  Matrix A(m, n, 0.99);
  Vector x(m, 1.99), Atx(n, 1.50);
  step_size = 2.99;
  weight = 3.99;
  forward_grad_for_log_loss<Matrix> p4(&A, &x, &Atx, step_size, weight);
  EXPECT_DOUBLE_EQ(p4.step_size, step_size);
  EXPECT_DOUBLE_EQ(p4.weight, weight);

  //test for class SpMat
  forward_grad_for_log_loss<SpMat> p1_Sp;
  EXPECT_DOUBLE_EQ(p1_Sp.step_size, 0.);
  EXPECT_DOUBLE_EQ(p1_Sp.weight, 1.);

  // constructor initialized by step_size
  step_size = 10.;
  forward_grad_for_log_loss<SpMat> p2_Sp(step_size);
  EXPECT_DOUBLE_EQ(p2_Sp.step_size, 10.);
  EXPECT_DOUBLE_EQ(p2_Sp.weight, 1.);

  // constructor initialized by step_size and weight
  weight = 9.99;
  forward_grad_for_log_loss<SpMat> p3_Sp(step_size, weight);
  EXPECT_DOUBLE_EQ(p3_Sp.step_size, 10.);
  EXPECT_DOUBLE_EQ(p3_Sp.weight, weight);

  // initialize with data
  m = 10;
  n = 20;
  SpMat A_Sp(m, n);
  step_size = 2.99;
  weight = 3.99;
  forward_grad_for_log_loss<SpMat> p4_Sp(&A_Sp, &x, &Atx, step_size, weight);
  EXPECT_DOUBLE_EQ(p4_Sp.step_size, step_size);
  EXPECT_DOUBLE_EQ(p4_Sp.weight, weight);
 
}

// test the update_step_size member function
TEST(ForwardLogLossTest, UpdateStepSize) {
  double step_size = 1.;
  // test for Matrix
  forward_grad_for_log_loss<Matrix> p(step_size);
  p.update_step_size(100.);
  EXPECT_DOUBLE_EQ(p.step_size, 100.);
  EXPECT_DOUBLE_EQ(p.weight, 1.);  

  // test for SpMat
  step_size = 1.;
  forward_grad_for_log_loss<SpMat> p_Sp(step_size);
  p_Sp.update_step_size(100.);
  EXPECT_DOUBLE_EQ(p_Sp.step_size, 100.);
  EXPECT_DOUBLE_EQ(p_Sp.weight, 1.);  
}


// test the scalar operator
TEST(ForwardLogLossTest, ScalarOperator) {
  
}


// test the vector operator 
TEST(ForwardLogLossTest, VectorOperator) {

}

// test the full update operator 
TEST(ForwardLogLossTest, FullOperator) {
 
}


// test for update_cache_vars()
TEST(ForwardLogLossTest, UpdateCacheVars) {
  
}

// test the forward gradient operator for least square

TEST(ForwardGradLS, Constructor) {
  // test for class Matrix
  // default constructor
  forward_grad_for_square_loss<Matrix> forward1;
  EXPECT_DOUBLE_EQ(forward1.step_size, 0.);
  EXPECT_DOUBLE_EQ(forward1.weight, 1.);
  
  // initialize with step_size and weight
  forward_grad_for_square_loss<Matrix> forward2(9.99, 8.99);
  EXPECT_DOUBLE_EQ(forward2.step_size, 9.99);
  EXPECT_DOUBLE_EQ(forward2.weight, 8.99);
  
  // initialize with data
  int m = 10;
  int n = 20;
  Matrix A(m, n, 0.99);
  Vector b(n, 1.99);
  Vector Atx(n, 0.);
  double step_size = 2.99;
  double weight = 3.99;
  forward_grad_for_square_loss<Matrix> forward3(&A, &b, &Atx, step_size, weight);
  EXPECT_DOUBLE_EQ(forward3.step_size, step_size);
  EXPECT_DOUBLE_EQ(forward3.weight, weight);

  // test for class SpMat
  // default constructor
  forward_grad_for_square_loss<SpMat> forward1_Sp;
  EXPECT_DOUBLE_EQ(forward1_Sp.step_size, 0.);
  EXPECT_DOUBLE_EQ(forward1_Sp.weight, 1.);
  
  // initialize with step_size and weight
  forward_grad_for_square_loss<SpMat> forward2_Sp(9.99, 8.99);
  EXPECT_DOUBLE_EQ(forward2_Sp.step_size, 9.99);
  EXPECT_DOUBLE_EQ(forward2_Sp.weight, 8.99);
  
  // initialize with data
  std::vector<Eigen::Triplet<double> > triplets(m*n);
  int k = 0, i, j;
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      triplets[k++] = Eigen::Triplet<double> (i, j, 0.99);
    }
  }

  SpMat A_Sp(m, n);
  A_Sp.setFromTriplets(triplets.begin(), triplets.end());

  step_size = 2.99;
  weight = 3.99;
  forward_grad_for_square_loss<SpMat> forward3_Sp(&A_Sp, &b, &Atx, step_size, weight);
  EXPECT_DOUBLE_EQ(forward3_Sp.step_size, step_size);
  EXPECT_DOUBLE_EQ(forward3_Sp.weight, weight);
}



// test the update_step_size member function
TEST(ForwardGradLS, UpdateStepSize) {
  double step_size = 1.;

  // test for class Matrix
  forward_grad_for_square_loss<Matrix> p(step_size);
  p.update_step_size(100.);
  EXPECT_DOUBLE_EQ(p.step_size, 100.);
  EXPECT_DOUBLE_EQ(p.weight, 1.);  

  // test for class SpMat
  forward_grad_for_square_loss<SpMat> p_Sp(step_size);
  p_Sp.update_step_size(100.);
  EXPECT_DOUBLE_EQ(p_Sp.step_size, 100.);
  EXPECT_DOUBLE_EQ(p_Sp.weight, 1.);  
}

// test the scalar operator
TEST(ForwardGradLS, ScalarOperator) {

  double step_size = 1.;

  
  forward_grad_for_square_loss<Matrix> p(step_size);
  // test numbers bigger than step_size
  EXPECT_DOUBLE_EQ(p(1.), DOUBLE_MARKER);
  // test numbers in [-step_size, step_size]
  EXPECT_DOUBLE_EQ(p(-1.), DOUBLE_MARKER);
  // test numbers smaller than -step_size
  EXPECT_DOUBLE_EQ(p(0.), DOUBLE_MARKER);
  // test with a given index
  EXPECT_DOUBLE_EQ(p(0., 100), DOUBLE_MARKER);  

  // test for class SpMat
  forward_grad_for_square_loss<SpMat> p_Sp(step_size);
  // test numbers bigger than step_size
  EXPECT_DOUBLE_EQ(p_Sp(1.), DOUBLE_MARKER);
  // test numbers in [-step_size, step_size]
  EXPECT_DOUBLE_EQ(p_Sp(-1.), DOUBLE_MARKER);
  // test numbers smaller than -step_size
  EXPECT_DOUBLE_EQ(p_Sp(0.), DOUBLE_MARKER);
  // test with a given index
  EXPECT_DOUBLE_EQ(p_Sp(0., 100), DOUBLE_MARKER);  
}

// test the vector operator 
TEST(ForwardGradLS, VectorOperator) {
  
  int n = 2;
  int m = 2;

  // test for class Matrix
  Matrix A(m, n, 1.);
  Vector b(n, 2.);
  Vector x(m, 0.);
  Vector Atx(n, 0.);
  double weight = 1.;
  double step_size = 10.;  
  forward_grad_for_square_loss<Matrix> forward(&A, &b, &Atx, step_size, weight);

  EXPECT_DOUBLE_EQ(forward(&x, 0), 40.);
  EXPECT_DOUBLE_EQ(x[0], 0.);
  
  forward.update_step_size(1.);
  EXPECT_DOUBLE_EQ(forward.step_size, 1.);
  EXPECT_DOUBLE_EQ(forward(&x, 1), 4.);

  // test for class SpMat
  std::vector<Eigen::Triplet<double> > triplets(m*n);
  int k = 0, i, j;
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      triplets[k++] = Eigen::Triplet<double> (i, j, 1.);
    }
  }

  SpMat A_Sp(m, n);
  A_Sp.setFromTriplets(triplets.begin(), triplets.end());

  std::fill(b.begin(), b.end(), 2.);
  std::fill(x.begin(), x.end(), 0.);
  std::fill(Atx.begin(), Atx.end(), 0.);
  weight = 1.;
  step_size = 10.;  
   
  forward_grad_for_square_loss<SpMat> forward_Sp(&A_Sp, &b, &Atx, step_size, weight);
  EXPECT_DOUBLE_EQ(forward_Sp(&x, 0), 40.);
  EXPECT_DOUBLE_EQ(x[0], 0.);
  
  forward_Sp.update_step_size(1.);
  EXPECT_DOUBLE_EQ(forward_Sp.step_size, 1.);
  EXPECT_DOUBLE_EQ(forward_Sp(&x, 1), 4.);
}



// test the full update operator 
TEST(ForwardGradLS, FullOperator) {

  int n = 2;
  int m = 2;

  // test for class Matrix
  Matrix A(m, n, 1.);
  Vector b(n, 2.);
  Vector v_in(n, 0.), v_out(n, 1.);  
  Vector Atx(n, 0.);
  double weight = 1.;
  double step_size = 10.;  
  forward_grad_for_square_loss<Matrix> forward(&A, &b, &Atx, step_size, weight);

  forward(&v_in, &v_out);
  
  EXPECT_DOUBLE_EQ(v_out[0], 40.);
  EXPECT_DOUBLE_EQ(v_out[1], 40.);  

  forward.update_step_size(1.);
  EXPECT_DOUBLE_EQ(forward.step_size, 1.);    
  forward(&v_in, &v_out);
  EXPECT_DOUBLE_EQ(v_out[0], 4.);
  EXPECT_DOUBLE_EQ(v_out[1], 4.); 

  //test for class SpMat 
  std::vector<Eigen::Triplet<double> > triplets(m*n);
  int k = 0, i, j;
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      triplets[k++] = Eigen::Triplet<double> (i, j, 1.);
    }
  }

  SpMat A_Sp(m, n);
  A_Sp.setFromTriplets(triplets.begin(), triplets.end());

  std::fill(b.begin(), b.end(), 2.);
  std::fill(v_in.begin(), v_in.end(), 0.);
  std::fill(v_out.begin(), v_out.end(), 1.);  
  std::fill(Atx.begin(), Atx.end(), 0.);
  weight = 1.;
  step_size = 10.; 
  forward_grad_for_square_loss<SpMat> forward_Sp(&A_Sp, &b, &Atx, step_size, weight);
  forward_Sp(&v_in, &v_out);
  
  EXPECT_DOUBLE_EQ(v_out[0], 40.);
  EXPECT_DOUBLE_EQ(v_out[1], 40.);  

  forward_Sp.update_step_size(1.);
  EXPECT_DOUBLE_EQ(forward_Sp.step_size, 1.);    
  v_out[0] = 0.;
  v_out[1] = 0.;
  forward_Sp(&v_in, &v_out);
  EXPECT_DOUBLE_EQ(v_out[0], 4.);
  EXPECT_DOUBLE_EQ(v_out[1], 4.); 
}

// test for update_cache_vars()
TEST(ForwardGradLS, UpdateCacheVars) {
  
}

// unit tests for the square hinge loss
// test for constructors
TEST(SquareHingeLossTest, Constructor) {
  // test for class Matrix
  // default constructor
  forward_grad_for_square_hinge_loss<Matrix> p1;
  EXPECT_DOUBLE_EQ(p1.step_size, 0.);
  EXPECT_DOUBLE_EQ(p1.weight, 1.);

  double step_size = 10.;
  
  // constructor initialized by step_size and weight
  double weight = 9.99;
  forward_grad_for_square_hinge_loss<Matrix> p3(step_size, weight);
  EXPECT_DOUBLE_EQ(p3.step_size, 10.);
  EXPECT_DOUBLE_EQ(p3.weight, weight);
  
  // initialize with data
  int m = 10;
  int n = 20;
  Matrix A(m, n, 0.99);
  Vector b(n, 1.);
  Vector Atx(m, 1.5);
  step_size = 2.99;
  weight = 3.99;
  forward_grad_for_square_hinge_loss<Matrix> p4(&A, &b, &Atx, step_size, weight);
  EXPECT_DOUBLE_EQ(p4.step_size, step_size);
  EXPECT_DOUBLE_EQ(p4.weight, weight);

  //test for class SpMat
 
  forward_grad_for_square_hinge_loss<SpMat> p1_Sp;
  EXPECT_DOUBLE_EQ(p1_Sp.step_size, 0.);
  EXPECT_DOUBLE_EQ(p1_Sp.weight, 1.);

  step_size = 10.;

  // constructor initialized by step_size and weight
  weight = 9.99;
  forward_grad_for_square_hinge_loss<SpMat> p3_Sp(step_size, weight);
  p3_Sp.update_step_size(10.0);
  EXPECT_DOUBLE_EQ(p3_Sp.step_size, 10.);
  EXPECT_DOUBLE_EQ(p3_Sp.weight, weight);

  // initialize with data
  m = 10;
  n = 20;
  SpMat A_Sp(m, n);
  step_size = 2.99;
  weight = 3.99;
  forward_grad_for_square_hinge_loss<SpMat> p4_Sp(&A_Sp, &b, &Atx, step_size, weight);
  EXPECT_DOUBLE_EQ(p4_Sp.step_size, step_size);
  EXPECT_DOUBLE_EQ(p4_Sp.weight, weight);
}

// test the update_step_size member function
TEST(SquareHingeLossTest, UpdateStepSize) {
  double step_size = 1.;
  // test for Matrix
  forward_grad_for_square_hinge_loss<Matrix> p;
  p.update_step_size(100.);
  EXPECT_DOUBLE_EQ(p.step_size, 100.);
  EXPECT_DOUBLE_EQ(p.weight, 1.);  

  // test for SpMat
  forward_grad_for_square_hinge_loss<SpMat> p_Sp;
  p_Sp.update_step_size(100.);
  EXPECT_DOUBLE_EQ(p_Sp.step_size, 100.);
  EXPECT_DOUBLE_EQ(p_Sp.weight, 1.);  
}

// test the scalar operator
TEST(SquareHingeLossTest, ScalarOperator) {

  double step_size = 1.;

  // test for Matrix
  forward_grad_for_square_hinge_loss<Matrix> p;
  // test numbers bigger than step_size
  EXPECT_DOUBLE_EQ(p(1.), DOUBLE_MARKER);
  // test numbers in [-step_size, step_size]
  EXPECT_DOUBLE_EQ(p(-1.), DOUBLE_MARKER);
  // test numbers smaller than -step_size
  EXPECT_DOUBLE_EQ(p(0.), DOUBLE_MARKER);
  // test with a given index
  EXPECT_DOUBLE_EQ(p(0., 100.), DOUBLE_MARKER);  

  // test for SpMat
  step_size = 1.;
  forward_grad_for_square_hinge_loss<SpMat> p_Sp;
  // test numbers bigger than step_size
  EXPECT_DOUBLE_EQ(p_Sp(1.), DOUBLE_MARKER);
  // test numbers in [-step_size, step_size]
  EXPECT_DOUBLE_EQ(p_Sp(-1.), DOUBLE_MARKER);
  // test numbers smaller than -step_size
  EXPECT_DOUBLE_EQ(p_Sp(0.), DOUBLE_MARKER);
  // test with a given index
  EXPECT_DOUBLE_EQ(p_Sp(0., 100.), DOUBLE_MARKER);  
}

// test the vector operator 
TEST(SquareHingeLossTest, VectorOperator) {
  
  int n = 2;
  int m = 2;

  // test for class Matrix
  Matrix A(m, n, 1.);
  Vector b(n, 1.);
  Vector x(m, 0.);
  Vector Atx(n, 0.);
  double weight = 1.;
  double step_size = 10.;  
  forward_grad_for_square_hinge_loss<Matrix> forward(&A, &b, &Atx, step_size, weight);

  EXPECT_DOUBLE_EQ(forward(&x, 0), 20.);
  EXPECT_DOUBLE_EQ(x[0], 0.);
  
  forward.update_step_size(1.);
  EXPECT_DOUBLE_EQ(forward.step_size, 1.);
  EXPECT_DOUBLE_EQ(forward(&x, 1), 2.);

  // test for class SpMat
  std::vector<Eigen::Triplet<double> > triplets(m*n);
  int k = 0, i, j;
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      triplets[k++] = Eigen::Triplet<double> (i, j, 1.);
    }
  }
  
  SpMat A_Sp(m, n);
  A_Sp.setFromTriplets(triplets.begin(), triplets.end());

  std::fill(b.begin(), b.end(), 1.);
  b[0] = -1.;
  std::fill(x.begin(), x.end(), 0.);
  std::fill(Atx.begin(), Atx.end(), 0.);
  weight = 1.;
  step_size = 10.; 
  forward_grad_for_square_hinge_loss<SpMat> forward_Sp(&A_Sp, &b, &Atx, step_size, weight);

  EXPECT_DOUBLE_EQ(forward_Sp(&x, 0), 0.);
  EXPECT_DOUBLE_EQ(x[0], 0.);
  
  forward_Sp.update_step_size(1.);
  EXPECT_DOUBLE_EQ(forward_Sp.step_size, 1.);
  EXPECT_DOUBLE_EQ(forward_Sp(&x, 1), 0.);
}



// test the full update operator 
TEST(SquareHingeLossTest, FullOperator) {

  int n = 2;
  int m = 2;

  // test for class Matrix
  Matrix A(m, n, 1.);
  Vector b(n, 1.);

  Vector v_in(n, 0.), v_out(n, 1.);  
  Vector Atx(n, 0.);
  double weight = 1.;
  double step_size = 10.;  
  forward_grad_for_square_hinge_loss<Matrix> forward(&A, &b, &Atx, step_size, weight);

  forward(&v_in, &v_out);
  
  EXPECT_DOUBLE_EQ(v_out[0], 20.);
  EXPECT_DOUBLE_EQ(v_out[1], 20.);  

  forward.update_step_size(1.);
  EXPECT_DOUBLE_EQ(forward.step_size, 1.);    
  forward(&v_in, &v_out);
  EXPECT_DOUBLE_EQ(v_out[0], 2.);
  EXPECT_DOUBLE_EQ(v_out[1], 2.);  

  // test for class SpMat
  std::vector<Eigen::Triplet<double> > triplets(m*n);
  int k = 0, i, j;
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      triplets[k++] = Eigen::Triplet<double> (i, j, 1.);
    }
  }
  
  SpMat A_Sp(m, n);
  A_Sp.setFromTriplets(triplets.begin(), triplets.end());

  std:fill(b.begin(), b.end(), 1.);
  std::fill(v_in.begin(), v_in.end(), 0.);
  std::fill(v_out.begin(), v_out.end(), 1.);  
  std::fill(Atx.begin(), Atx.end(), 0.);
  weight = 1.;
  step_size = 10.;  

  forward_grad_for_square_hinge_loss<SpMat> forward_Sp(&A_Sp, &b, &Atx, step_size, weight);

  forward_Sp(&v_in, &v_out);
  EXPECT_DOUBLE_EQ(v_out[0], 20.);
  EXPECT_DOUBLE_EQ(v_out[1], 20.);  

  forward_Sp.update_step_size(1.);
  EXPECT_DOUBLE_EQ(forward_Sp.step_size, 1.);    
  forward_Sp(&v_in, &v_out);
  EXPECT_DOUBLE_EQ(v_out[0], 2.);
  EXPECT_DOUBLE_EQ(v_out[1], 2.);    
}


// unit tests for the Huber loss
// test for constructors
TEST(HuberLossTest, Constructor) {
  //test for class Matrix
  // default constructor
  forward_grad_for_huber_loss<Matrix> p1;
  EXPECT_DOUBLE_EQ(p1.step_size, 0.);
  EXPECT_DOUBLE_EQ(p1.weight, 1.);
 
  // initialize with data
  int m = 10;
  int n = 20;
  Matrix A(m, n, 0.99);
  Vector b(n, 1.99);
  Vector Atx(m, 1.5);
  double step_size = 2.99;
  double weight = 3.99;
  forward_grad_for_huber_loss<Matrix> p4(&A, &b, &Atx, step_size, weight);
  EXPECT_DOUBLE_EQ(p4.step_size, step_size);
  EXPECT_DOUBLE_EQ(p4.weight, weight);

  //test for class SpMat
 
  forward_grad_for_huber_loss<SpMat> p1_Sp;
  EXPECT_DOUBLE_EQ(p1_Sp.step_size, 0.);
  EXPECT_DOUBLE_EQ(p1_Sp.weight, 1.);

  // initialize with data
  m = 10;
  n = 20;
  SpMat A_Sp(m, n);
  step_size = 2.99;
  weight = 3.99;
  forward_grad_for_huber_loss<SpMat> p4_Sp(&A_Sp, &b, &Atx, step_size, weight);
  EXPECT_DOUBLE_EQ(p4_Sp.step_size, step_size);
  EXPECT_DOUBLE_EQ(p4_Sp.weight, weight);
}

// test the update_step_size member function
TEST(HuberLossTest, UpdateStepSize) {
  double step_size = 1.;
  // test for Matrix
  forward_grad_for_huber_loss<Matrix> p;
  p.update_step_size(100.);
  EXPECT_DOUBLE_EQ(p.step_size, 100.);
  EXPECT_DOUBLE_EQ(p.weight, 1.);  

  // test for SpMat
  forward_grad_for_huber_loss<SpMat> p_Sp;
  p_Sp.update_step_size(100.);
  EXPECT_DOUBLE_EQ(p_Sp.step_size, 100.);
  EXPECT_DOUBLE_EQ(p_Sp.weight, 1.);  

}

// test the scalar operator
TEST(HuberLossTest, ScalarOperator) {

  double step_size = 1.;

  // test for Matrix
   forward_grad_for_huber_loss<Matrix> p;
  // test numbers bigger than step_size
  EXPECT_DOUBLE_EQ(p(1.), DOUBLE_MARKER);
  // test numbers in [-step_size, step_size]
  EXPECT_DOUBLE_EQ(p(-1.), DOUBLE_MARKER);
  // test numbers smaller than -step_size
  EXPECT_DOUBLE_EQ(p(0.), DOUBLE_MARKER);
  // test with a given index
  EXPECT_DOUBLE_EQ(p(0., 100.), DOUBLE_MARKER);  

  // test for SpMat
  step_size = 1.;
  forward_grad_for_huber_loss<SpMat> p_Sp;
  // test numbers bigger than step_size
  EXPECT_DOUBLE_EQ(p_Sp(1.), DOUBLE_MARKER);
  // test numbers in [-step_size, step_size]
  EXPECT_DOUBLE_EQ(p_Sp(-1.), DOUBLE_MARKER);
  // test numbers smaller than -step_size
  EXPECT_DOUBLE_EQ(p_Sp(0.), DOUBLE_MARKER);
  // test with a given index
  EXPECT_DOUBLE_EQ(p_Sp(0., 100.), DOUBLE_MARKER);  
}
// test the vector operator 
TEST(HuberLossTest, VectorOperator) {
  
  int n = 2;
  int m = 2;

  // test for class Matrix
  Matrix A(m, n, 1.);
  Vector b(n, 2.);
  b[0] = -2.;
  Vector x(m, 0.);
  Vector Atx(n, 0.);
  double weight = 1.;
  double step_size = 10.;
  double delta = 1.;
  forward_grad_for_huber_loss<Matrix> forward(&A, &b, &Atx, step_size, weight, delta);

  EXPECT_DOUBLE_EQ(forward(&x, 0), 0.);
  EXPECT_DOUBLE_EQ(x[0], 0.);
  
  forward.update_step_size(1.);
  EXPECT_DOUBLE_EQ(forward.step_size, 1.);
  EXPECT_DOUBLE_EQ(forward(&x, 1), 0.);

  // test for class SpMat
  std::vector<Eigen::Triplet<double> > triplets(m*n);
  int k = 0, i, j;
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      triplets[k++] = Eigen::Triplet<double> (i, j, 1.);
    }
  }

  SpMat A_Sp(m, n);
  A_Sp.setFromTriplets(triplets.begin(), triplets.end());
  
  std::fill(b.begin(), b.end(), 2.);
  b[0] = -2.;
  std::fill(x.begin(), x.end(), 0.);
  std::fill(Atx.begin(), Atx.end(), 0.);
  weight = 1.;
  step_size = 10.;
  delta = 1.;
  forward_grad_for_huber_loss<SpMat> forward_Sp(&A_Sp, &b, &Atx, step_size, weight, delta);

  EXPECT_DOUBLE_EQ(forward(&x, 0), 0.);
  EXPECT_DOUBLE_EQ(x[0], 0.);
}



// test the full update operator 
TEST(HuberLossTest, FullOperator) {
  
  int n = 2;
  int m = 2;

  // test for class Matrix
  Matrix A(m, n, 1.);
  Vector b(n, 2.);
  Vector v_in(n, 0.), v_out(n, 1.);  
  Vector Atx(n, 0.);
  double weight = 1.;
  double step_size = 10.;
  double delta = 1.;
  forward_grad_for_huber_loss<Matrix> forward(&A, &b, &Atx, step_size, weight);

  forward(&v_in, &v_out);
  
  EXPECT_DOUBLE_EQ(v_out[0], 20.);
  EXPECT_DOUBLE_EQ(v_out[1], 20.);  

  forward.update_step_size(1.);
  EXPECT_DOUBLE_EQ(forward.step_size, 1.);    
  forward(&v_in, &v_out);
  EXPECT_DOUBLE_EQ(v_out[0], 2.);
  EXPECT_DOUBLE_EQ(v_out[1], 2.);  

  // test for class SpMat
  std::vector<Eigen::Triplet<double> > triplets(m*n);
  int k = 0, i, j;
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      triplets[k++] = Eigen::Triplet<double> (i, j, 1.);
    }
  }

  SpMat A_Sp(m, n);
  A_Sp.setFromTriplets(triplets.begin(), triplets.end());

  std::fill(v_in.begin(), v_in.end(), 0.);
  std::fill(v_out.begin(), v_out.end(), 1.);  
  std::fill(Atx.begin() , Atx.end(), 0.);
  weight = 1.;
  step_size = 10.;
  delta = 1.;
  forward_grad_for_huber_loss<SpMat> forward_Sp(&A_Sp, &b, &Atx, step_size, weight);

  forward_Sp(&v_in, &v_out);
  
  EXPECT_DOUBLE_EQ(v_out[0], 20.);
  EXPECT_DOUBLE_EQ(v_out[1], 20.);  

  forward_Sp.update_step_size(1.);
  EXPECT_DOUBLE_EQ(forward_Sp.step_size, 1.);    
  forward_Sp(&v_in, &v_out);
  EXPECT_DOUBLE_EQ(v_out[0], 2.);
  EXPECT_DOUBLE_EQ(v_out[1], 2.);  
}

// test for update_cache_vars()
TEST(HuberLossTest, UpdateCacheVars) {

}


// unit test for forward operator for quadratic function
TEST(PortfolioTest, Constructor) {
  //test for class Matrix
  // default constructor
  portfolio_3s<Matrix> p1;
  EXPECT_DOUBLE_EQ(p1.step_size, 0.);
  EXPECT_DOUBLE_EQ(p1.weight, 1.);

  // constructor initialized by step_size

  double step_size = 10.;
  double weight = 9.99;  

  portfolio_3s<Matrix> p2(step_size);
  EXPECT_DOUBLE_EQ(p2.step_size, 10.);
  EXPECT_DOUBLE_EQ(p2.weight, 1.);

  // constructor initialized by step_size and weight

  portfolio_3s<Matrix> p3(step_size, weight);
  EXPECT_DOUBLE_EQ(p3.step_size, 10.);
  EXPECT_DOUBLE_EQ(p3.weight, weight);

  //test for class SpMat
  portfolio_3s<SpMat> p1_Sp;
  EXPECT_DOUBLE_EQ(p1_Sp.step_size, 0.);
  EXPECT_DOUBLE_EQ(p1_Sp.weight, 1.);

  // constructor initialized by step_size
  step_size = 10.;
  portfolio_3s<SpMat> p2_Sp(step_size);
  EXPECT_DOUBLE_EQ(p2_Sp.step_size, 10.);
  EXPECT_DOUBLE_EQ(p2_Sp.weight, 1.);

  // constructor initialized by step_size and weight
  weight = 9.99;
  portfolio_3s<SpMat> p3_Sp(step_size, weight);
  EXPECT_DOUBLE_EQ(p3_Sp.step_size, 10.);
  EXPECT_DOUBLE_EQ(p3_Sp.weight, weight);

  // initialize with data
  int n = 2;
  Matrix Q(n, n, 0.99);
  Vector epsilon(n, 1.99);
  double c = 0.1;
  step_size = 2.99;
  weight = 3.99;
  portfolio_3s<Matrix> p4(&Q, &epsilon, c, step_size, weight);
  EXPECT_DOUBLE_EQ(p4.step_size, step_size);
  EXPECT_DOUBLE_EQ(p4.weight, weight);

  // initialize with data

  SpMat Q_Sp(n, n);
  std::vector<Eigen::Triplet<double> > triplets(1);
  triplets.push_back(Eigen::Triplet<double> (0, 0, 1.));
  Q_Sp.setFromTriplets(triplets.begin(), triplets.end());  
  
  step_size = 2.99;
  weight = 3.99;
  portfolio_3s<SpMat> p4_Sp(&Q_Sp, &epsilon, c, step_size, weight);
  EXPECT_DOUBLE_EQ(p4_Sp.step_size, step_size);
  EXPECT_DOUBLE_EQ(p4_Sp.weight, weight);
  
}

// test the update_step_size member function
TEST(PortfolioTest, UpdateStepSize) {
  double step_size = 1.;
  // test for Matrix
  portfolio_3s<Matrix> p(step_size);
  p.update_step_size(100.);
  EXPECT_DOUBLE_EQ(p.step_size, 100.);
  EXPECT_DOUBLE_EQ(p.weight, 1.);  

  // test for SpMat
  step_size = 1.;
  portfolio_3s<SpMat> p_Sp(step_size);
  
  p_Sp.update_step_size(100.);
  EXPECT_DOUBLE_EQ(p_Sp.step_size, 100.);
  EXPECT_DOUBLE_EQ(p_Sp.weight, 1.);  
}

// test the scalar operator
TEST(PortfolioTest, ScalarOperator) {

  double step_size = 1.;

  // test for Matrix
  portfolio_3s<Matrix> p(step_size);
  // test numbers bigger than step_size
  EXPECT_DOUBLE_EQ(p(1.), DOUBLE_MARKER);
  // test numbers in [-step_size, step_size]
  EXPECT_DOUBLE_EQ(p(-1.), DOUBLE_MARKER);
  // test numbers smaller than -step_size
  EXPECT_DOUBLE_EQ(p(0.), DOUBLE_MARKER);
  // test with a given index
  EXPECT_DOUBLE_EQ(p(0., 100.), DOUBLE_MARKER);  

  // test for SpMat
  step_size = 1.;
  portfolio_3s<SpMat> p_Sp(step_size);
  // test numbers bigger than step_size
  EXPECT_DOUBLE_EQ(p_Sp(1.), DOUBLE_MARKER);
  // test numbers in [-step_size, step_size]
  EXPECT_DOUBLE_EQ(p_Sp(-1.), DOUBLE_MARKER);
  // test numbers smaller than -step_size
  EXPECT_DOUBLE_EQ(p_Sp(0.), DOUBLE_MARKER);
  // test with a given index
  EXPECT_DOUBLE_EQ(p_Sp(0., 100.), DOUBLE_MARKER);  
}

// test the vector operator 
TEST(PortfolioTest, VectorOperator) {
  
  int n = 2;
  Vector v(n, 0.);
  double step_size = 1.;
  int m = 2;
  
  // test for Matrix
  Matrix Q(m, n, 1.);
  Vector c(m, 1.);
  
  portfolio_3s<Matrix> p(&Q, &c, step_size);
  EXPECT_DOUBLE_EQ(p(&v, 0),-0.5);
  EXPECT_DOUBLE_EQ(v[0], 0.);
  
  p.update_step_size(10.);
  EXPECT_DOUBLE_EQ(p.step_size, 10.);
  EXPECT_DOUBLE_EQ(p(&v, 1), -0.5);
 

  // test for SpMat
  std::vector<Eigen::Triplet<double> > triplets(m*n);
  int k = 0, i, j;
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      triplets[k++] = Eigen::Triplet<double> (i, j, 1.);
    }
  }

  SpMat Q_Sp(m, n);
  Q_Sp.setFromTriplets(triplets.begin(), triplets.end());

  portfolio_3s<SpMat> p_Sp(&Q_Sp, &c, step_size);
  EXPECT_DOUBLE_EQ(p_Sp(&v, 0),-0.5);
  EXPECT_DOUBLE_EQ(v[0], 0.);
  
  p_Sp.update_step_size(1.);
  EXPECT_DOUBLE_EQ(p_Sp.step_size, 1.);
  EXPECT_DOUBLE_EQ(p_Sp(&v, 1), -0.5);
}

// test the full update operator 
TEST(PortfolioTest, FullOperator) {
  

}

// test for update_cache_vars()
TEST(PortfolioTest, UpdateCacheVars) {
  
}
