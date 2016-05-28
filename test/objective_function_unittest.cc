#include "gtest/gtest.h"
#include "matrices.h"
#include "algebra.h"
#include "operators.h"
#include "constants.h"
#include "util.h"
#include "math.h"
#include "algebra_namespace_switcher.h"
// used for DEBUG
#include<iostream>


/**************************************
   Unit test for objective function calculators
***************************************/

// l2_loss_log_objective
TEST(L2_loss_log_objective, Constructor) {
  int n = 10;
  Vector b(n, -1.);
  Vector Atx(n, log(exp(2) - 1));
  Vector x(n, -2.);
  double lambda = 2.0;

  EXPECT_DOUBLE_EQ(l2_log_loss_objective(b, x, Atx, lambda), 6 * n);
}

// l1_loss_log_objective
TEST(L1_loss_log_objective, Constructor) {
  int n = 10;
  Vector b(n, -1.);
  Vector Atx(n, log(exp(2) - 1));
  Vector x(n, -2.);
  double lambda = 3.0;

  EXPECT_DOUBLE_EQ(l1_log_loss_objective(b, x, Atx, lambda), 8 * n);
}

// logistic loss
TEST(ObjectiveFunctionOperators, LogisticLoss) {
  int n = 10;
  Vector b(n, -1.);
  Vector Atx(n, log(exp(2) - 1));
  Vector x(n, -2.);

  EXPECT_DOUBLE_EQ(log_loss(x, Atx, b), 2 * n);
}


// square loss
TEST(ObjectiveFunctionOperators, SquareLoss) {
  int n = 10;
  Vector b(n, -1.);
  Vector Atx(n, 2.);
  Vector x(n, -2.);

  EXPECT_DOUBLE_EQ(square_loss(x, Atx, b), 4.5 * n);
}

// square hinge loss
TEST(ObjectiveFunctionOperators, SquareHingeLoss) {
  int n = 10;
  Vector b(n, -1.);
  Vector Atx(n, 2.);
  Vector x(n, -2.);

  EXPECT_DOUBLE_EQ(square_hinge_loss(x, Atx, b), 4 * n);
  int n2 = 10;
  Vector b2(n, -1.);
  Vector Atx2(n, .5);
  Vector x2(n, -2.);

  EXPECT_DOUBLE_EQ(square_hinge_loss(x2, Atx2, b2),  n2);
}

// huber loss
TEST(ObjectiveFunctionOperators, HuberLoss) {
  int n = 10;
  Vector b(n, -1.);
  Vector Atx(n, 2.);
  Vector x(n, -2.);
  double delta = 1.5;

  EXPECT_DOUBLE_EQ(huber_loss(x, Atx, b, delta), 1.5 * 2.25 * n);
  int n2 = 10;
  Vector b2(n, -.5);
  Vector Atx2(n, .5);
  Vector x2(n, -2.);

  EXPECT_DOUBLE_EQ(huber_loss(x2, Atx2, b2, delta), .5 * n2);
}

// quadratic function
TEST(ObjectiveFunctionOperators, QuadraticFunction) {
  int n = 10;
  double d = 2.;
  Vector Qx(n, 3.);
  Vector x(n, 3.);
  Vector c(n, -1.);

  EXPECT_DOUBLE_EQ(quad_func(x, Qx, c, d), 1.5 * n + d);
}

// huber norm
TEST(ObjectiveFunctionOperators, HuberNorm) {
  
}

