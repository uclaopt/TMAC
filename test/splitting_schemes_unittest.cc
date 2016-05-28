#include <iostream>
#include "matrices.h"
#include "algebra.h"
#include "operators.h"
#include "gtest/gtest.h"
#include "constants.h"
#include "algebra_namespace_switcher.h"
#include "parameters.h"
#include "splitting_schemes.h"

/**************************************
   Unit test for splitting schemes
***************************************/


// test on PPA ; initialization
TEST(ProximalPointAlgorithmTest, Initialization) {

    int n = 3;
    Vector v(n, 2.99);
    double lambda = 0.1;
    // initialize the proximal operator
    prox_l1 forward(lambda);
    // initialize the proximal point algorithm
    ProximalPointAlgorithm<prox_l1> ppa(&v, forward, lambda);
    Params params;
    ppa.update_params(&params);
    // test the initializtion
    EXPECT_DOUBLE_EQ(2.99, (*ppa.x)[0]);
    EXPECT_DOUBLE_EQ(0.1, lambda);
    EXPECT_DOUBLE_EQ(0.618, ppa.relaxation_step_size);

}


// test on PPA : Functor
TEST(ProximalPointAlgorithmTest, Functor) {

    int n = 3;
    Vector v(n, 2.99);
    double lambda = 0.1;
    // initialize the proximal operator
    prox_l1 forward(lambda);
    // initialize the proximal point algorithm
    ProximalPointAlgorithm<prox_l1> ppa(&v, forward, lambda);
    Params params;
    params.step_size = lambda;
    ppa.update_params(&params);
    // test update element 0;
    ppa(0);
    EXPECT_DOUBLE_EQ(2.9282, (*ppa.x)[0]);
    
}


// test on Gradient Descent : Initialization 
TEST(GradientDescentAlgorithm, Initialization) {

    // Initialize the square_loss
    Matrix A(3, 3, 0);
    Vector x, b, Atx;
    double data_A[9] = {1., 2., 3., \
                        4., 5., 6., \
                        7., 8., 9.};
    double data_b[3] = {0., 0., 1.};
    double data_x[3] = {0., 1., -1.};
    double data_Atx[3] = {-3., -3., -3.};
    A.assign(data_A, data_A + 9);
    b.assign(data_b, data_b + 3);
    x.assign(data_x, data_x + 3);
    Atx.assign(data_Atx, data_Atx+3);
    forward_grad_for_square_loss<Matrix> forward(&A, &x, &Atx);
    // gradient descent algorithm
    GradientDescentAlgorithm<forward_grad_for_square_loss<Matrix>> gda(&x, forward);
    Params params;
    gda.update_params(&params);
    EXPECT_DOUBLE_EQ(0., (*gda.x)[0]);
    EXPECT_DOUBLE_EQ(1., (*gda.forward.A)(0, 0));
    EXPECT_DOUBLE_EQ(0.618, gda.relaxation_step_size);

}


// test on Gradient Descent : Functor
TEST(GradientDescentAlgorithm, Functor) {

    // Initialize the forward operator
    Matrix A(3, 3, 0);
    Vector x, b, Atx;
    double data_A[9] = {1., 2., 3., \
                        4., 5., 6., \
                        7., 8., 9.};
    double data_b[3] = {0., 0., 1.};
    double data_x[3] = {0., 1., -1.};
    double data_Atx[3] = {-3., -3., -3.};
    A.assign(data_A, data_A + 9);
    b.assign(data_b, data_b + 3);
    x.assign(data_x, data_x + 3);
    Atx.assign(data_Atx, data_Atx + 3);
    forward_grad_for_square_loss<Matrix> forward(&A, &b, &Atx);
    // gradient descent algorithm
    GradientDescentAlgorithm<forward_grad_for_square_loss<Matrix>> gda(&x, forward);
    Params params;
    params.step_size = 1.;
    gda.update_params(&params);
    gda(0);
    EXPECT_DOUBLE_EQ(12.978, x[0]);

}


// test on forward-backward splitting : Initialization
TEST(ForwardBackwardSplitting, Initialization) {

    // Initialize the forward and backward operator
    // In this test, the forward operator is the gradient operator for squared loss
    //               the backward operator is the prximal operator for L1 norm
    Matrix A(3, 3, 0);
    Vector x, b, Atx;
    double data_x[3] = {0., 1., -1.}; 
    double data_A[9] = {1., 0., 0., \
                        0., 1., 0., \
                        0., 0., 1.};   
    double data_b[3] = {0., 0., 0.};
    double data_Atx[3] = {0., 1., -1.};
    A.assign(data_A, data_A + 9);
    b.assign(data_b, data_b + 3);
    x.assign(data_x, data_x + 3);
    Atx.assign(data_Atx, data_Atx + 3);
    forward_grad_for_square_loss<Matrix> forward(&A, &b, &Atx);
    prox_l1 backward(1., 1.);
    ForwardBackwardSplitting<forward_grad_for_square_loss<Matrix>, prox_l1> fbS(&x, forward, backward);
    Params params;
    params.step_size = 1.;
    fbS.update_params(&params);
    EXPECT_DOUBLE_EQ(-1., (*fbS.x)[2]);
    EXPECT_DOUBLE_EQ(-1., (*fbS.forward.Atx)[2]);
    EXPECT_DOUBLE_EQ(1., fbS.backward.weight);
    EXPECT_DOUBLE_EQ(1., fbS.forward.step_size);
    EXPECT_DOUBLE_EQ(1., fbS.backward.step_size);
    EXPECT_DOUBLE_EQ(0.618, fbS.relaxation_step_size);

}


// test on forward-backward splitting : Functor
TEST(ForwardBackwardSplitting, Functor) {

    // Initialize the forward and backward operator
    // In this test, the forward operator is the gradient operator for squared loss
    //               the backward operator is the prximal operator for L1 norm
    Matrix A(3, 3, 0);
    Vector x, b, Atx;
    double data_x[3] = {0., 1., -1.}; 
    double data_A[9] = {1., 0., 0., \
                        0., 1., 0., \
                        0., 0., 1.};   
    double data_b[3] = {0., 0., 0.};
    double data_Atx[3] = {0., 1., -1.};
    A.assign(data_A, data_A + 9);
    b.assign(data_b, data_b + 3);
    x.assign(data_x, data_x + 3);
    Atx.assign(data_Atx, data_Atx + 3);
    forward_grad_for_square_loss<Matrix> forward(&A, &b, &Atx);
    prox_l1 backward(1., 1.);
    ForwardBackwardSplitting<forward_grad_for_square_loss<Matrix>, prox_l1> fbS(&x, forward, backward);
    Params params;
    params.step_size = 1.;
    fbS.update_params(&params);
    fbS(2);
    EXPECT_DOUBLE_EQ(-0.382, x[2]);

}


// test on backward-forward splitting : Initialization
TEST(BackwardForwardSplitting, Initialization) {

    // Initialize the forward and backward operator
    // In this test, the forward operator is the gradient operator for quadratic loss
    //               the backward operator is the prximal operator for L1 norm
    Matrix Q(3, 3, 0);
    Vector c, x;
    double data_Q[9] = {1., 0., 0., \
                        0., 1., 0., \
                        0., 0., 1.};   
    double data_c[3] = {0., 0., 0.};
    double data_x[3] = {0., 1., -3.};
    Q.assign(data_Q, data_Q + 9);
    c.assign(data_c, data_c + 3);
    x.assign(data_x, data_x + 3);
    forward_grad_for_qp<Matrix> forward(&Q, &c);
    prox_l1 backward(1., 1.);
    BackwardForwardSplitting<prox_l1, forward_grad_for_qp<Matrix>> bfS(&x, backward, forward);
    // Initialize the forward and backward operator
    Params params;
    params.step_size = 1.;
    bfS.update_params(&params);

    EXPECT_DOUBLE_EQ(-3., (*bfS.x)[2]);
    EXPECT_DOUBLE_EQ(1., bfS.backward.weight);
    EXPECT_DOUBLE_EQ(1., (*bfS.forward.Q)(1, 1));
    EXPECT_DOUBLE_EQ(1., bfS.backward.step_size);
    EXPECT_DOUBLE_EQ(1., bfS.forward.step_size);
    EXPECT_DOUBLE_EQ(0.618, bfS.relaxation_step_size);

}


// test on backward-forward splitting : Functor
TEST(BackwardForwardSplitting, Functor) {

    // Initialize the forward and backward operator
    // In this test, the backward operator is the prximal operator for L1 norm
    //               the forward operator is the gradient operator for quadratic loss
    Matrix Q(3, 3, 0);
    Vector c, x;
    double data_Q[9] = {1., 0., 0., \
                        0., 1., 0., \
                        0., 0., 1.};   
    double data_c[3] = {0., 0., 0.};
    double data_x[3] = {0., 1., -3.};
    Q.assign(data_Q, data_Q + 9);
    c.assign(data_c, data_c + 3);
    x.assign(data_x, data_x + 3);
    forward_grad_for_qp<Matrix> forward(&Q, &c);
    prox_l1 backward(1., 1.);
    BackwardForwardSplitting<prox_l1, forward_grad_for_qp<Matrix>> bfS(&x, backward, forward);
    // Initialize the forward and backward operator
    Params params;
    params.step_size = 1.;
    bfS.update_params(&params);
    bfS(2);
    EXPECT_DOUBLE_EQ(-1.146, x[2]);
}


// test on Peaceman-Rachford splitting : Initialization
TEST(PeacemanRachford, Initialization) {

    // Initialize the Peaceman-Rachford operator
    // In this test, two operators are both prximal operator for L1 norm
    prox_l1 first(2., 1.);
    prox_l1 second(2., 2.);
    Vector x;
    double data_x[3] = {0., 1., -1.}; 
    x.assign(data_x, data_x + 3);

    PeacemanRachfordSplitting<prox_l1, prox_l1> PRs(&x, first, second);
    // Initialize the Peaceman-Rachford operator
    Params params;
    params.step_size = 1.;
    PRs.update_params(&params);

    EXPECT_DOUBLE_EQ(-1., (*PRs.x)[2]);
    EXPECT_DOUBLE_EQ(1., PRs.op1.step_size);
    EXPECT_DOUBLE_EQ(1., PRs.op2.step_size);
    EXPECT_DOUBLE_EQ(1., PRs.op1.weight);
    EXPECT_DOUBLE_EQ(2., PRs.op2.weight);
    EXPECT_DOUBLE_EQ(0.618, PRs.relaxation_step_size);

}


// test on Peaceman-Rachford splitting : Functor
TEST(PeacemanRachford, Functor) {

    // Initialize the Peaceman-Rachford operator
    // In this test, two operators are both prximal operator for L1 norm
    prox_l1 first(2., 1.);
    prox_l1 second(2., 2.);
    Vector x;
    double data_x[3] = {0., 3., -1.}; 
    x.assign(data_x, data_x + 3);

    PeacemanRachfordSplitting<prox_l1, prox_l1> PRs(&x, first, second);
    // Initialize the Peaceman-Rachford operator
    Params params;
    params.step_size = 1.;
    PRs.update_params(&params);

    PRs(1);
    EXPECT_DOUBLE_EQ(1.764, x[1]);

}