#include "gtest/gtest.h"
#include "matrices.h"
#include "algebra.h"
#include "constants.h"
#include "algebra_namespace_switcher.h"

/******************************************
 *  Unit test for linear algebra functions
 ******************************************/


// test the shrinkage function
TEST(ShrinkTest, DataSetI) {

  EXPECT_DOUBLE_EQ(9., shrink(10., 1.));
  EXPECT_DOUBLE_EQ(0., shrink(0.5, 1.));
  EXPECT_DOUBLE_EQ(0, shrink(-0.5, 1.));
  EXPECT_DOUBLE_EQ(-1.5, shrink(-2.5, 1.));
}


// test the norm function
TEST(NormTest,DataSetI) {
  
  Vector tmp_vector;
  double tmp_data[] = {177.6, 0, 1.4, -27.3};
  tmp_vector.assign(tmp_data, tmp_data + 4);
  EXPECT_DOUBLE_EQ(3, norm(tmp_vector, 0));
  EXPECT_DOUBLE_EQ(177.6 + 0 + 1.4 + 27.3, norm(tmp_vector, 1));
  EXPECT_DOUBLE_EQ(sqrt(177.6 * 177.6 + 0 * 0 + 1.4 * 1.4 + (-27.3) * (-27.3)), norm(tmp_vector, 2));
}


// test the calculate_column_norm : Non-sparse Matrix
TEST(Calculate_Column_Norm_Test, Regular) {
  
  Matrix tmp_matrix(3, 3, 0);
  Vector nrm(3, 0);
  double tmp_data[] = {5., -1., -6., -3.1, 4.2, 5.7, 6.4, 7.0, -7.9};
  double result[] = {0., 0., 0.};
  tmp_matrix.assign(tmp_data, tmp_data + 3 * 3);
  for (int i = 0; i < 3; i++)
    result[i] = sqrt(tmp_matrix(0,i)*tmp_matrix(0,i)+tmp_matrix(1,i)*tmp_matrix(1,i)+tmp_matrix(2,i)*tmp_matrix(2,i));
  calculate_column_norm(tmp_matrix, nrm);
  for (int i = 0; i < 3; i++)
    EXPECT_DOUBLE_EQ(result[i], nrm[i]);
}

TEST(TransposeTest, Dense) {
  
  Matrix A(1, 2, 0);
  double tmp_data[] = {1, 2};
  A.assign(tmp_data, tmp_data + 1 * 2);
  Matrix At(2, 1, 0);
  transpose(A, At);
  EXPECT_DOUBLE_EQ(At(0, 0), 1);
  EXPECT_DOUBLE_EQ(At(1, 0), 2);  

}


TEST(TransposeTest, Sparse) {
  
  SpMat A(2, 3);
  A.insert(0, 0) = 1.;
  A.insert(0, 1) = 2.;  
  
  SpMat At(3, 2);
  transpose(A, At);
  EXPECT_DOUBLE_EQ(At(0, 0), 1.);
  EXPECT_DOUBLE_EQ(At(1, 0), 2.);
}


// test the calculate_column_norm : Sparse Matrix
TEST(Calculate_Column_Norm_Test, Sparse) {
  
  /*
    1.,   -3.,  0,   0,   1., 0
    0,    0,    0,   0,   9., 0
    0,    0,    0,   0,   0,  0
   -0.5,  0,    0 ,  0,   0,  0
    0,    0,    0,   -4., 0,  0
    0,    0,    7.,  0,   0,  0
   */
  SpMat A(6, 6);
  Vector result(6, 0);
  double compare[6] ={0.};
  A.insert(0, 0) = 1.;
  A.insert(0, 1) = -3.;
  A.insert(0, 4) = 1.;
  A.insert(1, 4) = 9.;
  A.insert(3, 0) = -.5;
  A.insert(4, 3) = -4.;
  A.insert(5, 2) = 7.;
  A.makeCompressed();
  calculate_column_norm(A, result);
  EXPECT_DOUBLE_EQ(sqrt(1 * 1 + (-.5) * (-.5)), result[0]);
  EXPECT_DOUBLE_EQ(sqrt(3. * 3.), result[1]);
  EXPECT_DOUBLE_EQ(sqrt(7. * 7.), result[2]);
  EXPECT_DOUBLE_EQ(sqrt(4. * 4.), result[3]);
  EXPECT_DOUBLE_EQ(sqrt(1. * 1. + 9. * 9.), result[4]);
  EXPECT_DOUBLE_EQ(sqrt(0.), result[5]);
}


// test the a = a - b and a = a + b
TEST(SUB_ADD, DataSetI) {
  
  Vector a2, b;
  double tmp_a[] = {177.6, 0, 1.4, -27.3};
  double tmp_b[] = {13., -10., 4., 64.};
  a2.assign(tmp_a, tmp_a + 4);
  b.assign(tmp_b, tmp_b + 4);
  add(a2, b);
  for (int i=0; i < 4; i++)
  {
    EXPECT_DOUBLE_EQ(a2[i], tmp_a[i] + tmp_b[i]);
  }
}


TEST(SUM, DATA) {
  Vector a;
  double data_a[3] = {1., 3., -7.5};
  double result = 0.;
  a.assign(data_a, data_a + 3);
  result = sum(a);
  EXPECT_DOUBLE_EQ( 1. + 3. + (-7.5), result);
  
}


// test a = a - scalar * A(row, :) : Dense_Dense
TEST(SUB, Dense_Dense) {
  
  Matrix tmp_matrix(3, 3, 0);
  Vector a;
  double v_data[] = {1., 2., 3.};
  double m_data[] = {5., -1., -6., -3.1, 4.2, 5.7, 6.4, 7.0, -7.9};
  double result[3]= {0., 0., 0.};
  tmp_matrix.assign(m_data, m_data + 3 * 3);
  //1st row
  a.assign(v_data, v_data + 3);
  add(a, tmp_matrix, 0, -2.5);
  result[0]= -11.5, result[1]= 4.5, result[2]= 18.;
  for (int i = 0; i < 3; i++)
    EXPECT_DOUBLE_EQ(result[i], a[i]);
  //2nd row
  a.assign(v_data, v_data + 3);
  add(a, tmp_matrix, 1, -7.);
  result[0]= 22.7, result[1]= -27.4, result[2]= -36.9;
  for (int i = 0; i < 3; i++)
    EXPECT_DOUBLE_EQ(result[i], a[i]);
  //3rd row
  a.assign(v_data, v_data + 3);
  add(a, tmp_matrix, 2, 3.);
  result[0]= 20.2, result[1]= 23., result[2]= -20.7;
  for (int i = 0; i < 3; i++)
    EXPECT_DOUBLE_EQ(result[i], a[i]);
  
//test the pointer version
  Vector* pa = &a;
  Matrix* pm = &tmp_matrix;
  //1st row
  a.assign(v_data, v_data + 3);
  add(pa, pm, 0, -2.5);
  result[0]= -11.5, result[1]= 4.5, result[2]= 18.;
  for (int i = 0; i < 3; i++)
    EXPECT_DOUBLE_EQ(result[i],(*pa)[i]);
  //2nd row
  a.assign(v_data, v_data + 3);
  add(a, tmp_matrix, 1, -7.);
  result[0]= 22.7, result[1]= -27.4, result[2]= -36.9;
  for (int i = 0; i < 3; i++)
    EXPECT_DOUBLE_EQ(result[i], (*pa)[i]);
  //3rd row
  a.assign(v_data, v_data + 3);
  add(a, tmp_matrix, 2, 3.);
  result[0]= 20.2, result[1]= 23., result[2]= -20.7;
  for (int i = 0; i < 3; i++)
    EXPECT_DOUBLE_EQ(result[i], (*pa)[i]);
}


// test a = a - scalar * A(row, :) : Dense_Sparse
TEST(SUB, Dense_Sparse) {
  /*
   1., -3., 0, 0, 1., 0
   0, 0, 0, 0, 9., 0
   0, 0, 0, 0, 0, 0
   -0.5, 0, 0, 0, 0, 0
   0, 0, 0, -4., 0, 0
   0, 0, 7., 0, 0, 0
   */
  SpMat tmp_matrix(6, 6);
  tmp_matrix.insert(0, 0) = 1.;
  tmp_matrix.insert(0, 1) = -3.;
  tmp_matrix.insert(0, 4) = 1.;
  tmp_matrix.insert(1, 4) = 9.;
  tmp_matrix.insert(3, 0) = -.5;
  tmp_matrix.insert(4, 3) = -4.;
  tmp_matrix.insert(5, 2) = 7.;
  tmp_matrix.makeCompressed();
  Vector a;
  double v_data[] = {1., 2., 2., -3., -1., 1.};
  double result[6] = {0., 0., 0., 0., 0., 0.};
  //1st row
  a.assign(v_data, v_data + 6);
  add(a, tmp_matrix, 0, -2.5);
  result[0]= -1.5, result[1]= 9.5, result[2]= 2.,result[3]= -3., result[4]= -3.5, result[5]= 1.;
  for (int i = 0; i < 6; i++)
    EXPECT_DOUBLE_EQ(result[i], a[i]);
  //3rd row
  a.assign(v_data, v_data + 6);
  add(a, tmp_matrix, 2, -7.);
  result[0]= 1., result[1]= 2., result[2]= 2.,result[3]= -3., result[4]= -1., result[5]= 1.;
  for (int i = 0; i < 6; i++)
    EXPECT_DOUBLE_EQ(result[i], a[i]);
  //5th row
  a.assign(v_data, v_data + 6);
  add(a, tmp_matrix, 4, 3.);
  result[0]= 1., result[1]= 2., result[2]= 2.,result[3]= -15., result[4]= -1., result[5]= 1.;
  for (int i = 0; i < 6; i++)
    EXPECT_DOUBLE_EQ(result[i], a[i]);

//test the pointer version
  Vector* pa = &a;
  SpMat* pm =&tmp_matrix;
  //1st row
  a.assign(v_data, v_data + 6);
  add(pa, pm, 0, -2.5);
  result[0]= -1.5, result[1]= 9.5, result[2]= 2.,result[3]= -3., result[4]= -3.5, result[5]= 1.;
  for (int i = 0; i < 6; i++)
    EXPECT_DOUBLE_EQ(result[i], (*pa)[i]);
  //3rd row
  a.assign(v_data, v_data + 6);
  add(pa, pm, 2, -7.);
  result[0]= 1., result[1]= 2., result[2]= 2.,result[3]= -3., result[4]= -1., result[5]= 1.;
  for (int i = 0; i < 6; i++)
    EXPECT_DOUBLE_EQ(result[i], (*pa)[i]);
  //5th row
  a.assign(v_data, v_data + 6);
  add(pa, pm, 4, 3.);
  result[0]= 1., result[1]= 2., result[2]= 2.,result[3]= -15., result[4]= -1., result[5]= 1.;
  for (int i = 0; i < 6; i++)
    EXPECT_DOUBLE_EQ(result[i], (*pa)[i]);
}


// test a = a - scalar * A(row, :) : Sparse_Density
TEST(SUB, Sparse_Dense) {

  SpVec a(6); //{-5., 0., 0., 1.5, 0., 0.};
  Matrix m(6, 6, 0);
  double m_data[] = {1., -3., 0, 0, 1., 0, 0, 0, 0, 0, 9., 0, 0, 0, 0, 0, 0, 0, -0.5, 0, \
    0, 0, 0, 0, 0, 0, 0, -4., 0, 0, 0, 0, 7., 0, 0, 0};
  double result[6] = {0., 0., 0., 0., 0., 0.};
  m.assign(m_data, m_data + 6 * 6);
  //1st row
  a.insert(0) = -5., a.insert(3) = 1.5;
  add(a, m, 0, -2.5);
  result[0]= -7.5, result[1]= 7.5, result[2]= 0.,result[3]= 1.5, result[4]= -2.5, result[5]= 0.;//{-7.5, 7.5, 0., 1.5, -2.5, 0.}
  for (int i = 0; i < 6; i++)
    EXPECT_DOUBLE_EQ(result[i], a.coeffRef(i));
  //3rd row
    a.coeffRef(0) = -5., a.coeffRef(1) = 0., a.coeffRef(2) = 0., \
    a.coeffRef(3) = 1.5, a.coeffRef(4) = 0., a.coeffRef(5) = 0.;
  add(a, m, 2, -7.);
  result[0]= -5., result[1]= 0., result[2]= 0.,result[3]= 1.5, result[4]= 0., result[5]= 0.;//-5., 0., 0., 1.5, 0., 0.
  for (int i = 0; i < 6; i++)
    EXPECT_DOUBLE_EQ(result[i], a.coeffRef(i));
  //5th row
  a.coeffRef(0) = -5., a.coeffRef(1) = 0., a.coeffRef(2) = 0., \
  a.coeffRef(3) = 1.5, a.coeffRef(4) = 0., a.coeffRef(5) = 0.;
  add(a, m, 4, 3.);
  result[0]= -5., result[1]= 0., result[2]= 0.,result[3]= -10.5, result[4]= 0., result[5]= 0.;//{-5., 0., 0., -10.5, 0., 0.}
  for (int i = 0; i < 6; i++)
    EXPECT_DOUBLE_EQ(result[i], a.coeffRef(i));
}


// test a = a - scalar * A(row, :) : Sparse_Sparse
TEST(SUB, Sparse_Sparse) {
  
  SpVec a(6); //{-5., 0., 0., 1.5, 0., 0.}
  /*
   1., -3., 0, 0, 1., 0
   0, 0, 0, 0, 9., 0
   0, 0, 0, 0, 0, 0
   -0.5, 0, 0, 0, 0, 0
   0, 0, 0, -4., 0, 0
   0, 0, 7., 0, 0, 0
   */
  SpMat tmp_matrix(6, 6);
  tmp_matrix.insert(0, 0) = 1.;
  tmp_matrix.insert(0, 1) = -3.;
  tmp_matrix.insert(0, 4) = 1.;
  tmp_matrix.insert(1, 4) = 9.;
  tmp_matrix.insert(3, 0) = -.5;
  tmp_matrix.insert(4, 3) = -4.;
  tmp_matrix.insert(5, 2) = 7.;
  tmp_matrix.makeCompressed();
  double result[6] = {0., 0., 0., 0., 0., 0.};
  //1st row
  a.insert(0) = -5., a.insert(3) = 1.5;
  add(a, tmp_matrix, 0, -2.5);
  result[0]= -7.5, result[1]= 7.5, result[2]= 0.,result[3]= 1.5, result[4]= -2.5, result[5]= 0.;//{-7.5, 7.5, 0., 1.5, -2.5, 0.}
  for (int i = 0; i < 6; i++)
    EXPECT_DOUBLE_EQ(result[i], a.coeffRef(i));
  //3rd row
  a.coeffRef(0) = -5., a.coeffRef(1) = 0., a.coeffRef(2) = 0., \
  a.coeffRef(3) = 1.5, a.coeffRef(4) = 0., a.coeffRef(5) = 0.;
  add(a, tmp_matrix, 2, -7.);
  result[0]= -5., result[1]= 0., result[2]= 0.,result[3]= 1.5, result[4]= 0., result[5]= 0.;//-5., 0., 0., 1.5, 0., 0.
  for (int i = 0; i < 6; i++)
    EXPECT_DOUBLE_EQ(result[i], a.coeffRef(i));
  //5th row
  a.coeffRef(0) = -5., a.coeffRef(1) = 0., a.coeffRef(2) = 0., \
  a.coeffRef(3) = 1.5, a.coeffRef(4) = 0., a.coeffRef(5) = 0.;
  add(a, tmp_matrix, 4, 3.);
  result[0]= -5., result[1]= 0., result[2]= 0.,result[3]= -10.5, result[4]= 0., result[5]= 0.;//{-5., 0., 0., -10.5, 0., 0.}
  for (int i = 0; i < 6; i++)
    EXPECT_DOUBLE_EQ(result[i], a.coeffRef(i));
}


// test a = a + lambda * b
TEST(ADD, Dense_Dense) {
  Vector a(3), b(3);
  double data_a[3] = {1., 3., -7.};
  double data_b[3] = {2., -5.5, 6.7};
  a.assign(data_a, data_a + 3), b.assign(data_b, data_b + 3);
  add(a, b, 9.);
  EXPECT_DOUBLE_EQ(19., a[0]);
  EXPECT_DOUBLE_EQ(-46.5, a[1]);
  EXPECT_DOUBLE_EQ(53.3, a[2]);
}


// test the scaling function
TEST(SCALING, Dense) {
  Vector a(3);
  double data_a[3] = {1., 3., -7.};
  a.assign(data_a, data_a + 3);
  scale(a, 5.);
  EXPECT_DOUBLE_EQ(5., a[0]);
  EXPECT_DOUBLE_EQ(15., a[1]);
  EXPECT_DOUBLE_EQ(-35., a[2]);
}


// test the addition of sparse vectors
TEST(ADD, Sparse_Sparse) {
  Vector a(6);
  SpVec b(6);
  a[0] =-.5, a[1] = 0., a[2] = 0., a[3] = 1.5, a[4] = 0., a[5] = 0.;//{-.5, 0., 0., 1.5, 0., 0.}
  b.insert(3) = 3., b.insert(5) = -2.; //{0., 0., 0., 3., 0., -2.};
  add(a, b);
  EXPECT_DOUBLE_EQ(-.5,a[0]);
  EXPECT_DOUBLE_EQ(0., a[1]);
  EXPECT_DOUBLE_EQ(0., a[2]);
  EXPECT_DOUBLE_EQ(4.5, a[3]);
  EXPECT_DOUBLE_EQ(0., a[4]);
  EXPECT_DOUBLE_EQ(-2., a[5]);
}


// test the inner product of two dense vectors
TEST(DOT, Vec_Vec) {
  Vector a(3, 0.), b(3, 0.);
  double data_a[3] = {1., 2., 3.};
  double data_b[3] = {7., -5., 1.5};
  double result = 0.;
  a.assign(data_a, data_a + 3), b.assign(data_b, data_b + 3);
  result = dot(a, b);
  EXPECT_DOUBLE_EQ(result, 1.5);
}


// test the inner product of a sparse matrix row and a vector
TEST(DOT, SpRow_Vec) {
  /*
   1., -3., 0, 0, 1., 0
   0, 0, 0, 0, 9., 0
   0, 0, 0, 0, 0, 0
   -0.5, 0, 0, 0, 0, 0
   0, 0, 0, -4., 0, 0
   0, 0, 7., 0, 0, 0
   */
  SpMat tmp_matrix(6, 6);
  tmp_matrix.insert(0, 0) = 1.;
  tmp_matrix.insert(0, 1) = -3.;
  tmp_matrix.insert(0, 4) = 1.;
  tmp_matrix.insert(1, 4) = 9.;
  tmp_matrix.insert(3, 0) = -.5;
  tmp_matrix.insert(4, 3) = -4.;
  tmp_matrix.insert(5, 2) = 7.;
  tmp_matrix.makeCompressed();
  Vector v(6, 0.);
  double result = 0.;
  double data_v[6] = {1., -1., 2.5, 0., -3., 0.5};
  v.assign(data_v, data_v+6);
  
  result = dot(tmp_matrix, v, 0);
  EXPECT_DOUBLE_EQ(1., result);
  
  result = dot(tmp_matrix, v, 1);
  EXPECT_DOUBLE_EQ(-27., result);
  
  result = dot(tmp_matrix, v, 2);
  EXPECT_DOUBLE_EQ(0., result);
  
  result = dot(tmp_matrix, v, 3);
  EXPECT_DOUBLE_EQ(-0.5, result);
  
  result = dot(tmp_matrix, v, 4);
  EXPECT_DOUBLE_EQ(0., result);
  
  result = dot(tmp_matrix, v, 5);
  EXPECT_DOUBLE_EQ(17.5, result);
  
  // test the pointer version
  SpMat* mp = &tmp_matrix;
  Vector* vp = & v;
  result = dot((*mp), *(vp), 0);
  EXPECT_DOUBLE_EQ(1., result);
  
  result = dot((*mp), *(vp), 1);
  EXPECT_DOUBLE_EQ(-27., result);
  
  result = dot((*mp), *(vp), 2);
  EXPECT_DOUBLE_EQ(0., result);
  
  result = dot((*mp), *(vp), 3);
  EXPECT_DOUBLE_EQ(-0.5, result);
  
  result = dot((*mp), *(vp), 4);
  EXPECT_DOUBLE_EQ(0., result);
  
  result = dot((*mp), *(vp), 5);
  EXPECT_DOUBLE_EQ(17.5, result);
}


// test the inner product of a matrix row and a vector
TEST(DOT, Row_Vec) {
  Matrix m(6, 6, 0);
  double m_data[] = {1., -3., 0, 0, 1., 0, 0, 0, 0, 0, 9., 0, 0, 0, 0, 0, 0, 0, -0.5, 0, \
    0, 0, 0, 0, 0, 0, 0, -4., 0, 0, 0, 0, 7., 0, 0, 0};
  m.assign(m_data, m_data + 6 * 6);
  Vector v(6, 0.);
  double data_v[6] = {1., -1., 2.5, 0., -3., 0.5};
  v.assign(data_v, data_v+6);
  
  double result = 0.;
  result = dot(m, v, 0);
  EXPECT_DOUBLE_EQ(1., result);
  
  result = dot(m, v, 1);
  EXPECT_DOUBLE_EQ(-27., result);
  
  result = dot(m, v, 2);
  EXPECT_DOUBLE_EQ(0., result);
  
  result = dot(m, v, 3);
  EXPECT_DOUBLE_EQ(-0.5, result);
  
  result = dot(m, v, 4);
  EXPECT_DOUBLE_EQ(0., result);
  
  result = dot(m, v, 5);
  EXPECT_DOUBLE_EQ(17.5, result);
  
  // test the pointer version
  Matrix* mp = &m;
  Vector* vp = & v;
  result = dot((*mp), *(vp), 0);
  EXPECT_DOUBLE_EQ(1., result);
  
  result = dot((*mp), *(vp), 1);
  EXPECT_DOUBLE_EQ(-27., result);
  
  result = dot((*mp), *(vp), 2);
  EXPECT_DOUBLE_EQ(0., result);
  
  result = dot((*mp), *(vp), 3);
  EXPECT_DOUBLE_EQ(-0.5, result);
  
  result = dot((*mp), *(vp), 4);
  EXPECT_DOUBLE_EQ(0., result);
  
  result = dot((*mp), *(vp), 5);
  EXPECT_DOUBLE_EQ(17.5, result);
}

// test the multiplication of the transpose of a sparse matrix and a vector
TEST(TransMultiplyTest, SpMatVec) {

  int rows = 3, cols = 2;
  SpMat A(rows, cols);
  A.coeffRef(0, 0) = 1;
  A.coeffRef(2, 1) = 2;  
  Vector x(rows, 1.);
  Vector Atx(cols, 0.);
  trans_multiply(A, x, Atx);
  EXPECT_DOUBLE_EQ(1., Atx[0]);
  EXPECT_DOUBLE_EQ(2., Atx[1]);  

}


// test the multiplication of the transpose of a matrix and a vector
TEST(TransMultiplyTest, MatrixVec) {
  Matrix m(6, 6, 0);
  double m_data[] = {1., -3., 0, 0, 1., 0, 0, 0, 0, 0, 9., 0, 0, 0, 0, 0, 0, 0, -0.5, 0, \
    0, 0, 0, 0, 0, 0, 0, -4., 0, 0, 0, 0, 7., 0, 0, 0};
  m.assign(m_data, m_data + 6 * 6);
  Vector v(6, 0.), Atx(6, 0.);
  double data_v[6] = {1., -1., 2.5, 0., -3., 0.5};
  v.assign(data_v, data_v+6);
  trans_multiply(m, v, Atx);
  // { 1. ,  -3. ,   3.5,  12. ,  -8. ,   0. }
  EXPECT_DOUBLE_EQ(1.,  Atx[0]);
  EXPECT_DOUBLE_EQ(-3., Atx[1]);
  EXPECT_DOUBLE_EQ(3.5, Atx[2]);
  EXPECT_DOUBLE_EQ(12., Atx[3]);
  EXPECT_DOUBLE_EQ(-8., Atx[4]);
  EXPECT_DOUBLE_EQ(0.,  Atx[5]);
}


// test the multiplication of a sparse matrix and a vector
TEST(MultiplyTest, SpMatVec) {

  int rows = 3, cols = 2;
  SpMat A(rows, cols);
  A.coeffRef(0, 0) = 1;
  A.coeffRef(2, 1) = 2;  
  Vector x(cols, 1.);
  Vector Ax(rows, 0.);
  multiply(A, x, Ax);
  EXPECT_DOUBLE_EQ(1., Ax[0]);
  EXPECT_DOUBLE_EQ(0., Ax[1]);
  EXPECT_DOUBLE_EQ(2., Ax[2]);  

}


// test the multiplication of a matrix and a vector
TEST(MultiplyTest, MatrixVec) {
  Matrix m(6, 6, 0);
  double m_data[] = {1., -3., 0, 0, 1., 0, 0, 0, 0, 0, 9., 0, 0, 0, 0, 0, 0, 0, -0.5, 0, \
    0, 0, 0, 0, 0, 0, 0, -4., 0, 0, 0, 0, 7., 0, 0, 0};
  m.assign(m_data, m_data + 6 * 6);
  Vector v(6, 0.), Ax(6, 0.);
  double data_v[6] = {1., -1., 2.5, 0., -3., 0.5};
  v.assign(data_v, data_v+6);
  multiply(m, v, Ax);
  // { 1. , -27. ,   0. ,  -0.5,   0. ,  17.5}
  EXPECT_DOUBLE_EQ(1.,   Ax[0]);
  EXPECT_DOUBLE_EQ(-27., Ax[1]);
  EXPECT_DOUBLE_EQ(0.,   Ax[2]);
  EXPECT_DOUBLE_EQ(-0.5, Ax[3]);
  EXPECT_DOUBLE_EQ(0.,   Ax[4]);
  EXPECT_DOUBLE_EQ(17.5, Ax[5]);
}

