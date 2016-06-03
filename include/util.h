/*
 * Helper functions header file
 *
 *
 */

#ifndef TMAC_INCLUDE_UTIL_H
#define TMAC_INCLUDE_UTIL_H
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include "constants.h"
#include "matrices.h"
#include "parameters.h"
#include "algebra.h"
#include <string>
#include <vector>
#include "algebra_namespace_switcher.h"



// function for printing help message;
void exit_with_help(string, string, bool);



void set_default_settings(Params* para);


void print_result_info();


void print_parameters(Params& para);

void parse_input_argv_demo(Params* para,
                           int argc,
                           char *argv[]);


void parse_input_argv_mm(Params* para,
                         int argc,
                         char *argv[],
                         std::string& data_file_name,
                         std::string& label_file_name);

void parse_input_argv_mm(Params* para,
                         int argc,
                         char *argv[],
                         std::string& data_file_name,
                         std::string& label_file_name,
                         double& lambda);


void parse_input_argv_libsvm(Params* para,
                             int argc,
                             char *argv[],
                             std::string& data_file_name);


void parse_input_argv_libsvm(Params* para,
                             int argc,
                             char *argv[],
                             std::string& data_file_name,
                             double& lambda);



/*************************************************
 *
 * Logistic regression related utility functions
 *
 *************************************************/


double log_loss_gradient_at_idx ( Matrix& A, Vector& b, Vector& Atx, int idx );


double log_loss_gradient_at_idx ( SpMat& A, Vector& b, Vector& Atx, int idx );


/********************************************************************
 *  calculates l2 regularized objective
 *  Input:
 *     A:      the data matrix with size num_features x num_samples
 *     (type T, it can be sparse matrix SpMat or Matrix)
 *     b:      the observation labels for each sample
 *     (Vector)
 *     x:      the unknowns
 *     (Vector, size is the number of features)
 *     Atx:    A'*x, which is stored in shared memory for efficient
 *             computation
 *     lambda: regularization parameter
 *      
 *  Output: objective value
 *     (double)
 *******************************************************************/
double l2_log_loss_objective (Vector& b, Vector& x, Vector& Atx, double lambda);


/********************************************************************
 *  calculates l2 regularized objective
 *  Input:
 *     A:      the data matrix with size num_features x num_samples
 *     (type T, it can be sparse matrix SpMat or Matrix)
 *     b:      the observation labels for each sample
 *     (Vector)
 *     x:      the unknowns
 *     (Vector, size is the number of features)
 *     Atx:    A'*x, which is stored in shared memory for efficient
 *             computation
 *     lambda: regularization parameter
 *      
 *  Output: objective value
 *     (double)
 *******************************************************************/

double l1_log_loss_objective (Vector& b, Vector& x, Vector& Atx, double lambda);
double log_loss(Vector& x, Vector& Atx, Vector& b);
double square_loss(Vector& x, Vector& Atx, Vector& b);
double square_hinge_loss(Vector& x, Vector& Atx, Vector& b);
double huber_loss(Vector& x, Vector& Atx, Vector& b, double delta);
double quad_func(Vector& x, Vector& Qx, Vector& c, double d);
double huber_norm(Vector& x, double delta);


/******************************
 * Load data in libsvm format
 ******************************/
template<typename SparseMatrixType>
bool loadLibSVM(SparseMatrixType& mat, Vector& v, const std::string& filename);


/******************************
 * Load data in MATLAB format (.m file)
 ******************************/
template<typename SparseMatrixType>
bool loadMatlabSparse(SparseMatrixType& mat, const std::string& filename);

double get_wall_time();
double get_cpu_time();

#endif
