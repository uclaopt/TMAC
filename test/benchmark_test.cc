#include <iostream>
#include "matrices.h"
#include "algebra.h"
#include "parameters.h"
#include "util.h"
#include <thread>
using namespace std;

void dot_ds_vec_ds_vec_test(Vector* x, Vector* y, int max_itr, int num_thread, int type) {

  int rows = x->size();
  for (int i = 0; i < max_itr / num_thread; ++i) {
    if (type == 0) {
      MyAlgebra::dot(*x, *y);
    } else {
      BLASAlgebra::dot(*x, *y);
    }
  }
  return;
  
}


void logistic_test(int max_itr, int num_thread) {

  double x = 88.88;
  for (int i = 0; i < max_itr / num_thread; ++i) {
    // auto x = (double)rand() / (double)RAND_MAX;      
    log(1 + exp(-x));
  }

  return;
  
}


void dot_ds_mtx_ds_vec_test(Matrix* A, Vector* x, int max_itr, int num_thread, int type) {

  int rows = A->rows();
  int cols = A->cols();
  for (int i = 0; i < max_itr / num_thread; ++i) {
    for (int j = 0; j < rows; j++) {
      // int idx = rand() % rows;
      if (type == 0) {
        MyAlgebra::dot(A, x, j);
      } else {
        BLASAlgebra::dot(A, x, j);        
      }
    }
  }
  return;
}



void dot_sp_mtx_ds_vec_test(SpMat* A, Vector* x, int max_itr, int num_thread, int type) {

  int rows = A->rows();
  int cols = A->cols();
  for (int i = 0; i < max_itr / num_thread; ++i) {
    for (int j = 0; j < rows; j++) {
      if (type == 0) {
        MyAlgebra::dot(A, x, j);
      } else {
        BLASAlgebra::dot(A, x, j);        
      }
    }
  }
  return;
}





void scale_test(Vector* x, int max_itr, int num_thread, int type) {

  int rows = x->size();
  for (int i = 0; i < max_itr / num_thread; ++i) {
    if (type == 0) {
      MyAlgebra::scale(*x, 0.88);
    } else {
      BLASAlgebra::scale(*x, 0.88);
    }
  }
  return;
}


void logistic_speedup_test(Params params) {

  vector<std::thread> mythreads;
  int num_tests = 1;
  int max_itrs = params.max_itrs;
  int max_num_thraeds = params.total_num_threads;
  cout << "A" << " = [... " << endl;
  double blas_time = 0., my_time = 0.;

  for (int test = 1; test <= num_tests; test++) {

    for (int runs = 0; runs < 10; runs++) {      

      for (int num_threads = 1; num_threads <= max_num_thraeds; num_threads*=2) {

        double start_time = get_wall_time();    
        for (int i = 0; i < num_threads; i++) {
          mythreads.push_back(std::thread(logistic_test, max_itrs, num_threads));
        }
        for (size_t i = 0; i < num_threads; i++) {
          mythreads[i].join();
        }
        double end_time = get_wall_time();
        my_time = end_time - start_time;
        mythreads.clear();        
        
        cout << setw(15) << setprecision(2) << scientific << num_threads;
        cout << setw(15) << setprecision(2) << my_time << endl;
      }
    }
  }
  cout << "];" << endl;
  
}

void dot_prod_speedup_test(Params params) {

  vector<std::thread> mythreads;
  int num_tests = 4;
  int max_itrs = params.max_itrs;
  int max_num_thraeds = params.total_num_threads;
  int len = 10000;
  cout << "A" << " = [... " << endl;
  double blas_time = 0., my_time = 0.;

  for (int test = 1; test <= num_tests; test++) {

    for (int runs = 0; runs < 10; runs++) {      
      Vector x(len), y(len);
      for (int num_threads = 1; num_threads <= max_num_thraeds; num_threads*=2) {
        for (int type = 0; type <= 1; type++) {
          double start_time = get_wall_time();    
          for (int i = 0; i < num_threads; i++) {
            mythreads.push_back(std::thread(dot_ds_vec_ds_vec_test, &x, &y, max_itrs, num_threads, type));
          }
          for (size_t i = 0; i < num_threads; i++) {
            mythreads[i].join();
          }
          double end_time = get_wall_time();
          if (type == 0) {
            my_time = end_time - start_time;
          } else {
            blas_time = end_time - start_time;
          }
          mythreads.clear();        
        }
      
        cout << setw(15) << setprecision(2) << scientific << num_threads;
        cout << setw(15) << setprecision(2) << my_time;
        cout << setw(15) << setprecision(2) << blas_time;
        cout << setw(15) << setprecision(2) << len << endl;
      }

    }
    len *= 10;
  }
  cout << "];" << endl;
  
}


void dot_ds_mtx_ds_vec_speedup_test(Params params) {

  vector<std::thread> mythreads;
  int num_tests = 4;
  int max_itrs = params.max_itrs;
  int max_num_thraeds = params.total_num_threads;
  int len = 1000;
  cout << "ds_mtx_ds_vec" << " = [... " << endl;
  double blas_time = 0., my_time = 0.;
  
  for (int test = 1; test <= num_tests; test++) {

    Matrix A(len/100, len, 1.);
    Vector x(len, 0.);

    for (int runs = 0; runs < 10; runs++) {      
      for (int num_threads = 1; num_threads <= max_num_thraeds; num_threads*=2) {
        for (int type = 0; type <= 1; type++) {
          double start_time = get_wall_time();    
          for (int i = 0; i < num_threads; i++) {
            mythreads.push_back(std::thread(dot_ds_mtx_ds_vec_test, &A, &x, max_itrs, num_threads, type));
          }
          for (size_t i = 0; i < num_threads; i++) {
            mythreads[i].join();
          }
          double end_time = get_wall_time();
          if (type == 0) {
            my_time = end_time - start_time;
          } else {
            blas_time = end_time - start_time;
          }
          mythreads.clear();        
        }
        
        cout << setw(15) << setprecision(2) << scientific << num_threads;
        cout << setw(15) << setprecision(2) << my_time;
        cout << setw(15) << setprecision(2) << blas_time;
        cout << setw(15) << setprecision(2) << len << endl;

      }

    }
    len *= 5;
  }

  cout << "];" << endl;
  
}

void dot_sp_mtx_ds_vec_speedup_test(Params params) {

  vector<std::thread> mythreads;
  int num_tests = 4;
  int max_itrs = params.max_itrs;
  int max_num_thraeds = params.total_num_threads;
  int len = 1000;
  cout << "sp_mtx_ds_vec" << " = [... " << endl;
  double blas_time = 0., my_time = 0.;
  
  for (int test = 1; test <= num_tests; test++) {
    int rows=len/100;
    int cols=len;

    srand(0);
    std::vector<Eigen::Triplet<double> > tripletList;
    for(int i=0;i<rows;++i) {
      for(int j=0;j<cols;++j) {
	auto v_ij = (double)rand() / (double)RAND_MAX;
	if(v_ij < 0.05) { // about 5% of nonzeros
	  tripletList.push_back(Eigen::Triplet<double>(i,j,v_ij));      //if larger than treshold, insert it
	}
      }
    }
    
    SpMat A(rows, cols);
    A.setFromTriplets(tripletList.begin(), tripletList.end());   //create the matrix
    Vector x(len, 0.);

    for (int runs = 0; runs < 10; runs++) {      
      for (int num_threads = 1; num_threads <= max_num_thraeds; num_threads*=2) {
        for (int type = 0; type <= 1; type++) {
          double start_time = get_wall_time();    
          for (int i = 0; i < num_threads; i++) {
            mythreads.push_back(std::thread(dot_sp_mtx_ds_vec_test, &A, &x, max_itrs, num_threads, type));
          }
          for (size_t i = 0; i < num_threads; i++) {
            mythreads[i].join();
          }
          double end_time = get_wall_time();
          if (type == 0) {
            my_time = end_time - start_time;
          } else {
            blas_time = end_time - start_time;
          }
          mythreads.clear();        
        }
        
        cout << setw(15) << setprecision(2) << scientific << num_threads;
        cout << setw(15) << setprecision(2) << my_time;
        cout << setw(15) << setprecision(2) << blas_time;
        cout << setw(15) << setprecision(2) << len << endl;

      }

    }
    len *= 5;
  }
  cout << "];" << endl;
  
}



void scale_speedup_test(Params params) {

  vector<std::thread> mythreads;
  int num_tests = 4;
  int max_itrs = params.max_itrs;
  int max_num_thraeds = params.total_num_threads;
  int len = 1000;
  cout << "A" << " = [... " << endl;
  double blas_time = 0., my_time = 0.;
  
  len = 1000;
  for (int test = 1; test <= num_tests; test++) {
    len *= 10;
    for (int runs = 0; runs < 10; runs++) {      
      Vector x(len);
      for (int num_threads = 1; num_threads <= max_num_thraeds; num_threads*=2) {
        for (int type = 0; type <= 1; type++) {
          double start_time = get_wall_time();    
          for (int i = 0; i < num_threads; i++) {
            mythreads.push_back(std::thread(scale_test, &x, max_itrs, num_threads, type));
          }
          for (size_t i = 0; i < num_threads; i++) {
            mythreads[i].join();
          }
          double end_time = get_wall_time();
          if (type == 0) {
            my_time = end_time - start_time;
          } else {
            blas_time = end_time - start_time;
          }
          mythreads.clear();        
        }
      
        cout << setw(15) << setprecision(2) << scientific << num_threads;
        cout << setw(15) << setprecision(2) << my_time;
        cout << setw(15) << setprecision(2) << blas_time;
        cout << setw(15) << setprecision(2) << len << endl;
      }
    }
  }
  cout << "];" << endl;
}






int main( int argc, char** argv ){

  // Step 0: Define the parameters and input file names  
  Params params;
  double lambda = 1.;
  
  // Step 1. Parse the input argument
  parse_input_argv_demo(&params, argc, argv);
  // scale_speedup_test(params);
  // dot_prod_speedup_test(params);
  dot_ds_mtx_ds_vec_speedup_test(params);

  dot_sp_mtx_ds_vec_speedup_test(params);  
  // logistic_speedup_test(params);
  return 0;
}
