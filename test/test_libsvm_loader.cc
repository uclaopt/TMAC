#include <iostream>
#include <fstream>
#include "matrices.h"
#include "util.h"
using namespace std;



int main(int argc, char *argv[]) {
  
  Params params;
  string data_file_name;
  string label_file_name;

  SpMat A;
  Vector b;
  
  parse_input_argv_libsvm(&params, argc, argv, data_file_name);
  loadLibSVM(A, b, data_file_name);
  cout << A.rows() << " " << A.cols() << endl;

}
