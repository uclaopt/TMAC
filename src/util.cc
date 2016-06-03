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
#include "util.h"

void exit_with_help(string app_name, string data_type = "matrix_market", bool use_regularization = false) {

  std::cout << "--------------------------------------------------------------" << std::endl;
  std::cout << "The usage for " << app_name << " is: " << std::endl;
  std::cout << "--------------------------------------------------------------" << std::endl;  
  std::cout << app_name << " [options] " << std::endl;

  if (data_type == "matrix_market") {
    std::cout << " -data      <data file with matrix market format, size: features x samples> \n";
    std::cout << " -label     <label file with matrix market format> \n";
    std::cout << " -nthread   <total number of threads, default is set to 2> \n";
    std::cout << " -epoch     <total number of epoch, default is set to 10> \n";
    if (use_regularization) {
      std::cout << " -lambda    <regularization parameter, default is set to 1> \n";
    }
  } else if (data_type == "libsvm") {
    std::cout << " -data      <data file with libsvm format> \n";
    std::cout << " -nthread   <total number of threads, default is set to 2> \n";
    std::cout << " -epoch     <total number of epoch, default is set to 10> \n";
    if (use_regularization) {
      std::cout << " -lambda    <regularization parameter, default is set to 1.> \n";
    }
  } else if (data_type == "no_data") {
    std::cout << " -epoch            <total number of epoch, default is set to 10> \n";
    std::cout << " -nthread          <total number of threads, default is set to 2> \n";    
    std::cout << " -problem_size     <the size of the problem, default is set to 0> \n";    
    if (use_regularization) {
    std::cout << " -lambda           <regularization parameter, default is set to 1.> \n";
    }
  }
  std::cout << "--------------------------------------------------------------" << std::endl;  
  abort();
}


void set_default_settings(Params* para) {
  para->max_itrs = MAX_ITRS;
  para->use_controller = USE_CONTROLLER;
  para->total_num_threads = TOTAL_NUM_THREADS;
}


void print_result_info() {


}


void print_parameters(Params& para) {
  cout << "Parameter settings:" << endl; 
  cout << "---------------------------------" << endl;
  cout.width(28);
  cout << left << "Problem size: " << right << setw(3) << para.problem_size << endl;
  cout.width(28);  
  cout << left << "TMAC step size: " << right << para.tmac_step_size << endl;
  cout.width(28);  
  cout << left << "Operator step size: " << right << para.step_size << endl;
  cout.width(28);  
  cout << left << "Use controller: " << right << ((para.use_controller) ? "true" : "false") << endl;
  cout << "---------------------------------" << endl;  
}


void parse_input_argv_mm(Params* para,
                         int argc,
                         char *argv[],
                         std::string& data_file_name,
                         std::string& label_file_name) {


  // set the values in para to default parameters
  set_default_settings(para);

  if (argc < 2) {
    exit_with_help(argv[0], "matrix_market", false);
  }
  
  for ( int i = 1; i < argc; ++i ) {
    if ( argv[i][0] != '-') {
      break;
    }
    if ( ++i >= argc ) {
      exit_with_help(argv[0], "matrix_market", false);
    }
    else if ( std::string ( argv[i-1] ) == "-epoch" ) {
      para->max_itrs = atoi ( argv[i] );
    }
    else if ( std::string ( argv[i-1] ) == "-data" ){
      data_file_name = std::string ( argv[i] );
    }
    else if ( std::string ( argv[i-1] ) == "-label" ) {
      label_file_name = std::string ( argv[i] );
    }
    else if ( std::string ( argv[i-1] ) == "-nthread" ) {
      para->total_num_threads = atoi ( argv[i] ) ;
    }
    else {
      exit_with_help(argv[0], "matrix_market", false);      
    }
  }
  return;
}



void parse_input_argv_mm(Params* para,
                         int argc,
                         char *argv[],
                         std::string& data_file_name,
                         std::string& label_file_name,
                         double& lambda) {

  // set the values in para to default parameters
  set_default_settings(para);
  
  if (argc < 2) {
    exit_with_help(argv[0], "matrix_market", true);
  }
  
  for ( int i = 1; i < argc; ++i ) {
    if ( argv[i][0] != '-' ) {
      break;
    }
    if ( ++i >= argc ) {
      exit_with_help(argv[0], "matrix_market", true);      
    }
    else if ( std::string ( argv[i-1] ) == "-epoch" ) {
      para->max_itrs = atoi ( argv[i] );
    }
    else if ( std::string ( argv[i-1] ) == "-data" ){
      data_file_name = std::string ( argv[i] );
    }
    else if ( std::string ( argv[i-1] ) == "-label" ) {
      label_file_name = std::string ( argv[i] );
    }
    else if ( std::string ( argv[i-1] ) == "-nthread" ) {
      para->total_num_threads = atoi ( argv[i] ) ;
    }
    else if ( std::string ( argv[i-1] ) == "-lambda" ) {
      lambda = atof ( argv[i] ) ;
    }
    else {
      exit_with_help(argv[0], "matrix_market", true);
    }
  }
  return;
}



void parse_input_argv_libsvm(Params* para,
                             int argc,
                             char *argv[],
                             std::string& data_file_name) {

  // set the values in para to default parameters
  set_default_settings(para);

  if (argc < 2) {
    exit_with_help(argv[0], "libsvm", false);
  }
  
  for ( int i = 1; i < argc; ++i ) {
    if ( argv[i][0] != '-') {
      break;
    }
    if ( ++i >= argc ) {
      exit_with_help(argv[0], "libsvm", false);
    }
    else if ( std::string ( argv[i-1] ) == "-epoch" ) {
      para->max_itrs = atoi ( argv[i] );
    }
    else if ( std::string ( argv[i-1] ) == "-data" ){
      data_file_name = std::string ( argv[i] );
    }
    else if ( std::string ( argv[i-1] ) == "-nthread" ) {
      para->total_num_threads = atoi ( argv[i] ) ;
    }
    else {
      exit_with_help(argv[0], "libsvm", false);      
    }
  }
  return;
}


void parse_input_argv_libsvm(Params* para,
                             int argc,
                             char *argv[],
                             std::string& data_file_name,
                             double& lambda) {

  // set the values in para to default parameters
  set_default_settings(para);

  if (argc < 2) {
    exit_with_help(argv[0], "libsvm", true);
  }
  
  for ( int i = 1; i < argc; ++i ) {
    if ( argv[i][0] != '-' ) {
      break;
    }
    if ( ++i >= argc ) {
      exit_with_help(argv[0], "libsvm", true);      
    }
    else if ( std::string ( argv[i-1] ) == "-epoch" ) {
      para->max_itrs = atoi ( argv[i] );
    }
    else if ( std::string ( argv[i-1] ) == "-data" ){
      data_file_name = std::string ( argv[i] );
    }
    else if ( std::string ( argv[i-1] ) == "-nthread" ) {
      para->total_num_threads = atoi ( argv[i] ) ;
    }
    else if ( std::string ( argv[i-1] ) == "-lambda" ) {
      lambda = atof ( argv[i] ) ;
    }
    else if ( std::string ( argv[i-1] ) == "-use_controller" ) {
      para->use_controller = atoi(argv[i]);
    }
    else {
      exit_with_help(argv[0], "libsvm", true);
    }
  }
  return;
}


void parse_input_argv_demo(Params* para,
                         int argc,
                         char *argv[]) {

  // set the values in para to default parameters
  set_default_settings(para);

  if (argc < 2) {
    exit_with_help(argv[0], "no_data", false);
  }
  
  for ( int i = 1; i < argc; ++i ) {
    if ( argv[i][0] != '-') {
      break;
    }
    if ( ++i >= argc ) {
      exit_with_help(argv[0], "no_data", false);
    }
    else if ( std::string ( argv[i-1] ) == "-epoch" ) {
      para->max_itrs = atoi ( argv[i] );
    }
    else if ( std::string ( argv[i-1] ) == "-nthread" ) {
      para->total_num_threads = atoi ( argv[i] ) ;
    }
    else if ( std::string ( argv[i-1] ) == "-problem_size" ) {
      para->problem_size = atoi ( argv[i] ) ;
    }
    else {
      exit_with_help(argv[0], "no_data", false);      
    }
  }
  return;
}



double log_loss_gradient_at_idx ( Matrix& A, Vector& b, Vector& Atx, int idx ) {
  double result = 0.;
  for ( unsigned i = 0; i < A.cols(); ++i )
    result -= A ( idx, i ) * b[i] / ( 1.+exp ( b[i] * Atx[i] ) );
  return result;
}


double log_loss_gradient_at_idx ( SpMat& A, Vector& b, Vector& Atx, int idx ) {

  double result = 0.;
  int i;
  for ( SpMat::InnerIterator it ( A, idx ); it; ++it ) {
    i = it.index();
    result -= it.value() * b[i] / ( 1.+exp ( b[i] * Atx[i] ) );
  }
  return result;
}


double l2_log_loss_objective (Vector& b, Vector& x, Vector& Atx, double lambda) {
  double tmp = 0.;
    for ( unsigned i = 0; i < b.size(); ++i ) {
      tmp += log ( 1. + exp ( -b[i] * Atx[i] ) );
    }
  double nrm = norm(x, 2);
  return 0.5 * lambda * nrm * nrm + tmp;
}

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

double l1_log_loss_objective (Vector& b, Vector& x, Vector& Atx, double lambda) {
  double tmp = 0.;
    for ( unsigned i = 0; i < b.size(); ++i ) {
      tmp += log ( 1. + exp ( -b[i] * Atx[i] ) );
    }
  double nrm = norm(x, 1);
  return lambda * nrm + tmp;
}


double log_loss(Vector& x, Vector& Atx, Vector& b) { 
  double tmp; 
  for (unsigned i = 0; i< b.size(); ++i) {
    tmp += log(1 + exp(-b[i] * Atx[i]));
  }
  return tmp;
}

double square_loss(Vector& x, Vector& Atx, Vector& b) {
  Vector v(b.size());
  for (unsigned i = 0; i< v.size(); ++i) {
    v[i] = Atx[i];
  }
  add(v, b, -1.);
  double tmp = norm(v), w = 1.;
  return tmp * tmp / 2.0;
}

double square_hinge_loss(Vector& x, Vector& Atx, Vector& b) {
  double tmp, M; 
  for (unsigned i = 0; i< b.size(); ++i) {
    M = abs(1 + b[i] * Atx[i]) + 1 - b[i] * Atx[i];
    tmp +=  M * M;
  }
  return tmp / 4.0;
}

double huber_loss(Vector& x, Vector& Atx, Vector& b, double delta) {
  double tmp, tmp2;
  for (unsigned i = 0; i< b.size(); ++i) {
    tmp2 = Atx[i] - b[i];
    tmp += (abs(tmp2) < delta)? (tmp2 * tmp2) : (delta * ( 2. *abs(tmp2) - delta));
  }
  return tmp / 2.;
}

double quad_func(Vector& x, Vector& Qx, Vector& c, double d) {
  return (dot(x, c) + dot(Qx, x) / 2. + d);
}

double huber_norm(Vector& x, double delta) {
  double res = 0.;
  for (int i = 0; i < x.size(); ++i) {
    res += (abs(x[i]) < delta) ? 0.5 * x[i] * x[i] : delta * (abs(x[i]) - delta / 2);
  }
  return res;
}

/******************************
 * Load data in libsvm format
 ******************************/
template<typename SparseMatrixType>
bool loadLibSVM(SparseMatrixType& mat, Vector& v, const std::string& filename) {
  typedef typename SparseMatrixType::Scalar Scalar;
  std::ifstream input(filename.c_str(),std::ios::in);
  
  int sym = 0;
  bool iscomplex = 0;
  bool isvector = 0;

  int M(0), N(0), NNZ(0);
  int row_idx, col_idx, loc;
  string line, str1, str2;
  double val;
  typedef Eigen::Triplet<Scalar,int> T;
  std::vector<T> elements;
  
  while(getline(input, line)) {
    istringstream iss(line);
    int a;
    string term;
    int cnt = 0;
    while (iss >> term) {
      if (cnt == 0) {
        v.push_back(stoi(term));
        cnt = 1;
      } else {
        ++NNZ;
        loc = term.find(':');
        str1 = term.substr(0, loc);
        str2 = term.substr(loc + 1);
        row_idx = atoi(str1.c_str()) - 1;
        N = max(row_idx + 1, N);
        col_idx = M;
        val = atof(str2.c_str());
        elements.push_back(T(row_idx, col_idx, val));
      }
    }
    ++M;  // number of rows in the file
  }
  // set up the matrix
  mat.resize(N, M);
  mat.reserve(NNZ);
  mat.setFromTriplets(elements.begin(), elements.end());
  input.close();
  return true;
}


template bool loadLibSVM <SpMat>(SpMat& mat, Vector& v, const std::string& filename);


/******************************
 * Load data in MATLAB format
 ******************************/

template<typename SparseMatrixType>
bool loadMatlabSparse(SparseMatrixType& mat, const std::string& filename) {
  typedef typename SparseMatrixType::Scalar Scalar;
  std::ifstream input(filename.c_str(),std::ios::in);
  if (!input) { // if the file is not open correctly
    cout << "error opening file \'" << filename << "\'\n";
    exit(1);
  }

  string line;
  istringstream iss;
  string term;
  // deal with first line:
  // expected: % ---------...--------
  while (getline(input, line)) { 
    iss.str(line);
    if (iss >> term) { // a nonempty line appears
      break;
    }
    else {
      if(iss.fail()) { // we have to clear the string stream once empty line is confronted
        iss.clear();
      }
    }
  }
  if (0 == term.compare("%")) {
    iss >> term;
    int i = 0;
    for(;i < term.length(); i++) { // expected to get  ---------...--------  
      if ('-' != term[i]) {
        cout << "error first line: not standard since \'%\' is followed by characters other than \'-\'\n";
        exit(1);            
      }
    }
  }
  else {
    cout << "error first line beginning not with %";
    exit(1);
  }
  getline(input, line); // get second line
  getline(input, line); // get third line

  // deal with fourth line
  // expected: % ---------...--------
  istringstream iss2;// we have to introduce a new object of istringstream -- just fix a bug and do not know why
  while (getline(input, line)) { 
    iss2.str(line);
    if (iss2 >> term) { // a nonempty line appears
      break;
    }
    else {
      if(iss2.fail()) { // we have to clear the string stream once empty line is confronted
        iss2.clear();
      }
    }
  }
  if (0 == term.compare("%")) {
    iss2 >> term;
    int i = 0;
    for(;i < term.length(); i++) { // expected to get  ---------...--------  
      if ('-' != term[i]) {
        cout << "error fourth line: not standard since \'%\' is followed by characters other than \'-\'\n";
        exit(1);            
      }
    }
  }
  else {
    cout << "error forth line beginning not with \'%\'";
    exit(1);
  }
  // deal with empty line
  while (getline(input, line)) { 
    iss.str(line);
    if (iss >> term) { // a nonempty line appears
      break;
    }
    else {
      if(iss.fail()) { // we have to clear the string stream once empty line is confronted
        iss.clear();
      }
    }
  }
  
  int M(0), N(0), NNZ(0);
  int row_idx, col_idx, location_break_string_term;
  
  double val;
  typedef Eigen::Triplet<Scalar,int> T;
  std::vector<T> elements;
  string str_row_idx, str_col_idx, str_rest_of_term;
  // expect to get LINE to be "A = zeros(5,6);"
  iss >> term;// expect to get "A"s
  iss >> term;// expect to get "="
  iss >> term;// expect to get "zero(5,6);"
  location_break_string_term = term.find('(');
  str_rest_of_term = term.substr(location_break_string_term + 1);
  location_break_string_term = str_rest_of_term.find(',');
  // read row index of dimension
  str_row_idx = str_rest_of_term.substr(0, location_break_string_term);
  M = stoi(str_row_idx);
  str_rest_of_term = str_rest_of_term.substr(location_break_string_term + 1);
  location_break_string_term = term.find(')');
  // read column index of dimension
  str_col_idx = str_rest_of_term.substr(0, location_break_string_term);
  N = stoi(str_col_idx);

  string str_val;
  // intializer to make sure the loop does not break due to previous errors
  while(getline(input, line)) {
    istringstream iss(line); // we have to introduce new object(s) of istringstream -- just fix a bug and do not know why
    iss >> term;// expect to get "A(1,1)"
    location_break_string_term = term.find('(');
    if (-1 == location_break_string_term) {
      break;
    }
    str_rest_of_term = term.substr(location_break_string_term + 1);
    location_break_string_term = str_rest_of_term.find(',');
    // read row index of the entry
    str_row_idx = str_rest_of_term.substr(0, location_break_string_term);
    row_idx = atoi(str_row_idx.c_str()) - 1;
    str_rest_of_term = str_rest_of_term.substr(location_break_string_term + 1);
    location_break_string_term = str_rest_of_term.find(')');
    // read column index of the entry
    str_col_idx = str_rest_of_term.substr(0, location_break_string_term);
    col_idx = atoi(str_col_idx.c_str()) - 1;
    str_rest_of_term = str_rest_of_term.substr(location_break_string_term + 1);
    iss >> term;// expect to get "="
    iss >> term;// expect to get "0.3403857266661332;"      
    location_break_string_term = term.find(';');
    // read value of the entry
    str_val = term.substr(0,location_break_string_term);
    val = atof(str_val.c_str());

    // insert the sparse entry
    elements.push_back(T(row_idx, col_idx, val));
    // adjust dimension of matrix if it is changed
    M = max(M, row_idx + 1);
    N = max(N, col_idx + 1);
  }
  // set up the matrix
  mat.resize(M, N);
  mat.setFromTriplets(elements.begin(), elements.end());

  input.close();
  return true;
}

template bool loadMatlabSparse <SpMat>(SpMat& mat, const std::string& filename);


//  Windows
#ifdef _WIN32
#include <Windows.h>
double get_wall_time(){
  LARGE_INTEGER time,freq;
  if (!QueryPerformanceFrequency(&freq)){
    //  Handle error
    return 0;
  }
  if (!QueryPerformanceCounter(&time)){
    //  Handle error
    return 0;
  }
  return (double)time.QuadPart / freq.QuadPart;
}
double get_cpu_time(){
  FILETIME a,b,c,d;
  if (GetProcessTimes(GetCurrentProcess(),&a,&b,&c,&d) != 0){
    //  Returns total user time.
    //  Can be tweaked to include kernel times as well.
    return
        (double)(d.dwLowDateTime |
                 ((unsigned long long)d.dwHighDateTime << 32)) * 0.0000001;
  }else{
    //  Handle error
    return 0;
  }
}

//  Posix/Linux
#else
#include <time.h>
#include <sys/time.h>
double get_wall_time(){
  struct timeval time;
  if (gettimeofday(&time,NULL)){
    //  Handle error
    return 0;
  }
  return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
double get_cpu_time(){
  return (double)clock() / CLOCKS_PER_SEC;
}
#endif
