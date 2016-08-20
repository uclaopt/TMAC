#ifndef TMAC_INCLUDE_TMAC_H
#define TMAC_INCLUDE_TMAC_H
#include "matrices.h"
#include "algebra.h"
#include "worker.h"
#include "controller.h"
#include "parameters.h"
#include "result.h"
#include "driver.h"
#include <thread>

template<typename Splitting>
void TMAC(Splitting op, Params parameters, Result res, Controller<Splitting> controller = Controller<Splitting>()) {
  if (parameters.async == true) {
    //construct a controller based on params if use failed to specify one
    if(parameters.use_controller == true && controller.def == true){
      Controller<Splitting> param_cont(parameters);
	  cout << parameters.total_num_threads << endl;
      AROCK (op, parameters, res, param_cont);
    }
    else{
    AROCK (op, parameters, res, controller);
    }
  } else {
    SYNC(op, parameters);
  }
}

#endif
