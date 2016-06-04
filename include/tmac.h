#ifndef TMAC_INCLUDE_TMAC_H
#define TMAC_INCLUDE_TMAC_H
#include "matrices.h"
#include "algebra.h"
#include "worker.h"
#include "controller.h"
#include "parameters.h"
#include "driver.h"
#include <thread>

template<typename Splitting>
void TMAC (Splitting op, Params parameters, Controller<Splitting> controller = Controller<Splitting>()) {
  if (parameters.async == true) {
    AROCK (op, parameters,controller);
  } else {
    SYNC(op, parameters);
  }
}

#endif
