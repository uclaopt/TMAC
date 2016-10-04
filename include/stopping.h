#ifndef TMAC_INCLUDE_STOPPING_H
#define TMAC_INCLUDE_STOPPING_H

#include "parameters.h"
#include "splitting_schemes.h"
#include "util.h"
#include <atomic> //protect some variables from race conditions
#include <unordered_set> //might need this in the future
#include <unordered_map>
#include <cmath>
#include <limits>
#include <thread>
#include <vector>
#include<condition_variable>
#include<mutex>
#include<thread>

/*
 * The following is experimental functionality
 * it is kept separate from the current controller until
 * a decisions is made in regard to merging the functionality
 */


template<typename SplittingScheme> struct Stopper{
  SplittingScheme* scheme;
  Stopper(SplittingScheme* scheme);
  void loop();
  std::thread spawn_loop();
  //no op functions to interface with code at the moment
  void process_update (std::thread::id id, size_t coord_idx, double residual) {}
  void add_worker(std::thread::id id, SplittingScheme& scheme, Params* p) {}
  void remove_worker( std::thread::id  id) {}
};

template<typename SplittingScheme>
Stopper<SplittingScheme>::Stopper(SplittingScheme* scheme):scheme(scheme){};

template<typename SplittingScheme>
void Stopper<SplittingScheme>::loop(){
  while( !scheme->eval_stop_crit() ) {
  }
}
template<typename SplittingScheme>
std::thread Stopper<SplittingScheme>::spawn_loop() {
  return std::thread(&Stopper<SplittingScheme>::loop,this);
}


#endif
