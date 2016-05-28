#ifndef AROCK_INCLUDE_CONTROLLER_H
#define AROCK_INCLUDE_CONTROLLER_H

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

std::mutex g_set_schemes_mutex; //mutex for thread safe interaction with set
std::mutex g_calculate_residual_mutex;
std::condition_variable g_condition_calculate_fpr;
bool g_calculate_average_fpr = false;


/*
 * Controller's knowledge of a worker thread
 */
template<typename SplittingScheme>
struct worker_info{
  Params *p;
  SplittingScheme* scheme;
  std::atomic<bool> active;
  size_t last_updated; //last time the worker processed an update
  size_t  max_delay;
  std::vector<int> delays;
  worker_info(Params *_p, SplittingScheme* _scheme): p(_p), scheme(_scheme), active(true), last_updated(0), max_delay(0) {
    delays.reserve(_p->max_itrs);
  };
  worker_info(const worker_info& w): p(w.p), scheme(w.scheme), active(w.active.load() ), last_updated(w.last_updated), max_delay(w.max_delay), delays(w.delays) {}; 
};



template<typename T>
struct Controller {
  Vector average_fpr; //stores an approximate fixed point residual: currently a moving window average
  double window_size; //constant for sliding window (1-window_size)*FPR[i]+window_size*new_residual
  double old_fpr_norm; //stores old fpr norm for purposes of comparison
  
  std::atomic<size_t> epoch_iter; //number of iterations bewteen controller calculations
  std::atomic<size_t> total_iter; //total number of iterations
  std::atomic<int> max_delay; //maximum delay experienced
  size_t update_epoch; //how many updates to a coordinate before we determine a new stepsize
  size_t num_worker; //num of worker threads
  double arock_max_stepsize; 
  double arock_min_stepsize;
  std::vector<size_t> coord_num_updates; //stores the last time a coordinate was updated
  unordered_map< std::thread::id ,worker_info<T> > worker_state; //map from worker threads to their state

  //this constructor is here for legacy reasons, deprecate soon
  Controller () { }
  Controller (const Controller& c) : epoch_iter(c.epoch_iter.load()), total_iter(c.total_iter.load()), average_fpr(c.average_fpr), window_size(c.window_size), old_fpr_norm(c.old_fpr_norm), update_epoch(c.update_epoch), coord_num_updates(c.coord_num_updates), worker_state(c.worker_state), num_worker(c.num_worker), arock_max_stepsize(c.arock_max_stepsize), arock_min_stepsize(c.arock_min_stepsize), max_delay(c.max_delay.load()) { }
  
  //note arock_max_stepsize/min_stepsize should become part of the struct PARAM in future versions along with window_size
  Controller (Params p) : epoch_iter(0), total_iter(0), window_size(.8), old_fpr_norm(std::numeric_limits<double>::max()), update_epoch(p.problem_size / 4), num_worker(0), arock_max_stepsize(100), arock_min_stepsize(.000001), average_fpr(p.problem_size,0), coord_num_updates(p.problem_size,0), max_delay(0) {
    worker_state.reserve(p.total_num_threads+10);
  }
  
  void process_update (std::thread::id id, size_t coord_idx, double residual) {
    
    //update worker statistics;
    auto& w = worker_state.at(id); 
    int delay = total_iter.load() - w.last_updated;
    w.last_updated =  total_iter.load();
    if ( delay > max_delay.load() ) {
      max_delay.store(delay);
    }
    if ( delay > w.max_delay ) {
      w.max_delay = delay;
    }
    w.delays.push_back(delay);
    //average in the new residual data
    average_fpr[coord_idx] = ( 1 - window_size ) * average_fpr[coord_idx] + (window_size) * residual;
    
    //update the counters
    ++coord_num_updates[coord_idx];
    total_iter.fetch_add(1); //ATOMIC UPDATE
    epoch_iter.fetch_add(1); //ATOMIC UPDATE

    //if epoch_iter goes above the update_epoch, alert controller to do calculations
    if ( epoch_iter.load() > update_epoch ) { 
      epoch_iter.store(0);
      g_calculate_average_fpr = true; //signal the fpr norm should be calculated
      g_condition_calculate_fpr.notify_one();
    }

  }
  
  void update_average_fpr (std::thread::id id, size_t coord_idx, double residual) {
    
    //average in the new residual data
    average_fpr[coord_idx] = ( 1 - window_size ) * average_fpr[coord_idx] + (window_size) * residual;
    
    //update the counters
    ++coord_num_updates[coord_idx];
    total_iter.fetch_add(1); //ATOMIC UPDATE
    epoch_iter.fetch_add(1); //ATOMIC UPDATE
    //if epoch_iter goes above the update_epoch, alert controller to do calculations
    
    if ( epoch_iter.load() > update_epoch ) { 
      epoch_iter.store(0);
      g_calculate_average_fpr = true; //signal the fpr norm should be calculated
      g_condition_calculate_fpr.notify_one();
    }
    
  }

  // remove worker
  void remove_worker( std::thread::id  id) {
    //remove members of the set schemes in a thread safe manner
    std::unique_lock<std::mutex> lock(g_set_schemes_mutex);
    
    auto loc = worker_state.find(id);
    if( loc != worker_state.end() ){
      loc->second.active.store(false);
      --num_worker;
    }
    lock.unlock();
    
    //notify controller thread that it might have nothing to control
    g_condition_calculate_fpr.notify_one();
  }

  // add worker 
  void add_worker(std::thread::id id, T& scheme, Params* p) {
    //add members of the set schemes in a thread safe manner
    std::unique_lock<std::mutex> lock(g_set_schemes_mutex);
    if( worker_state.find(id) == worker_state.end() ) {
      worker_state.insert( std::make_pair< std::thread::id , worker_info<T> >(std::move(id) , worker_info<T>(p, &scheme)));
      ++num_worker;
    }
    lock.unlock();
  } 
};

//runs while there are still active workers
template<typename T> void Controller_loop(Controller<T>& controller) {

  size_t* p_num_worker=&controller.num_worker;
  
  //this is a lambda function that is used to wake up the controller thread. Thread activates when wait_func returns true
  auto wait_func = [=]{return g_calculate_average_fpr || *p_num_worker == 0 ;};
  
  //controller idles until a set startup time passes, or a worker contacts controller
  double start_time = get_wall_time();
  while ( controller.num_worker == 0 && get_wall_time() - start_time < 30) {
  }
 
  while ( controller.num_worker > 0 ) { 
   
    //instantiate lock for condition variable
    std::unique_lock<std::mutex> condition_lock(g_calculate_residual_mutex);
   
    //idle until wait_func returns true, then acquire lock on condition_lock
    g_condition_calculate_fpr.wait(condition_lock, wait_func);
    g_calculate_average_fpr = false;
    
    auto newnorm = norm(controller.average_fpr, 2);
    //if newnorm is less than old norm, increase stepsize, else decrease stepsize
    //note we back off harder than we increase, to prevent endless oscillations
    std::cout << newnorm << std::endl;
    if ( newnorm < 0.9 * controller.old_fpr_norm ) {
      
      //when reading from set, no one insert/delete
      std::unique_lock<std::mutex> lock(g_set_schemes_mutex);
      //update each schemes stepsize
      int asdf = 0;
      for ( auto& it : controller.worker_state ) {
	if( it.second.active.load() == true ) {
          
	  auto& step = (it.second).scheme -> relaxation_step_size;
          if(!asdf) { std::cout << "step size " << step << std::endl; ++asdf;}
          step = (1.5*step < controller.arock_max_stepsize) ? (1.5)*step : controller.arock_max_stepsize;
	}
      }
      lock.unlock();
    }
    else {
      //when reading from set, no one insert/delete
      std::unique_lock<std::mutex> lock(g_set_schemes_mutex);
      int asdf = 0;
      for ( auto& it : controller.worker_state ) {
	if( it.second.active.load() == true ) {
	  auto& step = (it.second).scheme -> relaxation_step_size;
          if(!asdf) { std::cout << "step size " << step << std::endl; ++asdf;}          
	  step = (step/2 > controller.arock_min_stepsize) ? step/2 : controller.arock_min_stepsize;
	}
      } 
      lock.unlock();
      
    }
    controller.old_fpr_norm = newnorm; //new FPR to considerx
    //unlock and wait again
    condition_lock.unlock();    
  }
  std::cout << "ending approx residual is " << controller.old_fpr_norm << std::endl;
  std::cout << "Max delay per thread: \n";
  for(auto& it : controller.worker_state) {
    std::cout << "Thread id: " << it.first << " max_delay: " << it.second.max_delay <<std::endl;
  }
}


#endif
