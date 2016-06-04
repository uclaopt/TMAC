#ifndef BARRIER_H
#define BARRIER_H

#include<atomic>
#include<mutex>

class Barrier{
 public:
  void wait();
  Barrier(const Barrier&)= delete; //dangerous, so disallow
  Barrier(int thread_count=0);
  void init(int nthreads){thread_count = nthreads;}
 private:
  std::atomic<int> counter[2];
  int lock[2];
  int cur_idx;
  int thread_count;

};

Barrier::Barrier(int thread_count):thread_count(thread_count){
    counter[0].store(0);counter[1].store(0);
    lock[0]=lock[1]=cur_idx=0;
  }

void Barrier::wait(){
  //int idx= cur_idx.load();
  int idx= cur_idx;
  //if the current lock is not locked, lock it
  //if(lock[idx].load() == 0){
  if(lock[idx] == 0){
    //lock[idx].store(1);
    lock[idx] = 1;
  }
  //val is equal to the value of the counter before increment
  int val= counter[idx].fetch_add(1);
  
  //if wait has been called at least thread_count many times
  if( val >= thread_count-1){
    counter[idx].store(0); //counter for current idx has reached its max
    lock[idx]=0;
    cur_idx= cur_idx ^ 1;
    //cur_idx.fetch_xor(1); //replace 0 with 1 and 1 with 0
    //lock[idx].store(0); //deactive lock associated with the old active idx
  }
while(lock[idx] == 1);
//while(lock[idx].load() == 1 ); //spin while waiting for other threads to hit barrier
}

#endif
