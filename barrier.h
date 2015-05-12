#ifndef _BARRIER_H
#define _BARRIER_H

#include <mutex>
#include <cassert>
#include <condition_variable>
#include <vector>
#include <atomic>

using namespace std;

#ifdef LOG_BARRIER

struct lock_control{
  atomic_flag parent_lock = ATOMIC_FLAG_INIT;
  atomic_flag left_lock = ATOMIC_FLAG_INIT;
  atomic_flag right_lock = ATOMIC_FLAG_INIT;
  bool left = false;
  bool right = false;
  bool parent = false;
};

class barrier
{
  vector<int> v_a (cb.NT, 0);
  vector<lock_control> v_lc (cb.NT);

public:
    // Fill in your definition of the barrier constructor
    barrier(int NT = 2) {
      for(unsigned int i = 0; i < cb.NT; i++){
        (v_lc.at(i)).parent_lock.set_and_test(memory_order_acquire);
        (v_lc.at(i)).left_lock.set_and_test(memory_order_acquire);
        (v_lc.at(i)).right_lock.set_and_test(memory_order_acquire);
      }
    }
    // Fill in your definition of the barrier synchronization function
    void bsync(int TID) {
      v_a.at(TID) = TID;
      int parent = (TID-1)/2;
      int left_child = TID*2+1;
      int right_child = TID*2+2;
      bool pass = false;
      bool left_c = true;
      bool right_c = true;

      if(parent < 0){
        parent = 0;
      }

      if(!pass){
        if(left_child > cb.NT-1 || (v_lc.at(TID)).left == false){
          (v_lc.at(TID)).left_lock.clear(memory_order_release);
          (v_lc.at(TID)).left = true;
        }
        else
          while((v_lc.at(TID)).left != left_c);

        if(right_child > cb.NT-1 || (v_lc.at(TID)).right == false){
          (v_lc.at(TID)).right_lock.clear(memory_order_release);
          (v_lc.at(TID)).right = true;
        }
        else
          while((v_lc.at(TID)).right != right_c);

        if((v_lc.at(TID)).left && (v_lc.at(TID)).right){
          (v_lc.at(TID)).parent_lock.clear(memory_order_release);
          (v_lc.at(TID)).parent = true;
          pass = true;
        }
      }
      else{
        if((TID-1)%2 == 0)
          (v_lc.at(parent)).right = true;
        else if((TID-1)%2 == 1)
          (v_lc.at(parent)).left = true;
      }

      /*bool left_c = (v_lc.at(left_child)).parent;
      bool right_c = (v_lc.at(left_child)).parent;


      if(left_c && right_c){
        v_lc.at(TID)).parent_lock.clear(memory_order_release);
         (v_lc.at(TID)).parent = true;
        return;
      }
      else
        (v_lc.at(TID)).parent_lock.set_and_test(memory_order_acquire);

      if(!left_c)
        (v_lc.at(TID)).left_lock.set_and_test(memory_order_acquire);
      else
        (v_lc.at(TID)).left_lock.clear(memory_order_release);

      if(!right_c)
         (v_lc.at(TID)).right_lock.set_and_test(memory_order_acquire);
      else
        (v_lc.at(TID)).right_lock.clear(memory_order_release);

      */

    }

};
#else
class barrier
{
public:
  explicit barrier(int a_numThreads):m_nt(a_numThreads), m_ndone(0){ };
  inline void bsync(int dummy);
  
private:
  barrier(const barrier&);  // disable copying barriers.
  barrier(); //no null constructor.
  barrier& operator=(const barrier&); // disable assignment
  int m_nt;
  int m_ndone;
  std::mutex m_mutex;
  std::condition_variable m_cv;
};

inline void barrier::bsync(int dummy)
{
  std::unique_lock<std::mutex> lk(m_mutex,std::defer_lock);
  lk.lock();
  if(++m_ndone < m_nt)
     m_cv.wait(lk);
  else {
      m_ndone = 0;
      m_cv.notify_all();
      return;
    }
}


#endif
#endif

