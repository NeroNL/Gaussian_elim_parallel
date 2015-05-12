#ifndef _BARRIER_H
#define _BARRIER_H

#include <mutex>
#include <cassert>
#include <condition_variable>
#include <vector>
#include <atomic>

using namespace std;

#ifdef LOG_BARRIER
class barrier
{
public:
    // Fill in your definition of the barrier constructor
    barrier(int NT = 2) { }
    // Fill in your definition of the barrier synchronization function
    void bsync(int TID) {
    }
};
#else
class barrier
{
public:
  explicit barrier(int a_numThreads):m_nt(a_numThreads), m_ndone(0){ };
  inline void bsync(int dummy);
  
private:
  //barrier(const barrier&);  // disable copying barriers.
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

