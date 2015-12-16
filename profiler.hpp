#ifndef I4I_PROFILER_HPP
#define I4I_PROFILER_HPP

#include "papi.h"
#include <iostream>
#include <vector>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#include <cstdlib>
#include <cstdio>
#else
#error "This compiler does not understand OPENMP"
#endif

using std::printf;
using std::cerr;
using std::cout;
using std::endl;

class Profiler {

  struct Thread {
    float real_time;
    float proc_time;
    float mflops;
    long long flpops;
    bool is_running;
  };

  static std::vector<Thread> data;

public:
 
  static float real_time(int tid)  {return data.at(tid).real_time;       }
  static float proc_time(int tid)  {return data.at(tid).proc_time;       }
  static float mflops(int tid)     {return data.at(tid).mflops;          }
  static float gflops(int tid)     {return data.at(tid).mflops * 0.001;  }
  static long long flpops(int tid) {return data.at(tid).flpops;          }

  //Profiler() {}

  // must call this before everything.
  static void initialize(int num_threads) {

    data.resize(num_threads, Thread{0.,0.,0.,0u,false});

    int retval = PAPI_is_initialized();
    //if (retval == PAPI_NOT_INITED) {
		  retval = PAPI_library_init( PAPI_VER_CURRENT );
	    if ( retval != PAPI_VER_CURRENT )
        throw;
    //}
    retval =	PAPI_thread_init( ( unsigned	long ( * )( void ) ) ( omp_get_thread_num ) );
   	if ( retval != PAPI_OK ) {
       throw;
   	}
  }

  // must call this everytime you change the number of threads
  static void change_num_threads(int num_threads)
  {
    data.resize(num_threads, Thread{0.,0.,0.,0u,false}); 
  }

  // start must be inside the parallel region
  static void start(int thread_id) {

    float real_time, proc_time,mflops;
    long long flpops;
    int retval; 

    if (data.at(thread_id).is_running) {
      cerr << "ERROR: Profiler::start: profiling already started\n";
      throw;
    }      
    data.at(thread_id).is_running = true;

    if((retval=PAPI_flops(&real_time,&proc_time,&flpops,&mflops)) < PAPI_OK)
    {
      printf("Could not initialise PAPI_flops \n");
      printf("Your platform may not support floating point operation event.\n");
      printf("retval: %d\n", retval);
      exit(1);
    }

  }


  // stop also must be inside the parallel region
  static void stop(int thread_id) {
    if (!data.at(thread_id).is_running ) {
      cerr << "ERROR: Profiler::stop: profiling not started \n";
      throw;
    }      
 
    float& real_time = data.at(thread_id).real_time;
    float& proc_time = data.at(thread_id).proc_time;
    float& mflops = data.at(thread_id).mflops;
    long long& flpops = data.at(thread_id).flpops;
    int retval; 

    if((retval=PAPI_flops( &real_time, &proc_time, &flpops, &mflops))<PAPI_OK)
    {
      printf("retval: %d\n", retval);
      exit(1);
    }

    data.at(thread_id).is_running = false;

    PAPI_unregister_thread();
    long long values [] = {0};
    PAPI_stop_counters( values, 0);  
  }

  // outside parallel
  static float total_mflops()
  {
    float a = 0, b = 0;
    for (auto d : data) {
      a += d.flpops;
      b = std::max(d.proc_time, b);
    }
    return a / b * 0.000001;
  }
 
  // outside parallel
  static float total_real_time()
  {
    float t = 0;
    for (auto d : data)
      t = std::max(d.real_time, t);
    return t;
  }

  // must be called outside the parallel region
  static void check() {
    for (unsigned i = 0; i<data.size(); ++i) {
      if (data.at(i).is_running) {
        printf("Thread race condition at %d. Terminating\n", (int)i);
        throw;
      }
    }
  }


}; // end class


std::vector<Profiler::Thread> Profiler::data;




///
/// EXAMPLE
///



// the function f() does some time-consuming work
//void f()
//{
//    volatile double d = 0;
//    for(int n=0; n<10000; ++n)
//       for(int m=0; m<10000; ++m)
//           d += d*n*m;
//}
//
//void Thread( unsigned long n )
//{
//  int id = omp_get_thread_num();
//	printf( "Thread %d started\n", id );
//
//  Profiler::start( id );
//  for(int j =0 ; j<n; j++)
//    f();
//  Profiler::stop( id );
//
//	printf( "Thread %d finished: ", id );
//  printf( "time = %f, mflops = %f, flops = %lld\n", Profiler::real_time(id), Profiler::mflops(id), Profiler::flpops(id));
//}
//
//int main( int argc, char **argv )
//{
//
//  int n_threads;
//
//#pragma omp parallel
//#pragma omp single
//  n_threads = omp_get_num_threads();
// 
//  Profiler::initialize(n_threads);
//
//  cout << "running for " << n_threads << " threads" << endl;
//#pragma omp parallel
//	{
//		Thread(  ( 0*omp_get_thread_num(  ) + 1 ) );
//	}
//
//  cout << "total MFLOPS = " << Profiler::total_mflops() << endl;
//  cout << "total real time = " << Profiler::total_real_time() << endl;
//
//  Profiler::check();
//
//  cout << "running for " << 2 << " threads" << endl;
//  omp_set_num_threads(2);
//  Profiler::change_num_threads(2);
//#pragma omp parallel
//	{
//		Thread(  ( omp_get_thread_num(  ) + 1 ) );
//	}
//  
//  cout << "total MFLOPS = " << Profiler::total_mflops() << endl;
//  cout << "total real time = " << Profiler::total_real_time() << endl;
//
//	exit( 0 );
//}


#endif

