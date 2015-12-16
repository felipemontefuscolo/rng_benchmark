#include <trng/config.hpp>
#include <cstdlib>
#include <iostream>
#include <omp.h>
#include <vector>
#include <iostream>
#include <cstdint>
#include <boost/preprocessor.hpp>
#include <random>
#include <boost/random.hpp>
#include <boost/math/distributions.hpp>
#include <pcg_random.hpp>

// Random123
#include <boost/random/counter_based_engine.hpp>
#include <boost/random/threefry.hpp>

// trng
#include <trng/yarn2.hpp>
#include <trng/uniform_dist.hpp>
#include <trng/normal_dist.hpp>
#include <trng/beta_dist.hpp>
#include <trng/gamma_dist.hpp>
#include <trng/exponential_dist.hpp>
#include <trng/poisson_dist.hpp>
#include <trng/bernoulli_dist.hpp>

#include "options.hpp"
#include "profiler.hpp"

#include <fstream>

using std::cout;
using std::endl;
using std::setprecision;

int set_w = 30;

Parameters pars;

// some sugar syntax
#define TEST        pars.test_reproducibility
#define n_samples   pars.num_samples
#define seed        pars.random_seed

int get_num_threads() {
  int n; 
#pragma omp parallel
#pragma omp single
  n = omp_get_num_threads();
  return n;
}


enum Dist_e {

  // continuous
  UNIFORM,
  NORMAL,
  BETA,
  GAMMA,
  EXPONENTIAL,
  
  // discrete
  POISSON,
  BERNOULLI

};


template<class Dist_t, Dist_e>
struct Params;

template<class Dist_t> struct Params<Dist_t, UNIFORM    > { static Dist_t createDist() { return Dist_t(get<0>(pars.uniform_params),get<1>(pars.uniform_params)); }  }; 
template<class Dist_t> struct Params<Dist_t, NORMAL     > { static Dist_t createDist() { return Dist_t(get<0>(pars.normal_params),get<1>(pars.normal_params)); }  }; 
template<class Dist_t> struct Params<Dist_t, BETA       > { static Dist_t createDist() { return Dist_t(get<0>(pars.beta_params),get<1>(pars.beta_params)); }  }; 
template<class Dist_t> struct Params<Dist_t, GAMMA      > { static Dist_t createDist() { return Dist_t(get<0>(pars.gamma_params),get<1>(pars.gamma_params)); }  }; 
template<class Dist_t> struct Params<Dist_t, EXPONENTIAL> { static Dist_t createDist() { return Dist_t(pars.exponential_param); }  }; 
template<class Dist_t> struct Params<Dist_t, POISSON    > { static Dist_t createDist() { return Dist_t(pars.poisson_param); }  }; 
template<class Dist_t> struct Params<Dist_t, BERNOULLI  > { static Dist_t createDist() { return Dist_t(pars.bernoulli_param); }  }; 

template<> struct Params<trng::bernoulli_dist<double>, BERNOULLI  >
                  { static trng::bernoulli_dist<double> createDist() { return trng::bernoulli_dist<double>(0.5, 0, 1); }  }; 

template<class E, class D>
void discard(E& e, D& d, size_t n) {
  for (size_t i =0; i<n; ++i)
    d(e);
}

// *****************************************************************************************
// -----------------------------------------------------------------------------------------
// *****************************************************************************************

template<class Dist_t, Dist_e dist_e >
double test_trng_yarn2() {

  vector<double> nums;
  bool const test = pars.test_reproducibility && get_num_threads()>1;

  if (test)
    nums.resize(n_samples);

  #pragma omp parallel
  {
    size_t size = omp_get_num_threads();     // get total number of processes
    size_t rank = omp_get_thread_num();      // get rank of current process
    Profiler::start(rank);
    trng::yarn2 r(seed);

    // split PRN sequences by leapfrog method
    r.split(size, rank);      // choose sub-stream no. rank out of size streams

    Dist_t u ( Params<Dist_t,dist_e>::createDist() );                 // random number distribution
    
    for (size_t i = rank; i < n_samples; i+=size) {
       auto x = u(r);
       if (test)
         nums[i] = x;
     }
     Profiler::stop(rank);
  }

  cout << std::fixed << std::showpoint;
  cout <<                    std::left << std::setw(set_w) <<  __func__;
  cout << setprecision(4) << std::left << std::setw(set_w) << Profiler::total_real_time();
  cout << setprecision(4) << std::left << std::setw(set_w) << Profiler::total_mflops();

  // compare to single-threaded
  if (test) {
    bool failed = false;   
 
    trng::yarn2 r(seed);
    Dist_t u ( Params<Dist_t,dist_e>::createDist() ); 
    for (size_t i = 0; i < n_samples; i++) {
       auto x = u(r);
       if (nums[i] != x)
         failed = true;
     }
    
     cout  << (failed ? "FAILED" : "PASSED") << endl;
  }
  else
    cout << endl;

  return Profiler::total_real_time();
}

// *****************************************************************************************
// -----------------------------------------------------------------------------------------
// *****************************************************************************************

template<class Dist_t, Dist_e dist_e >
double test_pcg32() {

  vector<double> nums;
  bool const test = pars.test_reproducibility && get_num_threads()>1;

  if (test)
    nums.resize(n_samples);

  #pragma omp parallel
  {
    size_t size = omp_get_num_threads();     // get total number of processes
    size_t rank = omp_get_thread_num();      // get rank of current process
    Profiler::start(rank);
    
    #pragma omp for 
    for (size_t i=0; i<n_samples; i++) {
       pcg32 r(seed, i);
       Dist_t u ( Params<Dist_t,dist_e>::createDist() );  // random number distribution
       auto x = u(r);
       if (test)
         nums[i] = x;
     }
    Profiler::stop(rank);
  }

  cout << std::fixed << std::showpoint;
  cout <<                    std::left << std::setw(set_w) <<  __func__;
  cout << setprecision(4) << std::left << std::setw(set_w) << Profiler::total_real_time();
  cout << setprecision(4) << std::left << std::setw(set_w) << Profiler::total_mflops();

   // compare to single-threaded
  if (test) {
    bool failed = false;   
     
    for (size_t i = 0; i < n_samples; i++) {
       pcg32 r(seed, i);
       Dist_t u ( Params<Dist_t,dist_e>::createDist() );
       auto x = u(r);
       if (nums[i] != x)
         failed = true;
     }
    
     cout  << (failed ? "FAILED" : "PASSED") << endl;
  }
  else
    cout << endl;

 

  return Profiler::total_real_time();
}

// *****************************************************************************************
// -----------------------------------------------------------------------------------------
// *****************************************************************************************

template<class Dist_t, Dist_e dist_e >
double test_std_mt19937_leap() {

  vector<double> nums;
  bool const test = pars.test_reproducibility && get_num_threads()>1;

  if (test)
    nums.resize(n_samples);

  #pragma omp parallel
  {
    size_t size = omp_get_num_threads();     // get total number of processes
    size_t rank = omp_get_thread_num();      // get rank of current process
    Profiler::start(rank);
    std::mt19937 r(seed);
   
    Dist_t u ( Params<Dist_t,dist_e>::createDist() );  // random number distribution

    discard(r, u, rank);

    for (size_t i = rank; i < n_samples; i+=size) {
       auto x = u(r);
       discard(r, u, size-1);
       if (test)
         nums[i] = x;
     }
     Profiler::stop(rank);
  }

  cout << std::fixed << std::showpoint;
  cout <<                    std::left << std::setw(set_w) <<  __func__;
  cout << setprecision(4) << std::left << std::setw(set_w) << Profiler::total_real_time();
  cout << setprecision(4) << std::left << std::setw(set_w) << Profiler::total_mflops();

   // compare to single-threaded
  if (test) {
    bool failed = false;   
    
    std::mt19937 r(seed);
    Dist_t u ( Params<Dist_t,dist_e>::createDist() );  // random number distribution
    for (size_t i = 0; i < n_samples; i++) {
       auto x = u(r);
       if (nums[i] != x)
         failed = true;
     }
    
     cout  << (failed ? "FAILED" : "PASSED") << endl;
  }
  else
    cout << endl;

  return Profiler::total_real_time();
}

// *****************************************************************************************
// -----------------------------------------------------------------------------------------
// *****************************************************************************************

template<class Dist_t, Dist_e dist_e >
double test_boost_mt19937_leap() {

  vector<double> nums;
  bool const test = pars.test_reproducibility && get_num_threads()>1;

  if (test)
    nums.resize(n_samples);

  #pragma omp parallel
  {
    size_t size = omp_get_num_threads();     // get total number of processes
    size_t rank = omp_get_thread_num();      // get rank of current process
    Profiler::start(rank);
    boost::mt19937 r(seed);
   
    Dist_t u ( Params<Dist_t,dist_e>::createDist() );  // random number distribution

    discard(r, u, rank);

    for (size_t i = rank; i < n_samples; i+=size) {
       auto x = u(r);
       discard(r, u, size-1);
       if (test)
         nums[i] = x;
     }
     Profiler::stop(rank);
  }

  cout << std::fixed << std::showpoint;
  cout <<                    std::left << std::setw(set_w) <<  __func__;
  cout << setprecision(4) << std::left << std::setw(set_w) << Profiler::total_real_time();
  cout << setprecision(4) << std::left << std::setw(set_w) << Profiler::total_mflops();

   // compare to single-threaded
  if (test) {
    bool failed = false;   
    
    std::mt19937 r(seed);
    Dist_t u ( Params<Dist_t,dist_e>::createDist() );  // random number distribution
    for (size_t i = 0; i < n_samples; i++) {
       auto x = u(r);
       if (nums[i] != x)
         failed = true;
     }
    
     cout  << (failed ? "FAILED" : "PASSED") << endl;
  }
  else
    cout << endl;

  return Profiler::total_real_time();
}




// *****************************************************************************************
// -----------------------------------------------------------------------------------------
// *****************************************************************************************


template<class Dist_t, Dist_e dist_e >
double test_random123() {

  vector<double> nums;
  bool const test = pars.test_reproducibility && get_num_threads()>1;

  if (test)
    nums.resize(n_samples);

  #pragma omp parallel
  {
    size_t size = omp_get_num_threads();     // get total number of processes
    size_t rank = omp_get_thread_num();      // get rank of current process
    Profiler::start(rank);
    typedef boost::random::threefry<4, uint32_t> Prf;
    boost::random::counter_based_engine<uint32_t, Prf, 32> r(seed);

    Dist_t u ( Params<Dist_t,dist_e>::createDist() ); 

    #pragma omp for
    for (size_t i = 0; i < n_samples; ++i) {
       r.restart({i, size_t(0), size_t(0)});
       auto x = u(r);
       if (test)
         nums[i] = x;
     }
    Profiler::stop(rank);
  }

  cout << std::fixed << std::showpoint;
  cout <<                    std::left << std::setw(set_w) <<  __func__;
  cout << setprecision(4) << std::left << std::setw(set_w) << Profiler::total_real_time();
  cout << setprecision(4) << std::left << std::setw(set_w) << Profiler::total_mflops();

   // compare to single-threaded
  if (test) {
    bool failed = false;   
    typedef boost::random::threefry<4, uint32_t> Prf;
    boost::random::counter_based_engine<uint32_t, Prf, 32> r(seed);

    Dist_t u ( Params<Dist_t,dist_e>::createDist() ); 
    for (size_t i = 0; i < n_samples; i++) {
       r.restart({i, size_t(0), size_t(0)});
       auto x = u(r);
       if (nums[i] != x)
         failed = true;
     }
    
     cout  << (failed ? "FAILED" : "PASSED") << endl;
  }
  else
    cout << endl;

  return Profiler::total_real_time();
}


template<class T>
void print_to_file(vector<T>& v) {
  
}

int main(int argc, char* argv[]) {

  readParameters(pars, argc, argv);
  pars.print_values();

  Profiler::initialize(1); 

  for (auto n_threads : pars.num_threads) {

    omp_set_num_threads(n_threads);

    Profiler::change_num_threads(n_threads);

    bool const test = pars.test_reproducibility && n_threads>1;

    #pragma omp parallel
    #pragma omp single
    {
      cout << "Running with  " << omp_get_num_threads() << " threads" << endl ;
      cout << "-------------------------------------------------------------------------------------------------------------- " << endl ;
    }
    cout <<  std::left << std::setw(set_w) << " ";
    cout <<  std::left << std::setw(set_w) << "Real time (s)";
    cout <<  std::left << std::setw(set_w) << "MFLOP/s ";
    if (test)
      cout <<  std::left << std::setw(set_w) <<  "reproducibility test" << endl;
    else
      cout <<  endl;
    
 
 
    cout << "Uniform : " << "\n";
    test_trng_yarn2             <trng::uniform_dist<double>, UNIFORM >(); 
    test_pcg32                  <std::uniform_real_distribution<double>, UNIFORM>();
    test_std_mt19937_leap       <std::uniform_real_distribution<double>, UNIFORM>();
    test_boost_mt19937_leap     <boost::random::uniform_real_distribution<double>, UNIFORM>();
    test_random123              <boost::random::uniform_real_distribution<double>, UNIFORM>();
    cout << "\n";
  
    cout << "Normal : " << "\n";
    test_trng_yarn2             <trng::normal_dist<double>, NORMAL >(); 
    test_pcg32                  <std::normal_distribution<double>, NORMAL>();
    test_std_mt19937_leap       <std::normal_distribution<double>, NORMAL>();
    test_boost_mt19937_leap     <boost::random::normal_distribution<double>, NORMAL>();
    test_random123              <boost::random::normal_distribution<double>, NORMAL>();
    cout << "\n";
  
  //  cout << "Beta : " << "\n";
  //  test_trng_yarn2             <trng::uniform_dist<double>, BETA >(); 
  //  test_pcg32                  <std::uniform_real_distribution<double>, BETA>();
  //  test_std_mt19937_leap       <std::uniform_real_distribution<double>, BETA>();
  //  test_boost_mt19937_leap     <boost::random::uniform_real_distribution<double>, BETA>();
  //  test_random123              <boost::random::uniform_real_distribution<double>, BETA>();
  //  cout << "\n";
  
    cout << "Gamma : " << "\n";
    test_trng_yarn2             <trng::gamma_dist<double>, GAMMA >(); 
    test_pcg32                  <std::gamma_distribution<double>, GAMMA>();
    test_std_mt19937_leap       <std::gamma_distribution<double>, GAMMA>();
    test_boost_mt19937_leap     <boost::random::gamma_distribution<double>, GAMMA>();
    test_random123              <boost::random::gamma_distribution<double>, GAMMA>();
    cout << "\n";
  
    cout << "Exponential : " << "\n";
    test_trng_yarn2             <trng::exponential_dist<double>, EXPONENTIAL >(); 
    test_pcg32                  <std::exponential_distribution<double>, EXPONENTIAL>();
    test_std_mt19937_leap       <std::exponential_distribution<double>, EXPONENTIAL>();
    test_boost_mt19937_leap     <boost::random::exponential_distribution<double>, EXPONENTIAL>();
    test_random123              <boost::random::exponential_distribution<double>, EXPONENTIAL>();
    cout << "\n";
  
    cout << "Poisson : " << "\n";
    test_trng_yarn2             <trng::poisson_dist, POISSON >(); 
    test_pcg32                  <std::poisson_distribution<int>, POISSON>();
    test_std_mt19937_leap       <std::poisson_distribution<int>, POISSON>();
    test_boost_mt19937_leap     <boost::random::poisson_distribution<int>, POISSON>();
    test_random123              <boost::random::poisson_distribution<int>, POISSON>();
    cout << "\n";
  
    cout << "Bernoulli : " << "\n";
    test_trng_yarn2             <trng::bernoulli_dist<double>, BERNOULLI >(); 
    test_pcg32                  <std::bernoulli_distribution, BERNOULLI>();
    test_std_mt19937_leap       <std::bernoulli_distribution, BERNOULLI>();
    test_boost_mt19937_leap     <boost::random::bernoulli_distribution<double>, BERNOULLI>();
    test_random123              <boost::random::bernoulli_distribution<double>, BERNOULLI>();
    cout << "\n";

  }

}




