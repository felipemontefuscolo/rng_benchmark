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

using std::cout;
using std::endl;
using std::setprecision;

unsigned long seed = 5489u;
size_t n_samples = 1000;
int set_w = 30;

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

template<class Dist_t> struct Params<Dist_t, UNIFORM    > { static Dist_t createDist() { return Dist_t(0.,1.); }  }; 
template<class Dist_t> struct Params<Dist_t, NORMAL     > { static Dist_t createDist() { return Dist_t(0.,1.); }  }; 
template<class Dist_t> struct Params<Dist_t, BETA       > { static Dist_t createDist() { return Dist_t(1.,2.); }  }; 
template<class Dist_t> struct Params<Dist_t, GAMMA      > { static Dist_t createDist() { return Dist_t(1.,2.); }  }; 
template<class Dist_t> struct Params<Dist_t, EXPONENTIAL> { static Dist_t createDist() { return Dist_t(1.  ); }  }; 
template<class Dist_t> struct Params<Dist_t, POISSON    > { static Dist_t createDist() { return Dist_t(4  ); }  }; 
template<class Dist_t> struct Params<Dist_t, BERNOULLI  > { static Dist_t createDist() { return Dist_t(0.5); }  }; 

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

  double start = omp_get_wtime();

  #pragma omp parallel
  {
    size_t size = omp_get_num_threads();     // get total number of processes
    size_t rank = omp_get_thread_num();      // get rank of current process
    trng::yarn2 r(seed);

    // split PRN sequences by leapfrog method
    r.split(size, rank);      // choose sub-stream no. rank out of size streams

    Dist_t u ( Params<Dist_t,dist_e>::createDist() );                 // random number distribution
    
    for (size_t i = rank; i < n_samples; i+=size) {
       auto x = u(r);
     }
  }

  double end = omp_get_wtime();
  cout << std::fixed << std::showpoint;
  cout << setprecision(4) << std::left << std::setw(set_w) <<  __func__ <<  (end-start) << "\n";
  return end-start;
}

// *****************************************************************************************
// -----------------------------------------------------------------------------------------
// *****************************************************************************************

template<class Dist_t, Dist_e dist_e >
double test_pcg32() {

  double start = omp_get_wtime();

  #pragma omp parallel
  {
    size_t size = omp_get_num_threads();     // get total number of processes
    size_t rank = omp_get_thread_num();      // get rank of current process
    pcg32 r(seed, rank);

    Dist_t u ( Params<Dist_t,dist_e>::createDist() );  // random number distribution
    
    #pragma omp for 
    for (size_t i=0; i<n_samples; i++) {
       auto x = u(r);
     }
  }

  double end = omp_get_wtime();
  cout << std::fixed << std::showpoint;
  cout << std::left << std::setw(set_w) <<  __func__ <<  (end-start) << "\n";
  return end-start;
}

// *****************************************************************************************
// -----------------------------------------------------------------------------------------
// *****************************************************************************************

template<class Dist_t, Dist_e dist_e >
double test_std_mt19937_leap() {

  double start = omp_get_wtime();

  #pragma omp parallel
  {
    size_t size = omp_get_num_threads();     // get total number of processes
    size_t rank = omp_get_thread_num();      // get rank of current process
    std::mt19937 r(seed);
   
    Dist_t u ( Params<Dist_t,dist_e>::createDist() );  // random number distribution

    discard(r, u, rank);

    for (size_t i = rank; i < n_samples; i+=size) {
       auto x = u(r);
       discard(r, u, size-1);
     }
  }

  double end = omp_get_wtime();
  cout << std::fixed << std::showpoint;
  cout << setprecision(4) << std::left << std::setw(set_w) <<  __func__ <<  (end-start) << "\n";
  return end-start;
}

// *****************************************************************************************
// -----------------------------------------------------------------------------------------
// *****************************************************************************************

template<class Dist_t, Dist_e dist_e >
double test_boost_mt19937_leap() {

  double start = omp_get_wtime();

  #pragma omp parallel
  {
    size_t size = omp_get_num_threads();     // get total number of processes
    size_t rank = omp_get_thread_num();      // get rank of current process
    boost::mt19937 r(seed);
   
    Dist_t u ( Params<Dist_t,dist_e>::createDist() );  // random number distribution

    discard(r, u, rank);

    for (size_t i = rank; i < n_samples; i+=size) {
       auto x = u(r);
       discard(r, u, size-1);
     }
  }

  double end = omp_get_wtime();
  cout << std::fixed << std::showpoint;
  cout << setprecision(4) << std::left << std::setw(set_w) <<  __func__ <<  (end-start) << "\n";
  return end-start;
}



// *****************************************************************************************
// -----------------------------------------------------------------------------------------
// *****************************************************************************************

template<class Dist_t, Dist_e dist_e >
double test_std_mt19937_critical() {

  double start = omp_get_wtime();
  std::mt19937 r(seed);
  Dist_t u ( Params<Dist_t,dist_e>::createDist() );  // random number distribution

  #pragma omp parallel
  {
    size_t size = omp_get_num_threads();     // get total number of processes
    size_t rank = omp_get_thread_num();      // get rank of current process

    for (size_t i = rank; i < n_samples; i+=size) {
       #pragma omp critical
       auto x = u(r);
     }
  }

  double end = omp_get_wtime();
  cout << std::fixed << std::showpoint;
  cout << setprecision(4) << std::left << std::setw(set_w) <<  __func__ <<  (end-start) << "\n";
  return end-start;
}

// *****************************************************************************************
// -----------------------------------------------------------------------------------------
// *****************************************************************************************

template<class Dist_t, Dist_e dist_e >
double test_boost_mt19937_critical() {

  double start = omp_get_wtime();
  boost::mt19937 r(seed);
  Dist_t u ( Params<Dist_t,dist_e>::createDist() );  // random number distribution

  #pragma omp parallel
  {
    size_t size = omp_get_num_threads();     // get total number of processes
    size_t rank = omp_get_thread_num();      // get rank of current process

    for (size_t i = rank; i < n_samples; i+=size) {
       #pragma omp critical
       auto x = u(r);
     }
  }

  double end = omp_get_wtime();
  cout << std::fixed << std::showpoint;
  cout << setprecision(4) << std::left << std::setw(set_w) <<  __func__ <<  (end-start) << "\n";
  return end-start;
}



// *****************************************************************************************
// -----------------------------------------------------------------------------------------
// *****************************************************************************************

template<class Dist_t, Dist_e dist_e >
double test_random123() {

  double start = omp_get_wtime();

  #pragma omp parallel
  {
    size_t size = omp_get_num_threads();     // get total number of processes
    size_t rank = omp_get_thread_num();      // get rank of current process
    typedef boost::random::threefry<4, uint32_t> Prf;
    boost::random::counter_based_engine<uint32_t, Prf, 32> r(seed);

    Dist_t u ( Params<Dist_t,dist_e>::createDist() ); 

    #pragma omp for
    for (size_t i = 0; i < n_samples; ++i) {
       r.restart({i, size_t(0), size_t(0)});
       auto x = u(r);
     }
  }

  double end = omp_get_wtime();
  cout << std::fixed << std::showpoint;
  cout << setprecision(4) << std::left << std::setw(set_w) <<  __func__ <<  (end-start) << "\n";
  return end-start;
}




int main(int argc, char* argv[]) {


  cout.precision(5);
  if (argc > 1) {
    n_samples = std::stoul(std::string(argv[1]));
  }

  #pragma omp parallel
  #pragma omp single
  {
    cout << "Num of threads " << omp_get_num_threads() << endl;
    cout << "Num of samples: " << n_samples << endl << endl;
  }

  cout << "Uniform : " << "\n";
//  test_trng_yarn2             <trng::uniform_dist<double>, UNIFORM >(); 
//  test_pcg32                  <std::uniform_real_distribution<double>, UNIFORM>();
//  test_std_mt19937_leap       <std::uniform_real_distribution<double>, UNIFORM>();
//  test_boost_mt19937_leap     <boost::random::uniform_real_distribution<double>, UNIFORM>();
//  test_std_mt19937_critical   <std::uniform_real_distribution<double>, UNIFORM>();
//  test_boost_mt19937_critical <boost::random::uniform_real_distribution<double>, UNIFORM>();
  test_random123              <boost::random::uniform_real_distribution<double>, UNIFORM>();
  cout << "\n";

  cout << "Normal : " << "\n";
//  test_trng_yarn2             <trng::normal_dist<double>, NORMAL >(); 
//  test_pcg32                  <std::normal_distribution<double>, NORMAL>();
//  test_std_mt19937_leap       <std::normal_distribution<double>, NORMAL>();
//  test_boost_mt19937_leap     <boost::random::normal_distribution<double>, NORMAL>();
//  test_std_mt19937_critical   <std::normal_distribution<double>, NORMAL>();
//  test_boost_mt19937_critical <boost::random::normal_distribution<double>, NORMAL>();
  test_random123              <boost::random::normal_distribution<double>, NORMAL>();
  cout << "\n";

//  cout << "Beta : " << "\n";
//  test_trng_yarn2             <trng::uniform_dist<double>, BETA >(); 
//  test_pcg32                  <std::uniform_real_distribution<double>, BETA>();
//  test_std_mt19937_leap       <std::uniform_real_distribution<double>, BETA>();
//  test_boost_mt19937_leap     <boost::random::uniform_real_distribution<double>, BETA>();
//  test_std_mt19937_critical   <std::uniform_real_distribution<double>, BETA>();
//  test_boost_mt19937_critical <boost::random::uniform_real_distribution<double>, BETA>();
//  test_random123              <boost::random::uniform_real_distribution<double>, BETA>();
//  cout << "\n";

  cout << "Gamma : " << "\n";
//  test_trng_yarn2             <trng::gamma_dist<double>, GAMMA >(); 
//  test_pcg32                  <std::gamma_distribution<double>, GAMMA>();
//  test_std_mt19937_leap       <std::gamma_distribution<double>, GAMMA>();
//  test_boost_mt19937_leap     <boost::random::gamma_distribution<double>, GAMMA>();
//  test_std_mt19937_critical   <std::gamma_distribution<double>, GAMMA>();
//  test_boost_mt19937_critical <boost::random::gamma_distribution<double>, GAMMA>();
  test_random123              <boost::random::gamma_distribution<double>, GAMMA>();
  cout << "\n";

  cout << "Exponential : " << "\n";
//  test_trng_yarn2             <trng::exponential_dist<double>, EXPONENTIAL >(); 
//  test_pcg32                  <std::exponential_distribution<double>, EXPONENTIAL>();
//  test_std_mt19937_leap       <std::exponential_distribution<double>, EXPONENTIAL>();
//  test_boost_mt19937_leap     <boost::random::exponential_distribution<double>, EXPONENTIAL>();
//  test_std_mt19937_critical   <std::exponential_distribution<double>, EXPONENTIAL>();
//  test_boost_mt19937_critical <boost::random::exponential_distribution<double>, EXPONENTIAL>();
  test_random123              <boost::random::exponential_distribution<double>, EXPONENTIAL>();
  cout << "\n";

  cout << "Poisson : " << "\n";
//  test_trng_yarn2             <trng::poisson_dist, POISSON >(); 
//  test_pcg32                  <std::poisson_distribution<int>, POISSON>();
//  test_std_mt19937_leap       <std::poisson_distribution<int>, POISSON>();
//  test_boost_mt19937_leap     <boost::random::poisson_distribution<int>, POISSON>();
//  test_std_mt19937_critical   <std::poisson_distribution<int>, POISSON>();
//  test_boost_mt19937_critical <boost::random::poisson_distribution<int>, POISSON>();
  test_random123              <boost::random::poisson_distribution<int>, POISSON>();
  cout << "\n";

  cout << "Bernoulli : " << "\n";
//  test_trng_yarn2             <trng::bernoulli_dist<double>, BERNOULLI >(); 
//  test_pcg32                  <std::bernoulli_distribution, BERNOULLI>();
//  test_std_mt19937_leap       <std::bernoulli_distribution, BERNOULLI>();
//  test_boost_mt19937_leap     <boost::random::bernoulli_distribution<double>, BERNOULLI>();
//  test_std_mt19937_critical   <std::bernoulli_distribution, BERNOULLI>();
//  test_boost_mt19937_critical <boost::random::bernoulli_distribution<double>, BERNOULLI>();
  test_random123              <boost::random::bernoulli_distribution<double>, BERNOULLI>();
  cout << "\n";


}



