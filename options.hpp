// Copyright Vladimir Prus 2002-2004.
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

/* Shows how to use both command line and config file. */

#ifndef BENCH_OPTIONS_HPP
#define BENCH_OPTIONS_HPP

#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
namespace po = boost::program_options;

#include <sstream>
#include <iostream>
#include <fstream>
#include <iterator>
#include <tuple>
using namespace std;

using std::cout;
using std::endl;
using std::setprecision;

// A helper function to simplify the main part.
template<class T>
ostream& operator<<(ostream& os, const vector<T>& v)
{
    copy(v.begin(), v.end(), ostream_iterator<T>(os, " ")); 
    return os;
}
template<class T, class U>
ostream& operator<<(ostream& os, const tuple<T,U>& v)
{
    os << get<0>(v) << " " << get<1>(v);
    return os;
}
// end helper function


struct Parameters {

  typedef tuple<double,double> pair_t;

  bool test_reproducibility = false;
  size_t num_samples = 2000000 ;
  vector<int> num_threads = vector<int>{1, 4, 8, 16, 32} ;
  unsigned long random_seed = 23948u;
  pair_t uniform_params = pair_t{0., 1.0};
  pair_t normal_params = pair_t{0., 1.0};
  pair_t beta_params = pair_t{1.0, 2.0};
  pair_t gamma_params = pair_t{1.0, 2.0};
  double exponential_param = 1.0;
  int poisson_param = 4;
  double bernoulli_param = 0.5;

  void print_values() const {

    cout << "This code is running with the following options:\n\n";

    cout << "  test_reproducibility:  " <<  test_reproducibility << endl;
    cout << "  num_samples:           " <<  num_samples       << endl;   
    cout << "  num_threads:           " <<  num_threads       << endl;   
    cout << "  random_seed:           " <<  random_seed       << endl;  
    cout << "  uniform_params:        " <<  uniform_params    << endl;    
    cout << "  normal_params:         " <<  normal_params     << endl;      
    cout << "  beta_params:           " <<  beta_params       << endl;  
    cout << "  gamma_params:          " <<  gamma_params      << endl;  
    cout << "  exponential_param:     " <<  exponential_param << endl;  
    cout << "  poisson_param:         " <<  poisson_param     << endl;  
    cout << "  bernoulli_param:       " <<  bernoulli_param   << endl; 

    cout << endl;
  }


};

template<class T>
struct VecT : public vector<T>
{
  using vector<T>::vector;
};
// This is for reading vectors with boost::program_options
template<class T>
void validate(boost::any& v, const std::vector<std::string>& values, VecT<T>*, int)
{
  VecT<T> dvalues;
  for(vector<string>::const_iterator it = values.begin(); it != values.end(); ++it) {
    stringstream ss(*it);
    copy(istream_iterator<T>(ss), istream_iterator<T>(), back_inserter(dvalues));
    if(!ss.eof()) {
      throw ("Invalid parameter(s) specification");
    }
  }
  v = dvalues;
}
// This is for reading tuples2 with boost::program_options
template<class U>
struct TupleT : public tuple<U,U>
{
  using tuple<U,U>::tuple;
};
template<class T>
void validate(boost::any& v, const std::vector<std::string>& values, TupleT<T>*, int) 
{
  VecT<T> dvalues;
  for(vector<string>::const_iterator it = values.begin(); it != values.end(); ++it) {
    stringstream ss(*it);
    copy(istream_iterator<T>(ss), istream_iterator<T>(), back_inserter(dvalues));
    if(!ss.eof()) {
      throw ("Invalid parameter(s) specification");
    }
  }
  v = TupleT<T>{dvalues.at(0), dvalues.at(1)};
}

int readParameters(Parameters &opt, int ac, char* av[]) {
    try {
        string config_file;
    
        // Declare a group of options that will be 
        // allowed only on command line
        po::options_description generic("Generic options");
        generic.add_options()
            ("version,v", "print version string")
            ("help", "produce help message")
            ("config,c", po::value<string>(&config_file)->default_value("parameters.cfg"),
                  "name of a file of a configuration.")
            ;
    
        // Declare a group of options that will be 
        // allowed both on command line and in
        // config file
        po::options_description config("Configuration");
        config.add_options()
            ("include-path,I",    po::value<vector<string>>()->composing(),  "include path")
            ("test_reproducibility", po::value<bool       >(), "check if the result is thread independent (slow)")      
            ("num_samples",       po::value<size_t        >(), "number of random numbers generated")      
            ("num_threads",       po::value<VecT<int>     >()->multitoken(), "list of number of threads separated by space")        
            ("random_seed",       po::value<unsigned long >(), "seed")           
            ("uniform_params",    po::value<TupleT<double>>()->multitoken(), "two floats separated by space")        
            ("normal_params",     po::value<TupleT<double>>()->multitoken(), "two floats separated by space")      
            ("beta_params",       po::value<TupleT<double>>()->multitoken(), "two floats separated by space")          
            ("gamma_params",      po::value<TupleT<double>>()->multitoken(), "two floats separated by space")          
            ("exponential_param", po::value<double        >(), "a float")         
            ("poisson_param",     po::value<int           >(), "an integer")         
            ("bernoulli_param",   po::value<double        >(), "a float")  
            ;

        // Hidden options, will be allowed both on command line and
        // in config file, but will not be shown to the user.
        po::options_description hidden("Hidden options");
        hidden.add_options()
            ("input-file", po::value< vector<string> >(), "input file")
            ;

        
        po::options_description cmdline_options;
        cmdline_options.add(generic).add(config).add(hidden);

        po::options_description config_file_options;
        config_file_options.add(config).add(hidden);

        po::options_description visible("Allowed options");
        visible.add(generic).add(config);
        
        po::positional_options_description p;
        p.add("input-file", -1);
        
        po::variables_map vm;
        store(po::command_line_parser(ac, av).
              options(cmdline_options).positional(p).run(), vm);
        notify(vm);
        
        ifstream ifs(config_file.c_str());
        if (!ifs)
        {
            cout << "can not open config file: " << config_file << "\n";
            throw;
        }
        else
        {
            store(parse_config_file(ifs, config_file_options), vm);
            notify(vm);
        }
    
        if (vm.count("help")) {
            cout << visible << "\n";
            return 0;
        }

        if (vm.count("version")) {
            cout << "Multiple sources example, version 1.0\n";
            return 0;
        }

        if (vm.count("include-path"))
        {
            vector<string> ss = vm["include-path"].as< vector<string> >();
            cout << "Include paths are: " 
                 //<< vm["include-path"].as< vector<string> >() << "\n";
                 << ss << "\n";
        }

        if (vm.count("input-file"))
        {
            cout << "Input files are: " 
                 << vm["input-file"].as< vector<string> >() << "\n";
        }

       if (vm.count("test_reproducibility")) opt.test_reproducibility = vm["test_reproducibility"].as<bool>();      
       if (vm.count("num_samples"      )) opt.num_samples       = vm["num_samples"      ].as< size_t        >();     
       if (vm.count("num_threads"      )) opt.num_threads       = vm["num_threads"      ].as< VecT<int>     >();     
       if (vm.count("random_seed"      )) opt.random_seed       = vm["random_seed"      ].as< unsigned long >();     
       if (vm.count("uniform_params"   )) opt.uniform_params    = vm["uniform_params"   ].as< TupleT<double>>();      
       if (vm.count("normal_params"    )) opt.normal_params     = vm["normal_params"    ].as< TupleT<double>>();     
       if (vm.count("beta_params"      )) opt.beta_params       = vm["beta_params"      ].as< TupleT<double>>();   
       if (vm.count("gamma_params"     )) opt.gamma_params      = vm["gamma_params"     ].as< TupleT<double>>();    
       if (vm.count("exponential_param")) opt.exponential_param = vm["exponential_param"].as< double        >();       
       if (vm.count("poisson_param"    )) opt.poisson_param     = vm["poisson_param"    ].as< int           >();      
       if (vm.count("bernoulli_param"  )) opt.bernoulli_param   = vm["bernoulli_param"  ].as< double        >();         

    }
    catch(exception& e)
    {
        cout << e.what() << " fuck \n";
        return 1;
    }    

    return 0;
}




#endif
