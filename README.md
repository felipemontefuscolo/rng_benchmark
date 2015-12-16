Random Number Generators benchmark
==================================



Dependencies:

  - g++ C++11 compiler
  - PAPI (version >= 3)
  - OpenMP
  - Boost

The tested libraries are:

  - [Trng](http://numbercrunch.de/trng/)
  - [PCG](http://www.pcg-random.org/)
  - [C++11 Standard Library](http://en.cppreference.com/w/cpp/numeric/random)
  - [Boost.Random](http://www.boost.org/doc/libs/1_59_0/doc/html/boost_random.html)
  - [Random123](https://github.com/DEShawResearch/Random123-Boost)

In the file `parameters.cfg` you can choose the seed and the parameters of the distributions. If `test_reproducibility` is 1, the code will print the files:

    Bernoulli.dat  Exponential.dat  Gamma.dat  Normal.dat  Poisson.dat  Uniform.dat

Each file contains the corresponding distribution for each library. It is written row-wise and NOT col-wise.

You can use octave to check the distribution, for instance:

    m = dlmread('Normal.dat')
    hist(m(1,:))
    hist(m(2,:))
    hist(m(3,:))
    ...

The order obey the order you see in the code output.

Please Edit the Makefile and adapt these parameters to your needs.






