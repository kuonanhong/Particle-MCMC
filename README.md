# SSM
There are two main purposes for this code. One is to provide base classes for different filters (e.g. the Kalman Filter, Sequential Importance Sampling with Resampling (SISR), the Auxiliary Particle Filter (APF), etc.). The other is to provide base classes for particle Markov chain Monte Carlo algorithms. This second category of base classes make use of the first category. Please see the different examples for demonstrations.

## Installation
You have to build this yourself. Make sure to compile with C++11 enabled (`-std=c++11`), to include the linker option `-lpthread`, and to include the `include` directory. Note, also, that this code all makes use of the [Eigen library](http://eigen.tuxfamily.org/).

## Organization of `src`
1. `distributions` is all the code related to evaluating densities, mass functions, or pseudo-random number generation.
2. `examples` are all of the self-contained examples that implement specific models and algorithms.
3. `mcmc_algos` are classes performing MCMC techniques for specific models. These classes inherit from the classes in `mcmc_bases`.
4. `mcmc_bases` are all of the base classes for certain types of MCMC algorithms. Inheriting from these allows the user to skip writing a bunch of code.
5. `models` are particle filtering classes for specific state space models.
6. `filter_bases` are Kalman filtering and particle filtering base classes. Inheriting from these abstracts away all of the details of how this is performed. Particle filter base classes are pure virtual, while the closed form filtering classes (fshmm and lgssm) are not.
7. `utilities` functions that might be useful and have no other home (e.g. reading or logging data, transforming parameters with common functions, etc.).

## Documentation
[More details on everything can be found here.](https://tbrown122387.github.io/ssm/)

## To-Do List
- Make a parameter class that organizes which are transformed, constrained, etc.
- Implement more resampling options (different algorithms, different schedules, etc.)
- Make sure all resamplers and classes are log-weights only
- unit testing...
