#ifndef _INDIVIDUAL_HPP_
#define _INDIVIDUAL_HPP_

#include <random>
#include "parameters.hpp"

class
{
    public:
        // learn from parent
        // first locus is mom, 2nd locus is dad
        double pp[2] = {0.0,0.0};
        
        // learn from individual connected to parent
        // first locus is female rel, 2nd locus is male rel
        double pc[2] = {0.0,0.0};
        
        // learn from random individual 
        // first locus is female, 2nd locus is male
        double pr[2] = {0.0,0.0};

        // default constructor
        Individual();

        // copy constructor
        Individual(Individual const &other);

        // assignment operator
        void operator=(Individual const &other);

        // birth constructor: 
        // - mother
        // - father
        // - parameter object containing mutation rates
        // - random number generator
        Individual(Individual const &mother
                    ,Individual const &father
                    ,Parameters const &params
                    ,std::mt19937 &rng);
};

#endif
