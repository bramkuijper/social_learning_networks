#ifndef _INDIVIDUAL_HPP_
#define _INDIVIDUAL_HPP_

#include <random>
#include "parameters.hpp"

// the individual object
// containing the traits and states 
// of an individual

class Individual
{
    public:
        // prob to learn as an individual learner
        // social learning is 1 - il
        double il = 0.0;

        // prob to make network connection with parent
        // first locus is mom, 2nd locus is dad
        double pp[2] = {0.0,0.0};
        
        // prob to make network connection with non-parent relative
        // first locus is female rel, 2nd locus is male rel
        double pc[2] = {0.0,0.0};
        
        // prob to make network connection with random individual
        // first locus is female, 2nd locus is male
        double pr[2] = {0.0,0.0};

        // variable specifying individual fitness
        double Wi = 1.0;
       
        // the repertoire as in Smolla & Akcay (see p8)
        std::vector < int > repertoire;

        // default constructor
        Individual(int const repertoire_size, Parameters const &params);

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

        double update_W(bool const envt_is_2);
};

#endif
