#ifndef _PATCH_HPP
#define _PATCH_HPP

#include <vector>
#include <random>
#include "individual.hpp"

// blueprint of the patch class

class Patch
{
    public:
        std::vector < std::vector < Individual > > breeders;
        std::vector < std::vector < Individual > > juveniles;
        std::vector < std::vector < bool > > network;

        // number of female and male breeders
        int nf = 5;
        int nm = 5;

        double Wtot = 0;


        // local environmental state
        bool envt2;

        // Patch constructor: make a patch with nf females and nm males
        Patch(Parameters const &params);
        
        // copy constructor - makes a patch out of another patch
        Patch(Patch const &other);

        void operator=(Patch const &other);

        // Calculate total fitness of the patch.
        // This depends on the environment
        // but this is a  class member variable
        double calculate_W();

        // create a random network (typically used for initialization purposes)
        void make_random_network(double const p_connect
                ,std::mt19937 &rng);
};


#endif 
