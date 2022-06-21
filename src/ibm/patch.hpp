#ifndef _PATCH_HPP
#define _PATCH_HPP

#include <vector>
#include "individual.hpp"

// blueprint of the patch class

class Patch
{
    public:
        std::vector < std::vector < Individual > > breeders;
        std::vector < std::vector < Individual > > juveniles;
        std::vector < std::vector < bool > > network;

        // number of female and male breeders
        int nf;
        int nm;


        // local environmental state
        bool envt2;

        // Patch constructor: make a patch with nf females and nm males
        Patch(Parameters const &params);
        
        // copy constructor - makes a patch out of another patch
        Patch(Patch const &other);

        void operator=(Patch const &other);
};


#endif 
