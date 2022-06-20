#ifndef _PATCH_HPP
#define _PATCH_HPP

#include "individual.hpp"

// blueprint of the patch class

class Patch
{
    private:
        std::vector < std::vector < Individual > > breeders;
        std::vector < std::vector < Individual > > juveniles;

    public:
        // Patch constructor: make a patch with nf females and nm males
        Patch(int const nf
                ,int const nm);
        
        // copy constructor - makes a patch out of another patch
        Patch(Patch const &other);

        void operator=(Patch const &other);
};


#endif 
