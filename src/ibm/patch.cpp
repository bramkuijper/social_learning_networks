#include "patch.hpp"
#include "individual.hpp"

        // make a patch with nf females and nm males
Patch::Patch(
        int const nf
        ,int const nm) :
    
{
    // make two vectors of individuals, one for females
    // make two vectors of individuals, the other for males
    std::vector<Individual> stack_females(nf, Individual());
    std::vector<Individual> stack_males(nm, Individual());

    breeders.push_back(stack_females);
    breeders.push_back(stack_males);
    
    assert(breeders.size() == 2);
    
    std::vector<Individual> juvsf(0, Individual());
    std::vector<Individual> juvsm(0, Individual());

    juveniles.push_back(juvsf, juvsm);

    assert(juveniles.size() == 2);
} // Patch::Patch()
