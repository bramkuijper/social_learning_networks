#include "patch.hpp"
#include "individual.hpp"

// make a patch with nf females and nm males
Patch::Patch(
        int const nf
        ,int const nm) :
{
    // initialize the population of male and female breeders
    
    // make two vectors of individuals, one for females
    // make two vectors of individuals, the other for males
    std::vector < Individual > stack_females(nf, Individual());
    std::vector < Individual > stack_males(nm, Individual());

    // add the females and then the males, resulting in a 
    // 2 dimensional vector of breeders
    breeders.push_back(stack_females);
    breeders.push_back(stack_males);
    
    assert(breeders.size() == 2);

    // make the juvenile population, following the same recipe 
    // as for the female and male breeders
    std::vector<Individual> juvsf(0, Individual());
    std::vector<Individual> juvsm(0, Individual());

    juveniles.push_back(juvsf, juvsm);
    assert(juveniles.size() == 2);

} // end Patch::Patch()
 
// make a copy constructor. This constructor makes a copy of
// an object. This is particularly useful if you store patches 
// in vectors (which we do, to create the metapopulation as a whole)
// as that requires the successive copying of a patch object
Patch::Patch(Patch const &other) :
    breeders{other.breeders[0],other.breeders[1]}
    ,juveniles{other.juveniles[0],other.juveniles[1]}
{
} // end Patch::Patch(Patch const &other)

// overload the assignment operator so that we can assign
// one patch variable to another patch variable, without losing
// some of its contents
void Patch::operator=(Patch const &other)
{
    for (int sex_idx = 0; sex_idx < 2; ++sex_idx)
    {
        breeders[sex_idx] = other.breeders[sex_idx];
        juveniles[sex_idx] = other.juveniles[sex_idx];
    }
} // end Patch::Patch(Patch const &other)
