#include <random>
#include <algorithm>
#include "individual.hpp"
#include "parameters.hpp"

Individual::Individual(int const repertoire_size) : // data member initializer list
    il{0.0}
    ,pp{0.0,0.0}
    ,pc{0.0,0.0}
    ,pr{0.0,0.0}
    ,repertoire(repertoire_size,0)
    ,Wi{0.0}
{}

Individual::Individual(Individual const &other) : // data member initializer list
    il{other.il}
    ,pp{other.pp[0],other.pp[1]}
    ,pc{other.pc[0],other.pc[1]}
    ,pr{other.pr[0],other.pr[1]}
    ,repertoire(other.repertoire)
    ,Wi(other.Wi)
{
}

// assignment operator allows one to assign
// individual other to the current individual 
// while making sure
// all the data is neatly copied 
void Individual::operator=(Individual const &other)
{
    il = other.il;
    for (int sex_idx = 0; sex_idx < 2; ++sex_idx)
    {
        pp[sex_idx] = other.pp[sex_idx];
        pc[sex_idx] = other.pc[sex_idx];
        pr[sex_idx] = other.pr[sex_idx];
    }

    repertoire = other.repertoire;
    Wi = other.Wi;

} // end operator=()

Individual::Individual(
        Individual const &mother
        ,Individual const &father
        ,Parameters const &params
        ,std::mt19937 &rng) :
    Wi{0.0}
    ,repertoire(params.n_traits,0) // start with empty repertoire, of size n_traits
{
    // initialize a bunch of random number distributions
    std::uniform_real_distribution<> uniform{0,1.0};

    // normal distribution to change allelic values
    std::normal_distribution<> normal{0,params.sdmu};

    // bernoulli distribution to determine whether
    // allele will be inherited from mom vs dad
    std::bernoulli_distribution from_mother{0.5};
        

    // inherit allele from one or the other
    il = from_mother(rng) ? 
        mother.il : father.il;

    // if there is a mutation, update the allelic value
    if (uniform(rng) < params.mu_il)
    {
        il += normal(rng);
    }

    // restrict value between 0 and 1
    il = std::clamp(il, 0.0, 1.0);



    // inherit and mutate alleles  - traits are sex
    // specific (learning from male vs female)
    for (int sex_idx = 0; sex_idx < 2; ++sex_idx)
    {


        // inherit allele from mom or dad
        pp[sex_idx] = from_mother(rng) ? 
            mother.pp[sex_idx] : father.pp[sex_idx];

        // if there is a mutation, update the allelic value
        if (uniform(rng) < params.mu_pp)
        {
            pp[sex_idx] += normal(rng);
        }

        // restrict value between 0 and 1
        pp[sex_idx] = std::clamp(pp[sex_idx], 0.0, 1.0);


        // do the same re the other traits...
        pc[sex_idx] = from_mother(rng) ? 
            mother.pc[sex_idx] : father.pc[sex_idx];

        if (uniform(rng) < params.mu_pc)
        {
            pc[sex_idx] += normal(rng);
        }

        pc[sex_idx] = std::clamp(pc[sex_idx], 0.0, 1.0);



        pr[sex_idx] = from_mother(rng) ? 
            mother.pr[sex_idx] : father.pr[sex_idx];

        if (uniform(rng) < params.mu_pr)
        {
            pr[sex_idx] += normal(rng);
        }

        pr[sex_idx] = std::clamp(pr[sex_idx], 0.0, 1.0);
    } // end for int sex_idx
} // end birth constructor

// obtain fitness
double Individual::update_W(bool const is_envt2)
{
    if (is_envt2)
    {
        int ntraits = 0;
        for (int trait_idx = 0; 
                trait_idx < repertoire.size(); ++trait_idx)
        {
            if (repertoire[trait_idx] > 0)
            {
                ++ntraits;
            }
        }

        return (ntraits);
    }

    // this is never reached unless envt = 1

    // specialist envt: find iterator pointing to 
    // max element in list
    auto iterator_of_max_elmt = std::max_element(
            repertoire.begin()
            ,repertoire.end());

    // return value of max element
    return(*iterator_of_max_elmt);
} // update_W
