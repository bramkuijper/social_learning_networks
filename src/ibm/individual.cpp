#include "individual.hpp"

Individual::Individual() : // data member initializer list
    pp{0.0,0.0}
    ,pc{0.0,0.0}
    ,pr{0.0,0.0}
{}

Individual::Individual(Individual const &other) : // data member initializer list
    pp{other.pp[0],other.pp[1]}
    ,pc{other.pc[0],other.pc[1]}
    ,pr{other.pr[0],other.pr[1]}
{
}

void Individual::operator=(Individual const &other)
{
    for (int sex_idx = 0; sex_idx < 2; ++sex_idx)
    {
        pp[sex_idx] = other.pp[sex_idx];
        pc[sex_idx] = other.pc[sex_idx];
        pr[sex_idx] = other.pr[sex_idx];
    }
}

void Individual::Individual(
        Individual const &mother
        ,Individual const &father
        ,Parameters const &params
        ,std::mt19937 &rng)
{
    // initialize a bunch of random number distributions
    std::uniform_real_distribution<> uniform{0,1.0};

    // normal distribution to change allelic values
    std::normal_distribution<> normal{0,params.sdmu};

    // bernoulli distribution to determine whether
    // allele will be inherited from mom vs dad
    std::bernoulli_distribution from_mother{0.5};

    // inherit and mutate alleles
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
    }
} // end birth constructor
