#include <vector> // to work with vectors
#include <cassert> // to assert things - for debugging 
#include "patch.hpp"
#include "individual.hpp"



// make a patch with nf females and nm males
Patch::Patch(Parameters const &params) :
    nf{params.n[0]} // initialize number females
    ,nm{params.n[1]} // initialize number males
    ,envt2{false} // initialize local environmental state: is it 2 or 1?
    ,Wtot{0} // initialize group level fitness value
{
    // make a 'standard' individual with which we will
    // initialize the population. It has a repertoire size of
    // params.n_traits
    Individual standard_individual(
            params.n_traits
            ,params);
    
    // initialize the population of male and female breeders
    
    // make two vectors of individuals, one for females
    // make two vectors of individuals, the other for males


    std::vector < Individual > stack_females(nf, standard_individual);
    std::vector < Individual > stack_males(nm, standard_individual);

    // add the females and then the males, resulting in a 
    // 2 dimensional vector of breeders
    breeders.push_back(stack_females);
    breeders.push_back(stack_males);
    
    assert(breeders.size() == 2);

    // make the juvenile population, following the same recipe 
    // as for the female and male breeders
    std::vector<Individual> juvsf(0, standard_individual);
    std::vector<Individual> juvsm(0, standard_individual);

    // add the females and then the males, resulting in a 
    // 2 dimensional vector of juveniles 
    juveniles.push_back(juvsf);
    juveniles.push_back(juvsm);
    assert(juveniles.size() == 2);

    // make the social network which spans 
    // (nf + nm) x (nf + nm) dimensions
    //
    // make a single row first of size nf + nm
    std::vector < bool > network_row(nf + nm, false);
    // then assign that row nf + nm times to the network
    // creating a (nf + nm) x (nf + nm) dimensional array 
    // (aka a (nf+nm) dimensional square matrix)
    for (int row_idx = 0; row_idx < nf + nm; ++row_idx)
    {
        network.push_back(network_row);
    }

    // update fitness for this patch
    calculate_W();
} // end Patch::Patch()

// calculate fitness for each patch
double Patch::calculate_W()
{
    Wtot = 0;

    for (int breeder_f_idx = 0; 
            breeder_f_idx < breeders[female].size();
            ++breeder_f_idx)
    {
        Wtot += breeders[female][breeder_f_idx].update_W(envt2);
    }
    
    for (int breeder_m_idx = 0; 
            breeder_m_idx < breeders[male].size();
            ++breeder_m_idx)
    {
        Wtot += breeders[male][breeder_m_idx].update_W(envt2);
    }

    return(Wtot);
} // end Patch::calculate_W()

// make a copy constructor. This constructor makes a copy of
// an object. This is particularly useful if you store patches 
// in vectors (which we do, to create the metapopulation as a whole)
// as that requires the successive copying of a patch object
Patch::Patch(Patch const &other) :
    nf{other.nf}
    ,nm{other.nm}
    ,envt2{other.envt2}
    ,breeders{other.breeders[0],other.breeders[1]}
    ,juveniles{other.juveniles[0],other.juveniles[1]}
{
    for (int row_idx = 0; row_idx < nf + nm; ++row_idx)
    {
        network.push_back(other.network[row_idx]);
    }

    assert(network.size() == nf + nm);
} // end Patch::Patch(Patch const &other)




// overload the assignment operator so that we can assign
// one patch variable to another patch variable, without losing
// some of its contents
void Patch::operator=(Patch const &other)
{
    nf = other.nf;
    nm = other.nm;
    envt2 = other.envt2;

    for (int sex_idx = 0; sex_idx < 2; ++sex_idx)
    {
        breeders[sex_idx] = other.breeders[sex_idx];
        juveniles[sex_idx] = other.juveniles[sex_idx];
    }

    network.clear();
    
    for (int row_idx = 0; row_idx < nf + nm; ++row_idx)
    {
        network.push_back(other.network[row_idx]);
    }
    
    assert(network.size() == nf + nm);
} // end Patch::Patch(Patch const &other)

// make a random network in this patch with connection probability p
void Patch::make_random_network(double const p
        ,std::mt19937 &rng)
{
    std::uniform_real_distribution<double> uniform{0.0,1.0};
    for (int row_idx = 0; row_idx < network.size(); ++row_idx)
    {
        for (int col_idx = 0; col_idx < network.size(); ++col_idx)
        {
            network[col_idx][row_idx] = uniform(rng) < p;
        }
    }
} // end make_random_network()
