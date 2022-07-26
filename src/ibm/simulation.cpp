// the heart of the code:  Simulation::functions describing the 
// different parts of the life cycle



#include <random> // random number generating libraries
#include <vector>
#include <cassert>
#include <numeric>
#include "simulation.hpp"
#include "patch.hpp"
#include "individual.hpp"
#include "parameters.hpp"

// the simulation constructor which 
// initializes everything
Simulation::Simulation(Parameters const &par) :
    rd{} // initialize random device, see *.hpp file
    ,seed{rd()} // initialize seed
    ,rng_r{seed} // initialize the random number generator
    ,uniform{0.0,1.0} // initialize the uniform distribution
    ,patch_sampler{0,(int)par.n_patches - 1} // initialize the patch sampling distribution
    ,data_file{par.base_name.c_str()} // initialize the data file by giving it a name
    ,data_file_network{par.base_name_matrix_file.c_str()} // initialize network data file
                                                          // by giving it a name
    ,par{par} // initialize the parameter data member with the constructor argument
    ,metapop{par.n_patches, Patch(par)} // initialize a meta population each with n1 individuals of species 1 and n2 individuals of species 2
    ,W_global_total{0.0} // set total fitness to 0
{
    // initialize patch environments
    //
    // count the number of patches in environmental state 2
    // (number of patches in environmental state 1 = par.n_patches - n_patches_2)
    n_patches_2 = 0;

    // overall frequency of environment 2 is 
    // switch_rate[1->2]/(switch_rate[1->2] + switch_rate[2->1])
    double prob_envt2 = par.switch_rate[1]/(par.switch_rate[1] + 
            par.switch_rate[0]);

    
    // set the environment for each patch
    // initialize the network for each patch
    // as a random network
    for (int patch_idx = 0; patch_idx < metapop.size(); ++patch_idx)
    {
        // randomly set the environment as a function of the probability
        // of getting environment 2. The value can either be true (envt2)
        // or false (envt1)
        metapop[patch_idx].envt2 = uniform(rng_r) < prob_envt2;

        n_patches_2 += metapop[patch_idx].envt2;

        // create a random network in each patch
        metapop[patch_idx].make_random_network(
                par.p_network_init
                ,rng_r);

        W_global_total += metapop[patch_idx].calculate_W();
    }

} // end IBM_Mutualism constructor

// run the simulation
void Simulation::run()
{
    // write the headers to the output file
    // so that all columns have names in the output file
    write_data_headers();

    // we loop over all the time steps and perform a whole life cycle
    // each and every time step
    for (time_step = 0; time_step < par.max_time_steps; ++time_step)
    {
        // total rate of environmental change
        // is:
        // number of envt 1 patches * switch_rate[0]
        // number of envt 2 patches * switch_rate[1]
        double rate_of_change[2] = {
            ((int)metapop.size() - n_patches_2) * par.switch_rate[0] 
            ,n_patches_2 * par.switch_rate[1]
        };

        // the individual death rate is set 1 (i.e., a standard rate
        // relative to all other rates of things happening), which means
        // that the overall rate at which some dies 
        // in the population is simply N
        double death_rate = metapop.size() * (par.n[female] + par.n[male]);

        // make a probability distribution of 
        // the different events
        std::discrete_distribution<int> event_chooser{
            rate_of_change[0]
            ,rate_of_change[1]
            ,death_rate};

        // sample from the event chooser
        // and then perform actions according to what
        // was sampled
        switch(event_chooser(rng_r))
        {
            case 0:
                environmental_change(0); // a patch in envtal state 1 changes its envt
                break;
            case 1:
                environmental_change(1); // a patch in envtal state 2 changes its envt
                break;
            case 2:
                death_birth(); // a random individual dies and is replaced
                break;
            default:
                throw std::range_error("Something is going wrong in the event chooser. There are more event options than currently accommodated for.");
                break;
        } // end switch
   
        // output stats to file every n time steps
        if (time_step % par.data_output_interval == 0)
        {
            std::cout << "time step: " << time_step << std::endl;
            write_data();
        }
    } // end for time_step

    // write the networks to a file
    write_all_networks();

    // write parameters to output
    write_parameters();
} // end Simulation::run()

// change the environment of a local pod
// argument specifies current environmental state of the patch
// that needs to be changed
void Simulation::environmental_change(bool const is_envt2)
{
    // auxiliary variable to contain the fitness
    // of the current patch
    double W_this_patch;

    // sample a random patch that will change
    int random_patch = patch_sampler(rng_r);

    for (int attempt_idx = 0; 
            attempt_idx < par.max_search_attempts; 
            ++attempt_idx)
    {
        // if one finds an environment that matches the local 
        // patch change the local environment
        if (metapop[random_patch].envt2 == is_envt2)
        {
            // get current fitness on this patch
            W_this_patch = metapop[random_patch].calculate_W();

            // subtract this from the local sum
            W_global_total -= W_this_patch;

            metapop[random_patch].envt2 = !metapop[random_patch].envt2;

            // new environment is 2
            //
            // update the counter of the patches in different environmental states
            if (metapop[random_patch].envt2)
            {
                ++n_patches_2;
            }
            else
            {
                --n_patches_2;
            }

            // calculate fitness for this changed patch - 
            // as the environment
            // has now changed its previous fitness value does not 
            // make much sense anymore
            // hence we recalculate
            W_this_patch = metapop[random_patch].calculate_W();

            W_global_total += W_this_patch;

            break;
        }
        else // ok the sampled patch was not in the desired state, sample another
        {
            random_patch = patch_sampler(rng_r);
        }
    } // end for
} // end void Simulation::environmental_change


// write parameters to a file
void Simulation::write_parameters()
{
    data_file << std::endl
        << std::endl;

    std::string sex_str;

    for (int sex_idx = 0; sex_idx < 2; ++sex_idx)
    {
        sex_str = sex_idx == female ? "f" : "m";
        data_file << "d" << sex_str << ";" << par.d[sex_idx] << std::endl;
        data_file << "n" << sex_str << ";" << par.n[sex_idx] << std::endl;
    }

    for (int envt_idx = 0; envt_idx < 2; ++envt_idx)
    {
        data_file << "switch_rate" 
            << (envt_idx + 1) << ";" 
            << par.switch_rate[envt_idx] << std::endl;
    }

    data_file << "random_seed;" << seed << std::endl;
    data_file << "n_traits;" << par.n_traits << std::endl;
    data_file << "n_patches;" << par.n_patches << std::endl;
    data_file << "n_patches;" << par.n_patches << std::endl;
    data_file << "n_learning_attempts;" << par.n_learning_attempts << std::endl;

    // see eq. (1) in Smolla & Akcay
    data_file << "gamma;" << par.gamma << std::endl;
    // see eq. (2) in Smolla & Akcay
    data_file << "sigma;" << par.sigma << std::endl;

    data_file << "mu_il;" << par.mu_il << std::endl;
    data_file << "mu_pp;" << par.mu_pp << std::endl;

    data_file << "mu_pc;" << par.mu_pc << std::endl;
    data_file << "mu_pr;" << par.mu_pr << std::endl;
} // Simulation::write_parameters()

// death followed by a birth
void Simulation::death_birth()
{
    // pick random patch to select individual to die
    //
    // effectively there are n_patches * (nf + nm) individuals
    int random_patch_idx = patch_sampler(rng_r);

    // OK this patch will change - let's subtract its fitness from 
    // the global fitness count
    double W_this_patch = metapop[random_patch_idx].calculate_W();

    W_global_total -= W_this_patch;

    // whether it is a male or a female that dies
    // if 0 (false) female, if 1 (true) male - see also 
    // the Sex enum in parameters.hpp
    Sex is_male = (Sex)(uniform(rng_r) < (double)par.n[male] / (par.n[female] + par.n[male]));

    // make a distribution to randomly sample the individual
    // of a particular sex who will die 
    std::uniform_int_distribution<int> 
        to_die_sampler(0, metapop[random_patch_idx].breeders[is_male].size() - 1);

    // sample the individual
    // that will die on this random patch
    int individual_idx = to_die_sampler(rng_r);

    // ask where ths individual resides in the network matrix
    int individual_location_in_network = individual_idx;

    // if this individual i is a male, note that its index
    // in the network is then nf + i
    // this is because the network is a (nf + nm) x (nf + nm) 
    // matrix (i.e., males have rows below those of females)
    if (is_male)
    {
        individual_location_in_network += par.n[female];
    }

    // get number of columns in network matrix
    // and store it in variable to speed things up
    // when we will loop through the matrix
    int ncols = metapop[random_patch_idx].network[
        individual_location_in_network].size();

    // bounds checking of dying individual's location in network
    assert(individual_location_in_network >= 0);
    assert(individual_location_in_network < ncols);
    assert(individual_location_in_network < metapop[random_patch_idx].nf 
            + metapop[random_patch_idx].nm);

    // set everything in the social learning network 
    // relating to this individual to 0 - as the indiv
    for (int col_idx = 0; col_idx < ncols; ++col_idx)
    {
        // set all connections in which 
        // individual serves as learner
        // to 0
        metapop[random_patch_idx].network[
            individual_location_in_network][col_idx] = false;

        // set all connections in which 
        // individual serves as model
        // to 0
        metapop[random_patch_idx].network[col_idx][
            individual_location_in_network] = false;
    }

    make_new_individual(
            random_patch_idx 
            ,is_male
            ,individual_idx
            );

    // add fitness of this newly changed patch to the global sum
    W_this_patch = metapop[random_patch_idx].calculate_W();

    W_global_total += W_this_patch;
} // end Simulation::death_birth()

void Simulation::make_new_individual(
        int const local_patch_idx // index of patch where offspring will establish as breeder
        ,bool const offspring_is_male // offspring will be male or female, this depends on who has died
        ,int const individual_idx // index in vector of breeders that this individual will take up
        )
{
    assert(local_patch_idx >= 0);
    assert(local_patch_idx < metapop.size());
    
    assert(individual_idx >= 0);
    assert(individual_idx < metapop[local_patch_idx].breeders[offspring_is_male].size());

    // fitness of the local patch
    double Wloc = metapop[local_patch_idx].Wtot;

    // to calculate probability that vacant breeding position 
    // will be taken my immigrant individual
    // calculate productivity of all remote patches (including the local one, 
    // which does not matter too much if the total # patches is large)
    // and more realistic as by accident individual can disperse and
    // still end up in original patch (as is also reflected in random
    // patch sampling below).
    double Wremote = W_global_total;

    // divide through patches, as 
    // any dispersing offspring may land on any of Npatches patches.
    // Hence probability of landing on this patch is 1/Npatches
    Wremote /= metapop.size();

    // by default, sample from local patch, unless..
    int patch_origin_idx = local_patch_idx;

    // ... we sample a remotely born offspring with a 
    // probability equal to  
    // the average productivity of a remote patch
    // times the dispersal rate 
    if (uniform(rng_r) < par.d[offspring_is_male] * Wremote / 
            (par.d[offspring_is_male] * Wremote + (1.0 - par.d[offspring_is_male]) * Wloc))
    {
        patch_origin_idx = patch_sampler(rng_r);
    }

    // now sample parents according to their fitnesses within that patch
    // and make offspring. This may include the 'dead' parent, implying that the
    // parent only dies after reproduction (cf. Smolla & Akcay)
    std::vector<double> Wifs;

    // set up the sampling distribution values of productivity for males and females
    for (int female_idx = 0; female_idx < par.n[female]; ++female_idx)
    {
        Wifs.push_back(metapop[patch_origin_idx].breeders[female][female_idx].Wi);
    }
    
    // vector with male fitness
    std::vector<double> Wims;

    for (int male_idx = 0; male_idx < par.n[male]; ++male_idx)
    {
        Wims.push_back(metapop[patch_origin_idx].breeders[male][male_idx].Wi);
    }


    // make sampling distributions of female and male fitness values
    // so that those individuals with larger values are more likely to be
    // chosen from the ditsitribution
    std::discrete_distribution<int> female_sampler(Wifs.begin(), Wifs.end());
    std::discrete_distribution<int> male_sampler(Wims.begin(), Wims.end());

    int mother_idx = female_sampler(rng_r);
    int father_idx = male_sampler(rng_r);

    // bounds checking
    assert(mother_idx >= 0);
    assert(mother_idx < par.n[female]);
    assert(mother_idx < metapop[patch_origin_idx].breeders[female].size());

    assert(father_idx >= 0);
    assert(father_idx < par.n[male]);
    assert(father_idx < metapop[patch_origin_idx].breeders[male].size());

    // use the birth constructor  in individual.cpp
    // to produce a new offspring
    Individual offspring(
        metapop[patch_origin_idx].breeders[female][mother_idx] // mom
        ,metapop[patch_origin_idx].breeders[male][father_idx] // dad
        ,par // Parameter object to get at mutation rates
        ,rng_r);


    // put offspring in place of the now dead breeder
    // and then perform learning
    metapop[local_patch_idx].
        breeders[offspring_is_male][individual_idx] = offspring;

    // make the network 
    // for this individual
    generate_network(
        local_patch_idx
        ,patch_origin_idx
        ,individual_idx
        ,(Sex)offspring_is_male
        ,mother_idx
        ,father_idx
        );

    // now perform learning
    learn(
        local_patch_idx // patch where this individual eventually lives
        ,patch_origin_idx  // patch where individual is born (i.e., where its parents are)
        ,(Sex)offspring_is_male // whether individual is male or female
        ,individual_idx // its index in the vector of breeders
        ,mother_idx // the idx of the mom
        ,father_idx // the idx of the dad
        );
}// end Simulation::sample_parents()

void Simulation::generate_network(
        int const local_patch_idx
        ,int const patch_origin_idx
        ,int const offspring_idx
        ,Sex const offspring_sex
        ,int const mother_idx
        ,int const father_idx
        )
{
    // bounds checking - once we set NDEBUG this is all turned off
    assert(local_patch_idx >= 0);
    assert(local_patch_idx < metapop.size());
    
    assert(patch_origin_idx >= 0);
    assert(patch_origin_idx < metapop.size());
    
    assert(offspring_idx >= 0);
    assert(offspring_idx < metapop[local_patch_idx].breeders[offspring_sex].size());

    // helper variable informing us about offspring's place in network matrix
    int individual_network_idx = offspring_idx;

    // as male columns and rows in the network matrix are placed 
    // below those of the female columns and rows, increase index
    // if this is a male, by the total columns/rows allocated to females
    if (offspring_sex == male)
    {
        individual_network_idx += par.n[female];
    }

    // similarly, make a help variable that gets dad's 
    // location in the network matrix
    int father_network_idx = father_idx + par.n[female];
    
    // helper variable informing one what total rows and columns 
    // are in network matrix
    int total_rows_cols_matrix = par.n[female] + par.n[male];

    // if debugging is on, check whether network entries 
    // have been set to 0 for this individual
    // when its predecessor died
#ifndef NDEBUG    
    for (int col_idx = 0; 
            col_idx < total_rows_cols_matrix; ++col_idx)
    {
        assert(!metapop[local_patch_idx].network[individual_network_idx][col_idx]);
        assert(!metapop[local_patch_idx].network[col_idx][individual_network_idx]);
    }
#endif


    // make a copy of the individual - easier to lookup trait values
    // rather than having to use umpteen array indices
    Individual focal = metapop[local_patch_idx].breeders[offspring_sex][offspring_idx];


    // individual always learns from the patch on arrival for now
    //
    // hence only learns from parents or individuals connected to parents
    // when parents are present on the local patch.
    if (local_patch_idx == patch_origin_idx)
    {
        // with probability pp[female] make a connection to mom
        if (uniform(rng_r) < focal.pp[female])
        {
            metapop[local_patch_idx].network[individual_network_idx][mother_idx] = true;
        }

        // with probability pp[male] make a connection to dad
        if (uniform(rng_r) < focal.pp[male])
        {
            metapop[local_patch_idx].
                network[individual_network_idx][father_network_idx] = true;
        }
    } // end if (local_patch_idx == patch_origin_idx)

    // with probability pc[female] make a connection to all the individuals
    // in mom's network
    // with probability pc[male] make a connection to all the individuals
    // in dad's network
    for (int col_idx = 0; 
            col_idx < total_rows_cols_matrix; ++col_idx)
    {
        // skip mother, father and self 
        if ((col_idx == individual_network_idx || 
                col_idx == mother_idx) ||
                col_idx == father_network_idx
                )
        {
            continue;
        }

        // if we are among mom's community (ie., patch of origin
        // is the patch of arrival) and mom has connection to other individual
        // then potentially copy that connection
        if (local_patch_idx == patch_origin_idx 
                && metapop[local_patch_idx].network[mother_idx][col_idx])
        {
            // copy it with probability pc[female]
            if (uniform(rng_r) < focal.pc[female])
            {
                // OK copy it to offspring
                metapop[local_patch_idx].
                    network[individual_network_idx][col_idx] = true;
            }
        } // if dad has connection to other individual
        else if (local_patch_idx == patch_origin_idx
                && metapop[local_patch_idx].network[father_network_idx][col_idx])
        {
            // copy it with probability pc[male]
            if (uniform(rng_r) < focal.pc[male])
            {
                // OK copy it to offspring
                metapop[local_patch_idx].
                    network[individual_network_idx][col_idx] = true;
            }
        } 
        else // other ind has neither connection to mom or dad, 
        {    // hence see whether we copy randomly
            
            // first determine sex of potential model
            Sex other_ind_sex = col_idx >= par.n[female] ? 
                male
                :
                female;
            
            // then evaluate probability of establishing connection
            // with that individual of that sex
            if (uniform(rng_r) < focal.pr[other_ind_sex])
            {
                metapop[local_patch_idx].
                    network[individual_network_idx][col_idx] = true;
            }
        } 
    }
} // end Simulation::generate_network


// writes the headers to the data files
void Simulation::write_data_headers()
{
    // headers of the network file
    data_file_network << "time_step;patch;from;to;" << std::endl;

    // headers of the normal data file
    data_file << "time_step;mean_il;var_il;";

    std::string sex_abbr[2] = {"f","m"};

    for (int sex_idx = 0; sex_idx < 2; ++sex_idx)
    {
        data_file << "mean_pp" << sex_abbr[sex_idx] 
            << ";var_pp" << sex_abbr[sex_idx] << ";";

        data_file << "mean_pc" << sex_abbr[sex_idx] 
            << ";var_pc" << sex_abbr[sex_idx] << ";";
        
        data_file << "mean_pr" << sex_abbr[sex_idx] 
            << ";var_pr" << sex_abbr[sex_idx] << ";";
    } // for (int sex_idx = 0; sex_idx < 2; ++sex_idx)

    data_file << "W_global_total;mean_repertoire_size;mean_pi_ts;" << std::endl;
} // end Simulation::write_data_headers()

// actually learn from others 
// and increase repertoire
void Simulation::learn(
    int const local_patch_idx // patch where this individual eventually lives
    ,int const patch_of_origin_idx  // patch where individual is born (i.e., where its parents are)
    ,Sex const offspring_sex // whether individual is male or female
    ,int const individual_idx // its index in the vector of breeders
    ,int const mother_id // the idx of the mom
    ,int const father_id // the idx of the dad
    )
{
    // individual should have an empty reportoire
    // check whether that is indeed the case
    int sum_rep = std::accumulate(
            metapop[local_patch_idx].breeders[offspring_sex][individual_idx].repertoire.begin()
            ,metapop[local_patch_idx].breeders[offspring_sex][individual_idx].repertoire.end()
            ,0);

    assert(sum_rep == 0);

    
    // go through all traits of this individual's network
    // and make the pi distribution
    int max_network_size = par.n[female] + par.n[male];

    // helper variable informing us about offspring's place in network matrix
    int individual_network_idx = individual_idx;

    // as male columns and rows in the network matrix are placed 
    // below those of the female columns and rows, increase index
    // if this is a male, by the total columns/rows allocated to females
    if (offspring_sex == male)
    {
        individual_network_idx += par.n[female];
    }

    // make a vector that describes the probability distribution 
    // pi_i as in Smolla & Akcay over all traits 0 <= i <= par.n_traits. 
    // Initially, we give each trait an equal nonzero weighting
    // of pi_i = 1. Occurrence of traits in the network then increase the
    // count of pi_i > 1. But that pi_i > 1 does not matter as later on we put 
    // this all in a std::discrete_distribution which samples according
    // to relative size. Hence we do not need to scale by n_traits or something
    // like that
    std::vector<int> pi_t(par.n_traits,1);

    // sum of the repertoire size among neighbours
    int sum_repertoire_size = 0;

    // implicit assumption in this pi_i distribution is that when 
    // all values are 0, then trait is chosen randomly
    
    // helper variable to translate
    // network position of the model to location in stack of breeders
    int model_breeder_idx;

    // helper variable to translate
    // network position of the model to sex
    Sex model_sex;

    // loop through all the potential connections of the network
    for (int col_idx = 0; col_idx < max_network_size; ++col_idx)
    {
        // update helper variables that help us to identify 
        // model in stack of breeders and what traits (s)he has
        model_breeder_idx = col_idx >= par.n[female] ? 
            col_idx - par.n[female] : col_idx;

        model_sex = col_idx >= par.n[female] ?
            male : female; 

        assert(model_breeder_idx >= 0);
        assert(model_breeder_idx < par.n[model_sex]);

        // OK, individual is in neighbourhood - network connection is true
        if (metapop[local_patch_idx].network[individual_network_idx][col_idx]) 
        {
            // look at the individual's traits and update pi_i distribution
            // as described on p 8 in Smolla & Akcay
            for (int trait_idx = 0; trait_idx < par.n_traits; ++trait_idx)
            {
                // increase count of this trait
                // we need to look at trait i in this model
                // where model is of sex model_sex
                // where model has the index model_breeder_idx;
                if (metapop[patch_of_origin_idx].breeders[model_sex][
                        model_breeder_idx].repertoire[trait_idx] > 0)
                {
                    // update count of the pi_vector with every 
                    // nonzero trait that we encounter
                    ++pi_t[trait_idx]; 
                    ++sum_repertoire_size;
                }
            } // end for (int trait_idx = 0; trait_idx < par.n_traits, ++trait_idx)

        } // end if (metapop[local_patch_idx].network[individual_network_idx][col_idx]) 

    }// end for (int col_idx = 0; col_idx < max_network_size; ++col_idx)
    
    // make a discrete distribution to choose the trait
    std::discrete_distribution<int> trait_chooser(
            pi_t.begin(),pi_t.end());

    // helper variable selecting the trait that is chosen
    int trait;

    // helper variables to specify probabilities of individual
    // and social learning given that a trait has been selected
    double smolla_akcay_eq_1 = par.gamma / par.n_traits;
    double smolla_akcay_eq_2; 

    // now go through the rounds of learning
    for (int learning_attempt_idx = 0; 
            learning_attempt_idx < par.n_learning_attempts;
            ++learning_attempt_idx)
    {
        // pick a trait according to pi_i
        trait = trait_chooser(rng_r);

        // check whether we perform individual learning
        if (uniform(rng_r) < 
                metapop[local_patch_idx].breeders[
                    offspring_sex][individual_idx].il
                    ) 
        {
            // successful individual learning 
            if (uniform(rng_r) < smolla_akcay_eq_1)
            {

                assert(trait >= 0);
                assert(trait < par.n_traits);


                // success - sample a  trait and increase one's proficiency
                // by one as on p8 in Smolla Akcay
                ++metapop[local_patch_idx].breeders[offspring_sex]
                    [individual_idx].repertoire[trait];
            } // if (uniform(rng_r) < smolla_akcay_eq_1)
        }
        else // ok we perform social learning
        {
            // sigma * pi_t^2
            smolla_akcay_eq_2 = par.sigma * (double)pi_t[trait] / sum_repertoire_size * (double)pi_t[trait] / sum_repertoire_size;

            if (uniform(rng_r) < smolla_akcay_eq_2)
            {
                ++metapop[local_patch_idx].breeders[offspring_sex]
                    [individual_idx].repertoire[trait];
            }
        }
    } // end for int learning attempt

    // and we are done - individual should have achieved a repertoire
    // fingers crossed...

} // end Simulation::learn()

// write statistics
void Simulation::write_data()
{
    // means
    double mean_pp[2] = {0.0,0.0};
    double mean_pc[2] = {0.0,0.0};
    double mean_pr[2] = {0.0,0.0};
    double mean_il = 0.0;
   
    // sums of squares for the variances 
    double ss_pp[2] = {0.0,0.0};
    double ss_pc[2] = {0.0,0.0};
    double ss_pr[2] = {0.0,0.0};
    double ss_il = 0.0;

    // some helper variables to store trait values
    double il,pp,pc,pr;

    for (int patch_idx = 0; patch_idx < metapop.size(); ++patch_idx)
    {
        for (int female_idx = 0; female_idx < par.n[female]; ++female_idx)
        {
            il = metapop[patch_idx].breeders[female][female_idx].il;

            mean_il += il;
            ss_il += il * il;

            for (int sex_idx = 0; sex_idx < 2; ++sex_idx)
            {
                pp = metapop[patch_idx].breeders[female][female_idx].pp[sex_idx];
                
                mean_pp[sex_idx] += pp;
                ss_pp[sex_idx] += pp * pp;
                
                pc = metapop[patch_idx].breeders[female][female_idx].pc[sex_idx];
                
                mean_pc[sex_idx] += pc;
                ss_pc[sex_idx] += pc * pc;
                
                pr = metapop[patch_idx].breeders[female][female_idx].pr[sex_idx];
                
                mean_pr[sex_idx] += pr;
                ss_pr[sex_idx] += pr * pr;
            }
        }
        
        for (int male_idx = 0; male_idx < par.n[male]; ++male_idx)
        {
            il = metapop[patch_idx].breeders[male][male_idx].il;

            mean_il += il;
            ss_il += il * il;

            for (int sex_idx = 0; sex_idx < 2; ++sex_idx)
            {
                pp = metapop[patch_idx].breeders[male][male_idx].pp[sex_idx];
                
                mean_pp[sex_idx] += pp;
                ss_pp[sex_idx] += pp * pp;
                
                pc = metapop[patch_idx].breeders[male][male_idx].pc[sex_idx];
                
                mean_pc[sex_idx] += pc;
                ss_pc[sex_idx] += pc * pc;
                
                pr = metapop[patch_idx].breeders[male][male_idx].pr[sex_idx];
                
                mean_pr[sex_idx] += pr;
                ss_pr[sex_idx] += pr * pr;
            }
        }
    } // for (int patch_idx = 0; patch_idx < metapop.size(); ++patch_idx)

    int n_tot = metapop.size() * (par.n[female] + par.n[male]);
    
    mean_il /= n_tot;
    double var_il = ss_il / n_tot - mean_il * mean_il;

    data_file << time_step << ";" << mean_il << ";" << var_il << ";";

    double var_pp[2] = {0.0,0.0};
    double var_pc[2] = {0.0,0.0};
    double var_pr[2] = {0.0,0.0};

    for (int sex_idx = 0; sex_idx < 2; ++sex_idx)
    {
        mean_pp[sex_idx] /= n_tot;

        var_pp[sex_idx] = ss_pp[sex_idx] / n_tot - 
            mean_pp[sex_idx] * mean_pp[sex_idx];

        data_file << mean_pp[sex_idx] << ";" << var_pp[sex_idx] << ";";

        // and now pr
        
        mean_pr[sex_idx] /= n_tot;

        var_pr[sex_idx] = ss_pr[sex_idx] / n_tot - 
            mean_pr[sex_idx] * mean_pr[sex_idx];

        data_file << mean_pr[sex_idx] << ";" << var_pr[sex_idx] << ";";

        // and now pc
        
        mean_pc[sex_idx] /= n_tot;

        var_pc[sex_idx] = ss_pc[sex_idx] / n_tot - 
            mean_pc[sex_idx] * mean_pc[sex_idx];

        data_file << mean_pc[sex_idx] << ";" << var_pc[sex_idx] << ";";
    } // for (int sex_idx = 0; sex_idx < 2; ++sex_idx)

    data_file << W_global_total << ";";

    double mean_repertoire_size = 0.0;
    double mean_latest_pi_ts = 0.0;

    for (std::vector < double >::iterator it = latest_repertoire_sizes.begin();
            it != latest_repertoire_sizes.end();
            ++it)
    {
        mean_repertoire_size += *it;
    }
    
    mean_repertoire_size /= latest_repertoire_sizes.size();
    
    for (std::vector < double >::iterator it = latest_pi_ts.begin();
            it != latest_pi_ts.end();
            ++it)
    {
        mean_latest_pi_ts += *it;
    }

    mean_latest_pi_ts /= latest_pi_ts.size();

    data_file << mean_repertoire_size << ";" 
        << mean_latest_pi_ts << ";"
        << std::endl;

} // end Simulation::write_data()

// print all nodes and possibly edges 
// in long format (data.frame format)
// To get this into R, see section 3.2 
// from the great tutorial here: 
// https://kateto.net/network-visualization 
void Simulation::write_all_networks()
{
    // helper variable to consider the number of
    // rows and columns of the matrix
    int sum_cols, sum_rows;

    int nf, nm;

    // loop through all patches and plot the respective networks
    for (int patch_idx = 0; patch_idx < metapop.size(); ++patch_idx)
    {
        nf = metapop[patch_idx].nf;
        nm = metapop[patch_idx].nm;
        sum_cols = nf + nm;
        sum_rows = sum_cols;

        // bounds checking
        assert(metapop[patch_idx].network.size() == sum_rows);

        // loop through all rows and columns of learning network
        // and print nonzero connections to file
        for (int row_idx = 0; row_idx < sum_rows; ++row_idx)
        {
            // bounds checking
            assert(metapop[patch_idx].network[row_idx].size() == sum_cols);

            for (int col_idx = 0; col_idx < sum_cols; ++col_idx)
            {
                // only print an entry if there is a connection
                if (metapop[patch_idx].network[row_idx][col_idx])
                {
                    // each row prints current time step and patch index 
                    data_file_network << time_step << ";" << patch_idx << ";";

                    if (row_idx < par.n[female])
                    {
                        // the learner is a female - print her ID
                        data_file_network << "f" << (row_idx + 1) << ";";
                    }
                    else
                    {
                        // the learner is a male  - print his ID
                        // we need to subtract nf from the row index to start
                        // counting males
                        data_file_network << "m" << (row_idx - par.n[female] + 1) << ";";
                    }
                    
                    if (col_idx < par.n[female])
                    {
                        // the model is a female - print her ID
                        data_file_network << "f" << (col_idx + 1) << ";";
                    }
                    else
                    {
                        // the model is a male  - print his ID
                        // we need to subtract nf from the col index to start
                        // counting males
                        data_file_network << "m" << (col_idx - par.n[female] + 1) << ";";
                    }

                    data_file_network << std::endl;
                } // end if connection yes
            } // end for col idx
        } // end for row_idx
    } // end for patch idx
} // end void Simulation::write_all_networks
