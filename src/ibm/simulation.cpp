// the heart of the code:  Simulation::functions describing the 
// different parts of the life cycle

#include <random> // random number generating libraries
#include <vector>
#include <cassert>
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
    ,patch_sampler{0,(int)par.npatches - 1} // initialize the patch sampling distribution
    ,data_file{par.base_name.c_str()} // initialize the data file by giving it a name
    ,data_file_network{par.base_name_matrix_file.c_str()} // initialize network data file
                                                          // by giving it a name
    ,par{par} // initialize the parameter data member with the constructor argument
    ,metapop{par.npatches, Patch(par)} // initialize a meta population each with n1 individuals of species 1 and n2 individuals of species 2
{
    // initialize patch environments
    //
    // count the number of patches in environmental state 2
    // (number of patches in environmental state 1 = par.npatches - n_patches_2)
    n_patches_2 = 0;

    // overall frequency of environment 2 is 
    // switch_rate[1->2]/(switch_rate[1->2] + switch_rate[2->1])
    double prob_envt2 = par.switch_rate[1]/(par.switch_rate[1] + 
            par.switch_rate[0]);

    for (int patch_idx = 0; patch_idx < metapop.size(); ++patch_idx)
    {
        metapop[patch_idx].envt2 = uniform(rng_r) < prob_envt2;

        n_patches_2 += metapop[patch_idx].envt2;
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

        double death_rate = metapop.size() * (par.n[0] + par.n[1]);

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
    } // end for time_step

    // write the networks to a file
    void write_all_networks();

    // write parameters to output
    void write_parameters();
} // end Simulation::run()

// change the environment of a local pod
// argument specifies current environmental state of the patch
// that needs to be changed
void Simulation::environmental_change(bool const is_envt2)
{
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
            metapop[random_patch].envt2 != metapop[random_patch].envt2;

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
            metapop[random_patch].calculate_W();

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
            << (envt_idx + 1) << ";" << 
            << par.switch_rate[envt_idx] << std::endl;
    }

    data_file << "random_seed;" << seed << std::endl;
    data_file << "n_traits;" << par.n_traits << std::endl;
    data_file << "npatches;" << par.npatches << std::endl;
    data_file << "n_learning_attempts;" << par.n_learning_attempts << std::endl;

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
    // effectively there are npatches * (nf + nm) individuals
    int random_patch_idx = patch_sampler(rng_r);

    // whether it is a male or a female that dies
    // if 0 (false) female, if 1 (true) male - see also 
    // the Sex enum in parameters.hpp
    Sex is_male = uniform(rng_r) < par.nm / (par.nf + par.nm);

    // make a distribution to randomly sample the individual
    // of a particular sex who will die 
    std::uniform_int_distribution<int> 
        to_die_sampler(0, metapop[random_patch_idx].breeders[is_male].size());

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
        individual_location_in_network += 
            metapop[random_patch_idx].nf;
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
} // end Simulation::death_birth()

void Simulation::make_new_individual(
        int const local_patch_idx // index of patch where offspring will establish as breeder
        ,bool const offspring_is_male // offspring will be male or female 
        ,int const individual_idx // index in vector of breeders that this individual will take up
        )
{
    assert(local_patch_idx >= 0);
    assert(local_patch_idx < metapop.size());

    // fitness of the local patch
    double Wloc = metapop[local_patch_idx].Wtot();

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

    // sample from local patch unless..
    int patch_to_sample = local_patch_idx;

    // we sample a remotely born offspring
    if (uniform(rng) < d[is_male] * Wremote / 
            (d[is_male] * Wremote + (1.0 - d[is_male]) * Wloc))
    {
        patch_to_sample = patch_sampler(rng);
    }

    // now sample parents according to their fitnesses within that patch
    // and make offspring. This may include the 'dead' parent, implying that the
    // parent only dies after reproduction (cf. Smolla & Akcay)
    std::vector<double> Wifs;

    for (int female_idx = 0; female_idx < par.n[female]; ++female_idx)
    {
        Wifs.push_back(metapop[patch_to_sample].breeders[female][female_idx].Wi);
    }
    
    std::vector<double> Wims;

    for (int male_idx = 0; male_idx < par.n[male]; ++male_idx)
    {
        Wims.push_back(metapop[patch_to_sample].breeders[male][male_idx].Wi);
    }


    std::discrete_distribution<int> female_sampler(Wifs.begin(), Wifs.end());
    std::discrete_distribution<int> male_sampler(Wims.begin(), Wims.end());


    // bounds checking
    assert(mother_idx >= 0);
    assert(mother_idx < par.n[female]);
    assert(mother_idx < metapop[patch_to_sample].breeders[female].size());

    assert(father_idx >= 0);
    assert(father_idx < par.n[male]);
    assert(father_idx < metapop[patch_to_sample].breeders[male].size());

    // use the birth constructor  in individual.cpp
    // to produce a new offspring
    Individual offspring(
        metapop[patch_to_sample].breeders[female][mother_idx] // mom
        ,metapop[patch_to_sample].breeders[male][father_idx] // dad
        ,par // Parameter object to get at mutation rates
        ,rng);


    // put offspring in place of the now dead breeder
    // and then perform learning
    metapop[local_patch_idx].breeders[offspring_is_male][individual_idx] = offspring;

    // now perform learning, which 
    // could happen either post or pre dispersal
}// end Simulation::sample_parents()

// writes the headers to the data files
void Simulation::write_data_headers()
{
    // headers of the network file
    data_file_network << "time_step;patch;from;to;" << std::endl;

    // headers of the normal data file
    data_file << "time_step;" << std::endl;;
}

void Simulation::write_data()
{
}

// print all nodes and possibly edges 
// in long format (data.frame format)
// To get this into R, see section 3.2 
// from the great tutorial here: 
// https://kateto.net/network-visualization 
void Simulation::write_all_networks()
{
    // auxiliary variable to consider the number of
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

                    if (row_idx < nf)
                    {
                        // the learner is a female - print her ID
                        data_file_network << "f" << (row_idx + 1) << ";";
                    }
                    else
                    {
                        // the learner is a male  - print his ID
                        // we need to subtract nf from the row index to start
                        // counting males
                        data_file_network << "m" << (row_idx - nf + 1) << ";";
                    }
                    
                    if (col_idx < nf)
                    {
                        // the model is a female - print her ID
                        data_file_network << "f" << (col_idx + 1) << ";";
                    }
                    else
                    {
                        // the model is a male  - print his ID
                        // we need to subtract nf from the col index to start
                        // counting males
                        data_file_network << "m" << (col_idx - nf + 1) << ";";
                    }

                    data_file_network << std::endl;
                } // end if connection yes
            } // end for col idx
        } // end for row_idx
    } // end for patch idx
} // end void Simulation::write_all_networks
