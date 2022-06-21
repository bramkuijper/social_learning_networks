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
    ,patch_sampler{0,(int)par.npatches-1} // initialize the patch sampling distribution
    ,data_file{par.base_name.c_str()} // initialize the data file by giving it a name
    ,data_file_network{par.base_name_matrix_file.c_str()} // initialize network data file
                                                          // by giving it a name
    ,par{par} // initialize the parameter data member with the constructor argument
    ,metapop{par.npatches, Patch(par)} // initialize a meta population each with n1 individuals of species 1 and n2 individuals of species 2
{
    // initialize patch environments
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

} // end IBM_Mutualism

// run the simulation
void Simulation::run()
{
    // write the headers to the output file
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

        switch(event_chooser(rng_r))
        {
            case 0:
                environmental_change(0);
                break;
            case 1:
                environmental_change(1);
                break;
            case 2:
                death_birth();
                break;
            default:
                throw std::range_error("Something is going wrong in the event chooser. Not good.");
                break;
        }
    }

    // write parameters to output
    void write_parameters();
} // end Simulation::run()

// change the environment of a local pod
// argument specifies current environmental state of the patch
// that needs to be changed
void Simulation::environmental_change(bool const is_envt2)
{
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

            break;
        }
    } // end for
}


// write parameters to a file
void Simulation::write_parameters()
{
    data_file << std::endl
        << std::endl;

    for (int sex_idx = 0; sex_idx < 2; ++sex_idx)
    {

    }
}

void Simulation::death_birth()
{
    // pick random individual to die
    // effectively there are npatches * (nf + nm) individuals
    int random_patch = patch_sampler(rng_r);

    // whether it is a male or a female that dies
    bool is_male = uniform(rng_r) < 0.5;

    // make a distribution to randomly sample the individual
    // of a particular sex who will die 
    std::uniform_int_distribution<int> 
        to_die_sampler(0, metapop[random_patch].breeders[is_male].size());

    // sample the individual
    int individual_id = to_die_sampler(rng_r);

    int individual_location_in_network = individual_id;

    // if this individual i is a male, note that its index
    // in the network is nf + i
    if (is_male)
    {
        individual_location_in_network += metapop[random_patch].nf;
    }

    // get number of columns in aux variable to speed things up
    int ncols = metapop[random_patch].network[
        individual_location_in_network].size();

    // bounds checking
    assert(individual_location_in_network >= 0);
    assert(individual_location_in_network < ncols);
    assert(individual_location_in_network < metapop[random_patch].nf 
            + metapop[random_patch].nm);

    // set everything in the social learning network 
    // relating to this individual to 0
    for (int col_idx = 0; col_idx < ncols; ++col_idx)
    {
        // set all connections in which 
        // individual serves as learner
        // to 0
        metapop[random_patch].network[
            individual_location_in_network][col_idx] = false;

        // set all connections in which 
        // individual serves as model
        // to 0
        metapop[random_patch].network[col_idx][
            individual_location_in_network] = false;
    }
} // end death_birth



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
