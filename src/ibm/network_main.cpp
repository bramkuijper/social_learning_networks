#include "simulation.hpp"
#include "parameters.hpp"

// the main part of the code
int main(int argc, char **argv)
{
    Parameters params;

    std::string temp_file_name_str(argv[1]);

    params.base_name = temp_file_name_str + "_stats";
    params.base_name_matrix_file = temp_file_name_str + "_data";

    Simulation sim(params);

    sim.run();
}
