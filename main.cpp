/*! \file main.hpp
 *  \brief Source code for the executable
 */

// Includes wrapper for all DSMC modules
#include "dsmc.hpp"

// Includes headers for file input
#include <fstream>
#include <string>

int main()
{

/*!
 * The name for the configuration file is read from a input file; the name (string)
 * is stored in input_files/read_input_file.txt
 */
std::ifstream file_stream("input_files/read_input_file.txt");
std::string input_file_name;
file_stream >> input_file_name;
file_stream.close();

/*!
 * DSMC class requires the name for the input configuration file; the DSMC loop is
 * called in the constructor (i.e. the following executes the DSMC procedure defined
 * in the conf. file)
 */
try
{
  DSMC my_dsmc_instanziation(input_file_name);
}
catch(const char* ex)
{
  std::cout << ex << std::endl;
}

return 0;

}
