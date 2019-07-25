#include "dsmc.hpp"

int main()
{

try
{
  DSMC my_dsmc_instanziation("input_files/drop_par.dat.txt");
}
catch(const char* ex)
{
  std::cout << ex << std::endl;
}

return 0;

}
