#include "dsmc.hpp"

int main()
{

try
{
  DSMC my_dsmc_instanziation("input_files/uniform02.dat.txt");
}
catch(const char* ex)
{
  std::cout << ex << std::endl;
}

return 0;

}
