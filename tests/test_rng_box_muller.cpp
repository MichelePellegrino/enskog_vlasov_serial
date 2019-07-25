#include <iostream>
#include <vector>

#include "random.hpp"
#include "types.hpp"

int main()
{

std::cout << "Insert the number of particles:" << std::endl;

int N;
std::cin >> N;

std::vector<double> vx(N, 0.0);
std::vector<double> vy(N, 0.0);
std::vector<double> vz(N, 0.0);

ev_random::CustomRngObject<RNG> rng;

for (int i = 0; i<N; ++i)
{

rng.sample_box_muller(
  1.0, 0.0, 0.0, 1.0,
  vx[i], vy[i], vz[i]);

std::cout << i << " particles" << std::endl;

}

return 0;

}
