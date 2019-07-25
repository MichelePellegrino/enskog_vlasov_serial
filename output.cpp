#include "output.hpp"
#include "collisions.hpp"
#include "density.hpp"
#include "force_field.hpp"

Output::Output(DSMC* dsmc):
  Motherbase(dsmc)
  { }

void
Output::output_kernel
(void)
{

  std::ofstream file1("output_files/kernel.txt");

  // Output potential kernel
  file1 << mean_field->get_kernel_matrix();
  file1.close();

}

void
Output::output_forces
(void)
{

  std::ofstream file1("output_files/forces_x.txt");
  std::ofstream file2("output_files/forces_y.txt");

  // Output Fx
  file1 << mean_field->get_force_x();
  file1.close();

  // Output Fy
  file2 << mean_field->get_force_y();
  file2.close();

}

void
Output::output_density
(void)
{

  std::ofstream file1("output_files/npc.txt");
  std::ofstream file2("output_files/aveta.txt");

  // Output no. particles per cell
  file1 << density->get_npc();
  file1.close();

  // Output averaged reduced density
  file2 << density->get_aveta();
  file2.close();

}

void
Output::output_majorants
(void)
{

  std::ofstream file1("output_files/A.txt");
  std::ofstream file2("output_files/C.txt");

  // Output A majorant
  file1 << collision_handler->get_a11();
  file1.close();

  // Output C majorant
  file2 << collision_handler->get_vrmax11();
  file2.close();

}

void
Output::output_collisions
(void)
{

  std::ofstream file1("output_files/n_coll.txt");

  // Output potential kernel
  file1 << collision_handler->get_n_coll_cell();
  file1.close();

}

// DEBUG

void
Output::output_fun_vec
(const std::vector<real_number>& x, const std::vector<real_number>& y)
{

  std::ofstream file1("output_files/fun_vec.txt");
  int n = x.size();
  for (int i = 0; i<n; ++i)
    file1 << x[i] << "\t" << y[i] << "\n";
  file1.close();

}
