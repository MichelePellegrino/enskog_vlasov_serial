#include "force_field.hpp"
#include "configuration.hpp"
#include "density.hpp"
#include "grid.hpp"

// DEBUG
// # # # # #
#include "output.hpp"
// # # # # #

ForceField::ForceField(DSMC* dsmc):
  Motherbase(dsmc),
  n_cutoff_x( density->get_n_cutoff_x() ),
  n_cutoff_y( density->get_n_cutoff_y() ),
  diamol( species->get_diam_fluid() ),
  dx( grid->get_dx() ),
  dy( grid->get_dy() ),
  kernel_function( potential->get_pot_kernel() ),
  finite_integrator( psi, DUMMY_A, DUMMY_B ),
  infinite_integrator( psi, DUMMY_A, DUMMY_B ),
  kernel_matrix(n_cutoff_x, n_cutoff_y, 0.0),
  force_x_matrix(0, grid->get_n_cells_x(), 0, grid->get_n_cells_y(), 0.0),
  force_y_matrix(0, grid->get_n_cells_x(), 0, grid->get_n_cells_y(), 0.0),
  force_x_convolutioner(force_x_matrix, kernel_matrix,
    density->get_num_dens_cell(), 0.0),
  force_y_convolutioner(force_y_matrix, kernel_matrix,
    density->get_num_dens_cell(), 0.0)
{
  std::cout << "### COMPUTING POTENTIAL KERNEL MATRIX ###" << std::endl;
  compute_kernel_matrix();
  // TESTING MEAN FIELD
  // # # # # #
  /*
  std::cout << "### TEST: computing force field ###" << std::endl;
  compute_force_field();
  */
  // # # # # #
}

void
ForceField::compute_kernel_matrix
(void)
{
  real_number pot_int(0.0), distx(0.0), disty(0.0);
  for (int i = -n_cutoff_x; i<=n_cutoff_x; ++i)
  {
    for (int j = -n_cutoff_y; j<=n_cutoff_y; ++j)
    {
      distx = i*dx;
      disty = j*dy;
      dist2 = distx*distx + disty*disty;
      pot_int = compute_integral();
      kernel_matrix(i,j) = pot_int;
    }
  }
}

real_number
ForceField::compute_integral
(void)
{
  real_number tmp = dist2 - diamol*diamol;
  if ( tmp > 0.0 )
  {
    return (
      infinite_integrator.integrate(psi, ev_const::minfty, -sqrt(tmp))
      + finite_integrator.integrate(psi, -sqrt(tmp), sqrt(tmp))
      + infinite_integrator.integrate(psi, sqrt(tmp), ev_const::pinfty)   );
  }
  else if ( tmp < 0.0 )
  {
    return (
      infinite_integrator.integrate(psi, ev_const::minfty, -sqrt(-tmp))
      + infinite_integrator.integrate(psi, sqrt(-tmp), ev_const::pinfty)  );
  }
  else
  {
    return (
      infinite_integrator.integrate(psi, ev_const::minfty, -zero_threshold)
      + infinite_integrator.integrate(psi, zero_threshold, ev_const::pinfty) );
  }
}

void
ForceField::compute_force_field
(void)
{
  force_x_convolutioner.convolute();
  force_x_matrix /= (dx*dy);
  force_y_convolutioner.convolute();
  force_y_matrix /= (dx*dy);
}


// TESTING
void
ForceField::testing_output_kernel_function
(int n_values)
{
  real_number pot_int(0.0), dist(0.0);
  real_number step = (n_cutoff_x*dx)/(n_values);
  std::vector<real_number> distance_values;
  std::vector<real_number> kernel_values;
  for (int i = 0; i<=n_values; ++i)
  {
    dist = i*step;
    dist2 = dist*dist;
    pot_int = compute_integral();
    distance_values.push_back(dist);
    kernel_values.push_back(pot_int);
  }
  output->output_fun_vec(distance_values, kernel_values);
}
