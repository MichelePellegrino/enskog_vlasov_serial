#include "density.hpp"
#include "grid.hpp"
#include "species.hpp"
#include "particles.hpp"
#include "configuration.hpp"

// DEBUG
// # # # # #
#include <algorithm>
#include <set>
// # # # # #

DensityKernel::DensityKernel
(DSMC* dsmc):
  Motherbase(dsmc),
  reduce_factor( ( ev_const::pi/6.0 ) * ev_utility::power<3>(species->get_diam_fluid()) ),
  ns_x( (int)( species->get_hdiam_fluid() / ( grid->get_dx() * ev_const::sqrt2 ) ) ),
  ns_y( (int)( species->get_hdiam_fluid() / ( grid->get_dy() * ev_const::sqrt2 ) ) ),
  stencil_x(ns_x, ns_y, 0.0),
  stencil_y(ns_x, ns_y, 0.0),
  weights(ns_x, ns_y, 0.0),
  x_min(grid->get_x_min()),
  y_min(grid->get_y_min()),
  n_cutoff_x( (int)(conf->get_x_extra()/grid->get_dx())+1 ),
  n_cutoff_y( (int)(conf->get_x_extra()/grid->get_dx())+1 ),
  n_part_cell( 0, grid->get_n_cells_x(), 0, grid->get_n_cells_y(), 0 ),
  num_dens_cell( ev_matrix::MaskMatrix<real_number> (
    -n_cutoff_x, grid->get_n_cells_x()+n_cutoff_x,
    -n_cutoff_y, grid->get_n_cells_y()+n_cutoff_y ),
    n_cutoff_x, n_cutoff_y ),
  reduced_density( ev_matrix::MaskMatrix<real_number> (
    -n_cutoff_x, grid->get_n_cells_x()+n_cutoff_x,
    -n_cutoff_y, grid->get_n_cells_y()+n_cutoff_y ),
    n_cutoff_x, n_cutoff_y ),
  average_reduced_density( 0, grid->get_n_cells_x(), 0, grid->get_n_cells_y(), 0.0 ),
  avg_convolutioner( average_reduced_density, weights, reduced_density, 0.0 ),
  idx_cell( ensemble->get_n_particles(), 0 ),
  idx_map( ensemble->get_n_particles(), 0 ),
  cum_num( grid->get_n_cells()+1, 0 ),
  raw_num( grid->get_n_cells(), 0 )
  {

    // Initialize weights
    real_number dx = grid->get_dx();
    real_number dy = grid->get_dy();
    real_number sigma = species->get_diam_fluid();
    real_number hsigma = sigma/2.0;
    real_number sum_w = 0.0;
    real_number ns_x = (int)( hsigma / ( dx * ev_const::sqrt2 ) ) ;
    real_number ns_y = (int)( hsigma / ( dy * ev_const::sqrt2 ) ) ;
    std::cout << "### COMPUTING AVERAGING WEIGHTS ###" << std::endl;

    for (int i =-ns_x; i<=ns_x; ++i)
    {
      for (int j =-ns_y; j<=ns_y; ++j)
      {
        stencil_x(i, j) = i * dx;
        stencil_y(i, j) = j * dy;
        weights(i,j) = 12.0 / ( ev_const::pi * ev_utility::power<3>(sigma) )
          * sqrt( hsigma * hsigma - stencil_x(i,j) * stencil_x(i,j)
          - stencil_y(i,j) * stencil_y(i,j) ) * dx * dy;
        sum_w = sum_w + weights(i,j);
      }
    }

    weights /= sum_w;
    std::cout << "### COMPUTING PARTICLE MAP ###" << std::endl;
    compute_ind_map_part();

  }

void
DensityKernel::binning
(void)
{
  n_part_cell = 0;
  int NP = ensemble->get_n_particles();
  idx_cell.assign(NP, 0);
  int i,j;
  for (int k = 0; k<NP; ++k)
  {
    i = ensemble->get_cx(k);
    j = ensemble->get_cy(k);
    n_part_cell(i,j) += 1;
    idx_cell[k] = grid->lexico(i,j);
  }
  // Setting cumulate density
  int NC = grid->get_n_cells();
  cum_num.assign(NC+1, 0);
  for (int k = 1; k<NC+1; ++k)
    cum_num[k] = cum_num[k-1] + n_part_cell(grid->lexico_inv(k-1).first, grid->lexico_inv(k-1).second);
}

void
DensityKernel::compute_ind_map_part
(void)
{
  int NP = ensemble->get_n_particles();
  int jc, k;
  idx_map.assign(NP, 0);
  raw_num.assign(grid->get_n_cells(), 0);
  for (int jp = 0; jp<NP; ++jp)
  {
    jc = idx_cell[jp];
    k = cum_num[jc] + raw_num[jc];
    raw_num[jc] += 1;
    idx_map[k] = jp;
  }
}

void
DensityKernel::fill_dummy_field
(void)
{
  // Copying inner data
  num_dens_cell.copy_patch<int>(n_part_cell, 0, 0);
  // CASE: periodic B.C. on each and every edge...
  for(int r = 0; r < N_BUF; ++r )
    num_dens_cell.set_outer_block( r, num_dens_cell.get_inner_block( ev_matrix::reflect_idx(r) ) );
  num_dens_cell /= grid->get_cell_volume();
}

void
DensityKernel::compute_reduced_density
(void)
{
  reduced_density = num_dens_cell;
  reduced_density *= reduce_factor;
}

void
DensityKernel::compute_avg_density
(void)
{
  avg_convolutioner.convolute();
}

void
DensityKernel::perform_density_kernel
(void)
{
  binning();
  fill_dummy_field();
  compute_reduced_density();
  compute_avg_density();
  // DEBUG
  // # # # # #
  // print_binned_particles();
  // print_reduced_numdens();
  // print_reduced_aveta();
  // # # # # #
}

// TESTING
void
DensityKernel::print_binned_particles
(void) const
{
  std::cout << " >> binned particles = " << n_part_cell.sum() << std::endl;
}

void
DensityKernel::print_reduced_numdens
(void) const
{
  std::cout << " >> max numdens = " << num_dens_cell.maxCoeff() << std::endl;
}

void
DensityKernel::print_reduced_aveta
(void) const
{
  std::cout << " >> max aveta = " << average_reduced_density.maxCoeff() << std::endl;
}

void
DensityKernel::test_particle_cell_map(void)
const
{
  std::cout << "### TESTING PARTICLE-CELL MAP ###" << std::endl;
  std::cout << " >> iof(end) = " << cum_num[grid->get_n_cells()-1] << std::endl;
  int jp;
  std::vector<int> pt;
  for ( int n = 0; n<grid->get_n_cells(); ++n )
  {
    for ( int k=cum_num[n]; k<cum_num[n+1]; ++k )
    {
      jp = idx_map[k];
      pt.push_back( ensemble->get_p_tag(jp) );
    }
  }
  std::sort(pt.begin(), pt.end());
  std::cout << " >> pt(beg) = " << pt[0] << std::endl;
  for ( int p = 1; p<ensemble->get_n_particles(); ++p)
  {
    // assert( pt[p]==(pt[p-1]+1) && "Something's wrong with the particle-cell map" );
    if ( pt[p]!=(pt[p-1]+1) )
      std::cout << "(" << pt[p-1] << "->" << pt[p] << ")" << std::endl;
  }
}
