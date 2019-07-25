#include "sampling.hpp"
#include "grid.hpp"
#include "particles.hpp"
#include "density.hpp"

Sampler::Sampler(DSMC* dsmc):

  Motherbase(dsmc),

  inner_counter( 0, grid->get_n_cells_x(), 0, grid->get_n_cells_y(), 0 ),
  inner_counter_cast( 0, grid->get_n_cells_x(), 0, grid->get_n_cells_y(), 0.0 ),
  
  vx_avg( 0, grid->get_n_cells_x(), 0, grid->get_n_cells_y(), 0.0 ),
  vy_avg( 0, grid->get_n_cells_x(), 0, grid->get_n_cells_y(), 0.0 ),
  vz_avg( 0, grid->get_n_cells_x(), 0, grid->get_n_cells_y(), 0.0 ),
  temp_avg( 0, grid->get_n_cells_x(), 0, grid->get_n_cells_y(), 0.0 ),
  pxx_avg( 0, grid->get_n_cells_x(), 0, grid->get_n_cells_y(), 0.0 ),
  pyy_avg( 0, grid->get_n_cells_x(), 0, grid->get_n_cells_y(), 0.0 ),
  pzz_avg( 0, grid->get_n_cells_x(), 0, grid->get_n_cells_y(), 0.0 ),
  pxy_avg( 0, grid->get_n_cells_x(), 0, grid->get_n_cells_y(), 0.0 ),
  pxz_avg( 0, grid->get_n_cells_x(), 0, grid->get_n_cells_y(), 0.0 ),
  pyz_avg( 0, grid->get_n_cells_x(), 0, grid->get_n_cells_y(), 0.0 ),
  qx_avg( 0, grid->get_n_cells_x(), 0, grid->get_n_cells_y(), 0.0 ),
  qy_avg( 0, grid->get_n_cells_x(), 0, grid->get_n_cells_y(), 0.0 ),
  qz_avg( 0, grid->get_n_cells_x(), 0, grid->get_n_cells_y(), 0.0 )

  { }

void
Sampler::reset
(void)
{
  inner_counter = 0;
  inner_counter_cast = 0.0;
  vx_avg = 0.0;
  vy_avg = 0.0;
  vz_avg = 0.0;
  temp_avg = 0.0;
  pxx_avg = 0.0;
  pyy_avg = 0.0;
  pzz_avg = 0.0;
  pxy_avg = 0.0;
  pxz_avg = 0.0;
  pyz_avg = 0.0;
  qx_avg = 0.0;
  qy_avg = 0.0;
  qz_avg = 0.0;
}

void
Sampler::sample
(void)
{
  int nc = grid->get_n_cells();
  int i, j, idx_p;
  real_number vx, vy, vz, e_kin;
  outer_counter++;
  for(int idx_c = 0; idx_c<nc; ++idx_c)
  {
    i = grid->lexico_inv(idx_c).first;
    j = grid->lexico_inv(idx_c).second;
    for(int k = density->iof(idx_c); k<density->iof(idx_c+1); ++k)
    {
      idx_p = density->ind(k);
      vx = ensemble->get_vx(idx_p);
      vy = ensemble->get_vy(idx_p);
      vz = ensemble->get_vz(idx_p);
      inner_counter(i,j)++;
      vx_avg(i,j) += vx;
      vy_avg(i,j) += vy;
      vz_avg(i,j) += vz;
      pxx_avg(i,j) += vx*vx;
      pyy_avg(i,j) += vy*vy;
      pzz_avg(i,j) += vz*vz;
      pxy_avg(i,j) += vx*vy;
      pxz_avg(i,j) += vx*vz;
      pyz_avg(i,j) += vy*vz;
      e_kin = vx*vx + vy*vy + vz*vz;
      temp_avg(i,j) += e_kin;
      qz_avg(i,j) += vz*e_kin;
      qx_avg(i,j) += vx*e_kin;
      qy_avg(i,j) += vy*e_kin;
    }
  }
}

void
Sampler::average
(void)
{
  inner_counter_cast.copy_cast<int>(inner_counter);
  vx_avg /= inner_counter_cast;
  vy_avg /= inner_counter_cast;
  vz_avg /= inner_counter_cast;
  pxx_avg /= inner_counter_cast;  pxx_avg -= vx_avg*vx_avg;
  pyy_avg /= inner_counter_cast;  pyy_avg -= vy_avg*vy_avg;
  pzz_avg /= inner_counter_cast;  pzz_avg -= vz_avg*vz_avg;
  pxy_avg /= inner_counter_cast;  pxy_avg -= vx_avg*vy_avg;
  pxz_avg /= inner_counter_cast;  pxz_avg -= vx_avg*vz_avg;
  pyz_avg /= inner_counter_cast;  pyz_avg -= vy_avg*vz_avg;
  qz_avg /= 2.0*inner_counter_cast;
  qx_avg /= 2.0*inner_counter_cast;
  qy_avg /= 2.0*inner_counter_cast;
  temp_avg /= inner_counter_cast;
  temp_avg = ( temp_avg -
    vx_avg*vx_avg - vy_avg*vy_avg - vz_avg*vz_avg ) / 3.0;
}
