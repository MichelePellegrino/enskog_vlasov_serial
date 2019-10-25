/*! \file sampling.hpp
 *  \brief Header containing class for collecting and averaging samples
 */

#ifndef EV_SAMPLING_HPP
#define EV_SAMPLING_HPP

#include "motherbase.hpp"
#include "matrix.hpp"

/*! \class Sampler
 *  \brief Class for collecting and averaging samples
 *
 *  In its current implementation the class samples:
 *    - number of particles per cell
 *    - number density
 *    - averaged reduced density
 *    - velocity field
 *    - temperature
 *    - streaming pressure tensor
 *    - heat flux
 */
class Sampler : protected Motherbase
{

private:

  int outer_counter = 0;

  ev_matrix::MaskMatrix<int> inner_counter;
  ev_matrix::MaskMatrix<real_number> inner_counter_cast;
  ev_matrix::MaskMatrix<real_number> dt_factor;

  ev_matrix::MaskMatrix<real_number> vx_avg;
  ev_matrix::MaskMatrix<real_number> vy_avg;
  ev_matrix::MaskMatrix<real_number> vz_avg;
  ev_matrix::MaskMatrix<real_number> temp_avg;
  ev_matrix::MaskMatrix<real_number> pxx_avg;
  ev_matrix::MaskMatrix<real_number> pyy_avg;
  ev_matrix::MaskMatrix<real_number> pzz_avg;
  ev_matrix::MaskMatrix<real_number> pxy_avg;
  ev_matrix::MaskMatrix<real_number> pxz_avg;
  ev_matrix::MaskMatrix<real_number> pyz_avg;
  ev_matrix::MaskMatrix<real_number> qx_avg;
  ev_matrix::MaskMatrix<real_number> qy_avg;
  ev_matrix::MaskMatrix<real_number> qz_avg;

  ev_matrix::MaskMatrix<real_number> numdens_avg;
  ev_matrix::MaskMatrix<real_number> aveta_avg;
  ev_matrix::MaskMatrix<real_number> fx_avg;
  ev_matrix::MaskMatrix<real_number> fy_avg;

public:

  Sampler(DSMC*);
  ~Sampler() = default;

  void reset(void);
  void sample(void);
  void average(void);

  inline const ev_matrix::MaskMatrix<real_number>& get_vx_avg(void) const { return vx_avg; }
  inline ev_matrix::MaskMatrix<real_number>& get_vx_avg(void) { return vx_avg; }
  inline const ev_matrix::MaskMatrix<real_number>& get_vy_avg(void) const { return vy_avg; }
  inline ev_matrix::MaskMatrix<real_number>& get_vy_avg(void) { return vy_avg; }
  inline const ev_matrix::MaskMatrix<real_number>& get_vz_avg(void) const { return vz_avg; }
  inline ev_matrix::MaskMatrix<real_number>& get_vz_avg(void) { return vz_avg; }
  inline const ev_matrix::MaskMatrix<real_number>& get_temp_avg(void) const { return temp_avg; }
  inline ev_matrix::MaskMatrix<real_number>& get_temp_avg(void) { return temp_avg; }
  inline const ev_matrix::MaskMatrix<real_number>& get_pxx_avg(void) const { return pxx_avg; }
  inline ev_matrix::MaskMatrix<real_number>& get_pxx_avg(void) { return pxx_avg; }
  inline const ev_matrix::MaskMatrix<real_number>& get_pyy_avg(void) const { return pyy_avg; }
  inline ev_matrix::MaskMatrix<real_number>& get_pyy_avg(void) { return pyy_avg; }
  inline const ev_matrix::MaskMatrix<real_number>& get_pzz_avg(void) const { return pzz_avg; }
  inline ev_matrix::MaskMatrix<real_number>& get_pzz_avg(void) { return pzz_avg; }
  inline const ev_matrix::MaskMatrix<real_number>& get_pxy_avg(void) const { return pxy_avg; }
  inline ev_matrix::MaskMatrix<real_number>& get_pxy_avg(void) { return pxy_avg; }
  inline const ev_matrix::MaskMatrix<real_number>& get_pxz_avg(void) const { return pxz_avg; }
  inline ev_matrix::MaskMatrix<real_number>& get_pxz_avg(void) { return pxz_avg; }
  inline const ev_matrix::MaskMatrix<real_number>& get_pyz_avg(void) const { return pyz_avg; }
  inline ev_matrix::MaskMatrix<real_number>& get_pyz_avg(void) { return pyz_avg; }
  inline const ev_matrix::MaskMatrix<real_number>& get_qx_avg(void) const { return qx_avg; }
  inline ev_matrix::MaskMatrix<real_number>& get_qx_avg(void) { return qx_avg; }
  inline const ev_matrix::MaskMatrix<real_number>& get_qy_avg(void) const { return qy_avg; }
  inline ev_matrix::MaskMatrix<real_number>& get_qy_avg(void) { return qy_avg; }
  inline const ev_matrix::MaskMatrix<real_number>& get_qz_avg(void) const { return qz_avg; }
  inline ev_matrix::MaskMatrix<real_number>& get_qz_avg(void) { return qz_avg; }
  inline const ev_matrix::MaskMatrix<real_number>& get_numdens_avg(void) const { return numdens_avg; }
  inline ev_matrix::MaskMatrix<real_number>& get_numdens_avg(void) { return numdens_avg; }
  inline const ev_matrix::MaskMatrix<real_number>& get_forces_x_avg(void) const { return fx_avg; }
  inline ev_matrix::MaskMatrix<real_number>& get_forces_x_avg(void) { return fx_avg; }
  inline const ev_matrix::MaskMatrix<real_number>& get_forces_y_avg(void) const { return fy_avg; }
  inline ev_matrix::MaskMatrix<real_number>& get_forces_y_avg(void) { return fy_avg; }
  inline const ev_matrix::MaskMatrix<real_number>& get_aveta_avg(void) const { return aveta_avg; }
  inline ev_matrix::MaskMatrix<real_number>& get_aveta_avg(void) { return aveta_avg; }

};

#endif /* EV_SAMPLING_HPP */
