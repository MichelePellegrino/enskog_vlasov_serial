#ifndef EV_SAMPLING_HPP
#define EV_SAMPLING_HPP

#include "motherbase.hpp"
#include "matrix.hpp"

class Sampler : protected Motherbase
{

private:

  int outer_counter;
  ev_matrix::MaskMatrix<int> inner_counter;
  ev_matrix::MaskMatrix<real_number> inner_counter_cast;
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

public:

  Sampler(DSMC*);
  ~Sampler() = default;
  void reset(void);
  void sample(void);
  void average(void);

  const ev_matrix::MaskMatrix<real_number>& get_temp_avg(void) const { return temp_avg; }
  ev_matrix::MaskMatrix<real_number>& get_temp_avg(void) { return temp_avg; }

};

#endif /* EV_SAMPLING_HPP */
