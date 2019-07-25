#ifndef EV_COLLISIONS_HPP
#define EV_COLLISIONS_HPP

#include "motherbase.hpp"
#include "matrix.hpp"
#include "utility.hpp"
#include <valarray>
#include <cmath>
#include <algorithm>
#include <vector>

#define TEST_COEFF_MULT 5

/*

  USE MEANINGFUL NAMES PLEASE !!!

*/

class CollisionHandler : protected Motherbase
{

private:

  // Global number of collisions
  int n_fake = 0;
  int n_real = 0;
  int n_total = 0;

  ev_matrix::MaskMatrix<real_number> a11;               // Init. = 0.0 (A_i); as many as no. cells
  ev_matrix::MaskMatrix<real_number> vrmax11;           // Max relative speed (C_i); as many as no. cells

  /* NEED TO BE IMPLEMENTED !!! */
  // ev_matrix::MaskMatrix<real_number> freq11;         // As many as no. cells

  ev_matrix::MaskMatrix<real_number> anew;              // Computed at next iteration

  /* THIS HAS NOT BEEN USED YET! */
  // ev_matrix::MaskMatrix<real_number> vrmaxnew;          //    "   "   "   "

  int n_coll;
  ev_matrix::MaskMatrix<int> n_coll_cell;

  std::valarray<real_number> scaled_k;        // Idea: store the unit vector to generate each time
  std::valarray<real_number> rel_vel;
  std::valarray<real_number> delta;           // ???

  inline void gen_scaled_k(void)              // Creates the vector poiting to the cell jc2 (scaled by sigma)
  {
    rng->sample_unit_sphere(scaled_k[0], scaled_k[1], scaled_k[2]);
    scaled_k *= species->get_diam_fluid();
  }

  // Utilities
  real_number vr = 0.0, scalar_prod = 0.0;
  std::vector<int> cells_ind;

  // References from other classes
  const int& npart;
  const int& nx, ny;
  const real_number& xmin, xmax, ymin, ymax;
  const real_number& rdx, rdy;
  const real_number& sigma, delta_t;

  void setup_cell_ind(void);

public:

  CollisionHandler(DSMC*);
  ~CollisionHandler() = default;

  void compute_majorants(void);
  void compute_collision_number(void);
  void perform_collisions(void);

  void perform_collision_kernel(void);

  // Getters
  inline int get_n_fake(void) const { return n_fake; }
  inline int get_n_real(void) const { return n_real; }
  inline int get_n_total(void) const { return n_total; }

  inline int get_n_coll(void) const { return n_coll; }

  inline const ev_matrix::MaskMatrix<real_number>& get_a11(void) const { return a11; }
  inline const ev_matrix::MaskMatrix<real_number>& get_vrmax11(void) const { return vrmax11; }
  inline const ev_matrix::MaskMatrix<int>& get_n_coll_cell(void) const { return n_coll_cell; }

};

#endif /* EV_COLLISIONS_HPP */
