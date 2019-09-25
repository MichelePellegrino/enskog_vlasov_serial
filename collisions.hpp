/*! \file collisions.hpp
 *  \brief Header containing class for collisions simulation
 */

#ifndef EV_COLLISIONS_HPP
#define EV_COLLISIONS_HPP

#include "motherbase.hpp"
#include "matrix.hpp"
#include "utility.hpp"

/* Valarray turns out to be handy when dealing to fixed-dimension array behaving
    like mathematical vectors */
#include <valarray>
#include <cmath>
#include <algorithm>
#include <vector>

/*! \def TEST_COEFF_MULT
    \brief Default value for the initial number of evaluation to estimate A_i and C_i
*/
#ifndef TEST_COEFF_MULT
#define TEST_COEFF_MULT 5
#endif

/*! \def DEFAULT_ALPHA_1
    \brief Parameter tuning out-of-bound/real ratio
*/
#ifndef DEFAULT_ALPHA_1
#define DEFAULT_ALPHA_1 1e-2
#endif

/*! \def DEFAULT_ALPHA_2
    \brief Parameter tuning the 'lowering' of estimated majorant
*/
#ifndef DEFAULT_ALPHA_2
#define DEFAULT_ALPHA_2 0.99
#endif

// The macros below are unused:
/* #ifndef DEFAULT_ALPHA_3
   #define DEFAULT_ALPHA_3 1e2
   #endif */
/* #ifndef DEFAULT_ALPHA_4
   #define DEFAULT_ALPHA_4 0.99
   #endif */

// NB: the variable names used in this class are not at all intuitive (more
//     meaningful names may be employed).

/*! \class CollisionHandler
 *  \brief Class for collisions simulation
 *
 *  CollisionHandler wraps collisions counters and all methods to intialize majorants,
 *  to compute the expected no. of collisions and to simulate collisions
 */
class CollisionHandler : protected Motherbase
{

private:

  // Global number of collisions
  int n_fake = 0;       /*!< Number of false collisions                 */
  int n_real = 0;       /*!< Number of true collisions                  */
  int n_total = 0;      /*!< Total number  collisions                   */
  int n_fake_idx = 0;   /*!< True collisions for which probability > 1  */

  // Vectors for storage
  std::vector<int> n_fake_store;    /*!< Number of false collisions for each time-step        */
  std::vector<int> n_real_store;    /*!< Number of true collisions for each time-step         */
  std::vector<int> n_total_store;   /*!< Number of total collisions for each time-step        */
  std::vector<int> n_out_store;     /*!< Number of out-of-bound collisions for each time-step */

  // Collisional majorants
  ev_matrix::MaskMatrix<real_number> a11;       /*!< Majorant A_i (density * correlation) */
  ev_matrix::MaskMatrix<real_number> vrmax11;   /*!< Majorant C_i (relative speed)        */

  // New values for the majorants
  ev_matrix::MaskMatrix<real_number> anew;      /*!< Updated values for A_i */
  ev_matrix::MaskMatrix<real_number> vrmaxnew;  /*!< Updated values for C_i */

  // UNUSED
  // ev_matrix::MaskMatrix<real_number> freq11;

  int n_coll;                               /*!< Total number of collisions         */
  ev_matrix::MaskMatrix<int> n_coll_cell;   /*!< Number of collisions for each cell */

  std::valarray<real_number> scaled_k;      /*!< Scaled vector pointing to target cell (sigma*k)  */
  std::valarray<real_number> rel_vel;       /*!< Relative velocity between particles              */
  std::valarray<real_number> delta;         /*!< Vector to be added/subtracted to update velocity */

  /*! \fn void CollisionHandler::gen_scaled_k(void)
      \brief Creates the vector poiting to the cell jc2 (scaled by sigma)
  */
  inline void gen_scaled_k(void)
  {
    rng->sample_unit_sphere(scaled_k[0], scaled_k[1], scaled_k[2]);
    scaled_k *= species->get_diam_fluid();
  }

  // Utilities
  real_number vr = 0.0, scalar_prod = 0.0;
  std::vector<int> cells_ind;

  // Controlling number of collisions
  const real_number alpha_1 = DEFAULT_ALPHA_1;  /*!< First coefficient for collision number control   */
  const real_number alpha_2 = DEFAULT_ALPHA_2;  /*!< Second coefficient for collision number control  */

  // UNUSED
  // const real_number alpha_3 = DEFAULT_ALPHA_3;
  // const real_number alpha_4 = DEFAULT_ALPHA_4;

  /*! \fn void CollisionHandler::update_majorants(void)
      \brief Updates majorants according to the predefined coefficients alpha_1, alpha_1
  */
  void update_majorants(void);

  // References from other classes
  const int& npart;
  const int& nx, ny;
  const real_number& xmin, xmax, ymin, ymax;
  const real_number& rdx, rdy;
  const real_number& sigma, delta_t;

  /*! \fn void CollisionHandler::setup_cell_ind(void)
      \brief Set-up the vector containing the index for each cell

      Cells indices are randomly selected from this array; once the cell is selected
      the index is erased from the array
  */
  void setup_cell_ind(void);

public:

  CollisionHandler(DSMC*);
  ~CollisionHandler() = default;

  // Each step of collisional stage (including initialization)
  void compute_majorants(void);
  void compute_collision_number(void);
  void perform_collisions(void);

  // Collisional stage in a packet
  void perform_collision_kernel(void);

  // GETTERS
  inline int get_n_fake(void) const { return n_fake; }
  inline int get_n_real(void) const { return n_real; }
  inline int get_n_total(void) const { return n_total; }
  inline int get_n_coll(void) const { return n_coll; }
  inline const ev_matrix::MaskMatrix<real_number>& get_a11(void) const { return a11; }
  inline const ev_matrix::MaskMatrix<real_number>& get_vrmax11(void) const { return vrmax11; }
  inline const ev_matrix::MaskMatrix<int>& get_n_coll_cell(void) const { return n_coll_cell; }
  inline std::vector<int>& get_n_fake_store(void) { return n_fake_store; }
  inline const std::vector<int>& get_n_fake_store(void) const { return n_fake_store; }
  inline std::vector<int>& get_n_real_store(void) { return n_real_store; }
  inline const std::vector<int>& get_n_real_store(void) const { return n_real_store; }
  inline std::vector<int>& get_n_total_store(void) { return n_total_store; }
  inline const std::vector<int>& get_n_total_store(void) const { return n_total_store; }
  inline std::vector<int>& get_n_out_store(void) { return n_out_store; }
  inline const std::vector<int>& get_n_out_store(void) const { return n_out_store; }

};

#endif /* EV_COLLISIONS_HPP */
