/*! \file particles.hpp
 *  \brief Header containing the structures to store particle and ensemble properties
 */

#ifndef EV_PARTICLES_HPP
#define EV_PARTICLES_HPP

#include "motherbase.hpp"

#include <vector>

/*! \struct Particle
 *  \brief A struct for a single particle
 */
struct Particle
{
  real_number xp, yp;
  real_number vx, vy, vz;
  int cell_x, cell_y;
  int p_tag;
};

/*! \class Ensemble
 *  \brief A class for an ensemble of particles
 */
class Ensemble : protected Motherbase
{

private:

  int n_particles;                    // Total number of particles
  real_number T_ini;
  real_number mass;

  std::vector<Particle> particles;    // Local vector of particles

  real_number barycentre_x;
  real_number barycentre_y;

  real_number total_speed_x;
  real_number total_speed_y;

  /*! \fn populate
   *  \brief Randomly populate the phase space
   */
  void populate(void);

public:

  Ensemble(DSMC*);
  ~Ensemble() = default;

  /* *** TESTING *** */
  /*! \fn test_stream
   *  \brief Tests particle streaming: periodic b.c. on every edge
   */
  void test_stream(real_number dt = 0.1);

  /* UTILITIES (unused) */
  // void save_to_file(const DefaultString&) const;

  void compute_baricentre(void);
  void compute_total_speed(void);

  // Parameter getters
  inline const int& get_n_particles(void) const { return n_particles; }

  // Data getters
  inline std::vector<Particle>& data(void) { return particles; }
  inline const std::vector<Particle>& data(void) const { return particles; }
  inline Particle& data(int i) { return particles[i]; }
  inline const Particle& data(int i) const { return particles[i]; }

  inline const real_number get_xp(int k) const { return particles[k].xp; }
  inline real_number& get_xp(int k) { return particles[k].xp; }
  inline const real_number get_yp(int k) const { return particles[k].yp; }
  inline real_number& get_yp(int k) { return particles[k].yp; }

  inline const real_number get_vx(int k) const { return particles[k].vx; }
  inline real_number& get_vx(int k) { return particles[k].vx; }
  inline const real_number get_vy(int k) const { return particles[k].vy; }
  inline real_number& get_vy(int k) { return particles[k].vy; }
  inline const real_number get_vz(int k) const { return particles[k].vz; }
  inline real_number& get_vz(int k) { return particles[k].vz; }

  inline int get_cx(int k) const { return particles[k].cell_x; }
  inline int get_cy(int k) const { return particles[k].cell_y; }

  inline real_number get_bar_x () const { return barycentre_x; }
  inline real_number get_bar_y () const { return barycentre_y; }

  inline real_number get_tot_vel_x () const { return total_speed_x; }
  inline real_number get_tot_vel_y () const { return total_speed_y; }

};


#endif /* EV_PARTICLES_HPP */
