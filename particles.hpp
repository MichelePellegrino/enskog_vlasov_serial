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

// NB CLEAN THE CODE AND REMOVE UNUSED BUFFERS

/*! \class Ensemble
 *  \brief A class for an ensemble of particles
 */
class Ensemble : protected Motherbase
{

private:

  int n_particles;                    /*!< Total number of particles */
  real_number T_ini;
  real_number mass;

  std::vector<Particle> particles;    /*!< Local vector of particles */

  /*! \fn populate
   *  \brief Randomly populate the phase space
   */
  void populate(void);

public:

  Ensemble(DSMC*);
  ~Ensemble() = default;

  /* FOR TESTING */
  /*! \fn test_stream
   *  \brief Tests particle streaming: periodic b.c. on every edge
   */
  void test_stream(real_number dt = 0.1);

  /* UTILITIES */
  // void save_to_file(const std::string&) const;

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

};


#endif /* EV_PARTICLES_HPP */
