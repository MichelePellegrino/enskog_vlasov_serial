/*! \file dmsc.hpp
 *  \brief Header containing main class wrapper for DSMC algorithm
 */

#ifndef EV_DSMC_HPP
#define EV_DSMC_HPP

// Including headers for utility libraries (or classes without circular dependencies)
#include "types.hpp"
#include "random.hpp"
#include "stopwatch.hpp"
#include "correlations.hpp"
#include "species.hpp"
#include "times.hpp"
#include "potential.hpp"

// Include STL libraries for (a) pointers, (b) static assert, (c) CPU times map
#include <memory>
#include <cassert>
#include <map>

/*! \def DEFAULT_ITER_THERMO
    \brief Default thermostat iterations (only for testing purposes)
*/
#ifndef DEFAULT_ITER_THERMO
#define DEFAULT_ITER_THERMO 20
#endif

/*! \def DEFAULT_ITER_SAMPLE
    \brief Default sampling iterations (only for testing purposes)
*/
#ifndef DEFAULT_ITER_SAMPLE
#define DEFAULT_ITER_SAMPLE 200
#endif

/*! \def DEFAULT_DUMMY_ITER
    \brief Default total simulation iterations (only for testing purposes)
*/
#ifndef DEFAULT_DUMMY_ITER
#define DEFAULT_DUMMY_ITER 800
#endif

/*!
 *  Stopwatch tags for partial times have been defined with meaningful names
 */

#define DENSITY_TAG   0
#define FORCES_TAG    1
#define ADVECT_TAG    2
#define COLLISION_TAG 3
#define SAMPLING_TAG  4

#define N_TAGS        5

/*!
 *  Forward declarations are needed in order for the design to work: DSMC needs
 *  to encapsulate pointers to all modules, even if referenced classes are constructed
 *  in the source code
 */

class ConfigurationReader;
class Boundary;
class Grid;
class Thermostat;
class Ensemble;
class DensityKernel;
class ForceField;
class CollisionHandler;
class Sampler;
class Output;

template<MarchingType tm_type> class TimeMarching;

/*! \class DSMC
 *  \brief Class for the overall DSMC procedure
 *
 *  DSMC wraps references to all classes implementing modules; moreover, it defines
 *  the functions for the dsmc loop, as well as testing and output functionalities
 */
class DSMC
{

public:
  typedef ev_random::CustomRngObject<RNG> RandomEngine;
  typedef ev_correlations::CorrelationFunction<CORR> CorrelationFun;

private:

  DefaultPointer<ConfigurationReader> conf;               /*!< Configuration reader                     */
  DefaultPointer<RandomEngine> rng;                       /*!< Random numbers generator                 */
  DefaultPointer<Species> species;                        /*!< Parameters for particle species          */
  DefaultPointer<Times> times;                            /*!< Start/end times and intervals            */
  DefaultPointer<Boundary> boundary;                      /*!< Parameters defining boundary cond.       */
  DefaultPointer<Grid> grid;                              /*!< Numerical parameters of the grid         */
  DefaultPointer<Ensemble> ensemble;                      /*!< Particles (storage and population)       */
  DefaultPointer<Thermostat> thermostat;                  /*!< Thermostat (rescaling velocities)        */
  DefaultPointer<DensityKernel> density;                  /*!< Density kernel (storage and computation) */
  DefaultPointer<NondirectionalPairPotential> potential;  /*!< Expression of the long-range potential   */
  DefaultPointer<ForceField> mean_field;                  /*!< Forces kernel (storage and computation)  */
  DefaultPointer<TimeMarching<TM>> time_marching;         /*!< Advection scheme                         */
  DefaultPointer<CollisionHandler> collision_handler;     /*!< Collision simulator, majorants storage   */
  DefaultPointer<Sampler> sampler;                        /*!< Sampling of macroscopic quantities       */
  DefaultPointer<Output> output;                          /*!< Output functionalities                   */
  CorrelationFun correlation;                             /*!< Exp. the short-range correlation         */
  Stopwatch<DefaultWatchPrecision> stopwatch;             /*!< Partial CPU times counter                */

  // NEED TO STORE LOCALLY:
  int n_iter_thermo = DEFAULT_ITER_THERMO;                  /*!< Number of thermostat iterations (stored locally)   */
  int n_iter_sample = DEFAULT_ITER_SAMPLE;                  /*!< Number of sampling iterations (  "  "  )           */
  bool mean_field_gg;                                       /*!< Perform mean-field computation (yes = 1, no = 0)   */
  std::map < int, std::vector<int> > stored_elapsed_times;  /*!< Cumulative elapsed times (see #define tags above)  */

public:

  DSMC(const DefaultString&);
  ~DSMC() = default;

  // GETTERS
  inline DefaultPointer<ConfigurationReader>& get_conf() { return conf; }
  inline DefaultPointer<Species>& get_species() { return species; }
  inline DefaultPointer<Times>& get_times() { return times; }
  inline DefaultPointer<RandomEngine>& get_rng() { return rng; }
  inline DefaultPointer<Boundary>& get_boundary() { return boundary; }
  inline DefaultPointer<Grid>& get_grid() { return grid; }
  inline DefaultPointer<Ensemble>& get_ensemble() { return ensemble; }
  inline DefaultPointer<Thermostat>& get_thermostat() { return thermostat; }
  inline DefaultPointer<DensityKernel>& get_density() { return density; }
  inline DefaultPointer<NondirectionalPairPotential>& get_potential() { return potential; }
  inline DefaultPointer<ForceField>& get_mean_field() { return mean_field; }
  inline DefaultPointer<TimeMarching<TM>>& get_time_marching() { return time_marching; }
  inline DefaultPointer<CollisionHandler>& get_collision_handler() { return collision_handler; }
  inline DefaultPointer<Sampler>& get_sampler() { return sampler; }
  inline DefaultPointer<Output>& get_output() { return output; }
  inline CorrelationFun& get_correlation() { return correlation; }

  // TESTING FEATURES
  void test_boundary_info(void);
  void test_species_info(void);
  void test_thermostat(void);
  void test_density(void);
  void test_force_field(void);
  void test_time_marching(void);
  void test_collisions(void);
  void test_sampling(void);
  void test_output(void);
  void display_barycentre(void) const;
  void display_total_speed(void) const;

  // INITIALIZATION AND DSMC LOOP
  void initialize_simulation(void);
  void dsmc_iteration(void);
  void dsmc_loop(void); /* UNUSED */

  // OUTPUT FEATURES
  void output_all_samples(void);
  void output_all_samples(real_number);
  void output_collision_statistics(void);
  void output_elapsed_times(void);

};

#endif /* DSMC_HPP */
