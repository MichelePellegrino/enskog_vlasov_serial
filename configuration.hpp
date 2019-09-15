/*! \file configuration.hpp
 *  \brief Header containing the class that reads configuration data from file
 */

#ifndef EV_CONFIGURATION_HPP
#define EV_CONFIGURATION_HPP

// Includes the base class for the DSMC modules
#include "motherbase.hpp"

// Includes the headers for data structures
#include <string>
#include <array>

// Forward declaration for the DSMC wrapper (redundant: may be deleted)
class DSMC;

/*! \class ConfigurationReader
 *  \brief Class used for reading configuration file
 *
 *  ConfigurationReader reads the configuration files passed as constructor parameter
 *  to DSMC: it stores all conf. parameters and allows access from other modules
 *  by getter functions; thid class also set up the cell volume, given the chosen
 *  configuration of the fluid flow.
 */
class ConfigurationReader : protected Motherbase
{

private:

  // FILE NAME
  DefaultString conf_file_name;             /*!< Name of the configuration file */

  // TAGS
  DefaultString numtest;                    /*!< Test serial number */

  // MOLECULAR PROPERTIES
  real_number mass_fluid;                   /*!< Molecular mass of gas particles          */
  real_number diam_fluid;                   /*!< Cross-section of gas particles           */
  char pot_gas;                             /*!< Type of gas-gas potential                */
  char mean_f_gg;                           /*!< Allows ('Y', 'y') gas-gas mean-field     */
  real_number phi11;                        /*!< Well-depth of gas-gas potential          */
  real_number gamma11;                      /*!< Exponent of gas-gas potential            */
  // real_number alpha11;                   /*!< Parameter of gas-gas potential (UNUSED)  */
  real_number mass_solid;                   /*!< Molecular mass of wall particles         */
  real_number diam_solid;                   /*!< Cross-section of wall particles          */
  char pot_wall;                            /*!< Type of gas-wall potential               */
  char mean_f_gw;                           /*!< Allows ('Y', 'y') gas-wall mean-field    */
  real_number phi12;                        /*!< Well-depth of gas-wall potential         */
  real_number gamma12;                      /*!< Exponent of gas-wall potential           */
  // real_number alpha12;                   /*!< Parameter of gas-wall potential (UNUSED) */

  // PHYSICAL STATE OF FLUID-WALLS
  real_number T_ini;                        /*!< Initial fluid temperature                              */
  real_number T_ini1;                       /*!< Initial fluid temperature (UNUSED)                     */
  real_number eta_w0;                       /*!< Wall reduced density, for mean-field computation       */
  real_number eta_w1;                       /*!< Wall reduced density, for collisions simulation        */
  real_number eta_liq0;                     /*!< Reduced density of the gas phase                       */
  real_number eta_liq1;                     /*!< Reduced density of the liquid phase                    */
  int liq_interf;                           /*!< Type of initial gas-liquid configuration               */
  real_number x_liq_interf;		              /*!< Thickness of liquid slab (parallel to x)               */
  real_number y_liq_interf;			            /*!< Thickness of liquid slab (parallel to y)               */
  real_number r_liq_interf;				          /*!< Radius of liquid drop                                  */
  real_number liq_drop_offset;              /*!< Liquid drop offset                                     */
  int niter_thermo;			                    /*!< Number of thermostat iterations                        */
  real_number T_ref;					              /*!< Reference temperature for the thermostat               */
  bool fix_baricentre;				              /*!< Keep baricentre fixed (1) or not (0) (UNUSED)          */
  std::array<char, 4> wall_cond;            /*!< Nature of boundary conditions                          */
  bool set_px;                              /*!< Sets (1) or ignore (0) periodicity in mean-f x [1/0])  */
  bool set_py;                              /*!< Sets (1) or ignore (0) periodicity in mean-f y [1/0])  */
  bool set_eta_px;                          /*!< Sets (1) or ignore (0) periodicity in avdens x [1/0])  */
  bool set_eta_py;                          /*!< Sets (1) or ignore (0) periodicity in avdens y [1/0])  */
  std::array<real_number, 4> T_w;           /*!< Walls temperature                                      */
  real_number f_x_1;                        /*!< Flux coefficient x1-wall (UNUSED)                      */
  real_number f_x_2;                        /*!< Flux coefficient x2-wall (UNUSED)                      */
  real_number f_y_1;		                    /*!< Flux coefficient y1-wall (UNUSED)                      */
  real_number f_y_2;			                  /*!< Flux coefficient x2-wall (UNUSED)                      */
  std::array<real_number, 4> p_e;           /*!< Probability of Maxwellian walls                        */
  std::array<real_number, 4> U_wx;          /*!< X-component of walls velocity                          */
  std::array<real_number, 4> U_wy;          /*!< Y-component of walls velocity                          */
  real_number L_x_1;	                      /*!< Distance between walls and x1-wall                     */
  real_number L_x_2;		                    /*!< Distance between walls and x2-wall                     */
  real_number L_y_1;			                  /*!< Distance between walls and y1-wall                     */
  real_number L_y_2;				                /*!< Distance between walls and y2-wall                     */
  real_number F_x_ext;                      /*!< External field (x-direction)                           */
  // real_number F_y_ext;                   /*!< External field (y-direction)                           */

  // COMPUTATIONAL PARAMETERS
  // int n_Iw;                              /*!< Number of spline points for gw potential (UNUSED)  */
  // int n_If;                              /*!< Number of spline points for gg potential (UNUSED)  */
  real_number x_min;                        /*!< Inferior physical grid boundary, x direction       */
  real_number x_max;			                  /*!< Superior physical grid boundary, x direction       */
  real_number y_min;		                    /*!< Inferior physical grid boundary, y direction       */
  real_number y_max;                        /*!< Superior physical grid boundary, y direction       */
  real_number x_extra;                      /*!< Mean-field cut-off, in x direction                 */
  real_number y_extra;                      /*!< Mean-field cut-off, in y direction                 */
  int n_cells_x;                            /*!< Number of cells in x-direction                     */
  int n_cells_y;                            /*!< Number of cells in y-direction                     */
  int n_part;                               /*!< Number of particles                                */
  bool mean_vel;                            /*!< Fix mean velocity (UNUSED)                         */
  bool restart;                             /*!< Restart from previous configuration (UNUSED)       */
  real_number vx_ini;                       /*!< Prescribed initial velocity in x-direction         */
  real_number vy_ini;                       /*!< Prescribed initial velocity in y-direction         */
  real_number vz_ini;                       /*!< Prescribed initial velocity in z-direction         */
  real_number t_ini;		                    /*!< Initial time                                       */
  real_number t_max;	                      /*!< Final time                                         */
  real_number t_im;                         /*!< Initial time for sampling (UNUSED)                 */
  int nc_out;                               /*!< Number of intervals [t_ini, t_max] (UNUSED)        */
  int iter_out_ini;                         /*!< Start counter for outer-samp. iterations (UNUSED)  */
  real_number delta_t;                      /*!< Time step                                          */
  int qwrite;                               /*!< Tag for writing macro quantities (UNUSED)          */
  int seed;                                 /*!< Seed for RNG                                       */
  int Nv;                                   /*!< Number of nodes for velocity distribution (UNUSED) */
  int routine_choice;                       /*!< Routine for collisions computation (UNUSED)        */
  int split_type;                           /*!< Routine for parallel collisions (UNUSED)           */
  char c_med_comp_type;                     /*!< Routine for mean-field kernel computation (UNUSED) */
  bool collstat;                            /*!< Output (1) or not (0) collisions statistics        */
  int ndom;                                 /*!< Number of subdom. for stat. aggregation (UNUSED)   */
  int niter_sampling;                       /*!< Number of sampling iterations                      */

  // DERIVED PARAMETERS
  /*!
   *  Some numerical parameters have to be derived beforehand, and then passe to the
   *  classes which need them
   */
  real_number dx, dy;

  // INITIAL CONFIGURATION
  /*!
   *  ConfigurationReader defines the parameters characterizing the initial configuration
   *  (in particular: the channel section, cells volume, proportion between the number
   *  of liquid/gas particles)
   */
  real_number channel_area, channel_area0, channel_area1;
  real_number homogeneous_density, homogeneous_density0, homogeneous_density1;
  real_number channel_volume, channel_volume0, channel_volume1;
  real_number channel_section, cell_volume;
  int npart0, npart1;

  // PRESCRIBED DENSITY PROFILE IN X (Y) DIRECTION
  /*!
   *  It is possible to prescribe a desired density profile in x (y) direction;
   *  the profile received as input has to be coherent in its nature and its values
   *  with the values in the conf. file (i.e. liq_interf)
   */
  std::string density_profile_file_name;
  std::vector<real_number> density_profile;
  std::vector<real_number> npc_fraction;

public:

  // CONSTRUCTOR / DESTRUCTIOR

  ConfigurationReader(DSMC*, const DefaultString&);
  ~ConfigurationReader() = default;

  // GETTER METHODS

  inline char get_mean_f_gg() const { return mean_f_gg; }

  inline int get_seed() const { return seed; }
  inline int get_L_x_1() const { return L_x_1; }
  inline int get_L_x_2() const { return L_x_2; }
  inline int get_L_y_1() const { return L_y_1; }
  inline int get_L_y_2() const { return L_y_2; }
  inline const std::array<char, 4>& get_wall_cond() const { return wall_cond; }
  inline const std::array<real_number, 4>& get_p_e() const { return p_e; }

  inline int get_liq_interf() const { return liq_interf; }
  inline real_number get_x_liq_interf() const { return x_liq_interf; }
  inline real_number get_y_liq_interf() const { return y_liq_interf; }

  inline real_number get_diam_fluid() { return diam_fluid; }
  inline real_number get_diam_solid() { return diam_solid; }
  inline real_number get_mass_fluid() { return mass_fluid; }
  inline real_number get_mass_solid() { return mass_solid; }

  inline real_number get_phi11() { return phi11; }
  inline real_number get_gamma11() { return gamma11; }
  inline int get_n_cells_x() const { return n_cells_x; }
  inline int get_n_cells_y() const { return n_cells_y; }
  inline real_number get_x_min() const { return x_min; }
  inline real_number get_x_max() const { return x_max; }
  inline real_number get_y_min() const { return y_min; }
  inline real_number get_y_max() const { return y_max; }
  inline real_number get_x_extra() const { return x_extra; }
  inline real_number get_y_extra() const { return y_extra; }

  inline real_number get_eta_liq0() const { return eta_liq0; }
  inline real_number get_eta_liq1() const { return eta_liq1; }

  inline int get_n_part() const { return n_part; }

  inline int get_nc_out() const { return nc_out; }
  inline int get_restart() const { return restart; }
  inline real_number get_t_ini() const { return t_ini; }
  inline real_number get_t_max() const { return t_max; }
  inline real_number get_t_im() const { return t_im; }
  inline real_number get_delta_t() const { return delta_t; }

  inline int get_niter_thermo() const { return niter_thermo; }
  inline int get_niter_sampling() const { return niter_sampling; }
  inline real_number get_T_ref() const { return T_ref; }
  inline real_number get_T_ini() const { return T_ini; }

  inline real_number get_dx() const { return dx; }
  inline real_number get_dy() const { return dy; }

  inline real_number get_vx_ini() const { return vx_ini; }
  inline real_number get_vy_ini() const { return vy_ini; }
  inline real_number get_vz_ini() const { return vz_ini; }

  inline real_number get_channel_area() const { return channel_area; }
  inline real_number get_channel_area0() const { return channel_area0; }
  inline real_number get_channel_area1() const { return channel_area1; }
  inline real_number get_homogeneous_density() const { return homogeneous_density; }
  inline real_number get_homogeneous_density0() const { return homogeneous_density0; }
  inline real_number get_homogeneous_density1() const { return homogeneous_density1; }
  inline real_number get_channel_volume() const { return channel_volume; }
  inline real_number get_channel_volume0() const { return channel_volume0; }
  inline real_number get_channel_volume1() const { return channel_volume1; }
  inline real_number get_channel_section() const { return channel_section; }
  inline real_number get_cell_volume() const { return cell_volume; }
  inline int get_npart0() const { return npart0; }
  inline int get_npart1() const { return npart1; }

  inline const std::vector<real_number>& get_density_profile(void) { return density_profile; }
  inline real_number get_density_profile(int k) { return density_profile[k]; }
  inline const std::vector<real_number>& get_npc_fraction(void) { return npc_fraction; }
  inline real_number get_npc_fraction(int k) { return npc_fraction[k]; }

private:

  /*! \fn read_conf_file
   *  \brief Reads the conf. file and stores all parameters
   */
  void read_conf_file( void );

  /*! \fn setup_initial_configuration
   *  \brief Set-up the initial fluid configuration
   */
  void setup_initial_configuration( void );

};

#endif /* EV_CONFIGURATION_HPP */
