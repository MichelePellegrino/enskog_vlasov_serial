#include "dsmc.hpp"
#include "configuration.hpp"
#include "boundary.hpp"
#include "grid.hpp"
#include "particles.hpp"
#include "density.hpp"
#include "force_field.hpp"
#include "collisions.hpp"
#include "advection.hpp"
#include "thermostat.hpp"
#include "sampling.hpp"
#include "output.hpp"

DSMC::DSMC(const DefaultString& file_name):
conf (
  new ConfigurationReader(this, file_name)
),
rng (
  new RandomEngine(conf->get_seed())
),
species (
  new Species(conf->get_diam_fluid(), conf->get_diam_solid(),
    conf->get_mass_fluid(), conf->get_mass_solid())
),
times (
  new Times(conf->get_nc_out(), conf->get_restart(), conf->get_niter_thermo(),
    conf->get_t_ini(), conf->get_t_max(), conf->get_t_im(), conf->get_delta_t())
),
boundary (
  new Boundary(this)
),
grid (
  new Grid(this)
),
ensemble (
  new Ensemble(this)
),
thermostat (
  new Thermostat(this)
),
density (
  new DensityKernel(this)
),
potential (
  new SutherlandMie(conf->get_phi11(), conf->get_diam_fluid(), conf->get_gamma11())
),
mean_field (
  new ForceField(this)
),
time_marching (
  new TimeMarching<TM>(this)
),
collision_handler (
  new CollisionHandler(this)
),
sampler (
  new Sampler(this)
),
output (
  new Output(this)
),
correlation (),
stopwatch(4, "millisecond"),
n_iter_thermo ( conf->get_niter_thermo() ),
n_iter_sample ( conf->get_niter_sampling() )
{

  if ( conf->get_mean_f_gg() == 'y' || conf->get_mean_f_gg() == 'Y' )
  {
    mean_field_gg = true;
    std::cout << "### MEAN-FIELD ON ###" << std::endl;
  }
  else
  {
    mean_field_gg = false;
    std::cout << "### MEAN-FIELD OFF ###" << std::endl;
  }

  // *** PRELIMINARY TESTING ***
  // # # # # #
  // test_species_info();
  // mean_field->testing_output_kernel_function(500);
  // test_boundary_info();
  // test_density();
  // test_force_field();
  // test_thermostat();
  // test_time_marching();
  // test_density();
  // test_collisions();
  // test_sampling();
  // test_output();
  // # # # # #

  std::cout << "### INITIALIZE DSMC SIMULATION ###" << std::endl;
  initialize_simulation();
  test_output();
  display_barycentre();
  display_total_speed();

  std::cout << "### TESTING DSMC ITERATIONS ###" << std::endl;
  int dummy_max_iter = DEFAULT_DUMMY_ITER;
  for (int t = 0; t <= dummy_max_iter; ++t)
  {
    std::cout << " >> iter = " << t << std::endl;
    dsmc_iteration();
    if (t%n_iter_thermo==0)
    {
      std::cout << "    applying thermostat ..." << std::endl;
      thermostat->rescale_velocity();
    }
    if (t%n_iter_sample==0)
    {
      std::cout << "    averaging ..." << std::endl;
      sampler->average();
      output_all_samples(t);
      sampler->reset();
      // output_collision_statistics();
      // output_elapsed_times();
    }
    display_barycentre();
    display_total_speed();
  }
  output->output_collisions_stat();

  std::cout << "### FINALIZING DSMC SIMULATION ###" << std::endl;

}


// TESTING FUNCTIONALITIES

void
DSMC::test_boundary_info
(void)
{
  std::cout << "### TEST: boundary info ###" << std::endl;
  std::cout << "L_x1 = " << boundary->get_Lx1() << ";\tL_y1 = " << boundary->get_Ly1() << std::endl;
  std::cout << "L_x2 = " << boundary->get_Lx2() << ";\tL_y2 = " << boundary->get_Ly2() << std::endl;
  std::cout << "### TEST: grid info ###" << std::endl;
  std::cout << "dx = " << grid->get_dx() << ";\tdy = " << grid->get_dy() << std::endl;
  std::cout << "dx^-1 = " << grid->get_rdx() << ";\tdy^-1 = " << grid->get_rdy() << std::endl;
  std::cout << "n_x = " << grid->get_n_cells_x() << ";\tn_y = " << grid->get_n_cells_y() << std::endl;
  std::cout << "x_max = " << grid->get_x_max() << ";\ty_max = " << grid->get_y_max() << std::endl;
  std::cout << "x_min = " << grid->get_x_min() << ";\ty_min = " << grid->get_y_min() << std::endl;
  std::cout << "x_lp = " << grid->get_xlp() << ";\ty_lp = " << grid->get_ylp() << std::endl;
  std::cout << "x_lm = " << grid->get_xlm() << ";\ty_lm = " << grid->get_ylm() << std::endl;
  std::cout << "cell_volume = " << grid->get_cell_volume() << std::endl;
}

void
DSMC::test_species_info
(void)
{
  std::cout << "### TEST: boundary info ###" << std::endl;
  std::cout << "sigma_g = " << conf->get_diam_fluid() << std::endl;
  std::cout << "mass_g = " << conf->get_mass_fluid() << std::endl;
  std::cout << "phi_g = " << conf->get_phi11() << std::endl;
  std::cout << "gamma_g = " << conf->get_gamma11() << std::endl;
}

void
DSMC::test_thermostat
(void)
{
  std::cout << "### TEST: thermostat ###" << std::endl;
  thermostat->rescale_velocity();
}

void
DSMC::test_density
(void)
{
  std::cout << "### TEST: binning and number density ###" << std::endl;
  density->binning();
  density->fill_dummy_field();
  density->compute_reduced_density();
  std::cout << "### TEST: performing density convolution ###" << std::endl;
  density->compute_avg_density();
}

void
DSMC::test_force_field
(void)
{
  std::cout << "### TEST: computing force field ###" << std::endl;
  mean_field->compute_force_field();
}

void
DSMC::test_time_marching
(void)
{
  std::cout << "### TEST: advection ###" << std::endl;
  time_marching->update_ensemble();
}

void
DSMC::test_collisions
(void)
{
  std::cout << "### TEST: computing collision majorants ###" << std::endl;
  collision_handler->compute_majorants();
  std::cout << "### TEST: computing collsion number ###" << std::endl;
  collision_handler->compute_collision_number();
  std::cout << "### TEST: performing collisions ###" << std::endl;
  collision_handler->perform_collisions();
}

void
DSMC::test_sampling
(void)
{
  std::cout << "### TEST: sampling and averaging ###" << std::endl;
  sampler->sample();
  sampler->average();
}

void
DSMC::test_output
(void)
{
  std::cout << "### TEST: output ###" << std::endl;
  output->output_kernel();
  output->output_forces();
  output->output_density();
  output->output_majorants();
  output->output_collisions();
}

void
DSMC::display_barycentre
(void) const
{
  ensemble->compute_baricentre();
  std::cout << "### TEST: baricentre ###" << std::endl;
  std::cout << " >> x_b = " << ensemble->get_bar_x() << ";\ty_b = " << ensemble->get_bar_y() << std::endl;
}

void
DSMC::display_total_speed
(void) const
{
  ensemble->compute_total_speed();
  std::cout << "### TEST: total speed ###" << std::endl;
  std::cout << " >> v_x = " << ensemble->get_tot_vel_x() << ";\tv_y = " << ensemble->get_tot_vel_y() << std::endl;
}


// DSMC SIMULATION STEPS

void
DSMC::initialize_simulation
(void)
{
  std::cout << "### INITIALIZINING DENSITY FIELD ###" << std::endl;
  density->perform_density_kernel();
  std::cout << "### INITIALIZINING COLLISIONS MAJORANTS ###" << std::endl;
  collision_handler->compute_majorants();
}

void
DSMC::dsmc_iteration
(void)
{
  std::cout << "### PERFORMING DSMC ITERATION ###" << std::endl;
  if ( mean_field_gg )
  {
    std::cout << "    computing force field ..." << std::endl;
    stopwatch.local_start(FORCES_TAG);
    mean_field->compute_force_field();
    stopwatch.local_stop(FORCES_TAG);
    stored_elapsed_times[FORCES_TAG].push_back(stopwatch.get_local_elapsed(FORCES_TAG));
  }
  std::cout << "    propagating ensemble ..." << std::endl;
  stopwatch.local_start(ADVECT_TAG);
  time_marching->update_ensemble();
  stopwatch.local_stop(ADVECT_TAG);
  stored_elapsed_times[ADVECT_TAG].push_back(stopwatch.get_local_elapsed(ADVECT_TAG));
  std::cout << "    computing density ..." << std::endl;
  stopwatch.local_start(DENSITY_TAG);
  density->perform_density_kernel();
  stopwatch.local_stop(DENSITY_TAG);
  stored_elapsed_times[DENSITY_TAG].push_back(stopwatch.get_local_elapsed(DENSITY_TAG));
  std::cout << "    simulating collisions ..." << std::endl;
  stopwatch.local_start(COLLISION_TAG);
  collision_handler->perform_collision_kernel();
  stopwatch.local_stop(COLLISION_TAG);
  stored_elapsed_times[COLLISION_TAG].push_back(stopwatch.get_local_elapsed(COLLISION_TAG));
  std::cout << "    sampling ..." << std::endl;
  stopwatch.local_start(SAMPLING_TAG);
  sampler->sample();
  stopwatch.local_stop(SAMPLING_TAG);
  stored_elapsed_times[SAMPLING_TAG].push_back(stopwatch.get_local_elapsed(SAMPLING_TAG));
}

void
DSMC::dsmc_loop
(void)
{
  // *** TO BE CONTINUED ... ***
}

void
DSMC::output_all_samples
(void)
{
  std::cout << "### OUTPUT ALL SAMPLES ###" << std::endl;
  output->output_sample(sampler->get_vx_avg(), "output_files/samples/test_sample_vx.txt");
  output->output_sample(sampler->get_vy_avg(), "output_files/samples/test_sample_vy.txt");
  output->output_sample(sampler->get_vz_avg(), "output_files/samples/test_sample_vz.txt");
  output->output_sample(sampler->get_temp_avg(), "output_files/samples/test_temp_avg.txt");
  output->output_sample(sampler->get_pxx_avg(), "output_files/samples/test_sample_pxx.txt");
  output->output_sample(sampler->get_pyy_avg(), "output_files/samples/test_sample_pyy.txt");
  output->output_sample(sampler->get_pzz_avg(), "output_files/samples/test_sample_pzz.txt");
  output->output_sample(sampler->get_pxy_avg(), "output_files/samples/test_sample_pxy.txt");
  output->output_sample(sampler->get_pxz_avg(), "output_files/samples/test_sample_pxz.txt");
  output->output_sample(sampler->get_pyz_avg(), "output_files/samples/test_sample_pyz.txt");
  output->output_sample(sampler->get_qx_avg(), "output_files/samples/test_sample_qx.txt");
  output->output_sample(sampler->get_qy_avg(), "output_files/samples/test_sample_qy.txt");
  output->output_sample(sampler->get_qz_avg(), "output_files/samples/test_sample_qz.txt");
  output->output_sample(sampler->get_numdens_avg(), "output_files/samples/test_sample_numdens.txt");
  output->output_sample(sampler->get_forces_x_avg(), "output_files/samples/test_sample_fx.txt");
  output->output_sample(sampler->get_forces_y_avg(), "output_files/samples/test_sample_fy.txt");
}

void
DSMC::output_all_samples
(real_number t)
{
  std::cout << "### OUTPUT ALL SAMPLES ###" << std::endl;
  output->output_sample(sampler->get_vx_avg(), "output_files/samples/test_sample_vx", t);
  output->output_sample(sampler->get_vy_avg(), "output_files/samples/test_sample_vy", t);
  output->output_sample(sampler->get_vz_avg(), "output_files/samples/test_sample_vz", t);
  output->output_sample(sampler->get_temp_avg(), "output_files/samples/test_temp_avg", t);
  output->output_sample(sampler->get_pxx_avg(), "output_files/samples/test_sample_pxx", t);
  output->output_sample(sampler->get_pyy_avg(), "output_files/samples/test_sample_pyy", t);
  output->output_sample(sampler->get_pzz_avg(), "output_files/samples/test_sample_pzz", t);
  output->output_sample(sampler->get_pxy_avg(), "output_files/samples/test_sample_pxy", t);
  output->output_sample(sampler->get_pxz_avg(), "output_files/samples/test_sample_pxz", t);
  output->output_sample(sampler->get_pyz_avg(), "output_files/samples/test_sample_pyz", t);
  output->output_sample(sampler->get_qx_avg(), "output_files/samples/test_sample_qx", t);
  output->output_sample(sampler->get_qy_avg(), "output_files/samples/test_sample_qy", t);
  output->output_sample(sampler->get_qz_avg(), "output_files/samples/test_sample_qz", t);
  output->output_sample(sampler->get_numdens_avg(), "output_files/samples/test_sample_numdens", t);
  output->output_sample(sampler->get_forces_x_avg(), "output_files/samples/test_sample_fx", t);
  output->output_sample(sampler->get_forces_y_avg(), "output_files/samples/test_sample_fy", t);
}

void
DSMC::output_collision_statistics
(void)
{
  output->output_vector(collision_handler->get_n_fake_store(), "output_files/collisions_fake.txt");
  output->output_vector(collision_handler->get_n_real_store(), "output_files/collisions_real.txt");
  output->output_vector(collision_handler->get_n_total_store(), "output_files/collisions_total.txt");
  output->output_vector(collision_handler->get_n_out_store(), "output_files/collisions_out.txt");
}

void
DSMC::output_elapsed_times
(void)
{
  output->output_vector(stored_elapsed_times[DENSITY_TAG], "output_files/times_density.txt");
  output->output_vector(stored_elapsed_times[FORCES_TAG], "output_files/times_forces.txt");
  output->output_vector(stored_elapsed_times[ADVECT_TAG], "output_files/times_advection.txt");
  output->output_vector(stored_elapsed_times[COLLISION_TAG], "output_files/times_collision.txt");
  output->output_vector(stored_elapsed_times[SAMPLING_TAG], "output_files/times_sampling.txt");
}
