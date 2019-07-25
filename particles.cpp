#include "particles.hpp"
#include "grid.hpp"
#include "boundary.hpp"
#include "configuration.hpp"

Ensemble::Ensemble
(DSMC* dsmc):
  Motherbase(dsmc),
  n_particles(conf->get_n_part()),
  T_ini(conf->get_T_ini()),
  mass(species->get_mass_fluid()),
  particles(n_particles)
  {
    std::cout << "### POPULATING ENSEMBLE ###" << std::endl;
    populate();
  }

void
Ensemble::populate
(void)
{

  // Basic type
  for ( int i = 0; i<n_particles; ++i )
  {
    particles[i].xp = grid->get_x_min() + rng->sample_uniform() * ( grid->get_x_max() - grid->get_x_min() );
    particles[i].yp = grid->get_y_min() + rng->sample_uniform() * ( grid->get_y_max() - grid->get_y_min() );
    particles[i].cell_x = (int) ( (particles[i].xp-grid->get_x_min() ) / grid->get_dx() );
    particles[i].cell_y = (int) ( (particles[i].yp-grid->get_y_min() ) / grid->get_dy() );
    /*
    particles[i].vx = - 1.0 + 2.0 * rng->sample_uniform();
    particles[i].vy = - 1.0 + 2.0 * rng->sample_uniform();
    particles[i].vz = 0.0;
    */
    rng->sample_box_muller (
      mass, 0.0, 0.0, T_ini,
      particles[i].vx,
      particles[i].vy,
      particles[i].vz );
    particles[i].p_tag = i;
  }

}

void
Ensemble::test_stream
(real_number dt)
{
  // PERIODIC B.C. ON EVERY EDGE
  for ( int i = 0; i<n_particles; ++i )
  {
    particles[i].xp = particles[i].xp + particles[i].vx*dt;
    particles[i].yp = particles[i].yp + particles[i].vy*dt;
    if ( particles[i].xp > grid->get_x_max() )
      particles[i].xp = grid->get_x_min() + particles[i].xp - grid->get_x_max();
    else if ( particles[i].xp < grid->get_x_min() )
      particles[i].xp = grid->get_x_max() - grid->get_x_min() + particles[i].xp;
    if ( particles[i].yp > grid->get_y_max() )
      particles[i].yp = grid->get_y_min() + particles[i].yp - grid->get_y_max();
    else if ( particles[i].yp < grid->get_y_min() )
      particles[i].yp = grid->get_y_max() - grid->get_y_min() + particles[i].yp;
  }
}
