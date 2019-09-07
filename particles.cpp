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

  int nx = grid->get_n_cells_x();
  int ny = grid->get_n_cells_y();
  int npc = n_particles / (nx*ny);
  int k = 0;

  switch( conf->get_liq_interf() )
  {
    // See configuration
    case -1:
      for ( int i = 0; i<nx; ++i )
      {
        for ( int j = 0; j<ny; ++j )
        {
          for ( int k_loc = 0; k_loc<npc; ++k_loc )
          {
            particles[k].xp = grid->get_xc(i);
            particles[k].yp = grid->get_yc(j);
            k++;
          }
        }
      }
    break;
    case 0:
      for ( int i = 0; i<n_particles; ++i )
      {
        particles[i].xp = grid->get_x_min() + rng->sample_uniform() * ( grid->get_x_max() - grid->get_x_min() );
        particles[i].yp = grid->get_y_min() + rng->sample_uniform() * ( grid->get_y_max() - grid->get_y_min() );
      }
      break;
    case 5:
      for ( int i = 0; i<conf->get_npart1(); ++i )
      {
        particles[i].xp = grid->get_x_min() + rng->sample_uniform() * ( grid->get_x_max() - grid->get_x_min() );
        particles[i].yp = rng->sample_uniform() * conf->get_y_liq_interf() - conf->get_y_liq_interf()/2.0;
      }
      for ( int i = conf->get_npart1(); i<n_particles; ++i )
      {
        particles[i].xp = grid->get_x_min() + rng->sample_uniform() * ( grid->get_x_max() - grid->get_x_min() );
        do {
          particles[i].yp = grid->get_y_min() + rng->sample_uniform() * ( grid->get_y_max() - grid->get_y_min() );
        } while( abs( particles[i].yp ) < conf->get_y_liq_interf()/2.0 );
      }
      break;
    case 6:
      for ( int i = 0; i<conf->get_npart1(); ++i )
      {
        particles[i].xp = rng->sample_uniform() * conf->get_x_liq_interf() - conf->get_x_liq_interf()/2.0;
        particles[i].yp = grid->get_y_min() + rng->sample_uniform() * ( grid->get_y_max() - grid->get_y_min() );
      }
      for ( int i = conf->get_npart1(); i<n_particles; ++i )
      {
        particles[i].yp = grid->get_y_min() + rng->sample_uniform() * ( grid->get_y_max() - grid->get_y_min() );
        do {
          particles[i].xp = grid->get_x_min() + rng->sample_uniform() * ( grid->get_x_max() - grid->get_x_min() );
        }   while( abs( particles[i].xp ) < conf->get_x_liq_interf()/2.0 );
      }
      break;
    default:
      std::cerr << "Unrecognized configuration" << std::endl;
  }

  for ( int i = 0; i<n_particles; ++i )
  {
    particles[i].cell_x = (int) ( (particles[i].xp - grid->get_x_min() ) / grid->get_dx() );
    particles[i].cell_y = (int) ( (particles[i].yp - grid->get_y_min() ) / grid->get_dy() );
    rng->sample_box_muller (
      mass, 0.0, 0.0, T_ini,
      particles[i].vx,
      particles[i].vy,
      particles[i].vz );
    particles[i].p_tag = i;
  }

}

void
Ensemble::compute_baricentre
(void)
{

  barycentre_x = 0; barycentre_y = 0;

  for ( int i = 0; i<n_particles; ++i )
  {

    barycentre_x += particles[i].xp;
    barycentre_y += particles[i].yp;

  }

  barycentre_x /= n_particles;
  barycentre_y /= n_particles;

}

void
Ensemble::compute_total_speed
(void)
{

  total_speed_x = 0; total_speed_y = 0;

  for ( int i = 0; i<n_particles; ++i )
  {

    total_speed_x += particles[i].vx;
    total_speed_y += particles[i].vy;

  }

  total_speed_x /= n_particles;
  total_speed_y /= n_particles;

}


/*
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
*/
