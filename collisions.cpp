/*! \file collisions.cpp
 *  \brief Source code for the class implementing collisions simulation
 */

#include "collisions.hpp"
#include "grid.hpp"
#include "particles.hpp"
#include "density.hpp"
#include "times.hpp"

CollisionHandler::CollisionHandler(DSMC* dsmc):
  Motherbase(dsmc),
  n_fake_store(),
  n_real_store(),
  n_total_store(),
  n_out_store(),
  a11( 0, grid->get_n_cells_x(), 0, grid->get_n_cells_y(), 0.0 ),
  vrmax11( 0, grid->get_n_cells_x(), 0, grid->get_n_cells_y(), 0.0 ),
  anew( 0, grid->get_n_cells_x(), 0, grid->get_n_cells_y(), 0.0 ),
  vrmaxnew( 0, grid->get_n_cells_x(), 0, grid->get_n_cells_y(), 0.0 ),
  n_coll_cell( 0, grid->get_n_cells_x(), 0, grid->get_n_cells_y(), 0 ),
  scaled_k( 0.0, 3 ),
  rel_vel( 0.0, 3 ),
  delta( 0.0, 3 ),
  cells_ind( grid->get_n_cells(), 0 ),
  npart( ensemble->get_n_particles() ),
  nx( grid->get_n_cells_y() ),
  ny( grid->get_n_cells_y() ),
  xmin( grid->get_x_min() ),
  xmax( grid->get_x_max() ),
  ymin( grid->get_y_min() ),
  ymax( grid->get_y_max() ),
  rdx( grid->get_rdx() ),
  rdy( grid->get_rdy() ),
  sigma( species->get_diam_fluid() ),
  delta_t( times->get_delta_t() )
  { }

void
CollisionHandler::compute_majorants
(void)
{
  int idx_p1, ip1, jp1, idx_ck, ick, jck, ichk, jchk, jjp2, jp2;
  real_number xk, yk, xkh, ykh, chi11;
  int np2;
  int ntest = TEST_COEFF_MULT*npart;
  // Reset all majorants
  a11.fill(0.0);
  anew.fill(0.0);
  vrmax11.fill(0.0);
  vrmaxnew.fill(0.0);
  for (int itest = 0; itest<ntest; itest++)
  {
    idx_p1 = (int)( rng->sample_uniform() * npart );
    ip1 = (int)( ( ensemble->get_xp(idx_p1) - xmin ) * rdx );
    jp1 = (int)( ( ensemble->get_yp(idx_p1) - ymin ) * rdy );
    /*
    ip1 = ensemble->get_cx(idx_p1);
    jp1 = ensemble->get_cx(idx_p1);
    */
    gen_scaled_k();
    xk = ensemble->get_xp(idx_p1) - scaled_k[0];
    yk = ensemble->get_yp(idx_p1) - scaled_k[1];
    /* WHAT ABOUT PERIODIC B.C. ??? */
    if ( ( xk >= xmin && xk <= xmax ) && ( yk >= ymin && yk <= ymax ) )
    {
      xkh = xk + scaled_k[0]/2.0;
      ykh = yk + scaled_k[1]/2.0;
      ick = (int)( ( xk - xmin ) * rdx );
      jck = (int)( ( yk - ymin ) * rdy );
      idx_ck = grid->lexico(ick, jck);
      if ( density->get_npc(ick, jck) >= 1 )
      {
        ichk = (int)( (xkh - xmin) * rdx );
        jchk = (int)( (ykh - ymin) * rdy );
        chi11 = correlation( density->get_aveta(ichk, jchk) );
        a11(ip1, jp1) = std::max(a11(ip1, jp1), density->get_numdens(ip1, jp1) * chi11);
        a11(ick, jck) = std::max(a11(ick, jck), density->get_numdens(ick, jck) * chi11);
        anew(ip1, jp1) = a11(ip1, jp1);
        np2 = density->iof(idx_ck);
        jjp2 = np2 + (int)( rng->sample_uniform() * density->get_npc(ick, jck) );
        jp2 = density->ind(jjp2);
        vr = sqrt (
            ev_utility::power<2>( ensemble->get_vx(jp2) - ensemble->get_vx(idx_p1) )
          + ev_utility::power<2>( ensemble->get_vy(jp2) - ensemble->get_vy(idx_p1) )
          + ev_utility::power<2>( ensemble->get_vz(jp2) - ensemble->get_vz(idx_p1) )
        );
        vrmax11(ip1, jp1) = std::max( vrmax11(ip1, jp1), vr );
        vrmax11(ick,jck) = std::max( vrmax11(ick,jck), vr );
        vrmaxnew(ip1, jp1) = vrmax11(ip1, jp1);
        vrmaxnew(ick, jck) = vrmax11(ick, jck);
      }
    }
  }
}

void
CollisionHandler::compute_collision_number
(void)
{
  int nx = grid->get_n_cells_x(), ny = grid->get_n_cells_y();
  real_number pfnc, cnc;
  n_coll_cell = 0;
  n_coll = 0;
  for (int i = 0; i<nx; ++i)
  {
    for (int j = 0; j<ny; ++j)
    {
      // cnc = ev_const::pi2*sigma*sigma*a11(i,j)*vrmax11(i,j)*delta_t;
      cnc = ev_const::pi2*sigma*sigma*a11(i,j)*vrmax11(i,j)*delta_t;
      n_coll_cell(i,j) = (int)cnc;
      pfnc = cnc - n_coll_cell(i,j);
      if ( rng->sample_uniform() < pfnc ) n_coll_cell(i,j)++;
    }
  }
  n_coll = n_coll_cell.sum();
}

void
CollisionHandler::setup_cell_ind
(void)
{
  int nc = grid->get_n_cells();
  for (int i = 0; i<nc; ++i)
    cells_ind[i] = i;
}

void
CollisionHandler::perform_collisions
(void)
{

  int nc = grid->get_n_cells();
  int idx, idx_cell1, i_cell1, j_cell1, idx_p1;
  int idx_cell2, idx_hcell2, i_cell2, j_cell2, idx_p2;
  real_number xk, yk, xkh, ykh, aa, fk;
  // Set-up
  setup_cell_ind();
  compute_collision_number();
  // Reset number of collisions
  n_fake_idx = 0; n_fake = 0; n_real = 0; n_total = 0;
  while ( nc > 0 )
  {
    // (1) Select a cell at random with with equiprobability
    idx = (int)( rng->sample_uniform() * nc );
    idx_cell1 = cells_ind[idx];
    cells_ind[idx] = cells_ind[nc-1];
    nc -= 1;
    i_cell1 = grid->lexico_inv(idx_cell1).first;
    j_cell1 = grid->lexico_inv(idx_cell1).second;
    // (2) Select a particle belonging to cell at random with equiprobability
    for (int i1 = 0; i1 < n_coll_cell(i_cell1, j_cell1); i1++)
    {
      idx_p1 = density->iof(idx_cell1) + (int)( rng->sample_uniform() *
        density->get_npc(i_cell1, j_cell1) );
      idx_p1 = density->ind(idx_p1);
      // (3) Select, at random, the unit k-vector
      gen_scaled_k();
      xk = ensemble->get_xp(idx_p1) - scaled_k[0];
      yk = ensemble->get_yp(idx_p1) - scaled_k[1];
      xkh = xk + scaled_k[0]/2.0;
      ykh = yk + scaled_k[1]/2.0;
      /* BOUNDARY CONDITIONS: we suppose periodic b.c. */
      if ( xk <= xmin || xk >= xmax )
      {
        xk = xk - round( (xk-0.5*(xmax+xmin))/(xmax-xmin) ) * (xmax-xmin);
        xkh = xkh - round( (xkh-0.5*(xmax+xmin))/(xmax-xmin) ) * (xmax-xmin);
      }
      if ( yk <= ymin || yk >= ymax )
      {
        yk = yk - round( (yk-0.5*(ymax+ymin))/(ymax-ymin) ) * (ymax-ymin);
        ykh = ykh - round( (ykh-0.5*(ymax+ymin))/(ymax-ymin) ) * (ymax-ymin);
      }
      i_cell2 = (int)( (xk-xmin)*rdx );
      j_cell2 = (int)( (yk-ymin)*rdy );
      if( ( density->get_npc(i_cell2, j_cell2)>0 ) )
      {
        // (4) Select at random with equiprobability a particle in the cell where the k-vector points
        idx_cell2 = grid->lexico(i_cell2, j_cell2);
        idx_hcell2 = grid->lexico( (int)((xkh-xmin)*rdx), (int)((ykh-ymin)*rdy) );
        idx = density->iof(idx_cell2) + (int)( rng->sample_uniform() * density->get_npc(i_cell2, j_cell2) );
        idx_p2 = density->ind(idx);
        rel_vel[0] = ensemble->get_vx(idx_p2) - ensemble->get_vx(idx_p1);
        rel_vel[1] = ensemble->get_vy(idx_p2) - ensemble->get_vy(idx_p1);
        rel_vel[2] = ensemble->get_vz(idx_p2) - ensemble->get_vz(idx_p1);
        vr = sqrt( rel_vel[0]*rel_vel[0]+rel_vel[1]*rel_vel[1]+rel_vel[2]*rel_vel[2] );
        vrmaxnew(i_cell1, j_cell1) = std::max( vrmaxnew(i_cell1, j_cell1), vr );
        vrmaxnew(i_cell2, j_cell2) = std::max( vrmaxnew(i_cell2, j_cell2), vr );
        scalar_prod = rel_vel[0]*scaled_k[0] + rel_vel[1]*scaled_k[1] + rel_vel[2]*scaled_k[2];
        scalar_prod *= sigma;
        aa = density->get_numdens(i_cell2, j_cell2) * correlation(
          density->get_aveta(grid->lexico_inv(idx_hcell2).first, grid->lexico_inv(idx_hcell2).second) );
        anew(i_cell1, j_cell1) = std::max( anew(i_cell1, j_cell1), aa );
        anew(i_cell2, j_cell2) = std::max( anew(i_cell2, j_cell2),
          density->get_numdens(i_cell1, j_cell1) * aa / density->get_numdens(i_cell2, j_cell2) );
        if( a11(i_cell2, j_cell2) == 0.0 )
          a11(i_cell2, j_cell2) = anew(i_cell2, j_cell2);
        if( scalar_prod > 0.0 )
        {
          // vrmax11 HAS NOT BEEN UPDATED! WHY?
          fk = scalar_prod * aa / ( a11(i_cell1, j_cell1) * vrmax11(i_cell1, j_cell1) );
          if (fk > 1.0)
            n_fake_idx++;
          if ( rng->sample_uniform() < fk )
          {
            n_real++;
            delta = scaled_k*scalar_prod*species->get_diam_fluid();
            // delta = rel_vel*scalar_prod;
            ensemble->get_vx(idx_p1) += delta[0];
            ensemble->get_vy(idx_p1) += delta[1];
            ensemble->get_vz(idx_p1) += delta[2];
            ensemble->get_vx(idx_p2) -= delta[0];
            ensemble->get_vy(idx_p2) -= delta[1];
            ensemble->get_vz(idx_p2) -= delta[2];
          }
          else
          {
            n_fake++;
          }
          n_total++;
        }
      }
    }
  }
  // Display statistics
  std::cout << "collisions performed" << std::endl;
  std::cout << "total = " << n_total << "\t real = " << n_real << "\t fake = " << n_fake << "\t out-range = " << n_fake_idx << std::endl;
  n_fake_store.push_back(n_fake);
  n_real_store.push_back(n_real);
  n_total_store.push_back(n_total);
  n_out_store.push_back(n_fake_idx);

  update_majorants();

}

void
CollisionHandler::update_majorants
(void)
{
  if ( (double)n_fake_idx > alpha_1*(double)n_real )
  {
    a11 = anew;
    vrmax11 = vrmaxnew;
  }
  else
  {
    a11 = alpha_2*a11;
    vrmax11 = alpha_2*vrmax11;
  }
  /*
  if ( (double)n_fake > alpha_3*(double)n_real )
  {
    a11 = alpha_4*a11;
    vrmax11 = alpha_4*vrmax11;
  }
  */
}

void
CollisionHandler::perform_collision_kernel
(void)
{
  // compute_majorants();
  compute_collision_number();
  perform_collisions();
}
