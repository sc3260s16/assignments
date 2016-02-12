#include "energy_force.h"
#include "params.h"
#include "atoms.h"
#include "timer.h"

//************************************************************************
// compute_long_range_correction() function
//   - Calculates long range correction due to finite interaction cutoff.
//   - Arguments:
//       - len_jo: struct containing leonard jones interaction parameters.
//       - m_pars: struct containing misc. parameters.
//       - energy_long: long-range correction to energy.
//       - force_long: long-range correction to force.
//************************************************************************
void compute_long_range_correction(const lj_params * len_jo, misc_params * m_pars,
                                   float * energy_long, float * force_long )
{

   float ulongpre = m_pars->float_N * 8.0 * len_jo->eps * 
                    m_pars->pi * m_pars->density;
   *energy_long = ulongpre * ( len_jo->sig12 / ( 9.0 * len_jo->rcut9 ) - 
                               len_jo->sig6 / ( 6.0 * len_jo->rcut3 ) );

   float vlongpre = 96.0 * len_jo->eps * m_pars->pi * m_pars->density;
   *force_long = -1.0 * vlongpre * ( len_jo->sig12 / ( 9.0 * len_jo->rcut9 ) - 
                                     len_jo->sig6 / ( 6.0 * len_jo->rcut3 ) ); 

}

//************************************************************************
// compute_energy_and_force() function
//   - Calculates energy and force acting on each atom.
//   - Arguments:
//       - myatoms: struct containing all atomic information.
//       - len_jo: struct containing lennard jones interaction parameters.
//       - m_pars: struct containing misc. parameters.
//************************************************************************
void compute_energy_and_force( Atoms * myatoms, const lj_params * len_jo, 
                               misc_params * m_pars )
{

   timeit(1,0);
   int atomi, atomj;
   for (atomi=0; atomi < myatoms->N; atomi++)
   {
      myatoms->fx[atomi] = 0.0;
      myatoms->fy[atomi] = 0.0;
      myatoms->fz[atomi] = 0.0;
   }
   myatoms->pot_energy = 0.0;
   myatoms->virial = 0.0;
   
   for (atomi=0; atomi < myatoms->N; atomi++)
   {

      for (atomj=atomi+1 ; atomj < myatoms->N; atomj++)
      {

         float xxi = myatoms->xx[atomi] - myatoms->xx[atomj];
	 ++m_pars->flop; 
         xxi = minimum_image( xxi, m_pars->side, m_pars->sideh, &m_pars->flop );
         float yyi = myatoms->yy[atomi] - myatoms->yy[atomj];
	 ++m_pars->flop; 
         yyi = minimum_image( yyi, m_pars->side, m_pars->sideh, &m_pars->flop  );
         float zzi = myatoms->zz[atomi] - myatoms->zz[atomj];
	 ++m_pars->flop; 
         zzi = minimum_image( zzi, m_pars->side, m_pars->sideh, &m_pars->flop  );

         float dis2 = xxi*xxi + yyi*yyi + zzi*zzi;
	 ++m_pars->flop; 
	 ++m_pars->flop; 
	 ++m_pars->flop; 
	 ++m_pars->flop; 
	 ++m_pars->flop; 
         if ( dis2 <= len_jo->rcut2 )
         {
            float dis2i = 1.0 / dis2;
	    ++m_pars->flop; 
            float dis6i = dis2i * dis2i * dis2i;
	    ++m_pars->flop; 
	    ++m_pars->flop; 
            float dis12i = dis6i * dis6i;
	    ++m_pars->flop; 
            myatoms->pot_energy += len_jo->sig12 * dis12i - 
                                   len_jo->sig6 * dis6i;
	    ++m_pars->flop; 
	    ++m_pars->flop; 
	    ++m_pars->flop; 
	    ++m_pars->flop; 
            float fterm = dis2i * ( 2.0 * len_jo->sig12 * dis12i -
                                          len_jo->sig6 * dis6i );
	    ++m_pars->flop; 
	    ++m_pars->flop; 
	    ++m_pars->flop; 
	    ++m_pars->flop; 
	    ++m_pars->flop; 
            myatoms->virial -= fterm * dis2;
	    ++m_pars->flop; 
	    ++m_pars->flop; 
            
            myatoms->fx[atomi] += fterm * xxi;
	    ++m_pars->flop; 
	    ++m_pars->flop; 
            myatoms->fy[atomi] += fterm * yyi;
	    ++m_pars->flop; 
	    ++m_pars->flop; 
            myatoms->fz[atomi] += fterm * zzi;
	    ++m_pars->flop; 
	    ++m_pars->flop; 
            myatoms->fx[atomj] -= fterm * xxi;
	    ++m_pars->flop; 
	    ++m_pars->flop; 
            myatoms->fy[atomj] -= fterm * yyi;
	    ++m_pars->flop; 
	    ++m_pars->flop; 
            myatoms->fz[atomj] -= fterm * zzi;
	    ++m_pars->flop; 
	    ++m_pars->flop; 
         
         }

      } 

   }
   for (atomi=0; atomi < myatoms->N; atomi++)
   {
      myatoms->fx[atomi] *= 24.0 * len_jo->eps;
      ++m_pars->flop; 
      ++m_pars->flop; 
      myatoms->fy[atomi] *= 24.0 * len_jo->eps;
      ++m_pars->flop; 
      ++m_pars->flop; 
      myatoms->fz[atomi] *= 24.0 * len_jo->eps;
      ++m_pars->flop; 
      ++m_pars->flop; 
   }
   myatoms->pot_energy *= 4.0 * len_jo->eps;
   ++m_pars->flop; 
   ++m_pars->flop; 
   myatoms->virial *= 24.0 * len_jo->eps;
   ++m_pars->flop; 
   ++m_pars->flop; 
   timeit(1,1);

}

//**********************************************************************
// minimum_image() function
//   - Finds the nearest images of atom i and j, and returns distance.
//   - Arguments:
//       - dist: 1d distance between atoms i and j in central sim. cell.
//       - box_length: length of simulation cell.
//       - half_box_length: half of length of simulation cell.
//**********************************************************************
float minimum_image( const float dist, const float box_length, 
                     const float half_box_length,  long int *flop  )
{

   float min_dist = dist;
   if (dist > half_box_length ) {
     min_dist = dist - box_length; 
     ++flop; 
   }
   if (dist < -half_box_length ) {
     min_dist = dist + box_length;
     ++flop; 
   }
   return min_dist; 

}
