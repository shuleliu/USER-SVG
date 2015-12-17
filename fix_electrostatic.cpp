/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Koenraad Janssens and David Olmsted (SNL)
------------------------------------------------------------------------- */

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "mpi.h"
#include "fix_electrostatic.h"
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "comm.h"
#include "output.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

#define BIG 1000000000

/* ---------------------------------------------------------------------- */

FixElectrostatic::FixElectrostatic(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  MPI_Comm_rank(world,&me);
  atom->elec_pot_flag = 1;
  n_avogadro = 6.022e23;

  if (narg != 5) error->all(FLERR,"Illegal fix electrostatic command");

  scalar_flag = 1;
  global_freq = 1;
  extscalar = 1;

  peratom_flag = 1;
  size_peratom_cols = 2;
  peratom_freq = 1;

  alf = atof(arg[3]); // decay for wolf method for when I was using it.
  cut_coul = atof(arg[4]);

  //set comm size, added by S. Liu
 
  comm_forward=1;
  // initializations

  cut_coulsq = cut_coul * cut_coul;

}

/* ---------------------------------------------------------------------- */

int FixElectrostatic::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= THERMO_ENERGY;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixElectrostatic::init()
{

  int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

/* ---------------------------------------------------------------------- */

void FixElectrostatic::init_list(int id, NeighList *ptr)
{
  printf ("init_list being called\n");
  list = ptr;
}

/* ---------------------------------------------------------------------- */
void FixElectrostatic::setup_pre_force(int vflag)
{
  pre_force(vflag);
}


void FixElectrostatic::pre_force(int vflag)
{

  // THE EFLAG VFLAG STUFF IS LEFT IN HERE, BUT I SHOULD NEVER NEED IT
  // SINCE THIS PART DOESN'T AFFECT THE ENERGY DIRECTLY (YET!)
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double qtmp_i,qtmp_j,xtmp,ytmp,ztmp,delx,dely,delz,ecoul,fpair;
  double rsq,forcecoul,factor_coul;
  double prefactor;
  double r,rexp;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double erfcc,erfcd,v_sh,dvdrr,e_self,e_shift,f_shift,qisq;

  ecoul = 0.0;
  //if (eflag || vflag) ev_setup(eflag,vflag);
  //else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double *elec_pot = atom->elec_pot;
  double *prtn_conc = atom->prtn_conc;
  double *mass = atom->mass;
  double *rho = atom->rho;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  //double *special_coul = force->special_coul;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;
  double qelectron = force->qelectron;
  
  double r2inv, rinv;
  double print_elec, q_compare, mirror_elec_pot;

  // self and shifted coulombic energy

  v_sh = 0.0;

  //e_shift = erfc(alf*cut_coul)/cut_coul;
  //f_shift = -(e_shift+ 2.0*alf/MY_PIS * exp(-alf*alf*cut_coul*cut_coul)) /  cut_coul;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
//    qtmp_i = prtn_conc[i] * qelectron;
  //  for debug
    qtmp_i = prtn_conc[i]*qelectron*mass[itype]*n_avogadro/rho[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    // WOLF POTENTIAL STUFF, IGNORE
    //qisq = qtmp_i*qtmp_i;
    //e_self = -(e_shift/2.0 + alf/MY_PIS) * qtmp_i*qqrd2e;
    //printf ("e_shift %f alf %f MY_PIS %f qisq %f qqrd2e %f qtmp_i %f xtmp %f\n", e_shift, alf, MY_PIS, qisq, qqrd2e, qtmp_i, xtmp);
    //elec_pot[i] += e_self;
    //if (evflag) ev_tally(i,i,nlocal,0,0.0,e_self,0.0,0.0,0.0,0.0);

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      //factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];
  //    qtmp_j = prtn_conc[j]*qelectron;
  //    for debug
        qtmp_j = prtn_conc[j]*qelectron*mass[jtype]*n_avogadro/rho[j];
      //printf ("i %d qtmp_j %e prtn_conc[j] %e mass[jtype] %e qelectron %e n_avogadro %e\n",atom->tag[j], qtmp_j, prtn_conc[j] , mass[jtype] , qelectron , n_avogadro);

      if (rsq < cut_coulsq) {

        // WOLF POTENTIAL STUFF, IGNORE
        //r = sqrt(rsq);
        //prefactor = qqrd2e*qtmp_j/r;
        //erfcc = erfc(alf*r);
        //erfcd = exp(-alf*alf*r*r);
        //v_sh = (erfcc - e_shift*r) * prefactor;

        r = sqrt(rsq);
        rinv = 1/r;

        v_sh = qqrd2e * qtmp_j * (rinv - 1/cut_coul);
        elec_pot[i] += v_sh;
        print_elec = v_sh;
        if(abs(print_elec)>10000.0)
           {
            printf("abnormal electrostatic potential for i=%d, from j=%d, potential = %e\n",i,j,print_elec);
           }
        //printf ("prenewton i %d %f j %d %f added %f to i - %d %d\n", atom->tag[i] ,elec_pot[i], atom->tag[j], elec_pot[j], print_elec, i , j);

/*        if (newton_pair || j < nlocal) {
          //mirror_elec_pot = e_self;
          mirror_elec_pot = qqrd2e * qtmp_i * (rinv - 1/cut_coul);
          elec_pot[j] += mirror_elec_pot;
          if(abs(mirror_elec_pot)>10000.0)
            {
             printf("abnormal electrostatic potential for j=%d, from i=%d, potential = %e\n",j,i,mirror_elec_pot);
            } 
          //printf ("newton    i %d %f j %d %f added %f to j - %d %d\n", atom->tag[i] ,elec_pot[i], atom->tag[j], elec_pot[j], mirror_elec_pot, i , j);
        }
*/

        // This force related part is left over from copying from a normal electrostatic pair, can be deleted later.
        //dvdrr = (erfcc/rsq + 2.0*alf/MY_PIS * erfcd/r) + f_shift;
        //forcecoul = dvdrr*rsq*prefactor;
        //if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;
        //fpair = forcecoul / rsq;

        //f[i][0] += delx*fpair;
        //f[i][1] += dely*fpair;
        //f[i][2] += delz*fpair;
        //if (newton_pair || j < nlocal) {
        //  f[j][0] -= delx*fpair;
        //  f[j][1] -= dely*fpair;
        //  f[j][2] -= delz*fpair;
        //}

        //if (eflag) {
        //  if (rsq < cut_coulsq) {
        //    ecoul = v_sh;
        //    if (factor_coul < 1.0) ecoul -= (1.0-factor_coul)*prefactor;
        //    } else ecoul = 0.0;
        //}

        //if (evflag) ev_tally(i,j,nlocal,newton_pair,
        //                     0.0,ecoul,fpair,delx,dely,delz);
      }
    }
    //printf("Atom: %d Self %f Shift %f Elec pot: %f alf %f cut_coul %f\n", atom->tag[i] ,e_self, e_shift,  elec_pot[i], alf , cut_coul );
  }

//  if (newton_pair) comm->reverse_comm();
  comm->forward_comm_fix(this);

  //  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */
