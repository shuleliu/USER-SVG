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

#include "math.h"
#include "stdlib.h"
#include "pair_sph_protondiffusion.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "memory.h"
#include "error.h"
#include "neigh_list.h"
#include "domain.h"
#include "modify.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairSPHProtonDiffusion::PairSPHProtonDiffusion(LAMMPS *lmp) : Pair(lmp)
{
  restartinfo = 0;
  n_avogadro = 6.022e23;
  pre_factor = 0.5 * force->mvv2e;
  printf ("Pair being called\n");
}

/* ---------------------------------------------------------------------- */

PairSPHProtonDiffusion::~PairSPHProtonDiffusion() {
  printf ("Destroyer being called\n");
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
    memory->destroy(alpha);
    memory->destroy(dcoeff);
  }
}

/* ---------------------------------------------------------------------- */

void PairSPHProtonDiffusion::compute(int eflag, int vflag) {
  int i, j, ii, jj, inum, jnum, itype, jtype;
  int ipore;
  double xtmp, ytmp, ztmp,rcyltmp, delx, dely, delz,rtmp,xtmp_red;

  int *ilist, *jlist, *numneigh, **firstneigh;
  double imass, jmass, h, ih, ihsq;
  double rsq, wfd, D, deltaPC, delta_chem_pot;
  double Di,Dj;

  if (eflag || vflag)
    ev_setup(eflag, vflag);
  else
    evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double *prtn_conc = atom->prtn_conc;
  double *dprtn_conc = atom->dprtn_conc;
  double *q = atom->q;
  double *dq = atom->dq;
  double *rho = atom->rho;
  double *e = atom->e;
  double *mass = atom->mass;
  double *elec_pot = atom->elec_pot;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  int *tag = atom->tag;
  double **v = atom->v;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  double kt;
  double qelectron = force->qelectron;
  
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms and do proton dffusion

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    jlist = firstneigh[i];
    jnum = numneigh[i];

    imass = mass[itype];
    //    if (atom->elec_pot_flag == 1)
    //   {
    //     printf(" i= %d, potential = %e\n",i,elec_pot[i]); 
    //  }
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];
      jmass = mass[jtype];

      if (rsq < cutsq[itype][jtype]) {
        h = cut[itype][jtype];
        ih = 1.0 / h;
        ihsq = ih * ih;

        kt = sph_temp; //Needs to be manually input for now while the bead movement is still being validated.
        kt *= force->boltz;

        // kernel function
        wfd = h - sqrt(rsq);
        if (domain->dimension == 3) {
          // Lucy Kernel, 3d
          // Note that wfd, the derivative of the weight function with respect to r,
          // is lacking a factor of r.
          // The missing factor of r is recovered by
          // deltaE, which is missing a factor of 1/r
          wfd = -25.066903536973515383e0 * wfd * wfd * ihsq * ihsq * ihsq * ih;
        } else {
          // Lucy Kernel, 2d
          wfd = -19.098593171027440292e0 * wfd * wfd * ihsq * ihsq * ihsq;
        }

        jmass = mass[jtype];

        Di=0;
        Dj=0;
        //introducing position dependent diffusion coefficient
        //lamellar
        if(diff_geo==1){
        //   if(itype!=3)
           Di = alpha[itype][jtype] + dcoeff[itype][jtype]*fabs(ztmp);

        //   if(jtype!=3)  
           Dj = alpha[itype][jtype] + dcoeff[itype][jtype]*fabs(x[j][2]);  
           }
        //cylinder, the proton gradient is along the z direction
        else if(diff_geo==2){
           rcyltmp = sqrt(xtmp*xtmp + ytmp*ytmp);
        //   if(itype!=3)
           Di = alpha[itype][jtype] + dcoeff[itype][jtype]*rcyltmp;
           rcyltmp = sqrt(x[j][0]*x[j][0] + x[j][1]*x[j][1]);
        //   if(jtype!=3)  
           Dj = alpha[itype][jtype] + dcoeff[itype][jtype]*rcyltmp;
           }
        //cluster-channel, gradient along the x direction
        else if(diff_geo==3){
         //  if(itype!=3){
           ipore = int((xtmp - center0)/x_per);
           xtmp_red = xtmp - (center0 + x_per*double(ipore));
           rtmp = sqrt(xtmp_red*xtmp_red + ytmp*ytmp + ztmp*ztmp); 
           if(rtmp<=rpore){
              Di = alpha[itype][jtype] + dcoeff[itype][jtype]*rtmp;
              }
           else
              {
               //for channel beads, the diffusion coefficient is set to 
               //be the one for random morphology
               Di=alpha[itype][jtype]/2.0;
              }
         //  }
           //diffusion coefficient for j
        //   if(jtype!=3){ 
           ipore = int((x[j][0] - center0)/x_per);
           xtmp_red = x[j][0] - (center0 + x_per*double(ipore));
           rtmp = sqrt(xtmp_red*xtmp_red + x[j][1]*x[j][1] + x[j][2]*x[j][2]); 
           if(rtmp<=rpore){
              Dj = alpha[itype][jtype] + dcoeff[itype][jtype]*rtmp;
              }
           else
              {
               //for channel beads, the diffusion coefficient is set to 
               //be the one for random morphology
               Dj=alpha[itype][jtype]/2.0;
              }
          //   }
             } 
        else{
         //  if(itype!=3)  
           Di = alpha[itype][jtype]; // diffusion coefficient, included since this may become a function of position
         //  if(jtype!=3)  
           Dj = alpha[itype][jtype]; // diffusion coefficient, included since this may become a function of position
        }

        if (Di == 0)
          continue; // Skips everything following if there is going to be no diffusion anyway

        //debug, Di cannot be less than 0
        if(Di<0){
           printf("ztmp = %e\n",ztmp);
           printf("Di = %e\n",Di);
           error->one(FLERR,"Diffusion Coefficient less than zero");
           }

     //   if ((prtn_conc[i] == 0) || (prtn_conc[j] == 0))
       //   error->one(FLERR,"Proton concentration = 0"); //Can't have proton concentration of 0, doesn't make physical sense
       //diffusion term only applied to non-sidechain beads
  //     if((itype!=3)&&(jtype!=3)){
 //         delta_chem_pot = log(prtn_conc[i] / prtn_conc[j]);
   //       delta_chem_pot *= kt;
    //    }

    //    if ((atom->elec_pot_flag == 1)&&(with_elec==1))
    //      { 
    //        delta_chem_pot += qelectron * (elec_pot[i] - elec_pot[j]);
    //      } 

//        deltaPC = jmass / (kt * rho[j]);
//        deltaPC *= delta_chem_pot ;
      //  if((itype!=3)&&(jtype!=3)){
//            deltaPC *= (prtn_conc[i]*Di + prtn_conc[j]*Dj) * wfd;
        //  }
     /*   else if((itype!=3)&&(jtype==3)){
            deltaPC *= prtn_conc[i]*Di*wfd;
            }
        else if((itype==3)&&(jtype!=3)){
            deltaPC *= prtn_conc[j]*Dj*wfd;
            }
        else {
            deltaPC *= 0.0;
            } */
  //begin here with the new integration scheme
  //this shouldn't be deleted afterwards

      deltaPC = jmass/rho[j]*(Di + Dj)*(prtn_conc[i] - prtn_conc[j])*wfd \
        + jmass/(kt*rho[j])*qelectron*(elec_pot[i] - elec_pot[j])*(prtn_conc[i]*Di + prtn_conc[j]*Dj)*wfd;
//    add contribution from water flow, this is probably problematic
//    because  it lets the total proton concentration diverge
//        deltaPC += jmass/rho[j]*prtn_conc[j]*(v[j][0]*delx+v[j][1]*dely+v[j][2]*delz)*wfd;
//         if(abs(deltaPC)>100000)
  //         printf("Abnormal deltaPC for atom %d, deltaPC = %e\n",i,deltaPC);

        // Note here alpha ii not alpha ij. This means one can create "static" beads that don't change concentration
        // but can still transfer concentration to neighboring beads, useful for boundary conditions.
        if (alpha[itype][itype] != 0){ 
          dprtn_conc[i] += deltaPC;
          if(itype!=3){
             dq[i] += deltaPC*qelectron*mass[itype]*n_avogadro/rho[i];
             }
          else{
             dq[i]=0.0;
             }
          }
//        printf ("i: %i j: %i Protojust computed: %f elec_pot[i] %f elec_pot[j] %f kt %f delta_chem_pot %f\n",tag[i], tag[j], prtn_conc[i], dprtn_conc[i], deltaPC, elec_pot[i] , elec_pot[j], kt, delta_chem_pot); 
//        printf ("i: %i j: %i Proton Conc: %e Dproton conc: %e D just computed: %e jmass %e rho %e diff %e wfd %e kt %e delta_chem_pot %e\n",tag[i], tag[j], prtn_conc[i], dprtn_conc[i], deltaPC, jmass, rho[j], D, wfd, kt, delta_chem_pot); 
        if (newton_pair || j < nlocal) {
          if (alpha[jtype][jtype] != 0){
            dprtn_conc[j] -= deltaPC;
            if(jtype!=3){
               dq[j] -= deltaPC*qelectron*mass[jtype]*n_avogadro/rho[j];
               }
            else{
               dq[j]==0.0;
              }
            }
        }

      }
    }
  }
}

/* ----------------------------------------------------------------------
 allocate all arrays
 ------------------------------------------------------------------------- */

void PairSPHProtonDiffusion::allocate() {
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");
  memory->create(cut, n + 1, n + 1, "pair:cut");
  memory->create(alpha, n + 1, n + 1, "pair:alpha");
  memory->create(dcoeff, n + 1, n + 1, "pair:dcoeff");
}

/* ----------------------------------------------------------------------
 global settings
 ------------------------------------------------------------------------- */

void PairSPHProtonDiffusion::settings(int narg, char **arg) {
  if (narg != 6)
    error->all(FLERR,
        "Illegal number of setting arguments for pair_style sph/protondiffusion");
  sph_temp = force -> numeric(arg[0]);
  // variable added by S. Liu to determine whether electrostatics will 
  // be used in proton concentration calculation: 0 means no electrostatics
  // and 1 means with electrostatics
  with_elec = force -> inumeric(arg[1]); 
  // variable added by S. Liu to determin which type of position dependent diffusion
  // constant should be used, 1 represents lamellar and 2 represents cylinder
  // 3 for cluster-channel
  diff_geo = force -> inumeric(arg[2]); 
  //parameters below are for cluster-channel only, for cylinder and lamellar
  //they're set to be zero
  //the x coordinate of the center of first pore
  center0 = force -> numeric(arg[3]); 
  //periodicity in the x direction
  x_per = force -> numeric(arg[4]);
  //radius of the pore
  rpore = force-> numeric(arg[5]);

}

/* ----------------------------------------------------------------------
 set coeffs for one or more type pairs
 ------------------------------------------------------------------------- */

void PairSPHProtonDiffusion::coeff(int narg, char **arg) {
  if (narg != 5)
    error->all(FLERR,"Incorrect number of args for pair_style sph/protondiffusion coefficients");
  if (!allocated)
    allocate();

  int ilo, ihi, jlo, jhi;
  force->bounds(arg[0], atom->ntypes, ilo, ihi);
  force->bounds(arg[1], atom->ntypes, jlo, jhi);

  double alpha_one = force->numeric(arg[2]);
  double cut_one   = force->numeric(arg[3]);
  double dcoeff_one =force->numeric(arg[4]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      //printf("setting cut[%d][%d] = %f\n", i, j, cut_one);
      cut[i][j] = cut_one;
      alpha[i][j] = alpha_one;
      dcoeff[i][j] = dcoeff_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0)
    error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
 init for one type pair i,j and corresponding j,i
 ------------------------------------------------------------------------- */

double PairSPHProtonDiffusion::init_one(int i, int j) {
  printf ("init_one being called\n");
  if (setflag[i][j] == 0) {
    error->all(FLERR,"All pair sph/protondiffusion coeffs are not set");
  }

  cut[j][i] = cut[i][j];
  alpha[j][i] = alpha[i][j];
  dcoeff[j][i] = dcoeff[i][j];

  return cut[i][j];
}

/* ---------------------------------------------------------------------- */

double PairSPHProtonDiffusion::single(int i, int j, int itype, int jtype,
    double rsq, double factor_coul, double factor_lj, double &fforce) {
  fforce = 0.0;

  return 0.0;
}
