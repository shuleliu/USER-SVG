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
#include "pair_sph_nafion.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "domain.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairSPHNafion::PairSPHNafion(LAMMPS *lmp) : Pair(lmp)
{
  restartinfo = 0;
  int seed = 557438290; //Given seed now for debugging purposes, should be made an input setting to the pair command 

  random = new RanPark(lmp,seed);
  printf ("Random seed being initialised.\n");
}

/* ---------------------------------------------------------------------- */

PairSPHNafion::~PairSPHNafion() {
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(viscosity);
    memory->destroy(kappa_therm);
    memory->destroy(pressure_scale);
  }
  delete random;

}

/* ---------------------------------------------------------------------- */

void PairSPHNafion::compute(int eflag, int vflag) {
  // The form of the force evaluation comes from Vasquez-Quesada, Ellero, Espanol, JCP 2009, 130, 034901

  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, fpair;

  int *ilist, *jlist, *numneigh, **firstneigh;
  double vxtmp, vytmp, vztmp, imass, jmass, fi, fj, fvisc_v, fvisc_r, h, ih, ihsq, ihcub;
  double f_elec;
  double rsq, r, wfd, delVdotDelR, mu, deltaE, deltaH, ti, tj, lrc, random_amplitude, dyad_product, random_energy;
  double delx[3] , delv[3], f_rand[3];
  double dwein[3][3];
  double n_avogadro;

  //if (eflag || vflag)
  //  ev_setup(eflag, vflag);
  //else
  //  evflag = vflag_fdotr = 0;

  double **v = atom->vest;
  double **x = atom->x;
  double **f = atom->f;
  double *rho = atom->rho;
  double *mass = atom->mass;
  double *de = atom->de;
  double *e = atom->e;
  double *cv = atom->cv;
  double *drho = atom->drho;
  double *prtn_conc = atom->prtn_conc;
  double *elec_pot = atom->elec_pot;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  double dt = update-> dt;
  double boltz = force->boltz;
  double qelectron=force->qelectron;
  //double boltz = 0.;

  //double sqrdt = sqrt(dt * 0.5);

  n_avogadro=6.022e23;

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
    vxtmp = v[i][0];
    vytmp = v[i][1];
    vztmp = v[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    imass = mass[itype];

    // compute pressure of particle i with my EOS
    myEOS(rho[i], e[i], cv[i], &fi, &ti);
//printf ("%e\n",fi);
    if (ti < 1e-2)
      printf ("Low temp:Atomi %d temp %e\n", atom->tag[i], ti);
    fi /= (rho[i] * rho[i]);
//printf ("%e\n",fi);
    //printf("fi = %f\n", fi);
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx[0] = xtmp - x[j][0];
      delx[1] = ytmp - x[j][1];
      delx[2] = ztmp - x[j][2];
      rsq = delx[0] * delx[0] + delx[1] * delx[1] + delx[2] * delx[2];
      r = sqrt(rsq);
      jtype = type[j];
      jmass = mass[jtype];

      if (rsq < cutsq[itype][jtype]) {
        h = cut[itype][jtype];
        ih = 1.0 / h;
        ihsq = ih * ih;
        ihcub = ihsq * ih;

        delv[0] = vxtmp - v[j][0];
        delv[1] = vytmp - v[j][1];
        delv[2] = vztmp - v[j][2];
        wfd = h - r;
        if (domain->dimension == 3) {
          // Lucy Kernel, 3d
          // Note that wfd, the derivative of the weight function with respect to r,
          // is lacking a factor of r.
          // The missing factor of r is recovered by
          // (1) using delV . delX instead of delV . (delX/r) and
          // (2) using f[i][0] += delx * fpair instead of f[i][0] += (delx/r) * fpair
          wfd = -25.066903536973515383e0 * wfd * wfd * ihsq * ihsq * ihsq * ih;
        } else {
          // Lucy Kernel, 2d
          wfd = -19.098593171027440292e0 * wfd * wfd * ihsq * ihsq * ihsq;
        }

        // compute pressure of particle j with my EOS
        myEOS(rho[j], e[j], cv[j], &fj, &tj);
        if (tj < 1e-2)
          printf ("Low tempj:Atom %d temp %e\n", atom->tag[j], tj);
 
        fj /= (rho[j] * rho[j]);

        // apply long-range correction to model a LJ fluid with cutoff
        // this implies that the modelled LJ fluid has cutoff == SPH cutoff
        //lrc = - 11.1701 * (ihcub * ihcub * ihcub - 1.5 * ihcub);
        //fi += lrc;
        //fj += lrc;

        // Scale the pressure forces
        // This is an attempt to include differences between beads repulsion to create net attraction,
        // this is only in the very very early stages of testing out, but setting to 1 in input removes it.
        fi *= pressure_scale[itype][jtype];
        fj *= pressure_scale[itype][jtype];

        // dot product of velocity delta and distance vector
        delVdotDelR = delx[0] * delv[0] + delx[1] * delv[1] + delx[2] * delv[2];
        
        //modified by S. Liu for incompressible fluid, see Eq. 23 of
        //LAMMPS-SPH guide
        fvisc_v = imass * jmass * (viscosity[itype][jtype]+viscosity[jtype][itype])*wfd;
        fvisc_v /= (rho[i] * rho[j]);
/*
        if (delVdotDelR < 0.) {
          fvisc_r = 5 * viscosity[itype][jtype] / 3 ;
          fvisc_r *= delVdotDelR;
          fvisc_r /= rsq * rho[i] * rho[j];    
        } else {
          fvisc_r = 0.;
        }*/

        // total pair force & thermal energy increment
        fpair = -imass * jmass * (fi + fj) * wfd;
        //printf ("fi %e fj %e fvisc_r %e wfd %e\n", fi, fj, fvisc_r, wfd);

        //electrostatic force from other particles
        f_elec = qelectron*n_avogadro*prtn_conc[i]*jmass/rho[j]*elec_pot[j]*wfd;
        fpair += f_elec;
        // create symettric wiener matrix
        // This one still has trace, according to paper I'm using, Vasquez-Quesada, Ellero, Espanol, JCP 2009, 130, 034901
        for (int wii = 0 ; wii < 3 ; wii++) {
          for (int wjj = wii ; wjj < 3 ; wjj++) {
            dwein[wii][wjj] = random->gaussian();
            if (wii != wjj)
              dwein[wjj][wii] = dwein[wii][wjj];
          }
        }
        random_amplitude = (-40/3) * viscosity[itype][jtype] * boltz * ti * tj * imass * jmass * wfd/ (ti +tj);
        random_amplitude /= (rho[i] * rho[j] * rsq * dt) ; //rsq is here to cancel r from dwein.rij to give eij, dt is for dt^.5/dt = 1/dt^.5
        random_amplitude = sqrt(random_amplitude);
   //     random_amplitude = 0;  // <- Un commenting this turns off random fluctutions

//printf ("dwein:\n %e %e %e \n%e %e %e \n%e %e %e \nx y z %e %e %e\n", dwein[0][0],dwein[0][1],dwein[0][2],dwein[1][0],dwein[1][1],dwein[1][2],dwein[2][0],dwein[2][1], dwein[2][2],delx,dely,delz); 
//printf ("Amplitude v%e b%e ti%e tj%e mi%e mj%e wfd%e pi%e pj%e\n",viscosity[itype][jtype] , boltz , ti , tj , imass , jmass , wfd, rho[i] , rho[j]);
        f_rand[0] = dwein[0][0] * delx[0] + dwein[0][1] * delx[1] + dwein[0][2] * delx[2];
        f_rand[1] = dwein[1][0] * delx[0] + dwein[1][1] * delx[1] + dwein[1][2] * delx[2];
        f_rand[2] = dwein[2][0] * delx[0] + dwein[2][1] * delx[1] + dwein[2][2] * delx[2];
        //deltaE = -0.5 * fpair * delVdotDelR;

        f[i][0] += delx[0] * fpair + delv[0] * fvisc_v + f_rand[0] * random_amplitude;
//printf ("i %d j %d fpair %e fvisc_r %e fvisc_v %e f_rand0 %e f_rand1 %e f_rand2 %e force %e \n", atom->tag[i], atom->tag[j], fpair, fvisc_r, fvisc_v, f_rand[0] * random_amplitude,f_rand[1] * random_amplitude,f_rand[2] * random_amplitude, f[i][0]); 
        f[i][1] += delx[1] * fpair + delv[1] * fvisc_v + f_rand[1] * random_amplitude;
//printf ("mi%e mj%e fi%e wfd%e rhoi%e rhoj%e\n",imass, jmass, fi, wfd, rho[i],rho[j]);
//printf ("i %d j %d fpair %e fvisc_r %e fvisc_v %e fv %e fr %e fi %e force %e \n", atom->tag[i], atom->tag[j], fpair, fvisc_r, fvisc_v, delv[1] * fvisc_v,-imass * jmass * (-fvisc_r) * wfd * delx[1],-imass * jmass * (fi) * wfd * delx[1], f[i][1]); 
//printf ("f%e e%e dv%e ",f_rand[1] * random_amplitude, f_rand[1] * random_amplitude * delv[1] / 2, delv[1] * f_rand[1] );
        f[i][2] += delx[2] * fpair + delv[2] * fvisc_v + f_rand[2] * random_amplitude;

        // and change in density
        //drho[i] += jmass * delVdotDelR * wfd; //See NOTE below if this is commented out

        dyad_product = f_rand[0] * delv[0] + f_rand[1] * delv[1] + f_rand[2] * delv[2];
//printf ("dyad %e ", dyad_product);

/* NOTE Large problem with using the Verlet algorithm means the velocities used for the energy and density propagation
would be a half step behind the "actual" velocities of the step. The original SPH module gets around this by assuming the 
velocity is smoothly changing and that v(t+dt) ~ vest = v(t) + dt * a(t). This is a pretty good approximation for normal SPH
but unfortunately for SDPD, the velocity is not smoothly changing, and vest tends to be a terrible approximation for v(t+dt) giving 
terrible conservation of energy and density. We need an improved algorithm, but for now, having it run "canonical" and doing
a density sum every step (pair_sph_rhosum) does good enough, though it means viscous force isn't great.
What this means is anything involving de or drho is commented out in pair_sph_nafion
see Litvinov, Ellero, Hu, Adams J Comput Phys 229 (2010) 5457-5464 */

        random_energy = -4 * kappa_therm[itype][jtype] * boltz * ti * tj * wfd * imass * jmass / (rho[i] * rho[j] * dt);
        random_energy = sqrt(random_energy);
        random_energy *= random->gaussian();
        
        deltaE = delv[0] * delv[0] + delv[1] * delv[1] + delv[2] * delv[2];
        if (delVdotDelR < 0.) 
          deltaE += (delVdotDelR * delVdotDelR / rsq);
        deltaE *= -5 * viscosity[itype][jtype] / 6;
        deltaE *= imass * jmass * wfd;
        deltaE /= (rho[i] * rho[j]);
//printf ("%d %d E %e ", atom->tag[i] , atom->tag[j], deltaE);
        deltaE -= random_amplitude * dyad_product / 2;
//printf ("R %e d %e\n", random_amplitude * dyad_product / 2, dyad_product);
        //deltaE += random_energy;
//printf ("%e ",random_energy);
        deltaH = 2 * kappa_therm[itype][jtype] * (ti - tj);
        deltaH *= imass * jmass * wfd;
        deltaH /= (rho[i] * rho[j]);
//printf ("H %e\n",deltaH);
        //de[i] += deltaE;
        //de[i] += random_energy;
        //de[i] += deltaH;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx[0] * fpair + delv[0] * fvisc_v + f_rand[0] * random_amplitude;
          f[j][1] -= delx[1] * fpair + delv[1] * fvisc_v + f_rand[1] * random_amplitude;
          f[j][2] -= delx[2] * fpair + delv[2] * fvisc_v + f_rand[2] * random_amplitude;
          //de[j] += deltaE;
          //de[j] -= random_energy;
          //de[j] -= deltaH;
          //drho[j] += imass * delVdotDelR * wfd;
        }

        //if (evflag)
        //  ev_tally(i, j, nlocal, newton_pair, 0.0, 0.0, fpair, delx, dely, delz);
      }
    }
  }
  //if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
 allocate all arrays
 ------------------------------------------------------------------------- */

void PairSPHNafion::allocate() {
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");

  memory->create(cut, n + 1, n + 1, "pair:cut");
  memory->create(viscosity, n + 1, n + 1, "pair:viscosity");
  memory->create(kappa_therm, n + 1, n + 1, "pair:kappa_therm");
  memory->create(pressure_scale, n + 1, n + 1, "pair:pressure_scale");
}

/* ----------------------------------------------------------------------
 global settings
 ------------------------------------------------------------------------- */

void PairSPHNafion::settings(int narg, char **arg) {
  if (narg != 0)
    error->all(FLERR,
        "Illegal number of setting arguments for pair_style sph/nafion");
}

/* ----------------------------------------------------------------------
 set coeffs for one or more type pairs
 ------------------------------------------------------------------------- */

void PairSPHNafion::coeff(int narg, char **arg) {
  if (narg != 6)
    error->all(FLERR,
        "Incorrect args for pair_style sph/lj coefficients");
  if (!allocated)
    allocate();

  int ilo, ihi, jlo, jhi;
  force->bounds(arg[0], atom->ntypes, ilo, ihi);
  force->bounds(arg[1], atom->ntypes, jlo, jhi);

  double viscosity_one = force->numeric(arg[2]);
  double cut_one = force->numeric(arg[3]);
  double kappa_therm_one = force->numeric(arg[4]);
  double pressure_scale_one = force->numeric(arg[5]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      viscosity[i][j] = viscosity_one;
      printf("setting cut[%d][%d] = %f\n", i, j, cut_one);
      cut[i][j] = cut_one;
      kappa_therm[i][j] = kappa_therm_one;
      pressure_scale[i][j] = pressure_scale_one;
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

double PairSPHNafion::init_one(int i, int j) {

  if (setflag[i][j] == 0) {
    error->all(FLERR,"All pair sph/lj coeffs are not set");
  }

  cut[j][i] = cut[i][j];
  viscosity[j][i] = viscosity[i][j];
  kappa_therm[j][i] = kappa_therm[i][j];
  pressure_scale[j][i] = pressure_scale[i][j];

  return cut[i][j];
}

/* ---------------------------------------------------------------------- */

double PairSPHNafion::single(int i, int j, int itype, int jtype,
    double rsq, double factor_coul, double factor_lj, double &fforce) {
  fforce = 0.0;

  return 0.0;
}


/*double PairSPHNafion::LJEOS2(double rho, double e, double cv) {


  double T = e / cv;
  if (T < 1.e-2) T = 1.e-2;
  //printf("%f %f\n", T, rho);
  double iT = 0.1e1 / T;
  //double itpow1_4 = exp(0.25 * log(iT)); //pow(iT, 0.1e1 / 0.4e1);
  double itpow1_4 = pow(iT, 0.1e1 / 0.4e1);
  double x = rho * itpow1_4;
  double xsq = x * x;
  double xpow3 = xsq * x;
  double xpow4 = xsq * xsq;
  double xpow9 = xpow3 * xpow3 * xpow3;


  return (0.1e1 + rho * (0.3629e1 + 0.7264e1 * x + 0.104925e2 * xsq + 0.11460e2
      * xpow3 + 0.21760e1 * xpow9 - itpow1_4 * itpow1_4 * (0.5369e1 + 0.13160e2
      * x + 0.18525e2 * xsq - 0.17076e2 * xpow3 + 0.9320e1 * xpow4) + iT
      * (-0.3492e1 + 0.18698e2 * x - 0.35505e2 * xsq + 0.31816e2 * xpow3
          - 0.11195e2 * xpow4)) * itpow1_4) * rho * T;
}*/


/* --------------------------------------------------------------------------------------------- */
/* my very very simplistic EOS,
   suggested by Litvinov et al, phys rev E 77, 066703 (2008)
   P = Po (RHO/RHOo)^gamma + b
   where in this case, Po = 1e5 (atmosphereic pressure)
   RHOo = 1000, water density at stp
   b = 0 and gamma = 4, for no particular reason, but lower than 4 doesn't give stable results.
*/

void PairSPHNafion::myEOS(double rho, double e, double cv, double *p, double *t) {
  double T = e/cv; // e and cv are set in data file, I usually set them to give 300K 
  //double T = 300;
  *p = 1e5 * rho * rho * rho * rho / 1e12;
  *t = T;
}

/* ------------------------------------------------------------------------------ */

/* Jirí Kolafa, Ivo Nezbeda
 * "The Lennard-Jones fluid: an accurate analytic and theoretically-based equation of state",
 *  Fluid Phase Equilibria 100 pp. 1-34 (1994) */
/*double PairSPHNafion::LJEOS2(double rho, double e, double cv) {
 double T = e / cv;

 double sT = sqrt(T);
 double isT = 1.0 / sT;
 double dC = -0.063920968 * log(T) + 0.011117524 / T - 0.076383859 / sT
 + 1.080142248 + 0.000693129 * sT;
 double eta = 3.141592654 / 6. * rho * (dC * dC * dC);
 double zHS = (1 + eta * (1 + eta * (1 - eta / 1.5 * (1 + eta))))
 / ((1. - eta) * (1. - eta) * (1. - eta));
 double BC = (((((-0.58544978 * isT + 0.43102052) * isT + .87361369) * isT
 - 4.13749995) * isT + 2.90616279) * isT - 7.02181962) / T + 0.02459877;
 double gammaBH = 1.92907278;

 double sum = ((2.01546797 * 2 + rho * ((-28.17881636) * 3 + rho
 * (28.28313847 * 4 + rho * (-10.42402873) * 5))) + (-19.58371655 * 2
 + rho * (+75.62340289 * 3 + rho * ((-120.70586598) * 4 + rho
 * (+93.92740328 * 5 + rho * (-27.37737354) * 6)))) / sqrt(T)
 + ((29.34470520 * 2 + rho * ((-112.35356937) * 3 + rho * (+170.64908980
 * 4 + rho * ((-123.06669187) * 5 + rho * 34.42288969 * 6))))
 + ((-13.37031968) * 2 + rho * (65.38059570 * 3 + rho
 * ((-115.09233113) * 4 + rho * (88.91973082 * 5 + rho
 * (-25.62099890) * 6)))) / T) / T) * rho * rho;
 return ((zHS + BC / exp(gammaBH * rho * rho) * rho * (1 - 2 * gammaBH * rho
 * rho)) * T + sum) * rho;
 }
*/
