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

#ifdef PAIR_CLASS

PairStyle(sph/protondiffusion,PairSPHProtonDiffusion)

#else

#ifndef LMP_PAIR_SPH_PROTONDIFFUSION_H
#define LMP_PAIR_SPH_PROTONDIFFUSION_H

#include "pair.h"

namespace LAMMPS_NS {

class PairSPHProtonDiffusion : public Pair {
 public:
  PairSPHProtonDiffusion(class LAMMPS *);
  virtual ~PairSPHProtonDiffusion();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  virtual double init_one(int, int);
  virtual double single(int, int, int, int, double, double, double, double &);

 protected:
  double n_avogadro;
  double **cut, **alpha, **dcoeff;
  double sph_temp;
  double center0,x_per,rpore;
  int  with_elec;
  int  diff_geo;
  void allocate();

 private:
  double pre_factor;
};

}

#endif
#endif
