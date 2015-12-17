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

#ifdef COMPUTE_CLASS

ComputeStyle(meso_prtn_conc/atom,ComputeMesoPrtn_ConcAtom)

#else

#ifndef LMP_COMPUTE_MESO_PRTN_CONC_ATOM_H
#define LMP_COMPUTE_MESO_PRTN_CONC_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeMesoPrtn_ConcAtom : public Compute {
 public:
  ComputeMesoPrtn_ConcAtom(class LAMMPS *, int, char **);
  ~ComputeMesoPrtn_ConcAtom();
  void init();
  void compute_peratom();
  double memory_usage();

 private:
  int nmax;
  double *prtn_concVector;
};

}

#endif
#endif
