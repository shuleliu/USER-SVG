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

#ifdef FIX_CLASS

FixStyle(electrostatic,FixElectrostatic)

#else

#ifndef LMP_FIX_ELECTROSTATIC_H
#define LMP_FIX_ELECTROSTATIC_H

#include "fix.h"

namespace LAMMPS_NS {

class FixElectrostatic : public Fix {
 public:
  FixElectrostatic(class LAMMPS *, int, char **);
  //~FixElectrostatic();
  int setmask();
  void init();
  void init_list(int, class NeighList *);
  void setup_pre_force(int);
  void pre_force(int);
  //double compute_scalar();
  //int pack_comm(int, int *, double *, int, int *);
  //void unpack_comm(int, int, double *);
  //double memory_usage();
 
 protected:
  double cut_coul,cut_coulsq,alf; 

 private:
  int me;
  double n_avogadro;

  class NeighList *list;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix electrostatic file open failed

The fix electrostatic command could not open a specified file.

E: Fix electrostatic file read failed

The fix electrostatic command could not read the needed parameters from a
specified file.

E: Fix electrostatic found self twice

The neighbor lists used by fix electrostatic are messed up.  If this
error occurs, it is likely a bug, so send an email to the
"developers"_http://lammps.sandia.gov/authors.html.

*/
