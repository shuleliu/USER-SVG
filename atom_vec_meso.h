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

#ifdef ATOM_CLASS

AtomStyle(meso,AtomVecMeso)

#else

#ifndef LMP_ATOM_VEC_MESO_H
#define LMP_ATOM_VEC_MESO_H

#include "atom_vec.h"

namespace LAMMPS_NS {

class AtomVecMeso : public AtomVec {
 public:
  AtomVecMeso(class LAMMPS *, int, char **);
  ~AtomVecMeso() {}
  void grow(int);
  void grow_reset();
  void copy(int, int, int);
  int pack_comm(int, int *, double *, int, int *);
  int pack_comm_vel(int, int *, double *, int, int *);
  void unpack_comm(int, int, double *);
  void unpack_comm_vel(int, int, double *);
  int pack_reverse(int, int, double *);
  void unpack_reverse(int, int *, double *);

  int pack_comm_hybrid(int, int *, double *);
  int unpack_comm_hybrid(int, int, double *);

  int pack_border_hybrid(int, int *, double *);
  int unpack_border_hybrid(int, int, double *);

  int pack_reverse_hybrid(int, int, double *);
  int unpack_reverse_hybrid(int, int *, double *);

  int pack_border(int, int *, double *, int, int *);
  int pack_border_vel(int, int *, double *, int, int *);
  void unpack_border(int, int, double *);
  void unpack_border_vel(int, int, double *);
  int pack_exchange(int, double *);
  int unpack_exchange(double *);
  int size_restart();
  int pack_restart(int, double *);
  int unpack_restart(double *);
  void create_atom(int, double *);
  void data_atom(double *, tagint, char **);
  int data_atom_hybrid(int, char **);
  bigint memory_usage();

 private:
  int *tag,*type,*mask;
  tagint *image;
  double **x,**v,**f;
  double *q;  //charges added by S. Liu
  double *dq;
  int *molecule;
  int **nspecial,**special;
  int *num_bond;
  int **bond_type,**bond_atom;
  double *rho, *drho, *e, *de, *cv, *prtn_conc, *dprtn_conc, *elec_pot;
  double **vest; // estimated velocity during force computation
};

}

#endif
#endif
