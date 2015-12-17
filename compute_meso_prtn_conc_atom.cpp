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

#include "string.h"
#include "compute_meso_prtn_conc_atom.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeMesoPrtn_ConcAtom::ComputeMesoPrtn_ConcAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute meso_prtn_conc/atom command");
  if (atom->prtn_conc_flag != 1) error->all(FLERR,"compute meso_prtn_conc/atom command requires atom_style with proton concentration (e.g. meso)");

  peratom_flag = 1;
  size_peratom_cols = 0;

  nmax = 0;
  prtn_concVector = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeMesoPrtn_ConcAtom::~ComputeMesoPrtn_ConcAtom()
{
  memory->sfree(prtn_concVector);
}

/* ---------------------------------------------------------------------- */

void ComputeMesoPrtn_ConcAtom::init()
{

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"prtn_concVector/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute prtn_concVector/atom");
}

/* ---------------------------------------------------------------------- */

void ComputeMesoPrtn_ConcAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow prtn_concVector array if necessary

  if (atom->nlocal > nmax) {
    memory->sfree(prtn_concVector);
    nmax = atom->nmax;
    prtn_concVector = (double *) memory->smalloc(nmax*sizeof(double),"atom:prtn_concVector");
    vector_atom = prtn_concVector;
  }

  // compute proton concentration for each atom in group

  double *prtn_conc = atom->prtn_conc;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
              prtn_concVector[i] = prtn_conc[i];
              //printf ("atom number:%i proton concentration: %f\n" , i , prtn_conc[i]); 
      }
      else {
              prtn_concVector[i] = 0.0;
      }
    }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeMesoPrtn_ConcAtom::memory_usage()
{
  double bytes = nmax * sizeof(double);
  return bytes;
}
