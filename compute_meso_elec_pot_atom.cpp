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
#include "compute_meso_elec_pot_atom.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeMesoElec_PotAtom::ComputeMesoElec_PotAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute meso_elec_pot/atom command");
  if (atom->elec_pot_flag != 1) error->all(FLERR,"compute meso_elec_pot/atom command requires atom_style with electrostatic potential (e.g. meso)");

  peratom_flag = 1;
  size_peratom_cols = 0;

  nmax = 0;
  elec_potVector = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeMesoElec_PotAtom::~ComputeMesoElec_PotAtom()
{
  memory->sfree(elec_potVector);
}

/* ---------------------------------------------------------------------- */

void ComputeMesoElec_PotAtom::init()
{

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"elec_potVector/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute elec_potVector/atom");
}

/* ---------------------------------------------------------------------- */

void ComputeMesoElec_PotAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow elec_potVector array if necessary

  if (atom->nlocal > nmax) {
    memory->sfree(elec_potVector);
    nmax = atom->nmax;
    elec_potVector = (double *) memory->smalloc(nmax*sizeof(double),"atom:elec_potVector");
    vector_atom = elec_potVector;
  }

  // compute electrostatic potential for each atom in group

  double *elec_pot = atom->elec_pot;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
              elec_potVector[i] = elec_pot[i];
              //printf ("atom number:%i proton concentration: %f\n" , i , elec_pot[i]); 
      }
      else {
              elec_potVector[i] = 0.0;
      }
    }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeMesoElec_PotAtom::memory_usage()
{
  double bytes = nmax * sizeof(double);
  return bytes;
}
