// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: William C. Witt
------------------------------------------------------------------------- */

#include "system.hpp"
#include "fix_profess.h"

#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "update.h"

#include <cstring>
#include <iostream>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixProfess::FixProfess(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  scalar_flag = 1;
  global_freq = 1;
  extscalar = 1;
  energy_global_flag = 1;
  virial_global_flag = 1;
  thermo_energy = 1;
  thermo_virial = 1;
}

FixProfess::~FixProfess()
{
}

int FixProfess::setmask()
{
  int mask = 0;
  //mask |= PRE_REVERSE;
  mask |= POST_FORCE;
  //mask |= MIN_POST_FORCE;
  return mask;
}

void FixProfess::post_force(int vflag)
{
  auto system = profess::System::create(
        {{{domain->h[0],          0.0,          0.0},
          {domain->h[5], domain->h[1],          0.0},
          {domain->h[4], domain->h[3], domain->h[2]}}},
        400,
        {"a","ev"});
  std::vector<std::array<double,3>> xyz_coords(atom->nlocal);
  for (int i=0; i<atom->nlocal; i++) {
    xyz_coords[i][0] = atom->x[i][0];
    xyz_coords[i][1] = atom->x[i][1];
    xyz_coords[i][2] = atom->x[i][2];
  }
  system.add_ions("li.gga.recpot", xyz_coords, "a")
        .add_electrons()
        .add_wang_teter_functional()
        .add_hartree_functional()
        .add_ion_electron_functional()
        .add_ion_ion_interaction()
        .add_perdew_burke_ernzerhof_functional();
  profess_energy = system.minimize_energy().energy("ev");
  auto forces = system.forces("ev/a");
  for (int i=0; i<atom->nlocal; i++) {
    atom->f[i][0] = forces[i][0];
    atom->f[i][1] = forces[i][1];
    atom->f[i][2] = forces[i][2];
  }
}

double FixProfess::compute_scalar()
{
  return profess_energy;
}
