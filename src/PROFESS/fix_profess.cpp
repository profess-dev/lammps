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
  thermo_energy = thermo_virial = 1;
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
  std::cout << "hello from post_force" << std::endl;
}

double FixProfess::compute_scalar()
{
  return 1.0;
}
