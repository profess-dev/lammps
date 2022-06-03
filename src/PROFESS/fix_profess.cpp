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
   Contributing author: Christian Negre (LANL)
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

using namespace LAMMPS_NS;
using namespace FixConst;

extern "C" {
  void profess(int *, int *, double *, int *, int *,
             double *, double *, double *, double *,
             double *, double *, double *, int *,
             double *, double *, double *, double *, int * , bool *);
  int profess_abiversion();
}

// the ABIVERSION number here must be kept consistent
// with its counterpart in the PROFESS library and the
// prototype above. We want to catch mismatches with
// a meaningful error messages, as they can cause
// difficult to debug crashes or memory corruption.

#define PROFESS_ABIVERSION 20180622

/* ---------------------------------------------------------------------- */

FixProfess::FixProfess(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (strcmp(update->unit_style,"metal") != 0)
    error->all(FLERR,"Must use units metal with fix profess command");

  if (comm->nprocs != 1)
    error->all(FLERR,"Fix profess currently runs only in serial");

  if (PROFESS_ABIVERSION != profess_abiversion())
    error->all(FLERR,"LAMMPS is linked against incompatible PROFESS library");

  if (narg != 4) error->all(FLERR,"Illegal fix profess command");

  scalar_flag = 1;
  global_freq = 1;
  extscalar = 1;
  energy_global_flag = 1;
  virial_global_flag = 1;
  thermo_energy = thermo_virial = 1;

  // store ID of compute pe/atom used to generate Coulomb potential for PROFESS
  // null pointer means PROFESS will compute Coulombic potential

  coulomb = 0;
  id_pe = nullptr;

  if (strcmp(arg[3],"NULL") != 0) {
    coulomb = 1;
    error->all(FLERR,"Fix profess does not yet support a LAMMPS calculation of a Coulomb potential");

    id_pe = utils::strdup(arg[3]);
    c_pe = modify->get_compute_by_id(id_pe);
    if (!c_pe) error->all(FLERR,"Could not find fix profess compute ID {}", id_pe);
    if (c_pe->peatomflag == 0)
      error->all(FLERR,"Fix profess compute ID does not compute pe/atom");
  }

  // initializations

  nmax = 0;
  qpotential = nullptr;
  fprofess = nullptr;

  profess_energy = 0.0;
}

/* ---------------------------------------------------------------------- */

FixProfess::~FixProfess()
{
  delete[] id_pe;
  memory->destroy(qpotential);
  memory->destroy(fprofess);
}

/* ---------------------------------------------------------------------- */

int FixProfess::setmask()
{
  int mask = 0;
  mask |= PRE_REVERSE;
  mask |= POST_FORCE;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixProfess::init()
{
  // error checks

  if (domain->dimension == 2)
    error->all(FLERR,"Fix profess requires 3d problem");

  if (coulomb) {
    if (atom->q_flag == 0 || force->pair == nullptr || force->kspace == nullptr)
      error->all(FLERR,"Fix profess cannot compute Coulomb potential");

    c_pe = modify->get_compute_by_id(id_pe);
    if (!c_pe) error->all(FLERR,"Could not find fix profess compute ID {}", id_pe);
  }

  // must be fully periodic or fully non-periodic

  if (domain->nonperiodic == 0) pbcflag = 1;
  else if (!domain->xperiodic && !domain->yperiodic && !domain->zperiodic)
    pbcflag = 0;
  else error->all(FLERR,"Fix profess requires 3d simulation");

  // create qpotential & fprofess if needed
  // for now, assume nlocal will never change

  if (coulomb && qpotential == nullptr) {
    memory->create(qpotential,atom->nlocal,"profess:qpotential");
    memory->create(fprofess,atom->nlocal,3,"profess:fprofess");
  }
}

/* ---------------------------------------------------------------------- */

void FixProfess::init_list(int /*id*/, NeighList * /*ptr*/)
{
  // list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixProfess::setup(int vflag)
{
  newsystem = 1;
  post_force(vflag);
  newsystem = 0;
}

/* ---------------------------------------------------------------------- */

void FixProfess::min_setup(int vflag)
{
  newsystem = 1;
  post_force(vflag);
  newsystem = 0;
}

/* ---------------------------------------------------------------------- */

void FixProfess::setup_pre_reverse(int eflag, int vflag)
{
  pre_reverse(eflag,vflag);
}

/* ----------------------------------------------------------------------
   integrate electronic degrees of freedom
------------------------------------------------------------------------- */

void FixProfess::initial_integrate(int /*vflag*/) {}

/* ----------------------------------------------------------------------
   store eflag, so can use it in post_force to tally per-atom energies
------------------------------------------------------------------------- */

void FixProfess::pre_reverse(int eflag, int /*vflag*/)
{
  eflag_caller = eflag;
}

/* ---------------------------------------------------------------------- */

void FixProfess::post_force(int vflag)
{
  int eflag = eflag_caller;
  ev_init(eflag,vflag);

  // compute Coulombic potential = pe[i]/q[i]
  // invoke compute pe/atom
  // wrap with clear/add and trigger pe/atom calculation every step

  if (coulomb) {
    modify->clearstep_compute();

    if (!(c_pe->invoked_flag & Compute::INVOKED_PERATOM)) {
      c_pe->compute_peratom();
      c_pe->invoked_flag |= Compute::INVOKED_PERATOM;
    }

    modify->addstep_compute(update->ntimestep+1);

    double *pe = c_pe->vector_atom;
    double *q = atom->q;
    int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++)
      if (q[i]) qpotential[i] = pe[i]/q[i];
      else qpotential[i] = 0.0;
  }

  // hardwire these unsupported flags for now

  int coulombflag = 0;
  neighflag = 0;

  // set flags used by PROFESS
  // NOTE: PROFESS does not compute per-atom energies or virials

  int flags[6];

  flags[0] = pbcflag;         // 1 for fully periodic, 0 for fully non-periodic
  flags[1] = coulombflag;     // 1 for LAMMPS computes Coulombics, 0 for PROFESS
  flags[2] = eflag_atom;      // 1 to return per-atom energies, 0 for no
  flags[3] = vflag_global && thermo_virial;    // 1 to return global/per-atom
  flags[4] = vflag_atom && thermo_virial;      //   virial, 0 for no
  flags[5] = neighflag;       // 1 to pass neighbor list to PROFESS, 0 for no

  // setup PROFESS arguments

  int natoms = atom->nlocal;
  double *coords = &atom->x[0][0];
  int *type = atom->type;
  int ntypes = atom->ntypes;
  double *mass = &atom->mass[1];
  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;
  double *forces;
  bool professerror = false;
  if (coulomb) forces = &fprofess[0][0];
  else forces = &atom->f[0][0];
  int maxiter = -1;

  profess(flags,&natoms,coords,type,&ntypes,mass,boxlo,boxhi,&domain->xy,
        &domain->xz,&domain->yz,forces,&maxiter,&profess_energy,
        &atom->v[0][0],&update->dt,virial,&newsystem,&professerror);

  if (professerror) error->all(FLERR,"Internal PROFESS problem");

  // sum PROFESS forces to LAMMPS forces
  // e.g. LAMMPS may compute Coulombics at some point

  if (coulomb) {
    double **f = atom->f;
    int nlocal = atom->nlocal;
    for (int i = 0; i < nlocal; i++) {
      f[i][0] += fprofess[i][0];
      f[i][1] += fprofess[i][1];
      f[i][2] += fprofess[i][2];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixProfess::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   integrate electronic degrees of freedom
------------------------------------------------------------------------- */

void FixProfess::final_integrate() {}

/* ---------------------------------------------------------------------- */

void FixProfess::reset_dt()
{
  //dtv = update->dt;
  //dtf = 0.5 * update->dt * force->ftm2v;
}

/* ----------------------------------------------------------------------
   DFTB energy from PROFESS
------------------------------------------------------------------------- */

double FixProfess::compute_scalar()
{
  return profess_energy;
}

/* ----------------------------------------------------------------------
   memory usage of local arrays
------------------------------------------------------------------------- */

double FixProfess::memory_usage()
{
  double bytes = 0.0;
  if (coulomb) bytes += (double)nmax * sizeof(double);
  if (coulomb) bytes += (double)nmax*3 * sizeof(double);
  return bytes;
}
