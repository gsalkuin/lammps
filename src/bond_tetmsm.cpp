/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Gabriel Alkuino (Syracuse University)

   bond_style tetmsm

   Harmonic spring for tetrahedral mass-spring models.

   Energy:   U = (1/2) kappa * r0 * (1 - r/r0)^2
   Force:    F = -kappa * (r - r0) / r0

   r0 auto-computed from initial geometry on first timestep,
   cached by atom tag pair, serialized to binary restart files.

   Coefficients:  bond_coeff TYPE kappa
------------------------------------------------------------------------- */

#include "bond_tetmsm.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neighbor.h"

#include <cmath>
#include <cstring>
#include <vector>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

BondTetMSM::BondTetMSM(LAMMPS *_lmp) : Bond(_lmp)
{
  writedata = 1;
  max_r0 = 0.0;
}

/* ---------------------------------------------------------------------- */

BondTetMSM::~BondTetMSM()
{
  if (allocated && !copymode) {
    memory->destroy(setflag);
    memory->destroy(kappa);
  }
}

/* ---------------------------------------------------------------------- */

double BondTetMSM::get_r0(int64_t t1, int64_t t2, double r_current)
{
  int64_t lo = (t1 < t2) ? t1 : t2;
  int64_t hi = (t1 < t2) ? t2 : t1;
  BondKey key{lo, hi};

  auto it = r0_map.find(key);
  if (it != r0_map.end()) return it->second;

  if (r_current <= 0.0)
    error->one(FLERR, "Bond tetmsm: zero or negative initial bond length");
  r0_map[key] = r_current;
  if (r_current > max_r0) max_r0 = r_current;
  return r_current;
}

/* ---------------------------------------------------------------------- */

void BondTetMSM::compute(int eflag, int vflag)
{
  int i1, i2, n, type;
  double delx, dely, delz, ebond, fbond;
  double rsq, r, r0, strain;

  ebond = 0.0;
  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  tagint *tag = atom->tag;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (n = 0; n < nbondlist; n++) {
    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    type = bondlist[n][2];

    delx = x[i1][0] - x[i2][0];
    dely = x[i1][1] - x[i2][1];
    delz = x[i1][2] - x[i2][2];

    rsq = delx * delx + dely * dely + delz * delz;
    r = sqrt(rsq);

    r0 = get_r0(tag[i1], tag[i2], r);

    // strain = (r - r0) / r0
    // U      = (1/2) kappa * r0 * strain^2
    // F      = -kappa * strain
    // fbond  = F / r  (applied as f[i] += fbond * del)

    strain = (r - r0) / r0;

    if (eflag) ebond = 0.5 * kappa[type] * r0 * strain * strain;

    if (r > 0.0)
      fbond = -kappa[type] * strain / r;
    else
      fbond = 0.0;

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += fbond * delx;
      f[i1][1] += fbond * dely;
      f[i1][2] += fbond * delz;
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] -= fbond * delx;
      f[i2][1] -= fbond * dely;
      f[i2][2] -= fbond * delz;
    }

    if (evflag) ev_tally(i1, i2, nlocal, newton_bond, ebond, fbond,
                         delx, dely, delz);
  }
}

/* ---------------------------------------------------------------------- */

void BondTetMSM::allocate()
{
  allocated = 1;
  const int np1 = atom->nbondtypes + 1;

  memory->create(kappa, np1, "bond:kappa");
  memory->create(setflag, np1, "bond:setflag");
  for (int i = 1; i < np1; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   bond_coeff TYPE kappa
------------------------------------------------------------------------- */

void BondTetMSM::coeff(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR, "Incorrect args for bond coefficients: "
                            "expected 'bond_coeff TYPE kappa'");
  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->nbondtypes, ilo, ihi, error);

  double kappa_one = utils::numeric(FLERR, arg[1], false, lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    kappa[i] = kappa_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR, "Incorrect args for bond coefficients");
}

/* ---------------------------------------------------------------------- */

double BondTetMSM::equilibrium_distance(int /*type*/)
{
  // if map already populated (from compute or restart), use cached max
  if (max_r0 > 0.0) return max_r0;

  // map not yet populated — scan per-atom bond arrays
  // (available from read_data before neighbor lists are built)

  double **x = atom->x;
  int nlocal = atom->nlocal;
  int *num_bond = atom->num_bond;
  tagint **bond_atom = atom->bond_atom;

  double maxlen = 0.0;
  for (int i = 0; i < nlocal; i++) {
    for (int j = 0; j < num_bond[i]; j++) {
      int k = atom->map(bond_atom[i][j]);
      if (k < 0) continue;
      double dx = x[i][0] - x[k][0];
      double dy = x[i][1] - x[k][1];
      double dz = x[i][2] - x[k][2];
      maxlen = fmax(maxlen, sqrt(dx * dx + dy * dy + dz * dz));
    }
  }

  double maxall;
  MPI_Allreduce(&maxlen, &maxall, 1, MPI_DOUBLE, MPI_MAX, world);
  return maxall;
}

/* ---------------------------------------------------------------------- */

double BondTetMSM::single(int type, double rsq, int i, int j, double &fforce)
{
  double r = sqrt(rsq);
  double r0 = get_r0(atom->tag[i], atom->tag[j], r);
  double strain = (r - r0) / r0;

  fforce = 0.0;
  if (r > 0.0) fforce = -kappa[type] * strain / r;

  return 0.5 * kappa[type] * r0 * strain * strain;
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file:
     1. per-type kappa array
     2. per-bond r0 map (n entries, then [lo, hi, r0] triples)
------------------------------------------------------------------------- */

void BondTetMSM::write_restart(FILE *fp)
{
  fwrite(&kappa[1], sizeof(double), atom->nbondtypes, fp);

  int n = r0_map.size();
  fwrite(&n, sizeof(int), 1, fp);
  for (auto &[key, val] : r0_map) {
    fwrite(&key.lo, sizeof(int64_t), 1, fp);
    fwrite(&key.hi, sizeof(int64_t), 1, fp);
    fwrite(&val, sizeof(double), 1, fp);
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts to all procs
------------------------------------------------------------------------- */

void BondTetMSM::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0)
    utils::sfread(FLERR, &kappa[1], sizeof(double), atom->nbondtypes,
                  fp, nullptr, error);
  MPI_Bcast(&kappa[1], atom->nbondtypes, MPI_DOUBLE, 0, world);

  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;

  // read per-bond r0 map

  int n = 0;
  if (comm->me == 0)
    utils::sfread(FLERR, &n, sizeof(int), 1, fp, nullptr, error);
  MPI_Bcast(&n, 1, MPI_INT, 0, world);

  // each entry: lo (int64) + hi (int64) + r0 (double) = 24 bytes
  const int entry_bytes = 2 * sizeof(int64_t) + sizeof(double);
  std::vector<char> buf(n * entry_bytes);
  if (comm->me == 0 && n > 0)
    utils::sfread(FLERR, buf.data(), entry_bytes, n, fp, nullptr, error);
  if (n > 0)
    MPI_Bcast(buf.data(), n * entry_bytes, MPI_CHAR, 0, world);

  r0_map.clear();
  r0_map.reserve(n);
  char *ptr = buf.data();
  for (int i = 0; i < n; i++) {
    int64_t lo, hi;
    double r0;
    memcpy(&lo, ptr, sizeof(int64_t)); ptr += sizeof(int64_t);
    memcpy(&hi, ptr, sizeof(int64_t)); ptr += sizeof(int64_t);
    memcpy(&r0, ptr, sizeof(double));  ptr += sizeof(double);
    r0_map[{lo, hi}] = r0;
  }

  // rebuild max_r0 from loaded map
  max_r0 = 0.0;
  for (auto &[key, val] : r0_map)
    if (val > max_r0) max_r0 = val;
}

/* ---------------------------------------------------------------------- */

void BondTetMSM::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nbondtypes; i++)
    fprintf(fp, "%d %g\n", i, kappa[i]);
}

/* ---------------------------------------------------------------------- */

void *BondTetMSM::extract(const char *str, int &dim)
{
  dim = 1;
  if (strcmp(str, "kappa") == 0) return (void *) kappa;
  return nullptr;
}
