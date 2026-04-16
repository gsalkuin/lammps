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

   improper_style tetmsm

   Tetrahedral volume penalty for mass-spring models on tet meshes.

   Energy:   U = (1/2) kappa * V0 * (1 - V/V0)^2
   Force:    f_i = -kappa * (V - V0) / V0 * dV/dx_i

   V0 auto-computed from initial geometry on first timestep,
   cached by atom tag quadruplet, serialized to binary restart files.

   Coefficients:  improper_coeff TYPE kappa
------------------------------------------------------------------------- */

#include "improper_tetmsm.h"

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

ImproperTetMSM::ImproperTetMSM(LAMMPS *_lmp) : Improper(_lmp)
{
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

ImproperTetMSM::~ImproperTetMSM()
{
  if (allocated && !copymode) {
    memory->destroy(setflag);
    memory->destroy(kappa);
  }
}

/* ---------------------------------------------------------------------- */

double ImproperTetMSM::compute_tet_volume(double **x,
                                           int i1, int i2, int i3, int i4)
{
  double ax = x[i2][0] - x[i1][0];
  double ay = x[i2][1] - x[i1][1];
  double az = x[i2][2] - x[i1][2];

  double bx = x[i3][0] - x[i1][0];
  double by = x[i3][1] - x[i1][1];
  double bz = x[i3][2] - x[i1][2];

  double cx = x[i4][0] - x[i1][0];
  double cy = x[i4][1] - x[i1][1];
  double cz = x[i4][2] - x[i1][2];

  return (ax * (by * cz - bz * cy) +
          ay * (bz * cx - bx * cz) +
          az * (bx * cy - by * cx)) / 6.0;
}

/* ---------------------------------------------------------------------- */

double ImproperTetMSM::get_v0(int64_t t1, int64_t t2, int64_t t3, int64_t t4,
                               double **x, int i1, int i2, int i3, int i4)
{
  ImpKey key{t1, t2, t3, t4};
  auto it = v0_map.find(key);
  if (it != v0_map.end()) return it->second;

  double vol = compute_tet_volume(x, i1, i2, i3, i4);
  if (vol <= 0.0)
    error->one(FLERR, "Improper tetmsm: non-positive reference volume; "
               "check vertex ordering");
  v0_map[key] = vol;
  return vol;
}

/* ---------------------------------------------------------------------- */

void ImproperTetMSM::compute(int eflag, int vflag)
{
  int i1, i2, i3, i4, n, type;
  double eimproper, f1[3], f2[3], f3[3], f4[3];
  double vb1x, vb1y, vb1z, vb2x, vb2y, vb2z, vb3x, vb3y, vb3z;

  eimproper = 0.0;
  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  tagint *tag = atom->tag;
  int **improperlist = neighbor->improperlist;
  int nimproperlist = neighbor->nimproperlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (n = 0; n < nimproperlist; n++) {
    i1 = improperlist[n][0];
    i2 = improperlist[n][1];
    i3 = improperlist[n][2];
    i4 = improperlist[n][3];
    type = improperlist[n][4];

    double v0 = get_v0(tag[i1], tag[i2], tag[i3], tag[i4],
                        x, i1, i2, i3, i4);

    // Edge vectors from vertex 1 (i1)

    double ax = x[i2][0] - x[i1][0];
    double ay = x[i2][1] - x[i1][1];
    double az = x[i2][2] - x[i1][2];

    double bx = x[i3][0] - x[i1][0];
    double by = x[i3][1] - x[i1][1];
    double bz = x[i3][2] - x[i1][2];

    double cx = x[i4][0] - x[i1][0];
    double cy = x[i4][1] - x[i1][1];
    double cz = x[i4][2] - x[i1][2];

    // Cross products for volume gradients

    double bxc_x = by * cz - bz * cy;
    double bxc_y = bz * cx - bx * cz;
    double bxc_z = bx * cy - by * cx;

    double cxa_x = cy * az - cz * ay;
    double cxa_y = cz * ax - cx * az;
    double cxa_z = cx * ay - cy * ax;

    double axb_x = ay * bz - az * by;
    double axb_y = az * bx - ax * bz;
    double axb_z = ax * by - ay * bx;

    // V = (1/6) a . (b x c)

    double vol = (ax * bxc_x + ay * bxc_y + az * bxc_z) / 6.0;

    // strain_v = (V - V0) / V0
    // U        = (1/2) kappa * V0 * strain_v^2
    // f_i      = -kappa * strain_v * dV/dx_i

    double strain_v = (vol - v0) / v0;

    if (eflag) eimproper = 0.5 * kappa[type] * v0 * strain_v * strain_v;

    // prefactor absorbs 1/6 from dV/dx:
    //   f_i = -kappa * strain_v * (1/6)(cross product)

    double prefactor = -kappa[type] * strain_v / 6.0;

    f2[0] = prefactor * bxc_x;
    f2[1] = prefactor * bxc_y;
    f2[2] = prefactor * bxc_z;

    f3[0] = prefactor * cxa_x;
    f3[1] = prefactor * cxa_y;
    f3[2] = prefactor * cxa_z;

    f4[0] = prefactor * axb_x;
    f4[1] = prefactor * axb_y;
    f4[2] = prefactor * axb_z;

    f1[0] = -(f2[0] + f3[0] + f4[0]);
    f1[1] = -(f2[1] + f3[1] + f4[1]);
    f1[2] = -(f2[2] + f3[2] + f4[2]);

    // apply force to each of 4 atoms

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += f1[0];
      f[i1][1] += f1[1];
      f[i1][2] += f1[2];
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] += f2[0];
      f[i2][1] += f2[1];
      f[i2][2] += f2[2];
    }

    if (newton_bond || i3 < nlocal) {
      f[i3][0] += f3[0];
      f[i3][1] += f3[1];
      f[i3][2] += f3[2];
    }

    if (newton_bond || i4 < nlocal) {
      f[i4][0] += f4[0];
      f[i4][1] += f4[1];
      f[i4][2] += f4[2];
    }

    // virial: LAMMPS improper convention
    //   vb1 = x[i1] - x[i2] = -a
    //   vb2 = x[i3] - x[i2] = b - a
    //   vb3 = x[i4] - x[i3] = c - b

    if (evflag) {
      vb1x = -ax;
      vb1y = -ay;
      vb1z = -az;

      vb2x = bx - ax;
      vb2y = by - ay;
      vb2z = bz - az;

      vb3x = cx - bx;
      vb3y = cy - by;
      vb3z = cz - bz;

      ev_tally(i1, i2, i3, i4, nlocal, newton_bond, eimproper, f1, f3, f4,
               vb1x, vb1y, vb1z, vb2x, vb2y, vb2z, vb3x, vb3y, vb3z);
    }
  }
}

/* ---------------------------------------------------------------------- */

void ImproperTetMSM::allocate()
{
  allocated = 1;
  const int np1 = atom->nimpropertypes + 1;

  memory->create(kappa, np1, "improper:kappa");
  memory->create(setflag, np1, "improper:setflag");
  for (int i = 1; i < np1; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   improper_coeff TYPE kappa
------------------------------------------------------------------------- */

void ImproperTetMSM::coeff(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR, "Incorrect args for improper coefficients: "
                            "expected 'improper_coeff TYPE kappa'");
  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->nimpropertypes, ilo, ihi, error);

  double kappa_one = utils::numeric(FLERR, arg[1], false, lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    kappa[i] = kappa_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR, "Incorrect args for improper coefficients");
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file:
     1. per-type kappa array
     2. per-tet v0 map (n entries, then [t1, t2, t3, t4, v0] tuples)
------------------------------------------------------------------------- */

void ImproperTetMSM::write_restart(FILE *fp)
{
  fwrite(&kappa[1], sizeof(double), atom->nimpropertypes, fp);

  int n = v0_map.size();
  fwrite(&n, sizeof(int), 1, fp);
  for (auto &[key, val] : v0_map) {
    fwrite(&key.t1, sizeof(int64_t), 1, fp);
    fwrite(&key.t2, sizeof(int64_t), 1, fp);
    fwrite(&key.t3, sizeof(int64_t), 1, fp);
    fwrite(&key.t4, sizeof(int64_t), 1, fp);
    fwrite(&val, sizeof(double), 1, fp);
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts to all procs
------------------------------------------------------------------------- */

void ImproperTetMSM::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0)
    utils::sfread(FLERR, &kappa[1], sizeof(double), atom->nimpropertypes,
                  fp, nullptr, error);
  MPI_Bcast(&kappa[1], atom->nimpropertypes, MPI_DOUBLE, 0, world);

  for (int i = 1; i <= atom->nimpropertypes; i++) setflag[i] = 1;

  // read per-tet v0 map

  int n = 0;
  if (comm->me == 0)
    utils::sfread(FLERR, &n, sizeof(int), 1, fp, nullptr, error);
  MPI_Bcast(&n, 1, MPI_INT, 0, world);

  // each entry: 4 tags (int64) + v0 (double) = 40 bytes
  const int entry_bytes = 4 * sizeof(int64_t) + sizeof(double);
  std::vector<char> buf(n * entry_bytes);
  if (comm->me == 0 && n > 0)
    utils::sfread(FLERR, buf.data(), entry_bytes, n, fp, nullptr, error);
  if (n > 0)
    MPI_Bcast(buf.data(), n * entry_bytes, MPI_CHAR, 0, world);

  v0_map.clear();
  v0_map.reserve(n);
  char *ptr = buf.data();
  for (int i = 0; i < n; i++) {
    int64_t t1, t2, t3, t4;
    double v0;
    memcpy(&t1, ptr, sizeof(int64_t)); ptr += sizeof(int64_t);
    memcpy(&t2, ptr, sizeof(int64_t)); ptr += sizeof(int64_t);
    memcpy(&t3, ptr, sizeof(int64_t)); ptr += sizeof(int64_t);
    memcpy(&t4, ptr, sizeof(int64_t)); ptr += sizeof(int64_t);
    memcpy(&v0, ptr, sizeof(double));  ptr += sizeof(double);
    v0_map[{t1, t2, t3, t4}] = v0;
  }
}

/* ---------------------------------------------------------------------- */

void ImproperTetMSM::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nimpropertypes; i++)
    fprintf(fp, "%d %g\n", i, kappa[i]);
}

/* ---------------------------------------------------------------------- */

void *ImproperTetMSM::extract(const char *str, int &dim)
{
  dim = 1;
  if (strcmp(str, "kappa") == 0) return (void *) kappa;
  return nullptr;
}
