/* -*- c++ -*- ----------------------------------------------------------
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
   Reference length r0 is auto-computed from initial geometry and
   stored in binary restart files.

   Energy:   U = (1/2) kappa_E * r0 * (1 - r/r0)^2
   Force:    F = -kappa_E * (r - r0) / r0

   kappa_E has units of [energy/length] and is the spring stiffness
   (strain energy per unit reference length at unit strain).
   The conventional spring constant is kE = kappa_E / r0.

   Coeffs:   bond_coeff TYPE kappa_E
------------------------------------------------------------------------- */

#ifdef BOND_CLASS
// clang-format off
BondStyle(tetmsm,BondTetMSM);
// clang-format on
#else

#ifndef LMP_BOND_TETMSM_H
#define LMP_BOND_TETMSM_H

#include "bond.h"

#include <unordered_map>
#include <cstdint>

namespace LAMMPS_NS {

class BondTetMSM : public Bond {
 public:
  BondTetMSM(class LAMMPS *);
  ~BondTetMSM() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  double equilibrium_distance(int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_data(FILE *) override;
  double single(int, double, int, int, double &) override;
  void *extract(const char *, int &) override;

 protected:
  double *kappa;

  struct BondKey {
    int64_t lo, hi;
    bool operator==(const BondKey &o) const {
      return lo == o.lo && hi == o.hi;
    }
  };
  struct BondKeyHash {
    size_t operator()(const BondKey &k) const {
      size_t h = std::hash<int64_t>()(k.lo);
      h ^= std::hash<int64_t>()(k.hi) + 0x9e3779b9 + (h << 6) + (h >> 2);
      return h;
    }
  };
  std::unordered_map<BondKey, double, BondKeyHash> r0_map;
  double max_r0;

  double get_r0(int64_t t1, int64_t t2, double r_current);

  virtual void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif
