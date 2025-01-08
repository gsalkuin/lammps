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

#ifdef BOND_CLASS
// clang-format off
BondStyle(bpm/beam,BondBPMBeam);
// clang-format on
#else

#ifndef LMP_BOND_BPM_BEAM_H
#define LMP_BOND_BPM_BEAM_H

#include "bond_bpm.h"

namespace LAMMPS_NS {

class BondBPMBeam : public BondBPM {
 public:
  BondBPMBeam(class LAMMPS *);
  ~BondBPMBeam() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void init_style() override;
  void settings(int, char **) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;
  double single(int, double, int, int, double &) override;

 protected:
  double *E, *G;
  double *ecrit, *thetacrit, *psicrit;
  double *gnorm, *gslide, *groll, *gtwist;
  int smooth_flag, normalize_flag;

  double elastic_forces(int, int, int, double, double, double, double *, double *, double *,
                        double *, double *, double *, double &);
  void damping_forces(int, int, int, double *, double *, double *, double *, double *);

  void allocate();
  void store_data();
  double store_bond(int, int, int);
};

}    // namespace LAMMPS_NS

#endif
#endif
