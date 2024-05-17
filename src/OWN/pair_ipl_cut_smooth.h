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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(ipl/cut/smooth,PairIPLCutSmooth);
//PairStyle(lj/cut/sphere,PairLJCutSphere);
// clang-format on
#else

#ifndef LMP_PAIR_IPL_CUT_SMOOTH_H
#define LMP_PAIR_IPL_CUT_SMOOTH_H
//#ifndef LMP_PAIR_LJ_CUT_SPHERE_H
//#define LMP_PAIR_LJ_CUT_SPHERE_H

#include "pair.h"

namespace LAMMPS_NS {

class PairIPLCutSmooth : public Pair {
 public:
  PairIPLCutSmooth(class LAMMPS *);
  ~PairIPLCutSmooth() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;
  void write_data(FILE *) override;
  void write_data_all(FILE *) override;
  double single(int, int, int, int, double, double, double, double &) override;
  void *extract(const char *, int &) override;

 protected:
  double cut_global;
  int    powern;
  double c0, c2, c4;
  double *rmax;
  double **cut;
  double **epsilon;

  virtual void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif
