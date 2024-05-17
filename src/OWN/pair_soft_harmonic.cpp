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

#include "pair_soft_harmonic.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_special.h"
#include "memory.h"
#include "neigh_list.h"

#include <cstring>

using namespace LAMMPS_NS;
using MathSpecial::powint;
using MathSpecial::square;

/* ---------------------------------------------------------------------- */

PairSoftHarmonic::PairSoftHarmonic(LAMMPS *lmp) :
    Pair(lmp), rmax(nullptr), cut(nullptr), epsilon(nullptr)
{
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairSoftHarmonic::~PairSoftHarmonic()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(rmax);
    memory->destroy(cut);
    memory->destroy(epsilon);
  }
}

/* ---------------------------------------------------------------------- */

void PairSoftHarmonic::compute(int eflag, int vflag)
{
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, rtmp, delx, dely, delz, evdwl, sigma, fpair;
  double rcutsq, rsq, r2inv, ros, forcelj, factor_lj;
  //double rcutsq, rsq, r2inv, rosigma2, rosigma4, sigmaorn, forcelj, factor_lj;
  int *ilist, *jlist, *numneigh, **firstneigh;

  evdwl = 0.0;
  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  double *radius = atom->radius;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    rtmp = radius[i];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];

      if (true) {
      //if (rsq < cutsq[itype][jtype]) {

        // cutsq is maximum cutoff per type. Now compute and apply real cutoff

        sigma = 2.0 * mix_distance(rtmp, radius[j]);
        rcutsq = square(sigma);
        //rcutsq = square(cut[itype][jtype] * sigma);
        if (rsq < rcutsq) {

          //harmonic
          r2inv = 1.0 /rsq;
          ros   = sqrt(rsq)/sigma;
          forcelj = pow(1.0-ros, alpha-1.0)*ros;
          //r2inv = 1.0 / rsq;
          //rosigma2 = rsq/(sigma*sigma);
          //rosigma4 = rosigma2*rosigma2;
          //sigmaorn = powint(sqrt(1.0/rosigma2), powern);

          //forcelj = -1.0*(double)powern*sigmaorn 
          //        + 4.0*c4 * rosigma4
          //        + 2.0*c2 * rosigma2;
          fpair = factor_lj * forcelj * r2inv;
          

          f[i][0] += delx * fpair;
          f[i][1] += dely * fpair;
          f[i][2] += delz * fpair;
          if (newton_pair || j < nlocal) {
            f[j][0] -= delx * fpair;
            f[j][1] -= dely * fpair;
            f[j][2] -= delz * fpair;
          }

          if (eflag) {
            evdwl = pow(1.0-ros, alpha)/alpha;
            //evdwl = sigmaorn + c4*rosigma4 + c2*rosigma2+c0;
            evdwl *= factor_lj;
          }

          if (evflag) ev_tally(i, j, nlocal, newton_pair, evdwl, 0.0, fpair, delx, dely, delz);
        }
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairSoftHarmonic::allocate()
{
  allocated = 1;
  int n = atom->ntypes + 1;

  memory->create(setflag, n, n, "pair:setflag");
  for (int i = 1; i < n; i++)
    for (int j = i; j < n; j++) setflag[i][j] = 0;

  memory->create(cutsq, n, n, "pair:cutsq");

  memory->create(rmax, n, "pair:rmax");
  memory->create(cut, n, n, "pair:cut");
  memory->create(epsilon, n, n, "pair:epsilon");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairSoftHarmonic::settings(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR, "Illegal pair_style command");

  alpha = utils::numeric(FLERR, arg[0], false, lmp);

  cut_global = utils::numeric(FLERR, arg[1], false, lmp);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i, j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairSoftHarmonic::coeff(int narg, char **arg)
{
  if (narg != 3) error->all(FLERR, "Incorrect args for pair coefficients");
  //if (narg < 3 || narg > 4) error->all(FLERR, "Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  double epsilon_one = utils::numeric(FLERR, arg[2], false, lmp);
  double cut_one = cut_global;
  //if (narg == 4) cut_one = utils::numeric(FLERR, arg[3], false, lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      epsilon[i][j] = epsilon_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairSoftHarmonic::init_style()
{
  Pair::init_style();

  if (!atom->radius_flag)
    error->all(FLERR, "Pair style soft/harmonic requires atom attribute radius");
  if (mix_flag != ARITHMETIC)
    error->all(FLERR, "Pair_modify mix must be arithmetic for pair style soft/harmonic");

  // determine max radius per type

  int *type = atom->type;
  double *radius = atom->radius;
  rmax[0] = 0.0;
  for (int itype = 1; itype <= atom->ntypes; ++itype) {
    double rmax_one = 0.0;
    for (int i = 0; i < atom->nlocal; ++i) {
      if (type[i] == itype) rmax_one = MAX(rmax_one, radius[i]);
    }
    MPI_Allreduce(&rmax_one, &rmax[itype], 1, MPI_DOUBLE, MPI_MAX, world);
  }
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairSoftHarmonic::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    epsilon[i][j] = mix_energy(epsilon[i][i], epsilon[j][j], 0.0, 0.0);
    cut[i][j] = mix_distance(cut[i][i], cut[j][j]);
  }

  epsilon[j][i] = epsilon[i][j];
  cut[j][i] = cut[i][j];

  // since cut is a scaled by the mixed diameter, report maximum possible cutoff.

  return cut[i][j] * 2.0 * mix_distance(rmax[i], rmax[j]);
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairSoftHarmonic::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i, j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j], sizeof(int), 1, fp);
      if (setflag[i][j]) {
        fwrite(&epsilon[i][j], sizeof(double), 1, fp);
        fwrite(&cut[i][j], sizeof(double), 1, fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairSoftHarmonic::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i, j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR, &setflag[i][j], sizeof(int), 1, fp, nullptr, error);
      MPI_Bcast(&setflag[i][j], 1, MPI_INT, 0, world);
      if (setflag[i][j]) {
        if (me == 0) {
          utils::sfread(FLERR, &epsilon[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &cut[i][j], sizeof(double), 1, fp, nullptr, error);
        }
        MPI_Bcast(&epsilon[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&cut[i][j], 1, MPI_DOUBLE, 0, world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairSoftHarmonic::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global, sizeof(double), 1, fp);
  fwrite(&alpha,     sizeof(double), 1, fp);
  fwrite(&offset_flag, sizeof(int), 1, fp);
  fwrite(&mix_flag, sizeof(int), 1, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairSoftHarmonic::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    utils::sfread(FLERR, &cut_global, sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &alpha,     sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &offset_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &mix_flag, sizeof(int), 1, fp, nullptr, error);
  }
  MPI_Bcast(&cut_global, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&alpha, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&offset_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&mix_flag, 1, MPI_INT, 0, world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairSoftHarmonic::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++) fprintf(fp, "%d %g\n", i, epsilon[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairSoftHarmonic::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp, "%d %d %g %g\n", i, j, epsilon[i][j], cut[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairSoftHarmonic::single(int i, int j, int itype, int jtype, double rsq,
                               double /*factor_coul*/, double factor_lj, double &fforce)
{
  double r2inv, sigma, ros, forcelj, philj;
  sigma = 2.0 * mix_distance(atom->radius[i], atom->radius[j]);
  r2inv = 1.0 / rsq;
  //r     = sqrt(rsq);
  ros   = sqrt(rsq) / sigma;
  forcelj = pow(1.0 - ros, alpha-1.0) * ros;


  //double r2inv, sigma, rosigma2, rosigma4, sigmaorn, forcelj, philj;

  //sigma = 2.0 * mix_distance(atom->radius[i], atom->radius[j]);
  ////rcutsq = square(cut[itype][jtype] * sigma);
  //r2inv = 1.0 / rsq;
  //rosigma2 = rsq/(sigma*sigma);
  //rosigma4 = rosigma2*rosigma2;
  //sigmaorn = powint(sqrt(1.0/rosigma2), powern);

  //forcelj = -1.0*(double)powern*sigmaorn 
  //        + 4.0*c4 * rosigma4
  //        + 2.0*c2 * rosigma2;


  fforce = factor_lj * forcelj * r2inv;

  philj = pow(1.0 - ros, alpha) / alpha;
  return factor_lj * philj;
}

/* ---------------------------------------------------------------------- */

void *PairSoftHarmonic::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str, "epsilon") == 0) return (void *) epsilon;
  return nullptr;
}
