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
   Contributing author: Paul Crozier (SNL)
------------------------------------------------------------------------- */

#include "pair_coresoftened.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_special.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
//#include "respa.h"
//#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using MathSpecial::powint;

/* ---------------------------------------------------------------------- */

PairCoresoftenedCut::PairCoresoftenedCut(LAMMPS *lmp) : Pair(lmp)
{
  //respa_enable = 1;
  //born_matrix_enable = 1;
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairCoresoftenedCut::~PairCoresoftenedCut()
{
  if (copymode) return;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(epsilon);
    memory->destroy(sigma);
  }
}

/* ---------------------------------------------------------------------- */

void PairCoresoftenedCut::compute(int eflag, int vflag)
{
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz, evdwl, fpair;
  double rsq, r2inv, r6inv, r14inv, rk, thrks, forcelj, factor_lj;
  int *ilist, *jlist, *numneigh, **firstneigh;

  evdwl = 0.0;
  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
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

      // V(r) = (sigma/r)^14 + 0.5*(1-tanh(k(r-sigma1)))
      if (rsq < cutsq[itype][jtype]) {
        r2inv  = 1.0 / rsq;
        r6inv  = r2inv*r2inv*r2inv;
        r14inv = r6inv*r6inv*r2inv;///////////////////////
        //r14inv = r2inv*r2inv*r2inv*r2inv*r2inv*r2inv*r2inv;///////////////////////
        //r14inv = powint(r2inv,7);
        rk     = sqrt(rsq) * indk;
        thrks  = tanh(rk - ks1); 

        forcelj = 14.0 * r14inv + rk / 2.0 * (1.0 - thrks * thrks);
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
          evdwl = r14inv + 0.5 * (1.0 - thrks);
          evdwl *= factor_lj;
        }


        if (evflag) ev_tally(i, j, nlocal, newton_pair, evdwl, 0.0, fpair, delx, dely, delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}


/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairCoresoftenedCut::allocate()
{
  allocated = 1;
  int n = atom->ntypes + 1;

  memory->create(setflag, n, n, "pair:setflag");
  for (int i = 1; i < n; i++)
    for (int j = i; j < n; j++) setflag[i][j] = 0;

  memory->create(cutsq, n, n, "pair:cutsq");

  memory->create(cut, n, n, "pair:cut");
  memory->create(epsilon, n, n, "pair:epsilon");
  memory->create(sigma, n, n, "pair:sigma");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairCoresoftenedCut::settings(int narg, char **arg)
{
  if (narg != 3) error->all(FLERR, "Illegal pair_style command: cut/sigma1/k");

  cut_global = utils::numeric(FLERR, arg[0], false, lmp);
  sigma1     = utils::numeric(FLERR, arg[1], false, lmp);
  indk       = utils::numeric(FLERR, arg[2], false, lmp);

  ks1        = indk * sigma1;

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

void PairCoresoftenedCut::coeff(int narg, char **arg)
{
  if (narg < 4 || narg > 5) error->all(FLERR, "Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  double epsilon_one = utils::numeric(FLERR, arg[2], false, lmp);
  double sigma_one = utils::numeric(FLERR, arg[3], false, lmp);

  double cut_one = cut_global;
  if (narg == 5) cut_one = utils::numeric(FLERR, arg[4], false, lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      epsilon[i][j] = epsilon_one;
      sigma[i][j] = sigma_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

///* ----------------------------------------------------------------------
//   init specific to this pair style
//------------------------------------------------------------------------- */
//
//void PairCoresoftenedCut::init_style()
//{
//  Pair::init_style();
//
//  //int list_style = NeighConst::REQ_DEFAULT;//////////////////////////
//  //// request regular or rRESPA neighbor list
//
//  //int list_style = NeighConst::REQ_DEFAULT;
//
//  //if (update->whichflag == 1 && utils::strmatch(update->integrate_style, "^respa")) {
//  //  auto respa = dynamic_cast<Respa *>(update->integrate);
//  //  if (respa->level_inner >= 0) list_style = NeighConst::REQ_RESPA_INOUT;
//  //  if (respa->level_middle >= 0) list_style = NeighConst::REQ_RESPA_ALL;
//  //}
//  //neighbor->add_request(this, list_style);
//
//  //// set rRESPA cutoffs
//
//  //if (utils::strmatch(update->integrate_style, "^respa") &&
//  //    (dynamic_cast<Respa *>(update->integrate))->level_inner >= 0)
//  //  cut_respa = (dynamic_cast<Respa *>(update->integrate))->cutoff;
//  //else
//  //  cut_respa = nullptr;
//}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairCoresoftenedCut::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    epsilon[i][j] = mix_energy(epsilon[i][i], epsilon[j][j], sigma[i][i], sigma[j][j]);
    sigma[i][j] = mix_distance(sigma[i][i], sigma[j][j]);
    cut[i][j] = mix_distance(cut[i][i], cut[j][j]);
  }

  epsilon[j][i] = epsilon[i][j];
  sigma[j][i] = sigma[i][j];
  cut[j][i] = cut[i][j];

  //lj1[i][j] = 48.0 * epsilon[i][j] * pow(sigma[i][j], 12.0);
  //lj2[i][j] = 24.0 * epsilon[i][j] * pow(sigma[i][j], 6.0);
  //lj3[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j], 12.0);
  //lj4[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j], 6.0);

  //if (offset_flag && (cut[i][j] > 0.0)) {
  //  double ratio = sigma[i][j] / cut[i][j];
  //  offset[i][j] = 4.0 * epsilon[i][j] * (pow(ratio, 12.0) - pow(ratio, 6.0));
  //} else
  //  offset[i][j] = 0.0;

  //lj1[j][i] = lj1[i][j];
  //lj2[j][i] = lj2[i][j];
  //lj3[j][i] = lj3[i][j];
  //lj4[j][i] = lj4[i][j];
  //offset[j][i] = offset[i][j];

  //// check interior rRESPA cutoff

  //if (cut_respa && cut[i][j] < cut_respa[3])
  //  error->all(FLERR, "Pair cutoff < Respa interior cutoff");

  //// compute I,J contribution to long-range tail correction
  //// count total # of atoms of type I and J via Allreduce

  //if (tail_flag) {
  //  int *type = atom->type;
  //  int nlocal = atom->nlocal;

  //  double count[2], all[2];
  //  count[0] = count[1] = 0.0;
  //  for (int k = 0; k < nlocal; k++) {
  //    if (type[k] == i) count[0] += 1.0;
  //    if (type[k] == j) count[1] += 1.0;
  //  }
  //  MPI_Allreduce(count, all, 2, MPI_DOUBLE, MPI_SUM, world);

  //  double sig2 = sigma[i][j] * sigma[i][j];
  //  double sig6 = sig2 * sig2 * sig2;
  //  double rc3 = cut[i][j] * cut[i][j] * cut[i][j];
  //  double rc6 = rc3 * rc3;
  //  double rc9 = rc3 * rc6;
  //  double prefactor = 8.0 * MY_PI * all[0] * all[1] * epsilon[i][j] * sig6 / (9.0 * rc9);
  //  etail_ij = prefactor * (sig6 - 3.0 * rc6);
  //  ptail_ij = 2.0 * prefactor * (2.0 * sig6 - 3.0 * rc6);
  //}

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairCoresoftenedCut::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i, j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j], sizeof(int), 1, fp);
      if (setflag[i][j]) {
        fwrite(&epsilon[i][j], sizeof(double), 1, fp);
        fwrite(&sigma[i][j], sizeof(double), 1, fp);
        fwrite(&cut[i][j], sizeof(double), 1, fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairCoresoftenedCut::read_restart(FILE *fp)
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
          utils::sfread(FLERR, &sigma[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &cut[i][j], sizeof(double), 1, fp, nullptr, error);
        }
        MPI_Bcast(&epsilon[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&sigma[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&cut[i][j], 1, MPI_DOUBLE, 0, world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairCoresoftenedCut::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global, sizeof(double), 1, fp);
  fwrite(&sigma1, sizeof(double), 1, fp);
  fwrite(&indk, sizeof(double), 1, fp);
  //fwrite(&offset_flag, sizeof(int), 1, fp);
  fwrite(&mix_flag, sizeof(int), 1, fp);
  //fwrite(&tail_flag, sizeof(int), 1, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairCoresoftenedCut::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    utils::sfread(FLERR, &cut_global, sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &sigma1, sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &indk, sizeof(double), 1, fp, nullptr, error);
    //utils::sfread(FLERR, &offset_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &mix_flag, sizeof(int), 1, fp, nullptr, error);
    //utils::sfread(FLERR, &tail_flag, sizeof(int), 1, fp, nullptr, error);

    ks1 = indk * sigma1;
  }
  MPI_Bcast(&cut_global, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&sigma1, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&indk, 1, MPI_DOUBLE, 0, world);
  //MPI_Bcast(&offset_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&mix_flag, 1, MPI_INT, 0, world);
  //MPI_Bcast(&tail_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&ks1, 1, MPI_DOUBLE, 0, world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairCoresoftenedCut::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++) fprintf(fp, "%d %g %g\n", i, epsilon[i][i], sigma[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairCoresoftenedCut::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp, "%d %d %g %g %g\n", i, j, epsilon[i][j], sigma[i][j], cut[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairCoresoftenedCut::single(int /*i*/, int /*j*/, int itype, int jtype, double rsq,
                         double /*factor_coul*/, double factor_lj, double &fforce)
{
  double r2inv, r6inv, r14inv, rk, thrks, forcelj, philj;

  r2inv  = 1.0/rsq;
  r6inv  = r2inv*r2inv*r2inv;
  r14inv = r6inv*r6inv*r2inv;

  //r14inv = powint(r2inv,7);
  //r      = sqrt(rsq);
  //rk     = r * indk;
  rk       = sqrt(rsq) * indk;
  thrks  = tanh(rk - ks1); 

  forcelj = 14.0 * r14inv + rk / 2.0 * (1.0 - thrks * thrks);
  //fforce = forcelj * r2inv;//////////////////////////////////
  fforce = factor_lj * forcelj * r2inv;

  philj = r14inv + 0.5 * (1.0 - thrks);
  return factor_lj * philj;
}

/* ---------------------------------------------------------------------- */

void *PairCoresoftenedCut::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str, "epsilon") == 0) return (void *) epsilon;
  if (strcmp(str, "sigma") == 0) return (void *) sigma;
  return nullptr;
}