/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 * Copyright (c) 2015 Peter W. Draper (p.w.draper@durham.ac.uk)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>

#include "error.h"
#include "part.h"
#include "cell.h"
#include "tools.h"
#include "swift.h"

/**
 *  Factorize a given integer, attempts to keep larger pair of factors.
 */
void factor(int value, int *f1, int *f2) {
  int j;
  int i;

  j = (int)sqrt(value);
  for (i = j; i > 0; i--) {
    if ((value % i) == 0) {
      *f1 = i;
      *f2 = value / i;
      break;
    }
  }
}

/**
 * @brief Compute the average number of pairs per particle using
 *      a brute-force O(N^2) computation.
 *
 * @param dim The space dimensions.
 * @param parts The #part array.
 * @param N The number of parts.
 * @param periodic Periodic boundary conditions flag.
 */

void pairs_n2(double *dim, struct part *__restrict__ parts, int N,
              int periodic) {

  int i, j, k, count = 0;
  // int mj, mk;
  // double maxratio = 1.0;
  double r2, dx[3], rho = 0.0;
  double rho_max = 0.0, rho_min = 100;

  /* Loop over all particle pairs. */
  for (j = 0; j < N; j++) {
    if (j % 1000 == 0) {
      printf("pairs_n2: j=%i.\n", j);
      fflush(stdout);
    }
    for (k = j + 1; k < N; k++) {
      for (i = 0; i < 3; i++) {
        dx[i] = parts[j].x[i] - parts[k].x[i];
        if (periodic) {
          if (dx[i] < -dim[i] / 2)
            dx[i] += dim[i];
          else if (dx[i] > dim[i] / 2)
            dx[i] -= dim[i];
        }
      }
      r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
      if (r2 < parts[j].h * parts[j].h || r2 < parts[k].h * parts[k].h) {
        runner_iact_density(r2, NULL, parts[j].h, parts[k].h, &parts[j],
                            &parts[k]);
        /* if ( parts[j].h / parts[k].h > maxratio )
            {
            maxratio = parts[j].h / parts[k].h;
            mj = j; mk = k;
            }
        else if ( parts[k].h / parts[j].h > maxratio )
            {
            maxratio = parts[k].h / parts[j].h;
            mj = j; mk = k;
            } */
      }
    }
  }

  /* Aggregate the results. */
  for (k = 0; k < N; k++) {
    // count += parts[k].icount;
    rho += parts[k].density.wcount;
    rho_min = fmin(parts[k].density.wcount, rho_min);
    rho_min = fmax(parts[k].density.wcount, rho_max);
  }

  /* Dump the result. */
  printf("pairs_n2: avg. density per part is %.3f (nr. pairs %.3f).\n",
         rho / N + 32.0 / 3, ((double)count) / N);
  printf("pairs_n2: densities are in [ %e , %e ].\n", rho_min / N + 32.0 / 3,
         rho_max / N + 32.0 / 3);
  /* printf( "pairs_n2: maximum ratio between parts %i [%e,%e,%e] and %i
     [%e,%e,%e] is %.3f/%.3f\n" ,
      mj , parts[mj].x[0] , parts[mj].x[1] , parts[mj].x[2] ,
      mk , parts[mk].x[0] , parts[mk].x[1] , parts[mk].x[2] ,
      parts[mj].h , parts[mk].h ); fflush(stdout); */
  fflush(stdout);
}

void pairs_single_density(double *dim, long long int pid,
                          struct part *__restrict__ parts, int N,
                          int periodic) {

  int i, k;
  // int mj, mk;
  // double maxratio = 1.0;
  double r2, dx[3];
  float fdx[3];
  struct part p;
  // double ih = 12.0/6.25;

  /* Find "our" part. */
  for (k = 0; k < N && parts[k].id != pid; k++)
    ;
  if (k == N) error("Part not found.");
  p = parts[k];
  printf("pairs_single: part[%i].id == %lli.\n", k, pid);

  p.rho = 0.0;
  p.density.wcount = 0.0;
  // p.icount = 0;
  p.rho_dh = 0.0;

  /* Loop over all particle pairs. */
  for (k = 0; k < N; k++) {
    if (parts[k].id == p.id) continue;
    for (i = 0; i < 3; i++) {
      dx[i] = p.x[i] - parts[k].x[i];
      if (periodic) {
        if (dx[i] < -dim[i] / 2)
          dx[i] += dim[i];
        else if (dx[i] > dim[i] / 2)
          dx[i] -= dim[i];
      }
      fdx[i] = dx[i];
    }
    r2 = fdx[0] * fdx[0] + fdx[1] * fdx[1] + fdx[2] * fdx[2];
    if (r2 < p.h * p.h) {
      runner_iact_nonsym_density(r2, fdx, p.h, parts[k].h, &p, &parts[k]);
      /* printf( "pairs_simple: interacting particles %lli [%i,%i,%i] and %lli
         [%i,%i,%i], r=%e.\n" ,
          pid , (int)(p.x[0]*ih) , (int)(p.x[1]*ih) , (int)(p.x[2]*ih) ,
          parts[k].id , (int)(parts[k].x[0]*ih) , (int)(parts[k].x[1]*ih) ,
         (int)(parts[k].x[2]*ih) ,
          sqrtf(r2) ); */
    }
  }

  /* Dump the result. */
  printf("pairs_single: wcount of part %lli (h=%e) is %f.\n", p.id, p.h,
         p.density.wcount + 32.0 / 3);
  fflush(stdout);
}

void pairs_all_density(struct runner *r, struct cell *ci, struct cell *cj) {

  float r2, hi, hj, hig2, hjg2, dx[3];
  struct part *pi, *pj;

  /* Implements a double-for loop and checks every interaction */
  for (int i = 0; i < ci->count; ++i) {

    pi = &ci->parts[i];
    hi = pi->h;
    hig2 = hi * hi * kernel_gamma2;

    for (int j = 0; j < cj->count; ++j) {

      pj = &cj->parts[j];

      /* Pairwise distance */
      r2 = 0.0f;
      for (int k = 0; k < 3; k++) {
        dx[k] = ci->parts[i].x[k] - cj->parts[j].x[k];
        r2 += dx[k] * dx[k];
      }

      /* Hit or miss? */
      if (r2 < hig2) {

        /* Interact */
        runner_iact_nonsym_density(r2, dx, hi, pj->h, pi, pj);
      }
    }
  }

  /* Reverse double-for loop and checks every interaction */
  for (int j = 0; j < cj->count; ++j) {

    pj = &cj->parts[j];
    hj = pj->h;
    hjg2 = hj * hj * kernel_gamma2;

    for (int i = 0; i < ci->count; ++i) {

      pi = &ci->parts[i];

      /* Pairwise distance */
      r2 = 0.0f;
      for (int k = 0; k < 3; k++) {
        dx[k] = cj->parts[j].x[k] - ci->parts[i].x[k];
        r2 += dx[k] * dx[k];
      }

      /* Hit or miss? */
      if (r2 < hjg2) {

        /* Interact */
        runner_iact_nonsym_density(r2, dx, hj, pi->h, pj, pi);
      }
    }
  }
}

void pairs_single_grav(double *dim, long long int pid,
                       struct gpart *__restrict__ parts, int N, int periodic) {

  int i, k;
  // int mj, mk;
  // double maxratio = 1.0;
  double r2, dx[3];
  float fdx[3], a[3] = {0.0, 0.0, 0.0}, aabs[3] = {0.0, 0.0, 0.0};
  struct gpart pi, pj;
  // double ih = 12.0/6.25;

  /* Find "our" part. */
  for (k = 0; k < N; k++)
    if ((parts[k].id > 0 && parts[k].part->id == pid) || parts[k].id == -pid)
      break;
  if (k == N) error("Part not found.");
  pi = parts[k];
  pi.a[0] = 0.0f;
  pi.a[1] = 0.0f;
  pi.a[2] = 0.0f;

  /* Loop over all particle pairs. */
  for (k = 0; k < N; k++) {
    if (parts[k].id == pi.id) continue;
    pj = parts[k];
    for (i = 0; i < 3; i++) {
      dx[i] = pi.x[i] - pj.x[i];
      if (periodic) {
        if (dx[i] < -dim[i] / 2)
          dx[i] += dim[i];
        else if (dx[i] > dim[i] / 2)
          dx[i] -= dim[i];
      }
      fdx[i] = dx[i];
    }
    r2 = fdx[0] * fdx[0] + fdx[1] * fdx[1] + fdx[2] * fdx[2];
    runner_iact_grav(r2, fdx, &pi, &pj);
    a[0] += pi.a[0];
    a[1] += pi.a[1];
    a[2] += pi.a[2];
    aabs[0] += fabsf(pi.a[0]);
    aabs[1] += fabsf(pi.a[1]);
    aabs[2] += fabsf(pi.a[2]);
    pi.a[0] = 0.0f;
    pi.a[1] = 0.0f;
    pi.a[2] = 0.0f;
  }

  /* Dump the result. */
  message(
      "acceleration on gpart %lli is a=[ %e %e %e ], |a|=[ %.2e %.2e %.2e ].\n",
      pi.part->id, a[0], a[1], a[2], aabs[0], aabs[1], aabs[2]);
}

/**
 * @brief Test the density function by dumping it for two random parts.
 *
 * @param N number of intervals in [0,1].
 */

void density_dump(int N) {

  int k;
  float r2[4] = {0.0f, 0.0f, 0.0f, 0.0f}, hi[4], hj[4];
  struct part *pi[4], *pj[4], Pi[4], Pj[4];

  /* Init the interaction parameters. */
  for (k = 0; k < 4; k++) {
    Pi[k].mass = 1.0f;
    Pi[k].rho = 0.0f;
    Pi[k].density.wcount = 0.0f;
    Pj[k].mass = 1.0f;
    Pj[k].rho = 0.0f;
    Pj[k].density.wcount = 0.0f;
    hi[k] = 1.0;
    hj[k] = 1.0;
    pi[k] = &Pi[k];
    pj[k] = &Pj[k];
  }

  for (k = 0; k <= N; k++) {
    r2[3] = r2[2];
    r2[2] = r2[1];
    r2[1] = r2[0];
    r2[0] = ((float)k) / N;
    Pi[0].density.wcount = 0;
    Pj[0].density.wcount = 0;
    runner_iact_density(r2[0], NULL, hi[0], hj[0], &Pi[0], &Pj[0]);
    printf(" %e %e %e", r2[0], Pi[0].density.wcount, Pj[0].density.wcount);
    Pi[0].density.wcount = 0;
    Pj[0].density.wcount = 0;
    Pi[1].density.wcount = 0;
    Pj[1].density.wcount = 0;
    Pi[2].density.wcount = 0;
    Pj[2].density.wcount = 0;
    Pi[3].density.wcount = 0;
    Pj[3].density.wcount = 0;
    runner_iact_vec_density(r2, NULL, hi, hj, pi, pj);
    printf(" %e %e %e %e\n", Pi[0].density.wcount, Pi[1].density.wcount,
           Pi[2].density.wcount, Pi[3].density.wcount);
  }
}
