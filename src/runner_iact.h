/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
#ifndef SWIFT_RUNNER_IACT_H
#define SWIFT_RUNNER_IACT_H

/* Includes. */
#include "const.h"
#include "kernel.h"
#include "part.h"
#include "riemann.h"
#include "voronoi.h"

/**
 * @file  runner_iact.h
 * @brief SPH interaction functions following the Gadget-2 version of SPH.
 *
 * The interactions computed here are the ones presented in the Gadget-2 paper
 *and use the same
 * numerical coefficients as the Gadget-2 code. When used with the Spline-3
 *kernel, the results
 * should be equivalent to the ones obtained with Gadget-2 up to the rounding
 *errors and interactions
 * missed by the Gadget-2 tree-code neighbours search.
 *
 * The code uses internal energy instead of entropy as a thermodynamical
 *variable.
 */

/**
 * @brief Density loop
 */

__attribute__((always_inline)) INLINE static void runner_iact_density(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

  float xd[3];

  /* voronoi calculation, added as test */
  voronoi_intersect(dx, pi, pj);

  xd[0] = -dx[0];
  xd[1] = -dx[1];
  xd[2] = -dx[2];
  voronoi_intersect(xd, pj, pi);

}

/**
 * @brief Density loop (non-symmetric version)
 */

__attribute__((always_inline)) INLINE static void runner_iact_nonsym_density(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

  /* test voronoi calculation */
  voronoi_intersect(dx, pi, pj);

}

__attribute__((always_inline)) INLINE static void runner_iact_force_common(
    float r2, float *dx, struct part *pi, struct part *pj) {

  int idx;

  idx = voronoi_get_face_index(pi, pj);
  if(idx >= 0){

    float A, midface[3], vface[3];
    float WL[5], WR[5], Whalf[5];
    float n_unit[3], r, ehalf, v2, flux[5][3];

    A = pi->voronoi.face_areas[idx];
    midface[0] = pi->voronoi.face_midpoints[3*idx];
    midface[1] = pi->voronoi.face_midpoints[3*idx+1];
    midface[2] = pi->voronoi.face_midpoints[3*idx+2];

    voronoi_get_face_velocity(r2, dx, pi, pj, midface, vface);
    r = sqrtf(r2);

    WL[0] = pi->primitives.rho;
    WL[1] = pi->primitives.v[0] - vface[0];
    WL[2] = pi->primitives.v[1] - vface[1];
    WL[3] = pi->primitives.v[2] - vface[2];
    WL[4] = pi->primitives.P;

    WR[0] = pj->primitives.rho;
    WR[1] = pj->primitives.v[0] - vface[0];
    WR[2] = pj->primitives.v[1] - vface[1];
    WR[3] = pj->primitives.v[2] - vface[2];
    WR[4] = pj->primitives.P;

    /* dx = pi->x - pj->x */
    n_unit[0] = -dx[0]/r;
    n_unit[1] = -dx[1]/r;
    n_unit[2] = -dx[2]/r;

    riemann_solver_solve(WL, WR, Whalf, n_unit);

    Whalf[1] += vface[0];
    Whalf[2] += vface[1];
    Whalf[3] += vface[2];

    v2 = Whalf[1] * Whalf[1] + Whalf[2] * Whalf[2] + Whalf[3] * Whalf[3];
    ehalf = Whalf[4] / (const_hydro_gamma - 1.0f) / Whalf[0] + 0.5 * v2;

    flux[0][0] = Whalf[0] * (Whalf[1] - vface[0]);
    flux[0][1] = Whalf[0] * (Whalf[2] - vface[1]);
    flux[0][2] = Whalf[0] * (Whalf[3] - vface[2]);

    flux[1][0] = Whalf[0] * Whalf[1] * (Whalf[1] - vface[0]) + Whalf[4];
    flux[1][1] = Whalf[0] * Whalf[1] * (Whalf[2] - vface[1]);
    flux[1][2] = Whalf[0] * Whalf[1] * (Whalf[3] - vface[2]);

    flux[2][0] = Whalf[0] * Whalf[2] * (Whalf[1] - vface[0]);
    flux[2][1] = Whalf[0] * Whalf[2] * (Whalf[2] - vface[1]) + Whalf[4];
    flux[2][2] = Whalf[0] * Whalf[2] * (Whalf[3] - vface[2]);

    flux[3][0] = Whalf[0] * Whalf[3] * (Whalf[1] - vface[0]);
    flux[3][1] = Whalf[0] * Whalf[3] * (Whalf[2] - vface[1]);
    flux[3][2] = Whalf[0] * Whalf[3] * (Whalf[3] - vface[2]) + Whalf[4];

    flux[4][0] = Whalf[0] * ehalf * (Whalf[1] - vface[0]) + Whalf[4] * Whalf[1];
    flux[4][1] = Whalf[0] * ehalf * (Whalf[2] - vface[1]) + Whalf[4] * Whalf[2];
    flux[4][2] = Whalf[0] * ehalf * (Whalf[3] - vface[2]) + Whalf[4] * Whalf[3];

    pi->conserved.dm += A * (flux[0][0] * n_unit[0] + flux[0][1] * n_unit[1] + flux[0][2] * n_unit[2]);
    pi->conserved.dp[0] += A * (flux[1][0] * n_unit[0] + flux[1][1] * n_unit[1] + flux[1][2] * n_unit[2]);
    pi->conserved.dp[1] += A * (flux[2][0] * n_unit[0] + flux[2][1] * n_unit[1] + flux[2][2] * n_unit[2]);
    pi->conserved.dp[2] += A * (flux[3][0] * n_unit[0] + flux[3][1] * n_unit[1] + flux[3][2] * n_unit[2]);
    pi->conserved.de += A * (flux[4][0] * n_unit[0] + flux[4][1] * n_unit[1] + flux[4][2] * n_unit[2]);

    pj->conserved.dm -= A * (flux[0][0] * n_unit[0] + flux[0][1] * n_unit[1] + flux[0][2] * n_unit[2]);
    pj->conserved.dp[0] -= A * (flux[1][0] * n_unit[0] + flux[1][1] * n_unit[1] + flux[1][2] * n_unit[2]);
    pj->conserved.dp[1] -= A * (flux[2][0] * n_unit[0] + flux[2][1] * n_unit[1] + flux[2][2] * n_unit[2]);
    pj->conserved.dp[2] -= A * (flux[3][0] * n_unit[0] + flux[3][1] * n_unit[1] + flux[3][2] * n_unit[2]);
    pj->conserved.de -= A * (flux[4][0] * n_unit[0] + flux[4][1] * n_unit[1] + flux[4][2] * n_unit[2]);

  }

}

/**
 * @brief Force loop
 */

__attribute__((always_inline)) INLINE static void runner_iact_force(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

  /* Moving mesh stuff */
  runner_iact_force_common(r2, dx, pi, pj);

}

/**
 * @brief Force loop (non-symmetric version)
 */

__attribute__((always_inline)) INLINE static void runner_iact_nonsym_force(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

  /* Moving mesh stuff */
  runner_iact_force_common(r2, dx, pi, pj);

}

#endif /* SWIFT_RUNNER_IACT_H */
