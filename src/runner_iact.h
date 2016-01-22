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

  /* old density calculation */
  float r = sqrtf(r2), ri = 1.0f / r;
  float xi, xj;
  float h_inv;
  float wi, wj, wi_dx, wj_dx;
  float mi, mj;
  float dvdr;
  float dv[3], curlvr[3];
  int k;

  /* Get the masses. */
  mi = pi->mass;
  mj = pj->mass;

  /* Compute dv dot r */
  dv[0] = pi->v[0] - pj->v[0];
  dv[1] = pi->v[1] - pj->v[1];
  dv[2] = pi->v[2] - pj->v[2];
  dvdr = dv[0] * dx[0] + dv[1] * dx[1] + dv[2] * dx[2];
  dvdr *= ri;

  /* Compute dv cross r */
  curlvr[0] = dv[1] * dx[2] - dv[2] * dx[1];
  curlvr[1] = dv[2] * dx[0] - dv[0] * dx[2];
  curlvr[2] = dv[0] * dx[1] - dv[1] * dx[0];
  for (k = 0; k < 3; k++) curlvr[k] *= ri;

  /* Compute density of pi. */
  h_inv = 1.0 / hi;
  xi = r * h_inv;
  kernel_deval(xi, &wi, &wi_dx);

  pi->rho += mj * wi;
  pi->rho_dh -= mj * (3.0 * wi + xi * wi_dx);
  pi->density.wcount += wi;
  pi->density.wcount_dh -= xi * wi_dx;

  pi->density.div_v += mj * dvdr * wi_dx;
  for (k = 0; k < 3; k++) pi->density.curl_v[k] += mj * curlvr[k] * wi_dx;

  /* Compute density of pj. */
  h_inv = 1.0 / hj;
  xj = r * h_inv;
  kernel_deval(xj, &wj, &wj_dx);

  pj->rho += mi * wj;
  pj->rho_dh -= mi * (3.0 * wj + xj * wj_dx);
  pj->density.wcount += wj;
  pj->density.wcount_dh -= xj * wj_dx;

  pj->density.div_v += mi * dvdr * wj_dx;
  for (k = 0; k < 3; k++) pj->density.curl_v[k] += mi * curlvr[k] * wj_dx;
  
  /* voronoi calculation, added as test */
  voronoi_intersect(r2, dx, pi, pj);
  
  float xd[3];
  xd[0] = -dx[0];
  xd[1] = -dx[1];
  xd[2] = -dx[2];
  voronoi_intersect(r2, xd, pj, pi);

}

/**
 * @brief Density loop (non-symmetric version)
 */

__attribute__((always_inline)) INLINE static void runner_iact_nonsym_density(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

  /* old density calculation */
  float r, ri;
  float xi;
  float h_inv;
  float wi, wi_dx;
  float mj;
  float dvdr;
  float dv[3], curlvr[3];
  int k;

  /* Get the masses. */
  mj = pj->mass;

  /* Get r and r inverse. */
  r = sqrtf(r2);
  ri = 1.0f / r;

  /* Compute dv dot r */
  dv[0] = pi->v[0] - pj->v[0];
  dv[1] = pi->v[1] - pj->v[1];
  dv[2] = pi->v[2] - pj->v[2];
  dvdr = dv[0] * dx[0] + dv[1] * dx[1] + dv[2] * dx[2];
  dvdr *= ri;

  /* Compute dv cross r */
  curlvr[0] = dv[1] * dx[2] - dv[2] * dx[1];
  curlvr[1] = dv[2] * dx[0] - dv[0] * dx[2];
  curlvr[2] = dv[0] * dx[1] - dv[1] * dx[0];
  for (k = 0; k < 3; k++) curlvr[k] *= ri;

  h_inv = 1.0 / hi;
  xi = r * h_inv;
  kernel_deval(xi, &wi, &wi_dx);

  pi->rho += mj * wi;
  pi->rho_dh -= mj * (3.0 * wi + xi * wi_dx);
  pi->density.wcount += wi;
  pi->density.wcount_dh -= xi * wi_dx;

  pi->density.div_v += mj * dvdr * wi_dx;
  for (k = 0; k < 3; k++) pi->density.curl_v[k] += mj * curlvr[k] * wi_dx;

  /* test voronoi calculation */
  voronoi_intersect(r2, dx, pi, pj);
}

__attribute__((always_inline)) INLINE static void runner_iact_force_common(
    float r2, float *dx, struct part *pi, struct part *pj) {

  int idx = voronoi_get_face_index(pi, pj);
  if(idx >= 0){
    float A = pi->voronoi.face_areas[idx];
    float WL[5], WR[5], Whalf[5];
    float n_unit[3];
    float r = sqrtf(r2);
//    float vproj;
    float ehalf;
    float v2;
    float flux[5][3];

    WL[0] = pi->primitives.rho;
    WL[1] = pi->primitives.v[0];
    WL[2] = pi->primitives.v[1];
    WL[3] = pi->primitives.v[2];
    WL[4] = pi->primitives.P;

    WR[0] = pj->primitives.rho;
    WR[1] = pj->primitives.v[0];
    WR[2] = pj->primitives.v[1];
    WR[3] = pj->primitives.v[2];
    WR[4] = pj->primitives.P;

    n_unit[0] = dx[0]/r;
    n_unit[1] = dx[1]/r;
    n_unit[2] = dx[2]/r;

    riemann_solver_solve(WL, WR, Whalf, n_unit);

//    vproj = Whalf[1]*n_unit[0] + Whalf[2]*n_unit[1] + Whalf[3]*n_unit[2];
    v2 = Whalf[1] * Whalf[1] + Whalf[2] * Whalf[2] + Whalf[3] * Whalf[3];
    ehalf = Whalf[4] / (const_hydro_gamma - 1.0f) / Whalf[0] + 0.5 * v2;

    flux[0][0] = Whalf[0] * Whalf[1];
    flux[0][1] = Whalf[0] * Whalf[2];
    flux[0][2] = Whalf[0] * Whalf[3];

    flux[1][0] = Whalf[0] * Whalf[1] * Whalf[1] + Whalf[4];
    flux[1][1] = Whalf[0] * Whalf[1] * Whalf[2];
    flux[1][2] = Whalf[0] * Whalf[1] * Whalf[3];

    flux[2][0] = Whalf[0] * Whalf[2] * Whalf[1];
    flux[2][1] = Whalf[0] * Whalf[2] * Whalf[2] + Whalf[4];
    flux[2][2] = Whalf[0] * Whalf[2] * Whalf[3];

    flux[3][0] = Whalf[0] * Whalf[3] * Whalf[1];
    flux[3][1] = Whalf[0] * Whalf[3] * Whalf[2];
    flux[3][2] = Whalf[0] * Whalf[3] * Whalf[3] + Whalf[4];

    flux[4][0] = Whalf[0] * ehalf * Whalf[1] + Whalf[4] * Whalf[1];
    flux[4][1] = Whalf[0] * ehalf * Whalf[2] + Whalf[4] * Whalf[2];
    flux[4][2] = Whalf[0] * ehalf * Whalf[3] + Whalf[4] * Whalf[3];

//    pi->conserved.dm -= (Whalf[0] * vproj) * A;
//    pi->conserved.dp[0] -= (Whalf[0] * Whalf[1] * vproj + Whalf[4]*n_unit[0]) * A;
//    pi->conserved.dp[1] -= (Whalf[0] * Whalf[2] * vproj + Whalf[4]*n_unit[1]) * A;
//    pi->conserved.dp[2] -= (Whalf[0] * Whalf[3] * vproj + Whalf[4]*n_unit[2]) * A;
//    pi->conserved.de -= (Whalf[0] * ehalf * vproj + Whalf[4] * vproj) * A;

//    pj->conserved.dm += (Whalf[0] * vproj) * A;
//    pj->conserved.dp[0] += (Whalf[0] * Whalf[1] * vproj + Whalf[4]*n_unit[0]) * A;
//    pj->conserved.dp[1] += (Whalf[0] * Whalf[2] * vproj + Whalf[4]*n_unit[1]) * A;
//    pj->conserved.dp[2] += (Whalf[0] * Whalf[3] * vproj + Whalf[4]*n_unit[2]) * A;
//    pj->conserved.de += (Whalf[0] * ehalf * vproj + Whalf[4] * vproj) * A;
    pi->conserved.dm -= A * (flux[0][0] * n_unit[0] + flux[0][1] * n_unit[1] + flux[0][2] * n_unit[2]);
    pi->conserved.dp[0] -= A * (flux[1][0] * n_unit[0] + flux[1][1] * n_unit[1] + flux[1][2] * n_unit[2]);
    pi->conserved.dp[1] -= A * (flux[2][0] * n_unit[0] + flux[2][1] * n_unit[1] + flux[2][2] * n_unit[2]);
    pi->conserved.dp[2] -= A * (flux[3][0] * n_unit[0] + flux[3][1] * n_unit[1] + flux[3][2] * n_unit[2]);
    pi->conserved.de -= A * (flux[4][0] * n_unit[0] + flux[4][1] * n_unit[1] + flux[4][2] * n_unit[2]);

    pj->conserved.dm += A * (flux[0][0] * n_unit[0] + flux[0][1] * n_unit[1] + flux[0][2] * n_unit[2]);
    pj->conserved.dp[0] += A * (flux[1][0] * n_unit[0] + flux[1][1] * n_unit[1] + flux[1][2] * n_unit[2]);
    pj->conserved.dp[1] += A * (flux[2][0] * n_unit[0] + flux[2][1] * n_unit[1] + flux[2][2] * n_unit[2]);
    pj->conserved.dp[2] += A * (flux[3][0] * n_unit[0] + flux[3][1] * n_unit[1] + flux[3][2] * n_unit[2]);
    pj->conserved.de += A * (flux[4][0] * n_unit[0] + flux[4][1] * n_unit[1] + flux[4][2] * n_unit[2]);
  }

}

/**
 * @brief Force loop
 */

__attribute__((always_inline)) INLINE static void runner_iact_force(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

  float r = sqrtf(r2), ri = 1.0f / r;
  float xi, xj;
  float hi_inv, hi2_inv;
  float hj_inv, hj2_inv;
  float wi, wj, wi_dx, wj_dx, wi_dr, wj_dr, w, dvdr;
  float mi, mj, POrho2i, POrho2j, rhoi, rhoj;
  float v_sig, omega_ij, Pi_ij, alpha_ij, tc, v_sig_u;
  // float dt_max;
  float f;
  int k;

  /* Get some values in local variables. */
  mi = pi->mass;
  mj = pj->mass;
  rhoi = pi->rho;
  rhoj = pj->rho;
  POrho2i = pi->force.POrho2;
  POrho2j = pj->force.POrho2;

  /* Get the kernel for hi. */
  hi_inv = 1.0f / hi;
  hi2_inv = hi_inv * hi_inv;
  xi = r * hi_inv;
  kernel_deval(xi, &wi, &wi_dx);
  wi_dr = hi2_inv * hi2_inv * wi_dx;

  /* Get the kernel for hj. */
  hj_inv = 1.0f / hj;
  hj2_inv = hj_inv * hj_inv;
  xj = r * hj_inv;
  kernel_deval(xj, &wj, &wj_dx);
  wj_dr = hj2_inv * hj2_inv * wj_dx;

  /* Compute dv dot r. */
  dvdr = (pi->v[0] - pj->v[0]) * dx[0] + (pi->v[1] - pj->v[1]) * dx[1] +
         (pi->v[2] - pj->v[2]) * dx[2];
  dvdr *= ri;

  /* Compute the relative velocity. (This is 0 if the particles move away from
   * each other and negative otherwise) */
  omega_ij = fminf(dvdr, 0.f);

  /* Compute signal velocity */
  v_sig = pi->force.c + pj->force.c - 2.0f * omega_ij;

  /* Compute viscosity parameter */
  alpha_ij = -0.5f * (pi->alpha + pj->alpha);

  /* Compute viscosity tensor */
  Pi_ij = alpha_ij * v_sig * omega_ij / (rhoi + rhoj);

  /* Apply balsara switch */
  Pi_ij *= (pi->force.balsara + pj->force.balsara);

  /* Thermal conductivity */
  v_sig_u = sqrtf(2.f * (const_hydro_gamma - 1.f) *
                  fabs(rhoi * pi->u - rhoj * pj->u) / (rhoi + rhoj));
  tc = const_conductivity_alpha * v_sig_u / (rhoi + rhoj);
  tc *= (wi_dr + wj_dr);

  /* Get the common factor out. */
  w = ri *
      ((POrho2i * wi_dr + POrho2j * wj_dr) + 0.25f * Pi_ij * (wi_dr + wj_dr));

  /* Use the force, Luke! */
  for (k = 0; k < 3; k++) {
    f = dx[k] * w;
    pi->a[k] -= mj * f;
    pj->a[k] += mi * f;
  }

  /* Get the time derivative for u. */
  pi->force.u_dt +=
      mj * dvdr * (POrho2i * wi_dr + 0.125f * Pi_ij * (wi_dr + wj_dr));
  pj->force.u_dt +=
      mi * dvdr * (POrho2j * wj_dr + 0.125f * Pi_ij * (wi_dr + wj_dr));

  /* Add the thermal conductivity */
  pi->force.u_dt += mj * tc * (pi->u - pj->u);
  pj->force.u_dt += mi * tc * (pj->u - pi->u);

  /* Get the time derivative for h. */
  pi->force.h_dt -= mj * dvdr / rhoj * wi_dr;
  pj->force.h_dt -= mi * dvdr / rhoi * wj_dr;

  /* Update the signal velocity. */
  pi->force.v_sig = fmaxf(pi->force.v_sig, v_sig);
  pj->force.v_sig = fmaxf(pj->force.v_sig, v_sig);

  /* Moving mesh stuff */
  runner_iact_force_common(r2, dx, pi, pj);
}

/**
 * @brief Force loop (non-symmetric version)
 */

__attribute__((always_inline)) INLINE static void runner_iact_nonsym_force(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

  float r = sqrtf(r2), ri = 1.0f / r;
  float xi, xj;
  float hi_inv, hi2_inv;
  float hj_inv, hj2_inv;
  float wi, wj, wi_dx, wj_dx, wi_dr, wj_dr, w, dvdr;
  float /*mi,*/ mj, POrho2i, POrho2j, rhoi, rhoj;
  float v_sig, omega_ij, Pi_ij, alpha_ij, tc, v_sig_u;
  // float dt_max;
  float f;
  int k;

  /* Get some values in local variables. */
  // mi = pi->mass;
  mj = pj->mass;
  rhoi = pi->rho;
  rhoj = pj->rho;
  POrho2i = pi->force.POrho2;
  POrho2j = pj->force.POrho2;

  /* Get the kernel for hi. */
  hi_inv = 1.0f / hi;
  hi2_inv = hi_inv * hi_inv;
  xi = r * hi_inv;
  kernel_deval(xi, &wi, &wi_dx);
  wi_dr = hi2_inv * hi2_inv * wi_dx;

  /* Get the kernel for hj. */
  hj_inv = 1.0f / hj;
  hj2_inv = hj_inv * hj_inv;
  xj = r * hj_inv;
  kernel_deval(xj, &wj, &wj_dx);
  wj_dr = hj2_inv * hj2_inv * wj_dx;

  /* Compute dv dot r. */
  dvdr = (pi->v[0] - pj->v[0]) * dx[0] + (pi->v[1] - pj->v[1]) * dx[1] +
         (pi->v[2] - pj->v[2]) * dx[2];
  dvdr *= ri;

  /* Compute the relative velocity. (This is 0 if the particles move away from
   * each other and negative otherwise) */
  omega_ij = fminf(dvdr, 0.f);

  /* Compute signal velocity */
  v_sig = pi->force.c + pj->force.c - 2.0f * omega_ij;

  /* Compute viscosity parameter */
  alpha_ij = -0.5f * (pi->alpha + pj->alpha);

  /* Compute viscosity tensor */
  Pi_ij = alpha_ij * v_sig * omega_ij / (rhoi + rhoj);

  /* Apply balsara switch */
  Pi_ij *= (pi->force.balsara + pj->force.balsara);

  /* Thermal conductivity */
  v_sig_u = sqrtf(2.f * (const_hydro_gamma - 1.f) *
                  fabs(rhoi * pi->u - rhoj * pj->u) / (rhoi + rhoj));
  tc = const_conductivity_alpha * v_sig_u / (rhoi + rhoj);
  tc *= (wi_dr + wj_dr);

  /* Get the common factor out. */
  w = ri *
      ((POrho2i * wi_dr + POrho2j * wj_dr) + 0.25f * Pi_ij * (wi_dr + wj_dr));

  /* Use the force, Luke! */
  for (k = 0; k < 3; k++) {
    f = dx[k] * w;
    pi->a[k] -= mj * f;
  }

  /* Get the time derivative for u. */
  pi->force.u_dt +=
      mj * dvdr * (POrho2i * wi_dr + 0.125f * Pi_ij * (wi_dr + wj_dr));

  /* Add the thermal conductivity */
  pi->force.u_dt += mj * tc * (pi->u - pj->u);

  /* Get the time derivative for h. */
  pi->force.h_dt -= mj * dvdr / rhoj * wi_dr;

  /* Update the signal velocity. */
  pi->force.v_sig = fmaxf(pi->force.v_sig, v_sig);
  pj->force.v_sig = fmaxf(pj->force.v_sig, v_sig);

  /* Moving mesh stuff */
  runner_iact_force_common(r2, dx, pi, pj);
}

#endif /* SWIFT_RUNNER_IACT_H */
