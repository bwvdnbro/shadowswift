/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_RUNNER_IACT_LEGACY_H
#define SWIFT_RUNNER_IACT_LEGACY_H

/* Includes. */
#include "const.h"
#include "kernel.h"
#include "part.h"
#include "vector.h"

/**
 * @file  runner_iact_legacy.h
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
}

/**
 * @brief Density loop (Vectorized version)
 */
__attribute__((always_inline)) INLINE static void runner_iact_vec_density(
    float *R2, float *Dx, float *Hi, float *Hj, struct part **pi,
    struct part **pj) {

#ifdef VECTORIZE

  vector r, r2, ri, xi, xj, hi, hj, hi_inv, hj_inv, wi, wj, wi_dx, wj_dx;
  vector rhoi, rhoj, rhoi_dh, rhoj_dh, wcounti, wcountj, wcounti_dh, wcountj_dh;
  vector mi, mj;
  vector dx[3], dv[3];
  vector vi[3], vj[3];
  vector dvdr, div_vi, div_vj;
  vector curlvr[3], curl_vi[3], curl_vj[3];
  int k, j;

#if VEC_SIZE == 8
  mi.v = vec_set(pi[0]->mass, pi[1]->mass, pi[2]->mass, pi[3]->mass,
                 pi[4]->mass, pi[5]->mass, pi[6]->mass, pi[7]->mass);
  mj.v = vec_set(pj[0]->mass, pj[1]->mass, pj[2]->mass, pj[3]->mass,
                 pj[4]->mass, pj[5]->mass, pj[6]->mass, pj[7]->mass);
  for (k = 0; k < 3; k++) {
    vi[k].v = vec_set(pi[0]->v[k], pi[1]->v[k], pi[2]->v[k], pi[3]->v[k],
                      pi[4]->v[k], pi[5]->v[k], pi[6]->v[k], pi[7]->v[k]);
    vj[k].v = vec_set(pj[0]->v[k], pj[1]->v[k], pj[2]->v[k], pj[3]->v[k],
                      pj[4]->v[k], pj[5]->v[k], pj[6]->v[k], pj[7]->v[k]);
  }
  for (k = 0; k < 3; k++)
    dx[k].v = vec_set(Dx[0 + k], Dx[3 + k], Dx[6 + k], Dx[9 + k], Dx[12 + k],
                      Dx[15 + k], Dx[18 + k], Dx[21 + k]);
#elif VEC_SIZE == 4
  mi.v = vec_set(pi[0]->mass, pi[1]->mass, pi[2]->mass, pi[3]->mass);
  mj.v = vec_set(pj[0]->mass, pj[1]->mass, pj[2]->mass, pj[3]->mass);
  for (k = 0; k < 3; k++) {
    vi[k].v = vec_set(pi[0]->v[k], pi[1]->v[k], pi[2]->v[k], pi[3]->v[k]);
    vj[k].v = vec_set(pj[0]->v[k], pj[1]->v[k], pj[2]->v[k], pj[3]->v[k]);
  }
  for (k = 0; k < 3; k++)
    dx[k].v = vec_set(Dx[0 + k], Dx[3 + k], Dx[6 + k], Dx[9 + k]);
#endif

  /* Get the radius and inverse radius. */
  r2.v = vec_load(R2);
  ri.v = vec_rsqrt(r2.v);
  ri.v = ri.v - vec_set1(0.5f) * ri.v * (r2.v * ri.v * ri.v - vec_set1(1.0f));
  r.v = r2.v * ri.v;

  hi.v = vec_load(Hi);
  hi_inv.v = vec_rcp(hi.v);
  hi_inv.v = hi_inv.v - hi_inv.v * (hi_inv.v * hi.v - vec_set1(1.0f));
  xi.v = r.v * hi_inv.v;

  hj.v = vec_load(Hj);
  hj_inv.v = vec_rcp(hj.v);
  hj_inv.v = hj_inv.v - hj_inv.v * (hj_inv.v * hj.v - vec_set1(1.0f));
  xj.v = r.v * hj_inv.v;

  kernel_deval_vec(&xi, &wi, &wi_dx);
  kernel_deval_vec(&xj, &wj, &wj_dx);

  /* Compute dv. */
  dv[0].v = vi[0].v - vj[0].v;
  dv[1].v = vi[1].v - vj[1].v;
  dv[2].v = vi[2].v - vj[2].v;

  /* Compute dv dot r */
  dvdr.v = (dv[0].v * dx[0].v) + (dv[1].v * dx[1].v) + (dv[2].v * dx[2].v);
  dvdr.v = dvdr.v * ri.v;

  /* Compute dv cross r */
  curlvr[0].v = dv[1].v * dx[2].v - dv[2].v * dx[1].v;
  curlvr[1].v = dv[2].v * dx[0].v - dv[0].v * dx[2].v;
  curlvr[2].v = dv[0].v * dx[1].v - dv[1].v * dx[0].v;
  for (k = 0; k < 3; k++) curlvr[k].v *= ri.v;

  rhoi.v = mj.v * wi.v;
  rhoi_dh.v = mj.v * (vec_set1(3.0f) * wi.v + xi.v * wi_dx.v);
  wcounti.v = wi.v;
  wcounti_dh.v = xi.v * wi_dx.v;
  div_vi.v = mj.v * dvdr.v * wi_dx.v;
  for (k = 0; k < 3; k++) curl_vi[k].v = mj.v * curlvr[k].v * wi_dx.v;

  rhoj.v = mi.v * wj.v;
  rhoj_dh.v = mi.v * (vec_set1(3.0f) * wj.v + xj.v * wj_dx.v);
  wcountj.v = wj.v;
  wcountj_dh.v = xj.v * wj_dx.v;
  div_vj.v = mi.v * dvdr.v * wj_dx.v;
  for (k = 0; k < 3; k++) curl_vj[k].v = mi.v * curlvr[k].v * wj_dx.v;

  for (k = 0; k < VEC_SIZE; k++) {
    pi[k]->rho += rhoi.f[k];
    pi[k]->rho_dh -= rhoi_dh.f[k];
    pi[k]->density.wcount += wcounti.f[k];
    pi[k]->density.wcount_dh -= wcounti_dh.f[k];
    pi[k]->density.div_v += div_vi.f[k];
    for (j = 0; j < 3; j++) pi[k]->density.curl_v[j] += curl_vi[j].f[k];
    pj[k]->rho += rhoj.f[k];
    pj[k]->rho_dh -= rhoj_dh.f[k];
    pj[k]->density.wcount += wcountj.f[k];
    pj[k]->density.wcount_dh -= wcountj_dh.f[k];
    pj[k]->density.div_v += div_vj.f[k];
    for (j = 0; j < 3; j++) pj[k]->density.curl_v[j] += curl_vj[j].f[k];
  }

#else

  for (int k = 0; k < VEC_SIZE; k++)
    runner_iact_density(R2[k], &Dx[3 * k], Hi[k], Hj[k], pi[k], pj[k]);

#endif
}

/**
 * @brief Density loop (non-symmetric version)
 */

__attribute__((always_inline)) INLINE static void runner_iact_nonsym_density(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

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
}

/**
 * @brief Density loop (non-symmetric vectorized version)
 */

__attribute__((always_inline))
    INLINE static void runner_iact_nonsym_vec_density(float *R2, float *Dx,
                                                      float *Hi, float *Hj,
                                                      struct part **pi,
                                                      struct part **pj) {

#ifdef VECTORIZE

  vector r, r2, ri, xi, hi, hi_inv, wi, wi_dx;
  vector rhoi, rhoi_dh, wcounti, wcounti_dh, div_vi;
  vector mj;
  vector dx[3], dv[3];
  vector vi[3], vj[3];
  vector dvdr;
  vector curlvr[3], curl_vi[3];
  int k, j;

#if VEC_SIZE == 8
  mj.v = vec_set(pj[0]->mass, pj[1]->mass, pj[2]->mass, pj[3]->mass,
                 pj[4]->mass, pj[5]->mass, pj[6]->mass, pj[7]->mass);
  for (k = 0; k < 3; k++) {
    vi[k].v = vec_set(pi[0]->v[k], pi[1]->v[k], pi[2]->v[k], pi[3]->v[k],
                      pi[4]->v[k], pi[5]->v[k], pi[6]->v[k], pi[7]->v[k]);
    vj[k].v = vec_set(pj[0]->v[k], pj[1]->v[k], pj[2]->v[k], pj[3]->v[k],
                      pj[4]->v[k], pj[5]->v[k], pj[6]->v[k], pj[7]->v[k]);
  }
  for (k = 0; k < 3; k++)
    dx[k].v = vec_set(Dx[0 + k], Dx[3 + k], Dx[6 + k], Dx[9 + k], Dx[12 + k],
                      Dx[15 + k], Dx[18 + k], Dx[21 + k]);
#elif VEC_SIZE == 4
  mj.v = vec_set(pj[0]->mass, pj[1]->mass, pj[2]->mass, pj[3]->mass);
  for (k = 0; k < 3; k++) {
    vi[k].v = vec_set(pi[0]->v[k], pi[1]->v[k], pi[2]->v[k], pi[3]->v[k]);
    vj[k].v = vec_set(pj[0]->v[k], pj[1]->v[k], pj[2]->v[k], pj[3]->v[k]);
  }
  for (k = 0; k < 3; k++)
    dx[k].v = vec_set(Dx[0 + k], Dx[3 + k], Dx[6 + k], Dx[9 + k]);
#endif

  /* Get the radius and inverse radius. */
  r2.v = vec_load(R2);
  ri.v = vec_rsqrt(r2.v);
  ri.v = ri.v - vec_set1(0.5f) * ri.v * (r2.v * ri.v * ri.v - vec_set1(1.0f));
  r.v = r2.v * ri.v;

  hi.v = vec_load(Hi);
  hi_inv.v = vec_rcp(hi.v);
  hi_inv.v = hi_inv.v - hi_inv.v * (hi_inv.v * hi.v - vec_set1(1.0f));
  xi.v = r.v * hi_inv.v;

  kernel_deval_vec(&xi, &wi, &wi_dx);

  /* Compute dv. */
  dv[0].v = vi[0].v - vj[0].v;
  dv[1].v = vi[1].v - vj[1].v;
  dv[2].v = vi[2].v - vj[2].v;

  /* Compute dv dot r */
  dvdr.v = (dv[0].v * dx[0].v) + (dv[1].v * dx[1].v) + (dv[2].v * dx[2].v);
  dvdr.v = dvdr.v * ri.v;

  /* Compute dv cross r */
  curlvr[0].v = dv[1].v * dx[2].v - dv[2].v * dx[1].v;
  curlvr[1].v = dv[2].v * dx[0].v - dv[0].v * dx[2].v;
  curlvr[2].v = dv[0].v * dx[1].v - dv[1].v * dx[0].v;
  for (k = 0; k < 3; k++) curlvr[k].v *= ri.v;

  rhoi.v = mj.v * wi.v;
  rhoi_dh.v = mj.v * (vec_set1(3.0f) * wi.v + xi.v * wi_dx.v);
  wcounti.v = wi.v;
  wcounti_dh.v = xi.v * wi_dx.v;
  div_vi.v = mj.v * dvdr.v * wi_dx.v;
  for (k = 0; k < 3; k++) curl_vi[k].v = mj.v * curlvr[k].v * wi_dx.v;

  for (k = 0; k < VEC_SIZE; k++) {
    pi[k]->rho += rhoi.f[k];
    pi[k]->rho_dh -= rhoi_dh.f[k];
    pi[k]->density.wcount += wcounti.f[k];
    pi[k]->density.wcount_dh -= wcounti_dh.f[k];
    pi[k]->density.div_v += div_vi.f[k];
    for (j = 0; j < 3; j++) pi[k]->density.curl_v[j] += curl_vi[j].f[k];
  }

#else

  for (int k = 0; k < VEC_SIZE; k++)
    runner_iact_nonsym_density(R2[k], &Dx[3 * k], Hi[k], Hj[k], pi[k], pj[k]);

#endif
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
  float v_sig, omega_ij, Pi_ij;
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
  v_sig = pi->force.c + pj->force.c - 3.0f * omega_ij;

  /* Compute viscosity tensor */
  Pi_ij = -const_viscosity_alpha * v_sig * omega_ij / (rhoi + rhoj);

  /* Apply balsara switch */
  Pi_ij *= (pi->force.balsara + pj->force.balsara);

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

  /* Get the time derivative for h. */
  pi->force.h_dt -= mj * dvdr / rhoj * wi_dr;
  pj->force.h_dt -= mi * dvdr / rhoi * wj_dr;

  /* Update the signal velocity. */
  pi->force.v_sig = fmaxf(pi->force.v_sig, v_sig);
  pj->force.v_sig = fmaxf(pj->force.v_sig, v_sig);
}

/**
 * @brief Force loop (Vectorized version)
 */

__attribute__((always_inline)) INLINE static void runner_iact_vec_force(
    float *R2, float *Dx, float *Hi, float *Hj, struct part **pi,
    struct part **pj) {

#ifdef VECTORIZE

  vector r, r2, ri;
  vector xi, xj;
  vector hi, hj, hi_inv, hj_inv;
  vector hi2_inv, hj2_inv;
  vector wi, wj, wi_dx, wj_dx, wi_dr, wj_dr, dvdr;
  vector w;
  vector piPOrho2, pjPOrho2, pirho, pjrho;
  vector mi, mj;
  vector f;
  vector dx[3];
  vector vi[3], vj[3];
  vector pia[3], pja[3];
  vector piu_dt, pju_dt;
  vector pih_dt, pjh_dt;
  vector ci, cj, v_sig, vi_sig, vj_sig;
  vector omega_ij, Pi_ij, balsara;
  int j, k;

/* Load stuff. */
#if VEC_SIZE == 8
  mi.v = vec_set(pi[0]->mass, pi[1]->mass, pi[2]->mass, pi[3]->mass,
                 pi[4]->mass, pi[5]->mass, pi[6]->mass, pi[7]->mass);
  mj.v = vec_set(pj[0]->mass, pj[1]->mass, pj[2]->mass, pj[3]->mass,
                 pj[4]->mass, pj[5]->mass, pj[6]->mass, pj[7]->mass);
  piPOrho2.v =
      vec_set(pi[0]->force.POrho2, pi[1]->force.POrho2, pi[2]->force.POrho2,
              pi[3]->force.POrho2, pi[4]->force.POrho2, pi[5]->force.POrho2,
              pi[6]->force.POrho2, pi[7]->force.POrho2);
  pjPOrho2.v =
      vec_set(pj[0]->force.POrho2, pj[1]->force.POrho2, pj[2]->force.POrho2,
              pj[3]->force.POrho2, pj[4]->force.POrho2, pj[5]->force.POrho2,
              pj[6]->force.POrho2, pj[7]->force.POrho2);
  pirho.v = vec_set(pi[0]->rho, pi[1]->rho, pi[2]->rho, pi[3]->rho, pi[4]->rho,
                    pi[5]->rho, pi[6]->rho, pi[7]->rho);
  pjrho.v = vec_set(pj[0]->rho, pj[1]->rho, pj[2]->rho, pj[3]->rho, pj[4]->rho,
                    pj[5]->rho, pj[6]->rho, pj[7]->rho);
  ci.v =
      vec_set(pi[0]->force.c, pi[1]->force.c, pi[2]->force.c, pi[3]->force.c,
              pi[4]->force.c, pi[5]->force.c, pi[6]->force.c, pi[7]->force.c);
  cj.v =
      vec_set(pj[0]->force.c, pj[1]->force.c, pj[2]->force.c, pj[3]->force.c,
              pj[4]->force.c, pj[5]->force.c, pj[6]->force.c, pj[7]->force.c);
  vi_sig.v = vec_set(pi[0]->force.v_sig, pi[1]->force.v_sig, pi[2]->force.v_sig,
                     pi[3]->force.v_sig, pi[4]->force.v_sig, pi[5]->force.v_sig,
                     pi[6]->force.v_sig, pi[7]->force.v_sig);
  vj_sig.v = vec_set(pj[0]->force.v_sig, pj[1]->force.v_sig, pj[2]->force.v_sig,
                     pj[3]->force.v_sig, pj[4]->force.v_sig, pj[5]->force.v_sig,
                     pj[6]->force.v_sig, pj[7]->force.v_sig);
  for (k = 0; k < 3; k++) {
    vi[k].v = vec_set(pi[0]->v[k], pi[1]->v[k], pi[2]->v[k], pi[3]->v[k],
                      pi[4]->v[k], pi[5]->v[k], pi[6]->v[k], pi[7]->v[k]);
    vj[k].v = vec_set(pj[0]->v[k], pj[1]->v[k], pj[2]->v[k], pj[3]->v[k],
                      pj[4]->v[k], pj[5]->v[k], pj[6]->v[k], pj[7]->v[k]);
  }
  for (k = 0; k < 3; k++)
    dx[k].v = vec_set(Dx[0 + k], Dx[3 + k], Dx[6 + k], Dx[9 + k], Dx[12 + k],
                      Dx[15 + k], Dx[18 + k], Dx[21 + k]);
  balsara.v =
      vec_set(pi[0]->force.balsara, pi[1]->force.balsara, pi[2]->force.balsara,
              pi[3]->force.balsara, pi[4]->force.balsara, pi[5]->force.balsara,
              pi[6]->force.balsara, pi[7]->force.balsara) +
      vec_set(pj[0]->force.balsara, pj[1]->force.balsara, pj[2]->force.balsara,
              pj[3]->force.balsara, pj[4]->force.balsara, pj[5]->force.balsara,
              pj[6]->force.balsara, pj[7]->force.balsara);
#elif VEC_SIZE == 4
  mi.v = vec_set(pi[0]->mass, pi[1]->mass, pi[2]->mass, pi[3]->mass);
  mj.v = vec_set(pj[0]->mass, pj[1]->mass, pj[2]->mass, pj[3]->mass);
  piPOrho2.v = vec_set(pi[0]->force.POrho2, pi[1]->force.POrho2,
                       pi[2]->force.POrho2, pi[3]->force.POrho2);
  pjPOrho2.v = vec_set(pj[0]->force.POrho2, pj[1]->force.POrho2,
                       pj[2]->force.POrho2, pj[3]->force.POrho2);
  pirho.v = vec_set(pi[0]->rho, pi[1]->rho, pi[2]->rho, pi[3]->rho);
  pjrho.v = vec_set(pj[0]->rho, pj[1]->rho, pj[2]->rho, pj[3]->rho);
  ci.v =
      vec_set(pi[0]->force.c, pi[1]->force.c, pi[2]->force.c, pi[3]->force.c);
  cj.v =
      vec_set(pj[0]->force.c, pj[1]->force.c, pj[2]->force.c, pj[3]->force.c);
  vi_sig.v = vec_set(pi[0]->force.v_sig, pi[1]->force.v_sig, pi[2]->force.v_sig,
                     pi[3]->force.v_sig);
  vj_sig.v = vec_set(pj[0]->force.v_sig, pj[1]->force.v_sig, pj[2]->force.v_sig,
                     pj[3]->force.v_sig);
  for (k = 0; k < 3; k++) {
    vi[k].v = vec_set(pi[0]->v[k], pi[1]->v[k], pi[2]->v[k], pi[3]->v[k]);
    vj[k].v = vec_set(pj[0]->v[k], pj[1]->v[k], pj[2]->v[k], pj[3]->v[k]);
  }
  for (k = 0; k < 3; k++)
    dx[k].v = vec_set(Dx[0 + k], Dx[3 + k], Dx[6 + k], Dx[9 + k]);
  balsara.v = vec_set(pi[0]->force.balsara, pi[1]->force.balsara,
                      pi[2]->force.balsara, pi[3]->force.balsara) +
              vec_set(pj[0]->force.balsara, pj[1]->force.balsara,
                      pj[2]->force.balsara, pj[3]->force.balsara);
#else
#error
#endif

  /* Get the radius and inverse radius. */
  r2.v = vec_load(R2);
  ri.v = vec_rsqrt(r2.v);
  ri.v = ri.v - vec_set1(0.5f) * ri.v * (r2.v * ri.v * ri.v - vec_set1(1.0f));
  r.v = r2.v * ri.v;

  /* Get the kernel for hi. */
  hi.v = vec_load(Hi);
  hi_inv.v = vec_rcp(hi.v);
  hi_inv.v = hi_inv.v - hi_inv.v * (hi.v * hi_inv.v - vec_set1(1.0f));
  hi2_inv.v = hi_inv.v * hi_inv.v;
  xi.v = r.v * hi_inv.v;
  kernel_deval_vec(&xi, &wi, &wi_dx);
  wi_dr.v = hi2_inv.v * hi2_inv.v * wi_dx.v;

  /* Get the kernel for hj. */
  hj.v = vec_load(Hj);
  hj_inv.v = vec_rcp(hj.v);
  hj_inv.v = hj_inv.v - hj_inv.v * (hj.v * hj_inv.v - vec_set1(1.0f));
  hj2_inv.v = hj_inv.v * hj_inv.v;
  xj.v = r.v * hj_inv.v;
  kernel_deval_vec(&xj, &wj, &wj_dx);
  wj_dr.v = hj2_inv.v * hj2_inv.v * wj_dx.v;

  /* Compute dv dot r. */
  dvdr.v = ((vi[0].v - vj[0].v) * dx[0].v) + ((vi[1].v - vj[1].v) * dx[1].v) +
           ((vi[2].v - vj[2].v) * dx[2].v);
  dvdr.v = dvdr.v * ri.v;

  /* Get the time derivative for h. */
  pih_dt.v = mj.v / pjrho.v * dvdr.v * wi_dr.v;
  pjh_dt.v = mi.v / pirho.v * dvdr.v * wj_dr.v;

  /* Compute the relative velocity. (This is 0 if the particles move away from
   * each other and negative otherwise) */
  omega_ij.v = vec_fmin(dvdr.v, vec_set1(0.0f));

  /* Compute signal velocity */
  v_sig.v = ci.v + cj.v - vec_set1(3.0f) * omega_ij.v;

  /* Compute viscosity tensor */
  Pi_ij.v = -balsara.v * vec_set1(const_viscosity_alpha) * v_sig.v *
            omega_ij.v / (pirho.v + pjrho.v);
  Pi_ij.v *= (wi_dr.v + wj_dr.v);

  /* Get the common factor out. */
  w.v = ri.v * ((piPOrho2.v * wi_dr.v + pjPOrho2.v * wj_dr.v) +
                vec_set1(0.25f) * Pi_ij.v);

  /* Use the force, Luke! */
  for (k = 0; k < 3; k++) {
    f.v = dx[k].v * w.v;
    pia[k].v = mj.v * f.v;
    pja[k].v = mi.v * f.v;
  }

  /* Get the time derivative for u. */
  piu_dt.v =
      mj.v * dvdr.v * (piPOrho2.v * wi_dr.v + vec_set1(0.125f) * Pi_ij.v);
  pju_dt.v =
      mi.v * dvdr.v * (pjPOrho2.v * wj_dr.v + vec_set1(0.125f) * Pi_ij.v);

  /* compute the signal velocity (this is always symmetrical). */
  vi_sig.v = vec_fmax(vi_sig.v, v_sig.v);
  vj_sig.v = vec_fmax(vj_sig.v, v_sig.v);

  /* Store the forces back on the particles. */
  for (k = 0; k < VEC_SIZE; k++) {
    pi[k]->force.u_dt += piu_dt.f[k];
    pj[k]->force.u_dt += pju_dt.f[k];
    pi[k]->force.h_dt -= pih_dt.f[k];
    pj[k]->force.h_dt -= pjh_dt.f[k];
    pi[k]->force.v_sig = vi_sig.f[k];
    pj[k]->force.v_sig = vj_sig.f[k];
    for (j = 0; j < 3; j++) {
      pi[k]->a[j] -= pia[j].f[k];
      pj[k]->a[j] += pja[j].f[k];
    }
  }

#else

  for (int k = 0; k < VEC_SIZE; k++)
    runner_iact_force(R2[k], &Dx[3 * k], Hi[k], Hj[k], pi[k], pj[k]);

#endif
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
  float v_sig, omega_ij, Pi_ij;
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
  v_sig = pi->force.c + pj->force.c - 3.0f * omega_ij;

  /* Compute viscosity tensor */
  Pi_ij = -const_viscosity_alpha * v_sig * omega_ij / (rhoi + rhoj);

  /* Apply balsara switch */
  Pi_ij *= (pi->force.balsara + pj->force.balsara);

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

  /* Get the time derivative for h. */
  pi->force.h_dt -= mj * dvdr / rhoj * wi_dr;

  /* Update the signal velocity. */
  pi->force.v_sig = fmaxf(pi->force.v_sig, v_sig);
  pj->force.v_sig = fmaxf(pj->force.v_sig, v_sig);
}

/**
 * @brief Force loop (Vectorized non-symmetric version)
 */

__attribute__((always_inline)) INLINE static void runner_iact_nonsym_vec_force(
    float *R2, float *Dx, float *Hi, float *Hj, struct part **pi,
    struct part **pj) {

#ifdef VECTORIZE

  vector r, r2, ri;
  vector xi, xj;
  vector hi, hj, hi_inv, hj_inv;
  vector hi2_inv, hj2_inv;
  vector wi, wj, wi_dx, wj_dx, wi_dr, wj_dr, dvdr;
  vector w;
  vector piPOrho2, pjPOrho2, pirho, pjrho;
  vector mj;
  vector f;
  vector dx[3];
  vector vi[3], vj[3];
  vector pia[3];
  vector piu_dt;
  vector pih_dt;
  vector ci, cj, v_sig, vi_sig, vj_sig;
  vector omega_ij, Pi_ij, balsara;
  int j, k;

/* Load stuff. */
#if VEC_SIZE == 8
  mj.v = vec_set(pj[0]->mass, pj[1]->mass, pj[2]->mass, pj[3]->mass,
                 pj[4]->mass, pj[5]->mass, pj[6]->mass, pj[7]->mass);
  piPOrho2.v =
      vec_set(pi[0]->force.POrho2, pi[1]->force.POrho2, pi[2]->force.POrho2,
              pi[3]->force.POrho2, pi[4]->force.POrho2, pi[5]->force.POrho2,
              pi[6]->force.POrho2, pi[7]->force.POrho2);
  pjPOrho2.v =
      vec_set(pj[0]->force.POrho2, pj[1]->force.POrho2, pj[2]->force.POrho2,
              pj[3]->force.POrho2, pj[4]->force.POrho2, pj[5]->force.POrho2,
              pj[6]->force.POrho2, pj[7]->force.POrho2);
  pirho.v = vec_set(pi[0]->rho, pi[1]->rho, pi[2]->rho, pi[3]->rho, pi[4]->rho,
                    pi[5]->rho, pi[6]->rho, pi[7]->rho);
  pjrho.v = vec_set(pj[0]->rho, pj[1]->rho, pj[2]->rho, pj[3]->rho, pj[4]->rho,
                    pj[5]->rho, pj[6]->rho, pj[7]->rho);
  ci.v =
      vec_set(pi[0]->force.c, pi[1]->force.c, pi[2]->force.c, pi[3]->force.c,
              pi[4]->force.c, pi[5]->force.c, pi[6]->force.c, pi[7]->force.c);
  cj.v =
      vec_set(pj[0]->force.c, pj[1]->force.c, pj[2]->force.c, pj[3]->force.c,
              pj[4]->force.c, pj[5]->force.c, pj[6]->force.c, pj[7]->force.c);
  vi_sig.v = vec_set(pi[0]->force.v_sig, pi[1]->force.v_sig, pi[2]->force.v_sig,
                     pi[3]->force.v_sig, pi[4]->force.v_sig, pi[5]->force.v_sig,
                     pi[6]->force.v_sig, pi[7]->force.v_sig);
  vj_sig.v = vec_set(pj[0]->force.v_sig, pj[1]->force.v_sig, pj[2]->force.v_sig,
                     pj[3]->force.v_sig, pj[4]->force.v_sig, pj[5]->force.v_sig,
                     pj[6]->force.v_sig, pj[7]->force.v_sig);
  for (k = 0; k < 3; k++) {
    vi[k].v = vec_set(pi[0]->v[k], pi[1]->v[k], pi[2]->v[k], pi[3]->v[k],
                      pi[4]->v[k], pi[5]->v[k], pi[6]->v[k], pi[7]->v[k]);
    vj[k].v = vec_set(pj[0]->v[k], pj[1]->v[k], pj[2]->v[k], pj[3]->v[k],
                      pj[4]->v[k], pj[5]->v[k], pj[6]->v[k], pj[7]->v[k]);
  }
  for (k = 0; k < 3; k++)
    dx[k].v = vec_set(Dx[0 + k], Dx[3 + k], Dx[6 + k], Dx[9 + k], Dx[12 + k],
                      Dx[15 + k], Dx[18 + k], Dx[21 + k]);
  balsara.v =
      vec_set(pi[0]->force.balsara, pi[1]->force.balsara, pi[2]->force.balsara,
              pi[3]->force.balsara, pi[4]->force.balsara, pi[5]->force.balsara,
              pi[6]->force.balsara, pi[7]->force.balsara) +
      vec_set(pj[0]->force.balsara, pj[1]->force.balsara, pj[2]->force.balsara,
              pj[3]->force.balsara, pj[4]->force.balsara, pj[5]->force.balsara,
              pj[6]->force.balsara, pj[7]->force.balsara);
#elif VEC_SIZE == 4
  mj.v = vec_set(pj[0]->mass, pj[1]->mass, pj[2]->mass, pj[3]->mass);
  piPOrho2.v = vec_set(pi[0]->force.POrho2, pi[1]->force.POrho2,
                       pi[2]->force.POrho2, pi[3]->force.POrho2);
  pjPOrho2.v = vec_set(pj[0]->force.POrho2, pj[1]->force.POrho2,
                       pj[2]->force.POrho2, pj[3]->force.POrho2);
  pirho.v = vec_set(pi[0]->rho, pi[1]->rho, pi[2]->rho, pi[3]->rho);
  pjrho.v = vec_set(pj[0]->rho, pj[1]->rho, pj[2]->rho, pj[3]->rho);
  ci.v =
      vec_set(pi[0]->force.c, pi[1]->force.c, pi[2]->force.c, pi[3]->force.c);
  cj.v =
      vec_set(pj[0]->force.c, pj[1]->force.c, pj[2]->force.c, pj[3]->force.c);
  vi_sig.v = vec_set(pi[0]->force.v_sig, pi[1]->force.v_sig, pi[2]->force.v_sig,
                     pi[3]->force.v_sig);
  vj_sig.v = vec_set(pj[0]->force.v_sig, pj[1]->force.v_sig, pj[2]->force.v_sig,
                     pj[3]->force.v_sig);
  for (k = 0; k < 3; k++) {
    vi[k].v = vec_set(pi[0]->v[k], pi[1]->v[k], pi[2]->v[k], pi[3]->v[k]);
    vj[k].v = vec_set(pj[0]->v[k], pj[1]->v[k], pj[2]->v[k], pj[3]->v[k]);
  }
  for (k = 0; k < 3; k++)
    dx[k].v = vec_set(Dx[0 + k], Dx[3 + k], Dx[6 + k], Dx[9 + k]);
  balsara.v = vec_set(pi[0]->force.balsara, pi[1]->force.balsara,
                      pi[2]->force.balsara, pi[3]->force.balsara) +
              vec_set(pj[0]->force.balsara, pj[1]->force.balsara,
                      pj[2]->force.balsara, pj[3]->force.balsara);
#else
#error
#endif

  /* Get the radius and inverse radius. */
  r2.v = vec_load(R2);
  ri.v = vec_rsqrt(r2.v);
  ri.v = ri.v - vec_set1(0.5f) * ri.v * (r2.v * ri.v * ri.v - vec_set1(1.0f));
  r.v = r2.v * ri.v;

  /* Get the kernel for hi. */
  hi.v = vec_load(Hi);
  hi_inv.v = vec_rcp(hi.v);
  hi_inv.v = hi_inv.v - hi_inv.v * (hi.v * hi_inv.v - vec_set1(1.0f));
  hi2_inv.v = hi_inv.v * hi_inv.v;
  xi.v = r.v * hi_inv.v;
  kernel_deval_vec(&xi, &wi, &wi_dx);
  wi_dr.v = hi2_inv.v * hi2_inv.v * wi_dx.v;

  /* Get the kernel for hj. */
  hj.v = vec_load(Hj);
  hj_inv.v = vec_rcp(hj.v);
  hj_inv.v = hj_inv.v - hj_inv.v * (hj.v * hj_inv.v - vec_set1(1.0f));
  hj2_inv.v = hj_inv.v * hj_inv.v;
  xj.v = r.v * hj_inv.v;
  kernel_deval_vec(&xj, &wj, &wj_dx);
  wj_dr.v = hj2_inv.v * hj2_inv.v * wj_dx.v;

  /* Compute dv dot r. */
  dvdr.v = ((vi[0].v - vj[0].v) * dx[0].v) + ((vi[1].v - vj[1].v) * dx[1].v) +
           ((vi[2].v - vj[2].v) * dx[2].v);
  dvdr.v = dvdr.v * ri.v;

  /* Get the time derivative for h. */
  pih_dt.v = mj.v / pjrho.v * dvdr.v * wi_dr.v;

  /* Compute the relative velocity. (This is 0 if the particles move away from
   * each other and negative otherwise) */
  omega_ij.v = vec_fmin(dvdr.v, vec_set1(0.0f));

  /* Compute signal velocity */
  v_sig.v = ci.v + cj.v - vec_set1(3.0f) * omega_ij.v;

  /* Compute viscosity tensor */
  Pi_ij.v = -balsara.v * vec_set1(const_viscosity_alpha) * v_sig.v *
            omega_ij.v / (pirho.v + pjrho.v);
  Pi_ij.v *= (wi_dr.v + wj_dr.v);

  /* Get the common factor out. */
  w.v = ri.v * ((piPOrho2.v * wi_dr.v + pjPOrho2.v * wj_dr.v) +
                vec_set1(0.25f) * Pi_ij.v);

  /* Use the force, Luke! */
  for (k = 0; k < 3; k++) {
    f.v = dx[k].v * w.v;
    pia[k].v = mj.v * f.v;
  }

  /* Get the time derivative for u. */
  piu_dt.v =
      mj.v * dvdr.v * (piPOrho2.v * wi_dr.v + vec_set1(0.125f) * Pi_ij.v);

  /* compute the signal velocity (this is always symmetrical). */
  vi_sig.v = vec_fmax(vi_sig.v, v_sig.v);
  vj_sig.v = vec_fmax(vj_sig.v, v_sig.v);

  /* Store the forces back on the particles. */
  for (k = 0; k < VEC_SIZE; k++) {
    pi[k]->force.u_dt += piu_dt.f[k];
    pi[k]->force.h_dt -= pih_dt.f[k];
    pi[k]->force.v_sig = vi_sig.f[k];
    pj[k]->force.v_sig = vj_sig.f[k];
    for (j = 0; j < 3; j++) pi[k]->a[j] -= pia[j].f[k];
  }

#else

  for (int k = 0; k < VEC_SIZE; k++)
    runner_iact_nonsym_force(R2[k], &Dx[3 * k], Hi[k], Hj[k], pi[k], pj[k]);

#endif
}

#endif /* SWIFT_RUNNER_IACT_LEGACY_H */
