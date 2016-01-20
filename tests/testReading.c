/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2015 Matthieu Schaller (matthieu.schaller@durham.ac.uk).
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

#include "swift.h"
#include <stdlib.h>

int main() {

  int N = -1, periodic = -1;
  int i, j, k, n;
  double dim[3];
  struct part *parts = NULL;

  /* Properties of the ICs */
  const double boxSize = 1.;
  const int L = 4;
  const double rho = 2.;
  const double P = 1.;
  const double gamma = 5. / 3.;

  /* Read data */
  read_ic_single("input.hdf5", dim, &parts, &N, &periodic);

  /* Check global properties read are correct */
  assert(dim[0] == boxSize);
  assert(dim[1] == boxSize);
  assert(dim[2] == boxSize);
  assert(N == L * L * L);
  assert(periodic == 1);

  /* Check particles */
  for (n = 0; n < N; ++n) {

    /* Check that indices are in a reasonable range */
    unsigned long long index = parts[n].id;
    assert(index < N);

    /* Check masses */
    float mass = parts[n].mass;
    float correct_mass = boxSize * boxSize * boxSize * rho / N;
    assert(mass == correct_mass);

    /* Check smoothing length */
    float h = parts[n].h;
    float correct_h = 2.251 * boxSize / L;
    assert(h == correct_h);

    /* Check internal energy */
    float u = parts[n].u;
    float correct_u = P / ((gamma - 1.) * rho);
    assert(u == correct_u);

    /* Check velocity */
    assert(parts[n].v[0] == 0.);
    assert(parts[n].v[1] == 0.);
    assert(parts[n].v[2] == 0.);

    /* Check positions */
    k = index % 4;
    j = ((index - k) / 4) % 4;
    i = (index - k - 4 * j) / 16;
    double correct_x = i * boxSize / L + boxSize / (2 * L);
    double correct_y = j * boxSize / L + boxSize / (2 * L);
    double correct_z = k * boxSize / L + boxSize / (2 * L);
    assert(parts[n].x[0] == correct_x);
    assert(parts[n].x[1] == correct_y);
    assert(parts[n].x[2] == correct_z);

    /* Check accelerations */
    assert(parts[n].a[0] == 0.);
    assert(parts[n].a[1] == 0.);
    assert(parts[n].a[2] == 0.);
  }

  /* Clean-up */
  free(parts);

  return 0;
}
