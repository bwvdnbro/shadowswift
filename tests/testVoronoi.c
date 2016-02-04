/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
#include "part.h"
#include "voronoi.h"

/**
 * @brief Check if voronoi_volume_tetrahedron() works
 */
void test_voronoi_volume_tetrahedron(){
  float v1[3] = {0., 0., 0.};
  float v2[3] = {0., 0., 1.};
  float v3[3] = {0., 1., 0.};
  float v4[3] = {1., 0., 0.};
  
  float V = voronoi_volume_tetrahedron(v1, v2, v3, v4);
  assert(V == 1.0f/6.0f);
}

/**
 * @brief Check if voronoi_calculate_cell() works
 */
void test_calculate_cell(){
  struct part p;
  p.x[0] = 0.5;
  p.x[1] = 0.5;
  p.x[2] = 0.5;
  voronoi_initialize(&p);
  voronoi_calculate_cell(&p);
  message("volume: %g", p.voronoi.volume);

  assert(p.voronoi.volume == 1.0f);
  assert(p.voronoi.centroid[0] = 0.5f);
  assert(p.voronoi.centroid[1] = 0.5f);
  assert(p.voronoi.centroid[2] = 0.5f);
}

int main(int argc, char **argv){
  test_voronoi_volume_tetrahedron();
  test_calculate_cell();
/*  test_calculate_faces();*/
/*  test_calculate_centroid_tetrahedron();*/
/*  test_cell_intersections();*/
/*  test_degeneracies();*/
}
