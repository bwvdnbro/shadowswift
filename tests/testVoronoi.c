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
 * @brief Check if voronoi_centroid_tetrahedron() works
 */
void test_voronoi_centroid_tetrahedron(){
  float v1[3] = {0., 0., 0.};
  float v2[3] = {0., 0., 1.};
  float v3[3] = {0., 1., 0.};
  float v4[3] = {1., 0., 0.};
  
  float centroid[3];
  voronoi_centroid_tetrahedron(centroid, v1, v2, v3, v4);
  assert(centroid[0] == 0.25f);
  assert(centroid[1] == 0.25f);
  assert(centroid[2] == 0.25f);
}

/**
 * @brief Check if voronoi_calculate_cell() works
 */
void test_calculate_cell(){
  struct voronoi_box b;
  b.anchor[0] = 0.0f;
  b.anchor[1] = 0.0f;
  b.anchor[2] = 0.0f;
  b.sides[0] = 1.0f;
  b.sides[1] = 1.0f;
  b.sides[2] = 1.0f;
  struct part p;
  p.x[0] = 0.5;
  p.x[1] = 0.5;
  p.x[2] = 0.5;
  voronoi_initialize(&p, &b);
  voronoi_calculate_cell(&p);

  assert(p.voronoi.volume == 1.0f);
  assert(p.voronoi.centroid[0] = 0.5f);
  assert(p.voronoi.centroid[1] = 0.5f);
  assert(p.voronoi.centroid[2] = 0.5f);
}

/**
 * @brief Check if voronoi_calculate_faces() works
 */
void test_voronoi_calculate_faces(){
  struct voronoi_box b;
  b.anchor[0] = 0.0f;
  b.anchor[1] = 0.0f;
  b.anchor[2] = 0.0f;
  b.sides[0] = 1.0f;
  b.sides[1] = 1.0f;
  b.sides[2] = 1.0f;
  struct part p;
  p.x[0] = 0.5;
  p.x[1] = 0.5;
  p.x[2] = 0.5;
  voronoi_initialize(&p, &b);
  voronoi_calculate_faces(&p);

  assert(p.voronoi.face_areas[0] == 1.0f);
  assert(p.voronoi.face_areas[1] == 1.0f);
  assert(p.voronoi.face_areas[2] == 1.0f);
  assert(p.voronoi.face_areas[3] == 1.0f);
  assert(p.voronoi.face_areas[4] == 1.0f);
  assert(p.voronoi.face_areas[5] == 1.0f);

  assert(p.voronoi.face_midpoints[ 0] == 0.5f);
  assert(p.voronoi.face_midpoints[ 1] == 0.0f);
  assert(p.voronoi.face_midpoints[ 2] == 0.5f);

  assert(p.voronoi.face_midpoints[ 3] == 0.0f);
  assert(p.voronoi.face_midpoints[ 4] == 0.5f);
  assert(p.voronoi.face_midpoints[ 5] == 0.5f);

  assert(p.voronoi.face_midpoints[ 6] == 0.5f);
  assert(p.voronoi.face_midpoints[ 7] == 0.5f);
  assert(p.voronoi.face_midpoints[ 8] == 0.0f);

  assert(p.voronoi.face_midpoints[ 9] == 0.5f);
  assert(p.voronoi.face_midpoints[10] == 0.5f);
  assert(p.voronoi.face_midpoints[11] == 1.0f);

  assert(p.voronoi.face_midpoints[12] == 0.5f);
  assert(p.voronoi.face_midpoints[13] == 1.0f);
  assert(p.voronoi.face_midpoints[14] == 0.5f);

  assert(p.voronoi.face_midpoints[15] == 1.0f);
  assert(p.voronoi.face_midpoints[16] == 0.5f);
  assert(p.voronoi.face_midpoints[17] == 0.5f);

  assert(p.voronoi.ngbs[0] == VORONOI_BOX_FRONT);
  assert(p.voronoi.ngbs[1] == VORONOI_BOX_LEFT);
  assert(p.voronoi.ngbs[2] == VORONOI_BOX_BOTTOM);
  assert(p.voronoi.ngbs[3] == VORONOI_BOX_TOP);
  assert(p.voronoi.ngbs[4] == VORONOI_BOX_BACK);
  assert(p.voronoi.ngbs[5] == VORONOI_BOX_RIGHT);
}

/**
 * @brief Convert a float to an unsigned integer
 *
 * We use this to exactly compare float results with round off
 */
unsigned int get_bytes(float f){
  union{
    float f;
    unsigned int u;
  } u;
  u.f = f;
  return u.u;
}

/**
 * @brief Convert an unsigned integer to a float
 *
 * This can be used to exactly reproduce a floating point result with round off
 */
float get_float(unsigned int ui){
  union{
    float f;
    unsigned int u;
  } u;
  u.u = ui;
  return u.f;
}

/**
 * @brief Test if the contents of the cell corresponds to what we expect
 */
void test_cell_contents(struct part *p, int nv, float *v, int *o, int *e,
                        int *f, unsigned long long *n){
  assert(p->geometry.nvert == nv);
  int k = 0;
  for(int i = 0; i < nv; i++){
    fprintf(stderr, "%i == %i?\n", p->geometry.orders[i], o[i]);
    assert(p->geometry.orders[i] == o[i]);
    for(int j = 0; j < 3; j++){
      fprintf(stderr, "%g == %g?\n", p->geometry.vertices[3*i+j], v[3*i+j]);
      assert(p->geometry.vertices[3*i+j] == v[3*i+j]);
    }
    for(int j = 0; j < o[i]; j++){
      fprintf(stderr, "%i == %i?\n", p->geometry.edges[k+j], e[k+j]);
      assert(p->geometry.edges[k+j] == e[k+j]);

      fprintf(stderr, "%i == %i?\n", p->geometry.edgeindices[k+j], f[k+j]);
      assert(p->geometry.edgeindices[k+j] == f[k+j]);

      fprintf(stderr, "%llu == %llu?\n", p->geometry.ngbs[k+j], n[k+j]);
      assert(p->geometry.ngbs[o[i]+j] == n[o[i]+j]);
    }
    k += o[i];
  }
}

/**
 * @brief Test a single cell intersection and verify the results
 */
void test_cell_intersection(float x1, float y1, float z1, float x2,
                            float y2, float z2, int nv, float *v, int *o,
                            int *e, int *f, unsigned long long *n){
  struct voronoi_box b;
  b.anchor[0] = 0.0f;
  b.anchor[1] = 0.0f;
  b.anchor[2] = 0.0f;
  b.sides[0] = 1.0f;
  b.sides[1] = 1.0f;
  b.sides[2] = 1.0f;
  struct part particles[2];
  particles[0].x[0] = x1;
  particles[0].x[1] = y1;
  particles[0].x[2] = z1;
  particles[0].id = 0;
  particles[1].x[0] = x2;
  particles[1].x[1] = y2;
  particles[1].x[2] = z2;
  particles[1].id = 1;
  voronoi_initialize(&particles[0], &b);
  voronoi_initialize(&particles[1], &b);
  float dx[3];
  dx[0] = x1 - x2;
  dx[1] = y1 - y2;
  dx[2] = z1 - z2;
  voronoi_intersect(dx, &particles[0], &particles[1]);
  test_cell_contents(&particles[0], nv, v, o, e, f, n);
}

/**
 * @brief Test general cell intersections
 */
void test_cell_intersections(){
  {
    float r = 0.45f;
    float l = 1.0f - r;
    float v[24] = {-0.5f, -0.5f, -0.5f,
                   -0.5f, 0.5f, -0.5f,
                   0.5f, -0.5f, -0.5f,
                   0.5f, 0.5f, -0.5f,
                   -0.5f, -0.5f, -0.5f*r + 0.5f*l,
                   0.5f, -0.5f, -0.5f*r + 0.5f*l,
                   0.5f, 0.5f, -0.5f*r + 0.5f*l,
                   -0.5f, 0.5f, -0.5f*r + 0.5f*l};
    int o[8] = {3, 3, 3, 3, 3, 3, 3, 3};
    int e[24] = {4, 1, 2,
                 7, 3, 0,
                 0, 3, 5,
                 1, 6, 2,
                 7, 0, 5,
                 4, 2, 6,
                 5, 3, 7,
                 6, 1, 4};
    int f[24] = {1, 2, 0,
                 1, 0, 1,
                 2, 2, 1,
                 1, 1, 1,
                 2, 0, 0,
                 2, 2, 0,
                 2, 1, 0,
                 2, 0, 0};
    unsigned long long n[24] =
      {VORONOI_BOX_FRONT, VORONOI_BOX_LEFT, VORONOI_BOX_BOTTOM,
       VORONOI_BOX_LEFT, VORONOI_BOX_BACK, VORONOI_BOX_BOTTOM,
       VORONOI_BOX_FRONT, VORONOI_BOX_BOTTOM, VORONOI_BOX_RIGHT,
       VORONOI_BOX_BOTTOM, VORONOI_BOX_BACK, VORONOI_BOX_RIGHT,
       1, VORONOI_BOX_LEFT, VORONOI_BOX_FRONT,
       1, VORONOI_BOX_FRONT, VORONOI_BOX_RIGHT,
       1, VORONOI_BOX_RIGHT, VORONOI_BOX_BACK,
       1, VORONOI_BOX_BACK, VORONOI_BOX_LEFT};
    test_cell_intersection(0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.6f,
                           8, v, o, e, f, n);
  }
  /* add additional cases below */
}

/**
 * @brief Set the coordinates of a particle to the given coordinates and
 *  initialize the corresponding Voronoi cell
 */
void set_coordinates(struct part *p, double x, double y, double z,
                     unsigned int id){
  struct voronoi_box b;
  b.anchor[0] = 0.0f;
  b.anchor[1] = 0.0f;
  b.anchor[2] = 0.0f;
  b.sides[0] = 1.0f;
  b.sides[1] = 1.0f;
  b.sides[2] = 1.0f;
  /* set coordinates */
  p->x[0] = x;
  p->x[1] = y;
  p->x[2] = z;
  p->id = id;
  voronoi_initialize(p, &b);
}

/**
 * @brief Test degenerate cell intersections by subjecting a small cube to a
 *  number of consecutive intersections
 */
void test_degeneracies(){
  unsigned int idx = 0;
  /* make a small cube */
  struct part particles[100];
  set_coordinates(&particles[idx], 0.1, 0.1, 0.1, idx);
  idx++;
  set_coordinates(&particles[idx], 0.2, 0.1, 0.1, idx);
  idx++;
  set_coordinates(&particles[idx], 0.1, 0.2, 0.1, idx);
  idx++;
  set_coordinates(&particles[idx], 0.1, 0.1, 0.2, idx);
  idx++;
  /* corner on cutting plane */
  set_coordinates(&particles[idx], 0.2, 0.2, 0.2, idx);
  idx++;
  /* edge on cutting plane */
  set_coordinates(&particles[idx], 0.2, 0.1, 0.2, idx);
  idx++;
  set_coordinates(&particles[idx], 0.2, 0.2, 0.1, idx);
  idx++;
  /* cutting plane is diagonal */
  set_coordinates(&particles[idx], 0.05, 0.1, 0.05, idx);
  idx++;
  /* order 4 vertex (found after an impressive display of analytical geometry
     of which I'm rather proud) */
  float t = 0.5/0.0475;
  set_coordinates(&particles[idx], 0.0075*t+0.1, 0.0075*t+0.1, 0.1-0.0025*t,
                  idx);
  idx++;
  /* order 4 vertex with float edge */
  t = 0.35/0.06125;
  set_coordinates(&particles[idx], 0.0075*t+0.1, 0.015*t+0.1, 0.1-0.005*t,
                  idx);
  idx++;
  /* plane that was already encountered */
  t = 0.5/0.0475;
  set_coordinates(&particles[idx], 0.0075*t+0.1, 0.0075*t+0.1, 0.1-0.0025*t,
                  idx);
  idx++;
  /* no intersection (just to cover all code) */
  set_coordinates(&particles[idx], 0.3, 0.3, 0.3, idx);
  idx++;
  set_coordinates(&particles[idx], 0.3, 0.1, 0.3, idx);
  idx++;
  /* order 5 vertex */
  t = 0.04/0.0175;
  set_coordinates(&particles[idx], 0.1-0.0075*t, 0.1+0.00375*t,
                  0.1+0.00625*t, idx);
  idx++;
  /* plane with order 5 vertex */
  set_coordinates(&particles[idx], 0.1, 0.2, 0.1, idx);
  idx++;
  /* edge with order 5 vertex that looses an edge */
  t = -0.1/0.095;
  set_coordinates(&particles[idx], 0.1-0.015*t, 0.1+0.015*t, 0.1-0.005*t,
                  idx);
  idx++;

  for(int i = 1; i < idx; i++){
    float dx[3];
    dx[0] = particles[0].x[0] - particles[i].x[0];
    dx[1] = particles[0].x[1] - particles[i].x[1];
    dx[2] = particles[0].x[2] - particles[i].x[2];
    voronoi_intersect(dx, &particles[0], &particles[i]);
  }

  /* due to round off error, we have to store these values in binary format */
  float v[15] =
    {get_float(3156465441), get_float(1028443350), get_float(3173242620),
     get_float(3163877883), get_float(1028443350), get_float(3171046343),
     get_float(1028443350), get_float(1028443350), get_float(1028443350),
     get_float(3175926984), get_float(3179841669), get_float(841482240),
     get_float(1028443350), get_float(3006267392), get_float(3184315597)};
  int o[5] = {3, 3, 4, 3, 3};
  int e[16] = {2, 4, 1,
               0, 3, 2,
               4, 0, 1, 3,
               2, 1, 4,
               3, 0, 2};
  int f[16] = {1, 1, 0,
               2, 1, 2,
               2, 0, 2, 0,
               3, 1, 0,
               2, 1, 0};
  unsigned long long n[16] = {14, 9, 7,
                              14, 7, 13,
                              15, 9, 14, 13,
                              15, 13, 7,
                              15, 7, 9};
  test_cell_contents(&particles[0], 5, v, o, e, f, n);
}

/**
 * @brief Test the Voronoi cell algorithms
 */
int main(int argc, char **argv){
  test_voronoi_volume_tetrahedron();
  test_voronoi_centroid_tetrahedron();
  test_calculate_cell();
  test_voronoi_calculate_faces();
  test_cell_intersections();
  test_degeneracies();
}
