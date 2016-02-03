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
#ifndef SWIFT_VORONOI_H
#define SWIFT_VORONOI_H

/* Includes. */
#include "part.h"
#include <string.h>

/* Box boundary flags used to signal cells neighbouring the box boundary */
#define BOX_FRONT   0xfffffff0
#define BOX_BACK    0xfffffff1
#define BOX_TOP     0xfffffff2
#define BOX_BOTTOM  0xfffffff3
#define BOX_LEFT    0xfffffff4
#define BOX_RIGHT   0xfffffff5

/* Maximal number of vertices allowed for a Voronoi cell during the cell
   construction */
#define MAX_NUM_VERT 300
#define MAX_NUM_EDGE (3*MAX_NUM_VERT)

#define VORONOI_TOLERANCE 1.e-7

/* Internal representation of a Voronoi cell, which has more memory to store
   intermediate vertices, and which uses a more efficient memory layout */
struct voronoi_cell {

  /* Number of vertices */
  int nvert;

  /* Vertex coordinates */
  float vertices[3*MAX_NUM_VERT];

  /* Number of edges for every vertex */
  int orders[MAX_NUM_VERT];

  /* Offsets of the edges, edgeindices and neighbours corresponding to a
     particular vertex in the internal arrays */
  int offsets[MAX_NUM_VERT];

  /* Edge information */
  int edges[MAX_NUM_EDGE];

  /* Additional edge information */
  int edgeindices[MAX_NUM_EDGE];

  /* Neighbour information */
  unsigned long long ngbs[MAX_NUM_EDGE];

};

/**
 * @brief Method that loads geometry information from a #part into a
 *  #voronoi_cell
 */
__attribute__((always_inline)) INLINE static void voronoi_set_cell_values(
  struct part *p, struct voronoi_cell *c){

  int i;
  int numedge;

  c->nvert = p->geometry.nvert;
  for(i = 0; i < c->nvert; i++){
    if(i){
      c->offsets[i] = c->offsets[i-1] + p->geometry.orders[i-1];
    } else {
      c->offsets[i] = 0;
    }
  }
  numedge = c->offsets[i-1] + p->geometry.orders[i-1];

  memcpy(p->geometry.vertices, c->vertices, 3*c->nvert*sizeof(float));
  memcpy(p->geometry.orders, c->orders, c->nvert*sizeof(int));
  memcpy(p->geometry.edges, c->edges, numedge*sizeof(int));
  memcpy(p->geometry.edgeindices, c->edgeindices, numedge*sizeof(int));
  memcpy(p->geometry.ngbs, c->ngbs, numedge*sizeof(unsigned long long));

}

/**
 * @brief Method that loads a #voronoi_cell into the geometry information of a
 *  #part
 *
 * Since a #voronoi_cell can store more information than fits in a #part, we
 * need to check if this is possible.
 */
__attribute__((always_inline)) INLINE static void voronoi_set_particle_values(
  struct part *p, struct voronoi_cell *c){

  if(c->nvert > VORONOI_MAXVERT){
    error("Too many vertices!");
  }

  int i;
  int numedge;

  numedge = c->offsets[c->nvert-1]+c->orders[c->nvert-1];
  if(numedge > VORONOI_MAXEDGE){
    error("Too many edges!");
  }

  p->geometry.nvert = c->nvert;
  for(i = 0; i < c->nvert; i++){
    p->geometry.vertices[i] = c->vertices[i];
  }

  memcpy(c->orders, p->geometry.orders, c->nvert*sizeof(int));
  memcpy(c->edges, p->geometry.edges, numedge*sizeof(int));
  memcpy(c->edgeindices, p->geometry.edgeindices, numedge*sizeof(int));
  memcpy(c->ngbs, p->geometry.ngbs, numedge*sizeof(unsigned long long));

}

/**
 * @brief Get the index of the vertex pointed to by the given edge of the given
 *  vertex
 */
__attribute__((always_inline)) INLINE static int voronoi_get_edge(
  struct voronoi_cell *c, int vertex, int edge){
  return c->edges[ c->offsets[vertex] + edge ];
}

/**
 * @brief Get the index of vertex in the edge list of the vertex pointed to by
 *  the given edge of the given vertex
 */
__attribute__((always_inline)) INLINE static int voronoi_get_edgeindex(
  struct voronoi_cell *c, int vertex, int edge){
  return c->edgeindices[ c->offsets[vertex] + edge ];
}

/**
 * @brief Set the index of the vertex pointed to by the given edge of the given
 *  vertex
 */
__attribute__((always_inline)) INLINE static void voronoi_set_edge(
  struct voronoi_cell *c, int vertex, int edge, int value){
  c->edges[ c->offsets[vertex] + edge ] = value;
}

/**
 * @brief Set the index of vertex in the edge list of the vertex pointed to by
 *  the given edge of the given vertex
 */
__attribute__((always_inline)) INLINE static void voronoi_set_edgeindex(
  struct voronoi_cell *c, int vertex, int edge, int value){
  c->edgeindices[ c->offsets[vertex] + edge ] = value;
}

/**
 * @brief Get the neighbour for the given edge of the given vertex
 */
__attribute__((always_inline)) INLINE static int voronoi_get_ngb(
  struct voronoi_cell *c, int vertex, int edge){
  return c->ngbs[ c->offsets[vertex] + edge ];
}

/**
 * @brief Set the neighbour for the given edge of the given vertex
 */
__attribute__((always_inline)) INLINE static void voronoi_set_ngb(
  struct voronoi_cell *c, int vertex, int edge, int value){
  c->ngbs[ c->offsets[vertex] + edge ] = value;
}

/**
 * @brief Check if the given vertex is above, below or on the cutting plane
 */
__attribute__((always_inline)) INLINE static int voronoi_test_vertex(
    float *v, float *dx, float r2, float *test){

    *test = v[0]*dx[0] + v[1]*dx[1] + v[2]*dx[2] - r2;
    if(*test < -VORONOI_TOLERANCE){
        return -1;
    }
    if(*test > VORONOI_TOLERANCE){
        return 1;
    }
    return 0;

}

/**
 * @brief Initialize the cell as a cube with side 2*h, centered on the particle
 *  position
 */
__attribute__((always_inline)) INLINE static void voronoi_initialize(
    struct part *p){

  p->geometry.nvert = 8;

  double *origin = p->x;

  // (0, 0, 0) -- 0
  p->geometry.vertices[ 0] = 0. - origin[0];
  p->geometry.vertices[ 1] = 0. - origin[1];
  p->geometry.vertices[ 2] = 0. - origin[2];

  // (0, 0, 1)-- 1
  p->geometry.vertices[ 3] = 0. - origin[0];
  p->geometry.vertices[ 4] = 0. - origin[1];
  p->geometry.vertices[ 5] = 1. - origin[2];

  // (0, 1, 0) -- 2
  p->geometry.vertices[ 6] = 0. - origin[0];
  p->geometry.vertices[ 7] = 1. - origin[1];
  p->geometry.vertices[ 8] = 0. - origin[2];

  // (0, 1, 1) -- 3
  p->geometry.vertices[ 9] = 0. - origin[0];
  p->geometry.vertices[10] = 1. - origin[1];
  p->geometry.vertices[11] = 1. - origin[2];

  // (1, 0, 0) -- 4
  p->geometry.vertices[12] = 1. - origin[0];
  p->geometry.vertices[13] = 0. - origin[1];
  p->geometry.vertices[14] = 0. - origin[2];

  // (1, 0, 1) -- 5
  p->geometry.vertices[15] = 1. - origin[0];
  p->geometry.vertices[16] = 0. - origin[1];
  p->geometry.vertices[17] = 1. - origin[2];

  // (1, 1, 0) -- 6
  p->geometry.vertices[18] = 1. - origin[0];
  p->geometry.vertices[19] = 1. - origin[1];
  p->geometry.vertices[20] = 0. - origin[2];

  // (1, 1, 1) -- 7
  p->geometry.vertices[21] = 1. - origin[0];
  p->geometry.vertices[22] = 1. - origin[1];
  p->geometry.vertices[23] = 1. - origin[2];

  p->geometry.orders[0] = 3;
  p->geometry.orders[1] = 3;
  p->geometry.orders[2] = 3;
  p->geometry.orders[3] = 3;
  p->geometry.orders[4] = 3;
  p->geometry.orders[5] = 3;
  p->geometry.orders[6] = 3;
  p->geometry.orders[7] = 3;

  // edges are ordered counterclockwise w.r.t. a vector pointing from the
  // cell generator to the vertex
  // (0, 0, 0) corner
  p->geometry.edges[ 0] = 1;
  p->geometry.edges[ 1] = 2;
  p->geometry.edges[ 2] = 4;
  p->geometry.edgeindices[ 0] = 0;
  p->geometry.edgeindices[ 1] = 2;
  p->geometry.edgeindices[ 2] = 0;

  // (0, 0, 1) corner
  p->geometry.edges[ 3] = 0;
  p->geometry.edges[ 4] = 5;
  p->geometry.edges[ 5] = 3;
  p->geometry.edgeindices[ 3] = 0;
  p->geometry.edgeindices[ 4] = 2;
  p->geometry.edgeindices[ 5] = 1;

  // (0, 1, 0) corner
  p->geometry.edges[ 6] = 3;
  p->geometry.edges[ 7] = 6;
  p->geometry.edges[ 8] = 0;
  p->geometry.edgeindices[ 6] = 0;
  p->geometry.edgeindices[ 7] = 0;
  p->geometry.edgeindices[ 8] = 1;

  // (0, 1, 1) corner
  p->geometry.edges[ 9] = 2;
  p->geometry.edges[10] = 1;
  p->geometry.edges[11] = 7;
  p->geometry.edgeindices[ 9] = 0;
  p->geometry.edgeindices[10] = 2;
  p->geometry.edgeindices[11] = 0;

  // (1, 0, 0) corner
  p->geometry.edges[12] = 0;
  p->geometry.edges[13] = 6;
  p->geometry.edges[14] = 5;
  p->geometry.edgeindices[12] = 2;
  p->geometry.edgeindices[13] = 2;
  p->geometry.edgeindices[14] = 0;

  // (1, 0, 1) corner
  p->geometry.edges[15] = 4;
  p->geometry.edges[16] = 7;
  p->geometry.edges[17] = 1;
  p->geometry.edgeindices[15] = 2;
  p->geometry.edgeindices[16] = 1;
  p->geometry.edgeindices[17] = 1;

  // (1, 1, 0) corner
  p->geometry.edges[18] = 2;
  p->geometry.edges[19] = 7;
  p->geometry.edges[20] = 4;
  p->geometry.edgeindices[18] = 1;
  p->geometry.edgeindices[19] = 2;
  p->geometry.edgeindices[20] = 1;

  // (1, 1, 1) corner
  p->geometry.edges[21] = 3;
  p->geometry.edges[22] = 5;
  p->geometry.edges[23] = 6;
  p->geometry.edgeindices[21] = 2;
  p->geometry.edgeindices[22] = 1;
  p->geometry.edgeindices[23] = 1;

  // ngbs[3*i+j] is the neighbour corresponding to the plane clockwise of
  // edge j of vertex i
  // we do not need to initialize these fields...
  // however, we set them to a ridiculously large value to be able to track
  // faces without neighbour
  p->geometry.ngbs[ 0] = BOX_LEFT;   // (000) - (001)
  p->geometry.ngbs[ 1] = BOX_BOTTOM; // (000) - (010)
  p->geometry.ngbs[ 2] = BOX_FRONT;  // (000) - (100)

  p->geometry.ngbs[ 3] = BOX_FRONT;  // (001) - (000)
  p->geometry.ngbs[ 4] = BOX_TOP;    // (001) - (101)
  p->geometry.ngbs[ 5] = BOX_LEFT;   // (001) - (011)

  p->geometry.ngbs[ 6] = BOX_BACK;   // (010) - (011)
  p->geometry.ngbs[ 7] = BOX_BOTTOM; // (010) - (110)
  p->geometry.ngbs[ 8] = BOX_LEFT;   // (010) - (000)

  p->geometry.ngbs[ 9] = BOX_LEFT;   // (011) - (010)
  p->geometry.ngbs[10] = BOX_TOP;    // (011) - (001)
  p->geometry.ngbs[11] = BOX_BACK;   // (011) - (111)

  p->geometry.ngbs[12] = BOX_BOTTOM; // (100) - (000)
  p->geometry.ngbs[13] = BOX_RIGHT;  // (100) - (110)
  p->geometry.ngbs[14] = BOX_FRONT;  // (100) - (101)

  p->geometry.ngbs[15] = BOX_RIGHT;  // (101) - (100)
  p->geometry.ngbs[16] = BOX_TOP;    // (101) - (111)
  p->geometry.ngbs[17] = BOX_FRONT;  // (101) - (001)

  p->geometry.ngbs[18] = BOX_BACK;   // (110) - (010)
  p->geometry.ngbs[19] = BOX_RIGHT;  // (110) - (111)
  p->geometry.ngbs[20] = BOX_BOTTOM; // (110) - (100)

  p->geometry.ngbs[21] = BOX_TOP;    // (111) - (011)
  p->geometry.ngbs[22] = BOX_RIGHT;  // (111) - (101)
  p->geometry.ngbs[23] = BOX_BACK;   // (111) - (110)
}

/**
 * @brief Intersect particle pi with particle pj and adapt its Voronoi cell
 *  structure
 *
 * odx = x_i - x_j!!!
 */

__attribute__((always_inline)) INLINE static void voronoi_intersect(
    float *odx, struct part *pi, struct part *pj) {

  float dx[3];
  float r2;
  float u, l, q;
  int up, us, uw, lp, ls, lw, qp, qs, qw;
  int complicated;

  struct voronoi_cell c;
  voronoi_set_cell_values(pi, &c);

  dx[0] = -0.5f * odx[0];
  dx[1] = -0.5f * odx[1];
  dx[2] = -0.5f * odx[2];
  r2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];

  uw = voronoi_test_vertex(&c.vertices[0], dx, r2, &u);
  up = 0;
  complicated = 0;
  if(uw == 0){

    complicated = 1;

  } else {

    if(uw == 1){

      lp = voronoi_get_edge(&c, up, 0);
      lw = voronoi_test_vertex(&c.vertices[3*lp], dx, r2, &l);
      us = 1;
      while(us < c.orders[up] && l >= u){
        lp = voronoi_get_edge(&c, up, us);
        lw = voronoi_test_vertex(&c.vertices[3*lp], dx, r2, &l);
        us++;
      }
      us--;
      if(l >= u){
        error("Cell completely gone! This should not happen.");
      }
      ls = voronoi_get_edgeindex(&c, up, us);

      while(lw == 1){
        u = l;
        up = lp;
        us = 0;
        while(us < ls && l >= u){
          lp = voronoi_get_edge(&c, up, us);
          lw = voronoi_test_vertex(&c.vertices[3*lp], dx, r2, &l);
          us++;
        }
        if(l >= u){
          us++;
          while(us < c.orders[up] && l >= u){
            lp = voronoi_get_edge(&c, up, us);
            lw = voronoi_test_vertex(&c.vertices[3*lp], dx, r2, &l);
            us++;
          }
          if(l >= u){
            error("Cell completely gone! This should not happen.");
          }
        }
        us--;
        ls = voronoi_get_edgeindex(&c, up, us);
      }
      if(lw == 0){
        up = lp;
        complicated = 1;
      }

    } else { /* if(uw == 1) */

      qp = voronoi_get_edge(&c, up, 0);
      qw = voronoi_test_vertex(&c.vertices[3*qp], dx, r2, &q);
      us = 1;
      while(us < c.orders[up] && u >= q){
        qp = voronoi_get_edge(&c, up, us);
        qw = voronoi_test_vertex(&c.vertices[3*qp], dx, r2, &q);
        us++;
      }
      if(u >= q){
        /* cell unaltered */
        return;
      } else {
        us--;
      }

      while(qw == -1){
        qs = voronoi_get_edgeindex(&c, up, us);
        u = q;
        up = qp;
        us = 0;
        while(us < qs && u >= q){
          qp = voronoi_get_edge(&c, up, us);
          qw = voronoi_test_vertex(&c.vertices[3*qp], dx, r2, &q);
          us++;
        }
        if(u >= q){
          us++;
          while(us < c.orders[up] && u >= q){
            qp = voronoi_get_edge(&c, up, us);
            qw = voronoi_test_vertex(&c.vertices[3*qp], dx, r2, &q);
            us++;
          }
          if(u >= q){
            /* cell unaltered */
            return;
          }
        }
        us--;
      }
      if(qw == 1){
          lp = up;
          ls = us;
          l = u;
          up = qp;
          us = voronoi_get_edgeindex(&c, lp, ls);
          u = q;
      } else {
          up = qp;
          complicated = 1;
      }

    } /* if(uw == 1) */

  } /* if(uw == 0) */

  int vindex;
  int visitflags[MAX_NUM_VERT];
  int dstack[MAX_NUM_VERT];
/*  int dstack_size = 1;*/
  float r;
  int cs, rp;
  int double_edge = 0;
  int i, j, k;

  if(complicated){

    lp = voronoi_get_edge(&c, up, 0);
    lw = voronoi_test_vertex(&c.vertices[3*lp], dx, r2, &l);

    if(lw != -1){

      rp = lw;
      i = 1;
      lp = voronoi_get_edge(&c, up, i);
      lw = voronoi_test_vertex(&c.vertices[3*lp], dx, r2, &l);
      while(lw != -1){
        i++;
        if(i == c.orders[up]){
          error("Cell completely gone! This should not happen.");
        }
        lp = voronoi_get_edge(&c, up, i);
        lw = voronoi_test_vertex(&c.vertices[3*lp], dx, r2, &l);
      }

      j = i+1;
      while(j < c.orders[up]){
        lp = voronoi_get_edge(&c, up, j);
        lw = voronoi_test_vertex(&c.vertices[3*lp], dx, r2, &l);
        if(lw != -1){
          break;
        }
        j++;
      }

      if(j == c.orders[up] && i == 1 && rp == 0){
        k = c.orders[up];
        double_edge = 1;
      } else {
        k = j - i + 2;
      }

      /* create new order k vertex */
      vindex = c.nvert;
      c.nvert++;
      if(c.nvert == MAX_NUM_VERT){
        error("Too many vertices!");
      }
      c.orders[vindex] = k;
      c.offsets[vindex] = c.offsets[vindex-1] + c.orders[vindex-1];
      if(c.offsets[vindex] + k >= MAX_NUM_EDGE){
        error("Too many edges!");
      }

      k = 1;

      visitflags[vindex] = 0;
      c.vertices[3*vindex+0] = c.vertices[3*up+0];
      c.vertices[3*vindex+1] = c.vertices[3*up+1];
      c.vertices[3*vindex+2] = c.vertices[3*up+2];

      us = i-1;
      if(i < 0){
        i = c.orders[up]-1;
      }
      while(i < j){
          qp = voronoi_get_edge(&c, up, i);
          qs = voronoi_get_edgeindex(&c, up, i);
          voronoi_set_ngb(&c, vindex, k, voronoi_get_ngb(&c, up, i));
          voronoi_set_edge(&c, vindex, k, qp);
          voronoi_set_edgeindex(&c, vindex, k, qs);
          voronoi_set_edge(&c, qp, qs, vindex);
          voronoi_set_edgeindex(&c, qp, qs, k);
          voronoi_set_edge(&c, up, i, -1);
          i++;
          k++;
      }
      if(i == c.orders[up]){
          qs = 0;
      } else {
          qs = i;
      }

    } else { /* if(lw != -1) */

      i = c.orders[up]-1;
      lp = voronoi_get_edge(&c, up, i);
      lw = voronoi_test_vertex(&c.vertices[3*lp], dx, r2, &l);
      while(lw == -1){
        i--;
        if(i == 0){
          /* cell unaltered */
          return;
        }
        lp = voronoi_get_edge(&c, up, i);
        lw = voronoi_test_vertex(&c.vertices[3*lp], dx, r2, &l);
      }

      j = 1;
      qp = voronoi_get_edge(&c, up, j);
      qw = voronoi_test_vertex(&c.vertices[3*qp], dx, r2, &q);
      while(qw == -1){
        j++;
        qp = voronoi_get_edge(&c, up, j);
        qw = voronoi_test_vertex(&c.vertices[3*qp], dx, r2, &l);
      }

      if(i == j && qw == 0){
        double_edge = 1;
        k = c.orders[up];
      } else {
        k = c.orders[up] - i + j + 1;
      }

      /* create new order k vertex */
      vindex = c.nvert;
      c.nvert++;
      if(c.nvert == MAX_NUM_VERT){
        error("Too many vertices!");
      }
      c.orders[vindex] = k;
      c.offsets[vindex] = c.offsets[vindex-1] + c.orders[vindex-1];
      if(c.offsets[vindex] + k >= MAX_NUM_EDGE){
        error("Too many edges!");
      }
      k = 1;

      visitflags[vindex] = 0;
      c.vertices[3*vindex+0] = c.vertices[3*up+0];
      c.vertices[3*vindex+1] = c.vertices[3*up+1];
      c.vertices[3*vindex+2] = c.vertices[3*up+2];

      us = i;
      i++;
      while(i < c.orders[up]){
        qp = voronoi_get_edge(&c, up, i);
        qs = voronoi_get_edgeindex(&c, up, i);
        voronoi_set_ngb(&c, vindex, k, voronoi_get_ngb(&c, up, i));
        voronoi_set_edge(&c, vindex, k, qp);
        voronoi_set_edgeindex(&c, vindex, k, qs);
        voronoi_set_edge(&c, qp, qs, vindex);
        voronoi_set_edgeindex(&c, qp, qs, k);
        voronoi_set_edge(&c, up, i, -1);
        i++;
        k++;
      }
      i = 0;
      while(i < j){
        qp = voronoi_get_edge(&c, up, i);
        qs = voronoi_get_edgeindex(&c, up, i);
        voronoi_set_ngb(&c, vindex, k, voronoi_get_ngb(&c, up, i));
        voronoi_set_edge(&c, vindex, k, qp);
        voronoi_set_edgeindex(&c, vindex, k, qs);
        voronoi_set_edge(&c, qp, qs, vindex);
        voronoi_set_edgeindex(&c, qp, qs, k);
        voronoi_set_edge(&c, up, i, -1);
        i++;
        k++;
      }
      qs = j;
    }

    if(!double_edge){
        voronoi_set_ngb(&c, vindex, k, voronoi_get_ngb(&c, up, qs));
        voronoi_set_ngb(&c, vindex, 0, pj->id);
    } else {
        voronoi_set_ngb(&c, vindex, 0, voronoi_get_ngb(&c, up, qs));
    }

    dstack[0] = up;

    cs = k;
    qp = up;
    q = u;
    i = voronoi_get_edge(&c, up, us);
    us = voronoi_get_edgeindex(&c, up, us);
    up = i;
    visitflags[qp] = -vindex;

  } else { /* if(complicated) */

    r = u/(u-l);
    l = 1.0f - r;

    /* create new order 3 vertex */
    vindex = c.nvert;
    c.nvert++;
    if(c.nvert == MAX_NUM_VERT){
      error("Too many vertices!");
    }
    c.orders[vindex] = 3;
    c.offsets[vindex] = c.offsets[vindex-1] + c.orders[vindex-1];
    if(c.offsets[vindex] + 3 >= MAX_NUM_EDGE){
      error("Too many edges!");
    }

    visitflags[vindex] = 0;
    c.vertices[3*vindex+0] = c.vertices[3*lp+0] * r + c.vertices[3*up+0] * l;
    c.vertices[3*vindex+1] = c.vertices[3*lp+1] * r + c.vertices[3*up+1] * l;
    c.vertices[3*vindex+2] = c.vertices[3*lp+2] * r + c.vertices[3*up+2] * l;

    dstack[0] = up;

    voronoi_set_edge(&c, vindex, 1, lp);
    voronoi_set_edgeindex(&c, vindex, 1, ls);
    voronoi_set_edge(&c, lp, ls, vindex);
    voronoi_set_edgeindex(&c, lp, ls, 1);
    voronoi_set_edge(&c, up, us, -1);

    voronoi_set_ngb(&c, vindex, 0, pj->id);
    voronoi_set_ngb(&c, vindex, 1, voronoi_get_ngb(&c, up, us));
    voronoi_set_ngb(&c, vindex, 2, voronoi_get_ngb(&c, lp, ls));

    qs = us+1;
    if(qs == c.orders[up]){
      qs = 0;
    }
    qp = up;
    q = u;

    cs = 2;

  } /* if(complicated) */

  visitflags[0]++;
  dstack[0]++;
  cs++;

#ifdef DO_NOT_COMPILE
  int cp = vindex;
  int rp = vindex;
  // cp corresponds to the last added vertex
  // rp corresponds to the first added vertex
  vindex++;
  // loop around the edges until we arrive back at the original up vertex
  // we start the search from the next neighbour vertex counterclockwise from
  // neighbour lp
  // qp always holds the last vertex above the plane
  // lp always holds the last vertex below the plane
  unsigned int count = 0;
  while(qp != up || qs != us){
  //        cerr << "qp: " << qp << ", qs: " << qs << endl;
  //        cerr << "up: " << up << ", us: " << us << endl;
      count++;
      if(count == 1000){
  //            print_cell(c);
          cerr << "Endless loop" << endl;
          my_exit();
      }
      lp = edge(c, qp, qs);
      if(lp < 0){
          cerr << "Negative edge!" << endl;
          my_exit();
      }
      lw = test_vertex(&c->vertices[3*lp], x, r2, l);
      if(lw == 0){
          LOGR("complicated case");
          if(global_logging){
              print_cell(c);
          }
          LOGD("lp: %i, lw: %i", lp, lw);
          
          // if we end up here, the next vertex lies on the plane
          // there are several possibilities:
          //  - the previous vertex was also on the plane, but is not
          //    connected through a float edge with this one
          //  - the previous vertex was also on the plane and is connected
          //    with this one through a float edge
          //  - the previous vertex was above the plane
          
          // we need to figure out how many edges the vertex gets
          // for this, we count the number of neighbours below the plane
          // in the moving direction
          // if the previous vertex had a float edge, this vertex will have
          // one neighbour less (sic voro++)
          
          int k = 1;
          if(float_edge){
              k = 0;
          }
          qs = edgeindex(c, qp, qs);
          qp = lp;
          int iqs = qs;
          LOGD("qp: %i, qs: %i", qp, qs);
          
          do{
              k++;
              qs++;
              if(qs == c->vorders[qp]){
                  qs = 0;
              }
              lp = edge(c, qp, qs);
              LOGD("lp: %i, qp: %i, qs: %i", lp, qp, qs);
              lw = test_vertex(&c->vertices[3*lp], x, r2, l);
              LOGD("pi: %g %g %g", pi->x[0], pi->x[1], pi->x[2]);
              LOGD("pj: %g %g %g", pj->x[0], pj->x[1], pj->x[2]);
          } while(lw == -1);
          
          // turns out we need a way to signal if this vertex has been visited
          // before: visitflags
          int j = visitflags[qp];
          bool new_float_edge;
          int i;
          if(qp == up && qs == us){
              LOGR("qp == up && qs == us");
              new_float_edge = false;
              if(j > 0){
                  LOGR("k += c->vorders");
                  k += c->vorders[j];
              }
          } else {
              // NOTE TO SELF: a lot of different cases to unit test...
              if(j > 0){
                  LOGR("Vertex already visited");
                  k += c->vorders[j];
                  if(lw == 0){
                      LOGR("Vertex on plane");
                      i = -visitflags[lp];
                      if(i > 0){
                          LOGR("i > 0");
                          if(edge(c, i, c->vorders[i]-1) == j){
                              LOGR("new float edge");
                              new_float_edge = true;
                              k--;
                          } else {
                              LOGR("no new float edge");
                              new_float_edge = false;
                          }
                      } else {
                          LOGR("i <= 0");
                          if(j == rp && lp == up && edge(c, qp, qs) == us){
                              LOGR("new float edge");
                              new_float_edge = true;
                              k--;
                          } else {
                              LOGR("no new float edge");
                              new_float_edge = false;
                          }
                      }
                  } else {
                      LOGR("Vertex not on plane");
                      new_float_edge = false;
                  }
              } else {
                  LOGR("Vertex not yet visited");
                  if(lw == 0){
                      LOGR("Vertex on plane");
                      i = -visitflags[lp];
                      if(i == cp){
                          LOGR("new float edge");
                          new_float_edge = true;
                          k--;
                      } else {
                          LOGR("no new float edge");
                          new_float_edge = false;
                      }
                  } else {
                      LOGR("Vertex not on plane");
                      new_float_edge = false;
                  }
              }
          }
          
          LOGD("k: %i, new_float_edge: %i", k, new_float_edge);
          
          // create new vertex with order k
          // no idea why it matters if the vertex has already been visited...
          // oh no, I see: if the vertex has already been visited, it is
          // already a new vertex. And we don't want to add new vertices
          // twice, do we?
          // (actually, in our paradigm, we do)
          // problem: we will have to figure everything out ourselves...
          // consult the voro++ source code to implement this case
          if(j > 0){
              cerr << "You have ended up in a very special situation, that is"
                   << "not handled by this code.\nPlease contact your local"
                   << "code developer, tell him exactly what you did to get"
                   << "here, and ask him/her to solve this issue.\nThank you!"
                   << endl;
              my_exit();
          }
          
          LOGR("Vertex not yet visited");
          // allocate a new vertex
          vindex = grow_arrays(c, k);
          
          int *mvisitflags = (int*) realloc(visitflags, c->nvert*sizeof(int));
          visitflags = mvisitflags;
          visitflags[vindex] = 0;
          
          c->vertices[3*vindex+0] = c->vertices[3*qp+0];
          c->vertices[3*vindex+1] = c->vertices[3*qp+1];
          c->vertices[3*vindex+2] = c->vertices[3*qp+2];
  //                cerr << "Creating vertex" << endl;
  //                cerr << c->vertices[3*vindex] << "\t" << c->vertices[3*vindex+1] << "\t" << c->vertices[3*vindex+2] << endl;
          visitflags[qp] = -vindex;
          dstack_size++;
          int *more_dstack = (int*) realloc(dstack, dstack_size*sizeof(int));
          dstack = more_dstack;
          dstack[dstack_size-1] = qp;
          j = vindex;
          i = 0;
          
          if(!float_edge){
              LOGR("No float edge");
              set_ngb(c, j, i, pj->id);
              set_edge(c, j, i, cp);
              set_edgeindex(c, j, i, cs);
              set_edge(c, cp, cs, j);
              set_edgeindex(c, cp, cs, i);
              i++;
          }
          
          // copy the edges of the old vertex (one less if it was a new float edge)
          qs = iqs;
          int lim = k-1;
          if(new_float_edge){
              lim = k;
          }
          while(i < lim){
              qs++;
              if(qs == c->vorders[qp]){
                  qs = 0;
              }
              lp = edge(c, qp, qs);
              ls = edgeindex(c, qp, qs);
              set_ngb(c, j, i, ngb(c, qp, qs));
              set_edge(c, j, i, lp);
              set_edgeindex(c, j, i, ls);
              set_edge(c, lp, ls, j);
              set_edgeindex(c, lp, ls, i);
              set_edge(c, qp, qs, -1);
              i++;
          }
          qs++;
          if(qs == c->vorders[qp]){
              qs = 0;
          }
          cs = i;
          cp = j;
          
          if(new_float_edge){
              set_ngb(c, j, 0, ngb(c, qp, qs));
          } else {
              set_ngb(c, j, cs, ngb(c, qp, qs));
          }
          
          float_edge = new_float_edge;
          if(global_logging){
              print_cell(c);
  //                my_exit();
          }
      } else {
          LOGR("Normal case");
          if(lw == 1){
              LOGR("Above plane");
  //            cout << "Vertex above plane" << endl;
              // delete lp, it is still above the plane
              // continue with the next vertex of lp
              qs = edgeindex(c, qp, qs) + 1;
              if(qs == c->vorders[lp]){
                  qs = 0;
              }
              qp = lp;
              q = l;
              dstack_size++;
              int *more_dstack = (int*) realloc(dstack, dstack_size*sizeof(int));
              dstack = more_dstack;
              dstack[dstack_size-1] = qp;
          } else {
              LOGR("normal vertex creation");
  //            cout << "Vertex below plane" << endl;
              // since the previous vertex was above the plane, we have found an
              // intersecting edge
              // add a new vertex
              r = q/(q-l);
              l = 1. - r;

              vindex = grow_arrays(c, 3);

              c->vertices[3*vindex+0] = c->vertices[3*lp+0]*r +
                                        c->vertices[3*qp+0]*l;
              c->vertices[3*vindex+1] = c->vertices[3*lp+1]*r +
                                        c->vertices[3*qp+1]*l;
              c->vertices[3*vindex+2] = c->vertices[3*lp+2]*r +
                                        c->vertices[3*qp+2]*l;

              // get the index of the edge lp in qp
              ls = edgeindex(c, qp, qs);
              set_edge(c, vindex, 0, cp);
              set_edge(c, vindex, 1, lp);
              set_edgeindex(c, vindex, 0, cs);
              set_edgeindex(c, vindex, 1, ls);
              set_edge(c, lp, ls, vindex);
              set_edgeindex(c, lp, ls, 1);
              set_edge(c, cp, cs, vindex);
              set_edgeindex(c, cp, cs, 0);
              set_edge(c, qp, qs, -1);
              // ngbs
              set_ngb(c, vindex, 0, pj->id);
              set_ngb(c, vindex, 1, ngb(c, qp, qs));
              set_ngb(c, vindex, 2, ngb(c, lp, ls));
              // continue with the next edge of this vertex
              qs = qs+1;
              if(qs == c->vorders[qp]){
                  qs = 0;
              }
              cp = vindex;
              vindex++;
              cs = 2;
          }
      }
  }

  set_edge(c, cp, cs, rp);
  set_edge(c, rp, 0, cp);
  set_edgeindex(c, cp, cs, 0);
  set_edgeindex(c, rp, 0, cs);

  if(global_logging){
      print_cell(c);
  }
  //    cout << "delete stack:" << endl;
  //    for(int i = 0; i < dstack_size; i++){
  //        cout << dstack[i] << endl;
  //    }

  // add vertices connected to deleted ones to the delete stack
  // using dstack_size (which changes inside the loop is not a mistake:
  //  the newly deleted vertices can have connections that need to be
  //  deleted too, and this way, we make sure they are deleted as well
  //  this is also the reason why we reset the edges here, since otherwise
  //  we might get stuck in an endless loop
  // other interesting thought: we duplicate some vertices... We cannot simply
  // remove the vertices connected to duplicated vertices, as these might be
  // important. Or aren't they... Not sure...
  for(int i = 0; i < dstack_size; i++){
      for(int j = 0; j < c->vorders[dstack[i]]; j++){
          if(edge(c, dstack[i], j) >= 0){
              dstack_size++;
              int *more_dstack = (int*) realloc(dstack, dstack_size*sizeof(int));
              dstack = more_dstack;
              dstack[dstack_size-1] = edge(c, dstack[i], j);
              set_edge(c, dstack[i], j, -1);
              set_edgeindex(c, dstack[i], j, -1);
          }
      }
  }

  maxnumdel = max(dstack_size, maxnumdel);

  //    cout << "edges after:" << endl;
  //    for(int i = 0; i < 6*c->nvert; i++){
  //        cout << i << ": " << c->edges[i] << endl;
  //    }

  for(int i = 0; i < c->nvert; i++){
      for(int j = 0; j < c->vorders[i]; j++){
          int k = edge(c, i, j);
          if(k >= 0 && edge(c, k, 0) < 0){
              cerr << "Inconsistency!" << endl;
              cerr << "i: " << i << ", j: " << j << endl;
              print_cell(c);
              my_exit();
          }
      }
  }

  //    cout << "edges before:" << endl;
  //    for(int i = 0; i < 6*c->nvert; i++){
  //        cout << i << ": " << c->edges[i] << endl;
  //    }

  // loop through the edges. If a deactivated edge is detected, remove the
  // corresponding vertex by moving the next vertex (if it exists) to the
  // current position. Also move its edges and change the value at the
  // other endpoint of its edges
  // the way we do it now is a bit tricky, since copying in place is maybe
  // not safe if not all vertices have the same number of edges...
  // better make new arrays
  int *new_vorders = (int*) malloc(c->nvert*sizeof(int));
  float *new_vertices = (float*) malloc(3*c->nvert*sizeof(float));
  int *new_edgeoffsets = (int*) malloc(c->nvert*sizeof(int));
  int *new_edges = (int*) malloc((c->edgeoffsets[c->nvert-1] + 2*c->vorders[c->nvert-1])*sizeof(int));
  int *new_ngboffsets = (int*) malloc(c->nvert*sizeof(int));
  unsigned long long *new_ngbs = (unsigned long long*) malloc((c->ngboffsets[c->nvert-1]+c->vorders[c->nvert-1])*sizeof(unsigned long long));
  int newnvert;
  for(int i = 0; i < c->nvert; i++){
      if(edge(c, i, 0) < 0){
          // find the next valid vertex
          int j = i+1;
          while(j < c->nvert && edge(c, j, 0) < 0){
              j++;
          }
          if(j == c->nvert){
              // no more valid other vertex found, we can stop
              newnvert = i;
              break;
          }
          // copy the vertex position
          new_vertices[3*i+0] = c->vertices[3*j+0];
          new_vertices[3*i+1] = c->vertices[3*j+1];
          new_vertices[3*i+2] = c->vertices[3*j+2];
          
          // copy the vorder
          new_vorders[i] = c->vorders[j];
          
          // copy the edgeoffset
          if(i){
              new_edgeoffsets[i] = new_edgeoffsets[i-1] + 2*new_vorders[i-1];
          } else {
              new_edgeoffsets[i] = 0;
          }
          
          // copy the edges
          for(int k = 0; k < c->vorders[j]; k++){
              new_edges[new_edgeoffsets[i]+k] = edge(c, j, k);
              new_edges[new_edgeoffsets[i]+new_vorders[i]+k] = edgeindex(c, j, k);
          }
          
          // copy the ngbs
          if(i){
              new_ngboffsets[i] = new_ngboffsets[i-1] + new_vorders[i-1];
          } else {
              new_ngboffsets[i] = 0;
          }
          for(int k = 0; k < c->vorders[j]; k++){
              new_ngbs[new_ngboffsets[i]+k] = ngb(c, j, k);
          }
          
          for(int k = 0; k < c->vorders[j]; k++){
              int m = edge(c, j, k);
              int n = edgeindex(c, j, k);
              // < j would also work...
              if(m < i){
                  new_edges[new_edgeoffsets[m]+n] = i;
              } else {
                  set_edge(c, m, n, i);
              }
          }
          
          // deactivate the old edges
          set_edge(c, j, 0, -1);
      } else {
          // simply copy the original values
          new_vorders[i] = c->vorders[i];
          new_vertices[3*i+0] = c->vertices[3*i+0];
          new_vertices[3*i+1] = c->vertices[3*i+1];
          new_vertices[3*i+2] = c->vertices[3*i+2];
          new_edgeoffsets[i] = c->edgeoffsets[i];
          for(int j = 0; j < c->vorders[i]; j++){
              new_edges[c->edgeoffsets[i]+j] = edge(c, i, j);
              new_edges[c->edgeoffsets[i]+c->vorders[i]+j] = edgeindex(c, i, j);
          }
          new_ngboffsets[i] = c->ngboffsets[i];
          for(int j = 0; j < c->vorders[i]; j++){
              new_ngbs[c->ngboffsets[i]+j] = ngb(c, i, j);
          }
      }
  }

  free(c->vorders);
  free(c->vertices);
  free(c->edgeoffsets);
  free(c->edges);
  free(c->ngboffsets);
  free(c->ngbs);
  c->nvert = newnvert;
  int *mvorder = (int*) realloc(new_vorders, c->nvert*sizeof(int));
  c->vorders = mvorder;
  float *mvert = (float*) realloc(new_vertices, 3*c->nvert*sizeof(float));
  c->vertices = mvert;
  int *medgeof = (int*) realloc(new_edgeoffsets, c->nvert*sizeof(int));
  c->edgeoffsets = medgeof;
  int *medge = (int*) realloc(new_edges, (c->edgeoffsets[c->nvert-1] + 2*c->vorders[c->nvert-1])*sizeof(int));
  c->edges = medge;
  int *mngbof = (int*) realloc(new_ngboffsets, c->nvert*sizeof(int));
  c->ngboffsets = mngbof;
  unsigned long long *mngb = (unsigned long long*) realloc(new_ngbs, (c->ngboffsets[c->nvert-1]+c->vorders[c->nvert-1])*sizeof(unsigned long long));
  c->ngbs = mngb;

  maxnumvert = max(c->nvert, maxnumvert);
  maxnumedge = max(6*c->nvert, maxnumedge);

  //    for(int i = 0; i < c->nvert; i++){
  //        for(int j = 0; j < c->vorders[i]; j++){
  //            for(int k = 0; k < c->vorders[i]; k++){
  //                if(j != k){
  //                    if(c->edges[c->edgeoffsets[i]+j] == c->edges[c->edgeoffsets[i]+k]){
  //                        cerr << "float edge!" << endl;
  //                        print_cell(c);
  //                        my_exit();
  //                    }
  //                }
  //            }
  //        }
  //    }

  if(global_logging){
      print_cell(c);
  }

  // need to implement this!!
  //    collapse_order2(c);

  //    cout << "edges after:" << endl;
  //    for(int i = 0; i < 6*c->nvert; i++){
  //        cout << i << ": " << c->edges[i] << endl;
  //    }
  //    exit(1);

  //    for(int i = 0; i < c->nvert; i++){
  //        cout << i << ":" << endl;
  //        cout << c->vertices[3*i] << "\t" << c->vertices[3*i+1] << "\t"
  //             << c->vertices[3*i+2] << endl;
  //    }
  //    exit(1);

  free(visitflags);

  set_particle(pi, c);
#endif // DO_NOT_COMPILE

}

__attribute__((always_inline)) INLINE static float voronoi_volume_tetrahedron(
    float *v1, float *v2, float *v3, float *v4){

    float V;
    float r1[3], r2[3], r3[3];

    r1[0] = v2[0] - v1[0];
    r1[1] = v2[1] - v1[1];
    r1[2] = v2[2] - v1[2];
    r2[0] = v3[0] - v1[0];
    r2[1] = v3[1] - v1[1];
    r2[2] = v3[2] - v1[2];
    r3[0] = v4[0] - v1[0];
    r3[1] = v4[1] - v1[1];
    r3[2] = v4[2] - v1[2];
    V = fabs(r1[0]*r2[1]*r3[2] + r1[1]*r2[2]*r3[0] + r1[2]*r2[0]*r3[1] -
             r1[2]*r2[1]*r3[0] - r2[2]*r3[1]*r1[0] - r3[2]*r1[1]*r2[0]);
    V /= 6.;
    return V;

}

__attribute__((always_inline)) INLINE static void voronoi_centroid_tetrahedron(
    float *centroid, float *v1, float *v2, float *v3, float *v4){

    centroid[0] = 0.25f*(v1[0] + v2[0] + v3[0] + v4[0]);
    centroid[1] = 0.25f*(v1[1] + v2[1] + v3[1] + v4[1]);
    centroid[2] = 0.25f*(v1[2] + v2[2] + v3[2] + v4[2]);

}

__attribute__((always_inline)) INLINE static void voronoi_calculate_cell(
    struct part *p){

    float v1[3], v2[3], v3[3], v4[3];
    int i, j, k, l, m, n;
    float tcentroid[3];
    float tvol;

    /* we need to calculate the volume of the tetrahedra formed by the first
       vertex and the triangles that make up the other faces
       since we do not store faces explicitly, this means keeping track of the
       edges that have been processed somehow
       we follow the method used in voro++ and "flip" processed edges to
       negative values
       this also means that we need to process all triangles corresponding to
       an edge at once */
    p->voronoi.volume = 0.0f;
    v1[0] = p->voronoi.vertices[0];
    v1[1] = p->voronoi.vertices[1];
    v1[2] = p->voronoi.vertices[2];
    p->voronoi.centroid[0] = 0.0f;
    p->voronoi.centroid[1] = 0.0f;
    p->voronoi.centroid[2] = 0.0f;
    
    /* loop over all vertices (except the first one) */
    for(i = 1; i < p->voronoi.nvert; i++){
        v2[0] = p->voronoi.vertices[3*i];
        v2[1] = p->voronoi.vertices[3*i+1];
        v2[2] = p->voronoi.vertices[3*i+2];
        
        /*  loop over the edges of the vertex*/
        for(j = 0; j < 3; j++){
            k = p->voronoi.edges[6*i+j];
            /* check if the edge has already been processed */
            if(k >= 0){
                /* mark the edge as processed */
                p->voronoi.edges[6*i+j] = -k-1;
                
                /* do some magic
                   code below brainlessly copied from voro++ */
                l = p->voronoi.edges[6*i+3+j];
                if(l == 2){
                    l = 0;
                } else {
                    l++;
                }
                v3[0] = p->voronoi.vertices[3*k];
                v3[1] = p->voronoi.vertices[3*k+1];
                v3[2] = p->voronoi.vertices[3*k+2];
                m = p->voronoi.edges[6*k+l];
                p->voronoi.edges[6*k+l] = -1-m;
                while(m != i){
                    n = p->voronoi.edges[6*k+3+l];
                    if(n == 2){
                        n = 0;
                    } else {
                        n++;
                    }
                    v4[0] = p->voronoi.vertices[3*m];
                    v4[1] = p->voronoi.vertices[3*m+1];
                    v4[2] = p->voronoi.vertices[3*m+2];
                    tvol = voronoi_volume_tetrahedron(v1, v2, v3, v4);
                    p->voronoi.volume += tvol;
                    voronoi_centroid_tetrahedron(tcentroid, v1, v2, v3, v4);
                    p->voronoi.centroid[0] += tcentroid[0]*tvol;
                    p->voronoi.centroid[1] += tcentroid[1]*tvol;
                    p->voronoi.centroid[2] += tcentroid[2]*tvol;
                    k = m;
                    l = n;
                    v3[0] = v4[0];
                    v3[1] = v4[1];
                    v3[2] = v4[2];
                    m = p->voronoi.edges[6*k+l];
                    p->voronoi.edges[6*k+l] = -1-m;
                }
            }
        }
    }
    
    p->voronoi.centroid[0] /= p->voronoi.volume;
    p->voronoi.centroid[1] /= p->voronoi.volume;
    p->voronoi.centroid[2] /= p->voronoi.volume;
    
    /* centroid was calculated relative w.r.t. particle position */
    p->voronoi.centroid[0] += p->x[0];
    p->voronoi.centroid[1] += p->x[1];
    p->voronoi.centroid[2] += p->x[2];
    
    // unmark edges
    for(i = 0; i < p->voronoi.nvert; i++){
        for(j = 0; j < 3; j++){
            if(p->voronoi.edges[6*i+j] < 0){
                p->voronoi.edges[6*i+j] = -p->voronoi.edges[6*i+j]-1;
            } else {
                error("edge inconsistency");
            }
        }
    }

}

__attribute__((always_inline)) INLINE static void voronoi_calculate_faces(
    struct part *p){

    unsigned long long newngbs[300];
    int i, j, k, l, m, n;
    float area;
    float midpoint[3];
    float u[3], v[3], w[3];
    float loc_area;

    p->voronoi.nface = 0;
    for(i = 0; i < p->voronoi.nvert; i++){
        for(j = 0; j < 3; j++){
            k = p->voronoi.edges[6*i+j];
            if(k >= 0){
                newngbs[p->voronoi.nface] = p->voronoi.ngbs[3*i+j];
                area = 0.;
                midpoint[0] = 0.;
                midpoint[1] = 0.;
                midpoint[2] = 0.;
                p->voronoi.edges[6*i+j] = -1 - k;
                l = p->voronoi.edges[6*i+3+j] + 1;
                if(l == 3){
                    l = 0;
                }
                m = p->voronoi.edges[6*k+l];
                p->voronoi.edges[6*k+l] = -1 - m;
                while(m != i){
                    n = p->voronoi.edges[6*k+3+l] + 1;
                    if(n == 3){
                        n = 0;
                    }
                    u[0] = p->voronoi.vertices[3*k+0] - p->voronoi.vertices[3*i+0];
                    u[1] = p->voronoi.vertices[3*k+1] - p->voronoi.vertices[3*i+1];
                    u[2] = p->voronoi.vertices[3*k+2] - p->voronoi.vertices[3*i+2];
                    v[0] = p->voronoi.vertices[3*m+0] - p->voronoi.vertices[3*i+0];
                    v[1] = p->voronoi.vertices[3*m+1] - p->voronoi.vertices[3*i+1];
                    v[2] = p->voronoi.vertices[3*m+2] - p->voronoi.vertices[3*i+2];
                    w[0] = u[1]*v[2] - u[2]*v[1];
                    w[1] = u[2]*v[0] - u[0]*v[2];
                    w[2] = u[0]*v[1] - u[1]*v[0];
                    loc_area = sqrtf(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);
                    area += loc_area;
                    midpoint[0] += loc_area*(p->voronoi.vertices[3*k+0] +
                                             p->voronoi.vertices[3*i+0] +
                                             p->voronoi.vertices[3*m+0]);
                    midpoint[1] += loc_area*(p->voronoi.vertices[3*k+1] +
                                             p->voronoi.vertices[3*i+1] +
                                             p->voronoi.vertices[3*m+1]);
                    midpoint[2] += loc_area*(p->voronoi.vertices[3*k+2] +
                                             p->voronoi.vertices[3*i+2] +
                                             p->voronoi.vertices[3*m+2]);
                    k = m;
                    l = n;
                    m = p->voronoi.edges[6*k+l];
                    p->voronoi.edges[6*k+l] = -1 - m;
                }
                p->voronoi.face_areas[p->voronoi.nface] = 0.5f*area;
                p->voronoi.face_midpoints[3*p->voronoi.nface+0] = midpoint[0]/area/3.0f;
                p->voronoi.face_midpoints[3*p->voronoi.nface+1] = midpoint[1]/area/3.0f;
                p->voronoi.face_midpoints[3*p->voronoi.nface+2] = midpoint[2]/area/3.0f;
                /* face midpoint was calculated relative to particle position */
                p->voronoi.face_midpoints[3*p->voronoi.nface+0] += p->x[0];
                p->voronoi.face_midpoints[3*p->voronoi.nface+1] += p->x[1];
                p->voronoi.face_midpoints[3*p->voronoi.nface+2] += p->x[2];
                p->voronoi.nface++;
            }
        }
    }

    /* update ngbs */
    for(i = 0; i < 300; i++){
        p->voronoi.ngbs[i] = newngbs[i];
    }

    /* unmark edges */
    for(i = 0; i < p->voronoi.nvert; i++){
        for(j = 0; j < 3; j++){
            if(p->voronoi.edges[6*i+j] < 0){
                p->voronoi.edges[6*i+j] = -p->voronoi.edges[6*i+j]-1;
            } else {
                error("edge inconsistency");
            }
        }
    }

}

/**
 * @brief Get the index of the face between pi and pj
 *
 * If the face does not exist, we return -1.
 */

__attribute__((always_inline)) INLINE static int voronoi_get_face_index(
    struct part *pi, struct part *pj){

    /* alternative (less precise) neighbour finding */
/*    float ri2, rj2;*/
/*    for(int i = 0; i < pi->voronoi.nface; i++){*/
/*        ri2 = 0.0f;*/
/*        for(int j = 0; j < 3; j++){*/
/*            float x = pi->voronoi.face_midpoints[3*i+j] - pi->x[j];*/
/*            if(x < -0.5f){*/
/*                x += 1.0f;*/
/*            }*/
/*            if(x > 0.5f){*/
/*                x -= 1.0f;*/
/*            }*/
/*            ri2 += x * x;*/
/*        }*/
/*        rj2 = 0.0f;*/
/*        for(int j = 0; j < 3; j++){*/
/*            float x = pi->voronoi.face_midpoints[3*i+j] - pj->x[j];*/
/*            if(x < -0.5f){*/
/*                x += 1.0f;*/
/*            }*/
/*            if(x > 0.5f){*/
/*                x -= 1.0f;*/
/*            }*/
/*            rj2 += x * x;*/
/*        }*/
/*        if(fabs(ri2-rj2) < 1.e-8){*/
/*            return i;*/
/*        }*/
/*    }*/
/*    return -1;*/

    /* explicit neighbour tracking */
    int i;
    for(i = 0; i < pi->voronoi.nface; i++){
        if(pi->voronoi.ngbs[i] == pj->id){
            return i;
        }
    }
    /* particles are not cell neighbours */
    return -1;

}

/**
 * @brief Calculate the velocity of the face between pi and pj
 */

__attribute__((always_inline)) INLINE static void voronoi_get_face_velocity(
    float r2, float *dx, struct part *pi, struct part *pj, float *midface, float *vface){

    float xd[3];
    float vproj;

    vface[0] = 0.5f * (pi->primitives.v[0] + pj->primitives.v[0]);
    vface[1] = 0.5f * (pi->primitives.v[1] + pj->primitives.v[1]);
    vface[2] = 0.5f * (pi->primitives.v[2] + pj->primitives.v[2]);

    /* we cannot simply calculate xd as 0.5*(xi+xj), since this does not take
       the periodic corrections into account */
    xd[0] = pi->x[0] - 0.5f*dx[0];
    xd[1] = pi->x[1] - 0.5f*dx[1];
    xd[2] = pi->x[2] - 0.5f*dx[2];

    vproj = (pi->primitives.v[0] - pj->primitives.v[0]) * (midface[0] - xd[0])
              + (pi->primitives.v[1] - pj->primitives.v[1]) * (midface[1] - xd[1])
              + (pi->primitives.v[2] - pj->primitives.v[2]) * (midface[2] - xd[2]);

    /* minus sign due to the reverse definition of dx w.r.t. Springel 2010 */
    vface[0] -= vproj * dx[0] / r2;
    vface[1] -= vproj * dx[1] / r2;
    vface[2] -= vproj * dx[2] / r2;

}

#endif /* SWIFT_VORONOI_H */
