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
#define VORONOI_BOX_FRONT   18446744073709551600llu
#define VORONOI_BOX_BACK    18446744073709551601llu
#define VORONOI_BOX_TOP     18446744073709551602llu
#define VORONOI_BOX_BOTTOM  18446744073709551603llu
#define VORONOI_BOX_LEFT    18446744073709551604llu
#define VORONOI_BOX_RIGHT   18446744073709551605llu

struct voronoi_box {

  /* anchor of the box */
  float anchor[3];

  /* side lenghs of the box */
  float sides[3];

};

/* Maximal number of vertices allowed for a Voronoi cell during the cell
   construction */
#define VORONOI_MAX_NUM_VERT 300
#define VORONOI_MAX_NUM_EDGE (3*VORONOI_MAX_NUM_VERT)

#define VORONOI_TOLERANCE 1.e-7

/* Internal representation of a Voronoi cell, which has more memory to store
   intermediate vertices, and which uses a more efficient memory layout */
struct voronoi_cell {

  /* Number of vertices */
  int nvert;

  /* Vertex coordinates */
  float vertices[3*VORONOI_MAX_NUM_VERT];

  /* Number of edges for every vertex */
  int orders[VORONOI_MAX_NUM_VERT];

  /* Offsets of the edges, edgeindices and neighbours corresponding to a
     particular vertex in the internal arrays */
  int offsets[VORONOI_MAX_NUM_VERT];

  /* Edge information */
  int edges[VORONOI_MAX_NUM_EDGE];

  /* Additional edge information */
  int edgeindices[VORONOI_MAX_NUM_EDGE];

  /* Neighbour information */
  unsigned long long ngbs[VORONOI_MAX_NUM_EDGE];

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

  memcpy(c->vertices, p->geometry.vertices, 3*c->nvert*sizeof(float));
  memcpy(c->orders, p->geometry.orders, c->nvert*sizeof(int));
  memcpy(c->edges, p->geometry.edges, numedge*sizeof(int));
  memcpy(c->edgeindices, p->geometry.edgeindices, numedge*sizeof(int));
  memcpy(c->ngbs, p->geometry.ngbs, numedge*sizeof(unsigned long long));

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

  int numedge;

  numedge = c->offsets[c->nvert-1]+c->orders[c->nvert-1];
  if(numedge > VORONOI_MAXEDGE){
    error("Too many edges!");
  }

  p->geometry.nvert = c->nvert;

  memcpy(p->geometry.vertices, c->vertices, 3*c->nvert*sizeof(float));
  memcpy(p->geometry.orders, c->orders, c->nvert*sizeof(int));
  memcpy(p->geometry.edges, c->edges, numedge*sizeof(int));
  memcpy(p->geometry.edgeindices, c->edgeindices, numedge*sizeof(int));
  memcpy(p->geometry.ngbs, c->ngbs, numedge*sizeof(unsigned long long));

}

/**
 * @brief For debugging purposes
 */
__attribute__((always_inline)) INLINE static void voronoi_print_cell(
  struct part *p){

  int i, j;
  struct voronoi_cell c;

  voronoi_set_cell_values(p, &c);

  for(i = 0; i < c.nvert; i++){
    fprintf(stderr, "%i: %g %g %g (%i)\n", i, c.vertices[3*i],
                                           c.vertices[3*i+1], c.vertices[3*i+2],
                                           c.orders[i]);
    for(j = 0; j < c.orders[i]; j++){
      fprintf(stderr, "%i (%i)", c.edges[c.offsets[i]+j],
                                 c.edgeindices[c.offsets[i]+j]);
      if(j < c.orders[i]-1){
        fprintf(stderr, "\t");
      } else {
        fprintf(stderr, "\n");
      }
    }
  }
  fprintf(stderr, "\n");
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
    struct part *p, struct voronoi_box *b){

  p->geometry.nvert = 8;

/*  b->anchor[0] = -1.0f;*/
/*  b->anchor[1] = -1.0f;*/
/*  b->anchor[2] = -1.0f;*/
/*  b->sides[0] = 3.0f;*/
/*  b->sides[1] = 3.0f;*/
/*  b->sides[2] = 3.0f;*/

  /* (0, 0, 0) -- 0 */
  p->geometry.vertices[ 0] = b->anchor[0] - p->x[0];
  p->geometry.vertices[ 1] = b->anchor[1] - p->x[1];
  p->geometry.vertices[ 2] = b->anchor[2] - p->x[2];

  /* (0, 0, 1)-- 1 */
  p->geometry.vertices[ 3] = b->anchor[0] - p->x[0];
  p->geometry.vertices[ 4] = b->anchor[1] - p->x[1];
  p->geometry.vertices[ 5] = b->anchor[2] + b->sides[2] - p->x[2];

  /* (0, 1, 0) -- 2 */
  p->geometry.vertices[ 6] = b->anchor[0] - p->x[0];
  p->geometry.vertices[ 7] = b->anchor[1] + b->sides[1] - p->x[1];
  p->geometry.vertices[ 8] = b->anchor[2] - p->x[2];

  /* (0, 1, 1) -- 3 */
  p->geometry.vertices[ 9] = b->anchor[0] - p->x[0];
  p->geometry.vertices[10] = b->anchor[1] + b->sides[1] - p->x[1];
  p->geometry.vertices[11] = b->anchor[2] + b->sides[2] - p->x[2];

  /* (1, 0, 0) -- 4 */
  p->geometry.vertices[12] = b->anchor[0] + b->sides[0] - p->x[0];
  p->geometry.vertices[13] = b->anchor[1] - p->x[1];
  p->geometry.vertices[14] = b->anchor[2] - p->x[2];

  /* (1, 0, 1) -- 5 */
  p->geometry.vertices[15] = b->anchor[0] + b->sides[0] - p->x[0];
  p->geometry.vertices[16] = b->anchor[1] - p->x[1];
  p->geometry.vertices[17] = b->anchor[2] + b->sides[2] - p->x[2];

  /* (1, 1, 0) -- 6 */
  p->geometry.vertices[18] = b->anchor[0] + b->sides[0] - p->x[0];
  p->geometry.vertices[19] = b->anchor[1] + b->sides[1] - p->x[1];
  p->geometry.vertices[20] = b->anchor[2] - p->x[2];

  /* (1, 1, 1) -- 7 */
  p->geometry.vertices[21] = b->anchor[0] + b->sides[0] - p->x[0];
  p->geometry.vertices[22] = b->anchor[1] + b->sides[1] - p->x[1];
  p->geometry.vertices[23] = b->anchor[2] + b->sides[2] - p->x[2];

  p->geometry.orders[0] = 3;
  p->geometry.orders[1] = 3;
  p->geometry.orders[2] = 3;
  p->geometry.orders[3] = 3;
  p->geometry.orders[4] = 3;
  p->geometry.orders[5] = 3;
  p->geometry.orders[6] = 3;
  p->geometry.orders[7] = 3;

  /* edges are ordered counterclockwise w.r.t. a vector pointing from the
     cell generator to the vertex
     (0, 0, 0) corner */
  p->geometry.edges[ 0] = 1;
  p->geometry.edges[ 1] = 2;
  p->geometry.edges[ 2] = 4;
  p->geometry.edgeindices[ 0] = 0;
  p->geometry.edgeindices[ 1] = 2;
  p->geometry.edgeindices[ 2] = 0;

  /* (0, 0, 1) corner */
  p->geometry.edges[ 3] = 0;
  p->geometry.edges[ 4] = 5;
  p->geometry.edges[ 5] = 3;
  p->geometry.edgeindices[ 3] = 0;
  p->geometry.edgeindices[ 4] = 2;
  p->geometry.edgeindices[ 5] = 1;

  /* (0, 1, 0) corner */
  p->geometry.edges[ 6] = 3;
  p->geometry.edges[ 7] = 6;
  p->geometry.edges[ 8] = 0;
  p->geometry.edgeindices[ 6] = 0;
  p->geometry.edgeindices[ 7] = 0;
  p->geometry.edgeindices[ 8] = 1;

  /* (0, 1, 1) corner */
  p->geometry.edges[ 9] = 2;
  p->geometry.edges[10] = 1;
  p->geometry.edges[11] = 7;
  p->geometry.edgeindices[ 9] = 0;
  p->geometry.edgeindices[10] = 2;
  p->geometry.edgeindices[11] = 0;

  /* (1, 0, 0) corner */
  p->geometry.edges[12] = 0;
  p->geometry.edges[13] = 6;
  p->geometry.edges[14] = 5;
  p->geometry.edgeindices[12] = 2;
  p->geometry.edgeindices[13] = 2;
  p->geometry.edgeindices[14] = 0;

  /* (1, 0, 1) corner */
  p->geometry.edges[15] = 4;
  p->geometry.edges[16] = 7;
  p->geometry.edges[17] = 1;
  p->geometry.edgeindices[15] = 2;
  p->geometry.edgeindices[16] = 1;
  p->geometry.edgeindices[17] = 1;

  /* (1, 1, 0) corner */
  p->geometry.edges[18] = 2;
  p->geometry.edges[19] = 7;
  p->geometry.edges[20] = 4;
  p->geometry.edgeindices[18] = 1;
  p->geometry.edgeindices[19] = 2;
  p->geometry.edgeindices[20] = 1;

  /* (1, 1, 1) corner */
  p->geometry.edges[21] = 3;
  p->geometry.edges[22] = 5;
  p->geometry.edges[23] = 6;
  p->geometry.edgeindices[21] = 2;
  p->geometry.edgeindices[22] = 1;
  p->geometry.edgeindices[23] = 1;

  /* ngbs[3*i+j] is the neighbour corresponding to the plane clockwise of
     edge j of vertex i (when going from edge j to vertex i)
     we set them to a ridiculously large value to be able to track faces without
     neighbour */
  p->geometry.ngbs[ 0] = VORONOI_BOX_FRONT;  /* (000) - (001) */
  p->geometry.ngbs[ 1] = VORONOI_BOX_LEFT;   /* (000) - (010) */
  p->geometry.ngbs[ 2] = VORONOI_BOX_BOTTOM; /* (000) - (100) */

  p->geometry.ngbs[ 3] = VORONOI_BOX_LEFT;   /* (001) - (000) */
  p->geometry.ngbs[ 4] = VORONOI_BOX_FRONT;  /* (001) - (101) */
  p->geometry.ngbs[ 5] = VORONOI_BOX_TOP;    /* (001) - (011) */

  p->geometry.ngbs[ 6] = VORONOI_BOX_LEFT;   /* (010) - (011) */
  p->geometry.ngbs[ 7] = VORONOI_BOX_BACK;   /* (010) - (110) */
  p->geometry.ngbs[ 8] = VORONOI_BOX_BOTTOM; /* (010) - (000) */

  p->geometry.ngbs[ 9] = VORONOI_BOX_BACK;   /* (011) - (010) */
  p->geometry.ngbs[10] = VORONOI_BOX_LEFT;   /* (011) - (001) */
  p->geometry.ngbs[11] = VORONOI_BOX_TOP;    /* (011) - (111) */

  p->geometry.ngbs[12] = VORONOI_BOX_FRONT;  /* (100) - (000) */
  p->geometry.ngbs[13] = VORONOI_BOX_BOTTOM; /* (100) - (110) */
  p->geometry.ngbs[14] = VORONOI_BOX_RIGHT;  /* (100) - (101) */

  p->geometry.ngbs[15] = VORONOI_BOX_FRONT;  /* (101) - (100) */
  p->geometry.ngbs[16] = VORONOI_BOX_RIGHT;  /* (101) - (111) */
  p->geometry.ngbs[17] = VORONOI_BOX_TOP;    /* (101) - (001) */

  p->geometry.ngbs[18] = VORONOI_BOX_BOTTOM; /* (110) - (010) */
  p->geometry.ngbs[19] = VORONOI_BOX_BACK;   /* (110) - (111) */
  p->geometry.ngbs[20] = VORONOI_BOX_RIGHT;  /* (110) - (100) */

  p->geometry.ngbs[21] = VORONOI_BOX_BACK;   /* (111) - (011) */
  p->geometry.ngbs[22] = VORONOI_BOX_TOP;    /* (111) - (101) */
  p->geometry.ngbs[23] = VORONOI_BOX_RIGHT;  /* (111) - (110) */
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
  int visitflags[VORONOI_MAX_NUM_VERT];
  int dstack[VORONOI_MAX_NUM_VERT];
  int dstack_size = 1;
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
      if(c.nvert == VORONOI_MAX_NUM_VERT){
        error("Too many vertices!");
      }
      c.orders[vindex] = k;
      c.offsets[vindex] = c.offsets[vindex-1] + c.orders[vindex-1];
      if(c.offsets[vindex] + k >= VORONOI_MAX_NUM_EDGE){
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
      if(c.nvert == VORONOI_MAX_NUM_VERT){
        error("Too many vertices!");
      }
      c.orders[vindex] = k;
      c.offsets[vindex] = c.offsets[vindex-1] + c.orders[vindex-1];
      if(c.offsets[vindex] + k >= VORONOI_MAX_NUM_EDGE){
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
    if(c.nvert == VORONOI_MAX_NUM_VERT){
      error("Too many vertices!");
    }
    c.orders[vindex] = 3;
    c.offsets[vindex] = c.offsets[vindex-1] + c.orders[vindex-1];
    if(c.offsets[vindex] + 3 >= VORONOI_MAX_NUM_EDGE){
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

  int cp;
  int iqs;
  int new_double_edge;

  cp = vindex;
  rp = vindex;
  while(qp != up || qs != us){

    lp = voronoi_get_edge(&c, qp, qs);
    lw = voronoi_test_vertex(&c.vertices[3*lp], dx, r2, &l);
    if(lw == 0){

      k = 1;
      if(double_edge){
        k = 0;
      }
      qs = voronoi_get_edgeindex(&c, qp, qs);
      qp = lp;
      iqs = qs;

      k++;
      qs++;
      if(qs == c.orders[qp]){
          qs = 0;
      }
      lp = voronoi_get_edge(&c, qp, qs);
      lw = voronoi_test_vertex(&c.vertices[3*lp], dx, r2, &l);
      while(lw == -1){
        k++;
        qs++;
        if(qs == c.orders[qp]){
            qs = 0;
        }
        lp = voronoi_get_edge(&c, qp, qs);
        lw = voronoi_test_vertex(&c.vertices[3*lp], dx, r2, &l);
      }

      j = visitflags[qp];
      if(qp == up && qs == us){
        new_double_edge = 0;
        if(j > 0){
          k += c.orders[j];
        }
      } else {
        if(j > 0){
            k += c.orders[j];
            if(lw == 0){
              i = -visitflags[lp];
              if(i > 0){
                if(voronoi_get_edge(&c, i, c.orders[i]-1) == j){
                  new_double_edge = 1;
                  k--;
                } else {
                  new_double_edge = 0;
                }
              } else {
                if(j == rp && lp == up && voronoi_get_edge(&c, qp, qs) == us){
                  new_double_edge = 1;
                  k--;
                } else {
                  new_double_edge = 0;
                }
              }
            } else {
              new_double_edge = 0;
            }
        } else {
          if(lw == 0){
            i = -visitflags[lp];
            if(i == cp){
              new_double_edge = 1;
              k--;
            } else {
              new_double_edge = 0;
            }
          } else {
            new_double_edge = 0;
          }
        }
      }

      if(j > 0){
        voronoi_print_cell(pi);
        voronoi_print_cell(pj);
        error("Case not handled!");
      }

      /* create new order k vertex */
      vindex = c.nvert;
      c.nvert++;
      if(c.nvert == VORONOI_MAX_NUM_VERT){
        error("Too many vertices!");
      }
      c.orders[vindex] = k;
      c.offsets[vindex] = c.offsets[vindex-1] + c.orders[vindex-1];
      if(c.offsets[vindex] + k >= VORONOI_MAX_NUM_EDGE){
        error("Too many edges!");
      }

      visitflags[vindex] = 0;
      c.vertices[3*vindex+0] = c.vertices[3*qp+0];
      c.vertices[3*vindex+1] = c.vertices[3*qp+1];
      c.vertices[3*vindex+2] = c.vertices[3*qp+2];
      visitflags[qp] = -vindex;
      dstack[dstack_size] = qp;
      dstack_size++;
      j = vindex;
      i = 0;

      if(!double_edge){
        voronoi_set_ngb(&c, j, i, pj->id);
        voronoi_set_edge(&c, j, i, cp);
        voronoi_set_edgeindex(&c, j, i, cs);
        voronoi_set_edge(&c, cp, cs, j);
        voronoi_set_edgeindex(&c, cp, cs, i);
        i++;
      }

      qs = iqs;
      iqs = k-1;
      if(new_double_edge){
        iqs = k;
      }
      while(i < iqs){
        qs++;
        if(qs == c.orders[qp]){
            qs = 0;
        }
        lp = voronoi_get_edge(&c, qp, qs);
        ls = voronoi_get_edgeindex(&c, qp, qs);
        voronoi_set_ngb(&c, j, i, voronoi_get_ngb(&c, qp, qs));
        voronoi_set_edge(&c, j, i, lp);
        voronoi_set_edgeindex(&c, j, i, ls);
        voronoi_set_edge(&c, lp, ls, j);
        voronoi_set_edgeindex(&c, lp, ls, i);
        voronoi_set_edge(&c, qp, qs, -1);
        i++;
      }
      qs++;
      if(qs == c.orders[qp]){
        qs = 0;
      }
      cs = i;
      cp = j;

      if(new_double_edge){
        voronoi_set_ngb(&c, j, 0, voronoi_get_ngb(&c, qp, qs));
      } else {
        voronoi_set_ngb(&c, j, cs, voronoi_get_ngb(&c, qp, qs));
      }

      double_edge = new_double_edge;

    } else { /* if(lw == 0) */

      if(lw == 1){
        qs = voronoi_get_edgeindex(&c, qp, qs) + 1;
        if(qs == c.orders[lp]){
            qs = 0;
        }
        qp = lp;
        q = l;
        dstack[dstack_size] = qp;
        dstack_size++;
      } else {

        r = q/(q-l);
        l = 1.0f - r;

        /* create new order 3 vertex */
        vindex = c.nvert;
        c.nvert++;
        if(c.nvert == VORONOI_MAX_NUM_VERT){
          error("Too many vertices!");
        }
        c.orders[vindex] = 3;
        c.offsets[vindex] = c.offsets[vindex-1] + c.orders[vindex-1];
        if(c.offsets[vindex] + 3 >= VORONOI_MAX_NUM_EDGE){
          error("Too many edges!");
        }

        c.vertices[3*vindex+0] = c.vertices[3*lp+0] * r +
                                 c.vertices[3*qp+0] * l;
        c.vertices[3*vindex+1] = c.vertices[3*lp+1] * r +
                                 c.vertices[3*qp+1] * l;
        c.vertices[3*vindex+2] = c.vertices[3*lp+2] * r +
                                 c.vertices[3*qp+2] * l;

        ls = voronoi_get_edgeindex(&c, qp, qs);
        voronoi_set_edge(&c, vindex, 0, cp);
        voronoi_set_edge(&c, vindex, 1, lp);
        voronoi_set_edgeindex(&c, vindex, 0, cs);
        voronoi_set_edgeindex(&c, vindex, 1, ls);
        voronoi_set_edge(&c, lp, ls, vindex);
        voronoi_set_edgeindex(&c, lp, ls, 1);
        voronoi_set_edge(&c, cp, cs, vindex);
        voronoi_set_edgeindex(&c, cp, cs, 0);
        voronoi_set_edge(&c, qp, qs, -1);

        voronoi_set_ngb(&c, vindex, 0, pj->id);
        voronoi_set_ngb(&c, vindex, 1, voronoi_get_ngb(&c, qp, qs));
        voronoi_set_ngb(&c, vindex, 2, voronoi_get_ngb(&c, lp, ls));

        qs++;
        if(qs == c.orders[qp]){
          qs = 0;
        }
        cp = vindex;
        cs = 2;
      } /* if(lw == 1) */

    } /* if(lw == 0) */

  } /* while() */

  voronoi_set_edge(&c, cp, cs, rp);
  voronoi_set_edge(&c, rp, 0, cp);
  voronoi_set_edgeindex(&c, cp, cs, 0);
  voronoi_set_edgeindex(&c, rp, 0, cs);

  for(i = 0; i < dstack_size; i++){
    for(j = 0; j < c.orders[dstack[i]]; j++){
      if(voronoi_get_edge(&c, dstack[i], j) >= 0){
        dstack[dstack_size] = voronoi_get_edge(&c, dstack[i], j);
        dstack_size++;
        voronoi_set_edge(&c, dstack[i], j, -1);
        voronoi_set_edgeindex(&c, dstack[i], j, -1);
      }
    }
  }

  /* remove deleted vertices from all arrays */
  struct voronoi_cell new_cell;
  int m, n;
  for(vindex = 0; vindex < c.nvert; vindex++){
    j = vindex;
    /* find next edge that is not deleted */
    while(j < c.nvert && voronoi_get_edge(&c, j, 0) < 0){
      j++;
    }

    if(j == c.nvert){
      /* ready */
      break;
    }

    /* copy vertices */
    new_cell.vertices[3*vindex+0] = c.vertices[3*j+0];
    new_cell.vertices[3*vindex+1] = c.vertices[3*j+1];
    new_cell.vertices[3*vindex+2] = c.vertices[3*j+2];

    /* copy order */
    new_cell.orders[vindex] = c.orders[j];

    /* set offset */
    if(vindex){
      new_cell.offsets[vindex] = new_cell.offsets[vindex-1]
                                   + new_cell.orders[vindex-1];
    } else {
      new_cell.offsets[vindex] = 0;
    }

    /* copy edges, edgeindices and ngbs */
    for(k = 0; k < c.orders[j]; k++){
      voronoi_set_edge(&new_cell, vindex, k, voronoi_get_edge(&c, j, k));
      voronoi_set_edgeindex(&new_cell, vindex, k,
                            voronoi_get_edgeindex(&c, j, k));
      voronoi_set_ngb(&new_cell, vindex, k, voronoi_get_ngb(&c, j, k));
    }

    /* update other edges */
    for(k = 0; k < c.orders[j]; k++){
      m = voronoi_get_edge(&c, j, k);
      n = voronoi_get_edgeindex(&c, j, k);
      if(m < vindex){
        voronoi_set_edge(&new_cell, m, n, vindex);
      } else {
        voronoi_set_edge(&c, m, n, vindex);
      }
    }

    /* deactivate edge */
    voronoi_set_edge(&c, j, 0, -1);
  }
  new_cell.nvert = vindex;

  voronoi_set_particle_values(pi, &new_cell);

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
  struct voronoi_cell c;

  voronoi_set_cell_values(p, &c);

  /* we need to calculate the volume of the tetrahedra formed by the first
     vertex and the triangles that make up the other faces
     since we do not store faces explicitly, this means keeping track of the
     edges that have been processed somehow
     we follow the method used in voro++ and "flip" processed edges to
     negative values
     this also means that we need to process all triangles corresponding to
     an edge at once */
  p->voronoi.volume = 0.0f;
  v1[0] = c.vertices[0];
  v1[1] = c.vertices[1];
  v1[2] = c.vertices[2];
  p->voronoi.centroid[0] = 0.0f;
  p->voronoi.centroid[1] = 0.0f;
  p->voronoi.centroid[2] = 0.0f;
  
  /* loop over all vertices (except the first one) */
  for(i = 1; i < c.nvert; i++){

    v2[0] = c.vertices[3*i+0];
    v2[1] = c.vertices[3*i+1];
    v2[2] = c.vertices[3*i+2];
    
    /*  loop over the edges of the vertex*/
    for(j = 0; j < c.orders[i]; j++){

      k = voronoi_get_edge(&c, i, j);

      if(k >= 0){

        /* mark the edge as processed */
        voronoi_set_edge(&c, i, j, -k-1);

        l = voronoi_get_edgeindex(&c, i, j) + 1;
        if(l == c.orders[k]){
          l = 0;
        }
        v3[0] = c.vertices[3*k+0];
        v3[1] = c.vertices[3*k+1];
        v3[2] = c.vertices[3*k+2];
        m = voronoi_get_edge(&c, k, l);
        voronoi_set_edge(&c, k, l, -1-m);

        while(m != i){
          n = voronoi_get_edgeindex(&c, k, l) + 1;
          if(n == c.orders[m]){
            n = 0;
          }
          v4[0] = c.vertices[3*m+0];
          v4[1] = c.vertices[3*m+1];
          v4[2] = c.vertices[3*m+2];
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
          m = voronoi_get_edge(&c, k, l);
          voronoi_set_edge(&c, k, l, -1-m);
        } /* while() */

      } /* if(k >= 0) */

    } /* for(j) */

  } /* for(i) */

  p->voronoi.centroid[0] /= p->voronoi.volume;
  p->voronoi.centroid[1] /= p->voronoi.volume;
  p->voronoi.centroid[2] /= p->voronoi.volume;

  /* centroid was calculated relative w.r.t. particle position */
  p->voronoi.centroid[0] += p->x[0];
  p->voronoi.centroid[1] += p->x[1];
  p->voronoi.centroid[2] += p->x[2];

}

__attribute__((always_inline)) INLINE static void voronoi_calculate_faces(
    struct part *p){

  int i, j, k, l, m, n;
  float area;
  float midpoint[3];
  float u[3], v[3], w[3];
  float loc_area;
  struct voronoi_cell c;

  voronoi_set_cell_values(p, &c);

  p->voronoi.nface = 0;
  for(i = 0; i < c.nvert; i++){

    for(j = 0; j < c.orders[i]; j++){

      k = voronoi_get_edge(&c, i, j);

      if(k >= 0){

        p->voronoi.ngbs[p->voronoi.nface] = voronoi_get_ngb(&c, i, j);
        area = 0.;
        midpoint[0] = 0.;
        midpoint[1] = 0.;
        midpoint[2] = 0.;
        voronoi_set_edge(&c, i, j, -1-k);
        l = voronoi_get_edgeindex(&c, i, j) + 1;
        if(l == c.orders[k]){
          l = 0;
        }
        m = voronoi_get_edge(&c, k, l);
        voronoi_set_edge(&c, k, l, -1-m);

        while(m != i){
          n = voronoi_get_edgeindex(&c, k, l) + 1;
          if(n == c.orders[m]){
            n = 0;
          }
          u[0] = c.vertices[3*k+0] - c.vertices[3*i+0];
          u[1] = c.vertices[3*k+1] - c.vertices[3*i+1];
          u[2] = c.vertices[3*k+2] - c.vertices[3*i+2];
          v[0] = c.vertices[3*m+0] - c.vertices[3*i+0];
          v[1] = c.vertices[3*m+1] - c.vertices[3*i+1];
          v[2] = c.vertices[3*m+2] - c.vertices[3*i+2];
          w[0] = u[1]*v[2] - u[2]*v[1];
          w[1] = u[2]*v[0] - u[0]*v[2];
          w[2] = u[0]*v[1] - u[1]*v[0];
          loc_area = sqrtf(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);
          area += loc_area;
          midpoint[0] += loc_area*(c.vertices[3*k+0] +
                                   c.vertices[3*i+0] +
                                   c.vertices[3*m+0]);
          midpoint[1] += loc_area*(c.vertices[3*k+1] +
                                   c.vertices[3*i+1] +
                                   c.vertices[3*m+1]);
          midpoint[2] += loc_area*(c.vertices[3*k+2] +
                                   c.vertices[3*i+2] +
                                   c.vertices[3*m+2]);
          k = m;
          l = n;
          m = voronoi_get_edge(&c, k, l);
          voronoi_set_edge(&c, k, l, -1-m);
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

        if(p->voronoi.nface == VORONOI_MAXFACE){
          error("Too many faces!");
        }

      } /* if(k >= 0) */

    } /* for(j) */

  } /* for(i) */

}

/**
 * @brief Get the index of the face between pi and pj
 *
 * If the face does not exist, we return -1.
 */
__attribute__((always_inline)) INLINE static int voronoi_get_face_index(
    struct part *pi, struct part *pj){

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
    float r2, float *dx, struct part *pi, struct part *pj, float *midface,
    float *vface){

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
