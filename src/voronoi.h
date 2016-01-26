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

/* Box boundary flags used to signal cells neighbouring the box boundary */
#define BOX_FRONT   0xfffffff0
#define BOX_BACK    0xfffffff1
#define BOX_TOP     0xfffffff2
#define BOX_BOTTOM  0xfffffff3
#define BOX_LEFT    0xfffffff4
#define BOX_RIGHT   0xfffffff5

/**
 * @brief Check if the given vertex is above or below the cutting plane
 */

__attribute__((always_inline)) INLINE static int test_vertex(
    float *v, float *dx, float r2, float *test){

    *test = v[0]*dx[0] + v[1]*dx[1] + v[2]*dx[2] - r2;
    if(*test <= 0.){
        return -1;
    }
    if(*test > 0.){
        return 1;
    }
    return 0;

}

/**
 * @brief Initialize the cell as a cube with side 2*h, centered on the particle
 *  position
 */

__attribute__((always_inline)) INLINE static void voronoi_initialize(
    struct part *p, float h){

    p->voronoi.nvert = 8;

    /* (0, 0, 0) -- 0 */
    p->voronoi.vertices[ 0] = -h;
    p->voronoi.vertices[ 1] = -h;
    p->voronoi.vertices[ 2] = -h;

    /* (0, 0, 1)-- 1 */
    p->voronoi.vertices[ 3] = -h;
    p->voronoi.vertices[ 4] = -h;
    p->voronoi.vertices[ 5] =  h;

    /* (0, 1, 0) -- 2*/
    p->voronoi.vertices[ 6] = -h;
    p->voronoi.vertices[ 7] =  h;
    p->voronoi.vertices[ 8] = -h;

    /* (0, 1, 1) -- 3 */
    p->voronoi.vertices[ 9] = -h;
    p->voronoi.vertices[10] =  h;
    p->voronoi.vertices[11] =  h;

    /* (1, 0, 0) -- 4 */
    p->voronoi.vertices[12] =  h;
    p->voronoi.vertices[13] = -h;
    p->voronoi.vertices[14] = -h;

    /* (1, 0, 1) -- 5 */
    p->voronoi.vertices[15] =  h;
    p->voronoi.vertices[16] = -h;
    p->voronoi.vertices[17] =  h;

    /* (1, 1, 0) -- 6 */
    p->voronoi.vertices[18] =  h;
    p->voronoi.vertices[19] =  h;
    p->voronoi.vertices[20] = -h;

    /* (1, 1, 1) -- 7 */
    p->voronoi.vertices[21] =  h;
    p->voronoi.vertices[22] =  h;
    p->voronoi.vertices[23] =  h;

    /* edges are ordered counterclockwise w.r.t. a vector pointing from the
       cell generator to the vertex */
    /* (0, 0, 0) corner */
    p->voronoi.edges[ 0] = 1;
    p->voronoi.edges[ 1] = 2;
    p->voronoi.edges[ 2] = 4;
    p->voronoi.edges[ 3] = 0;
    p->voronoi.edges[ 4] = 2;
    p->voronoi.edges[ 5] = 0;

    /* (0, 0, 1) corner */
    p->voronoi.edges[ 6] = 0;
    p->voronoi.edges[ 7] = 5;
    p->voronoi.edges[ 8] = 3;
    p->voronoi.edges[ 9] = 0;
    p->voronoi.edges[10] = 2;
    p->voronoi.edges[11] = 1;

    /* (0, 1, 0) corner */
    p->voronoi.edges[12] = 3;
    p->voronoi.edges[13] = 6;
    p->voronoi.edges[14] = 0;
    p->voronoi.edges[15] = 0;
    p->voronoi.edges[16] = 0;
    p->voronoi.edges[17] = 1;

    /* (0, 1, 1) corner */
    p->voronoi.edges[18] = 2;
    p->voronoi.edges[19] = 1;
    p->voronoi.edges[20] = 7;
    p->voronoi.edges[21] = 0;
    p->voronoi.edges[22] = 2;
    p->voronoi.edges[23] = 0;

    /* (1, 0, 0) corner */
    p->voronoi.edges[24] = 0;
    p->voronoi.edges[25] = 6;
    p->voronoi.edges[26] = 5;
    p->voronoi.edges[27] = 2;
    p->voronoi.edges[28] = 2;
    p->voronoi.edges[29] = 0;

    /* (1, 0, 1) corner */
    p->voronoi.edges[30] = 4;
    p->voronoi.edges[31] = 7;
    p->voronoi.edges[32] = 1;
    p->voronoi.edges[33] = 2;
    p->voronoi.edges[34] = 1;
    p->voronoi.edges[35] = 1;

    /* (1, 1, 0) corner */
    p->voronoi.edges[36] = 2;
    p->voronoi.edges[37] = 7;
    p->voronoi.edges[38] = 4;
    p->voronoi.edges[39] = 1;
    p->voronoi.edges[40] = 2;
    p->voronoi.edges[41] = 1;

    /* (1, 1, 1) corner */
    p->voronoi.edges[42] = 3;
    p->voronoi.edges[43] = 5;
    p->voronoi.edges[44] = 6;
    p->voronoi.edges[45] = 2;
    p->voronoi.edges[46] = 1;
    p->voronoi.edges[47] = 1;

    /* ngbs[3*i+j] is the neighbour corresponding to the plane clockwise of
       edge j of vertex i */
    p->voronoi.ngbs[ 0] = BOX_LEFT;     /* (000) - (001) */
    p->voronoi.ngbs[ 1] = BOX_BOTTOM;   /* (000) - (010) */
    p->voronoi.ngbs[ 2] = BOX_FRONT;    /* (000) - (100) */

    p->voronoi.ngbs[ 3] = BOX_FRONT;    /* (001) - (000) */
    p->voronoi.ngbs[ 4] = BOX_TOP;      /* (001) - (101) */
    p->voronoi.ngbs[ 5] = BOX_LEFT;     /* (001) - (011) */

    p->voronoi.ngbs[ 6] = BOX_BACK;     /* (010) - (011) */
    p->voronoi.ngbs[ 7] = BOX_BOTTOM;   /* (010) - (110) */
    p->voronoi.ngbs[ 8] = BOX_LEFT;     /* (010) - (000) */

    p->voronoi.ngbs[ 9] = BOX_LEFT;     /* (011) - (010) */
    p->voronoi.ngbs[10] = BOX_TOP;      /* (011) - (001) */
    p->voronoi.ngbs[11] = BOX_BACK;     /* (011) - (111) */

    p->voronoi.ngbs[12] = BOX_BOTTOM;   /* (100) - (000) */
    p->voronoi.ngbs[13] = BOX_RIGHT;    /* (100) - (110) */
    p->voronoi.ngbs[14] = BOX_FRONT;    /* (100) - (101) */

    p->voronoi.ngbs[15] = BOX_RIGHT;    /* (101) - (100) */
    p->voronoi.ngbs[16] = BOX_TOP;      /* (101) - (111) */
    p->voronoi.ngbs[17] = BOX_FRONT;    /* (101) - (001) */

    p->voronoi.ngbs[18] = BOX_BACK;     /* (110) - (010) */
    p->voronoi.ngbs[19] = BOX_RIGHT;    /* (110) - (111) */
    p->voronoi.ngbs[20] = BOX_BOTTOM;   /* (110) - (100) */

    p->voronoi.ngbs[21] = BOX_TOP;      /* (111) - (011) */
    p->voronoi.ngbs[22] = BOX_RIGHT;    /* (111) - (101) */
    p->voronoi.ngbs[23] = BOX_BACK;     /* (111) - (110) */
}

/**
 * @brief Intersect particle pi with particle pj and adapt its Voronoi cell
 *  structure
 *
 * odx = x_i - x_j!!!
 */

__attribute__((always_inline)) INLINE static void voronoi_intersect(
    float *odx, struct part *pi, struct part *pj) {

    /* check if the cell is still valid */
    /* it could happen that a previous interaction silently exited, leading to
       a cell that has to be reinitialized in a next iteration and that is no
       longer valid. Using the cell in a consecutive interaction will lead to
       SEGFAULTs */
    if(!pi->voronoi.nvert){
        return;
    }

    float r2;
    float dx[3];
    
    float u, l, q;
    int uw, up, us, lw, lp, ls, qw, qp, qs, cs, cp, rp;
    
    /* delete stack for vertices that should be removed */
    int dstack[1000];
    int dstack_size = 1;
    
    float r;
    int vindex;
    
    int i, j, k;
    int newnvert;
    
    dx[0] = -0.5f*odx[0];
    dx[1] = -0.5f*odx[1];
    dx[2] = -0.5f*odx[2];
    
    r2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];

    uw = test_vertex(&pi->voronoi.vertices[0], dx, r2, &u);
    up = 0;
    if(uw == 1){
        /* the first vertex lies outside the plane: u is positive
           we try to find a connected vertex that is closer to the plane */
        lp = pi->voronoi.edges[6*up];
        lw = test_vertex(&pi->voronoi.vertices[3*lp], dx, r2, &l);
        us = 1;
        while(us < 3 && l >= u){
            lp = pi->voronoi.edges[6*up+us];
            lw = test_vertex(&pi->voronoi.vertices[3*lp], dx, r2, &l);
            us++;
        }
        /* the last loop increased the counter by one too much, decrease it
           again */
        us--;
        if(l >= u){
            /* emergency exit */
            pi->voronoi.nvert = 0;
            return;
        }
        /* ok, we found a closer vertex
           get the edge of this vertex that corresponds to our initial vertex
           this is the definition of the 3 last elements of edges */
        ls = pi->voronoi.edges[6*up+3+us];
        /* now try to find an edge that intersects with the plane */
        while(lw == 1){
            /* use the vertex lp as new reference point */
            u = l;
            up = lp;
            us = 0;
            while(us < ls && l >= u){
                lp = pi->voronoi.edges[6*up+us];
                lw = test_vertex(&pi->voronoi.vertices[3*lp], dx, r2, &l);
                us++;
            }
            if(l >= u){
                /* skip initial vertex */
                us++;
                while(us < 3 && l >= u){
                    lp = pi->voronoi.edges[6*up+us];
                    lw = test_vertex(&pi->voronoi.vertices[3*lp], dx, r2, &l);
                    us++;
                }
                if(l >= u){
                    /* emergency exit */
                    pi->voronoi.nvert = 0;
                    return;
                }
            }
            us--;
            
            /* we found a closer vertex
               if it lies on the other side of the plane, we have found our
               intersected edge. If not, use this vertex as a new reference
               point and try again. */
            ls = pi->voronoi.edges[6*up+3+us];
        }
    } else {
        /* the initial vertex lies below the plane: u is negative
           we try to find a closer vertex */
        qp = pi->voronoi.edges[6*up];
        qw = test_vertex(&pi->voronoi.vertices[3*qp], dx, r2, &q);
        us = 1;
        while(us < 3 && u >= q){
            qp = pi->voronoi.edges[6*up+us];
            qw = test_vertex(&pi->voronoi.vertices[3*qp], dx, r2, &q);
            us++;
        }
        /* decrease us to correct for the extra increase */
        if(u >= q){
            /* no closer vertex is found. Since the Voronoi cell is convex, this
               means all vertices lie below the plane. The plane does not
               intersect the current cell and can be discarded. */
            return;
        } else {
            us--;
        }
        /* now try to find an edge that intersects with the plane */
        while(qw == -1){
            /* make the closer vertex the new initial guess */
            qs = pi->voronoi.edges[6*up+3+us];
            u = q;
            up = qp;
            us = 0;
            while(us < qs && u >= q){
                qp = pi->voronoi.edges[6*up+us];
                qw = test_vertex(&pi->voronoi.vertices[3*qp], dx, r2, &q);
                us++;
            }
            if(u >= q){
                /* skip initial vertex */
                us++;
                while(us < 3 && u >= q){
                    qp = pi->voronoi.edges[6*up+us];
                    qw = test_vertex(&pi->voronoi.vertices[3*qp], dx, r2, &q);
                    us++;
                }
                if(u >= q){
                    /* no closer vertex was found: the entire cell lies below
                       the plane and the plane does not interact with it */
                    return;
                }
            }
            us--;
        }
        /* store the information about the point below and above the plane */
        lp = up;
        ls = us;
        l = u;
        up = qp;
        us = pi->voronoi.edges[6*lp+3+ls];
        u = q;
    }
    
    /* at this point, lp, ls and l correspond to respectively the vertex, edge
       and distance to the plane for the point below the plane
       up, us and u correspond to the same values for the point above the plane */
    
    /* note: this division could go wrong... */
    if(u-l == 0.0f){
        error("Division by zero!");
    }
    r = u/(u-l);
    l = 1. - r;
    /* add new vertex */
    vindex = pi->voronoi.nvert;
    pi->voronoi.nvert++;
    
    if(vindex == 100){
        /* too many vertices, emergency exit */
        pi->voronoi.nvert = 0;
        return;
    }
    
    pi->voronoi.vertices[3*vindex+0] = pi->voronoi.vertices[3*lp+0]*r + pi->voronoi.vertices[3*up+0]*l;
    pi->voronoi.vertices[3*vindex+1] = pi->voronoi.vertices[3*lp+1]*r + pi->voronoi.vertices[3*up+1]*l;
    pi->voronoi.vertices[3*vindex+2] = pi->voronoi.vertices[3*lp+2]*r + pi->voronoi.vertices[3*up+2]*l;
    /* remove up by adding it to the delete stack */
    dstack[0] = up;
    /* connect lp to the newly created vertex */
    pi->voronoi.edges[6*vindex+1] = lp;
    pi->voronoi.edges[6*vindex+4] = ls;
    pi->voronoi.edges[6*lp+ls] = vindex;
    /* 1 because we choose 1 above */
    pi->voronoi.edges[6*lp+3+ls] = 1;
    pi->voronoi.edges[6*up+us] = -1;
    /* ngbs */
    pi->voronoi.ngbs[3*vindex] = pj->id;
    pi->voronoi.ngbs[3*vindex+1] = pi->voronoi.ngbs[3*up+us];
    pi->voronoi.ngbs[3*vindex+2] = pi->voronoi.ngbs[3*lp+ls];
    
    qs = us+1;
    if(qs == 3){
        qs = 0;
    }
    qp = up;
    q = u;
    
    cs = 2;
    cp = vindex;
    rp = vindex;
    /* cp corresponds to the last added vertex
       rp corresponds to the first added vertex */
    vindex++;
    /* loop around the edges until we arrive back at the original up vertex
       we start the search from the next neighbour vertex counterclockwise from
       neighbour lp
       qp always holds the last vertex above the plane
       lp always holds the last vertex below the plane */
    while(qp != up || qs != us){
        lp = pi->voronoi.edges[6*qp+qs];
        if(lp < 0){
            /* emergency exit */
            pi->voronoi.nvert = 0;
            return;
        }
        lw = test_vertex(&pi->voronoi.vertices[3*lp], dx, r2, &l);
        if(lw == 0){
            /* emergency exit */
            pi->voronoi.nvert = 0;
            return;
        }
        if(lw == 1){
            /* delete lp, it is still above the plane
               continue with the next vertex of lp */
            qs = pi->voronoi.edges[6*qp+3+qs]+1;
            if(qs == 3){
                qs = 0;
            }
            qp = lp;
            q = l;
            dstack_size++;
            if(dstack_size == 1001){
                /* delete stack too small, emergency exit */
                pi->voronoi.nvert = 0;
                return;
            }
            dstack[dstack_size-1] = qp;
        } else {
            /* since the previous vertex was above the plane, we have found an
               intersecting edge
               add a new vertex */
            if(q-l == 0.0f){
                error("Division by zero!");
            }
            r = q/(q-l);
            l = 1. - r;
            vindex = pi->voronoi.nvert;
            pi->voronoi.nvert++;
            if(vindex == 100){
                /* too many vertices, emergency exit */
                pi->voronoi.nvert = 0;
                return;
            }
            pi->voronoi.vertices[3*vindex+0] = pi->voronoi.vertices[3*lp+0]*r +
                                               pi->voronoi.vertices[3*qp+0]*l;
            pi->voronoi.vertices[3*vindex+1] = pi->voronoi.vertices[3*lp+1]*r +
                                               pi->voronoi.vertices[3*qp+1]*l;
            pi->voronoi.vertices[3*vindex+2] = pi->voronoi.vertices[3*lp+2]*r +
                                               pi->voronoi.vertices[3*qp+2]*l;
            /* get the index of the edge lp in qp */
            ls = pi->voronoi.edges[6*qp+3+qs];
            pi->voronoi.edges[6*vindex+0] = cp;
            pi->voronoi.edges[6*vindex+1] = lp;
            pi->voronoi.edges[6*vindex+3] = cs;
            pi->voronoi.edges[6*vindex+4] = ls;
            pi->voronoi.edges[6*lp+ls] = vindex;
            pi->voronoi.edges[6*lp+3+ls] = 1;
            pi->voronoi.edges[6*cp+cs] = vindex;
            pi->voronoi.edges[6*cp+3+cs] = 0;
            pi->voronoi.edges[6*qp+qs] = -1;
            /* ngbs */
            pi->voronoi.ngbs[3*vindex] = pj->id;
            pi->voronoi.ngbs[3*vindex+1] = pi->voronoi.ngbs[3*qp+qs];
            pi->voronoi.ngbs[3*vindex+2] = pi->voronoi.ngbs[3*lp+ls];
            /* continue with the next edge of this vertex */
            qs = qs+1;
            if(qs == 3){
                qs = 0;
            }
            cp = vindex;
            vindex++;
            cs = 2;
        }
    }
    
    pi->voronoi.edges[6*cp+cs] = rp;
    pi->voronoi.edges[6*rp] = cp;
    pi->voronoi.edges[6*cp+3+cs] = 0;
    pi->voronoi.edges[6*rp+3] = cs;
    
    /* add vertices connected to deleted ones to the delete stack
       using dstack_size (which changes inside the loop) is not a mistake:
        the newly deleted vertices can have connections that need to be
        deleted too, and this way, we make sure they are deleted as well
        this is also the reason why we reset the edges here, since otherwise
        we might get stuck in an endless loop */
    for(i = 0; i < dstack_size; i++){
        for(j = 0; j < 3; j++){
            if(pi->voronoi.edges[6*dstack[i]+j] >= 0){
                dstack_size++;
                if(dstack_size == 1001){
                    /* delete stack too small, emergency exit */
                    pi->voronoi.nvert = 0;
                    return;
                }
                dstack[dstack_size-1] = pi->voronoi.edges[6*dstack[i]+j];
                pi->voronoi.edges[6*dstack[i]+j] = -1;
                pi->voronoi.edges[6*dstack[i]+3+j] = -1;
            }
        }
    }
    
    /* check edge consistency */
    for(i = 0; i < pi->voronoi.nvert; i++){
        for(j = 0; j < 3; j++){
            k = pi->voronoi.edges[6*i+j];
            if(k >= 0 && pi->voronoi.edges[6*k] < 0){
                error("edge inconsistency!");
            }
        }
    }
    
    /* loop through the edges. If a deactivated edge is detected, remove the
       corresponding vertex by moving the next vertex (if it exists) to the
       current position. Also move its edges and change the value at the
       other endpoint of its edges */
    newnvert = pi->voronoi.nvert;
    for(i = 0; i < pi->voronoi.nvert; i++){
        if(pi->voronoi.edges[6*i] < 0){
            /* find the next valid vertex */
            j = i+1;
            while(j < pi->voronoi.nvert && pi->voronoi.edges[6*j] < 0){
                j++;
            }
            if(j == pi->voronoi.nvert){
                /* no more valid other vertex found, we can stop */
                newnvert = i;
                break;
            }
            /* copy the vertex position */
            pi->voronoi.vertices[3*i+0] = pi->voronoi.vertices[3*j+0];
            pi->voronoi.vertices[3*i+1] = pi->voronoi.vertices[3*j+1];
            pi->voronoi.vertices[3*i+2] = pi->voronoi.vertices[3*j+2];
            
            /* copy the edges */
            pi->voronoi.edges[6*i+0] = pi->voronoi.edges[6*j+0];
            pi->voronoi.edges[6*i+1] = pi->voronoi.edges[6*j+1];
            pi->voronoi.edges[6*i+2] = pi->voronoi.edges[6*j+2];
            pi->voronoi.edges[6*i+3] = pi->voronoi.edges[6*j+3];
            pi->voronoi.edges[6*i+4] = pi->voronoi.edges[6*j+4];
            pi->voronoi.edges[6*i+5] = pi->voronoi.edges[6*j+5];

            /* copy the ngbs */
            pi->voronoi.ngbs[3*i+0] = pi->voronoi.ngbs[3*j+0];
            pi->voronoi.ngbs[3*i+1] = pi->voronoi.ngbs[3*j+1];
            pi->voronoi.ngbs[3*i+2] = pi->voronoi.ngbs[3*j+2];

            /* deactivate the old edges */
            pi->voronoi.edges[6*j] = -1;
            
            /* update the vertex index in the edges for other vertices */
            pi->voronoi.edges[6*pi->voronoi.edges[6*i+0]+pi->voronoi.edges[6*i+3]] = i;
            pi->voronoi.edges[6*pi->voronoi.edges[6*i+1]+pi->voronoi.edges[6*i+4]] = i;
            pi->voronoi.edges[6*pi->voronoi.edges[6*i+2]+pi->voronoi.edges[6*i+5]] = i;
        }
    }
    
    /* resize the internal arrays */
    pi->voronoi.nvert = newnvert;

}

__attribute__((always_inline)) INLINE static float calculate_volume_tetrahedron(
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

__attribute__((always_inline)) INLINE static void calculate_centroid_tetrahedron(
    float *centroid, float *v1, float *v2, float *v3, float *v4){

    centroid[0] = 0.25f*(v1[0] + v2[0] + v3[0] + v4[0]);
    centroid[1] = 0.25f*(v1[1] + v2[1] + v3[1] + v4[1]);
    centroid[2] = 0.25f*(v1[2] + v2[2] + v3[2] + v4[2]);

}

__attribute__((always_inline)) INLINE static void calculate_cell(
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
                    tvol = calculate_volume_tetrahedron(v1, v2, v3, v4);
                    p->voronoi.volume += tvol;
                    calculate_centroid_tetrahedron(tcentroid, v1, v2, v3, v4);
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

__attribute__((always_inline)) INLINE static void calculate_faces(
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
