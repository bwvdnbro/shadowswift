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

/**
 * @brief Check if the given vertex is above or below the cutting plane
 */

__attribute__((always_inline)) INLINE static int test_vertex(
    float *v, float *dx, float r2, float *test){

    *test = -v[0]*dx[0] - v[1]*dx[1] - v[2]*dx[2] - r2;
    if(*test <= 0.){
        return -1;
    }
    if(*test > 0.){
        return 1;
    }
    return 0;

}

/**
 * @brief Initialize the cell as a cube
 */

__attribute__((always_inline)) INLINE static void voronoi_initialize(
    struct part *p){
    
    double *origin = p->x;
    
    p->voronoi.nvert = 8;

    // (0, 0, 0) -- 0
    p->voronoi.vertices[0] = -1. - origin[0];
    p->voronoi.vertices[1] = -1. - origin[1];
    p->voronoi.vertices[2] = -1. - origin[2];
    
    // (0, 0, 1)-- 1
    p->voronoi.vertices[3] = -1. - origin[0];
    p->voronoi.vertices[4] = -1. - origin[1];
    p->voronoi.vertices[5] = 2. - origin[2];
    
    // (0, 1, 0) -- 2
    p->voronoi.vertices[6] = -1. - origin[0];
    p->voronoi.vertices[7] = 2. - origin[1];
    p->voronoi.vertices[8] = -1. - origin[2];
    
    // (0, 1, 1) -- 3
    p->voronoi.vertices[9] = -1. - origin[0];
    p->voronoi.vertices[10] = 2. - origin[1];
    p->voronoi.vertices[11] = 2. - origin[2];
    
    // (1, 0, 0) -- 4
    p->voronoi.vertices[12] = 2. - origin[0];
    p->voronoi.vertices[13] = -1. - origin[1];
    p->voronoi.vertices[14] = -1. - origin[2];
    
    // (1, 0, 1) -- 5
    p->voronoi.vertices[15] = 2. - origin[0];
    p->voronoi.vertices[16] = -1. - origin[1];
    p->voronoi.vertices[17] = 2. - origin[2];
    
    // (1, 1, 0) -- 6
    p->voronoi.vertices[18] = 2. - origin[0];
    p->voronoi.vertices[19] = 2. - origin[1];
    p->voronoi.vertices[20] = -1. - origin[2];
    
    // (1, 1, 1) -- 7
    p->voronoi.vertices[21] = 2. - origin[0];
    p->voronoi.vertices[22] = 2. - origin[1];
    p->voronoi.vertices[23] = 2. - origin[2];
    
    // edges are ordered counterclockwise w.r.t. a vector pointing from the
    // cell generator to the vertex
    // (0, 0, 0) corner
    p->voronoi.edges[0] = 1;
    p->voronoi.edges[1] = 2;
    p->voronoi.edges[2] = 4;
    p->voronoi.edges[3] = 0;
    p->voronoi.edges[4] = 2;
    p->voronoi.edges[5] = 0;
    
    // (0, 0, 1) corner
    p->voronoi.edges[6] = 0;
    p->voronoi.edges[7] = 5;
    p->voronoi.edges[8] = 3;
    p->voronoi.edges[9] = 0;
    p->voronoi.edges[10] = 2;
    p->voronoi.edges[11] = 1;
    
    // (0, 1, 0) corner
    p->voronoi.edges[12] = 3;
    p->voronoi.edges[13] = 6;
    p->voronoi.edges[14] = 0;
    p->voronoi.edges[15] = 0;
    p->voronoi.edges[16] = 0;
    p->voronoi.edges[17] = 1;
    
    // (0, 1, 1) corner
    p->voronoi.edges[18] = 2;
    p->voronoi.edges[19] = 1;
    p->voronoi.edges[20] = 7;
    p->voronoi.edges[21] = 0;
    p->voronoi.edges[22] = 2;
    p->voronoi.edges[23] = 0;
    
    // (1, 0, 0) corner
    p->voronoi.edges[24] = 0;
    p->voronoi.edges[25] = 6;
    p->voronoi.edges[26] = 5;
    p->voronoi.edges[27] = 2;
    p->voronoi.edges[28] = 2;
    p->voronoi.edges[29] = 0;
    
    // (1, 0, 1) corner
    p->voronoi.edges[30] = 4;
    p->voronoi.edges[31] = 7;
    p->voronoi.edges[32] = 1;
    p->voronoi.edges[33] = 2;
    p->voronoi.edges[34] = 1;
    p->voronoi.edges[35] = 1;
    
    // (1, 1, 0) corner
    p->voronoi.edges[36] = 2;
    p->voronoi.edges[37] = 7;
    p->voronoi.edges[38] = 4;
    p->voronoi.edges[39] = 1;
    p->voronoi.edges[40] = 2;
    p->voronoi.edges[41] = 1;
    
    // (1, 1, 1) corner
    p->voronoi.edges[42] = 3;
    p->voronoi.edges[43] = 5;
    p->voronoi.edges[44] = 6;
    p->voronoi.edges[45] = 2;
    p->voronoi.edges[46] = 1;
    p->voronoi.edges[47] = 1;
}

/**
 * @brief Intersect particle pi with particle pj and adapt its Voronoi cell
 *  structure
 *
 * dx = x_i - x_j!!!
 */

__attribute__((always_inline)) INLINE static void voronoi_intersect(
    float or2, float *odx, struct part *pi, struct part *pj) {

    float r2;
    float dx[3];
    dx[0] = (pj->x[0] - pi->x[0]);
    dx[1] = (pj->x[1] - pi->x[1]);
    dx[2] = (pj->x[2] - pi->x[2]);
    
    if(dx[0] > 0.5){
        dx[0] -= 1.;
    }
    if(dx[0] < -0.5){
        dx[0] += 1.;
    }
    if(dx[1] > 0.5){
        dx[1] -= 1.;
    }
    if(dx[1] < -0.5){
        dx[1] += 1.;
    }
    if(dx[2] > 0.5){
        dx[2] -= 1.;
    }
    if(dx[2] < -0.5){
        dx[2] += 1.;
    }
    
    dx[0] *= 0.5f;
    dx[1] *= 0.5f;
    dx[2] *= 0.5f;
    
    r2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];

    float u;
    int uw = test_vertex(&pi->voronoi.vertices[0], dx, r2, &u);
    
    int up = 0;
    int lp, ls, lw, us, qp, qw, qs;
    float l, q;
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
            error("impossible case 1");
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
                    error("impossible case 2");
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
    
    /* create a delete stack for vertices that should be removed */
    int dstack[100];
    int dstack_size = 1;
    
    float r = u/(u-l);
    l = 1. - r;
    /* add new vertex */
    int vindex = pi->voronoi.nvert;
    pi->voronoi.nvert++;
    
    if(vindex == 100){
        error("too many vertices!");
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
    
    qs = us+1;
    if(qs == 3){
        qs = 0;
    }
    qp = up;
    q = u;
    
    int cs = 2;
    int cp = vindex;
    int rp = vindex;
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
            error("negative edge!");
        }
        lw = test_vertex(&pi->voronoi.vertices[3*lp], dx, r2, &l);
        if(lw == 0){
            error("case not covered 3");
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
            if(dstack_size == 101){
                error("delete stack too small");
            }
            dstack[dstack_size-1] = qp;
        } else {
            /* since the previous vertex was above the plane, we have found an
               intersecting edge
               add a new vertex */
            r = q/(q-l);
            l = 1. - r;
            vindex = pi->voronoi.nvert;
            pi->voronoi.nvert++;
            if(vindex == 100){
                error("too many vertices");
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
       using dstack_size (which changes inside the loop is not a mistake:
        the newly deleted vertices can have connections that need to be
        deleted too, and this way, we make sure they are deleted as well
        this is also the reason why we reset the edges here, since otherwise
        we might get stuck in an endless loop */
    for(int i = 0; i < dstack_size; i++){
        for(int j = 0; j < 3; j++){
            if(pi->voronoi.edges[6*dstack[i]+j] >= 0){
                dstack_size++;
                if(dstack_size == 101){
                    error("delete stack too small");
                }
                dstack[dstack_size-1] = pi->voronoi.edges[6*dstack[i]+j];
                pi->voronoi.edges[6*dstack[i]+j] = -1;
                pi->voronoi.edges[6*dstack[i]+3+j] = -1;
            }
        }
    }
    
    /* check edge consistency */
    for(int i = 0; i < pi->voronoi.nvert; i++){
        for(int j = 0; j < 3; j++){
            int k = pi->voronoi.edges[6*i+j];
            if(k >= 0 && pi->voronoi.edges[6*k] < 0){
                error("edge inconsistency!");
            }
        }
    }
    
    /* loop through the edges. If a deactivated edge is detected, remove the
       corresponding vertex by moving the next vertex (if it exists) to the
       current position. Also move its edges and change the value at the
       other endpoint of its edges */
    int newnvert = pi->voronoi.nvert;
    for(int i = 0; i < pi->voronoi.nvert; i++){
        if(pi->voronoi.edges[6*i] < 0){
            /* find the next valid vertex */
            int j = i+1;
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

#endif /* SWIFT_VORONOI_H */
