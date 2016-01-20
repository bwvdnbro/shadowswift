/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <float.h>
#include <limits.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* This object's header. */
#include "space.h"

/* Local headers. */
#include "atomic.h"
#include "engine.h"
#include "error.h"
#include "kernel.h"
#include "lock.h"
#include "runner.h"

/* Shared sort structure. */
struct parallel_sort space_sort_struct;

/* Split size. */
int space_splitsize = space_splitsize_default;
int space_subsize = space_subsize_default;
int space_maxsize = space_maxsize_default;

/* Map shift vector to sortlist. */
const int sortlistID[27] = {
    /* ( -1 , -1 , -1 ) */ 0,
    /* ( -1 , -1 ,  0 ) */ 1,
    /* ( -1 , -1 ,  1 ) */ 2,
    /* ( -1 ,  0 , -1 ) */ 3,
    /* ( -1 ,  0 ,  0 ) */ 4,
    /* ( -1 ,  0 ,  1 ) */ 5,
    /* ( -1 ,  1 , -1 ) */ 6,
    /* ( -1 ,  1 ,  0 ) */ 7,
    /* ( -1 ,  1 ,  1 ) */ 8,
    /* (  0 , -1 , -1 ) */ 9,
    /* (  0 , -1 ,  0 ) */ 10,
    /* (  0 , -1 ,  1 ) */ 11,
    /* (  0 ,  0 , -1 ) */ 12,
    /* (  0 ,  0 ,  0 ) */ 0,
    /* (  0 ,  0 ,  1 ) */ 12,
    /* (  0 ,  1 , -1 ) */ 11,
    /* (  0 ,  1 ,  0 ) */ 10,
    /* (  0 ,  1 ,  1 ) */ 9,
    /* (  1 , -1 , -1 ) */ 8,
    /* (  1 , -1 ,  0 ) */ 7,
    /* (  1 , -1 ,  1 ) */ 6,
    /* (  1 ,  0 , -1 ) */ 5,
    /* (  1 ,  0 ,  0 ) */ 4,
    /* (  1 ,  0 ,  1 ) */ 3,
    /* (  1 ,  1 , -1 ) */ 2,
    /* (  1 ,  1 ,  0 ) */ 1,
    /* (  1 ,  1 ,  1 ) */ 0};

/**
 * @brief Get the shift-id of the given pair of cells, swapping them
 *      if need be.
 *
 * @param s The space
 * @param ci Pointer to first #cell.
 * @param cj Pointer second #cell.
 * @param shift Vector from ci to cj.
 *
 * @return The shift ID and set shift, may or may not swap ci and cj.
 */

int space_getsid(struct space *s, struct cell **ci, struct cell **cj,
                 double *shift) {

  int k, sid = 0, periodic = s->periodic;
  struct cell *temp;
  double dx[3];

  /* Get the relative distance between the pairs, wrapping. */
  for (k = 0; k < 3; k++) {
    dx[k] = (*cj)->loc[k] - (*ci)->loc[k];
    if (periodic && dx[k] < -s->dim[k] / 2)
      shift[k] = s->dim[k];
    else if (periodic && dx[k] > s->dim[k] / 2)
      shift[k] = -s->dim[k];
    else
      shift[k] = 0.0;
    dx[k] += shift[k];
  }

  /* Get the sorting index. */
  for (k = 0; k < 3; k++)
    sid = 3 * sid + ((dx[k] < 0.0) ? 0 : ((dx[k] > 0.0) ? 2 : 1));

  /* Switch the cells around? */
  if (runner_flip[sid]) {
    temp = *ci;
    *ci = *cj;
    *cj = temp;
    for (k = 0; k < 3; k++) shift[k] = -shift[k];
  }
  sid = sortlistID[sid];

  /* Return the sort ID. */
  return sid;
}

/**
 * @brief Recursively dismantle a cell tree.
 *
 */

void space_rebuild_recycle(struct space *s, struct cell *c) {

  int k;

  if (c->split)
    for (k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) {
        space_rebuild_recycle(s, c->progeny[k]);
        space_recycle(s, c->progeny[k]);
        c->progeny[k] = NULL;
      }
}

/**
 * @brief Re-build the cell grid.
 *
 * @param s The #space.
 * @param cell_max Maximum cell edge length.
 * @param verbose Print messages to stdout or not.
 */

void space_regrid(struct space *s, double cell_max, int verbose) {

  float h_max = s->cell_min / kernel_gamma / space_stretch, dmin;
  int i, j, k, cdim[3], nr_parts = s->nr_parts;
  struct cell *restrict c;
  // ticks tic;

  /* Run through the parts and get the current h_max. */
  // tic = getticks();
  if (s->cells != NULL) {
    for (k = 0; k < s->nr_cells; k++) {
      if (s->cells[k].h_max > h_max) h_max = s->cells[k].h_max;
    }
  } else {
    for (k = 0; k < nr_parts; k++) {
      if (s->parts[k].h > h_max) h_max = s->parts[k].h;
    }
    s->h_max = h_max;
  }

/* If we are running in parallel, make sure everybody agrees on
   how large the largest cell should be. */
#ifdef WITH_MPI
  {
    float buff;
    if (MPI_Allreduce(&h_max, &buff, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD) !=
        MPI_SUCCESS)
      error("Failed to aggregate the rebuild flag across nodes.");
    h_max = buff;
  }
#endif
  if (verbose) message("h_max is %.3e (cell_max=%.3e).", h_max, cell_max);

  /* Get the new putative cell dimensions. */
  for (k = 0; k < 3; k++)
    cdim[k] =
        floor(s->dim[k] / fmax(h_max * kernel_gamma * space_stretch, cell_max));

  /* Check if we have enough cells for periodicity. */
  if (s->periodic && (cdim[0] < 3 || cdim[1] < 3 || cdim[2] < 3))
    error(
        "Must have at least 3 cells in each spatial dimension when periodicity "
        "is switched on.");

/* In MPI-Land, we're not allowed to change the top-level cell size. */
#ifdef WITH_MPI
  if (cdim[0] < s->cdim[0] || cdim[1] < s->cdim[1] || cdim[2] < s->cdim[2])
    error("Root-level change of cell size not allowed.");
#endif

  /* Do we need to re-build the upper-level cells? */
  // tic = getticks();
  if (s->cells == NULL || cdim[0] < s->cdim[0] || cdim[1] < s->cdim[1] ||
      cdim[2] < s->cdim[2]) {

    /* Free the old cells, if they were allocated. */
    if (s->cells != NULL) {
      for (k = 0; k < s->nr_cells; k++) {
        space_rebuild_recycle(s, &s->cells[k]);
        if (s->cells[k].sort != NULL) free(s->cells[k].sort);
      }
      free(s->cells);
      s->maxdepth = 0;
    }

    /* Set the new cell dimensions only if smaller. */
    for (k = 0; k < 3; k++) {
      s->cdim[k] = cdim[k];
      s->h[k] = s->dim[k] / cdim[k];
      s->ih[k] = 1.0 / s->h[k];
    }
    dmin = fminf(s->h[0], fminf(s->h[1], s->h[2]));

    /* Allocate the highest level of cells. */
    s->tot_cells = s->nr_cells = cdim[0] * cdim[1] * cdim[2];
    if (posix_memalign((void *)&s->cells, 64,
                       s->nr_cells * sizeof(struct cell)) != 0)
      error("Failed to allocate cells.");
    bzero(s->cells, s->nr_cells * sizeof(struct cell));
    for (k = 0; k < s->nr_cells; k++)
      if (lock_init(&s->cells[k].lock) != 0) error("Failed to init spinlock.");

    /* Set the cell location and sizes. */
    for (i = 0; i < cdim[0]; i++)
      for (j = 0; j < cdim[1]; j++)
        for (k = 0; k < cdim[2]; k++) {
          c = &s->cells[cell_getid(cdim, i, j, k)];
          c->loc[0] = i * s->h[0];
          c->loc[1] = j * s->h[1];
          c->loc[2] = k * s->h[2];
          c->h[0] = s->h[0];
          c->h[1] = s->h[1];
          c->h[2] = s->h[2];
          c->dmin = dmin;
          c->depth = 0;
          c->count = 0;
          c->gcount = 0;
          c->super = c;
          lock_init(&c->lock);
        }

    /* Be verbose about the change. */
    if (verbose)
      message("set cell dimensions to [ %i %i %i ].", cdim[0], cdim[1],
              cdim[2]);
    fflush(stdout);

  } /* re-build upper-level cells? */
  // message( "rebuilding upper-level cells took %.3f ms." , (double)(getticks()
  // - tic) / CPU_TPS * 1000 );

  /* Otherwise, just clean up the cells. */
  else {

    /* Free the old cells, if they were allocated. */
    for (k = 0; k < s->nr_cells; k++) {
      space_rebuild_recycle(s, &s->cells[k]);
      s->cells[k].sorts = NULL;
      s->cells[k].nr_tasks = 0;
      s->cells[k].nr_density = 0;
      s->cells[k].nr_force = 0;
      s->cells[k].density = NULL;
      s->cells[k].force = NULL;
      s->cells[k].dx_max = 0.0f;
      s->cells[k].sorted = 0;
      s->cells[k].count = 0;
      s->cells[k].gcount = 0;
      s->cells[k].kick1 = NULL;
      s->cells[k].kick2 = NULL;
      s->cells[k].super = &s->cells[k];
    }
    s->maxdepth = 0;
  }
}

/**
 * @brief Re-build the cells as well as the tasks.
 *
 * @param s The #space in which to update the cells.
 * @param cell_max Maximal cell size.
 * @param verbose Print messages to stdout or not
 *
 */

void space_rebuild(struct space *s, double cell_max, int verbose) {

  int j, k, cdim[3], nr_parts = s->nr_parts, nr_gparts = s->nr_gparts;
  struct cell *restrict c, *restrict cells;
  struct part *restrict p;
  int *ind;
  double ih[3], dim[3];
  // ticks tic;

  /* Be verbose about this. */
  // message( "re)building space..." ); fflush(stdout);

  /* Re-grid if necessary, or just re-set the cell data. */
  space_regrid(s, cell_max, verbose);
  cells = s->cells;

  /* Run through the particles and get their cell index. */
  // tic = getticks();
  const int ind_size = s->size_parts;
  if ((ind = (int *)malloc(sizeof(int) * ind_size)) == NULL)
    error("Failed to allocate temporary particle indices.");
  ih[0] = s->ih[0];
  ih[1] = s->ih[1];
  ih[2] = s->ih[2];
  dim[0] = s->dim[0];
  dim[1] = s->dim[1];
  dim[2] = s->dim[2];
  cdim[0] = s->cdim[0];
  cdim[1] = s->cdim[1];
  cdim[2] = s->cdim[2];
  for (k = 0; k < nr_parts; k++) {
    p = &s->parts[k];
    for (j = 0; j < 3; j++)
      if (p->x[j] < 0.0)
        p->x[j] += dim[j];
      else if (p->x[j] >= dim[j])
        p->x[j] -= dim[j];
    ind[k] =
        cell_getid(cdim, p->x[0] * ih[0], p->x[1] * ih[1], p->x[2] * ih[2]);
    cells[ind[k]].count++;
  }
// message( "getting particle indices took %.3f ms." , (double)(getticks() -
// tic) / CPU_TPS * 1000 );

#ifdef WITH_MPI
  /* Move non-local parts to the end of the list. */
  int nodeID = s->e->nodeID;
  for (k = 0; k < nr_parts; k++)
    if (cells[ind[k]].nodeID != nodeID) {
      cells[ind[k]].count -= 1;
      nr_parts -= 1;
      struct part tp = s->parts[k];
      s->parts[k] = s->parts[nr_parts];
      s->parts[nr_parts] = tp;
      struct xpart txp = s->xparts[k];
      s->xparts[k] = s->xparts[nr_parts];
      s->xparts[nr_parts] = txp;
      int t = ind[k];
      ind[k] = ind[nr_parts];
      ind[nr_parts] = t;
    }

  /* Exchange the strays, note that this potentially re-allocates
     the parts arrays. */
  s->nr_parts =
      nr_parts + engine_exchange_strays(s->e, nr_parts, &ind[nr_parts],
                                        s->nr_parts - nr_parts);

  /* Re-allocate the index array if needed.. */
  if (s->nr_parts > ind_size) {
    int *ind_new;
    if ((ind_new = (int *)malloc(sizeof(int) * s->nr_parts)) == NULL)
      error("Failed to allocate temporary particle indices.");
    memcpy(ind_new, ind, sizeof(int) * nr_parts);
    free(ind);
    ind = ind_new;
  }

  /* Assign each particle to its cell. */
  for (k = nr_parts; k < s->nr_parts; k++) {
    p = &s->parts[k];
    ind[k] =
        cell_getid(cdim, p->x[0] * ih[0], p->x[1] * ih[1], p->x[2] * ih[2]);
    cells[ind[k]].count += 1;
    /* if ( cells[ ind[k] ].nodeID != nodeID )
        error( "Received part that does not belong to me (nodeID=%i)." , cells[
       ind[k] ].nodeID ); */
  }
  nr_parts = s->nr_parts;
#endif

  /* Sort the parts according to their cells. */
  // tic = getticks();
  space_parts_sort(s, ind, nr_parts, 0, s->nr_cells - 1);
  // message( "parts_sort took %.3f ms." , (double)(getticks() - tic) / CPU_TPS
  // * 1000 );

  /* Re-link the gparts. */
  for (k = 0; k < nr_parts; k++)
    if (s->parts[k].gpart != NULL) s->parts[k].gpart->part = &s->parts[k];

  /* Verify space_sort_struct. */
  /* for ( k = 1 ; k < nr_parts ; k++ ) {
      if ( ind[k-1] > ind[k] ) {
          error( "Sort failed!" );
          }
      else if ( ind[k] != cell_getid( cdim , parts[k].x[0]*ih[0] ,
     parts[k].x[1]*ih[1] , parts[k].x[2]*ih[2] ) )
          error( "Incorrect indices!" );
      } */

  /* We no longer need the indices as of here. */
  free(ind);

  /* Run through the gravity particles and get their cell index. */
  // tic = getticks();
  if ((ind = (int *)malloc(sizeof(int) * s->size_gparts)) == NULL)
    error("Failed to allocate temporary particle indices.");
  for (k = 0; k < nr_gparts; k++) {
    struct gpart *gp = &s->gparts[k];
    for (j = 0; j < 3; j++)
      if (gp->x[j] < 0.0)
        gp->x[j] += dim[j];
      else if (gp->x[j] >= dim[j])
        gp->x[j] -= dim[j];
    ind[k] =
        cell_getid(cdim, gp->x[0] * ih[0], gp->x[1] * ih[1], gp->x[2] * ih[2]);
    cells[ind[k]].gcount++;
  }
  // message( "getting particle indices took %.3f ms." , (double)(getticks() -
  // tic) / CPU_TPS * 1000 );

  /* TODO: Here we should exchange the gparts as well! */

  /* Sort the parts according to their cells. */
  // tic = getticks();
  gparts_sort(s->gparts, ind, nr_gparts, 0, s->nr_cells - 1);
  // message( "gparts_sort took %.3f ms." , (double)(getticks() - tic) / CPU_TPS
  // * 1000 );

  /* Re-link the parts. */
  for (k = 0; k < nr_gparts; k++)
    if (s->gparts[k].id > 0) s->gparts[k].part->gpart = &s->gparts[k];

  /* We no longer need the indices as of here. */
  free(ind);

  /* Hook the cells up to the parts. */
  // tic = getticks();
  struct part *finger = s->parts;
  struct xpart *xfinger = s->xparts;
  struct gpart *gfinger = s->gparts;
  for (k = 0; k < s->nr_cells; k++) {
    c = &cells[k];
    c->parts = finger;
    c->xparts = xfinger;
    c->gparts = gfinger;
    finger = &finger[c->count];
    xfinger = &xfinger[c->count];
    gfinger = &gfinger[c->gcount];
  }
  // message( "hooking up cells took %.3f ms." , (double)(getticks() - tic) /
  // CPU_TPS * 1000 );

  /* At this point, we have the upper-level cells, old or new. Now make
     sure that the parts in each cell are ok. */
  // tic = getticks();
  // for (k = 0; k < s->nr_cells; k++) space_split(s, &cells[k]);
  for (k = 0; k < s->nr_cells; k++)
    scheduler_addtask(&s->e->sched, task_type_split_cell, task_subtype_none,
                      k, 0, &cells[k], NULL, 0);
  engine_launch(s->e, s->e->nr_threads, 1 << task_type_split_cell);

  // message( "space_split took %.3f ms." , (double)(getticks() - tic) / CPU_TPS
  // * 1000 );
}

/**
 * @brief Sort the particles and condensed particles according to the given
 *indices.
 *
 * @param s The #space.
 * @param ind The indices with respect to which the parts are sorted.
 * @param N The number of parts
 * @param min Lowest index.
 * @param max highest index.
 */

void space_parts_sort(struct space *s, int *ind, int N, int min, int max) {
  // Populate the global parallel_sort structure with the input data.
  space_sort_struct.parts = s->parts;
  space_sort_struct.xparts = s->xparts;
  space_sort_struct.ind = ind;
  space_sort_struct.stack_size = 2 * (max - min + 1) + 10 + s->e->nr_threads;
  if ((space_sort_struct.stack = malloc(sizeof(struct qstack) *
                                        space_sort_struct.stack_size)) == NULL)
    error("Failed to allocate sorting stack.");
  for (int i = 0; i < space_sort_struct.stack_size; i++)
    space_sort_struct.stack[i].ready = 0;

  // Add the first interval.
  space_sort_struct.stack[0].i = 0;
  space_sort_struct.stack[0].j = N - 1;
  space_sort_struct.stack[0].min = min;
  space_sort_struct.stack[0].max = max;
  space_sort_struct.stack[0].ready = 1;
  space_sort_struct.first = 0;
  space_sort_struct.last = 1;
  space_sort_struct.waiting = 1;

  // Launch the sorting tasks.
  engine_launch(s->e, s->e->nr_threads, (1 << task_type_psort));

  /* Verify space_sort_struct. */
  /* for (int i = 1; i < N; i++)
    if (ind[i - 1] > ind[i])
      error("Sorting failed (ind[%i]=%i,ind[%i]=%i), min=%i, max=%i.", i - 1, ind[i - 1], i,
            ind[i], min, max);
  message("Sorting succeeded."); */

  // Clean up.
  free(space_sort_struct.stack);
}

void space_do_parts_sort() {

  /* Pointers to the sorting data. */
  int *ind = space_sort_struct.ind;
  struct part *parts = space_sort_struct.parts;
  struct xpart *xparts = space_sort_struct.xparts;

  /* Main loop. */
  while (space_sort_struct.waiting) {

    /* Grab an interval off the queue. */
    int qid =
        atomic_inc(&space_sort_struct.first) % space_sort_struct.stack_size;

    /* Wait for the entry to be ready, or for the sorting do be done. */
    while (!space_sort_struct.stack[qid].ready)
      if (!space_sort_struct.waiting) return;
      
    /* Get the stack entry. */
    int i = space_sort_struct.stack[qid].i;
    int j = space_sort_struct.stack[qid].j;
    int min = space_sort_struct.stack[qid].min;
    int max = space_sort_struct.stack[qid].max;
    space_sort_struct.stack[qid].ready = 0;

    /* Loop over sub-intervals. */
    while (1) {

      /* Bring beer. */
      const int pivot = (min + max) / 2;
      /* message("Working on interval [%i,%i] with min=%i, max=%i, pivot=%i.",
              i, j, min, max, pivot); */

      /* One pass of QuickSort's partitioning. */
      int ii = i;
      int jj = j;
      while (ii < jj) {
        while (ii <= j && ind[ii] <= pivot) ii++;
        while (jj >= i && ind[jj] > pivot) jj--;
        if (ii < jj) {
          int temp_i = ind[ii];
          ind[ii] = ind[jj];
          ind[jj] = temp_i;
          struct part temp_p = parts[ii];
          parts[ii] = parts[jj];
          parts[jj] = temp_p;
          struct xpart temp_xp = xparts[ii];
          xparts[ii] = xparts[jj];
          xparts[jj] = temp_xp;
        }
      }

      /* Verify space_sort_struct. */
      /* for (int k = i; k <= jj; k++)
        if (ind[k] > pivot) {
          message("sorting failed at k=%i, ind[k]=%i, pivot=%i, i=%i, j=%i.", k,
                  ind[k], pivot, i, j);
          error("Partition failed (<=pivot).");
        }
      for (int k = jj + 1; k <= j; k++)
        if (ind[k] <= pivot) {
          message("sorting failed at k=%i, ind[k]=%i, pivot=%i, i=%i, j=%i.", k,
                  ind[k], pivot, i, j);
          error("Partition failed (>pivot).");
        } */

      /* Split-off largest interval. */
      if (jj - i > j - jj + 1) {

        /* Recurse on the left? */
        if (jj > i && pivot > min) {
          qid = atomic_inc(&space_sort_struct.last) %
                space_sort_struct.stack_size;
          while (space_sort_struct.stack[qid].ready);
          space_sort_struct.stack[qid].i = i;
          space_sort_struct.stack[qid].j = jj;
          space_sort_struct.stack[qid].min = min;
          space_sort_struct.stack[qid].max = pivot;
          if (atomic_inc(&space_sort_struct.waiting) >=
              space_sort_struct.stack_size)
            error("Qstack overflow.");
          space_sort_struct.stack[qid].ready = 1;
        }

        /* Recurse on the right? */
        if (jj + 1 < j && pivot + 1 < max) {
          i = jj + 1;
          min = pivot + 1;
        } else
          break;

      } else {

        /* Recurse on the right? */
        if (pivot + 1 < max) {
          qid = atomic_inc(&space_sort_struct.last) %
                space_sort_struct.stack_size;
          while (space_sort_struct.stack[qid].ready);
          space_sort_struct.stack[qid].i = jj + 1;
          space_sort_struct.stack[qid].j = j;
          space_sort_struct.stack[qid].min = pivot + 1;
          space_sort_struct.stack[qid].max = max;
          if (atomic_inc(&space_sort_struct.waiting) >=
              space_sort_struct.stack_size)
            error("Qstack overflow.");
          space_sort_struct.stack[qid].ready = 1;
        }

        /* Recurse on the left? */
        if (jj > i && pivot > min) {
          j = jj;
          max = pivot;
        } else
          break;
      }

    } /* loop over sub-intervals. */

    atomic_dec(&space_sort_struct.waiting);

  } /* main loop. */
}

void gparts_sort(struct gpart *gparts, int *ind, int N, int min, int max) {

  struct qstack {
    volatile int i, j, min, max;
    volatile int ready;
  };
  struct qstack *qstack;
  int qstack_size = 2 * (max - min) + 10;
  volatile unsigned int first, last, waiting;

  int pivot;
  int i, ii, j, jj, temp_i, qid;
  struct gpart temp_p;

  /* for ( int k = 0 ; k < N ; k++ )
      if ( ind[k] > max || ind[k] < min )
          error( "ind[%i]=%i is not in [%i,%i]." , k , ind[k] , min , max ); */

  /* Allocate the stack. */
  if ((qstack = malloc(sizeof(struct qstack) * qstack_size)) == NULL)
    error("Failed to allocate qstack.");

  /* Init the interval stack. */
  qstack[0].i = 0;
  qstack[0].j = N - 1;
  qstack[0].min = min;
  qstack[0].max = max;
  qstack[0].ready = 1;
  for (i = 1; i < qstack_size; i++) qstack[i].ready = 0;
  first = 0;
  last = 1;
  waiting = 1;

  /* Main loop. */
  while (waiting > 0) {

    /* Grab an interval off the queue. */
    qid = (first++) % qstack_size;

    /* Get the stack entry. */
    i = qstack[qid].i;
    j = qstack[qid].j;
    min = qstack[qid].min;
    max = qstack[qid].max;
    qstack[qid].ready = 0;

    /* Loop over sub-intervals. */
    while (1) {

      /* Bring beer. */
      pivot = (min + max) / 2;

      /* One pass of QuickSort's partitioning. */
      ii = i;
      jj = j;
      while (ii < jj) {
        while (ii <= j && ind[ii] <= pivot) ii++;
        while (jj >= i && ind[jj] > pivot) jj--;
        if (ii < jj) {
          temp_i = ind[ii];
          ind[ii] = ind[jj];
          ind[jj] = temp_i;
          temp_p = gparts[ii];
          gparts[ii] = gparts[jj];
          gparts[jj] = temp_p;
        }
      }

      /* Verify space_sort_struct. */
      /* for ( int k = i ; k <= jj ; k++ )
         if ( ind[k] > pivot ) {
         message( "sorting failed at k=%i, ind[k]=%i, pivot=%i, i=%i, j=%i,
         N=%i." , k , ind[k] , pivot , i , j , N );
         error( "Partition failed (<=pivot)." );
         }
         for ( int k = jj+1 ; k <= j ; k++ )
         if ( ind[k] <= pivot ) {
         message( "sorting failed at k=%i, ind[k]=%i, pivot=%i, i=%i, j=%i,
         N=%i." , k , ind[k] , pivot , i , j , N );
         error( "Partition failed (>pivot)." );
         } */

      /* Split-off largest interval. */
      if (jj - i > j - jj + 1) {

        /* Recurse on the left? */
        if (jj > i && pivot > min) {
          qid = (last++) % qstack_size;
          qstack[qid].i = i;
          qstack[qid].j = jj;
          qstack[qid].min = min;
          qstack[qid].max = pivot;
          qstack[qid].ready = 1;
          if ((waiting++) >= qstack_size) error("Qstack overflow.");
        }

        /* Recurse on the right? */
        if (jj + 1 < j && pivot + 1 < max) {
          i = jj + 1;
          min = pivot + 1;
        } else
          break;

      } else {

        /* Recurse on the right? */
        if (pivot + 1 < max) {
          qid = (last++) % qstack_size;
          qstack[qid].i = jj + 1;
          qstack[qid].j = j;
          qstack[qid].min = pivot + 1;
          qstack[qid].max = max;
          qstack[qid].ready = 1;
          if ((waiting++) >= qstack_size) error("Qstack overflow.");
        }

        /* Recurse on the left? */
        if (jj > i && pivot > min) {
          j = jj;
          max = pivot;
        } else
          break;
      }

    } /* loop over sub-intervals. */

    waiting--;

  } /* main loop. */

  /* Verify space_sort_struct. */
  /* for ( i = 1 ; i < N ; i++ )
      if ( ind[i-1] > ind[i] )
          error( "Sorting failed (ind[%i]=%i,ind[%i]=%i)." , i-1 , ind[i-1] , i
     , ind[i] ); */

  /* Clean up. */
  free(qstack);
}

/**
 * @brief Mapping function to free the sorted indices buffers.
 */

void space_map_clearsort(struct cell *c, void *data) {

  if (c->sort != NULL) {
    free(c->sort);
    c->sort = NULL;
  }
}

/**
 * @brief Map a function to all particles in a cell recursively.
 *
 * @param c The #cell we are working in.
 * @param fun Function pointer to apply on the cells.
 * @param data Data passed to the function fun.
 */

static void rec_map_parts(struct cell *c,
                          void (*fun)(struct part *p, struct cell *c,
                                      void *data),
                          void *data) {

  int k;

  /* No progeny? */
  if (!c->split)
    for (k = 0; k < c->count; k++) fun(&c->parts[k], c, data);

  /* Otherwise, recurse. */
  else
    for (k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) rec_map_parts(c->progeny[k], fun, data);
}

/**
 * @brief Map a function to all particles in a space.
 *
 * @param s The #space we are working in.
 * @param fun Function pointer to apply on the cells.
 * @param data Data passed to the function fun.
 */

void space_map_parts(struct space *s,
                     void (*fun)(struct part *p, struct cell *c, void *data),
                     void *data) {

  int cid = 0;

  /* Call the recursive function on all higher-level cells. */
  for (cid = 0; cid < s->nr_cells; cid++)
    rec_map_parts(&s->cells[cid], fun, data);
}

/**
 * @brief Map a function to all particles in a cell recursively.
 *
 * @param c The #cell we are working in.
 * @param full Map to all cells, including cells with sub-cells.
 * @param fun Function pointer to apply on the cells.
 * @param data Data passed to the function fun.
 */

static void rec_map_cells_post(struct cell *c, int full,
                               void (*fun)(struct cell *c, void *data),
                               void *data) {

  int k;

  /* Recurse. */
  if (c->split)
    for (k = 0; k < 8; k++)
      if (c->progeny[k] != NULL)
        rec_map_cells_post(c->progeny[k], full, fun, data);

  /* No progeny? */
  if (full || !c->split) fun(c, data);
}

/**
 * @brief Map a function to all particles in a aspace.
 *
 * @param s The #space we are working in.
 * @param full Map to all cells, including cells with sub-cells.
 * @param fun Function pointer to apply on the cells.
 * @param data Data passed to the function fun.
 */

void space_map_cells_post(struct space *s, int full,
                          void (*fun)(struct cell *c, void *data), void *data) {

  int cid = 0;

  /* Call the recursive function on all higher-level cells. */
  for (cid = 0; cid < s->nr_cells; cid++)
    rec_map_cells_post(&s->cells[cid], full, fun, data);
}

static void rec_map_cells_pre(struct cell *c, int full,
                              void (*fun)(struct cell *c, void *data),
                              void *data) {

  int k;

  /* No progeny? */
  if (full || !c->split) fun(c, data);

  /* Recurse. */
  if (c->split)
    for (k = 0; k < 8; k++)
      if (c->progeny[k] != NULL)
        rec_map_cells_pre(c->progeny[k], full, fun, data);
}

void space_map_cells_pre(struct space *s, int full,
                         void (*fun)(struct cell *c, void *data), void *data) {

  int cid = 0;

  /* Call the recursive function on all higher-level cells. */
  for (cid = 0; cid < s->nr_cells; cid++)
    rec_map_cells_pre(&s->cells[cid], full, fun, data);
}

/**
 * @brief Split cells that contain too many particles.
 *
 * @param s The #space we are working in.
 * @param c The #cell under consideration.
 */

void space_split(struct space *s, struct cell *c) {

  int k, count = c->count, gcount = c->gcount, maxdepth = 0;
  float h, h_max = 0.0f, dt, dt_min = c->parts[0].dt, dt_max = dt_min;
  struct cell *temp;
  struct part *p, *parts = c->parts;
  struct xpart *xp, *xparts = c->xparts;

  /* Check the depth. */
  if (c->depth > s->maxdepth) s->maxdepth = c->depth;

  /* Split or let it be? */
  if (count > space_splitsize || gcount > space_splitsize) {

    /* No longer just a leaf. */
    c->split = 1;

    /* Create the cell's progeny. */
    for (k = 0; k < 8; k++) {
      temp = space_getcell(s);
      temp->count = 0;
      temp->gcount = 0;
      temp->loc[0] = c->loc[0];
      temp->loc[1] = c->loc[1];
      temp->loc[2] = c->loc[2];
      temp->h[0] = c->h[0] / 2;
      temp->h[1] = c->h[1] / 2;
      temp->h[2] = c->h[2] / 2;
      temp->dmin = c->dmin / 2;
      if (k & 4) temp->loc[0] += temp->h[0];
      if (k & 2) temp->loc[1] += temp->h[1];
      if (k & 1) temp->loc[2] += temp->h[2];
      temp->depth = c->depth + 1;
      temp->split = 0;
      temp->h_max = 0.0;
      temp->dx_max = 0.0;
      temp->nodeID = c->nodeID;
      temp->parent = c;
      c->progeny[k] = temp;
    }

    /* Split the cell data. */
    cell_split(c);

    /* Remove any progeny with zero parts. */
    for (k = 0; k < 8; k++)
      if (c->progeny[k]->count == 0 && c->progeny[k]->gcount == 0) {
        space_recycle(s, c->progeny[k]);
        c->progeny[k] = NULL;
      } else {
        space_split(s, c->progeny[k]);
        h_max = fmaxf(h_max, c->progeny[k]->h_max);
        dt_min = fminf(dt_min, c->progeny[k]->dt_min);
        dt_max = fmaxf(dt_max, c->progeny[k]->dt_max);
        if (c->progeny[k]->maxdepth > maxdepth)
          maxdepth = c->progeny[k]->maxdepth;
      }

    /* Set the values for this cell. */
    c->h_max = h_max;
    c->dt_min = dt_min;
    c->dt_max = dt_max;
    c->maxdepth = maxdepth;

  }

  /* Otherwise, collect the data for this cell. */
  else {

    /* Clear the progeny. */
    bzero(c->progeny, sizeof(struct cell *) * 8);
    c->split = 0;
    c->maxdepth = c->depth;

    /* Get dt_min/dt_max. */

    for (k = 0; k < count; k++) {
      p = &parts[k];
      xp = &xparts[k];
      xp->x_old[0] = p->x[0];
      xp->x_old[1] = p->x[1];
      xp->x_old[2] = p->x[2];
      dt = p->dt;
      h = p->h;
      if (h > h_max) h_max = h;
      if (dt < dt_min) dt_min = dt;
      if (dt > dt_max) dt_max = dt;
    }
    c->h_max = h_max;
    c->dt_min = dt_min;
    c->dt_max = dt_max;
  }

  /* Set ownership according to the start of the parts array. */
  c->owner = ((c->parts - s->parts) % s->nr_parts) * s->nr_queues / s->nr_parts;
}

/**
 * @brief Return a used cell to the cell buffer.
 *
 * @param s The #space.
 * @param c The #cell.
 */

void space_recycle(struct space *s, struct cell *c) {

  /* Lock the space. */
  lock_lock(&s->lock);

  /* Clear the cell. */
  if (lock_destroy(&c->lock) != 0) error("Failed to destroy spinlock.");

  /* Clear this cell's sort arrays. */
  if (c->sort != NULL) free(c->sort);

  /* Clear the cell data. */
  bzero(c, sizeof(struct cell));

  /* Hook this cell into the buffer. */
  c->next = s->cells_new;
  s->cells_new = c;
  s->tot_cells -= 1;

  /* Unlock the space. */
  lock_unlock_blind(&s->lock);
}

/**
 * @brief Get a new empty cell.
 *
 * @param s The #space.
 */

struct cell *space_getcell(struct space *s) {

  struct cell *c;
  int k;

  /* Lock the space. */
  lock_lock(&s->lock);

  /* Is the buffer empty? */
  if (s->cells_new == NULL) {
    if (posix_memalign((void *)&s->cells_new, 64,
                       space_cellallocchunk * sizeof(struct cell)) != 0)
      error("Failed to allocate more cells.");
    bzero(s->cells_new, space_cellallocchunk * sizeof(struct cell));
    for (k = 0; k < space_cellallocchunk - 1; k++)
      s->cells_new[k].next = &s->cells_new[k + 1];
    s->cells_new[space_cellallocchunk - 1].next = NULL;
  }

  /* Pick off the next cell. */
  c = s->cells_new;
  s->cells_new = c->next;
  s->tot_cells += 1;

  /* Unlock the space. */
  lock_unlock_blind(&s->lock);

  /* Init some things in the cell. */
  bzero(c, sizeof(struct cell));
  c->nodeID = -1;
  if (lock_init(&c->lock) != 0 || lock_init(&c->glock) != 0)
    error("Failed to initialize cell spinlocks.");

  return c;
}

/**
 * @brief Split the space into cells given the array of particles.
 *
 * @param s The #space to initialize.
 * @param dim Spatial dimensions of the domain.
 * @param parts Pointer to an array of #part.
 * @param N The number of parts in the space.
 * @param periodic flag whether the domain is periodic or not.
 * @param h_max The maximal interaction radius.
 * @param verbose Print messages to stdout or not
 *
 * Makes a grid of edge length > r_max and fills the particles
 * into the respective cells. Cells containing more than #space_splitsize
 * parts with a cutoff below half the cell width are then split
 * recursively.
 */

void space_init(struct space *s, double dim[3], struct part *parts, int N,
                int periodic, double h_max, int verbose) {

  /* Store everything in the space. */
  s->dim[0] = dim[0];
  s->dim[1] = dim[1];
  s->dim[2] = dim[2];
  s->periodic = periodic;
  s->nr_parts = N;
  s->size_parts = N;
  s->parts = parts;
  s->cell_min = h_max;
  s->nr_queues = 1;
  s->size_parts_foreign = 0;

  /* Check that all the particle positions are reasonable, wrap if periodic. */
  if (periodic) {
    for (int k = 0; k < N; k++)
      for (int j = 0; j < 3; j++) {
        while (parts[k].x[j] < 0) parts[k].x[j] += dim[j];
        while (parts[k].x[j] >= dim[j]) parts[k].x[j] -= dim[j];
      }
  } else {
    for (int k = 0; k < N; k++)
      for (int j = 0; j < 3; j++)
        if (parts[k].x[j] < 0 || parts[k].x[j] >= dim[j])
          error("Not all particles are within the specified domain.");
  }

  /* Allocate the xtra parts array. */
  if (posix_memalign((void *)&s->xparts, part_align,
                     N * sizeof(struct xpart)) != 0)
    error("Failed to allocate xparts.");
  bzero(s->xparts, N * sizeof(struct xpart));

  /* Initialize the velocities and internal energies. */
  for (int k = 0; k < N; k++) {
    struct part *p = &parts[k];
    struct xpart *xp = &s->xparts[k];
    xp->v_hdt[0] = p->v[0];
    xp->v_hdt[1] = p->v[1];
    xp->v_hdt[2] = p->v[2];
    xp->u_hdt = p->u;
  }

  /* For now, clone the parts to make gparts. */
  if (posix_memalign((void *)&s->gparts, part_align,
                     N * sizeof(struct gpart)) != 0)
    error("Failed to allocate gparts.");
  bzero(s->gparts, N * sizeof(struct gpart));
  /* for ( int k = 0 ; k < N ; k++ ) {
      s->gparts[k].x[0] = s->parts[k].x[0];
      s->gparts[k].x[1] = s->parts[k].x[1];
      s->gparts[k].x[2] = s->parts[k].x[2];
      s->gparts[k].v[0] = s->parts[k].v[0];
      s->gparts[k].v[1] = s->parts[k].v[1];
      s->gparts[k].v[2] = s->parts[k].v[2];
      s->gparts[k].mass = s->parts[k].mass;
      s->gparts[k].dt = s->parts[k].dt;
      s->gparts[k].id = s->parts[k].id;
      s->gparts[k].part = &s->parts[k];
      s->parts[k].gpart = &s->gparts[k];
      }
  s->nr_gparts = s->nr_parts; */
  s->nr_gparts = 0;
  s->size_gparts = s->size_parts;

  /* Init the space lock. */
  if (lock_init(&s->lock) != 0) error("Failed to create space spin-lock.");

  /* Build the cells and the tasks. */
  space_regrid(s, h_max, verbose);
}
