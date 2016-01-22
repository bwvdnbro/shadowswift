/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2015 Peter W. Draper (p.w.draper@durham.ac.uk)
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

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <pthread.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <fenv.h>

/* Conditional headers. */
#ifdef HAVE_LIBZ
#include <zlib.h>
#endif

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* Local headers. */
#include "swift.h"
#include "voronoi.h"

/* Ticks per second on this machine. */
#ifndef CPU_TPS
#define CPU_TPS 2.40e9
#endif

/* Engine policy flags. */
#ifndef ENGINE_POLICY
#define ENGINE_POLICY engine_policy_none
#endif




/**
 * @brief Main routine that loads a few particles and generates some output.
 *
 */

int main(int argc, char *argv[]) {

  int c, icount, j, k, N, periodic = 1;
  long long N_total = -1;
  int nr_threads = 1, nr_queues = -1, runs = INT_MAX;
  int dump_tasks = 0;
  int data[2];
  double dim[3] = {1.0, 1.0, 1.0}, shift[3] = {0.0, 0.0, 0.0};
  double h_max = -1.0, scaling = 1.0;
  double clock = DBL_MAX;
  struct part *parts = NULL;
  struct space s;
  struct engine e;
  struct UnitSystem us;
  char ICfileName[200] = "";
  char dumpfile[30];
  float dt_max = 0.0f;
  float dt_snap = 0.0f;
  unsigned int numsnap = 1;
  ticks tic;
  int nr_nodes = 1, myrank = 0, grid[3] = {1, 1, 1};
  FILE *file_thread;
  int with_outputs = 1;

/* Choke on FP-exceptions. */
// feenableexcept( FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW );

#ifdef WITH_MPI
  /* Start by initializing MPI. */
  int res, prov;
  if ((res = MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &prov)) != MPI_SUCCESS)
    error("Call to MPI_Init failed with error %i.", res);
  if (prov != MPI_THREAD_MULTIPLE)
    error("MPI does not provide the level of threading required "
          "(MPI_THREAD_MULTIPLE).");
  if ((res = MPI_Comm_size(MPI_COMM_WORLD, &nr_nodes)) != MPI_SUCCESS)
    error("MPI_Comm_size failed with error %i.", res);
  if ((res = MPI_Comm_rank(MPI_COMM_WORLD, &myrank)) != MPI_SUCCESS)
    error("Call to MPI_Comm_rank failed with error %i.", res);
  if ((res = MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN)) !=
      MPI_SUCCESS)
    error("Call to MPI_Comm_set_errhandler failed with error %i.", res);
  if (myrank == 0) message("MPI is up and running with %i node(s).", nr_nodes);
  fflush(stdout);

  /* Set a default grid so that grid[0]*grid[1]*grid[2] == nr_nodes. */
  factor(nr_nodes, &grid[0], &grid[1]);
  factor(nr_nodes / grid[1], &grid[0], &grid[2]);
  factor(grid[0] * grid[1], &grid[1], &grid[0]);
#endif

  /* Greeting message */
  if (myrank == 0) greetings();

  /* Init the space. */
  bzero(&s, sizeof(struct space));

  /* Parse the options */
  while ((c = getopt(argc, argv, "a:c:d:f:g:i:m:q:r:s:t:w:y:z:")) != -1)
    switch (c) {
      case 'a':
        if (sscanf(optarg, "%lf", &scaling) != 1)
          error("Error parsing cutoff scaling.");
        if (myrank == 0) message("scaling cutoff by %.3f.", scaling);
        fflush(stdout);
        break;
      case 'c':
        if (sscanf(optarg, "%lf", &clock) != 1) error("Error parsing clock.");
        if (myrank == 0) message("clock set to %.3e.", clock);
        fflush(stdout);
        break;
      case 'd':
        if (sscanf(optarg, "%f", &dt_max) != 1)
          error("Error parsing timestep.");
        if (myrank == 0) message("dt set to %e.", dt_max);
        fflush(stdout);
        break;
      case 'f':
        if (!strcpy(ICfileName, optarg)) error("Error parsing IC file name.");
        break;
      case 'g':
        if (sscanf(optarg, "%i %i %i", &grid[0], &grid[1], &grid[2]) != 3)
          error("Error parsing grid.");
        break;
      case 'i':
        if (sscanf(optarg, "%f", &dt_snap) != 1)
          error("Error parsing snapshot interval.");
        break;
      case 'm':
        if (sscanf(optarg, "%lf", &h_max) != 1) error("Error parsing h_max.");
        if (myrank == 0) message("maximum h set to %e.", h_max);
        fflush(stdout);
        break;
      case 'o':
        with_outputs = 0;
        break;
      case 'q':
        if (sscanf(optarg, "%d", &nr_queues) != 1)
          error("Error parsing number of queues.");
        break;
      case 'r':
        if (sscanf(optarg, "%d", &runs) != 1)
          error("Error parsing number of runs.");
        break;
      case 's':
        if (sscanf(optarg, "%lf %lf %lf", &shift[0], &shift[1], &shift[2]) != 3)
          error("Error parsing shift.");
        if (myrank == 0)
          message("will shift parts by [ %.3f %.3f %.3f ].", shift[0], shift[1],
                  shift[2]);
        break;
      case 't':
        if (sscanf(optarg, "%d", &nr_threads) != 1)
          error("Error parsing number of threads.");
        break;
      case 'w':
        if (sscanf(optarg, "%d", &space_subsize) != 1)
          error("Error parsing sub size.");
        if (myrank == 0) message("sub size set to %i.", space_subsize);
        break;
      case 'y':
        if(sscanf(optarg, "%d", &dump_tasks) != 1)
          error("Error parsing dump_tasks (-y)");
        break;
      case 'z':
        if (sscanf(optarg, "%d", &space_splitsize) != 1)
          error("Error parsing split size.");
        if (myrank == 0) message("split size set to %i.", space_splitsize);
        break;
      case '?':
        error("Unknown option.");
        break;
    }

#if defined(WITH_MPI)
  if (myrank == 0) {
    message("Running with %i thread(s) per node.", nr_threads);
    message("grid set to [ %i %i %i ].", grid[0], grid[1], grid[2]);

    if (nr_nodes == 1) {
      message("WARNING: you are running with one MPI rank.");
      message("WARNING: you should use the non-MPI version of this program." );
    }
    fflush(stdout);
  }
#else
  if (myrank == 0) message("Running with %i thread(s).", nr_threads);
#endif

  /* How large are the parts? */
  if (myrank == 0) {
    message("sizeof(struct part) is %li bytes.", (long int)sizeof(struct part));
    message("sizeof(struct gpart) is %li bytes.",
            (long int)sizeof(struct gpart));
  }

  /* Initialize unit system */
  initUnitSystem(&us);
  if (myrank == 0) {
    message("Unit system: U_M = %e g.", us.UnitMass_in_cgs);
    message("Unit system: U_L = %e cm.", us.UnitLength_in_cgs);
    message("Unit system: U_t = %e s.", us.UnitTime_in_cgs);
    message("Unit system: U_I = %e A.", us.UnitCurrent_in_cgs);
    message("Unit system: U_T = %e K.", us.UnitTemperature_in_cgs);
    message("Density units: %e a^%f h^%f.",
            conversionFactor(&us, UNIT_CONV_DENSITY),
            aFactor(&us, UNIT_CONV_DENSITY), hFactor(&us, UNIT_CONV_DENSITY));
    message("Entropy units: %e a^%f h^%f.",
            conversionFactor(&us, UNIT_CONV_ENTROPY),
            aFactor(&us, UNIT_CONV_ENTROPY), hFactor(&us, UNIT_CONV_ENTROPY));
  }

  /* Check whether an IC file has been provided */
  if (strcmp(ICfileName, "") == 0)
    error("An IC file name must be provided via the option -f");

  /* Read particles and space information from (GADGET) IC */
  tic = getticks();
#if defined(WITH_MPI)
#if defined(HAVE_PARALLEL_HDF5)
  read_ic_parallel(ICfileName, dim, &parts, &N, &periodic, myrank, nr_nodes,
                   MPI_COMM_WORLD, MPI_INFO_NULL);
#else
  read_ic_serial(ICfileName, dim, &parts, &N, &periodic, myrank, nr_nodes,
                 MPI_COMM_WORLD, MPI_INFO_NULL);
#endif
#else
  read_ic_single(ICfileName, dim, &parts, &N, &periodic);
#endif

  if (myrank == 0)
    message("reading particle properties took %.3f ms.",
            ((double)(getticks() - tic)) / CPU_TPS * 1000);
  fflush(stdout);

#if defined(WITH_MPI)
  long long N_long = N;
  MPI_Reduce(&N_long, &N_total, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
#else
  N_total = N;
#endif
  if (myrank == 0) message("Read %lld particles from the ICs", N_total);

  /* Apply h scaling */
  if (scaling != 1.0)
    for (k = 0; k < N; k++) parts[k].h *= scaling;

  /* Apply shift */
  if (shift[0] != 0 || shift[1] != 0 || shift[2] != 0)
    for (k = 0; k < N; k++) {
      parts[k].x[0] += shift[0];
      parts[k].x[1] += shift[1];
      parts[k].x[2] += shift[2];
    }

  /* Moving mesh initialization */
  for (k = 0; k < N; k++) {
    /* Initialize Voronoi cells */
    voronoi_initialize(&parts[k], 2.0f*parts[k].h);

    /* Convert thermal energy to pressure */
    parts[k].primitives.P *= (const_hydro_gamma - 1.0f) *
                               parts[k].primitives.rho;
    /* Initialize flow velocity */
    parts[k].primitives.v[0] = parts[k].v[0];
    parts[k].primitives.v[1] = parts[k].v[1];
    parts[k].primitives.v[2] = parts[k].v[2];
    /* Initialize mass to 0 */
    parts[k].conserved.m = 0.0f;
  }

  /* Set default number of queues. */
  if (nr_queues < 0) nr_queues = nr_threads;

  /* Initialize the space with this data. */
  tic = getticks();
  space_init(&s, dim, parts, N, periodic, h_max, myrank == 0);
  if (myrank == 0)
    message("space_init took %.3f ms.",
            ((double)(getticks() - tic)) / CPU_TPS * 1000);
  fflush(stdout);

  /* Set the default time step to 1.0f. */
  if (myrank == 0) message("dt_max is %e.", dt_max);

  /* Say a few nice things about the space we just created. */
  if (myrank == 0) {
    message("space dimensions are [ %.3f %.3f %.3f ].", s.dim[0], s.dim[1],
            s.dim[2]);
    message("space %s periodic.", s.periodic ? "is" : "isn't");
    message("highest-level cell dimensions are [ %i %i %i ].", s.cdim[0],
            s.cdim[1], s.cdim[2]);
    message("%i parts in %i cells.", s.nr_parts, s.tot_cells);
    message("maximum depth is %d.", s.maxdepth);
    // message( "cutoffs in [ %g %g ]." , s.h_min , s.h_max ); fflush(stdout);
  }

  /* Verify that each particle is in it's proper cell. */
  if (myrank == 0) {
    icount = 0;
    space_map_cells_pre(&s, 0, &map_cellcheck, &icount);
    message("map_cellcheck picked up %i parts.", icount);
  }

  if (myrank == 0) {
    data[0] = s.maxdepth;
    data[1] = 0;
    space_map_cells_pre(&s, 0, &map_maxdepth, data);
    message("nr of cells at depth %i is %i.", data[0], data[1]);
  }

  /* Dump the particle positions. */
  // space_map_parts( &s , &map_dump , shift );

  /* Initialize the engine with this space. */
  tic = getticks();
  if (myrank == 0) message("nr_nodes is %i.", nr_nodes);
  engine_init(&e, &s, dt_max, nr_threads, nr_queues, nr_nodes, myrank,
              ENGINE_POLICY | engine_policy_steal);
  if (myrank == 0)
    message("engine_init took %.3f ms.",
            ((double)(getticks() - tic)) / CPU_TPS * 1000);
  fflush(stdout);

#ifdef WITH_MPI
  /* Split the space. */
  engine_split(&e, grid);
  engine_redistribute(&e);
#endif

  if (with_outputs) {
    /* Write the state of the system as it is before starting time integration.
     */
    tic = getticks();
#if defined(WITH_MPI)
#if defined(HAVE_PARALLEL_HDF5)
    write_output_parallel(&e, &us, myrank, nr_nodes, MPI_COMM_WORLD,
                          MPI_INFO_NULL);
#else
    write_output_serial(&e, &us, myrank, nr_nodes, MPI_COMM_WORLD,
                        MPI_INFO_NULL);
#endif
#else
    write_output_single(&e, &us);
#endif
    if (myrank == 0)
      message("writing particle properties took %.3f ms.",
              ((double)(getticks() - tic)) / CPU_TPS * 1000);
    fflush(stdout);
  }

/* Init the runner history. */
#ifdef HIST
  for (k = 0; k < runner_hist_N; k++) runner_hist_bins[k] = 0;
#endif

  if (myrank == 0) {
    /* Inauguration speech. */
    if (runs < INT_MAX) {
      message(
          "Running on %lld particles for %i steps with %i threads and %i "
          "queues...",
          N_total, runs, e.nr_threads, e.sched.nr_queues);
    } else {
      message(
          "Running on %lld particles until t=%.3e with %i threads and %i "
          "queues...",
          N_total, clock, e.nr_threads, e.sched.nr_queues);
      if(!dt_snap){
        dt_snap = 0.1*clock;
      }
    }
    message("Writing a snapshot approximately every Delta_t=%.3e", dt_snap);
    fflush(stdout);
  }

  /* Legend */
  if (myrank == 0)
    printf("# Step  Time  time-step  CPU Wall-clock time [ms]\n");

  /* Let loose a runner on the space. */
  for (j = 0; j < runs && e.time < clock; j++) {

/* Repartition the space amongst the nodes? */
#if defined(WITH_MPI) && defined(HAVE_METIS)
    if (j % 100 == 2) e.forcerepart = 1;
#endif

    timers_reset(timers_mask_all);
#ifdef COUNTER
    for (k = 0; k < runner_counter_count; k++) runner_counter[k] = 0;
#endif

    /* Take a step. */
    engine_step(&e);

    if (with_outputs && numsnap*dt_snap <= e.time) {

#if defined(WITH_MPI)
#if defined(HAVE_PARALLEL_HDF5)
      write_output_parallel(&e, &us, myrank, nr_nodes, MPI_COMM_WORLD,
                            MPI_INFO_NULL);
#else
      write_output_serial(&e, &us, myrank, nr_nodes, MPI_COMM_WORLD,
                          MPI_INFO_NULL);
#endif
#else
      write_output_single(&e, &us);
#endif
      numsnap++;
    }

  /* Dump the task data using the given frequency. */
    if (dump_tasks && (dump_tasks == 1 || j % dump_tasks == 1)) {
#ifdef WITH_MPI

      /* Make sure output file is empty, only on one rank. */
      sprintf(dumpfile, "thread_info_MPI-step%d.dat", j);
      if (myrank == 0) {
          file_thread = fopen(dumpfile, "w");
          fclose(file_thread);
      }
      MPI_Barrier(MPI_COMM_WORLD);

      for (int i = 0; i < nr_nodes; i++) {

        /* Rank 0 decides the index of writing node, this happens one-by-one. */
        int kk = i;
        MPI_Bcast(&kk, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (i == myrank) {

          /* Open file and position at end. */
          file_thread = fopen(dumpfile, "a");

          fprintf(file_thread, " %03i 0 0 0 0 %lli 0 0 0 0\n", myrank,
                  e.tic_step);
          int count = 0;
          for (int l = 0; l < e.sched.nr_tasks; l++)
            if (!e.sched.tasks[l].skip && !e.sched.tasks[l].implicit) {
              fprintf(
                  file_thread, " %03i %i %i %i %i %lli %lli %i %i %i\n", myrank,
                  e.sched.tasks[l].rid, e.sched.tasks[l].type,
                  e.sched.tasks[l].subtype, (e.sched.tasks[l].cj == NULL),
                  e.sched.tasks[l].tic, e.sched.tasks[l].toc,
                  (e.sched.tasks[l].ci != NULL) ? e.sched.tasks[l].ci->count : 0,
                  (e.sched.tasks[l].cj != NULL) ? e.sched.tasks[l].cj->count : 0,
                  e.sched.tasks[l].flags);
              fflush(stdout);
              count++;
            }
          message("rank %d counted %d tasks", myrank, count);

          fclose(file_thread);
        }

        /* And we wait for all to synchronize. */
        MPI_Barrier(MPI_COMM_WORLD);
      }

#else
      sprintf(dumpfile, "thread_info-step%d.dat", j);
      file_thread = fopen(dumpfile, "w");
      for (int l = 0; l < e.sched.nr_tasks; l++)
        if (!e.sched.tasks[l].skip && !e.sched.tasks[l].implicit)
          fprintf(file_thread, " %i %i %i %i %lli %lli %i %i\n",
                  e.sched.tasks[l].rid, e.sched.tasks[l].type,
                  e.sched.tasks[l].subtype, (e.sched.tasks[l].cj == NULL),
                  e.sched.tasks[l].tic, e.sched.tasks[l].toc,
                  (e.sched.tasks[l].ci == NULL) ? 0 : e.sched.tasks[l].ci->count,
                  (e.sched.tasks[l].cj == NULL) ? 0 : e.sched.tasks[l].cj->count);
      fclose(file_thread);
#endif
    }

    /* Dump a line of aggregate output. */
    /*     if (myrank == 0) { */
    /*       printf("%i %e %.16e %.16e %.16e %.3e %.3e %i %.3e %.3e", j, e.time,
     */
    /*              e.ekin + e.epot, e.ekin, e.epot, e.dt, e.dt_step,
     * e.count_step, */
    /*              e.dt_min, e.dt_max); */
    /*       for (k = 0; k < timer_count; k++) */
    /*         printf(" %.3f", ((double)timers[k]) / CPU_TPS * 1000); */
    /*       printf("\n"); */
    /*       fflush(stdout); */
    /*     } */
    if (myrank == 0) {
      printf("%i %e %.3e", j, e.time, e.dt);
      printf(" %.3f", ((double)timers[timer_count - 1]) / CPU_TPS * 1000);
      printf("\n");
      fflush(stdout);
    }
  }

/* Print the values of the runner histogram. */
#ifdef HIST
  printf("main: runner histogram data:\n");
  for (k = 0; k < runner_hist_N; k++)
    printf(" %e %e %e\n",
           runner_hist_a + k * (runner_hist_b - runner_hist_a) / runner_hist_N,
           runner_hist_a +
               (k + 1) * (runner_hist_b - runner_hist_a) / runner_hist_N,
           (double)runner_hist_bins[k]);
#endif


  if (with_outputs) {
/* Write final output. */
#if defined(WITH_MPI)
#if defined(HAVE_PARALLEL_HDF5)
    write_output_parallel(&e, &us, myrank, nr_nodes, MPI_COMM_WORLD,
                          MPI_INFO_NULL);
#else
    write_output_serial(&e, &us, myrank, nr_nodes, MPI_COMM_WORLD,
                        MPI_INFO_NULL);
#endif
#else
    write_output_single(&e, &us);
#endif
  }

#ifdef WITH_MPI
  if (MPI_Finalize() != MPI_SUCCESS)
    error("call to MPI_Finalize failed with error %i.", res);
#endif

  /* Say goodbye. */
  if (myrank == 0) message("done.");

  /* All is calm, all is bright. */
  return 0;
}
