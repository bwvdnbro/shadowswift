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
#ifndef SWIFT_RUNNER_H
#define SWIFT_RUNNER_H

/* Includes. */
#include "cell.h"
#include "inline.h"

extern const float runner_shift[13 * 3];
extern const char runner_flip[27];

/* A struct representing a runner's thread and its data. */
struct runner {

  /* The id of this thread. */
  int id;

  /* The thread which it is running. */
  pthread_t thread;

  /* The queue to use to get tasks. */
  int cpuid, qid;

  /* The underlying runner. */
  struct engine *e;
};

/* Function prototypes. */
void runner_doghost(struct runner *r, struct cell *c);
void runner_dopair_density(struct runner *r, struct cell *ci, struct cell *cj);
void runner_doself_density(struct runner *r, struct cell *c);
void runner_dosub_density(struct runner *r, struct cell *ci, struct cell *cj,
                          int flags);
void runner_dosort(struct runner *r, struct cell *c, int flag, int clock);
void runner_dogsort(struct runner *r, struct cell *c, int flag, int clock);
void runner_dokick(struct runner *r, struct cell *c, int timer);
void runner_dokick1(struct runner *r, struct cell *c);
void runner_dokick2(struct runner *r, struct cell *c);
void *runner_main(void *data);

#endif /* SWIFT_RUNNER_H */
