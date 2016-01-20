/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2013 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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
#ifndef SWIFT_PROXY_H
#define SWIFT_PROXY_H

/* Includes. */
#include "cell.h"
#include "part.h"

/* Some constants. */
#define proxy_buffgrow 1.5
#define proxy_buffinit 100

/* Proxy tag arithmetic. */
#define proxy_tag_shift 8
#define proxy_tag_count 0
#define proxy_tag_parts 1
#define proxy_tag_xparts 2
#define proxy_tag_cells 3

/* Data structure for the proxy. */
struct proxy {

  /* ID of the node this proxy represents. */
  int mynodeID, nodeID;

  /* Incoming cells. */
  struct cell **cells_in;
  struct pcell *pcells_in;
  int nr_cells_in, size_cells_in, size_pcells_in;

  /* Outgoing cells. */
  struct cell **cells_out;
  struct pcell *pcells_out;
  int nr_cells_out, size_cells_out, size_pcells_out;

  /* The parts and xparts buffers for input and output. */
  struct part *parts_in, *parts_out;
  struct xpart *xparts_in, *xparts_out;
  int size_parts_in, size_parts_out;
  int nr_parts_in, nr_parts_out;

/* MPI request handles. */
#ifdef WITH_MPI
  MPI_Request req_parts_count_out, req_parts_count_in;
  MPI_Request req_parts_out, req_parts_in;
  MPI_Request req_xparts_out, req_xparts_in;
  MPI_Request req_cells_count_out, req_cells_count_in;
  MPI_Request req_cells_out, req_cells_in;
#endif
};

/* Function prototypes. */
void proxy_init(struct proxy *p, int mynodeID, int nodeID);
void proxy_parts_load(struct proxy *p, struct part *parts, struct xpart *xparts,
                      int N);
void proxy_parts_exch1(struct proxy *p);
void proxy_parts_exch2(struct proxy *p);
void proxy_addcell_in(struct proxy *p, struct cell *c);
void proxy_addcell_out(struct proxy *p, struct cell *c);
void proxy_cells_exch1(struct proxy *p);
void proxy_cells_exch2(struct proxy *p);

#endif /* SWIFT_PROXY_H */
