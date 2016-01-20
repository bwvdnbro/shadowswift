/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Matthieu Schaller (matthieu.schaller@durham.ac.uk).
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
#ifndef SWIFT_DEBUG_H
#define SWIFT_DEBUG_H

/* Includes. */
#include "cell.h"
#include "part.h"

void printParticle(struct part *parts, long long int i, int N);
void printgParticle(struct gpart *parts, long long int i, int N);
void printParticle_single(struct part *p);

#ifdef HAVE_METIS
#include "metis.h"
void dumpMETISGraph(const char *prefix, idx_t nvtxs, idx_t ncon,
                    idx_t *xadj, idx_t *adjncy, idx_t *vwgt, idx_t *vsize,
                    idx_t *adjwgt);

#endif
#endif /* SWIFT_DEBUG_H */
