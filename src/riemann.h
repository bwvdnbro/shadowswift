/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2015 Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
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
#ifndef RIEMANN_H
#define RIEMANN_H

#define EXACT_SOLVER

#ifdef EXACT_SOLVER
#include "riemann_exact.h"
#endif

#ifdef TRRS_SOLVER
#include "riemann_trrs.h"
#endif

#endif //RIEMANN_H
