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
#ifndef SWIFT_SWIFT_H
#define SWIFT_SWIFT_H

/* Config parameters. */
#include "../config.h"

/* Local headers. */
#include "atomic.h"
#include "cell.h"
#include "const.h"
#include "const.h"
#include "cycle.h"
#include "debug.h"
#include "engine.h"
#include "error.h"
#include "lock.h"
#include "map.h"
#include "multipole.h"
#include "parallel_io.h"
#include "part.h"
#include "queue.h"
#include "runner.h"
#include "scheduler.h"
#include "serial_io.h"
#include "single_io.h"
#include "space.h"
#include "task.h"
#include "timers.h"
#include "units.h"
#include "tools.h"
#include "version.h"

#ifdef LEGACY_GADGET2_SPH
#include "runner_iact_legacy.h"
#else
#include "runner_iact.h"
#endif
#include "runner_iact_grav.h"

#endif /* SWIFT_SWIFT_H */
