
# This file is part of SWIFT.
# Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
#                    Matthieu Schaller (matthieu.schaller@durham.ac.uk).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Common flags
MYFLAGS = -DTIMER

# Add the source directory and debug to CFLAGS
AM_CFLAGS = -I../src -DCPU_TPS=2.67e9 $(HDF5_CPPFLAGS)

AM_LDFLAGS =

MPI_THREAD_LIBS = @MPI_THREAD_LIBS@
MPI_LIBS = $(METIS_LIBS) $(MPI_THREAD_LIBS)
MPI_FLAGS = -DWITH_MPI $(METIS_INCS)

# Set-up the library
bin_PROGRAMS = swift swift_fixdt swift_mindt

# Build MPI versions as well?
if HAVEMPI
bin_PROGRAMS += swift_mpi swift_fixdt_mpi swift_mindt_mpi
endif

# engine_policy_setaffinity is available?
if HAVESETAFFINITY
ENGINE_POLICY_SETAFFINITY=| engine_policy_setaffinity
else
ENGINE_POLICY_SETAFFINITY=
endif

# Sources for swift
swift_SOURCES = main.c
swift_CFLAGS = $(MYFLAGS) $(AM_CFLAGS) -DENGINE_POLICY="engine_policy_multistep | engine_policy_keep $(ENGINE_POLICY_SETAFFINITY)"
swift_LDADD =  ../src/.libs/libswiftsim.a $(HDF5_LDFLAGS) $(HDF5_LIBS)

# Sources for swift_fixdt
swift_fixdt_SOURCES = main.c
swift_fixdt_CFLAGS = $(MYFLAGS) $(AM_CFLAGS) -DENGINE_POLICY="engine_policy_fixdt | engine_policy_keep $(ENGINE_POLICY_SETAFFINITY)"
swift_fixdt_LDADD =  ../src/.libs/libswiftsim.a $(HDF5_LDFLAGS) $(HDF5_LIBS)

# Sources for swift_mindt
swift_mindt_SOURCES = main.c
swift_mindt_CFLAGS = $(MYFLAGS) $(AM_CFLAGS) -DENGINE_POLICY="engine_policy_keep $(ENGINE_POLICY_SETAFFINITY)"
swift_mindt_LDADD =  ../src/.libs/libswiftsim.a $(HDF5_LDFLAGS) $(HDF5_LIBS)

# Sources for swift_mpi
swift_mpi_SOURCES = main.c
swift_mpi_CFLAGS = $(MYFLAGS) $(AM_CFLAGS) $(MPI_FLAGS) -DENGINE_POLICY="engine_policy_multistep | engine_policy_keep"
swift_mpi_LDADD =  ../src/.libs/libswiftsim_mpi.a $(HDF5_LDFLAGS) $(HDF5_LIBS) $(MPI_LIBS)

# Sources for swift_fixdt_mpi
swift_fixdt_mpi_SOURCES = main.c
swift_fixdt_mpi_CFLAGS = $(MYFLAGS) $(AM_CFLAGS) $(MPI_FLAGS) -DENGINE_POLICY="engine_policy_fixdt | engine_policy_keep"
swift_fixdt_mpi_LDADD =  ../src/.libs/libswiftsim_mpi.a $(HDF5_LDFLAGS) $(HDF5_LIBS) $(MPI_LIBS)

# Sources for swift_mindt_mpi
swift_mindt_mpi_SOURCES = main.c
swift_mindt_mpi_CFLAGS = $(MYFLAGS) $(AM_CFLAGS) $(MPI_FLAGS) -DENGINE_POLICY="engine_policy_keep"
swift_mindt_mpi_LDADD =  ../src/.libs/libswiftsim_mpi.a $(HDF5_LDFLAGS) $(HDF5_LIBS) $(MPI_LIBS)

