
                                Building SWIFT
                                ==============

SWIFT is built from a clean source repository using the commands:

   ./autogen
   ./configure
   make

and from a distribution tarball using:

   ./configure
   make

The compiler choice is GCC by default, but that can be changed using the "CC"
environment variable. This can be just set, or passed on the ./configure
command line, i.e.:

   bash:
      export CC=icc
      ./configure

   [t]csh:
      setenv CC=icc
      ./configure

or:

   ./configure CC=icc

to use an Intel compiler. The main "programs" can be found in the "examples/"
directory.

SWIFT has been successfully built and tested with the following compilers:

  - GCC 4.8.x  
  - Intel ICC 15.0.x
  - clang 3.4.x 

More recent versions and slightly older ones should also be able to
built the software.

By default an attempt to choose suitable set of optimizing compiler flags
will be made, targeted for the host machine of the build. If this doesn't
work or the binaries will for another architecture then you can stop the
selection of flags using:

   ./configure --disable-optimization

and then supply your own flags using the "CFLAGS" environment variable, as for
CC.

Note that any CFLAGS that you supply will be added to those determined by
configure in all circumstances. To build SWIFT with debugging support you
can use:

    ./configure --enable-debug --disable-optimization

You could also add some additional flags:

    ./configure --enable-debug --disable-optimization CFLAGS="-O2"

for instance. GCC address sanitizer flags can be included using the 

    ./configure --enable-sanitizer

option. Note this requires a GCC compiler version of at least 4.8.


                                 Dependencies
                                 ============

SWIFT depends on a number of third party libraries that should be available
before you can build it.


HDF5: a HDF5 library is required to read and write particle data. One of the
commands "h5cc" or "h5pcc" should be available. If "h5pcc" is located them a
parallel HDF5 built for the version of MPI located should be provided. If
the command is not available then it can be located using the "--with-hfd5"
configure option. The value should be the full path to the "h5cc" or "h5pcc"
commands.


MPI: an optional MPI library that fully supports MPI_THREAD_MULTIPLE.  
Before running configure the "mpirun" command should be available in the
shell. If your command isn't called "mpirun" then define the "MPIRUN"
environment variable, either in the shell or when running configure.


METIS: a build of the METIS library can be optionally used to optimize the
load between MPI nodes (requires an MPI library). This should be found in the
standard installation directories, or pointed at using the "--with-metis"
configuration option.  In this case the top-level installation directory of
the METIS build should be given. Note to use METIS you should at least supply
"--with-metis".


DOXYGEN: the doxygen library is required to create the SWIFT API documentation.
