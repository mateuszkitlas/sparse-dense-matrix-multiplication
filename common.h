#ifndef COMMON_H
#define COMMON_H

#include <stdio.h>
#include <cstdio>
#include <mpi.h>

#ifdef DEBUG

#define debug_s(arg) { fprintf(stderr, "%s:%d in %s()  | %s=%s \n", __FILE__, __LINE__, __FUNCTION__, #arg, arg); }
#define debug_d(arg) { fprintf(stderr, "%s:%d in %s()  | %s=%d \n", __FILE__, __LINE__, __FUNCTION__, #arg, arg); }

#else

#define debug_s(arg) {;};
#define debug_d(arg) {;};

#endif

int MPI_Wait(MPI_Request* request);

extern int mpi_meta_init[2];
extern int mpi_meta_init_size;
extern int &split_row_no_max;
extern int &split_nnz_max;

#endif
