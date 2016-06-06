#ifndef COMMON_H
#define COMMON_H

#include <stdio.h>
#include <cstdio>
#include <mpi.h>

#ifdef DEBUG

#define debug_s(arg) { fprintf(stderr, "%d  %s:%d in %s()  | %s=%s \n", mpi_rank, __FILE__, __LINE__, __FUNCTION__, #arg, arg); }
#define debug_d(arg) { fprintf(stderr, "%d  %s:%d in %s()  | %s=%d \n", mpi_rank, __FILE__, __LINE__, __FUNCTION__, #arg, arg); }
#define debug_f(arg) { fprintf(stderr, "%d  %s:%d in %s()  | %s=%0.2lf \n", mpi_rank, __FILE__, __LINE__, __FUNCTION__, #arg, arg); }
#define debug_p(a,b,c) { fprintf(stderr, "%d  %s:%d in %s()  | [%d,%d]%0.2lf \n", mpi_rank, __FILE__, __LINE__, __FUNCTION__, a,b,c); }
#define debug(arg) { fprintf(stderr, "%d  %s:%d in %s()  | --- %s --- \n", mpi_rank, __FILE__, __LINE__, __FUNCTION__, arg); }

#else

#define debug_s(arg) {;};
#define debug_d(arg) {;};
#define debug(arg) {;};

#endif

int MPI_Wait(MPI_Request* request);

extern int mpi_meta_init[2];
extern int mpi_meta_init_size;
extern int &split_row_no_max;
extern int &split_nnz_max;
extern int mpi_rank;

#endif
