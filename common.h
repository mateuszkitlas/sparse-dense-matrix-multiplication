#ifndef COMMON_H
#define COMMON_H

#include <stdio.h>
#include <cstdio>
#include <mpi.h>
#include <assert.h>

#ifdef DEBUG

#define debug_s(arg) { fprintf(stderr, "%d  %s:%d in %s()  | %s=%s \n", mpi_rank, __FILE__, __LINE__, __FUNCTION__, #arg, arg); }
#define debug_d(arg) { fprintf(stderr, "%d  %s:%d in %s()  | %s=%d \n", mpi_rank, __FILE__, __LINE__, __FUNCTION__, #arg, arg); }
#define debug_f(arg) { fprintf(stderr, "%d  %s:%d in %s()  | %s=%0.2lf \n", mpi_rank, __FILE__, __LINE__, __FUNCTION__, #arg, arg); }
#define debug_p(a,b,c) { fprintf(stderr, "%d  %s:%d in %s()  | [%d,%d]%0.2lf \n", mpi_rank, __FILE__, __LINE__, __FUNCTION__, a,b,c); }
#define debug(arg) { fprintf(stderr, "%d  %s:%d in %s()  | --- %s --- \n", mpi_rank, __FILE__, __LINE__, __FUNCTION__, arg); }

#else

#define debug_s(arg) ;
#define debug_d(arg) ;
#define debug(arg) ;
#define NDEBUG

#endif

#define ASSERTS \
  assert(block_count > 0);


struct Mpi_meta_init{
  int row_no_max;
  int nnz_max;
  int side;
};
extern Mpi_meta_init mpi_meta_init;

int MPI_Wait(MPI_Request* request);

extern int mpi_rank;
extern int num_processes;
extern void init_block_count(bool by_col);
extern int block_count;
extern int block_size(int block_no);
extern int first_side(bool get_col, bool by_col, int block_no);
extern int repl_fact;
extern int which_block(int matrix_i);

int mpi_no(int);

#endif
