#ifndef COMMON_H
#define COMMON_H

#include <stdio.h>
#include <cstdio>
#include <mpi.h>
#include <assert.h>
#include "sparse.h"

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


#ifdef IDENTITY_MATRIX
#undef NDEBUG
#endif


struct Mpi_meta_init{
  int row_no_max;
  int nnz_max;
  int side;
};
extern void broadcast_metadata();
extern void compute_metadata();

extern int use_inner;

int MPI_Wait(MPI_Request* request);

extern int mpi_rank;
extern int num_processes;
extern int block_count;
extern int bigger_blocks_count;
extern int& row_no_max;
extern int& nnz_max;
extern int& side;
extern int block_size(int block_no);
extern int first_side(int block_no);
extern int repl_fact;
extern int which_block(int matrix_i);
extern int max_block_size;
extern int min_block_size;

int mpi_no(int);


extern int my_block_col_no;
//----------------------
//---- inner
//----------------------

extern MPI_Comm mpi_inner_group_comm;
extern int my_block_row_no;
extern void inner_replicate_sparse(Sparse* my_sparse);

//----------------------
//---- column
//----------------------


//----------------------

#endif
