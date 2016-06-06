#include "common.h"

int mpi_meta_init[2] = {-1,-1};
int mpi_meta_init_size = 2;
int &split_row_no_max = mpi_meta_init[0];
int &split_nnz_max = mpi_meta_init[1];
int mpi_rank = 0;
int num_processes = 1;
int &block_count = num_processes;

int mpi_no(int x){
  return (num_processes + mpi_rank + x) % num_processes;
}

int MPI_Wait(MPI_Request* request){
  return MPI_Wait(request, MPI_STATUS_IGNORE);
}
