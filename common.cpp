#include "common.h"

Mpi_meta_init mpi_meta_init = {-1, -1, -1};

int mpi_rank =  -1;
int num_processes = -1;
void init_block_count(bool by_col){
  block_count = by_col ? num_processes : num_processes / repl_fact;
};
int block_count = -1;
int repl_fact = 1;

int mpi_no(int x){
  return (num_processes + mpi_rank + x) % num_processes;
}

int MPI_Wait(MPI_Request* request){
  return MPI_Wait(request, MPI_STATUS_IGNORE);
}

//divides the most equally
int block_size(int block_no){
  ASSERTS;
  int side = mpi_meta_init.side;
  return side / block_count + ((side % block_count) >= (block_no + 1) );
}

//TODO remove first_row and first_col from Sparse
int first_side(bool get_col, bool by_col, int block_no){
  ASSERTS;
  int &side = mpi_meta_init.side;

  if(by_col) //column A
    if(!get_col) //first_row
      return 0;

  if(!by_col) //inner ABC
    if(get_col) //first_col
      return 0;

  int bigger_block_size = block_size(0);
  int smaller_block_size = bigger_block_size - 1;
  int smaller_blocks_count = block_no - side % block_count;
  if(smaller_blocks_count > 0)
    return smaller_blocks_count * smaller_block_size + bigger_block_size * (side % block_count);
  else
    return block_no * bigger_block_size;
}
