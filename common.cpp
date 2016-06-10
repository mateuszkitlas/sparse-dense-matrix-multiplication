#include "common.h"

Mpi_meta_init mpi_meta_init = {-1, -1, -1};
int bigger_blocks_count = -1;
int use_inner = 0;
int mpi_rank = -1;
int num_processes = -1;
int block_count = -1;
int min_block_size = -1;
int max_block_size = -1;
int &side = mpi_meta_init.side;
int &row_no_max = mpi_meta_init.row_no_max;
int &nnz_max = mpi_meta_init.nnz_max;
int repl_fact = 1;
bool inner_coordinator = false;



int my_block_col_no = -1;
//----------------------
//---- inner
//----------------------

MPI_Comm mpi_inner_group_comm;
int my_block_row_no = -1;

//----------------------
//---- column
//----------------------


//----------------------



void compute_metadata(){
  block_count = use_inner
    ? num_processes / repl_fact
    : num_processes;
  bigger_blocks_count = side % block_count;
  min_block_size = side / block_count;
  max_block_size = (bigger_blocks_count > 0) + min_block_size;
  if(use_inner){
    my_block_col_no = mpi_rank % block_count;
    my_block_row_no = mpi_rank / block_count;
  } else {
    my_block_col_no = mpi_rank;
  }
}

void broadcast_metadata(){
  MPI_Bcast(&mpi_meta_init, sizeof(mpi_meta_init), MPI_BYTE, 0, MPI_COMM_WORLD);
}

void inner_replicate_sparse(Sparse* my_sparse){
  MPI_Comm_split(MPI_COMM_WORLD, my_block_row_no, my_block_col_no, &mpi_inner_group_comm);
  int csr_size = (my_block_col_no == 0)
    ? my_sparse->csr_size()
    : Sparse::csr_alloc_size(row_no_max, nnz_max);
  MPI_Bcast(my_sparse->csr, csr_size, MPI_BYTE, 0, mpi_inner_group_comm);
}


int mpi_no(int x){
  return (num_processes + mpi_rank + x) % num_processes;
}
int mpi_inner_no(int x){
  int pc2 = num_processes / repl_fact / repl_fact; //p/(c^2)
  int anti_group_row = (x + pc2) % pc2;
  int row = (my_block_row_no + anti_group_row + block_count) % block_count;
  return row * block_count + my_block_col_no;
}

int MPI_Wait(MPI_Request* request){
  return MPI_Wait(request, MPI_STATUS_IGNORE);
}

//divides the most equally
int block_size(int block_no){
  ASSERTS;
  return side / block_count + ((side % block_count) >= (block_no + 1) );
}

//TODO remove first_row and first_col from Sparse
int first_side(int block_no){
  ASSERTS;
  int smaller_blocks_count = block_no - side % block_count;
  if(smaller_blocks_count > 0)
    return smaller_blocks_count * min_block_size + max_block_size * (side % block_count);
  else
    return block_no * max_block_size;
}


int which_block(int matrix_i){
  int i2 = matrix_i - max_block_size * bigger_blocks_count;
#ifdef DEBUG
  int x;
  if(i2 < 0)
     x = matrix_i / max_block_size;
  else
     x = bigger_blocks_count + i2 / min_block_size;
  //fprintf(stderr, 
  //    "side %d, bigger blocks count%d, block_count %d, max_block_size() %d, min_block_size() %d, matrix_i %d, i2=%d, result=%d\n",
  //    side, bigger_blocks_count, block_count, max_block_size(), min_block_size(), matrix_i, i2, x);
#endif
  if(i2 < 0)
    return matrix_i / max_block_size;
  else
    return bigger_blocks_count + i2 / min_block_size;
}

#ifdef DEBUG
void debug_print(const char* file, int line, const char* fun, const char* arg_name, const char* s_arg, int d_arg, char type){
  char printer[1000];
  if(type == 's')
    sprintf(printer, "---------  %s  --------", s_arg);
  else
    sprintf(printer, "%d", d_arg);
  fprintf(stderr, "[r%d,c%d]  %s:%d in %s()  | %s=%s \n",
      my_block_row_no,
      my_block_col_no,
      file,
      line,
      fun,
      arg_name,
      printer
      );
}
#endif
