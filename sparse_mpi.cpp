#include "sparse.h"
#include <stdlib.h>
#include <assert.h>


void Sparse::send(int rank){
  assert(block_no >= 0);
  MPI_Isend(csr, csr_size(), MPI_BYTE, rank, 1000+block_no, MPI_COMM_WORLD, &send_req);
}

void Sparse::recv(int rank, int block_no){
  assert(block_no >= 0);
  MPI_Irecv(
      csr,
      csr_alloc_size(row_no_max, nnz_max),
      MPI_BYTE,
      rank, 1000+block_no, MPI_COMM_WORLD, &recv_req);
  this->block_no = block_no;
}

bool Sparse::MPI_Test(MPI_Request* req){
  if(*req == MPI_REQUEST_NULL){
    debug_s("warning");
    return true;
  } else {
    int flag;
    ::MPI_Test(req, &flag, MPI_STATUS_IGNORE);
    return 1==flag;
  }
}

bool Sparse::send_ready(){
  return Sparse::MPI_Test(&send_req);
}

bool Sparse::recv_ready(){
  if(Sparse::MPI_Test(&recv_req)){
    update_refs();
    return true;
  } else {
    return false;
  }
}
