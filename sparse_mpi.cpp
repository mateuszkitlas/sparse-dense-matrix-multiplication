#include "sparse.h"
#include "common.h"
#include <stdlib.h>


void Sparse::send(){
  _send( mpi_no(-repl_fact) );
}

void Sparse::recv(){
  _recv( mpi_no(repl_fact), (block_no + repl_fact) % block_count );
}

void Sparse::_send(int rank){
#ifdef DEBUG
  //printf("send %d -> %d, block_no=%d\n", mpi_rank, rank, block_no);
#endif
  assert(block_no >= 0);
  MPI_Isend(csr, csr_size(), MPI_BYTE, rank, 1000+block_no, MPI_COMM_WORLD, &send_req);
  assert(send_req != MPI_REQUEST_NULL);
  send_counter++;
}

void Sparse::_recv(int rank, int block_no){
  assert(block_no >= 0);
  assert(Sparse::csr_alloc_size(row_no_max, nnz_max) > 0);
  MPI_Irecv(
      csr,
      Sparse::csr_alloc_size(row_no_max, nnz_max),
      MPI_BYTE,
      rank, 1000+block_no, MPI_COMM_WORLD, &recv_req);
  //assert(recv_req != MPI_REQUEST_NULL);
  post_async_recv(rank, block_no);
}

void Sparse::post_async_recv(int rank, int block_no){
  this->block_no = block_no;
  recv_counter++;
  done_multiplication = false;
}

bool Sparse::MPI_Test(MPI_Request* req){
  if(*req == MPI_REQUEST_NULL){
    //debug_s("warning");
    //debug_d(*req);
    //debug_d(MPI_REQUEST_NULL);
    return true;
  } else {
    int flag;
    ::MPI_Test(req, &flag, MPI_STATUS_IGNORE);
    return *req == MPI_REQUEST_NULL;
  }
}

bool Sparse::send_ready(){
  return Sparse::MPI_Test(&send_req);
}

bool Sparse::recv_ready(){
  if(Sparse::MPI_Test(&recv_req)){
    update_refs();
    //debug_d(block_no);
    if(first_col < 0){
      print();
    }
    assert(first_col >= 0);
    assert(first_row >= 0);
    return true;
  } else {
    return false;
  }
}

void Sparse::recv_wait(){
  if(!recv_ready()){
    MPI_Wait(&recv_req);
    recv_ready();
  }
}
void Sparse::send_wait(){
  if(!send_ready()){
    MPI_Wait(&send_req);
    send_ready();
  }
}
