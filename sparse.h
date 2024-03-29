#ifndef SPARSE_H
#define SPARSE_H

#include <stdio.h>
#include <cstdio>
#include <mpi.h>
#include "common.h"

#define SPFOR(sp) for(sp->begin(); !sp->end(); sp->next())

class Sparse {
  public:
  void print();
  void printA(char*);
  //void printJA();
  //void printIA();

  void* csr;
  int &first_row, &first_col, &row_no, &col_no, &nnz;
  int *IA, *JA;
  double *A;

  static size_t csr_alloc_size(
      int row_no_max,
      int nnz_max);
  static void* csr_alloc(
      int row_no_max,
      int nnz_max);

  int in_row;
  void begin();
  double it_val(){ return A[iterA]; }
  double val(){ return it_val(); }
  int it_col(){ return JA[iterA]; }
  int it_row();
  int row(){
    return it_row() + first_row;
  }
  int col(){
    int result = it_col() + first_col;
    if(result < 0){
      debug_d(0+result);
      debug_d(0+first_col);
      debug_d(0+it_col());
    }
    return it_col() + first_col;
  }
  void next(){ ++iterA; }
  bool end(){ return iterA == nnz; }


  void insert(double v, int g_col, int g_row);
  void done_insert();

  size_t csr_size();
  static Sparse* mpi_create(){
    return create(
        ::mpi_meta_init.row_no_max,
        ::mpi_meta_init.nnz_max,
        -1,
        -1,
        ::mpi_meta_init.row_no_max,
        -1,
        ::mpi_meta_init.nnz_max);
  }
  static Sparse* create(
      int row_no_max,
      int nnz_max,
      int first_row,
      int first_col,
      int row_no,
      int col_no,
      int nnz);
  Sparse(void* csr, int row_no_max, int nnz_max);
  void free_csr();
  int side(); //row_no; asserts row_no == col_no
  void test();

  //------------------------
  //---  MPI
  //------------------------
  int send_counter;
  int recv_counter;
  bool done_multiplication;
  MPI_Request send_req, recv_req;
  void _send(int rank);
  void _recv(int rank, int block_no);
  void send();
  void recv();
  static bool MPI_Test(MPI_Request* req);
  bool send_ready();
  bool recv_ready();
  void recv_wait();
  void send_wait();
  void post_async_recv(int rank, int block_no);
  //------------------------

  void update_refs();

  int block_no;

  int iterA, iterIA;
  protected:
  int nnz_max, row_no_max;

};

#endif
