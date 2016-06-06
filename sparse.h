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
  void printA();
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
  int it_col(){ return JA[iterA]; }
  int it_row();
  void next(){ ++iterA; }
  bool end(){ return iterA == nnz; }


  void insert(double v, int g_col, int g_row);

  size_t csr_size();
  static Sparse* mpi_create(int row_no_max, int nnz_max){
    return create(
        row_no_max,
        nnz_max,
        -1,
        -1,
        row_no_max,
        -1,
        nnz_max);
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

  int split_nnz_max(bool by_col, int block_count);
  int split_row_no_max(bool by_col, int block_count);
  Sparse** split(bool by_col, int block_count, int split_row_no_max, int split_nnz_max);


  //------------------------
  //---  MPI
  //------------------------
  MPI_Request send_req, recv_req;
  void send(int rank);
  void recv(int rank, int block_no);
  static bool MPI_Test(MPI_Request* req);
  bool send_ready();
  bool recv_ready();
  void recv_wait(){
    if(!recv_ready()){
      MPI_Wait(&recv_req);
      recv_ready();
    }
  }
  void send_wait(){
    if(!send_ready()){
      MPI_Wait(&send_req);
      send_ready();
    }
  }
  //------------------------

  void update_refs();

  int block_no;

  private:
  int iterA, iterIA;
  protected:
  int nnz_max, row_no_max;

};

#endif
