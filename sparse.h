#ifndef SPARSE_H
#define SPARSE_H

#include <stdio.h>
#include <cstdio>
#include <mpi.h>

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

  static size_t csr_alloc_size(int row_no, int nnz);
  static void* csr_alloc();
  static void* csr_alloc(int row_no_max, int nnz_max);

  int in_row;
  void begin();
  double it_val(){ return A[iterA]; }
  double val(){ return it_val(); }
  int it_col(){ return JA[iterA]; }
  int it_row();
  int row(){
    return it_row() + first_row;
  }
  int col();
  void next(){ ++iterA; }
  bool end(){ return iterA == nnz; }


  void insert(double v, int g_col, int g_row);
  void done_insert();

  size_t csr_size();
  static Sparse* mpi_create();
  static Sparse* create(
      int first_row,
      int first_col,
      int row_no,
      int col_no,
      int nnz);
  Sparse(void* csr);
  void free_csr();
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
};

#endif
