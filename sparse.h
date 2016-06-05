#ifndef SPARSE_H
#define SPARSE_H

#include <stdio.h>
#include <cstdio>
#include <mpi.h>

#ifdef DEBUG
#define debug_s(arg) { fprintf(stderr, "%s:%d in %s()  | %s=%s \n", __FILE__, __LINE__, __FUNCTION__, #arg, arg); }
#define debug_d(arg) { fprintf(stderr, "%s:%d in %s()  | %s=%d \n", __FILE__, __LINE__, __FUNCTION__, #arg, arg); }
#else
#define debug_s(arg) {;};
#define debug_d(arg) {;};
#endif


class Sparse {
  public:
  void print();

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
  Sparse** split(bool by_col, int block_count);

  void begin(){
    iterA = 0;
    iterIA = 0;
  }
  double it_val(){ return A[iterA]; }
  int it_col(){ return JA[iterA]; }
  int it_row(){
    while(iterA >= IA[iterIA])
      ++iterIA;
    return iterIA-1;
  }
  void next(){ ++iterA; }
  bool end(){ return iterA == nnz; }


  void insert(double v, int g_col, int g_row){ //global row / col - this is child matrix
    debug_s("inserting");
    JA[iterA+1] = g_col - first_col;
    A[iterA+1] = v;

    int last_row = it_row();
    int new_row = g_row - first_row;
    if(last_row == new_row){
      ++IA[iterIA];
    } else {
      IA[iterIA + 1] = IA[iterIA] + 1;
      iterIA++;
    }
    iterA++;


    debug_s("inserted");
  }

  size_t csr_size(){
    return sizeof(int)*(5 + row_no + 1 + nnz) + sizeof(double)*(nnz);
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



  MPI_Request send_req, recv_req;
  void send(int rank){
    send_req = MPI_Isend(csr, csr_size, MPI_BYTE, rank, 1410, MPI_COMM_WORLD);
  }
  MPI_Request Sparse::recv(int rank, void* csr, int row_no_max, int nnz_max){
    MPI_Request recv_req = MPI_Irecv(
        csr,
        csr_alloc_size(row_no_max, nnz_max),
        MPI_BYTE,
        rank, 1410, MPI_COMM_WORLD);
  }

  private:
  int iterA, iterIA;
  int nnz_max, row_no_max;
};

#endif
