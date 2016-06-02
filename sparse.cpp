#include "sparse.h"
#include <stdlib.h>

void* Sparse::csr_alloc(
      int row_no_max,
      int nnz_max){
  //meta - first_row, first_col, row_no, col_no, nnz - int[5]
  //IA - int[row_no_max+1]
  //JA - int[nnz_max]
  //A - double[nnz_max]
  return malloc(
      sizeof(int)*(5 + row_no_max + 1 + nnz_max)
      + sizeof(double)*(nnz_max));
}

Sparse::~Sparse(){
  free(csr);
}

Sparse::Sparse(
      int first_row,
      int first_col,
      int row_no_max,
      int col_no_max,
      int nnz_max){
  this->csr = Sparse::csr_alloc(row_no_max, nnz_max);
  this->first_row = ((int*)csr) + 0;
  this->first_col = ((int*)csr) + 1;
  this->row_no = ((int*)csr) + 2;
  this->col_no = ((int*)csr) + 3;
  this->nnz = ((int*)csr) + 4;
  this->IA = ((int*)csr) + 5;
  this->JA = ((int*)csr) + 5 + row_no_max + 1;
  this->A = (double*)(((int*)csr) + 5 + row_no_max + 1 + nnz_max);

  *this->first_row = first_row;
  *this->first_col = first_col;
}

