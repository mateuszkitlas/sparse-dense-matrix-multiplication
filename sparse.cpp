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

//divides the most equally
inline int part_size(int a, int b, int i){ return a / b + (a % b == i + 1) }

int Sparse::nnz_max_by_col(parts){
  int nnzs = new int[parts];
  int m = 0;
  for(int i=0; i<parts; ++i)
    nnzs[i]=0;
  for(int i=0; i<*this->nnz; ++i)
    m = max(m, ++nnzs[this->JA[i]]);
  return m;
}
int Sparse::nnz_max_by_row(parts){
  int m = 0;
  int last_row = 0; //excluded
  int this_row = 0; //included
  for(int p=0; p<parts; ++p){
    this_row += part_size(*this->row_no, parts, p);
    m = max(m, IA[this_row]-IA[last_row]);
    last_row = this_row;
  }
  return m;
}

Sparse* Sparse::split(bool by_row, int c, int p){
  //row_no == col_no
  int row_no = *this->row_no;
  int part_count = p/c + (p%c > 0);
  Sparse* result = new Sparse*[row_no];
  for(int row; row<row_no; ++row){
    result[row]
    first_row = part_size(row_no, rows, row);
  }
  return result;
}

