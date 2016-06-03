#ifndef SPARSE_H
#define SPARSE_H

#include <vector>

using namespace std;

class Sparse {
  public:
  typedef vector<pair<int,double> > Row;
  vector<Row> rows;
  void* csr;
  int *first_row, *first_col, *row_no, *col_no, *nnz, *IA, *JA;
  double *A;
  ~Sparse();
  static void* csr_alloc(
      int row_no_max,
      int nnz_max);
  Sparse(
      int first_row,
      int first_col,
      int row_no,
      int col_no,
      int nnz_max);
  Sparse** split(int rows, cols);
  inline int r(){ return *this->first_row

};

#endif
