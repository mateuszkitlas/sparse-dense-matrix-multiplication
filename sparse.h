#ifndef SPARSE_H
#define SPARSE_H

#include <vector>

using namespace std;

class Sparse {
  public:
  void* csr;
  int *first_row, *first_col, *row_no, *col_no, *nnz, *IA, *JA;
  double *A;

  int iterA, iterIA;


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
  Sparse* split(int rows, cols);
  void it_begin(){
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
  bool next(){ ++iterA; return !end(); }
  bool end(){ return iterA == nnz; }


};

#endif
