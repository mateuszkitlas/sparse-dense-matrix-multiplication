#ifndef SPARSE_H
#define SPARSE_H

#include <vector>

using namespace std;

class Sparse {
  public:

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


  void insert(double v, int g_col, int g_row){ //global row / col - this is child matrix
    int last_row = it_row();
    int new_row = g_row - first_row;
    if(last_row == new_row){
      ++IA[iterIA];
    } else {
      IA[iterIA + 1] = IA[iterIA] + 1;
      iterIA++;
    }
    iterA++;

    A[iterA] = v;
    JA[iterA] = g_col - first_col;

    nnz++;
  }

  static Sparse* create(
      int row_no_max,
      int nnz_max,
      int first_row,
      int first_col,
      int row_no,
      int col_no,
      int nnz);
  Sparse(void* csr);
  void free_csr();
  int side(); //row_no; asserts row_no == col_no

  private:
  int iterA, iterIA;
  int nnz_max, row_no_max;
};

#endif
