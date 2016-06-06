#ifndef FULLSPARSE_H
#define FULLSPARSE_H

#include "sparse.h"


int mpi_meta_init[2] = {-1, -1};
int mpi_meta_init_size = 2;
int &split_row_no_max = mpi_meta_init[0];
int &split_nnz_max = mpi_meta_init[1];

class FullSparse: public Sparse {
  public:

  //TODO init!
  bool by_col;
  int block_count;
  int split_nnz_max;
  int split_row_no_max;
  int* split_nnzs;

  int side(); //row_no; asserts row_no == col_no

  void init_split(bool by_col, int block_count);
  Sparse** split();
  FullSparse* create(int row_no, int col_no, int nnz);

};

#endif
