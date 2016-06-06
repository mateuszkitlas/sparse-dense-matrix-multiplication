#ifndef FULLSPARSE_H
#define FULLSPARSE_H

#include "sparse.h"


class FullSparse: public Sparse {
  public:

  bool by_col;
  int block_count;
  int* split_nnzs;

  int side(); //row_no; asserts row_no == col_no

  void init_split(bool by_col, int block_count);
  Sparse** split();
  static FullSparse* create(int row_no, int col_no, int nnz);

  using Sparse::Sparse;

};

#endif
