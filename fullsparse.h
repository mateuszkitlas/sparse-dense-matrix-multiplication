#ifndef FULLSPARSE_H
#define FULLSPARSE_H

#include "sparse.h"


class FullSparse: public Sparse {
  public:

  bool by_col;
  int* split_nnzs;

  void init_split(bool by_col);
  Sparse** split();
  static FullSparse* create(int row_no, int col_no, int nnz);

  using Sparse::Sparse;

  ~FullSparse();

};

#endif
