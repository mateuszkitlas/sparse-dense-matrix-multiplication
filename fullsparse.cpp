#include "fullsparse.h"
#include "common.h"
#include <stdlib.h>
#include <assert.h>


//divides the most equally
inline int block_size(int matrix_size, int block_count, int block_no){
  return matrix_size / block_count + (matrix_size % block_count == block_no + 1);
}

inline int max_block_size(int matrix_size, int block_count){
  return block_size(matrix_size, block_count, 0);
}
inline int which_block(int matrix_size, int block_count, int matrix_i){
  int bigger_blocks_count = matrix_size % block_count;
  int block_size = matrix_size / block_count; //+1 in first bigger_blocks_count blocks
  int i2 = matrix_i - (block_size + 1) * bigger_blocks_count;
  if(i2 < 0)
    return matrix_i / (block_size + 1);
  else
    return (bigger_blocks_count - 1) + i2 / block_size;
}

void FullSparse::init_split(bool by_col, int block_count){
  this->by_col = by_col;
  this->block_count = block_count;


  int* nnzs = new int[block_count];
  int nnz_max = 0;

  for(begin(); !end(); next()){
    int block_no = which_block(side(), block_count, by_col ? it_col() : it_row());
    nnzs[block_no]++;
  }
  for(int block_no=0; block_no<block_count; ++block_no)
    nnz_max = std::max(nnz_max, nnzs[block_no]);

  split_nnzs = nnzs;
  split_nnz_max = nnz_max;
  split_row_no_max = by_col ? side() : max_block_size(side(), block_count);
}

Sparse** FullSparse::split(){
  Sparse** children = new Sparse*[block_count];

  int first_incl = 0; //first col/row in block
  int this_block_size;

  for(int block_no=0; block_no<block_count; ++block_no){
    this_block_size = block_size(side(), block_count, block_no);

    children[block_no] = Sparse::create(
      split_row_no_max,
      split_nnz_max,
      by_col ? 0 : first_incl, //first row
      by_col ? first_incl : 0, //first col
      by_col ? side() : this_block_size, //row_no
      by_col ? this_block_size : side(), //col_no
      split_nnzs[block_no] //nnz
    );

    children[block_no]->block_no = block_no;

    first_incl += this_block_size;
  }

  debug_s("partition data");
  while(!this->end()){
    int block_no = which_block(side(), block_count, by_col ? it_col() : it_row());
    debug_d(block_no);
    children[block_no]->insert(it_val(), it_col(), it_row());
    next();
  }
  return children;
}

FullSparse* FullSparse::create(
    int row_no,
    int col_no,
    int nnz){
  void* csr = Sparse::csr_alloc(row_no, nnz);

  *( ((int*)csr) + 2 ) = row_no;
  *( ((int*)csr) + 4 ) = nnz;

  FullSparse* sp = new FullSparse(csr, row_no, nnz);
  sp->first_row = 0;
  sp->first_col = 0;
  sp->col_no = col_no;

  sp->row_no_max = row_no;
  sp->nnz_max = nnz;
  sp->IA[0]=0;
  sp->IA[1]=0;
  sp->JA[0]=0;
  sp->JA[1]=0;
  sp->A[0]=0;

  assert(sp->first_row == 0);
  assert(sp->first_col == 0);
  assert(sp->row_no == row_no);
  assert(sp->col_no == col_no);
  assert(sp->nnz == nnz);
  assert(sp->IA < sp->JA);
  assert((void*)sp->JA < (void*)sp->A);


  return sp;
}



int FullSparse::side(){
  assert(row_no == col_no);
  return col_no;
}
