#include "fullsparse.h"
#include "common.h"
#include <stdlib.h>


inline int max_block_size(){ return block_size(0); }

inline int min_block_size(){ return block_size(block_count -1 ); }

inline int which_block(int matrix_i){
  int &side = mpi_meta_init.side;
  int bigger_blocks_count = side % block_count;
  int bigger_block_size = block_size(0);
  int i2 = matrix_i - bigger_block_size * bigger_blocks_count;
  if(i2 < 0)
    return matrix_i / bigger_block_size;
  else
    return bigger_blocks_count + i2 / (bigger_block_size - 1);
}

void FullSparse::init_split(bool by_col_){
  ASSERTS;
  this->by_col = by_col_;
  int* nnzs = new int[block_count]();

  int nnz_max = 0;

  for(begin(); !end(); next()){
    int block_i = which_block(by_col ? it_col() : it_row());
    nnzs[block_i]++;
  }
  nnz_max = *std::max_element(nnzs, nnzs + block_count);
  debug_d(nnz_max);

  this->split_nnzs = nnzs;
  ::mpi_meta_init.nnz_max = nnz_max;
  ::mpi_meta_init.row_no_max = by_col ? side() : max_block_size();
}

Sparse** FullSparse::split(){
  Sparse** children = new Sparse*[block_count];
  Sparse* sp;

  int first_incl = 0; //first col/row in block

  for(int block_i=0; block_i<block_count; ++block_i){
    debug_d(block_count);
    debug_d(::block_size(block_i));
    debug_d(split_nnzs[block_i]);
    debug_d(mpi_meta_init.side);

    sp = children[block_i] = Sparse::create(
      ::mpi_meta_init.row_no_max,
      ::mpi_meta_init.nnz_max,
      by_col ? 0 : first_incl, //first row
      by_col ? first_incl : 0, //first col
      by_col ? side() : ::block_size(block_i), //row_no
      by_col ? ::block_size(block_i) : side(), //col_no
      split_nnzs[block_i] //nnz
    );

    sp->block_no = block_i;

    first_incl += ::block_size(block_i);
  }

  debug("partition data");
  SPFOR(this){
    int block_i = which_block(by_col ? it_col() : it_row());
    sp = children[block_i];
    //debug_d(block_i);
    sp->insert(it_val(), it_col(), it_row());

    //assert(sp->it_val() == it_val());
    //assert(sp->it_col() <= it_col());
    //assert(sp->it_row() <= it_row());

    //debug_d(sp->nnz);

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
  assert((void*)sp->JA <= (void*)sp->A);


  return sp;
}



int FullSparse::side(){
  assert(row_no == col_no);
  assert(row_no == mpi_meta_init.side);
  return col_no;
}

FullSparse::~FullSparse(){
  delete split_nnzs;
}
