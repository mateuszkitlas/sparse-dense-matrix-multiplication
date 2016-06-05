#include "sparse.h"
#include <stdlib.h>
#include <assert.h>

size_t Sparse::csr_alloc_size(
      int row_no_max,
      int nnz_max){//per row
  //meta - first_row, first_col, row_no, col_no, nnz - int[5]
  //IA - int[row_no_max+1]
  //JA - int[matrix_size]
  //A - double[matrix_size]
  int matrix_size = nnz_max * row_no_max;
  return
      sizeof(int)*( 5 + row_no_max + 1 + matrix_size ) +
      sizeof(double)*( matrix_size );
}
void* Sparse::csr_alloc(
      int row_no_max,
      int nnz_max){//per row
  return malloc(Sparse::csr_alloc_size(row_no_max, nnz_max));
}

Sparse* Sparse::create(
    int row_no_max,
    int nnz_max,
    int first_row,
    int first_col,
    int row_no,
    int col_no,
    int nnz){
  void* csr = Sparse::csr_alloc(row_no_max, nnz_max);

  *( ((int*)csr) + 2 ) = row_no;
  *( ((int*)csr) + 4 ) = nnz;

  Sparse* sp = new Sparse(csr, row_no_max, nnz_max);
  sp->first_row = first_row;
  sp->first_col = first_col;
  sp->col_no = col_no;

  sp->row_no_max = row_no_max;
  sp->nnz_max = nnz_max;
  sp->IA[0]=0;
  sp->IA[1]=0;
  sp->JA[0]=0;
  sp->JA[1]=0;
  sp->A[0]=0;

  assert(sp->first_row == first_row);
  assert(sp->first_col == first_col);
  assert(sp->row_no == row_no);
  assert(sp->col_no == col_no);
  assert(sp->nnz == nnz);
  assert(sp->IA < sp->JA);
  assert((void*)sp->JA < (void*)sp->A);


  return sp;
}

Sparse::Sparse(void* csr, int row_no_max, int nnz_max)
    : csr( csr )
    , first_row( *( ((int*)csr) + 0 ) )
    , first_col( *( ((int*)csr) + 1 ) )
    , row_no(    *( ((int*)csr) + 2 ) )
    , col_no(    *( ((int*)csr) + 3 ) )
    , nnz(       *( ((int*)csr) + 4 ) )
  {
  assert(csr == this->csr);
  this->row_no_max = row_no_max;
  this->nnz_max = nnz_max;
  IA = ((int*)csr) + 5;
  JA = IA + row_no + 1;
  A = (double*)(JA + nnz);
  begin();
}

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

//int Sparse::nnz_max_by_col(int block_count){
//  int* nnzs = new int[block_count];
//  int m = 0;
//  for(int block_no=0; i<block_count; ++block_no)
//    nnzs[block_no]=0;
//  for(int i=0; i<this->nnz; ++i)
//    m = max(m, ++(nnzs[this->JA[i]]));
//  return m;
//}
//
//int Sparse::nnz_max_by_row(int block_count){
//  int m = 0;
//  int last_row = 0; //excluded
//  int this_row = 0; //included
//  for(int block_no=0; block_no<block_count; ++block_no){
//    this_row += block_size(this->row_no, block_count, block_no);
//    m = max(m, IA[this_row]-IA[last_row]);
//    last_row = this_row;
//  }
//  return m;
//}


Sparse** Sparse::split(bool by_col, int block_count){ //by_col == true -> split column as in "colmn A alg"
  debug_s("start sparse");
  Sparse** children = new Sparse*[block_count];
  int* children_nnz = new int*[block_count];
  int child_nnz_max = 0;

  for(begin(); !end(); next()){
    int block_no = which_block(side(), block_count, by_col ? it_col() : it_row());
    children_nnz[block_no]++;
  }
  for(int block_no=0; block_no<block_count; ++block_no)
    child_nnz_max = max(child_nnz_max, children_nnz[block_no]);

  //TODO tutaj chybba może być min(max_block_size, nnz_max)
  int child_row_no_max = by_col ? side() : max_block_size(side(), block_count);
  int first_incl = 0; //first col/row in block
  int this_block_size;
  //TODO MPI nnz_max
  //alloc child Sparses
  for(int block_no=0; block_no<block_count; ++block_no){
    debug_d(block_no);
    debug_d(block_count);
    this_block_size = block_size(side(), block_count, block_no);
    debug_d(this_block_size);
    debug_d(child_row_no_max);
    debug_d(child_nnz_max);

    children[block_no] = Sparse::create(
      child_row_no_max,
      child_nnz_max,
      by_col ? 0 : first_incl, //first row
      by_col ? first_incl : 0, //first col
      by_col ? side() : this_block_size, //row_no
      by_col ? this_block_size : side(), //col_no
      children_nnz[block_no] //nnz
    );

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

void Sparse::free_csr(){ free(csr); }
int Sparse::side(){
  assert(row_no == col_no);
  return col_no;
}

void Sparse::print(){
#ifdef DEBUG
  printf("\
    first_row = %d\n\
    first_col = %d\n\
    row_no = %d\n\
    col_no = %d\n\
    nnz = %d\n\n", first_row, first_col, row_no, col_no, nnz);
#endif
}
