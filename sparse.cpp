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


static Sparse& create(
    int first_row,
    int first_col,
    int row_no,
    int col_no,
    int nnz_max){
  return Sparse(
    first_row,
    first_col,
    row_no,
    col_no,
    nnz_max,
    Sparse::csr_alloc(row_no_max, nnz_max)
  );
}




Sparse::Sparse(
    int first_row,
    int first_col,
    int row_no_max,
    int col_no_max,
    int nnz_max,
    void* csr)
    : first_row( *( ((int*)csr) + 0 ) )
    , first_col( *( ((int*)csr) + 1 ) )
    , row_no(    *( ((int*)csr) + 2 ) )
    , col_no(    *( ((int*)csr) + 3 ) )
    , nnz(       *( ((int*)csr) + 4 ) )
    , IA( ((int*)csr) + 5 )
    , JA( ((int*)csr) + row_no_max + 1 )
    , IA( ((int*)csr) + row_no_max + 1 + nnz_max )
  {
  this->first_row = first_row;
  this->first_col = first_col;

  it_begin();
}

int Sparse::nnz_max_by_col(int block_count){
  int nnzs = new int[block_count];
  int m = 0;
  for(int i=0; i<block_count; ++i)
    nnzs[i]=0;
  for(int i=0; i<this->nnz; ++i)
    m = max(m, ++nnzs[this->JA[i]]);
  return m;
}

int Sparse::nnz_max_by_row(int block_count){
  int m = 0;
  int last_row = 0; //excluded
  int this_row = 0; //included
  for(int block_no=0; block_no<block_count; ++block_no){
    this_row += block_size(this->row_no(), block_count, block_no);
    m = max(m, IA[this_row]-IA[last_row]);
    last_row = this_row;
  }
  return m;
}

//divides the most equally
inline int block_size(int a, int b, int i){ return a / b + (a % b == i + 1) }
inline int which_block(int a, int b, int i){
  int c = a % b;
  int i2 = i - (b + 1) * c;
  if(i2 < 0)
    return i / (b + 1);
  else
    return c + i2 / b;
}

static int Sparse::block_count(int p, int c){ return p/c + (p%c > 0); }

Sparse** Sparse::split(bool by_col, int block_count){ //by_col == true -> split column as in "colmn A alg"
  //row_no == col_no
  Sparse** result = new Sparse*[block_count];
  int child_nnz_max = by_col ? this->nnz_my_by_col(block_count) : this->nnz_my_by_row(block_count);
  int first_incl = 0; //first col/row in block
  int this_block_size;
  //TODO MPI nnz_max
  //alloc child Sparses
  for(int i; i<block_count; ++i){
    this_block_size = block_size(this->row_no(), block_count, i);

    result[i] = by_col
      ? new Sparse( 0, first_incl, this->row_no(), this_block_size, child_nnz_max )
      : new Sparse( first_incl, 0, this_block_size, this->row_no(), child_nnz_max );

    first_incl += this_block_size;
  }

  //spread data
  while(!this->end()){
    Sparse *sp = result[ which_part( by_col ? it_col() : it_row() ) ];
    sp->insert(it_val(), it_col(), it_row());
    next();
  }
  return result;
}

