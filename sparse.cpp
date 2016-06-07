#include "sparse.h"
#include <stdlib.h>
#include "common.h"

void Sparse::update_refs(){
  IA = ((int*)csr) + 5;
  JA = IA + row_no + 1;
  A = (double*)(JA + nnz + 1);
}
size_t Sparse::csr_size(){
  return sizeof(int)*(5 + (row_no + 1) + (nnz + 1)) + sizeof(double)*(nnz);
}
size_t Sparse::csr_alloc_size(
      int row_no_max,
      int nnz_max){//per row
  //meta - first_row, first_col, row_no, col_no, nnz - int[5]
  //IA - int[row_no_max + 1]
  //JA - int[matrix_size + 1]
  //A - double[matrix_size]
  return
      sizeof(int)*( 5 + (row_no_max + 1) + (nnz_max + 1) )+
      sizeof(double)*( nnz_max );
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
  //sp->A[0]=0;

  assert(sp->first_row == first_row);
  assert(sp->first_col == first_col);
  assert(sp->row_no == row_no);
  assert(sp->col_no == col_no);
  assert(sp->nnz == nnz);
  assert(&sp->nnz + 1 == sp->IA);
  assert(sp->IA < sp->JA);
  assert((void*)sp->JA <= (void*)sp->A);

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
  send_req = MPI_REQUEST_NULL;
  recv_req = MPI_REQUEST_NULL;
  assert(csr == this->csr);
  this->row_no_max = row_no_max;
  this->nnz_max = nnz_max;
  update_refs();
  begin();
  block_no = -1;
}

void Sparse::free_csr(){ free(csr); }

void Sparse::printA(){
  char x[10000];
  begin();
  for(int r=0; r<row_no; ++r){
    for(int c=0; c<col_no; ++c){
      if(it_row() == r && it_col() == c){
        sprintf(x, "%s%.2lf ", x, it_val());
        next();
      }
      else
        sprintf(x, "%s     ", x);
    }
    sprintf(x, "%s\n", x);
  }
  printf("%s", x);
  assert(end());
}

void Sparse::print(){
#ifdef DEBUG
  /*
  printf("\
%d   block_no = %d\n\
    first_row = %d\n\
    first_col = %d\n\
    row_no = %d\n\
    col_no = %d\n\
    nnz = %d\n\n"
mpi_rank, block_no,
    first_row,
    first_col,
    row_no,
    col_no,
    nnz,
    );
  */
  /*
  printf("JA: ");
  for(int i=0; i<=nnz; ++i)
    printf("%d ", JA[i]);
  printf("\n");
  printf("IA: ");
  for(int i=0; i<=row_no; ++i)
    printf("%d ", IA[i]);
  printf("\n");*/
  printA();
#endif
}

void Sparse::begin(){
  iterA = 0;
  iterIA = 0;
}

int Sparse::it_row(){
  while(iterA >= IA[iterIA])
    ++iterIA;
  return iterIA-1;
}
void Sparse::insert(double v, int g_col, int g_row){ //global row / col - this is child matrix
  //debug_s("inserting");
  JA[iterA] = g_col - first_col;
  A[iterA] = v;

  int last_row = iterIA-1;
  //debug_d(last_row);
  int new_row = g_row - first_row;
  if(last_row == new_row){
    ++IA[iterIA];
  } else {
    IA[iterIA + 1] = IA[iterIA] + 1;
    iterIA++;
  }
  iterA++;


  //debug_s("inserted");
}
