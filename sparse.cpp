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
  send_counter = 0;
  recv_counter = 0;
  assert(csr == this->csr);
  this->row_no_max = row_no_max;
  this->nnz_max = nnz_max;
  update_refs();
  begin();
  block_no = -1;
  done_multiplication = false;
}

void Sparse::free_csr(){ free(csr); }

void Sparse::printA(char* x){
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
  assert(end());
}

void Sparse::print(){
#ifdef DEBUG
  char x[100000];
  
  sprintf(x, "\
%d   block_no = %d\n\
    first_row = %d\n\
    first_col = %d\n\
    row_no = %d\n\
    col_no = %d\n\
    nnz = %d\n\n",
mpi_rank, block_no,
    first_row,
    first_col,
    row_no,
    col_no,
    nnz
    );
  
  
  sprintf(x, "%sJA: ", x);
  for(int i=0; i<=nnz; ++i)
    sprintf(x, "%s%d ", x, JA[i]);
  sprintf(x, "%s\nIA: ", x);
  for(int i=0; i<=row_no; ++i)
    sprintf(x, "%s%d ", x, IA[i]);
  sprintf(x, "%s\n", x);
  
  printA(x);
  fprintf(stderr, "\n%s\n", x);
#endif
}

void Sparse::begin(){
  iterA = 0;
  iterIA = 0;
}

void Sparse::test(){
#ifdef IDENTITY_MATRIX
  int *per_row = new int[row_no](), *per_col = new int[col_no]();
  SPFOR(this){
    per_row[it_row()]++;
    per_col[it_col()]++;
  }
  assert(IA[row_no-1]>=IA[0]);
  for(int ia=1; ia<row_no; ia++)
    assert(IA[ia] - IA[ia-1] == per_row[ia]);
  for(int ja=0; ja<nnz; ja++)
    per_col[JA[ja]]--;
  for(int c=0; c<col_no; ++c)
    assert(per_col[c] == 0);
#endif
}

int Sparse::it_row(){
  while(iterA >= IA[iterIA])
    ++iterIA;
  return iterIA;
}
void Sparse::done_insert(){
  for(; iterIA<row_no; ++iterIA){
    IA[iterIA + 1] = IA[iterIA];
  }
  begin();
}

void Sparse::insert(double v, int g_col, int g_row){ //global row / col - this is child matrix
  //debug_s("inserting");
  JA[iterA] = g_col - first_col;
  A[iterA] = v;

  int last_row = iterIA;
  //debug_d(last_row);
  int new_row = g_row - first_row;
  for(int i=last_row; i<new_row; ++i){
    IA[iterIA + 1] = IA[iterIA];
    iterIA++;
  }
  ++IA[iterIA];
  iterA++;


  //debug_s("inserted");
}
