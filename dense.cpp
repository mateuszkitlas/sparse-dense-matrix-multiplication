#include <stdio.h>
#include <cstdio>
#include "dense.h"
#include "densematgen.h"

Dense::Dense(int row_no, int col_no, int first_row, int first_col){
  data = new double[row_no*col_no]();
  this->row_no = row_no;
  this->col_no = col_no;
  this->first_row = first_row;
  this->first_col = first_col;
#ifdef DEBUG
  this->name = 'c';
#endif
}

int Dense::ge_elements(double ge){
  int n=0;
  for(int r=0; r<row_no; ++r)
    for(int c=0; c<col_no; ++c)
      n += *my_val(r,c) >= ge;
  return n;
}

Dense::Dense(int row_no, int col_no, int first_row, int first_col, int seed){
  data = new double[row_no*col_no];
  this->row_no = row_no;
  this->col_no = col_no;
  this->first_row = first_row;
  this->first_col = first_col;
  for(int r=0; r<row_no; ++r)
    for(int c=0; c<col_no; ++c){
#ifdef IDENTITY_MATRIX
      *my_val(r,c) = r+first_row==c+first_col;
#else
      *my_val(r,c) = generate_double(seed, r + first_row, c + first_col);
#endif
    }
#ifdef DEBUG
  this->name = 'b';
#endif
  //print();
}


Dense::Dense(bool by_col, int block_no) : Dense(
    by_col ? mpi_meta_init.side : block_size(block_no),
    by_col ? block_size(block_no) : mpi_meta_init.side,
    first_side(false, by_col, block_no),
    first_side(true, by_col, block_no)
  ){
  if(by_col){
    assert(first_row == 0);
    assert(row_no == mpi_meta_init.side);
  }
  else {
    assert(first_col == 0);
    assert(col_no == mpi_meta_init.side);
  }
}


Dense::~Dense(){
  delete data;
}

double* Dense::val(int g_row, int g_col){
  return my_val(my_row(g_row), my_col(g_col));
}

double* Dense::my_val(int row, int col){
  return data + row*col_no + col;
}


int Dense::my_row(int g_row){
  int result = g_row - first_row;
#ifdef DEBUG
  if(result >= row_no || first_row > g_row){
    debug_s(&name);
    debug_d(first_row);
    debug_d(g_row);
    debug_d(result);
    debug_d(row_no);
  }
#endif
  assert(result < row_no);
  assert(first_row <= g_row);
  return result;
}

int Dense::my_col(int g_col){
  int result = g_col - first_col;
  //if(result >= col_no)
  //  printf("%d %d %d\n", g_col, first_col, col_no);
  assert(result < col_no);
  assert(first_col <= g_col);
  return result;
}

void Dense::print(){
#ifdef DEBUG
  char x[10000];
  sprintf(x, "--| %c %d |-----\n", name, mpi_no(0));
  for(int r=0; r<row_no; ++r){
    for(int c=0; c<col_no; ++c){
      sprintf(x, "%s%.2lf ", x, *my_val(r,c));
    }
    sprintf(x, "%s\n", x);
  }
  sprintf(x, "%s-----------\n", x);
  fprintf(stderr, "%s", x);
#endif
}

void Dense::send(){
  MPI_Request r;
  MPI_Isend(data, row_no*col_no, MPI_DOUBLE, 0, 2000, MPI_COMM_WORLD, &r);
  MPI_Wait(&r);
  debug("done");
}

void Dense::recv(int rank, MPI_Request* mpi_recv){
  MPI_Irecv(data, row_no*col_no, MPI_DOUBLE, rank, 2000, MPI_COMM_WORLD, mpi_recv);
}
