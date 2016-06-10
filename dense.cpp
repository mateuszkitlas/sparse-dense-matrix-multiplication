#include <stdio.h>
#include <cstdio>
#include "dense.h"
#include "densematgen.h"

void Dense::alloc(){
  data = new double[row_no*col_no];
}

void Dense::zero(){
  data = new double[row_no*col_no]();
}

Dense::Dense(int row_no, int col_no, int first_row, int first_col){
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

void Dense::generate(int seed){
  alloc();
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

Dense::Dense(int block_row_no, int block_col_no) : Dense(
    block_size(block_row_no),
    block_size(block_col_no),
    first_side(block_row_no),
    first_side(block_col_no)
  ){
}


Dense::Dense(int block_no) : Dense(
    side,
    block_size(block_no),
    0,
    first_side(block_no)
  ){
  assert(first_row == 0);
  assert(row_no == side);
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
