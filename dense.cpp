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
}
Dense::Dense(int row_no, int col_no, int first_row, int first_col, int seed){
  data = new double[row_no*col_no];
  this->row_no = row_no;
  this->col_no = col_no;
  this->first_row = first_row;
  this->first_col = first_col;
  for(int r=0; r<row_no; ++r)
    for(int c=0; c<col_no; ++c)
      *my_val(r,c) = generate_double(seed, r + first_row, c + first_col);

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
  assert(result < row_no);
  return result;
}

int Dense::my_col(int g_col){
  int result = g_col - first_col;
  assert(result < col_no);
  return result;
}
