#ifndef DENSE_H
#define DENSE_H

#include <stdio.h>
#include <cstdio>
#include "common.h"


class Dense {
  public:

#ifdef DEBUG
    char name;
#endif

  double* data;
    int row_no;
    int col_no;
    int first_row;
    int first_col;
    Dense(int row_no, int col_no, int first_row, int first_col);
    Dense(int row_no, int col_no, int first_row, int first_col, int seed);
    ~Dense();
    void print();

    double* val(int g_row, int g_col);

  private:

    inline double* my_val(int row, int col);
    inline int my_row(int g_row);
    inline int my_col(int g_col);
};

#endif
