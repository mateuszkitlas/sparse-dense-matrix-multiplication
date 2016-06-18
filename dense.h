#ifndef DENSE_H
#define DENSE_H

#include <stdio.h>
#include <cstdio>
#include "common.h"
#include <mpi.h>

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
    int block_no;
    Dense(int row_no, int col_no, int first_row, int first_col);
    Dense(int block_row_no, int block_col_no);
    Dense(int block_no);
    ~Dense();
    void print();
    void generate(int seed);
    void zero();
    void alloc();
    void final_send(int x, MPI_Request *r);

    double* val(int g_row, int g_col);

    void final_send();
    void final_recv(int x, MPI_Request *r);
    void _send(int rank);
    void _recv(int rank);
    void wait();
    void recv(int rank, MPI_Request* mpi_recv);
    int ge_elements(double ge);


    MPI_Request send_req, recv_req;
  private:

    inline double* my_val(int row, int col);
    inline int my_row(int g_row);
    inline int my_col(int g_col);
};

#endif
