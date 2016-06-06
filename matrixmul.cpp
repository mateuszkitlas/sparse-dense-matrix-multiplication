#include <stdlib.h>
#include <assert.h>
#include <getopt.h>

#include "densematgen.h"
#include "sparse.h"
#include "fullsparse.h"
#include "common.h"

int main(int argc, char * argv[])
{
  int show_results = 0;
  int use_inner = 0;
  int gen_seed = -1;
  int repl_fact = 1;

  int option = -1;

  double comm_start = 0, comm_end = 0, comp_start = 0, comp_end = 0;
  int num_processes = 1;
  int exponent = 1;
  double ge_element = 0;
  int count_ge = 0;

  Sparse* sparse = NULL;
  FullSparse* full_sparse = NULL; //only for coordinator
  bool by_col = true; //by_col == true -> split column as in "colmn A alg"

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  debug_d(mpi_rank);


  while ((option = getopt(argc, argv, "vis:f:c:e:g:")) != -1) {
    switch (option) {
    case 'v': show_results = 1; 
      break;
    case 'i':
      use_inner = 1;
      by_col = false;
      break;
    case 'f':
      //------------------------
      //---  load full_sparse
      //------------------------
      if ((mpi_rank) == 0)
      {
        FILE *file_sparse = fopen(optarg, "r");
        int row_no, col_no, nnz, nnz_max;
        fscanf(file_sparse, "%d%d%d%d", &row_no, &col_no, &nnz, &nnz_max);
        assert(col_no == row_no);

        full_sparse = FullSparse::create(row_no, col_no, nnz);
        for(int i=0; i<nnz; ++i)
          fscanf(file_sparse, "%lf", &full_sparse->A[i]);
        for(int i=0; i<(row_no+1); ++i)
          fscanf(file_sparse, "%d", &full_sparse->IA[i]);
        for(int i=0; i<nnz; ++i)
          fscanf(file_sparse, "%d", &full_sparse->JA[i]);
        debug_s("full_sparse loaded");
        full_sparse->print();
      }
      break;
    case 'c': repl_fact = atoi(optarg);
      break;
    case 's': gen_seed = atoi(optarg);
      break;
    case 'e': exponent = atoi(optarg);
      break;
    case 'g': count_ge = 1; 
      ge_element = atof(optarg);
      break;
    default: fprintf(stderr, "error parsing argument %c exiting\n", option);
      MPI_Finalize();
      return 3;
    }
  }
  if ((gen_seed == -1) || ((mpi_rank == 0) && (full_sparse == NULL)))
  {
    fprintf(stderr, "error: missing seed or sparse matrix file; exiting\n");
    MPI_Finalize();
    return 3;
  }



  //------------------------
  //---  scatter
  //------------------------
  debug("scatter barrier");
  MPI_Barrier(MPI_COMM_WORLD);
  comm_start = MPI_Wtime();
  MPI_Request mpi_meta_init_req;
  if(mpi_rank == 0){
    int block_count = num_processes;

    full_sparse->init_split(by_col, block_count);

    debug("coordinator ibcast");
    MPI_Ibcast(mpi_meta_init, mpi_meta_init_size, MPI_INT, 0, MPI_COMM_WORLD, &mpi_meta_init_req);

    Sparse** mini_sparses = full_sparse->split();
    full_sparse->free_csr();
    delete full_sparse;

    sparse = mini_sparses[0];
    //sparse->print();
    Sparse *sp;
    for(int block_no=1; block_no<block_count; ++block_no){
      sp = mini_sparses[block_no];
      //sp->print();
      sp->send(block_no);
    }

    debug("coordinator broadcast wait");

    MPI_Wait(&mpi_meta_init_req);
    for(int block_no=1; block_no<block_count; ++block_no){
      sp = mini_sparses[block_no];
      sp->send_wait();
      sp->free_csr();
      delete sp;
    }
    delete mini_sparses;
  } else {
    debug("waiting for broadcast from coordinator");
    MPI_Ibcast(mpi_meta_init, mpi_meta_init_size, MPI_INT, 0, MPI_COMM_WORLD, &mpi_meta_init_req);
    MPI_Wait(&mpi_meta_init_req);
    debug("got broadcast from coordinator");
    sparse = Sparse::mpi_create(split_row_no_max, split_nnz_max);
    sparse->recv(0, mpi_rank);
    sparse->recv_wait();
  }
  MPI_Barrier(MPI_COMM_WORLD);
  comm_end = MPI_Wtime();

  //------------------------
  //---  compute
  //------------------------
  comp_start = MPI_Wtime();
  MPI_Barrier(MPI_COMM_WORLD);
  comp_end = MPI_Wtime();


  //------------------------
  //---  free A
  //------------------------
  sparse->print();
  sparse->free_csr();
  delete sparse;

  if (show_results) 
  {
    // FIXME: replace the following line: print the whole result matrix
    printf("1 1\n42\n");
  }
  if (count_ge)
  {
    // FIXME: replace the following line: count ge elements
    printf("54\n");
  }


  MPI_Finalize();
  return 0;
}
