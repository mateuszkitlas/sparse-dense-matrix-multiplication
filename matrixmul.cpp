#include <stdlib.h>
#include <getopt.h>

#include "densematgen.h"
#include "dense.h"
#include "sparse.h"
#include "fullsparse.h"
#include "common.h"
#include "compute.h"

int main(int argc, char * argv[])
{
  int show_results = 0;
  int use_inner = 0;
  int gen_seed = -1;

  int option = -1;

  double comm_start = 0, comm_end = 0, comp_start = 0, comp_end = 0;
  int exponent = 1;
  double ge_element = 0;
  int count_ge = 0;

  Sparse** sparses = NULL;
  Sparse* my_sparse = NULL;
  Dense* dense_b = NULL;
  Dense* dense_c = NULL;
  FullSparse* full_sparse = NULL; //only for coordinator
  bool by_col = true; //by_col == true -> split column as in "colmn A alg"

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  assert(mpi_rank == mpi_no(0));


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
        mpi_meta_init.side = col_no;

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
  init_block_count(by_col);


  //------------------------
  //---  scatter
  //------------------------
  debug("scatter barrier");
  MPI_Barrier(MPI_COMM_WORLD);
  comm_start = MPI_Wtime();
  MPI_Request mpi_meta_init_req;
  if(mpi_rank == 0){

    full_sparse->init_split(by_col);

    debug("coordinator ibcast");
    MPI_Ibcast(&mpi_meta_init, sizeof(mpi_meta_init), MPI_BYTE, 0, MPI_COMM_WORLD, &mpi_meta_init_req);

    Sparse** mini_sparses = full_sparse->split();
    full_sparse->free_csr();
    delete full_sparse;

    my_sparse = mini_sparses[0];
    //sparse->print();
    for(int block_no=1; block_no<block_count; ++block_no){
      Sparse *sp = mini_sparses[block_no];
      //sp->print();
      sp->send(block_no);
    }

    debug("coordinator broadcast wait");

    MPI_Wait(&mpi_meta_init_req);
    for(int block_no=1; block_no<block_count; ++block_no){
      Sparse *sp = mini_sparses[block_no];
      sp->send_wait();
      sp->free_csr();
      delete sp;
    }
    delete mini_sparses;
  } else {
    debug("waiting for broadcast from coordinator");
    MPI_Ibcast(&mpi_meta_init, sizeof(mpi_meta_init), MPI_BYTE, 0, MPI_COMM_WORLD, &mpi_meta_init_req);
    MPI_Wait(&mpi_meta_init_req);
    debug("got broadcast from coordinator");
    my_sparse = Sparse::mpi_create();
    my_sparse->recv(0, mpi_no(0));
  }


  sparses = new Sparse*[repl_fact];
  my_sparse->recv_wait();
  debug("has my_sparse");
  //my_sparse->print();
  sparses[0] = my_sparse;
  for(int ci = 1; ci<repl_fact; ci++){
    sparses[ci] = Sparse::mpi_create();
    sparses[ci]->recv( mpi_no(ci), mpi_no(ci) );
    my_sparse->send( mpi_no(-ci) );
  }

  //------------------------
  //---  dense inits
  //------------------------
  int my_block_no = mpi_no(0);

  {
    if(by_col)
      assert((num_processes % repl_fact) == 0);
    else
      assert((num_processes % (repl_fact*repl_fact)) == 0);

    //debug_d(my_block_no);
    //debug_d(block_size(my_block_no));
    int
      row_no = by_col ? mpi_meta_init.side : block_size(my_block_no),
      col_no = by_col ? block_size(my_block_no) : mpi_meta_init.side,
      first_row = first_side(false, by_col, my_block_no),
      first_col = first_side(true, by_col, my_block_no);
    if(by_col){
      assert(first_row == 0);
      assert(row_no == mpi_meta_init.side);
    }
    else {
      assert(first_col == 0);
      assert(col_no == mpi_meta_init.side);
    }
    dense_b = new Dense(row_no, col_no, first_row, first_col, gen_seed);
    dense_c = new Dense(row_no, col_no, first_row, first_col);
  }




  /*for(int ci = 0; ci<repl_fact; ci++){
    sparses[ci]->recv_wait();
    sparses[ci]->send_wait();
  }*/
  MPI_Barrier(MPI_COMM_WORLD);
  comm_end = MPI_Wtime();

  //------------------------
  //---  compute
  //------------------------
  comp_start = MPI_Wtime();
  MPI_Barrier(MPI_COMM_WORLD);


  int ci=0;
  int done_blocks=0;
  int done_nothing=0;

  while(done_blocks < block_count){
    Sparse *sp = sparses[ci];
    bool wait = false;
    if(done_nothing >= repl_fact*2)
      wait = true;
    if(!sp->done_multiplication && sp->send_counter < num_processes){
      sp->recv_wait();
      if(sp->recv_ready()){
        sp->send();
        multiply(sp, dense_b, dense_c);
        done_blocks++;
        debug_d(sp->block_no);
        done_nothing = 0;
      } else done_nothing++;
    } else done_nothing++;
    if(sp->done_multiplication && sp->send_ready() && sp->recv_counter < num_processes){
      sp->recv();
      done_nothing = 0;
    } else done_nothing++;
    ci = (ci+1) % repl_fact;
  }

  for(int ci=0; ci<repl_fact; ++ci){
    Sparse *sp = sparses[ci];
    debug_d(sp->block_no);
    if(!sp->send_ready())
      sp->send_wait();
  }


  comp_end = MPI_Wtime();


  //------------------------
  //---  free A
  //------------------------
  for(int ci = 0; ci<repl_fact; ci++)
    sparses[ci]->free_csr();
  //my_sparse->print();
  delete sparses;
  delete dense_b;

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

  delete dense_c;


  MPI_Finalize();
  return 0;
}
