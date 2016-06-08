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
        fscanf(file_sparse, "%*d");
        for(int i=0; i<row_no; ++i)
          fscanf(file_sparse, "%d", &full_sparse->IA[i]);
        for(int i=0; i<nnz; ++i)
          fscanf(file_sparse, "%d", &full_sparse->JA[i]);
        debug_s("full_sparse loaded");
        //full_sparse->print();
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
    full_sparse->print();
    full_sparse->test();


#ifndef IDENTITY_MATRIX
    full_sparse->free_csr();
    delete full_sparse;
#endif
    my_sparse = mini_sparses[0];
    my_sparse->post_async_recv(0, 0);
    //sparse->print();
    for(int block_no=1; block_no<block_count; ++block_no){
      Sparse *sp = mini_sparses[block_no];
      //sp->print();
      sp->_send(block_no);
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
    my_sparse->_recv(0, mpi_no(0));
  }


  //------------------------
  //---  dense inits
  //------------------------
  int my_block_no = mpi_no(0);

  if(by_col)
    assert((num_processes % repl_fact) == 0);
  else
    assert((num_processes % (repl_fact*repl_fact)) == 0);

  int
    dense_row_no = by_col ? mpi_meta_init.side : block_size(my_block_no),
    dense_col_no = by_col ? block_size(my_block_no) : mpi_meta_init.side,
    dense_first_row = first_side(false, by_col, my_block_no),
    dense_first_col = first_side(true, by_col, my_block_no);
  if(by_col){
    assert(dense_first_row == 0);
    assert(dense_row_no == mpi_meta_init.side);
  }
  else {
    assert(dense_first_col == 0);
    assert(dense_col_no == mpi_meta_init.side);
  }

  dense_b = new Dense(dense_row_no, dense_col_no, dense_first_row, dense_first_col, gen_seed);




  MPI_Barrier(MPI_COMM_WORLD);
  comm_end = MPI_Wtime();


  //------------------------
  //---  compute
  //------------------------


  sparses = new Sparse*[repl_fact];
  my_sparse->recv_wait();
  assert(my_sparse->recv_ready());
  debug("has my_sparse");
  //my_sparse->print();
  sparses[0] = my_sparse;
  for(int ci = 1; ci<repl_fact; ci++){
    sparses[ci] = Sparse::mpi_create();
    sparses[ci]->_recv( mpi_no(ci), mpi_no(ci) );
    my_sparse->_send(mpi_no(-ci));
  }


  for(int ci = 0; ci<repl_fact; ci++){
    sparses[ci]->recv_wait();
    sparses[ci]->send_wait();
  }
  debug("done scatter");

  comp_start = MPI_Wtime();
  MPI_Barrier(MPI_COMM_WORLD);

  //------------------------
  //---  compute
  //------------------------


  for(int e=0; e<exponent; ++e){
    dense_c = new Dense(by_col, my_block_no);
    int ci=0;
    int done_blocks=0;
    int done_nothing=0;
    int sparse_cycles = block_count / repl_fact;

    int *cycles_done = new int[repl_fact]();

    while(done_blocks < block_count){
      Sparse *sp = sparses[ci];
      sp->recv_wait();
      sp->send();
      multiply(sp, dense_b, dense_c);
      if(++cycles_done[ci] < sparse_cycles){
        sp->send_wait();
        sp->recv();
      } else {
        assert(sp->recv_ready());
      }
      ci = (ci+1) % repl_fact;
      done_blocks++;
    }
    delete cycles_done;
    delete dense_b;
    dense_b = dense_c;
  }

  comp_end = MPI_Wtime();

  if (show_results)
  {
    debug("show results");
    if(mpi_rank == 0){
      Dense **whole_c = new Dense*[block_count];
      whole_c[0] = dense_c;
      MPI_Request *reqs = new MPI_Request[block_count-1];
      for(int block_i=1; block_i < block_count; ++block_i){
        Dense *den = new Dense(by_col, block_i);
        whole_c[block_i] = den;
        den->recv(block_i, &reqs[block_i-1]);
      }
      MPI_Waitall(block_count-1, reqs, MPI_STATUSES_IGNORE);
      delete reqs;
      int &side = mpi_meta_init.side;

#ifndef IDENTITY_MATRIX
      printf("%d %d\n", side, side);
      for(int r=0; r<side; ++r){
        for(int c=0; c<side; ++c){
          int block_i = which_block(by_col ? c : r);
          Dense *den = whole_c[block_i];
          printf("%*.5lf ", *den->val(r,c), 10);
        }
        printf("\n");
      }
#else
      Sparse *s = full_sparse;
      SPFOR(s){
        int block_i = which_block(by_col ? s->col() : s->row());
        Dense *den = whole_c[block_i];
        //fprintf(stderr, "%lf %lf %d [%d,%d]\n", *den->val(s->row(),s->col()), s->val(), block_i, s->row(), s->col());
        assert(*den->val(s->row(),s->col()) - s->val() < 0.01);
      }
      //fprintf(stderr, "cokolwiek\n");
      delete s;
#endif
      delete whole_c;
    } else {
      dense_c->send();
      delete dense_c;
    }
  }
  if (count_ge)
  {
    int my_ge_elements = dense_c->ge_elements(ge_element);
    int ge_elements = 0;
    MPI_Reduce(&my_ge_elements, &ge_elements, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if(mpi_rank == 0)
      printf("%d\n", ge_elements);
  }
  //------------------------
  //---  free mem
  //------------------------
  for(int ci=0; ci<repl_fact; ++ci){
    Sparse *sp = sparses[ci];
    sp->test();
    if(!sp->send_ready())
      sp->send_wait();
    sp->free_csr();
  }
  delete sparses;
  //delete dense_b;
  //------------------------



  MPI_Finalize();
  return 0;
}
