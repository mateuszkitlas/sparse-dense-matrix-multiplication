#include "compute.h"

void multiply(Sparse* a, Dense* b, Dense* c){
  int last_col_excl = b->col_no + b->first_col;
  SPFOR(a){
    for(int col_i = b->first_col; col_i < last_col_excl; ++col_i){
      *c->val(a->row(), col_i) += a->val() * *b->val(a->col(), col_i);
    }
  }
}


