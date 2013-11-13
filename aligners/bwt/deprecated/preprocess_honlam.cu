#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sequence.h>
#include <thrust/sort.h>
#include <thrust/gather.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sys/time.h>

#include "BW_preprocess.h"
#include "BW_io.h"

struct cond_vec_inc
{
  unsigned int rank;

  __host__ __device__
  void operator()(unsigned int &x)
  {
    if (x>=rank) x = x + 1;
  }
};

int main(int argc, char **argv)
{

  byte_vector Xorig, X, B, Bi;
  vector C, C1;
  comp_vector S, Si;
  comp_matrix O, Oi;
  comp_vector Scomp, Scompi;

  int s_ratio;

  timevars();
  initReplaceTable();

  if (argc!=4) {
    printf("Syntax:\n\t%s ref_file output_dir s_ratio\n", argv[0]);
    return 1;
  }

  s_ratio = atoi(argv[3]);

  tic("Loading duplicated reference for BWT calculation");
  load_duplicated_reference(&Xorig, &X, argv[1]);
  toc();

  printUIntVector(X.vector, X.n);
  printUIntVector(X.vector, X.n*2+1);

  tic("Calculating BWT from reference string");

  size_t distance, size = X.n+1;
  unsigned int aux, rank = 0;
  char c;

  cond_vec_inc inc_op;

  thrust::host_vector<unsigned int> h_F(size,0);
  thrust::host_vector<unsigned int> h_H(nA + 1, 0);

  printUIntVector(h_H.data(), h_H.size());
  printUIntVector(h_F.data(), h_F.size());

  //printf("*********************************************************\n");

  size_t i = size;

  do {    
    i--;

    c = X.vector[i];

    //1 - Compute the rank of T'
    thrust::host_vector<unsigned int>::iterator it = thrust::lower_bound(h_F.begin() + h_H[c+1], h_F.begin() + h_H[c+2], rank);

    distance = std::distance(h_F.begin() + h_H[c+1], it);

    if (distance == 0)
      rank = h_H[c+1];
    else
      rank = std::distance(h_F.begin(), it);

    //printf("The rank is %u\n", rank);

    //2 - Update H->vector
    for(size_t j=0; j < h_H.size(); j++) {
      if (j > c+1) h_H[j]++;
    }
    printUIntVector(h_H.data(), h_H.size());

    //3 - Update F vector
    if (i != size-1) {
      aux = h_F[0];
      h_F[0] = rank;
      printUIntVector(h_F.data(), h_F.size());

      //thrust::copy(h_F.begin() + rank, h_F.begin() + size - i, h_F.begin() + rank + 1);
      h_F[rank] = aux;
      printUIntVector(h_F.data(), h_F.size());

      inc_op.rank=rank;
      //thrust::for_each(h_F.begin() + 1, h_F.begin() + size - i, inc_op);
      printUIntVector(h_F.data(), h_F.size());
    } else {
      h_F[0] = rank;
      printUIntVector(h_F.data(), h_F.size());
    }

    //printf("*********************************************************\n");
    //printf("%lu\n", i);
  } while(i!=0);

  toc();

  S.n = size;
  S.vector  = ( unsigned int *) malloc(S.n * sizeof(unsigned int));
  checkMalloc(S.vector, "main");

  B.n = size;
  B.vector  = ( char *) malloc(B.n * sizeof(char));
  checkMalloc(B.vector, "main");

  for(unsigned int k=0, p; k<S.n; k++) {
    p = h_F[p];
    S.vector[p] = k;
    //printf("k = %u, p = %u\n", k, p);
  }

  //thrust::gather(d_S, d_S + S->n, bwt_comp.BW_aux + X->n, d_B);
  
  printUIntVector(S.vector, S.n);
  /* printUIntVector(B.vector, B.n); */

  /* calculateScomp(&S, &Scomp, s_ratio); */
  /* printUIntVector(Scomp.vector, Scomp.siz); */

  /* saveUIntCompVector(&S, argv[2], "S"); */
  free(S.vector);
  /* saveUIntCompVector(&Scomp, argv[2], "Scomp"); */
  /* free(Scomp.vector); */

  /* tic("Calculating prefix-trie matrices C and O"); */
  /* calculateC(&C, &C1, &B, 0); */
  /* calculateO(&O, &B); */
  /* toc(); */

  /* printUIntVector(C.vector, C.n); */
  /* printUIntVector(C1.vector, C1.n); */
  /* printCompMatrix(O); */

  /* saveCharVector(&X, argv[2], "X"); */
  /* saveUIntVector(&C, argv[2], "C"); */
  /* free(C.vector); */
  /* saveUIntVector(&C1, argv[2], "C1"); */
  /* free(C1.vector); */
  /* saveCharVector(&B, argv[2], "B"); */
  /* free(B.vector); */
  /* saveCompMatrix(&O, argv[2], "O"); */
  /* freeCompMatrix(&O); */

  /* tic("Calculating BWT from inverted reference string"); */
  /* calculateBWT(&Bi, &Si, &X, 1); */
  /* toc(); */

  /* free(X.vector); */
  /* free(Xorig.vector); */

  /* printUIntVector(Bi.vector, Bi.n); */
  /* printUIntVector(Si.vector, Si.n); */

  /* calculateScomp(&Si, &Scompi, s_ratio); */
  /* printUIntVector(Scompi.vector, Scompi.siz); */

  /* saveUIntCompVector(&Si, argv[2], "Si"); */
  /* free(Si.vector); */
  /* saveUIntCompVector(&Scompi, argv[2], "Scompi"); */
  /* free(Scompi.vector); */

  /* tic("Calculating inverted prefix-trie matrix Oi"); */
  /* calculateO(&Oi, &Bi); */
  /* toc(); */

  /* printCompMatrix(Oi); */

  /* saveCharVector(&Bi, argv[2], "Bi"); */
  /* free(Bi.vector); */

  /* saveCompMatrix(&Oi, argv[2], "Oi"); */
  /* freeCompMatrix(&Oi); */

  /* return 0; */

}
