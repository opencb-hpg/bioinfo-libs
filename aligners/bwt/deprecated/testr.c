#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sys/time.h>

#include "commons/commons.h"
#include "commons/string_utils.h"
#include "BW_io.h"

int found;
unsigned int enW;
char search[MAXLINE];

unsigned int read_index=0;

int main(int argc, char **argv) {

  byte_vector B;
  vector C;
  comp_matrix O;
  comp_vector S, Scomp, R, Rcomp;

  if (argc!=3) {
    printf("Syntax:\n\t%s input_dir nucleotides\n", argv[0]);
    return 1;
  }

  timevars();
  initReplaceTable(argv[2]);

  tic("Loading data structures from input directory");

  readUIntVector(&C, argv[1], "C");
  readCompMatrix(&O, argv[1], "O");
  readCharVector(&B, argv[1], "B");
  readUIntCompVector(&S, argv[1], "S");
  readUIntCompVector(&Scomp, argv[1], "Scomp");
  readUIntCompVector(&R, argv[1], "R");
  readUIntCompVector(&Rcomp, argv[1], "Rcomp");

  toc();

  for (unsigned int i=0; i<S.n; i++) {

    if (getScompValue(i, &Scomp, &C, &O) != S.vector[i]) {
      printf("S diferent!\n");
    }

    if (getRcompValue(i, &Rcomp, &C, &O) != R.vector[i]) {
      printf("R diferent!\n");
    }

  }

  printf("Memory frees\n");

  free(B.vector);
  free(S.vector);
  free(Scomp.vector);
  free(R.vector);
  free(Rcomp.vector);

  return 0;

}
