#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sys/time.h>

#include "commons/string_utils.h"
#include "BW_io.h"

int found;
unsigned int enW;
char search[MAXLINE];

unsigned int read_index=0;

int main(int argc, char **argv) {

  comp_matrix O;

  if (argc!=3) {
    printf("Syntax:\n\t%s input_dir nucleotides\n", argv[0]);
    return 1;
  }

  initReplaceTable(argv[2]);

  readCompMatrix(&O, argv[1], "O");

#if defined VECTOR_O_32BIT_COMPRESSION || VECTOR_O_64BIT_COMPRESSION

  for(unsigned int i=0; i<65; i++) {
    printf("%u\n", getOcompValue(0, i, &O));
  }

#else

  for(unsigned int i=0; i<65; i++) {
    printf("%u\n", O.desp[0][i]);
  }

#endif
  printf("Memory frees\n");

  freeCompMatrix(&O);

  return 0;

}
