#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sys/time.h>

#include "commons/commons.h"
#include "BW_search.h"
#include "BW_io.h"

int found;
unsigned int enW;
char search[MAXLINE];

unsigned int read_index=0;

int main(int argc, char **argv) {

  byte_vector Worig, W;
  vector C, C1;
  comp_vector S, Si;
  comp_matrix O, Oi;

  size_t *vec_k=NULL, *vec_l=NULL, *vec_ki=NULL, *vec_li=NULL;

  results_list r_list;
  int start, end;

  int flag_simple;

  unsigned int RESULTS;
  exome ex;

  FILE *queries_file, *output_file;

  if (argc!=7) {
    printf("Syntax:\n\t%s search_file input_dir output_file results_buffer flag_simple nucleotides\n", argv[0]);
    return 1;
  }

  timevars();
  initReplaceTable(argv[6]);

  queries_file = fopen(argv[1], "r");
  checkFileOpen(queries_file, argv[1]);
  output_file = fopen(argv[3], "w");
  checkFileOpen(output_file, argv[3]);

  RESULTS =  atoi(argv[4]);
  flag_simple = atoi(argv[5]);

  tic("Loading data structures from input directory");

  readUIntVector(&C, argv[2], "C");
  printUIntVector(C.vector, C.n);
  readUIntVector(&C1, argv[2], "C1");
  printUIntVector(C1.vector, C1.n);
  readCompMatrix(&O, argv[2], "O");
  printCompMatrix(O);
  readCompMatrix(&Oi, argv[2], "Oi");
  printCompMatrix(Oi);
  readUIntCompVector(&S, argv[2], "Scomp");
  printUIntVector(S.vector, S.n);
  readUIntCompVector(&Si, argv[2], "Scompi");
  printUIntVector(Si.vector, Si.n);

  toc();

  tic("Preparing search space");

  load_exome_file(&ex, argv[2]);
  new_results_list(&r_list, RESULTS);

  W.vector = (char *) malloc( MAXLINE * sizeof(char) );
  checkMalloc(W.vector, "main");
  Worig.vector = (char *) malloc( MAXLINE * sizeof(char) );
  checkMalloc(Worig.vector, "main");

  toc();

  unsigned int nW_aux;
  result res;

  while(nextFASTAToken(queries_file, Worig.vector, W.vector, &nW_aux)) {

    if (!W.vector) {
      printf ("Error de lectura de la cadena a buscar\n");
      return -1;
    }

    Worig.n = nW_aux;
    W.n = nW_aux;

    printUIntVector(Worig.vector, Worig.n);
    printUIntVector(W.vector, W.n);

    start = 0;
    end   = W.n-1;

    //W.n=end-start+1;

    r_list.num_results = 0;
    r_list.read_index  = read_index;

    //tic("Sequence mapping");

    if (!flag_simple) {


      vec_k  = (size_t *) malloc(W.n * sizeof(size_t));
      vec_l  = (size_t *) malloc(W.n * sizeof(size_t));
      vec_ki = (size_t *) malloc(W.n * sizeof(size_t));
      vec_li = (size_t *) malloc(W.n * sizeof(size_t));

      BWExactSearchVectorBackward(W.vector, start, end, 0, O.siz-2, vec_k,  vec_l,  &C, &C1, &O);
      BWExactSearchVectorForward(W.vector, start, end, 0, O.siz-2, vec_ki, vec_li, &C, &C1, &Oi);

      /* for (size_t i=0; i< W.n; i++) */
      /* 	printf("%lu - %lu | ", vec_k[i], vec_l[i]); */
      /* printf("\n\n"); */

      /* for (size_t i=0; i< W.n; i++) */
      /*   printf("%lu - %lu | ", vec_ki[i], vec_li[i]); */
      /* printf("\n\n"); */

      //printf("%u - %u\n", vec_k[nW-1], vec_l[nW-1]);
      //printf("%u - %u\n", vec_ki[0], vec_li[0]);

    }

    tic("Sequence mapping");

    if (flag_simple==0) 
{
      BWSearch1GPUHelper(
		W.vector,
		start,
		end,
		vec_k,
		vec_l,
		vec_ki,
		vec_li,
		&C,
		&C1,
		&O,
		&Oi,
		&r_list
		);

    } else if (flag_simple==1) {

      init_result(&res, 1);
      bound_result(&res, start, end);
      change_result(&res, 0, O.siz-2, end);

      BWSimpleSearch1Backward(
			     W.vector,
			     &C,
			     &C1,
			     &O,
			     &res,
			     &r_list
			     );

    } else {

      init_result(&res, 0);
      bound_result(&res, start, end);
      change_result(&res, 0, O.siz-2, start);

      BWSimpleSearch1Forward(
  			     W.vector,
  			     &C,
  			     &C1,
  			     &Oi,
  			     &res,
  			     &r_list
  			     );

    }

    toc();

    write_results(&r_list, &ex, &S, &Si, &C, &O, &Oi, Worig.vector, Worig.n, 0, output_file);
    printf("\nFound %u results.\n", r_list.num_results);

    read_index++;

    if (!flag_simple) {
      free(vec_k);
      free(vec_l);
      free(vec_ki);
      free(vec_li);
    }

    /*
      for (unsigned int i=0;i<p;i++)
      printf("%u - %u\n", elementos[i][0], elementos[i][1]);
    */

  }

  printf("Memory frees\n");

  free(W.vector);
  free(Worig.vector);

  free(C.vector);
  free(C1.vector);

  freeCompMatrix(&O);
  freeCompMatrix(&Oi);

  free(S.vector);
  free(Si.vector);

  free(r_list.list);

  fclose(queries_file);
  fclose(output_file);

  return 0;

}
