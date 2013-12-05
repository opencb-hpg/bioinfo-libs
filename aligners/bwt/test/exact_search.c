#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#include "../bwt_commons.h"
#include "../search/search.h"
#include "../search/io.h"
#include "../csalib/csa.h"

#include "commons/commons.h"

int main(int argc, char **argv) {

  bwt_index backward, forward;
  bwt_config_t bwt_config;

  char *Worig;
  ref_vector W;

  results_list rl_prev, rl_next, rl_prev_i, rl_next_i, rl_final;
  uintmax_t read_index=0;

  uintmax_t RESULTS, FRAGSIZE;
  exome ex;

  FILE *queries_file, *output_file;

  check_syntax(argc, 6, "exact_search search_file input_dir output_file results_buffer frag_size");

  timevars();

	char *nucleotides = (char *) malloc(128 * sizeof(char));
  read_config(nucleotides, &(bwt_config.duplicate_strand), argv[2]);
  bwt_init_replace_table(&bwt_config, nucleotides);

  queries_file = fopen(argv[1], "r");
  check_file_open(queries_file, argv[1]);
  output_file = fopen(argv[3], "w");
  check_file_open(output_file, argv[3]);

  RESULTS  = atoi(argv[4]);
  FRAGSIZE = atoi(argv[5]);

  if (FRAGSIZE <= 0) {
    fprintf(stderr, "Fragsize must be greater than 0\n");
    exit(1);
  }

  tic("Loading FM-Index");
  load_bwt_index(NULL, &backward, argv[2], 1, true, bwt_config);
  load_bwt_index(NULL, &forward, argv[2], 0, true, bwt_config);
  toc();

  tic("Preparing search space");

  load_exome_file(&ex, argv[2]);
  new_results_list(&rl_prev, RESULTS); new_results_list(&rl_prev_i, RESULTS);
  new_results_list(&rl_next, RESULTS); new_results_list(&rl_next_i, RESULTS);
  new_results_list(&rl_final, RESULTS);

  Worig = (char *) malloc( MAXLINE * sizeof(char) );
  check_malloc(Worig, "main");
  W.vector = (uint8_t *) malloc( MAXLINE * sizeof(uint8_t) );
  check_malloc(W.vector, "main");

  toc();

  intmax_t *k = (intmax_t*)malloc(RESULTS * sizeof(intmax_t));
  intmax_t *l = (intmax_t*)malloc(RESULTS * sizeof(intmax_t));

  uintmax_t nW_aux;

  tic("Sequence mapping");

  while(nextFASTAToken(queries_file, Worig, W.vector, &nW_aux, &bwt_config)) {

    if (!W.vector) {
      printf ("Error de lectura de la cadena a buscar\n");
      return -1;
    }

    W.n = nW_aux;

    rl_prev.read_index = read_index; rl_prev_i.read_index = read_index;
    rl_next.read_index = read_index; rl_next_i.read_index = read_index;
    rl_final.read_index = read_index;

    BWSearchCPU(W.vector, W.n, &backward, &forward, &rl_prev, &rl_next, &rl_prev_i, &rl_next_i, &rl_final, NULL, FRAGSIZE, 1, bwt_config.nA);
    write_results(&rl_final, k, l, &ex, &backward, &forward, Worig, nW_aux, 2, output_file, &bwt_config);

    rl_final.num_results=0;

    read_index++;

  }

  toc();

  printf("Memory frees\n");

  bwt_free_replace_table(&bwt_config);

  free(W.vector);
  free(Worig);

  free(k);
  free(l);

  free_bwt_index(NULL, &backward, true);
  free_bwt_index(NULL, &forward, true);

  free(rl_prev.list);
  free(rl_next.list);
  free(rl_prev_i.list);
  free(rl_next_i.list);
  free(rl_final.list);

  fclose(queries_file);
  fclose(output_file);

  return 0;
}
