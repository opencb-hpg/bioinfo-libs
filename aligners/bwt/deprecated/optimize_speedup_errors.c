#include <cuda_runtime_api.h>

#include "BW_gpu.cuh"

//#include "cuPrintf.cuh"

#define MAX_READ_GPU 256000
#define TAM_BLOQUE_GPU 32

float t_cpu=0, t_gpu=0;
int MAX_BUS_GPU=4096;

exome ex;

int main(int argc, char **argv) {

  char *h_Worig, *h_We, *d_We;
  unsigned int *h_nWe;

  vector h_C, d_C, h_C1, d_C1;
  comp_matrix h_O, d_O, h_Oi, d_Oi;
  comp_vector h_S, h_Si, h_R, h_Ri;

  blocked_results_lists rl_prev_cpu, rl_next_cpu, rl_prev_i_cpu, rl_next_i_cpu, rl_final_cpu;
  blocked_results_lists rl_prev_gpu, rl_next_gpu, rl_prev_i_gpu, rl_next_i_gpu, rl_final_gpu;
  results_list rl_prev, rl_next, rl_prev_i, rl_next_i, rl_final;

  unsigned int read_index=0;

  FILE *queries_file, *output_file;
  int RESULTS, FRAGSIZE;

  cudaError_t error;

  if (argc!=7) {
    printf("Syntax:\n\t%s search_file input_dir output_file results_buffer frag_size nucleotides\n", argv[0]);
    return 1;
  }

  cudaSetDevice(0);

  timevars();
  initReplaceTable(argv[6]);

  queries_file = fopen(argv[1], "r");
  checkFileOpen(queries_file, argv[1]);
  output_file = fopen(argv[3], "w");
  checkFileOpen(output_file, argv[3]);

  RESULTS  = atoi(argv[4]);
  FRAGSIZE = atoi(argv[5]);

  if (FRAGSIZE <= 0) {
    fprintf(stderr, "Fragsize must be greater than 0\n");
    exit(1);
  }

  tic("Loading FM-Index");

  readUIntVector(&h_C, argv[2], "C");
  readUIntVector(&h_C1, argv[2], "C1");
  readCompMatrixGPU(&h_O, argv[2], "O");
  readCompMatrixGPU(&h_Oi, argv[2], "Oi");

  readUIntCompVector(&h_S,  argv[2], "Scomp");
  readUIntCompVector(&h_Si, argv[2], "Scompi");
  readUIntCompVector(&h_R,  argv[2], "Rcomp");
  readUIntCompVector(&h_Ri, argv[2], "Rcompi");

  copyVectorGPU(&d_C,   &h_C,   sizeof(unsigned int));
  copyVectorGPU(&d_C1,  &h_C1,  sizeof(unsigned int));
  copyCompMatrixGPU(&d_O, &h_O);
  copyCompMatrixGPU(&d_Oi, &h_Oi);

  toc();

  tic("Preparing search space");

  load_exome_file(&ex, argv[2]);

  declare_blocked_results_list_cpu(&rl_prev_cpu,   RESULTS, MAX_BUS_GPU);
  declare_blocked_results_list_cpu(&rl_prev_i_cpu, RESULTS, MAX_BUS_GPU);
  declare_blocked_results_list_cpu(&rl_next_cpu,   RESULTS, MAX_BUS_GPU);
  declare_blocked_results_list_cpu(&rl_next_i_cpu, RESULTS, MAX_BUS_GPU);
  declare_blocked_results_list_cpu(&rl_final_cpu,  RESULTS, MAX_BUS_GPU);

  declare_blocked_results_list_gpu(&rl_prev_gpu,   RESULTS, MAX_BUS_GPU);
  declare_blocked_results_list_gpu(&rl_prev_i_gpu, RESULTS, MAX_BUS_GPU);
  declare_blocked_results_list_gpu(&rl_next_gpu,   RESULTS, MAX_BUS_GPU);
  declare_blocked_results_list_gpu(&rl_next_i_gpu, RESULTS, MAX_BUS_GPU);
  declare_blocked_results_list_gpu(&rl_final_gpu,  RESULTS, MAX_BUS_GPU);

  toc();

  h_Worig = (char*)malloc(MAX_BUS_GPU * MAXLINE     * sizeof(char));

  cudaMallocHost((void**) &h_We, MAX_BUS_GPU * MAXLINE * sizeof(char));
  cudaMallocHost((void**) &h_nWe, MAX_BUS_GPU * sizeof(unsigned int));

  cudaMalloc((void**) &d_We,  MAX_READ_GPU * MAXLINE * sizeof(char));
  manageCudaError();

  int TAM_BUS_GPU=0, NUM_BLOQUES_GPU=0;

  tic("Read mappings from disk");

  while(nextFASTAToken(queries_file, h_Worig + TAM_BUS_GPU * MAXLINE, h_We + TAM_BUS_GPU * MAXLINE, h_nWe + TAM_BUS_GPU)) {

    TAM_BUS_GPU++;
    if (TAM_BUS_GPU == MAX_BUS_GPU) break;

  }

  toc();

  NUM_BLOQUES_GPU = (TAM_BUS_GPU / TAM_BLOQUE_GPU);
  if (TAM_BUS_GPU % TAM_BLOQUE_GPU) NUM_BLOQUES_GPU++;

  cudaThreadSynchronize();
  tic("CPU -> GPU");
  cudaMemcpy(d_We, h_We, TAM_BUS_GPU * MAXLINE * sizeof(char), cudaMemcpyHostToDevice);
  manageCudaError();
  cudaThreadSynchronize();
  toc();

  cudaThreadSynchronize();
  tic("GPU Kernel");

  //cudaPrintfInit();
  printf("%d %d\n", NUM_BLOQUES_GPU, TAM_BLOQUE_GPU);

  BWSearchGPU(NUM_BLOQUES_GPU, TAM_BLOQUE_GPU, d_We, h_We, h_nWe[0], &d_C, &h_C, &d_C1, &h_C1, &d_O, &h_O, &d_Oi, &h_Oi, &h_S, &h_R, &h_Si, &h_Ri, &rl_prev_cpu, &rl_next_cpu, &rl_prev_i_cpu, &rl_next_i_cpu, &rl_final_cpu, &rl_prev_gpu, &rl_next_gpu, &rl_prev_i_gpu, &rl_next_i_gpu, &rl_final_gpu, FRAGSIZE, RESULTS);

  /* cudaPrintfDisplay(stdout, true); */
  /* cudaPrintfEnd(); */

  cudaThreadSynchronize();
  toc();

  cudaThreadSynchronize();
  tic("GPU -> CPU");
  copy_blocked_results_list_cpu(&rl_final_cpu, &rl_final_gpu, RESULTS, TAM_BUS_GPU);
  cudaThreadSynchronize();
  toc();

  //printf("%d\n", rl_final_cpu.num_results[0]);
  write_blocked_results(&rl_final_cpu, &ex, &h_S, &h_Si, &h_C, &h_O, &h_Oi, h_Worig, h_nWe[0], 1, output_file, RESULTS, TAM_BUS_GPU, 0);

  new_results_list(&rl_prev, RESULTS); new_results_list(&rl_prev_i, RESULTS);
  new_results_list(&rl_next, RESULTS); new_results_list(&rl_next_i, RESULTS);
  new_results_list(&rl_final, RESULTS);

  tic("CPU kernel");
  for(int i=0; i<MAX_BUS_GPU; i++) {

    rl_prev.num_results = 0; rl_prev_i.num_results = 0;
    rl_next.num_results = 0; rl_next_i.num_results = 0;
    rl_final.num_results = 0;

    rl_prev.read_index = read_index; rl_prev_i.read_index = read_index;
    rl_next.read_index = read_index; rl_next_i.read_index = read_index;
    rl_final.read_index = read_index;
  
    BWSearchCPU(h_We + i*MAXLINE, h_nWe[0], &h_C, &h_C1, &h_O, &h_Oi, &h_S, &h_R, &h_Si, &h_Ri, &rl_prev, &rl_next, &rl_prev_i, &rl_next_i, &rl_final, FRAGSIZE, 1);

    read_index++;

  }
  toc();

  /*
  for (int i=0; i<TAM_BUS_GPU; i++) {

    for (int j=0; j<h_nWe[i]; j++) {

      if (h_k[i*MAXLINE + j] != h_k2[i*MAXLINE + j]) {
	printf("Diferente %d %d\n", i, j);
	goto salir;
      }

    }

  }

  salir:
  */

  /*
  for (int i=0; i<h_nWe[0]; i++) {
    printf("%u ", h_k[i]);
  }
  printf("\n");

  printf("\n");

  for (int i=0; i<h_nWe[0]; i++) {
    printf("%u ", h_k2[i]);
  }
  printf("\n");
  */

  cudaFreeHost(h_We);
  cudaFreeHost(h_nWe);
  cudaFree(d_We);

  free(h_C.vector);
  cudaFree(d_C.vector);
  free(h_C1.vector);
  cudaFree(d_C1.vector);

  freeCompMatrixGPUHost(&h_O);
  freeCompMatrixGPUDevice(&d_O);
  freeCompMatrixGPUHost(&h_Oi);
  freeCompMatrixGPUDevice(&d_Oi);

  fclose(queries_file);
  fclose(output_file);

  return 0;

}
