#include <pthread.h>
#include <cuda_runtime_api.h>

#include "../gpu/gpu.cuh"
#include "../search/search.h"
#include "../search/io.h"

#define MAX_READ_GPU 256000
#define TAM_BLOQUE_GPU 32
#define MAX_SEARCH 200
#define NUM_CARDS 2
#define RESULTS 2000

char *h_Worig, *h_Worig2, *h_Worig3;
uint8_t *h_We, *h_We2, *h_We3;
uint64_t *h_nWe, *h_nWe2, *h_nWe3;

char *read_Worig, *gpu_Worig, *store_Worig, *swap_Worig;
uint8_t *read_We, *gpu_We, *store_We, *swap_We;
uint64_t *read_nWe, *gpu_nWe, *store_nWe, *swap_nWe;

bwt_index backward, forward, backward_rev, forward_rev;

comp_matrix h_O, h_rO, h_Oi, h_rOi;

vector h_C, h_rC, h_C1, h_rC1;
comp_vector S, Si;

intmax_t *h_k, *h_l, *h_k2, *h_l2;
intmax_t *h_ki, *h_li, *h_ki2, *h_li2;

intmax_t *gpu_h_k, *gpu_h_l;
intmax_t *store_h_k, *store_h_l;
intmax_t *swap_h_k, *swap_h_l;

intmax_t *gpu_h_ki, *gpu_h_li;
intmax_t *store_h_ki, *store_h_li;
intmax_t *swap_h_ki, *swap_h_li;

intmax_t *k, *l;

cudaError_t error;

struct timeval t1, t2, t1_write, t2_write;
float t_gpu=0, t_total=0, t_write=0;

pthread_mutex_t gpu_time_thread, print_results;
pthread_cond_t start_time, start_time_write, stop_time, stop_time_write;

int stop, start[NUM_CARDS], end;
int start_write=0, stop_write[NUM_CARDS], end_write;
uintmax_t tam_read_gpu=0, tam_read_gpu2=0;
uintmax_t tam_write=0, num_write=0;

FILE *output_file, *notfound_file;

exome ex;

int num_errors;

void *writeResults(void *threadid) {

	results_list r_list;

	char search[201];
	char plusminus[] = "-+";
	uintmax_t index, key;

	intmax_t *h_kaux, *h_laux;
	intmax_t *h_kiaux, *h_liaux;

	uintmax_t w=0;

	uintmax_t contador=0;
	uintmax_t descartadas=0;
	uintmax_t total=0;

	new_results_list(&r_list, RESULTS);

	while(true) {

		pthread_mutex_lock(&print_results);

		//printf("W -> Para\n");

		while(start_write!=NUM_CARDS) {
			pthread_cond_wait(&start_time_write, &print_results);
		}

		start_write=0;

		//printf("W -> Sigue\n");

		pthread_mutex_unlock(&print_results);

		if (end_write) {

			free(r_list.list);

			//printf("W -> Saliendo\n");
			printf("%lu founds of %lu -> %.2f, discarded %lu -> %.2f\n", contador, total, contador * 100.0 / total, descartadas, contador * 100.0 / (total-descartadas));
			pthread_exit(NULL);

		}

		//printf("W -> Swap buffers\n");

		tam_write = tam_read_gpu2;

		num_write = MAX_READ_GPU*w;

		swap_Worig = gpu_Worig;
		gpu_Worig = store_Worig;
		store_Worig = swap_Worig;

		swap_We = gpu_We;
		gpu_We = store_We;
		store_We = swap_We;

		swap_h_k = gpu_h_k;
		gpu_h_k = store_h_k;
		store_h_k = swap_h_k;

		swap_h_l = gpu_h_l;
		gpu_h_l = store_h_l;
		store_h_l = swap_h_l;

		swap_h_ki = gpu_h_ki;
		gpu_h_ki = store_h_ki;
		store_h_ki = swap_h_ki;

		swap_h_li = gpu_h_li;
		gpu_h_li = store_h_li;
		store_h_li = swap_h_li;

		swap_nWe = gpu_nWe;
		gpu_nWe = store_nWe;
		store_nWe = swap_nWe;

		pthread_mutex_lock(&print_results);

		//printf("W -> Ordena seguir a los hilos de búsqueda GPU\n");
		for(int i=0; i<NUM_CARDS; i++) {
			stop_write[i] = 1;
		}
		pthread_cond_broadcast(&stop_time_write);

		pthread_mutex_unlock(&print_results);

		//printf("W -> Empieza a escribir %d resultados\n", tam_write);

		if (!num_errors) {

			gettimeofday(&t1_write, NULL);

			for (uintmax_t i=0; i < tam_write; i++) {

				bool found = false;

				search[0] = '\0';
				strncat(search, store_Worig + i * MAXLINE, store_nWe[i]);

				int n_count = 0;
				for (uint64_t nn=0; nn<store_nWe[i]; nn++) {
					if (search[nn] == 'N') n_count++;
					if (n_count>1) break; 
				}

				if (n_count>num_errors) {

					descartadas++;

				} else {

					for(int type=1; type>=0; type--) {

						h_kaux = store_h_k + MAX_READ_GPU * type;
						h_laux = store_h_l + MAX_READ_GPU * type;

						//fprintf(output_file, "%u - %u\n", h_kaux[i], h_laux[i]);
						//printf("%u - %u\n", h_kaux[i], h_laux[i]);

						for (intmax_t j=h_kaux[i]; j<=h_laux[i]; j++) {

							key = get_SA(j, &backward);

							index = binsearch(ex.offset, ex.size, key);

							if(key + store_nWe[i] <= ex.offset[index]) {
								found = true;
								fprintf(output_file, "read_%ju %c %s %ju 0 %s %s\n", num_write + i, plusminus[type], ex.chromosome + (index-1)*IDMAX, (uintmax_t) ex.start[index-1] + (key - ex.offset[index-1]), search, search);
							}

						}

					} //type

				}

				if (!found)
					fprintf(notfound_file, ">read_%ju\n%s\n", num_write + i, search);
				else
					contador++;

			}

			w++;

			gettimeofday(&t2_write, NULL);
			t_gpu = (t2_write.tv_sec-t1_write.tv_sec)*1e6+(t2_write.tv_usec-t1_write.tv_usec);
			t_write += t_gpu;

			//printf("W -> Termina de escribir hasta el resultado %d\n", w*MAX_READ_GPU);

		} else {

			gettimeofday(&t1_write, NULL);

			for (uintmax_t i=0; i < tam_write; i++) {

				bool found = false, found2 = false;

				int n_count = 0;
				for (uint64_t nn=0; nn<store_nWe[i]; nn++) {
					if (store_Worig[i*MAXLINE + nn] == 'N') n_count++;
					if (n_count>1) break; 
				}	

				if (n_count>num_errors) {

					descartadas++;

				} else {

					for(int type=1; type>=0; type--) {

						h_kaux = store_h_k + MAX_READ_GPU * MAXLINE * type;
						h_laux = store_h_l + MAX_READ_GPU * MAXLINE * type;

						h_kiaux = store_h_ki + MAX_READ_GPU * MAXLINE * type;
						h_liaux = store_h_li + MAX_READ_GPU * MAXLINE * type;

						//fprintf(output_file, "%u - %u\n", h_kaux[i*MAXLINE + store_nWe[i]-1], h_laux[i*MAXLINE + store_nWe[i]-1]);
						//fprintf(output_file, "%u - %u\n", h_kiaux[i*MAXLINE], h_liaux[i*MAXLINE]);

						r_list.num_results = 0;
						r_list.read_index = num_write + i;

						if (type) {

							BWSearch1GPUHelper(
									store_We + i * MAXLINE,
									0,
									store_nWe[i]-1,
									h_kaux  + i*MAXLINE,
									h_laux  + i*MAXLINE,
									h_kiaux + i*MAXLINE,
									h_liaux + i*MAXLINE,
									&backward,
									&forward,
									&r_list
									);

							found2 = false;
							found2 = write_results(&r_list, k, l, &ex, &backward, &forward, store_Worig + i*MAXLINE, store_nWe[i], type, output_file);
							found = found || found2;

						} else {

							BWSearch1GPUHelper(
									store_We + i * MAXLINE,
									0,
									store_nWe[i]-1,
									h_kiaux + i*MAXLINE,
									h_liaux + i*MAXLINE,
									h_kaux  + i*MAXLINE,
									h_laux  + i*MAXLINE,
									&forward_rev,
									&backward_rev,
									&r_list
									);

							found2 = false;
							found2 = write_results(&r_list, k, l, &ex, &backward, &forward, store_Worig + i*MAXLINE, store_nWe[i], type, output_file);
							found = found || found2;

						}

					} //type

				}

				if (!found) {
					search[0] = '\0';
					strncat(search, store_Worig + i * MAXLINE, store_nWe[i]);
					fprintf(notfound_file, ">read_%ju\n%s\n", num_write + i, search);
				} else
					contador++;

				total++;

			}

			w++;

			gettimeofday(&t2_write, NULL);
			t_gpu = (t2_write.tv_sec-t1_write.tv_sec)*1e6+(t2_write.tv_usec-t1_write.tv_usec);
			t_write += t_gpu;
			//printf("%lu founds of %u -> %.2f, discarded %lu -> %.2f\n", contador, num_write + tam_write, contador * 100.0 / (num_write+tam_write), descartadas, contador * 100.0 / (num_write+tam_write-descartadas));
			//printf("W -> Termina de escribir hasta el resultado %d\n", w*MAX_READ_GPU);

			}

		}

	}

	void *gpuSearch(void *threadid) {

		int tid;
		tid = (long)threadid;

		vector d_C, d_C1, d_rC, d_rC1;
		comp_matrix d_O, d_rO;

		uint8_t *d_We;
		uint64_t *d_nWe;
		intmax_t *d_k, *d_l;

		cudaSetDevice(tid);
		manageCudaError();

		//printf("%d -> Copiar vector C\n", tid);
		copy_vector_gpu(&d_C,   &h_C);
		copy_vector_gpu(&d_C1,  &h_C1);
		copy_vector_gpu(&d_rC,  &h_rC);
		copy_vector_gpu(&d_rC1, &h_rC1);
		//printf("%d -> Termina de copiar vector C\n", tid);

		//printf("%d -> Copiar vector O\n", tid);
		if (tid==1) {

			if(!num_errors) {
				copy_comp_matrix_gpu(&d_O, &h_O);
			} else {
				copy_comp_matrix_gpu(&d_O, &h_O);
				reverse_strand_gpu_O(&d_rO, &d_O);
			}

		} else {

			if(!num_errors) {
				copy_comp_matrix_gpu(&d_O, &h_rO);
			} else {
				copy_comp_matrix_gpu(&d_O, &h_Oi);
				reverse_strand_gpu_O(&d_rO, &d_O);
			}

		}

		//printf("%d -> Termina vector O\n", tid);

		//printf("%d -> Vector W\n", tid);
		cudaMalloc((void**) &d_We,  MAX_READ_GPU * MAXLINE * sizeof(uint8_t));
		manageCudaError();
		cudaMalloc((void**) &d_nWe, MAX_READ_GPU * sizeof(uint64_t));
		manageCudaError();
		//printf("%d -> Termina vector O\n", tid);

		//printf("%d -> Vector k l\n", tid);
		cudaMalloc((void**) &d_k,   MAX_READ_GPU * MAXLINE * sizeof(intmax_t));
		manageCudaError();
		cudaMalloc((void**) &d_l,   MAX_READ_GPU * MAXLINE * sizeof(intmax_t));
		manageCudaError();
		//printf("%d -> Termina k l\n", tid);

		//printf("%d -> Antes del bucle\n", tid);

		uintmax_t tam_bus_gpu, num_bloques_gpu;

		while(1) {

			pthread_mutex_lock(&gpu_time_thread);
			//printf("%d -> Espera lectura\n", tid);

			while(!start[tid]) {
				pthread_cond_wait(&start_time, &gpu_time_thread);
			}

			start[tid]=0;

			//printf("%d -> Sigue\n", tid);

			pthread_mutex_unlock(&gpu_time_thread);

			if (end) {

				//printf("%d -> Saliendo\n", tid);

				pthread_mutex_lock(&print_results);

				//printf("%d -> Ordena terminar al hilo de escritura\n", tid);

				start_write++;
				end_write=1;

				pthread_cond_signal(&start_time_write);

				pthread_mutex_unlock(&print_results);

				cudaFree(d_C.vector);
				cudaFree(d_C1.vector);
				cudaFree(d_rC.vector);
				cudaFree(d_rC1.vector);

				free_comp_matrix_gpu_device(&d_rO, &d_O);

				cudaFree(d_We);
				cudaFree(d_nWe);

				cudaFree(d_k);
				cudaFree(d_l);

				pthread_exit(NULL);

			}

			//printf("%d -> Pone valores\n", tid);

			uintmax_t aux_mod = tam_read_gpu % TAM_BLOQUE_GPU;

			if (aux_mod)
				tam_bus_gpu = tam_read_gpu + (TAM_BLOQUE_GPU - aux_mod);
			else
				tam_bus_gpu = tam_read_gpu;

			num_bloques_gpu = tam_bus_gpu / TAM_BLOQUE_GPU;

			if (tid) tam_read_gpu2 = tam_read_gpu;

			//printf("%d -> Copia en la gpu\n", tid);

			cudaMemcpy(d_We, read_We, tam_bus_gpu * MAXLINE * sizeof(uint8_t), cudaMemcpyHostToDevice);
			manageCudaError();
			cudaMemcpy(d_nWe, read_nWe, tam_bus_gpu * sizeof(uint64_t), cudaMemcpyHostToDevice);
			manageCudaError();

			//printf("%d -> Ordena seguir al hilo de lectura\n", tid);

			pthread_mutex_lock(&gpu_time_thread);

			stop++;
			if (stop==NUM_CARDS)
				pthread_cond_signal(&stop_time);

			pthread_mutex_unlock(&gpu_time_thread);

			if(!num_errors) {

				if (tid) {
					//printf("%d -> Lanzo kernel 0 - %lu\n", tid, d_O.siz-2);
					BWExactSearchBackwardGPUWrapper(num_bloques_gpu, TAM_BLOQUE_GPU, d_We, d_nWe, d_k, d_l, 0, d_O.siz-2, &d_C, &d_C1, &d_O);
					manageCudaError();
				} else {
					//printf("%d -> Lanzo kernel 0 - %lu\n", tid, d_O.siz-2);
					BWExactSearchForwardGPUWrapper(num_bloques_gpu, TAM_BLOQUE_GPU, d_We, d_nWe, d_k, d_l, 0, d_O.siz-2, &d_rC, &d_rC1, &d_O);
					manageCudaError();
				}

				//printf("%d -> Copia %d resultados a la CPU\n", tid, tam_read_gpu2);
				cudaMemcpy(gpu_h_k + MAX_READ_GPU * tid, d_k, sizeof(intmax_t) * tam_bus_gpu, cudaMemcpyDeviceToHost);
				manageCudaError();
				cudaMemcpy(gpu_h_l + MAX_READ_GPU * tid, d_l, sizeof(intmax_t) * tam_bus_gpu, cudaMemcpyDeviceToHost);
				manageCudaError();
				//printf("%d -> Termina de copiar %d resultados a la CPU\n", tid, tam_read_gpu2);

			} else {

				if (tid) {

					//printf("%d -> Lanzo kernel 0 - %lu\n", tid, d_O.siz-2);
					BWExactSearchBackwardVectorGPUWrapper(num_bloques_gpu, TAM_BLOQUE_GPU, d_We, d_nWe, d_k, d_l, 0, d_O.siz-2, &d_C, &d_C1, &d_O);
					manageCudaError();

					cudaMemcpy(gpu_h_k + MAX_READ_GPU * MAXLINE, d_k, sizeof(intmax_t) * tam_bus_gpu * MAXLINE, cudaMemcpyDeviceToHost);
					manageCudaError();
					cudaMemcpy(gpu_h_l + MAX_READ_GPU * MAXLINE, d_l, sizeof(intmax_t) * tam_bus_gpu * MAXLINE, cudaMemcpyDeviceToHost);
					manageCudaError();

					//printf("%d -> Lanzo kernel 0 - %lu\n", tid, d_rO.siz-2);
					BWExactSearchForwardVectorGPUWrapper(num_bloques_gpu, TAM_BLOQUE_GPU, d_We, d_nWe, d_k, d_l, 0, d_rO.siz-2, &d_rC, &d_rC1, &d_rO);
					manageCudaError();

					cudaMemcpy(gpu_h_k,                          d_k, sizeof(intmax_t) * tam_bus_gpu * MAXLINE, cudaMemcpyDeviceToHost);
					manageCudaError();
					cudaMemcpy(gpu_h_l,                          d_l, sizeof(intmax_t) * tam_bus_gpu * MAXLINE, cudaMemcpyDeviceToHost);
					manageCudaError();

				} else {

					//printf("%d -> Lanzo kernel 0 - %lu\n", tid, d_O.siz-2);
					BWExactSearchForwardVectorGPUWrapper(num_bloques_gpu, TAM_BLOQUE_GPU, d_We, d_nWe, d_k, d_l, 0, d_O.siz-2, &d_C, &d_C1, &d_O);
					manageCudaError();

					cudaMemcpy(gpu_h_ki + MAX_READ_GPU * MAXLINE, d_k, sizeof(intmax_t) * tam_bus_gpu * MAXLINE, cudaMemcpyDeviceToHost);
					manageCudaError();
					cudaMemcpy(gpu_h_li + MAX_READ_GPU * MAXLINE, d_l, sizeof(intmax_t) * tam_bus_gpu * MAXLINE, cudaMemcpyDeviceToHost);
					manageCudaError();

					//printf("%d -> Lanzo kernel 0 - %lu\n", tid, d_rO.siz-2);
					BWExactSearchBackwardVectorGPUWrapper(num_bloques_gpu, TAM_BLOQUE_GPU, d_We, d_nWe, d_k, d_l, 0, d_rO.siz-2, &d_rC, &d_rC1, &d_rO);
					manageCudaError();

					cudaMemcpy(gpu_h_ki,                          d_k, sizeof(intmax_t) * tam_bus_gpu * MAXLINE, cudaMemcpyDeviceToHost);
					manageCudaError();
					cudaMemcpy(gpu_h_li,                          d_l, sizeof(intmax_t) * tam_bus_gpu * MAXLINE, cudaMemcpyDeviceToHost);
					manageCudaError();

				}

			}

			pthread_mutex_lock(&print_results);

			//printf("%d -> Ordena continuar al hilo de escritura\n", tid);
			start_write++;

			pthread_cond_signal(&start_time_write);

			//printf("%d -> Espera escritura\n", tid);

			while(!stop_write[tid]) {
				pthread_cond_wait(&stop_time_write, &print_results);
			}

			stop_write[tid]=0;

			//printf("%d -> Sigue\n", tid);

			pthread_mutex_unlock(&print_results);

		}

	}

	int main(int argc, char **argv) {

		gettimeofday(&t1, NULL);

		FILE *queries_file;

		check_syntax(argc, 7, "search_gpu f_mappings d_transform f_output f_notfound num_errors nucleotides");

		num_errors = atoi(argv[5]);
		init_replace_table(argv[6]);

		read_vector(&h_C,  argv[2], "C");
		read_vector(&h_C1, argv[2], "C1");
		read_comp_matrix_gpu(&h_O,  argv[2], "O");
		read_comp_vector(&S, argv[2], "S");

		reverse_strand_C(&h_rC, &h_C, &h_rC1, &h_C1);
		reverse_strand_O(&h_rO, &h_O);

		if (num_errors) {
			read_comp_matrix_gpu(&h_Oi, argv[2], "Oi");
			read_comp_vector(&Si, argv[2], "Si");
			reverse_strand_O(&h_rOi, &h_Oi);
		}

		backward.C  = h_C;
		backward.C1 = h_C1;
		backward.O  = h_O;
		backward.S  = S;

		forward.C  = h_C;
		forward.C1 = h_C1;
		forward.O  = h_Oi;
		forward.S  = Si;

		backward_rev.C  = h_rC;
		backward_rev.C1 = h_rC1;
		backward_rev.O  = h_rO;
		backward_rev.S  = S;

		forward_rev.C  = h_rC;
		forward_rev.C1 = h_rC1;
		forward_rev.O  = h_rOi;
		forward_rev.S  = Si;

		h_Worig  = (char*)malloc(MAX_READ_GPU * MAXLINE * sizeof(char));
		check_malloc(h_Worig,  "main");
		h_Worig2 = (char*)malloc(MAX_READ_GPU * MAXLINE * sizeof(char));
		check_malloc(h_Worig2, "main");
		h_Worig3 = (char*)malloc(MAX_READ_GPU * MAXLINE * sizeof(char));
		check_malloc(h_Worig3, "main");

		h_We  = (uint8_t*)malloc(MAX_READ_GPU * MAXLINE * sizeof(uint8_t));
		check_malloc(h_We,  "main");
		h_We2 = (uint8_t*)malloc(MAX_READ_GPU * MAXLINE * sizeof(uint8_t));
		check_malloc(h_We2, "main");
		h_We3 = (uint8_t*)malloc(MAX_READ_GPU * MAXLINE * sizeof(uint8_t));
		check_malloc(h_We3, "main");

		h_nWe  = (uint64_t*)malloc(MAX_READ_GPU * sizeof(uint64_t));
		check_malloc(h_nWe, "main");
		h_nWe2 = (uint64_t*)malloc(MAX_READ_GPU * sizeof(uint64_t));
		check_malloc(h_nWe2, "main");
		h_nWe3 = (uint64_t*)malloc(MAX_READ_GPU * sizeof(uint64_t));
		check_malloc(h_nWe3, "main");

		if (!num_errors) {

			h_k  = (intmax_t*)malloc(MAX_READ_GPU * NUM_CARDS * sizeof(intmax_t));
			check_malloc(h_k, "main");
			h_l  = (intmax_t*)malloc(MAX_READ_GPU * NUM_CARDS * sizeof(intmax_t));
			check_malloc(h_l, "main");
			h_k2 = (intmax_t*)malloc(MAX_READ_GPU * NUM_CARDS * sizeof(intmax_t));
			check_malloc(h_k2, "main");
			h_l2 = (intmax_t*)malloc(MAX_READ_GPU * NUM_CARDS * sizeof(intmax_t));
			check_malloc(h_l2, "main");

		} else {

			h_k  = (intmax_t*)malloc(MAX_READ_GPU * MAXLINE * NUM_CARDS * sizeof(intmax_t));
			check_malloc(h_k, "main");
			h_l  = (intmax_t*)malloc(MAX_READ_GPU * MAXLINE * NUM_CARDS * sizeof(intmax_t));
			check_malloc(h_l, "main");
			h_k2 = (intmax_t*)malloc(MAX_READ_GPU * MAXLINE * NUM_CARDS * sizeof(intmax_t));
			check_malloc(h_k2, "main");
			h_l2 = (intmax_t*)malloc(MAX_READ_GPU * MAXLINE * NUM_CARDS * sizeof(intmax_t));
			check_malloc(h_l2, "main");

			h_ki  = (intmax_t*)malloc(MAX_READ_GPU * MAXLINE * NUM_CARDS * sizeof(intmax_t));
			check_malloc(h_ki, "main");
			h_li  = (intmax_t*)malloc(MAX_READ_GPU * MAXLINE * NUM_CARDS * sizeof(intmax_t));
			check_malloc(h_li, "main");
			h_ki2 = (intmax_t*)malloc(MAX_READ_GPU * MAXLINE * NUM_CARDS * sizeof(intmax_t));
			check_malloc(h_ki2, "main");
			h_li2 = (intmax_t*)malloc(MAX_READ_GPU * MAXLINE * NUM_CARDS * sizeof(intmax_t));
			check_malloc(h_li2, "main");

		}

		k = (intmax_t*)malloc(RESULTS * sizeof(intmax_t));
		l = (intmax_t*)malloc(RESULTS * sizeof(intmax_t));

		read_Worig = h_Worig; gpu_Worig = h_Worig2; store_Worig = h_Worig3;
		read_We = h_We; gpu_We = h_We2; store_We = h_We3;
		read_nWe = h_nWe; gpu_nWe = h_nWe2; store_nWe = h_nWe3;

		gpu_h_k = h_k; store_h_k = h_k2;
		gpu_h_l = h_l; store_h_l = h_l2;

		gpu_h_ki = h_ki; store_h_ki = h_ki2;
		gpu_h_li = h_li; store_h_li = h_li2;

		queries_file = fopen(argv[1], "r");
		check_file_open(queries_file, argv[1]);
		output_file = fopen(argv[3], "w");
		check_file_open(output_file, argv[3]);
		notfound_file = fopen(argv[4], "w");
		check_file_open(notfound_file, argv[4]);

		pthread_t gpu_threads[NUM_CARDS];
		pthread_t write_thread;
		pthread_attr_t attr;
		int rc;

		pthread_mutex_init(&gpu_time_thread,  NULL);
		pthread_mutex_init(&print_results,  NULL);

		pthread_cond_init(&start_time, NULL);
		pthread_cond_init(&stop_time,  NULL);

		pthread_cond_init(&start_time_write, NULL);
		pthread_cond_init(&stop_time_write,  NULL);

		pthread_attr_init(&attr);
		pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

		stop=0; end=0;

		start_write=0; end_write=0;

		for(long i=0; i<NUM_CARDS; i++) {

			start[i]=0;
			stop_write[i]=0;

			rc = pthread_create(&gpu_threads[i], &attr, gpuSearch, (void *)i);
			if (rc){
				fprintf(stderr, "ERROR; return code from pthread_create() on gpuSearch thread %ld is %d\n", i, rc);
				return 1;
			}

		}

		rc = pthread_create(&write_thread, &attr, writeResults, (void *)NUM_CARDS);
		if (rc){                                                                  
			fprintf(stderr, "ERROR; return code from pthread_create() on writeResults thread %d is %d\n", NUM_CARDS, rc);
			return 1;
		}

		//printf("L -> Antes de cargar el fichero del exoma\n");
		load_exome_file(&ex, argv[2]);

		//printf("L -> Antes del bucle de lecturas\n");

		bool exit = true;

		while (exit) {

			//printf("L -> Empieza a leer\n");

			while (exit) {

				exit = nextFASTAToken(queries_file, read_Worig + tam_read_gpu * MAXLINE, read_We + tam_read_gpu * MAXLINE, read_nWe + tam_read_gpu);

				tam_read_gpu++;

				if (tam_read_gpu == MAX_READ_GPU) break;

			}

			//printf("L -> Termina de leer tam_read_gpu = %d, status(%d)\n", tam_read_gpu, exit);

			if (!exit) {
				tam_read_gpu--;
				if (tam_read_gpu == 0) break;
			}

			pthread_mutex_lock(&gpu_time_thread);

			//printf("L -> Digo a los GPU que copien\n");

			for(int i=0; i<NUM_CARDS; i++) {
				start[i] = 1;
			}
			pthread_cond_broadcast(&start_time);

			//printf("L -> Antes de parar\n");

			while(stop!=NUM_CARDS) {
				pthread_cond_wait(&stop_time, &gpu_time_thread);
			}

			stop = 0;

			pthread_mutex_unlock(&gpu_time_thread);

			//printf("L -> Swap buffers\n");

			swap_Worig = read_Worig;
			read_Worig = gpu_Worig;
			gpu_Worig = swap_Worig;

			swap_We = read_We;
			read_We = gpu_We;
			gpu_We = swap_We;

			swap_nWe = read_nWe;
			read_nWe = gpu_nWe;
			gpu_nWe = swap_nWe;

			tam_read_gpu = 0;

		}

		//Terminar los hilos
		pthread_mutex_lock(&gpu_time_thread);

		//printf("L -> Ordena salir a los hilos\n");
		for(int i=0; i<NUM_CARDS; i++) {
			start[i]=1;
		}

		end=1;

		pthread_cond_broadcast(&start_time);

		pthread_mutex_unlock(&gpu_time_thread);

		for(int i=0; i<NUM_CARDS; i++)
			pthread_join(gpu_threads[i], NULL);
		pthread_join(write_thread, NULL);

		pthread_mutex_destroy(&gpu_time_thread);
		pthread_mutex_destroy(&print_results);

		pthread_cond_destroy(&start_time);
		pthread_cond_destroy(&start_time_write);

		pthread_cond_destroy(&stop_time);
		pthread_cond_destroy(&stop_time_write);

		pthread_attr_destroy(&attr);

		free(h_C.vector);
		free(h_C1.vector);
		free(h_rC.vector);
		free(h_rC1.vector);

		free_comp_matrix_gpu_host(&h_rO, &h_O);
		free(S.vector);

		if (num_errors) {
			free_comp_matrix_gpu_host(&h_rOi, &h_Oi);
			free(Si.vector);
		}

		cudaFreeHost(h_We);

		free(h_nWe);
		free(h_nWe2);
		free(h_nWe3);

		free(h_k);
		free(h_l);
		free(h_k2);
		free(h_l2);

		free(k);
		free(l);

		fclose(queries_file);
		fflush(output_file);
		fclose(output_file);
		fflush(notfound_file);
		fclose(notfound_file);

		gettimeofday(&t2, NULL);
		t_total = (t2.tv_sec-t1.tv_sec)*1e6+(t2.tv_usec-t1.tv_usec);

		printf("Write: %f\n", t_write / 1000000);
		printf("Total: %f\n", t_total / 1000000);

		//printf("L -> Saliendo\n");

		return 0;

	}
