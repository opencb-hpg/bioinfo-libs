#include "csafm.h"

void free_comp_matrix(comp_matrix *reverse, comp_matrix *strand) {

	for (SA_TYPE i=0; i<strand->n_desp; i++) {
		free(strand->desp[i]);
#if defined FM_COMP_32 || FM_COMP_64
		free(strand->count[i]);
#endif
	}

	free(strand->desp);
	if (reverse != NULL) free(reverse->desp);
#if defined FM_COMP_32 || FM_COMP_64
	free(strand->count);
	if (reverse != NULL) free(reverse->count);
#endif

}


void read_vector(vector *vector, const char *directory, const char *name) {

  size_t err=0;
  FILE *fp;

  char path[500]; //TODO: Change to dynamic allocation to avoid buffer overflow
  path[0]='\0';
  strcat(path, directory);
  strcat(path, "/");
  strcat(path, name);
  strcat(path, ".vec");

  fp  = fopen(path,  "rb+");
  check_file_open(fp, path);

  err = fread(&vector->n, sizeof(SA_TYPE),  1, fp);
  check_file_read(err, 1, path);

  vector->vector = (SA_TYPE *) malloc(vector->n * sizeof(SA_TYPE));
  check_malloc(vector->vector, path);

  err = fread(vector->vector, sizeof(SA_TYPE), vector->n, fp);
  check_file_read(err, vector->n, path);

  fclose(fp);
}

void read_comp_vector(comp_vector *vector, const char *directory, const char *name) {

  size_t err=0;
  FILE *fp;

	char path[500];
  path[0]='\0';
  strcat(path, directory);
  strcat(path, "/");
  strcat(path, name);
  strcat(path, ".vec");

	fp  = fopen(path,  "rb+");
  check_file_open(fp, path);

	err = fread(&vector->siz, sizeof(SA_TYPE),  1, fp);
  check_file_read(err, 1, path);

	err = fread(&vector->n, sizeof(SA_TYPE),  1, fp);
  check_file_read(err, 1, path);

  err = fread(&vector->ratio, sizeof(SA_TYPE),  1, fp);
  check_file_read(err, 1, path);

  vector->vector = (SA_TYPE *) malloc(vector->n * sizeof(SA_TYPE));
  check_malloc(vector->vector, path);

  err = fread(vector->vector, sizeof(SA_TYPE), vector->n, fp);
  check_file_read(err, vector->n, path);

  fclose(fp);
}

void read_comp_matrix(comp_matrix *matrix, const char *directory, const char *name) {

	size_t err=0;
  FILE *fp;

  char path[500];

  path[0]='\0';
  strcat(path, directory);
  strcat(path, "/");
  strcat(path, name);
  strcat(path, ".desp");

  fp  = fopen(path,  "rb+");
  check_file_open(fp, path);

  err = fread(&matrix->siz,      sizeof(SA_TYPE),  1, fp);
  check_file_read(err, 1, path);

  err = fread(&matrix->n_desp,   sizeof(SA_TYPE),  1, fp);
  check_file_read(err, 1, path);

  err = fread(&matrix->m_desp,   sizeof(SA_TYPE),  1, fp);
  check_file_read(err, 1, path);

	matrix->desp = (SA_TYPE **) malloc(matrix->n_desp * sizeof(SA_TYPE *));
	check_malloc(matrix->desp, path);

  for (SA_TYPE i=0; i<matrix->n_desp; i++) {
    matrix->desp[i] = (SA_TYPE *) malloc(matrix->m_desp * sizeof(SA_TYPE));
    check_malloc(matrix->desp[i], path);

    err = fread(matrix->desp[i], sizeof(SA_TYPE), matrix->m_desp, fp);
    check_file_read(err, matrix->m_desp, path);
  }

	fclose(fp);

#if defined FM_COMP_32 || FM_COMP_64

  path[0]='\0';
  strcat(path, directory);
  strcat(path, "/");
  strcat(path, name);
  strcat(path, ".count");

  fp  = fopen(path,  "rb+");
  check_file_open(fp, path);

  err = fread(&matrix->n_count,   sizeof(SA_TYPE),  1, fp);
  check_file_read(err, 1, path);

  err = fread(&matrix->m_count,   sizeof(SA_TYPE),  1, fp);
  check_file_read(err, 1, path);

	matrix->count = (FM_COMP_TYPE **) malloc(matrix->n_count * sizeof(FM_COMP_TYPE *));
	check_malloc(matrix->count, path);

  for (SA_TYPE i=0; i<matrix->n_count; i++){
    matrix->count[i] = (FM_COMP_TYPE *) malloc(matrix->m_count * sizeof(FM_COMP_TYPE));
    check_malloc(matrix->count[i], path);
    err = fread(matrix->count[i], sizeof(FM_COMP_TYPE), matrix->m_count, fp);
    check_file_read(err, matrix->m_count, path);
  }

  fclose(fp);

#endif

}

void save_vector(vector *vector, const char *directory, const char *name) {

  size_t err=0;
  FILE *fp;

  char path[500];

  path[0]='\0';
  strcat(path, directory);
  strcat(path, "/");
  strcat(path, name);
  strcat(path, ".vec");

  fp  = fopen(path,  "wb+");
  check_file_open(fp, path);

  err = fwrite(&vector->n,     sizeof(SA_TYPE), 1, fp);
  check_file_write(err, 1, path);

  err = fwrite(vector->vector, sizeof(SA_TYPE), vector->n, fp);
  check_file_write(err, vector->n, path);

  fclose(fp);

}

void save_comp_vector(comp_vector *vector, const char *directory, const char *name) {

  size_t err=0;
  FILE *fp;

  char path[500];

  path[0]='\0';
  strcat(path, directory);
  strcat(path, "/");
  strcat(path, name);
  strcat(path, ".vec");

  fp  = fopen(path,  "wb+");
  check_file_open(fp, path);

  err = fwrite(&vector->siz,   sizeof(SA_TYPE), 1, fp);
  check_file_write(err, 1, path);

  err = fwrite(&vector->n,     sizeof(SA_TYPE), 1, fp);
  check_file_write(err, 1, path);

  err = fwrite(&vector->ratio, sizeof(SA_TYPE), 1, fp);
  check_file_write(err, 1, path);

  err = fwrite(vector->vector, sizeof(SA_TYPE), vector->n, fp);
  check_file_write(err, vector->n, path);

  fclose(fp);

}

void save_comp_matrix(comp_matrix *matrix, const char *directory, const char *name) {

  size_t err=0;
  FILE *fp;

  char path[500];

  path[0]='\0';
  strcat(path, directory);
  strcat(path, "/");
  strcat(path, name);
  strcat(path, ".desp");

  fp  = fopen(path,  "wb+");
  check_file_open(fp, path);

  err = fwrite(&matrix->siz,    sizeof(SA_TYPE), 1, fp);
  check_file_write(err, 1, path);

  err = fwrite(&matrix->n_desp, sizeof(SA_TYPE), 1, fp);
  check_file_write(err, 1, path);

  err = fwrite(&matrix->m_desp, sizeof(SA_TYPE), 1, fp);
  check_file_write(err, 1, path);

  for (SA_TYPE i=0; i<matrix->n_desp; i++) {
    err = fwrite(matrix->desp[i], sizeof(SA_TYPE), matrix->m_desp, fp);
    check_file_write(err, matrix->m_desp, path);
  }

  fclose(fp);

#if defined FM_COMP_32 || FM_COMP_64

  path[0]='\0';
  strcat(path, directory);
  strcat(path, "/");
  strcat(path, name);
  strcat(path, ".count");

  fp  = fopen(path,  "wb+");
  check_file_open(fp, path);

  err = fwrite(&matrix->n_count, sizeof(SA_TYPE), 1, fp);
  check_file_write(err, 1, path);

  err = fwrite(&matrix->m_count, sizeof(SA_TYPE), 1, fp);
  check_file_write(err, 1, path);

  for (SA_TYPE i=0; i<matrix->n_count; i++) {
    err = fwrite(matrix->count[i], sizeof(FM_COMP_TYPE), matrix->m_count, fp);
    check_file_write(err, matrix->m_count, path);
  }

  fclose(fp);

#endif

}
