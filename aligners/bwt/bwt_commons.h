#ifndef BWT_COMMONS_H
#define BWT_COMMONS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <inttypes.h>

#define MAXLINE      500
#define MAXLINECOMP  125

//The meaning of Insertion and Deletion may be swapped to the meaning of the SAM file format
#define MATCH	  0
#define DELETION  1
#define MISMATCH  2
#define INSERTION 3

#ifndef min
#define min(x,y) (((x)<(y))?(x):(y))
#endif

#ifndef max
#define max(x,y) (((x)>(y))?(x):(y))
#endif

#define check_syntax(argc, n, params)\
  if ((argc)!=(n)) {\
    fprintf(stderr, "Syntax:\n\t%s\n", (params));\
    exit(1);\
  }

#define check_malloc(D, path)\
  if ((D)==NULL) {\
    fprintf(stderr, "Data structure " #D " in %s is too large\n", (path));\
    exit(1);\
  }

#define check_file_open(fp, path)\
  if (!(fp)) {\
    fprintf(stderr, "Error opening file: %s\n", (path));\
    exit(1);\
  }

#define check_file_read(err, nmemb, path)\
  if ((err) != (nmemb)) {\
    fprintf(stderr, "Error reading file '%s'\n", (path));\
    exit(1);\
  }

#define check_file_write(err, nmemb, path)\
  if ((err) != (nmemb)) {\
    fprintf(stderr, "Error writing file '%s'\n", (path));\
    exit(1);\
  }

#ifdef VERBOSE_DBG

#define print_matrix(M,n,m)\
{\
  printf("Matrix " #M ":\n");\
  for (uintmax_t i_=0; i_<((uintmax_t) (n)); i_++) {\
    printf("%ju: ", (uintmax_t) i_);\
    for (uintmax_t j_=0; j_<((uintmax_t) (m)); j_++) {\
      printf("%ju ", (uintmax_t) (M)[i_][j_]);\
    }\
    printf("\n");\
  }\
}

#define print_vector(V,n)\
{\
  printf("Vector " #V ":\n");\
  for (uintmax_t i_=0; i_<((uintmax_t) (n)); i_++) {\
		printf("%ju ", (uintmax_t) (V)[i_]);\
  }\
  printf("\n");\
}

#define print_string(S)\
	printf("String " #S ":\n%s\n", S);

#else

#define print_matrix(M,n,m);
#define print_vector(V,n);
#define print_string(S);

#endif

void *mymalloc(size_t n);
void *myrealloc(void *ptr, size_t next, size_t last);
void myfree(void *p, size_t s);
void report_mem(const char *s);

extern size_t cur_alloc, max_alloc;


//-----------------------------------------------------------------------------


extern uint8_t nA;
extern uint8_t AA, CC, GG, TT;

extern uint8_t table[128];
extern char rev_table[4];
extern char reserve[4];

/**
 *  @brief Inits table for nucleotide coding/decoding 
 *  @return void
 * 
 *  Inits table[128] for nucleotide coding/decoding 
 */
void init_table();
void init_replace_table(const char *str);

/**
 *  @brief Encodes a sequence of plain nucleotides
 *  @param dest pointer to destination with encoded nucleotides
 *  @param src pointer to char with plain nucleotides
 *  @param length length of the nucleotide sequence
 * 
 *  Encodes a sequence of plain nucleotides
 */
void encode_bases(uint8_t* dest, char* src, uintmax_t length);

/**
 *  @brief Decodes a sequence of plain nucleotides
 *  @param dest pointer to destination with plain nucleotides
 *  @param src pointer to char with encoded nucleotides
 *  @param length length of the nucleotide sequence
 * 
 *  Encodes a sequence of plain nucleotides
 */
void decode_bases(char* dest, uint8_t* src, uintmax_t length);

void revstring(uint8_t *X, uintmax_t nX);
void revstrand(uint8_t *X, uintmax_t nX);

void duplicate_reverse(uint8_t *X, uintmax_t nX);


#endif
