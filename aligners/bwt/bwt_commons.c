#include "bwt_commons.h"

size_t cur_alloc=0, max_alloc=0;

void *mymalloc(size_t n)
{
  void *p;

  p = malloc(n);
  if (p == NULL) {
    printf("malloc failed.\n");
    exit(1);
  }
  cur_alloc += n;
  if (cur_alloc > max_alloc) {
    max_alloc = cur_alloc;
    //printf("allocated %ld\n",max_alloc);
  }

  if (n == 0) {
    printf("warning: 0 bytes allocated p=%p\n",p);
  }

  return p;
}

void *myrealloc(void *ptr, size_t next, size_t last)
{
  void *p;

  p = realloc(ptr, next);
  if (next > 0 && p == NULL) {
    printf("realloc failed. ptr=%p new=%zd old=%zd\n",ptr,next,last);
    exit(1);
  }
  cur_alloc += next - last;
  if (cur_alloc > max_alloc) {
    max_alloc = cur_alloc;
    //printf("allocated %ld\n",max_alloc);
  }

  return p;
}

void report_mem(const char *s)
{
  puts(s);
  printf("allocated total %zu   max %zu\n",cur_alloc,max_alloc);
}


void myfree(void *p, size_t s)
{
  free(p);
  cur_alloc -= s;
}


//-----------------------------------------------------------------------------

uint8_t nA;
uint8_t AA = -1, CC = -1, GG = -1, TT = -1;

uint8_t table[128];
char rev_table[4];
uint8_t reverse[4];

/*
void init_replace_table(const char *str) {

  if (str == NULL) {

    nA = 4;
    AA = 0; CC = 1; GG = 2; TT = 3;


    table['a'] = AA; table['A'] = AA;
    table['c'] = CC; table['C'] = CC;
    table['t'] = TT; table['T'] = TT;
    table['g'] = GG; table['G'] = GG;
    table['n'] = AA; table['N'] = AA;

    rev_table[AA] = 'A';
    rev_table[CC] = 'C';
    rev_table[GG] = 'G';
    rev_table[TT] = 'T';

  } else {

    nA = strlen(str);

    for (uint8_t i = 0; i < nA; i++) {
      rev_table[i] = toupper(str[i]);

      table[toupper(str[i])] = i;
      table[tolower(str[i])] = i;

      if      (toupper(str[i]) == 'A') AA = i;
      else if (toupper(str[i]) == 'C') CC = i;
      else if (toupper(str[i]) == 'G') GG = i;
      else if (toupper(str[i]) == 'T') TT = i;
    }

  }

  reverse[AA] = TT;
  reverse[CC] = GG;
  reverse[GG] = CC;
  reverse[TT] = AA;

}

void init_table() {
  init_replace_table(NULL);
}
*/

void encode_bases(uint8_t* dest, char* src, uintmax_t length, uint8_t *table) {
  for (uintmax_t i=0; i<length; i++)
    dest[i] = table[(uintmax_t)src[i]];
}

void decode_bases(char* dest, uint8_t* src, uintmax_t length, char *rev_table) {
  for (uintmax_t i=0; i<length; i++)
    dest[i] = rev_table[(uintmax_t)src[i]];
  dest[length] = '\0';
}

void revstring(uint8_t *X, uintmax_t nX) {

  char tmp;
  uintmax_t i, j;

  if (nX <= 1) return;

  for (i=0, j=nX-1; i<=j; i++, j--) {
    tmp = X[i];
    X[i] = X[j];
    X[j] = tmp;
  }

}

void revstrand(uint8_t *X, uintmax_t nX) {

  char tmp;
  uintmax_t i, j;

  if (nX <= 1) return;

  for (i=0, j=nX-1; i<=j; i++, j--) {
    tmp = reverse[X[i]];
    X[i] = reverse[X[j]];
    X[j] = tmp;
  }

}

void duplicate_reverse(uint8_t *X, uintmax_t nX) {

  uintmax_t i, j;

	if (nX == 0) return;

	for (i=0, j=nX*2-1; i<nX; i++, j--) {
    X[j] = reverse[X[i]];
  }

}
