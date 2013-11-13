#include "dbwt.h"

void bwt_wrapper(char *fname)
{
  uchar *T;
  long n;
  FILE *fp;
//  MMAP *map;

///////////////////////////////////////////////
// File loading

#if 1
	size_t err;
  fp = fopen(fname,"rb");
	check_file_open(fp, fname);

  fseek(fp,0,SEEK_END);
  n = ftell(fp);
  fseek(fp,0,SEEK_SET);

  T = (uchar *) mymalloc((n+1)*sizeof(uchar));
  report_mem("read string");

	for (long i=0; i<n; i += (1<<20)) {
    size_t size;
    printf("%ld \r",i/(1<<20));
    fflush(stdout);
    size = min(n-i,1<<20);
		err = fread(&T[i],1,size,fp);
		check_file_read(err, size, fname);
  }
  fclose(fp);
#else
  map = mymmap(fname);
  n = map->len;
  T = map->addr;
#endif
///////////////////////////////////////////////

	direct_bwt(T, n, NULL, NULL, false);
}

int main(int argc, char *argv[])
{
  printf("sizeof(uchar *)=%lu sizeof(uint)=%lu sizeof(ulong)=%lu\n",sizeof(uchar *),sizeof(uint),sizeof(ulong));

	bwt_wrapper(argv[1]);
 
  return 0;
}

/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
