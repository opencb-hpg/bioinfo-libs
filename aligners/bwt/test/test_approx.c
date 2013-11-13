/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#if 1
 #include <sys/timeb.h>
#else
 #include <sys/time.h>
 #include <sys/resource.h>
#endif
#include "csa.h"
#include "approx.h"

#if 1
typedef struct timeb mytimestruct;
void mygettime(mytimestruct *t)
{
  ftime(t);
}
double mylaptime(mytimestruct *before,mytimestruct *after)
{
  double t;
  t = after->time - before->time;
  t += (double)(after->millitm - before->millitm)/1000;
  return t;
}
#else
typedef mytimestruct struct rusage;
void mygettime(mytimestruct *t)
{
  getrusage(RUSAGE_SELF,t);
}
double mylaptime(mytimestruct *before,mytimestruct *after)
{
  double t;
  t = after->ru_utime.tv_sec - before->ru_utime.tv_sec;
  t += (double)(after->ru_utime.tv_usec
                - before->ru_utime.tv_usec)/1000000;
  return t;
}
#endif

void test_approxsearch(CSA *SA)
{
   long i,n_pos;
   int keylen;
   double t;
   unsigned char key[4096];
   char *buf,*c_ptr;
   mytimestruct before,after;
   approx_list *list;
   int max_error;

   while (!feof(stdin)) {
//      printf("\ninput score & key ");  fflush(stdout);
//      scanf("%ld %s",&min_score,key);
      printf("\ninput error len pos ");  fflush(stdout);
      scanf(" %d %d %ld",&max_error, &keylen, &i);
//      i = rand() % (SA->n-keylen-1) + 1;
//      SA->substring(key, SA, i, keylen);
//      i = 200000;
      printf("error = %d len = %d pos = %ld\n", max_error, keylen, i);
      SA->text(key,SA,i,i+keylen-1);
      key[keylen] = 0;
      printf("key %s\n",key);

      mygettime(&before);
      list = csa_approxsearch2(key, keylen, max_error, SA);
      mygettime(&after);
      t = mylaptime(&before,&after);
      printf("search %f sec\n",t);

      printf("search results\n");
      approx_list_print(SA, list);

      approx_list_free(list);

   }
}


int main(int argc, char *argv[])
{
  CSA SA;

  csa_read(&SA,argc-1,argv+1);

  test_approxsearch(&SA);
}
