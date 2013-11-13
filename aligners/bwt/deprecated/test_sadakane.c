#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sys/time.h>

#include "../BW_csafm.h"

int main(int argc, char **argv) {

  ref_vector B;

  check_syntax(argc, 3, "test_sadakane input_dir nucleotides");

  init_replace_table(argv[3]);

	read_ref_vector(&B, argv[1], "B");

	printf("dollar -> %u, size -> %u\n", B.dollar, B.n);

  for (unsigned int i=0; i<B.n; i++) {
		printf("%u\t-> %c\n", i, rev_table[B.vector[i]]);
	}

	printf("Memory frees\n");

	free(B.vector);

	return 0;

}
