#include "../search/io.h"
#include "../bwt_commons.h"
#include "../dbwt/dbwt.h"
#include "../csalib/csa.h"

#include "commons/commons.h"

int main(int argc, char **argv)
{

	ref_vector X, Xi;

	exome ex;
	bwt_config_t bwt_config;

	check_syntax(argc, 5, "preprocess ref_file output_dir duplicate_reverse nucleotides");

	timevars();

	int duplicate_strand = atoi(argv[3]);
	save_config(argv[4], duplicate_strand, argv[2]);
	bwt_config.duplicate_strand = duplicate_strand;
	bwt_init_replace_table(&bwt_config, argv[4]);

	encode_reference(&X, &ex, argv[1], &bwt_config);
	save_ref_vector(&X, argv[2], "X");
	save_exome_file(&ex, duplicate_strand, argv[2]);

	tic("Calc. Backward Burrows-Wheeler Transform -> Sadakane direct SAIS");
	direct_bwt(X.vector, X.n, argv[2], "backward", false);
	toc();

	tic("Calc. Backward Suffix Array -> CSALIB DNA");
	csa_new_from_bwt_gnu_bwt_wrapper(argv[2], "backward");
	toc();

	read_ref_vector(&Xi, argv[2], "X");
	revstring(Xi.vector, Xi.n);
	save_ref_vector(&Xi, argv[2], "Xi");

	tic("Calc. Forward Burrows-Wheeler Transform -> Sadakane direct SAIS");
	direct_bwt(Xi.vector, Xi.n, argv[2], "forward", false);
	toc();

	tic("Calc. Forward Suffix Array -> CSALIB DNA");
	csa_new_from_bwt_gnu_bwt_wrapper(argv[2], "forward");
	toc();

	bwt_free_replace_table(&bwt_config);

	return 0;

}
