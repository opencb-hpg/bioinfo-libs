#include "BW_preprocess.h"

int main(int argc, char **argv)
{

  ref_vector X, B, Bi;
  vector C, C1;
  comp_vector S, Si, Scomp, Scompi;
  comp_vector R, Ri, Rcomp, Rcompi;
  comp_matrix O, Oi;

  int s_ratio;

  exome ex;

	check_syntax(argc, 5, "preprocess_debug ref_file output_dir s_ratio nucleotides");

  timevars();
	init_replace_table(argv[4]);

  s_ratio = atoi(argv[3]);

  encode_reference(&X, &ex, true, argv[1]);
  save_exome_file(&ex, argv[2]);

  tic("Calculating BWT");
  calculateBWTdebug(&B, &S, &X, 0);
  toc();

  save_ref_vector(&X, argv[2], "X");

  print_vector(S.vector, S.n);
  print_vector(B.vector, B.n);

  tic("Calculating prefix-trie matrices C and O");
  calculate_C(&C, &C1, &B);
  calculate_O(&O, &B);
  toc();

  print_vector(C.vector, C.n);
  print_vector(C1.vector, C1.n);
  print_comp_matrix(O);

  save_ref_vector(&B, argv[2], "B");
  free(B.vector);
  save_vector(&C, argv[2], "C");
  free(C.vector);
  save_vector(&C1, argv[2], "C1");
  free(C1.vector);
  save_comp_matrix(&O, argv[2], "O");
  free_comp_matrix(NULL, &O);

  tic("Calculating R");
  calculate_R(&R, &S);
  toc();
  print_vector(R.vector, R.n);

  tic("Calculating Scomp Rcomp");
  compress_SR(&S, &Scomp, s_ratio);
  print_vector(Scomp.vector, Scomp.n);
  compress_SR(&R, &Rcomp, s_ratio);
  print_vector(Rcomp.vector, Rcomp.n);
  toc();

  save_comp_vector(&S, argv[2], "S");
  free(S.vector);
  save_comp_vector(&R, argv[2], "R");
  free(R.vector);
  save_comp_vector(&Scomp, argv[2], "Scomp");
  free(Scomp.vector);
  save_comp_vector(&Rcomp, argv[2], "Rcomp");
  free(Rcomp.vector);

  tic("Calculating BWT of reverse reference");
  calculateBWTdebug(&Bi, &Si, &X, 1);
  toc();

  save_ref_vector(&X, argv[2], "Xi");

  print_vector(Bi.vector, Bi.n);
  print_vector(Si.vector, Si.n);

  tic("Calculating inverted prefix-trie matrix Oi");
  calculate_O(&Oi, &Bi);
  toc();

  free(X.vector);

	print_comp_matrix(Oi);

	save_ref_vector(&Bi, argv[2], "Bi");
  free(Bi.vector);

	save_comp_matrix(&Oi, argv[2], "Oi");
  free_comp_matrix(NULL, &Oi);

  tic("Calculating Ri");
  calculate_R(&Ri, &Si);
  toc();

  print_vector(Ri.vector, Ri.n);

  tic("Calculating Scompi Rcompi");
  compress_SR(&Si, &Scompi, s_ratio);
  print_vector(Scompi.vector, Scompi.n);
  compress_SR(&Ri, &Rcompi, s_ratio);
  print_vector(Rcompi.vector, Rcompi.n);
  toc();

  save_comp_vector(&Si, argv[2], "Si");
  free(Si.vector);
  save_comp_vector(&Ri, argv[2], "Ri");
  free(Ri.vector);
  save_comp_vector(&Scompi, argv[2], "Scompi");
  free(Scompi.vector);
  save_comp_vector(&Rcompi, argv[2], "Rcompi");
  free(Rcompi.vector);

  return 0;

}
