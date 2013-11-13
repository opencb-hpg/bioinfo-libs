#include "../BW_preprocess.h"

int main(int argc, char **argv)
{

  comp_vector S, Scomp, Si, Scompi;
  comp_vector R, Rcomp, Ri, Rcompi;
  unsigned int ratio;

  if (argc!=3) {
    printf("Syntax:\n\t%s ratio dir\n", argv[0]);
    return 1;
  }

  ratio = atoi(argv[1]);

  readUIntCompVector(&S, argv[2], "S");
  printUIntVector(S.vector, S.n);
  calculateSRcomp(&S, &Scomp, ratio);
  printUIntVector(Scomp.vector, Scomp.n);
  saveUIntCompVector(&Scomp, argv[2], "Scomp");
  free(S.vector);
  free(Scomp.vector);

  readUIntCompVector(&Si, argv[2], "Si");
  printUIntVector(Si.vector, Si.n);
  calculateSRcomp(&Si, &Scompi, ratio);
  printUIntVector(Scompi.vector, Scompi.n);
  saveUIntCompVector(&Scompi, argv[2], "Scompi");
  free(Si.vector);
  free(Scompi.vector);

  readUIntCompVector(&R, argv[2], "R");
  printUIntVector(R.vector, R.n);
  calculateSRcomp(&R, &Rcomp, ratio);
  printUIntVector(Rcomp.vector, Rcomp.n);
  saveUIntCompVector(&Rcomp, argv[2], "Rcomp");
  free(R.vector);
  free(Rcomp.vector);

  readUIntCompVector(&Ri, argv[2], "Ri");
  printUIntVector(Ri.vector, Ri.n);
  calculateSRcomp(&Ri, &Rcompi, ratio);
  printUIntVector(Rcompi.vector, Rcompi.n);
  saveUIntCompVector(&Rcompi, argv[2], "Rcompi");
  free(Ri.vector);
  free(Rcompi.vector);

  return 0;

}
