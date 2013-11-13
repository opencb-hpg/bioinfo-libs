/**
 *  @brief Encodes the sequence X putting 4 bases in each element of a byte vector 
 *  @param X Sequence to encode
 *  @param nX Length of the sequence to encode
 *  @param Y Encoded sequence
 *  @return The length of the encoded sequence
 * 
 */

uintmax_t comp4basesInByte(REF_TYPE *X, uintmax_t nX, uint8_t *Y) {

  uintmax_t desp=2, iaux=0;
  uint8_t aux=0;

  if (nX>0) aux = X[0];

  for (uintmax_t i=1; i<nX; i++) {

    if (desp==8) {
      desp=0;
      Y[iaux] = aux;
      iaux++ ;
      aux = aux & 0;
    }

    aux = aux | (X[i] << desp);
    desp = desp +2;

  }

  Y[iaux] = aux;

	uintmax_t nY = nX / 4;
  if (nX % 4) nY++;

  return nY;
}
