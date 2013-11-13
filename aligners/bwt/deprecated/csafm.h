inline SA_TYPE getScompValueB(SA_TYPE m, comp_vector *Scomp, vector *C, comp_matrix *O, ref_vector *B) {
  
  SA_TYPE i, j;
  uint8_t b_aux;
  
  i=m; j=0;

  while (i % Scomp->ratio) {
 
    if (i == B->dollar) {

      i=0;

		} else {

			if (i > B->dollar)
    		b_aux = B->vector[i-1];
			else
				b_aux = B->vector[i];

      i = C->vector[b_aux] + getO(b_aux, i+1/*0 is -1*/, O);

    }

    j++;

  }

  return (Scomp->vector[i / Scomp->ratio] + j) % (O->siz-1);

}

inline SA_TYPE getRcompValueB(SA_TYPE m, comp_vector *Rcomp, vector *C, comp_matrix *O, ref_vector *B) {
  SA_TYPE i, j, k;
  uint8_t b_aux;

  i = (Rcomp->ratio - (m % Rcomp->ratio)) % Rcomp->ratio;
  k = m + i;

  if(k < Rcomp->siz) {
    j = Rcomp->vector[k / Rcomp->ratio];
  } else {
    j = Rcomp->vector[0];
    i = Rcomp->siz - m;
  }

  while (i) {


    if (j == B->dollar) {

      j=0;

    } else {

			if (j > B->dollar)
				b_aux = B->vector[j-1];
			else
				b_aux = B->vector[j];

      j = C->vector[b_aux] + getO(b_aux, j+1/*0 is -1*/, O);

    }

    i--;

  }

	return j;

}
