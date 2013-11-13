__global__ void BWExactIterativeSearchGPU(REF_TYPE *W, SA_TYPE *nW, SA_TYPE *nWe, int *k, int *l, int k_ini, int l_ini, int *C,  int *O,  int sizO) {

	intmax_t i, b; //, pos;
	SA_TYPE k2, l2;
	REF_TYPE val1, val2, val3, val4;
	int siz1, siz2, siz3, siz4;

	int offset  = blockIdx.x * blockDim.x + threadIdx.x;

	__shared__ int Cshared[4];

	if (threadIdx.x<4) { // Minimum data is a 32 block
		Cshared[threadIdx.x] = C[threadIdx.x];
	}

	__syncthreads();

	k2 = k_ini;
	l2 = l_ini;

	// First block of 4 bases should not be fully filled

	int real_size = nWe[offset];
	int resto = real_size % 4;

	b = W[offset*MAXLINECOMP + nW[offset] - 1];

	switch(resto) {

		case 0:
			val4 = ( b >> 6 ) & 3;
			siz4 = val4*sizO;

			k2 = Cshared[val4] + O[siz4 + k2    ] + 1;
			l2 = Cshared[val4] + O[siz4 + l2 + 1];

		case 3:
			val3 = ( b >> 4 ) & 3;
			siz3 = val3*sizO;

			k2 = Cshared[val3] + O[siz3 + k2    ] + 1;
			l2 = Cshared[val3] + O[siz3 + l2 + 1];

		case 2:
			val2 = ( b >> 2 ) & 3;
			siz2 = val2*sizO;

			k2 = Cshared[val2] + O[siz2 + k2    ] + 1;
			l2 = Cshared[val2] + O[siz2 + l2 + 1];

		case 1:
			val1 = ( b      ) & 3;
			siz1 = val1*sizO;

			k2 = Cshared[val1] + O[siz1 + k2    ] + 1;
			l2 = Cshared[val1] + O[siz1 + l2 + 1];

	}

	__syncthreads();

	for (i=nW[offset]-2; (k2<=l2) && (i>=0); i--) {

		b = W[offset*MAXLINECOMP + i];

		val4 = ( b >> 6 ) & 3;
		siz4 = val4*sizO;

		k2 = Cshared[val4] + O[siz4 + k2    ] + 1;
		l2 = Cshared[val4] + O[siz4 + l2 + 1];

		val3 = ( b >> 4 ) & 3;
		siz3 = val3*sizO;

		k2 = Cshared[val3] + O[siz3 + k2    ] + 1;
		l2 = Cshared[val3] + O[siz3 + l2 + 1];

		val2 = ( b >> 2 ) & 3;
		siz2 = val2*sizO;

		k2 = Cshared[val2] + O[siz2 + k2    ] + 1;
		l2 = Cshared[val2] + O[siz2 + l2 + 1];

		val1 = ( b      ) & 3;
		siz1 = val1*sizO;

		k2 = Cshared[val1] + O[siz1 + k2    ] + 1;
		l2 = Cshared[val1] + O[siz1 + l2 + 1];

	}

	__syncthreads();

	k[offset] = k2;
	l[offset] = l2;

}

__global__ void BWExactIterativeSearchGPURev(REF_TYPE *W, SA_TYPE *nW, SA_TYPE *nWe, SA_TYPE *k, SA_TYPE *l, SA_TYPE k_ini, SA_TYPE l_ini, int *C,  int *O,  int sizO) {

	intmax_t i, b;//, pos;
	SA_TYPE k2, l2;
	char val1, val2, val3, val4;
	int siz1, siz2, siz3, siz4;

	int offset  = blockIdx.x * blockDim.x + threadIdx.x;

	__shared__ int Cshared[4];

	if (threadIdx.x<4) { // Minimum data is a 32 block
		Cshared[threadIdx.x] = C[threadIdx.x];
	}

	__syncthreads();

	k2 = k_ini;
	l2 = l_ini;

	for (i=0; (k2<=l2) && (i<nW[offset]-1); i++) {

		b = W[offset*MAXLINECOMP + i];

		val1 = ( b      ) & 3;
		siz1 = val1*sizO;

		k2 = Cshared[val1] + O[siz1 + k2    ] + 1;
		l2 = Cshared[val1] + O[siz1 + l2 + 1];

		val2 = ( b >> 2 ) & 3;
		siz2 = val2*sizO;

		k2 = Cshared[val2] + O[siz2 + k2    ] + 1;
		l2 = Cshared[val2] + O[siz2 + l2 + 1];

		val3 = ( b >> 4 ) & 3;
		siz3 = val3*sizO;

		k2 = Cshared[val3] + O[siz3 + k2    ] + 1;
		l2 = Cshared[val3] + O[siz3 + l2 + 1];

		val4 = ( b >> 6 ) & 3;
		siz4 = val4*sizO;

		k2 = Cshared[val4] + O[siz4 + k2    ] + 1;
		l2 = Cshared[val4] + O[siz4 + l2 + 1];

	}

	__syncthreads();

	if (k2<=l2) {

		// Last block of 4 bases should not be fully filled

		int real_size = nWe[offset];
		int resto = real_size % 4;

		b = W[offset*MAXLINECOMP + nW[offset] - 1];

		switch(resto) {

			case 0:

				val1 = ( b      ) & 3;
				siz1 = val1*sizO;

				k2 = Cshared[val1] + O[siz1 + k2    ] + 1;
				l2 = Cshared[val1] + O[siz1 + l2 + 1];

				val2 = ( b >> 2 ) & 3;
				siz2 = val2*sizO;

				k2 = Cshared[val2] + O[siz2 + k2    ] + 1;
				l2 = Cshared[val2] + O[siz2 + l2 + 1];

				val3 = ( b >> 4 ) & 3;
				siz3 = val3*sizO;

				k2 = Cshared[val3] + O[siz3 + k2    ] + 1;
				l2 = Cshared[val3] + O[siz3 + l2 + 1];

				val4 = ( b >> 6 ) & 3;
				siz4 = val4*sizO;

				k2 = Cshared[val4] + O[siz4 + k2    ] + 1;
				l2 = Cshared[val4] + O[siz4 + l2 + 1];

				break;

			case 3:

				val1 = ( b      ) & 3;
				siz1 = val1*sizO;

				k2 = Cshared[val1] + O[siz1 + k2    ] + 1;
				l2 = Cshared[val1] + O[siz1 + l2 + 1];

				val2 = ( b >> 2 ) & 3;
				siz2 = val2*sizO;

				k2 = Cshared[val2] + O[siz2 + k2    ] + 1;
				l2 = Cshared[val2] + O[siz2 + l2 + 1];

				val3 = ( b >> 4 ) & 3;
				siz3 = val3*sizO;

				k2 = Cshared[val3] + O[siz3 + k2    ] + 1;
				l2 = Cshared[val3] + O[siz3 + l2 + 1];

				break;

			case 2:

				val1 = ( b      ) & 3;
				siz1 = val1*sizO;

				k2 = Cshared[val1] + O[siz1 + k2    ] + 1;
				l2 = Cshared[val1] + O[siz1 + l2 + 1];

				val2 = ( b >> 2 ) & 3;
				siz2 = val2*sizO;

				k2 = Cshared[val2] + O[siz2 + k2    ] + 1;
				l2 = Cshared[val2] + O[siz2 + l2 + 1];

				break;

			case 1:
				val1 = ( b      ) & 3;
				siz1 = val1*sizO;

				k2 = Cshared[val1] + O[siz1 + k2    ] + 1;
				l2 = Cshared[val1] + O[siz1 + l2 + 1];

				break;

		}

	}

	__syncthreads();

	k[offset] = k2;
	l[offset] = l2;

}


