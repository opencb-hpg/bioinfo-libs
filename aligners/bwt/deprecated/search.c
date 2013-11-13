bool BWSimpleSearch1Backward(uint8_t *W, vector *C, vector *C1, comp_matrix *O, result *res, results_list *r_list) {

	SA_TYPE _k,_l, _k_next, _l_next, _k_aux, _l_aux;
	SA_TYPE results, results_next;
	int16_t start, end, i;

	result r;

	start   = res->start;
	end     = res->pos;

	_k_next = res->k;
	_l_next = res->l;
	results_next = _l_next - _k_next;

	init_result(&r, 0);
	bound_result(&r, start, end);
	add_mismatch(&r, MATCH, -1, start);

	for(i=end; i>=start; i--) {

		_k = _k_next;
		_l = _l_next;

		//printf("%d:\n", i);

		if (_k > _l) {
			change_result(res, _k, _l, -1);
			return false;
		}

		BWiteration(_k, _l, _k_next, _l_next, W[i], C, C1, O);
		results      = results_next;
		results_next = _l_next - _k_next;
		//printf("(%lu, %lu, %lu)\t", results, _k, _l);

		if (results == results_next) continue;

		//Deletion
		change_result(&r, _k, _l, i-1);
		BWExactSearchBackward(W, C, C1, O, &r);
		if (r.k<=r.l) {
			modify_last_mismatch3(&r, DELETION, -1, i);
			add_result(&r, r_list);
		}

		for (uint8_t b=0;b<nA;b++) {

			BWiteration(_k, _l, _k_aux, _l_aux, b, C, C1, O);

			//printf("W -> %d, %d - %d\n", b, _k_aux, _l_aux);

			if (_k_aux > _l_aux) continue;

			//Insertion
			change_result(&r, _k_aux, _l_aux, i);
			BWExactSearchBackward(W, C, C1, O, &r);
			if (r.k<=r.l) {
				modify_last_mismatch3(&r, INSERTION, b, i);
				add_result(&r, r_list);
			}

			//Mismatch
			if (b!=(int)W[i]) {
				change_result(&r, _k_aux, _l_aux, i-1);
				BWExactSearchBackward(W, C, C1, O, &r);
				if (r.k<=r.l) {
					modify_last_mismatch3(&r, MISMATCH, b, i);
					add_result(&r, r_list);
				}

			}

		}

	}

	//Match at exit in res
	change_result(res, _k_next, _l_next, -1);

	return false;

}

bool BWSimpleSearch1Forward(uint8_t *W, vector *C, vector *C1, comp_matrix *O, result *res, results_list *r_list) {

	SA_TYPE _k, _l, _k_next, _l_next, _k_aux, _l_aux;
	SA_TYPE results, results_next;
	int16_t start, end, i;

	result r;

	start   = res->pos;
	end     = res->end;

	_k_next = res->k;
	_l_next = res->l;
	results_next = _l_next - _k_next;

	init_result(&r, 1);
	bound_result(&r, start, end);
	add_mismatch(&r, MATCH, -1, start);

	for(i=start; i<=end; i++) {

		_k = _k_next;
		_l = _l_next;

		//printf("%d:\n", i);

		if (_k > _l) {
			change_result(res, _k, _l, -1);
			return false;
		}

		BWiteration(_k, _l, _k_next, _l_next, W[i], C, C1, O);
		results      = results_next;
		results_next = _l_next - _k_next;
		if (results == results_next) continue;

		//Deletion
		change_result(&r, _k, _l, i+1);
		BWExactSearchForward(W, C, C1, O, &r);
		if (r.k<=r.l) {
			modify_last_mismatch3(&r, DELETION, -1, i);
			add_result(&r, r_list);
		}

		for (uint8_t b=0;b<nA;b++) {

			BWiteration(_k, _l, _k_aux, _l_aux, b, C, C1, O);

			//printf("W -> %d, %d - %d\n", b, _k_aux, _l_aux);

			if (_k_aux > _l_aux) continue;

			//Insertion
			change_result(&r, _k_aux, _l_aux, i);
			BWExactSearchForward(W, C, C1, O, &r);
			if (r.k<=r.l) {
				modify_last_mismatch3(&r, INSERTION, b, i);
				add_result(&r, r_list);
			}

			//Mismatch
			if (b!=(int)W[i]) {
				change_result(&r, _k_aux, _l_aux, i+1);
				BWExactSearchForward(W, C, C1, O, &r);
				if (r.k<=r.l) {
					modify_last_mismatch3(&r, MISMATCH, b, i);
					add_result(&r, r_list);
				}
			}

		}

	}

	//Match at exit in res
	change_result(res, _k_next, _l_next, -1);

	return false;

}

bool BWSearch1CPU(uint8_t *W, vector *C, vector *C1, comp_matrix *O, comp_matrix *Oi, result *res, results_list *r_list) {

	int16_t start, end, half, n;;
	SA_TYPE _k, _l;

	_k = res->k;
	_l = res->l;

	start = res->start;
	end   = res->end;

	n = end - start + 1;
	half = n / 2;

	result r;

	init_result(&r, 0);
	bound_result(&r, half, end);
	change_result(&r, _k, _l, end);

	BWExactSearchBackward(W, C, C1, O, &r);

	if (r.k <= r.l) {
		r.start = start;
		r.pos = half-1;
		BWSimpleSearch1Backward(W, C, C1, O, &r, r_list);

		if (r.k <= r.l) add_result(&r, r_list); //Match
	}

	half--;

	init_result(&r, 1);
	bound_result(&r, start, half);
	change_result(&r, _k, _l, start);

	BWExactSearchForward(W, C, C1, Oi, &r);
	if (r.k <= r.l) {
		r.pos = half+1;
		r.end = end;
		BWSimpleSearch1Forward(W, C, C1, Oi, &r, r_list);
	}

	return false;

}
