/*
    bwa_gpu a set of tools which allow short sequence alignment using the Burrows-Wheeler
    transform usign both CPU and GPU approaches.
    Copyright (C) 2011  Jose Salavert Torres, Ignacio Blanquer Espert,
                        Andres Tomas Dominguez, Vicente Hernandez Garcia,
	 		Ignacio Medina Castello, Joaquin Tarraga Gimenez,
			Joaquin Dopazo Blazquez

    Contact e-mail: josator@fiv.upv.es, iblanque@dsic.upv.es, atomas@dsic.upv.es,
                    vhernand@dsic.upv.es, imedina@cipf.es, jtarraga@cipf.es,
                    jdopazo@cipf.es

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include "BW_search.h"

void BWExactFinalResultsBackward(char *W, vector *C, vector *C1, comp_matrix *O, results_list *rl_prev, results_list *rl_next) {

  BWiterationVariables();
  size_t k, l;
  int start, pos;

  result *r_iterator;

  for (size_t ii=0; ii < rl_prev->num_results; ii++) {

    r_iterator = &rl_prev->list[ii];

    start = r_iterator->start;
    pos   = r_iterator->pos;

    k = r_iterator->k;
    l = r_iterator->l;

    for(int i=pos; i>=start; i--) {
      BWiteration(k, l, k, l, W[i], C, C1, O);
      if (k > l) break;
    }

    if (k <= l) {
      change_result(r_iterator, k, l, start-1);
      add_result(r_iterator, rl_next);
    }

  } //r_prev

  rl_prev->num_results = 0;

}

void BWExactFinalResultsForward(char *W, vector *C, vector *C1, comp_matrix *O, results_list *rl_prev, results_list *rl_next) {

  BWiterationVariables();
  size_t k, l;
  int pos, end;

  result *r_iterator;

  for (size_t ii=0; ii < rl_prev->num_results; ii++) {

    r_iterator = &rl_prev->list[ii];

    pos   = r_iterator->pos;
    end   = r_iterator->end;

    k = r_iterator->k;
    l = r_iterator->l;

    for(int i=pos; i<=end; i++) {
      BWiteration(k, l, k, l, W[i], C, C1, O);
      if (k > l) break;
    }

    if (k <= l) {
      change_result(r_iterator, k, l, end+1);
      add_result(r_iterator, rl_next);
    }

  } //r_prev

  rl_prev->num_results = 0;

}

void BWExactPartialResultsBackward(char *W, vector *C, vector *C1, comp_matrix *O, results_list *rl_prev, results_list *rl_next) {

  BWiterationVariables();
  size_t k, l, k_next, l_next;
  int start, pos;

  result *r_iterator;
  size_t results, results_next;

  for (size_t ii=0; ii < rl_prev->num_results; ii++) {

    r_iterator = &rl_prev->list[ii];

    start = r_iterator->start;
    pos   = r_iterator->pos;

    k_next = r_iterator->k;
    l_next = r_iterator->l;
    results_next = l_next - k_next;

    for(int i=pos; i>=start; i--) {

      k = k_next;
      l = l_next;

      if (k > l) break;

      BWiteration(k, l, k_next, l_next, W[i], C, C1, O);
      results      = results_next;
      results_next = l_next - k_next;
      if (results == results_next) continue;

      change_result(r_iterator, k, l, i);
      add_result(r_iterator, rl_next);

    }

    if (k_next <= l_next) {
      change_result(r_iterator, k_next, l_next, start-1);
      add_result(r_iterator, rl_next);
    }

  } //r_prev

  rl_prev->num_results = 0;

}

void BWExactPartialResultsForward(char *W, vector *C, vector *C1, comp_matrix *O, results_list *rl_prev, results_list *rl_next) {

  BWiterationVariables();
  size_t k, l, k_next, l_next;
  int pos, end;

  result *r_iterator;
  size_t results, results_next;

  for (size_t ii=0; ii < rl_prev->num_results; ii++) {

    r_iterator = &rl_prev->list[ii];

    pos   = r_iterator->pos;
    end   = r_iterator->end;

    k_next = r_iterator->k;
    l_next = r_iterator->l;
    results_next = l_next - k_next;

    for(int i=pos; i<=end; i++) {

      k = k_next;
      l = l_next;

      if (k > l) break;

      BWiteration(k, l, k_next, l_next, W[i], C, C1, O);
      results      = results_next;
      results_next = l_next - k_next;
      if (results == results_next) continue;

      change_result(r_iterator, k, l, i);
      add_result(r_iterator, rl_next);

    }

    if (k_next <= l_next) {
      change_result(r_iterator, k_next, l_next, end+1);
      add_result(r_iterator, rl_next);
    }

  } //r_prev

  rl_prev->num_results = 0;

}

//TODO: Estudiar si evitar inserciones-delecciones-inserciones beneficia cuando partes de un fragmento con búsqueda exacta
void BWBranchPartialResultsBackward(char *W, vector *C, vector *C1, comp_matrix *O, results_list *rl_prev, results_list *rl_next) {

  BWiterationVariables();
  size_t k, l, k_aux, l_aux;
  int start, pos;
  int r_num_mismatches, no_previous;
  int last_err_pos;
  int last_err_kind = 0, last_err_base = -1;

  result *r_iterator;

  for (size_t ii=0; ii < rl_prev->num_results; ii++) {

    r_iterator = &rl_prev->list[ii];

    start = r_iterator->start;
    pos   = r_iterator->pos;

    if (pos < start) {
      add_result(r_iterator, rl_next);
      continue;
    }

    no_previous = 1;
    r_num_mismatches = r_iterator->num_mismatches-1;
    if (r_num_mismatches>-1) {
      last_err_pos  = r_iterator->err_pos[r_num_mismatches];
      last_err_kind = r_iterator->err_kind[r_num_mismatches];
      last_err_base = r_iterator->err_base[r_num_mismatches];
    } else {
      last_err_pos  = start-1;
    }

    k = r_iterator->k;
    l = r_iterator->l;

    add_mismatch(r_iterator, DELETION, -1, pos);

    if (last_err_pos == pos + 1) { //Previous MISMATCH or DELETION

      if (last_err_kind == MISMATCH) { //Previous MISMATCH
	
	//Deletion
	if (W[pos]!=last_err_base) {
	  change_result(r_iterator, k, l, pos-1);
	  add_result(r_iterator, rl_next);
	}

	for (int b=0;b<nA;b++) {

	  BWiteration(k, l, k_aux, l_aux, b, C, C1, O);
      
	  if (k_aux > l_aux) continue;
	  
	  //Insertion
	  if (b!=W[last_err_pos]) {
	    change_result(r_iterator, k_aux, l_aux, pos);
	    modify_last_mismatch2(r_iterator, INSERTION, b);
	    add_result(r_iterator, rl_next);
	  }

	  //Mismatch
	  if (b!=W[pos]) {
	    change_result(r_iterator, k_aux, l_aux, pos-1);
	    modify_last_mismatch2(r_iterator, MISMATCH, b);
	    add_result(r_iterator, rl_next);
	  }

	}

	no_previous = 0;

      } else if (last_err_kind == DELETION) { //Previous DELETION

	//Deletion
	change_result(r_iterator, k, l, pos-1);
	add_result(r_iterator, rl_next);

	for (int b=0;b<nA;b++) {

	  BWiteration(k, l, k_aux, l_aux, b, C, C1, O);

	  if (k_aux > l_aux) continue;

	  // NO INSERTION

	  if (b!=W[pos]) { //Mismatch

	    if (b!=W[last_err_pos]) {
	      change_result(r_iterator, k_aux, l_aux, pos-1);
	      modify_last_mismatch2(r_iterator, MISMATCH, b);
	      add_result(r_iterator, rl_next);
	    }

	  }

	}

	no_previous = 0;

      }

    } else if (last_err_pos == pos) { //Previous INSERTION

      //NO DELETION
      
      for (int b=0;b<nA;b++) {

	BWiteration(k, l, k_aux, l_aux, b, C, C1, O);
	
	if (k_aux > l_aux) continue;
      
	//Insertion
	change_result(r_iterator, k_aux, l_aux, pos);
	modify_last_mismatch2(r_iterator, INSERTION, b);
	add_result(r_iterator, rl_next);

	//Mismatch
	if (b!=W[pos]) {

	  if (W[pos]!=last_err_base) {
	    r_iterator->pos = pos-1;
	    modify_last_mismatch1(r_iterator, MISMATCH);
	    add_result(r_iterator, rl_next);
	  }

	}

      }

      no_previous = 0;

    }

    if (no_previous) { //Previous MATCH

      //Deletion
      change_result(r_iterator, k, l, pos-1);
      add_result(r_iterator, rl_next);

      for (int b=0;b<nA;b++) {

	BWiteration(k, l, k_aux, l_aux, b, C, C1, O);

	if (k_aux > l_aux) continue;

	//Insertion
	change_result(r_iterator, k_aux, l_aux, pos);
	modify_last_mismatch2(r_iterator, INSERTION, b);
	add_result(r_iterator, rl_next);

	if (b!=W[pos]) { //Mismatch
	  r_iterator->pos = pos-1;
	  modify_last_mismatch1(r_iterator, MISMATCH);
	  add_result(r_iterator, rl_next);
	}

      }

    }

  }

  rl_prev->num_results = 0;

}

void BWBranchPartialResultsForward(char *W, vector *C, vector *C1, comp_matrix *O, results_list *rl_prev, results_list *rl_next) {

  BWiterationVariables();
  size_t k, l, k_aux, l_aux;
  int end, pos;
  int r_num_mismatches, no_previous;
  int last_err_pos;
  int last_err_kind = 0, last_err_base = -1;

  result *r_iterator;

  for (size_t ii=0; ii < rl_prev->num_results; ii++) {

    r_iterator = &rl_prev->list[ii];

    end = r_iterator->end;
    pos   = r_iterator->pos;

    if (pos > end) {
      add_result(r_iterator, rl_next);
      continue;
    }

    no_previous = 1;
    r_num_mismatches = r_iterator->num_mismatches-1;
    if (r_num_mismatches>-1) {
      last_err_pos  = r_iterator->err_pos[r_num_mismatches];
      last_err_kind = r_iterator->err_kind[r_num_mismatches];
      last_err_base = r_iterator->err_base[r_num_mismatches];
    } else {
      last_err_pos  = end+1;
    }

    k = r_iterator->k;
    l = r_iterator->l;

    add_mismatch(r_iterator, DELETION, -1, pos);

    if (last_err_pos == pos - 1) { //Previous MISMATCH or DELETION

      if (last_err_kind == MISMATCH) { //Previous MISMATCH
	
	//Deletion
	if (W[pos]!=last_err_base) {
	  change_result(r_iterator, k, l, pos+1);
	  add_result(r_iterator, rl_next);
	}

	for (int b=0;b<nA;b++) {

	  BWiteration(k, l, k_aux, l_aux, b, C, C1, O);
      
	  if (k_aux > l_aux) continue;
	  
	  //Insertion
	  if (b!=W[last_err_pos]) {
	    change_result(r_iterator, k_aux, l_aux, pos);
	    modify_last_mismatch2(r_iterator, INSERTION, b);
	    add_result(r_iterator, rl_next);
	  }

	  //Mismatch
	  if (b!=W[pos]) {
	    change_result(r_iterator, k_aux, l_aux, pos+1);
	    modify_last_mismatch2(r_iterator, MISMATCH, b);
	    add_result(r_iterator, rl_next);
	  }

	}

	no_previous = 0;

      } else if (last_err_kind == DELETION) { //Previous DELETION

	//Deletion
	change_result(r_iterator, k, l, pos+1);
	add_result(r_iterator, rl_next);

	for (int b=0;b<nA;b++) {

	  BWiteration(k, l, k_aux, l_aux, b, C, C1, O);

	  if (k_aux > l_aux) continue;

	  // NO INSERTION

	  if (b!=W[pos]) { //Mismatch

	    if (b!=W[last_err_pos]) {
	      change_result(r_iterator, k_aux, l_aux, pos+1);
	      modify_last_mismatch2(r_iterator, MISMATCH, b);
	      add_result(r_iterator, rl_next);
	    }

	  }

	}

	no_previous = 0;

      }

    } else if (last_err_pos == pos) { //Previous INSERTION

      //NO DELETION
      
      for (int b=0;b<nA;b++) {

	BWiteration(k, l, k_aux, l_aux, b, C, C1, O);
	
	if (k_aux > l_aux) continue;
      
	//Insertion
	change_result(r_iterator, k_aux, l_aux, pos);
	modify_last_mismatch2(r_iterator, INSERTION, b);
	add_result(r_iterator, rl_next);

	//Mismatch
	if (b!=W[pos]) {

	  if (W[pos]!=last_err_base) {
	    r_iterator->pos = pos+1;
	    modify_last_mismatch1(r_iterator, MISMATCH);
	    add_result(r_iterator, rl_next);
	  }

	}

      }

      no_previous = 0;

    }

    if (no_previous) { //Previous MATCH

      //Deletion
      change_result(r_iterator, k, l, pos+1);
      add_result(r_iterator, rl_next);

      for (int b=0;b<nA;b++) {

	BWiteration(k, l, k_aux, l_aux, b, C, C1, O);

	if (k_aux > l_aux) continue;

	//Insertion
	change_result(r_iterator, k_aux, l_aux, pos);
	modify_last_mismatch2(r_iterator, INSERTION, b);
	add_result(r_iterator, rl_next);

	if (b!=W[pos]) { //Mismatch
	  r_iterator->pos = pos+1;
	  modify_last_mismatch1(r_iterator, MISMATCH);
	  add_result(r_iterator, rl_next);
	}

      }

    }

  }

  rl_prev->num_results = 0;

}

void BWExactSearchVectorBackward(char *W, int start, int end, size_t k, size_t l, size_t *vec_k, size_t *vec_l, vector *C, vector *C1, comp_matrix *O) {

  if (k > l)       return;
  if (start > end) return;

  BWiterationVariables();
  size_t k2, l2;
  int last;
  int i, j;

  last = end-start;

  k2 = k;
  l2 = l;

  for(i=end, j=last; i>=start; i--, j--) {

    BWiteration(k2, l2, k2, l2, W[i], C, C1, O);

    vec_k[j] = k2;
    vec_l[j] = l2;

    if (k2 > l2) {
      i--; j--;
      break;
    }

  }

  for(;i>=start; i--, j--) {
    vec_k[j] = k2;
    vec_l[j] = l2;
  }

}

void BWExactSearchVectorForward(char *W, int start, int end, size_t k, size_t l, size_t *vec_k, size_t *vec_l, vector *C, vector *C1, comp_matrix *O) {

  //printf("init search\n");
  if (k > l) return;
  if (start > end) return;

  //printf("init variables\n");
  BWiterationVariables();
  size_t k2, l2;
  //int last;
  int i, j;

  //last = end-start;

  k2 = k;
  l2 = l;
  //printf("init for\n");

  //printf("init for 1\n");
  for(i=start, j=0; i<=end; i++, j++) {
    //    printf("init variables %d\n", i);

    BWiteration(k2, l2, k2, l2, W[i], C, C1, O);

    //printf("update variables %d\n", i);
    vec_k[j] = k2;
    vec_l[j] = l2;

    if (k2 > l2) {
      i++; j++;
      break;
    }

  }
  //printf("end for\n");

  //printf("end for 1, init for 2\n");
  for(; i<=end; i++, j++) {
    vec_k[j] = k2;
    vec_l[j] = l2;
  }

  //printf("end for 2\n");
}

void BWSearchCPUBackward(char *W, vector *C, vector *C1, comp_matrix *O, result *res, results_list *rl_prev, results_list *rl_next, int num_errors) {

  init_result(res, 0);
  res->pos = res->end;
  add_result(res, rl_prev);

  while (num_errors > 0) {

    BWExactPartialResultsBackward(W, C, C1, O, rl_prev, rl_next);
    printf("next -> %u\n", rl_next->num_results);

    BWBranchPartialResultsBackward(W, C, C1, O, rl_next, rl_prev);
    printf("prev -> %u\n", rl_prev->num_results);

    num_errors--;

  }

  BWExactFinalResultsBackward(W, C, C1, O, rl_prev, rl_next);
  printf("next -> %u\n", rl_next->num_results);

}

void BWSearchCPUForward(char *W, vector *C, vector *C1, comp_matrix *O, result *res, results_list *rl_prev, results_list *rl_next, int num_errors) {

  init_result(res, 1);
  res->pos = res->start;
  add_result(res, rl_prev);

  while (num_errors > 0) {

    BWExactPartialResultsForward(W, C, C1, O, rl_prev, rl_next);
    printf("next -> %u\n", rl_next->num_results);

    BWBranchPartialResultsForward(W, C, C1, O, rl_next, rl_prev);
    printf("prev -> %u\n", rl_prev->num_results);

    num_errors--;

  }

  BWExactFinalResultsForward(W, C, C1, O, rl_prev, rl_next);
  printf("next -> %u\n", rl_next->num_results);

}

//TODO: Comprobar que funciona igual y borrar el codigo de esta función
void BWSearch1(char *W, int start, int end, size_t *vec_k, size_t *vec_l, size_t *vec_ki, size_t *vec_li, vector *C, vector *C1, comp_matrix *O, comp_matrix *Oi, results_list *r_list) {

  BWiterationVariables();

  size_t _k, _l, _ki, _li, _k_aux, _l_aux, _ki_aux, _li_aux;
  size_t results, results_last;

  unsigned int i, j, half, n;

  result r;

  n = end - start;
  half = n / 2;

  init_result(&r, 0);
  bound_result(&r, start, end);

  if (vec_k[0] <= vec_l[0]) {
    change_result(&r, vec_k[0], vec_l[0], -1);
    add_result(&r, r_list);
  }

  add_mismatch(&r, MATCH, -1, start);

  results = vec_l[0] - vec_k[0];

  results_last = results;
  _k  = vec_k[1];
  _l  = vec_l[1];
  results = _l  - _k;

  //printf("B -> %d: %d -> %d, %u, %u\n", 0, results, results_last, _k, _l);

  if (results != results_last) {

    //printf("*B -> %d: %d -> %d, %u, %u\n", 0, results, results_last, _k, _l);

    //printf("%d: %d -> %d\n", 0, results, results_last);

    //Deletion
    change_result(&r, _k, _l, -1);
    modify_last_mismatch3(&r, DELETION, -1, start);
    add_result(&r, r_list);

    for (int b=0;b<nA;b++) {

      BWiteration(_k, _l, _k_aux, _l_aux, b, C, C1, O);
      //printf("W -> %d, %d - %d\n", b, _k_aux, _l_aux);

      if (_k_aux > _l_aux) continue;
      //printf("*W -> %d, %d - %d\n", b, _k_aux, _l_aux);

      int b_w = (int) W[start];

      //Missmatch
      if (b!=b_w) {
    	change_result(&r, _k_aux, _l_aux, -1);
    	modify_last_mismatch2(&r, MISMATCH, b);
    	add_result(&r, r_list);
      }

      //Insertion
      BWiteration(_k_aux, _l_aux, _k_aux, _l_aux, b_w, C, C1, O);

      if (_k_aux <= _l_aux) {
    	change_result(&r, _k_aux, _l_aux, -1);
    	modify_last_mismatch2(&r, INSERTION, b);
    	add_result(&r, r_list);
      }

    }

  }

  for (i=start+2, j=2; j<=half; i++, j++) {

    results_last = results;
    _k = vec_k[j];
    _l = vec_l[j];
    results = _l  - _k;

    //printf("B -> %d: %d -> %d, %u, %u\n", j-1, results, results_last, _k, _l);

    if (results == results_last) continue;

    //printf("*B -> %d: %d -> %d, %u, %u\n", j-1, results, results_last, _k, _l);

    //Deletion
    change_result(&r, _k, _l, i-2);
    BWExactSearchBackward(W, C, C1, O, &r);
    if (r.k<=r.l) {
      modify_last_mismatch3(&r, DELETION, -1, i-1);
      add_result(&r, r_list);
    }

    for (int b=0;b<nA;b++) {
      /////////////////////////////////////////////////////////////////////////////
      BWiteration(_k, _l, _k_aux, _l_aux, b, C, C1, O);

      if (_k_aux > _l_aux) continue;

      //Insertion
      change_result(&r, _k_aux, _l_aux, i-1);
      BWExactSearchBackward(W, C, C1, O, &r);
      if (r.k<=r.l) {
	modify_last_mismatch3(&r, INSERTION, b, i-1);
	add_result(&r, r_list);
      }

      //Mismatch
      if (b!=(int)W[i-1]) {
	change_result(&r, _k_aux, _l_aux, i-2);
    	BWExactSearchBackward(W, C, C1, O, &r);
    	if (r.k<=r.l) {
	  modify_last_mismatch3(&r, MISMATCH, b, i-1);
	  add_result(&r, r_list);
	}

      }

    }

  }

  //printf("\n");

  half--;
  results = vec_li[n] - vec_ki[n];

  r.dir=1; //Change direction

  results_last = results;
  _ki  = vec_ki[n-1];
  _li  = vec_li[n-1];
  results = _li - _ki;

  //printf("F-> %d: %d -> %d, %u, %u\n", n, results, results_last, _ki, _li);

  if (results != results_last) {

    //printf("*F -> %d: %d -> %d, %u - %u\n", i+1, results, results_last, _ki, _li);

    //Deletion
    change_result(&r, _ki, _li, -1);
    modify_last_mismatch3(&r, DELETION, -1, end);
    add_result(&r, r_list);

    for (int b=0;b<nA;b++) {

      BWiteration(_ki, _li, _ki_aux, _li_aux, b, C, C1, Oi);

      if (_ki_aux > _li_aux) continue;

      int b_w = (int) W[end];

      //Mismatch
      if (b!=b_w) {
    	change_result(&r, _ki_aux, _li_aux, -1);
    	modify_last_mismatch2(&r, MISMATCH, b);
    	add_result(&r, r_list);
      }

      //Insertion
      BWiteration(_ki_aux, _li_aux, _ki_aux, _li_aux, b_w, C, C1, Oi);

      //printf("\tI -> %d - %d\n", _ki_aux, _li_aux);

      if (_ki_aux <= _li_aux){
    	change_result(&r, _ki_aux, _li_aux, -1);
    	modify_last_mismatch2(&r, INSERTION, b);
    	add_result(&r, r_list);
      }

    }

  }

  for(i=end-2,j=n-2; j>=half; i--, j--) {

    results_last = results;
    _ki  = vec_ki[j];
    _li  = vec_li[j];
    results = _li - _ki;

    //printf("F -> %d: %d -> %d, %u - %u\n", i+1, results, results_last, _ki, _li);

    if (results == results_last) continue;

    //printf("*F -> %d: %d -> %d, %u - %u\n", i+1, results, results_last, _ki, _li);

    //Deletion
    change_result(&r, _ki, _li, i+2);
    BWExactSearchForward(W, C, C1, Oi, &r);
    if (r.k<=r.l) {
      modify_last_mismatch3(&r, DELETION, -1, i+1);
      add_result(&r, r_list);
    }

    for (int b=0;b<nA;b++) {

      BWiteration(_ki, _li, _ki_aux, _li_aux, b, C, C1, Oi);

      //printf("W -> %d, %d - %d\n", b, _ki_aux, _li_aux);

      if (_ki_aux > _li_aux) continue;

      //Insertion
      change_result(&r, _ki_aux, _li_aux, i+1);
      BWExactSearchForward(W, C, C1, Oi, &r);
      if (r.k<=r.l) {
	modify_last_mismatch3(&r, INSERTION, b, i+1);
	add_result(&r, r_list);
      }

      //Mismatch
      if (b!= (int) W[i+1]) {
	change_result(&r, _ki_aux, _li_aux, i+2);
    	BWExactSearchForward(W, C, C1, Oi, &r);
    	if (r.k<=r.l) {
	  modify_last_mismatch3(&r, MISMATCH, b, i+1);
	  add_result(&r, r_list);
	}

      }

    }

  }

  //printf("\n");

}

void BWSearch1CPU(char *W, vector *C, vector *C1, comp_matrix *O, comp_matrix *Oi, result *res, results_list *r_list) {

  unsigned int half, n;
  int start, end;
  unsigned int _k, _l;

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

}

void BWSimpleSearch1Backward(char *W, vector *C, vector *C1, comp_matrix *O, result *res, results_list *r_list) {

  BWiterationVariables();
  size_t _k,_l, _k_next, _l_next, _k_aux, _l_aux;
  size_t results, results_next;
  int start, end, i;

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
      return;
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

    for (int b=0;b<nA;b++) {

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

}

void BWSimpleSearch1Forward(char *W, vector *C, vector *C1, comp_matrix *O, result *res, results_list *r_list) {

  BWiterationVariables();
  size_t _k, _l, _k_next, _l_next, _k_aux, _l_aux;
  size_t results, results_next;
  int start, end, i;

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
      return;
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

    for (int b=0;b<nA;b++) {

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

}
