#ifndef COMMON_H
#define COMMON_H

#include <cstdint>
#include <immintrin.h>
#include "params.h"

typedef struct alignas(32) vec { uint16_t data[64]; } vec;

void inline vec_add(vec* out, const vec* a, const vec *b);
void inline vec_sub(vec* out, const vec* a, const vec *b);
void inline vec_reduce(vec* v);
void inline vec_mul_add(vec* out, const vec* a, uint16_t b);
unsigned char vec_hw(const vec* v);
int vec_is_permutation(vec* a, vec* b);

void print_vec(const vec* V);
void print_matrix(const vec* matrix, unsigned char rows);
void random_gaussian_elimination(vec* matrix, unsigned char *I, unsigned char k, mt19937 *mt);
void read_codes(vec* matrix1, vec* matrix2, vec* perm, vec* diagonal);
void dual_code(vec* M, vec* Mperp);


static const __m256i zeros = _mm256_set1_epi16(0);
static const __m256i ones  = _mm256_set1_epi16(1);
static const __m256i Qminusones = _mm256_set1_epi16(Q-1);
static const __m256i Qs = _mm256_set1_epi16(Q);
static const __m256i inverseofQ = _mm256_set1_epi16(Qinv);
static const __m256i multiples = _mm256_set1_epi16(Q*49);

void inline vec_add_reduce(vec* out, const vec* a, const vec *b){
	__m256i a1,a2,a3,a4,b1,b2,b3,b4,out1,out2,out3,out4,cmp1,cmp2,cmp3,cmp4;

    a1 = _mm256_load_si256((__m256i *)a->data);
    a2 = _mm256_load_si256((__m256i *)a->data + 1);
    a3 = _mm256_load_si256((__m256i *)a->data + 2);
    a4 = _mm256_load_si256((__m256i *)a->data + 3);
    b1 = _mm256_load_si256((__m256i *)b->data);
    b2 = _mm256_load_si256((__m256i *)b->data + 1);
    b3 = _mm256_load_si256((__m256i *)b->data + 2);
    b4 = _mm256_load_si256((__m256i *)b->data + 3);

    out1 = _mm256_add_epi16(a1,b1);
    out2 = _mm256_add_epi16(a2,b2);
    out3 = _mm256_add_epi16(a3,b3);
    out4 = _mm256_add_epi16(a4,b4);

    cmp1 = _mm256_cmpgt_epi16(out1,Qminusones);
    cmp2 = _mm256_cmpgt_epi16(out2,Qminusones);
    cmp3 = _mm256_cmpgt_epi16(out3,Qminusones);
    cmp4 = _mm256_cmpgt_epi16(out4,Qminusones);

    out1 = _mm256_sub_epi16(out1,_mm256_and_si256(Qs,cmp1));
    out2 = _mm256_sub_epi16(out2,_mm256_and_si256(Qs,cmp2));
    out3 = _mm256_sub_epi16(out3,_mm256_and_si256(Qs,cmp3));
    out4 = _mm256_sub_epi16(out4,_mm256_and_si256(Qs,cmp4));

    _mm256_store_si256((__m256i *) out->data    , out1);
    _mm256_store_si256((__m256i *) out->data + 1, out2);
    _mm256_store_si256((__m256i *) out->data + 2, out3);
    _mm256_store_si256((__m256i *) out->data + 3, out4);
}

void inline vec_mul_add(vec* out, const vec* a, uint16_t b){
	__m256i a1,a2,a3,a4,tmp1,tmp2,tmp3,tmp4,bb,out1,out2,out3,out4;

	bb = _mm256_set1_epi16(b);

    a1 = _mm256_load_si256((__m256i *)a->data);
    a2 = _mm256_load_si256((__m256i *)a->data + 1);
    a3 = _mm256_load_si256((__m256i *)a->data + 2);
    a4 = _mm256_load_si256((__m256i *)a->data + 3);

    out1 = _mm256_load_si256((__m256i *)out->data);
    out2 = _mm256_load_si256((__m256i *)out->data + 1);
    out3 = _mm256_load_si256((__m256i *)out->data + 2);
    out4 = _mm256_load_si256((__m256i *)out->data + 3);

    tmp1 = _mm256_mullo_epi16(a1,bb);
    tmp2 = _mm256_mullo_epi16(a2,bb);
    tmp3 = _mm256_mullo_epi16(a3,bb);
    tmp4 = _mm256_mullo_epi16(a4,bb);

    out1 = _mm256_add_epi16(out1,tmp1);
    out2 = _mm256_add_epi16(out2,tmp2);
    out3 = _mm256_add_epi16(out3,tmp3);
    out4 = _mm256_add_epi16(out4,tmp4);

    _mm256_store_si256((__m256i *) out->data, out1);
    _mm256_store_si256((__m256i *) out->data + 1, out2);
    _mm256_store_si256((__m256i *) out->data + 2, out3);
    _mm256_store_si256((__m256i *) out->data + 3, out4);
}

void inline vec_reduce(vec* v){
	__m256i a1,a2,a3,a4,b1,b2,b3,b4,cmp1,cmp2,cmp3,cmp4;

    a1 = _mm256_load_si256((__m256i *)v->data);
    a2 = _mm256_load_si256((__m256i *)v->data + 1);
    a3 = _mm256_load_si256((__m256i *)v->data + 2);
    a4 = _mm256_load_si256((__m256i *)v->data + 3);

#ifdef LEP
    cmp1 = _mm256_cmpgt_epi16(a1,multiples);
    cmp2 = _mm256_cmpgt_epi16(a2,multiples);
    cmp3 = _mm256_cmpgt_epi16(a3,multiples);
    cmp4 = _mm256_cmpgt_epi16(a4,multiples);

    a1 = _mm256_sub_epi16(a1,_mm256_and_si256(multiples,cmp1));
    a2 = _mm256_sub_epi16(a2,_mm256_and_si256(multiples,cmp2));
    a3 = _mm256_sub_epi16(a3,_mm256_and_si256(multiples,cmp3));
    a4 = _mm256_sub_epi16(a4,_mm256_and_si256(multiples,cmp4));
#endif

    b1 = _mm256_mulhi_epu16(a1,inverseofQ);
    b2 = _mm256_mulhi_epu16(a2,inverseofQ);
    b3 = _mm256_mulhi_epu16(a3,inverseofQ);
    b4 = _mm256_mulhi_epu16(a4,inverseofQ);

    b1 = _mm256_mullo_epi16(b1,Qs);
    b2 = _mm256_mullo_epi16(b2,Qs);
    b3 = _mm256_mullo_epi16(b3,Qs);
    b4 = _mm256_mullo_epi16(b4,Qs);

    a1 = _mm256_sub_epi16(a1,b1);
    a2 = _mm256_sub_epi16(a2,b2);
    a3 = _mm256_sub_epi16(a3,b3);
    a4 = _mm256_sub_epi16(a4,b4);

    _mm256_store_si256((__m256i *) v->data,     a1);
    _mm256_store_si256((__m256i *) v->data + 1, a2);
    _mm256_store_si256((__m256i *) v->data + 2, a3);
    _mm256_store_si256((__m256i *) v->data + 3, a4);
}

unsigned char inline vec_hw(const vec* v){
	__m256i a1,a2,a3,a4,b,c;

	a1 = _mm256_load_si256((__m256i *)v->data);
    a2 = _mm256_load_si256((__m256i *)v->data + 1);
    a3 = _mm256_load_si256((__m256i *)v->data + 2);
    a4 = _mm256_load_si256((__m256i *)v->data + 3);

    b = _mm256_cmpgt_epi16(a1,zeros);
    c = _mm256_and_si256(b,ones);

    b = _mm256_cmpgt_epi16(a2,zeros);
    b = _mm256_and_si256(b,ones);
    c = _mm256_add_epi16(b,c);

    b = _mm256_cmpgt_epi16(a3,zeros);
    b = _mm256_and_si256(b,ones);
    c = _mm256_add_epi16(b,c);

    b = _mm256_cmpgt_epi16(a4,zeros);
    b = _mm256_and_si256(b,ones);
    c = _mm256_add_epi16(b,c);

    c = _mm256_hadd_epi16(c,zeros);
    c = _mm256_hadd_epi16(c,zeros);
    c = _mm256_hadd_epi16(c,zeros);

	return _mm256_extract_epi16(c,0) + _mm256_extract_epi16(c,8);
}

// size of support of span of 2 vectors
unsigned char inline space_hw(const vec* v, const vec* w){
	__m256i a1,a2,a3,a4,b,c;

	a1 = _mm256_load_si256((__m256i *)v->data);
    a2 = _mm256_load_si256((__m256i *)v->data + 1);
    a3 = _mm256_load_si256((__m256i *)v->data + 2);
    a4 = _mm256_load_si256((__m256i *)v->data + 3);

    a1 = _mm256_or_si256(a1,_mm256_load_si256((__m256i *)w->data));
    a2 = _mm256_or_si256(a2,_mm256_load_si256((__m256i *)w->data+1));
    a3 = _mm256_or_si256(a3,_mm256_load_si256((__m256i *)w->data+2));
    a4 = _mm256_or_si256(a4,_mm256_load_si256((__m256i *)w->data+3));

    b = _mm256_cmpgt_epi16(a1,zeros);
    c = _mm256_and_si256(b,ones);

    b = _mm256_cmpgt_epi16(a2,zeros);
    b = _mm256_and_si256(b,ones);
    c = _mm256_add_epi16(b,c);

    b = _mm256_cmpgt_epi16(a3,zeros);
    b = _mm256_and_si256(b,ones);
    c = _mm256_add_epi16(b,c);

    b = _mm256_cmpgt_epi16(a4,zeros);
    b = _mm256_and_si256(b,ones);
    c = _mm256_add_epi16(b,c);

    c = _mm256_hadd_epi16(c,zeros);
    c = _mm256_hadd_epi16(c,zeros);
    c = _mm256_hadd_epi16(c,zeros);

	return _mm256_extract_epi16(c,0) + _mm256_extract_epi16(c,8);
}

#endif