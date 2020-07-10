#include <stdio.h>
#include "common.h"
#include "params.h"
#include <thread>

void generate_sparse_codewords(vec* Min, int weight, vec* out, int vectors){
	int found = 0;
	unsigned char I[N];

	vec M[K];
	memcpy(M,Min,sizeof(vec[K]));

	uintptr_t outint = (uintptr_t) out;
	mt19937 mt((unsigned int) outint);

	long int GEs = 0;

	while(1){
		random_gaussian_elimination(M,I,K,&mt);
		GEs ++;
		for (int i = 0; i < K; ++i)
		{
			for (int j = i+1; j < K; ++j)
			{
				memcpy(out + 2*found, M+i, sizeof(vec));

				if(vec_hw(out+2*found) <= weight){
					found ++;
					if(found % 1000 == 0)
						printf("%d\n", found);
					if (found == vectors){
						printf("GEs: %ld \n", GEs);
						return;
					}
					break;
				}

				for (int k = 0; k < Q-1; ++k)
				{
					vec_add_reduce(out+2*found, out+2*found, M+j);

					if(vec_hw(out+2*found) <= weight){
						found ++;
						if(found % 10000 == 0)
							printf("%d\n", found);
						if (found == vectors){
							printf("GEs: %ld \n", GEs);
							return;
						}
						break;
					}

				}
			}
		}
	}
}

#define VECS 25000

#define DLOGS

void normal_form(vec* v){
	vec V2;
	memcpy(&V2,v,sizeof(vec));
	vec V = {0};
	memset(v+1,0,sizeof(vec));
	for (int i = 0; i < Q-1; ++i)
	{
		vec_add_reduce(&V,&V,&V2);
		uint16_t counts[Q] = {0};
		for (int j = 0; j < N; ++j)
		{
			counts[V.data[j]] ++;
		}
		if(memcmp(counts, v[1].data, sizeof(counts)) > 0){
			memcpy(v[1].data,counts,sizeof(counts));
			memcpy(v, &V, sizeof(vec));
		}
	}
}

int compare_vecs(const void* A , const void* B){
	vec *a = (vec*) A;
	vec *b = (vec*) B;

	int ans = 0;
	for (int i = 0; i < N; ++i)
	{
		ans += a->data[i] - b->data[i] + a->data[i]*a->data[i]*30 - b->data[i]*b->data[i]*30 + a->data[i]*a->data[i]*a->data[i]*10 - b->data[i]*b->data[i]*b->data[i]*10 ; 
	}
	return ans;
}

void update_perm(vec* a, vec* b, uint16_t* perm){
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			if(a->data[i] != b->data[j]){
				perm[i*N+j] ++;
			}
		}
	}
}

int compare (const void * a, const void * b)
{
	return memcmp((unsigned char*) a + sizeof(vec), (unsigned char *) b+sizeof(vec), sizeof(vec));
}

int main(int argc, char const *argv[])
{
	srand (time(NULL));

	vec M1[K];
	vec M2[K];
	vec solution;

	vec M1perp[N-K] = {0};
	vec M2perp[N-K] = {0};

	read_codes(M1,M2,&solution, NULL);

	vec *sparse1 = (vec *) aligned_alloc(32,sizeof(vec[2*VECS]));
	vec *sparse2 = (vec *) aligned_alloc(32,sizeof(vec[2*VECS]));

	#define W 30

	printf("computing sparse codewords... \n");
	thread th1(generate_sparse_codewords, M1, W, sparse1,  VECS/2);
	thread th2(generate_sparse_codewords, M1, W, sparse1 + VECS, VECS/2);
	thread th3(generate_sparse_codewords, M2, W, sparse2,  VECS/2);
	thread th4(generate_sparse_codewords, M2, W, sparse2 + VECS, VECS/2);

	th1.join();
	th2.join();
	th3.join();
	th4.join();

	printf("computing normal forms \n");

	for (int i = 0; i < VECS; ++i)
	{
		normal_form(sparse1 + i*2);
		normal_form(sparse2 + i*2);
	}

	printf("processing codewords...\n");

	qsort(sparse1, VECS, sizeof(vec[2]), compare);
	qsort(sparse2, VECS, sizeof(vec[2]), compare);

	int pos1 = 0;
	int pos2 = 0;

	uint16_t perm[N*N] = {0};

	int collisions = 0;
	int self_collisions1 = 0;
	int self_collisions2 = 0;

	while (pos1 < VECS && pos2 < VECS){
		if(pos1 < VECS-1 && compare(sparse1 + pos1*2, sparse1+2+pos1*2) == 0){
			if( memcmp(sparse1+pos1*2,sparse1+2+pos1*2,sizeof(vec)) != 0 ){
				printf("SELF COLLISION: \n");
				print_vec(sparse1+pos1*2);
				print_vec(sparse1+2+pos1*2);
				printf("\n");
				self_collisions1 ++;
			}
			pos1 ++;
		}
		if(pos2 < VECS-1 && compare(sparse2 + pos2*2, sparse2+2+pos2*2) == 0){
			if(memcmp(sparse2+pos2*2,sparse2+2+pos2*2,sizeof(vec)) != 0)
				self_collisions2 ++;
			pos2 ++;
		}
		int x = compare(sparse1+pos1*2, sparse2+pos2*2);
		if(x == 0){
			collisions ++;
			print_vec(sparse1+pos1*2);
			print_vec(sparse2+pos2*2);
			update_perm(sparse1+pos1*2,sparse2+pos2*2,perm);
			printf("\n");
			pos1 ++;
		}
		else if( x>0 ){
			pos2 ++;
		}
		else{
			pos1 ++;
		}
	}

	printf("collisions: %d \n", collisions);
	printf("self collisions 1: %d \n", self_collisions1);
	printf("self collisions 2: %d \n", self_collisions2);

	for (int i = 0; i < N; ++i)
	{
		int x = -1;
		for (int j = 0; j < N; ++j)
		{
			if (perm[j*N+i] == 0){
				if( x == -1 ){
					x = j;
				}
				else{
					x = -2;
				}
			}
		}
		if (x >= 0)
			printf("%2d ", x);
		if (x == -1)
			printf(" / ");
		if (x == -2)
			printf("   ");
	}

	printf("\n");
	for (int i = 0; i < N; ++i)
	{
		printf("%2d ", solution.data[i]);
	}
	printf("\n");

	return 0;
}