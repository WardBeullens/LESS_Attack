#include <stdio.h>
#include "common.h"
#include "params.h"
#include <thread>

#define SPACES 2800000

void generate_sparse_spaces(vec* Min, int weight, vec* out, int spaces){
	int found = 0;

	vec M[N-K];
	memcpy(M,Min,sizeof(vec[N-K]));

	uintptr_t outint = (uintptr_t) out;

	mt19937 mt((unsigned int) outint);

	int count = 0;

	unsigned char I[N];
	while(1){
		random_gaussian_elimination(M,I,N-K, &mt);
		count ++;
		for (int i = 0; i < N-K; ++i)
		{
			if(vec_hw(M+i) > weight-1){
				continue;
			}

			for (int j = i+1; j < N-K; ++j)
			{
				if(space_hw(M+i,M+j) <= weight){
					memcpy(out + found*3  , M+i, sizeof(vec));
					memcpy(out + found*3+1, M+j, sizeof(vec));
					found ++;
					if (found == spaces){
						printf("number of gaussian eliminations: %d \n", count);
						return;
					}
				}
			}
		}

	}
}

int list_sparsests(vec* space, vec* sparsests){
	int number_of_sparsests = 1;
	int weight = vec_hw(space);

	memcpy(sparsests, space, sizeof(vec));
	for (int i = 0; i < Q; ++i)
	{
		vec_mul_add(space+1, space, 1);
		vec_reduce(space+1);
		int w = vec_hw(space+1);
		if (w<weight){
			weight = w;
			memcpy(sparsests,space+1,sizeof(vec));
			number_of_sparsests = 1;
		}
		else if(w == weight){
			memcpy(sparsests+number_of_sparsests,space+1,sizeof(vec));
			number_of_sparsests ++;
		}
	}
	return number_of_sparsests;
}

void swap_cols(vec* space, int i, int j){
	uint16_t temp = space[0].data[i];
	space[0].data[i] = space[0].data[j];
	space[0].data[j] = temp;

	temp = space[1].data[i];
	space[1].data[i] = space[1].data[j];
	space[1].data[j] = temp;
}

void minimize_bottom_right(vec* t, vec*out, int start){
	uint16_t count[Q] = {0};

	int highest_count = 0;
	int second_highest_count = 0;

	for (int i = start; i < N; ++i)
	{
		count[t->data[i]] ++;
	}

	for (int i = 0; i < Q; ++i)
	{
		if(count[i] >= highest_count){
			second_highest_count = highest_count;
			highest_count = count[i];
		}
		else if (count[i] >= second_highest_count){
			second_highest_count = count[i];
		}
	}

	vec temp = {0};
	if(highest_count < out->data[0])
		return;
	if(highest_count == out->data[0] && second_highest_count < out->data[1])
		return;

	for (int i = 0; i < Q; ++i)
	{
		if(count[i] != highest_count)
			continue;

		for (int j = 0; j < Q; ++j)
		{
			if(i == j || count[j] != second_highest_count)
				continue;

			int better = 0;

			for (int k = 0; k < Q; ++k)
			{
				temp.data[k] = count[(i + k*(Q+j-i))%Q];

				if(better == 0 && temp.data[k] < out->data[k]){
					break;
				}
				if(temp.data[k] > out->data[k]){
					better = 1;
				}
			}

			if( better > 0){
				memcpy(out, &temp, sizeof(vec));
			}
		}
	}
}

void normal_form(vec* space, vec* normal){
	vec sparsests[Q+1] = {0};

	normal->data[0] = 0;

	int number_of_sparsests = list_sparsests(space,sparsests);

	vec temp[2];
	for (int i = 0; i < number_of_sparsests; ++i)
	{

		// make basis with sparsest vector in top
		memcpy(temp,sparsests + i,sizeof(vec));
		if(memcmp(temp,space,sizeof(vec)) == 0){
			memcpy(temp+1,space+1,sizeof(vec));
		}
		else{
			memcpy(temp+1,space,sizeof(vec));
		}

		// make first nonzero entry of each column equal to 1
		for (int j = 0; j < N; ++j)
		{
			if(temp[0].data[j] != 0){
				temp[1].data[j] *= inverse[temp[0].data[j]];
				temp[1].data[j] %= Q;
				temp[0].data[j] = 1;
			}
			else if(temp[1].data[j] != 0){
				temp[1].data[j] = 1;
			}
		}

		// swap all the 1's to the back
		int start = 0;
		int end = N-1;
		while(start<end){
			if(temp[0].data[start] == 0){
				start ++;
			}
			else if(temp[0].data[end] == 1){
				end --;
			}
			else{
				swap_cols(temp,start,end);
				start ++;
			}
		}

		// minimize the bottom right part
		minimize_bottom_right(temp+1,normal,start);
	}
}

int compare (const void * a, const void * b)
{
	return memcmp((unsigned char*) a + sizeof(vec[2]), (unsigned char *) b+sizeof(vec[2]), sizeof(vec));
}

int is_zero_col(vec* a, int col){
	return !(a[0].data[col] || a[1].data[col]);
}

void update_perm(vec* a, vec* b, uint16_t* perm){
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			if( is_zero_col(a,i) != (is_zero_col(b,j)) ){
				perm[i*N+j] = 1;
			}
		}
	}
}

void normal_form_worker(vec *x, int len)
{
	for (int i = 0; i < len; ++i)
	{
		normal_form(x+(i*3),x+2+(i*3));
	}
}

int main(int argc, char const *argv[])
{
	time_t current_time;
	current_time = time(NULL);

	srand (time(NULL));

	vec M1[K];
	vec M2[K];
	vec permutation;
	vec diagonal;

	read_codes(M1,M2,&permutation,&diagonal);

	vec *sparse1 = (vec*) aligned_alloc(32,sizeof(vec[SPACES*3]));
	vec *sparse2 = (vec*) aligned_alloc(32,sizeof(vec[SPACES*3]));

	#define W 28

	printf("computing sparse spaces... \n");
	thread th1(generate_sparse_spaces, M1, W, sparse1, SPACES/2);
	thread th2(generate_sparse_spaces, M1, W, sparse1 + SPACES/2*3, SPACES/2);
	thread th3(generate_sparse_spaces, M2, W, sparse2, SPACES/2);
	thread th4(generate_sparse_spaces, M2, W, sparse2 + SPACES/2*3, SPACES/2);

	th1.join();
	th2.join();
	th3.join();
	th4.join();

	printf("ISD seconds: %ld\n", (long int) (time(NULL)-current_time));

	printf("computing normal forms...\n");
	th1 = thread(normal_form_worker, sparse1, SPACES/2);
	th2 = thread(normal_form_worker, sparse1+(SPACES/2*3), SPACES/2);
	th3 = thread(normal_form_worker, sparse2, SPACES/2);
	th4 = thread(normal_form_worker, sparse2+(SPACES/2*3), SPACES/2);
	
	th1.join();
	th2.join();
	th3.join();
	th4.join();

	// sort spaces
	qsort(sparse1,SPACES,sizeof(vec[3]),compare);
	qsort(sparse2,SPACES,sizeof(vec[3]),compare);

	int pos1 = 0;
	int pos2 = 0;

	uint16_t perm[N*N] = {0};

	int collisions = 0;

	// find collisions by walking through sorted lists
	while (pos1 < SPACES && pos2 < SPACES){
		int x = compare(sparse1+pos1*3, sparse2+pos2*3);
		if(x == 0){ // for each collision constrain the possible permutations 
			update_perm(sparse1+(pos1*3),sparse2+(pos2*3),perm);
			collisions ++;
			pos2 ++;
		}
		else if( x>0 ){
			pos2 ++;
		}
		else{
			pos1 ++;
		}
	}

	printf("collisions: %d \n",collisions);

	printf("real diagonal: \n");
	print_vec(&diagonal);
	printf("\n");
	printf("found permutation: \n");

	// check which possiblities are left
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
		if (x >= 0) // this position of the permutation is known
			printf("%2d ", x);
		if (x == -1) // there does not exists a value for this position of the permutation that satisfies all the constrainst -> something went wrong!
			printf(" / ");
		if (x == -2) // this position of the permutation is not known
			printf("   ");
	}

	printf("\n");

	printf("real permutation:  \n");
	print_vec(&permutation);

	printf("total seconds: %ld\n", (long int) (time(NULL)-current_time));

	return 0;
}