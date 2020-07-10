#include "common.h"
#include <stdio.h>
#include <cstdint>
#include <cstdlib>
#include <cstring>

void print_vec(const vec* V){
	for (int j = 0; j < N; ++j)
	{
		printf("%2d ",(int) V->data[j]);
	}
	printf("\n");
}

void print_matrix(const vec* matrix, unsigned char rows){
	for (int i = 0; i < rows; ++i)
	{
		print_vec(matrix+i);
	}
	printf("\n");
}

void random_gaussian_elimination(vec* matrix, unsigned char *I, unsigned char k, mt19937 *mt){
	unsigned char row = 0;

	while(row < k){
		I[row] = rand() % N;
		while(matrix[row].data[I[row]] == 0){
			I[row] = (*mt)() % N;
		}

		for (int i = 0; i < k; ++i)
		{
			if(i==row){
				continue;
			}
			
			uint16_t scalar = (inverse[matrix[row].data[I[row]]] * matrix[i].data[I[row]]) % Q;
			scalar = (Q - scalar)%Q;

			vec_mul_add(matrix + i, matrix + row, scalar);
			vec_reduce(matrix+i);
		}

		row += 1;

	}

}

void read_codes(vec* matrix1, vec* matrix2, vec* perm, vec* diagonal){

#ifdef PEP
	ifstream in("permutation_equivalent_codes.txt");
#endif

#ifdef LEP
	ifstream in("linearly_equivalent_codes.txt");
#endif

	string line, entry;
	
	// read M1
	for (int i = 0; i < K; ++i)
	{
		getline(in, line);
		istringstream iline(line);
		for (int j = 0; j < N; ++j)
		{
			getline(iline, entry, ' ');
			matrix1[i].data[j] = (uint16_t) stoi(entry);
		}
	}

	// read M2
	for (int i = 0; i < K; ++i)
	{
		getline(in, line);
		istringstream iline(line);
		for (int j = 0; j < N; ++j)
		{
			getline(iline, entry, ' ');
			matrix2[i].data[j] = (uint16_t) stoi(entry);
		}
	}

	// read permutation
	getline(in, line);
	istringstream iline(line);
	for (int j = 0; j < N; ++j)
	{
		getline(iline, entry, ' ');
		perm->data[j] = (uint16_t) stoi(entry);
	}

#ifdef LEP
	// read diagonal
	getline(in, line);
	istringstream iline2(line);
	for (int j = 0; j < N; ++j)
	{
		getline(iline2, entry, ' ');
		diagonal->data[j] = (uint16_t) stoi(entry);
	}
#endif

	in.close();
}