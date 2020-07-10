#ifndef PARAMS_H
#define PARAMS_H

#include <time.h>
#include <fstream>
#include <cstring>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <string>
#include <random>
#include <cstdint>

using namespace std;

#ifdef PEP
	#define N 60
	#define K 25

	#define Q 31
	#define Qinv 2115 

	static const uint16_t inverse[Q] = {0, 1, 16, 21, 8, 25, 26, 9, 4, 7, 28, 17, 13, 12, 20, 29, 2, 11, 19, 18, 14, 3, 24, 27, 22, 5, 6, 23, 10, 15, 30}; 
#endif

#ifdef LEP
	#define N 54
	#define K 27

	#define Q 53
	#define Qinv 1237

	static const uint16_t inverse[Q] = {0, 1, 27, 18, 40, 32, 9, 38, 20, 6, 16, 29, 31, 49, 19, 46, 10, 25, 3, 14, 8, 48, 41, 30, 42, 17, 51, 2, 36, 11, 23, 12, 5, 45, 39, 50, 28, 43, 7, 34, 4, 22, 24, 37, 47, 33, 15, 44, 21, 13, 35, 26, 52}; 
#endif

#endif