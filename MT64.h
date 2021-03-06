/* This is a C++ program for a 64-bit version of Mersenne Twister pseudorandom number
   generator based on the C program available at
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/VERSIONS/C-LANG/mt19937-64.c
*/

#ifndef _MT64_H
#define _MT64_H

#include <iostream>

#define NN 312
#define MM 156
#define MATRIX_A_64 0xB5026F5AA96619E9ULL
#define UM 0xFFFFFFFF80000000ULL /* Most significant 33 bits */
#define LM 0x7FFFFFFFULL /* Least significant 31 bits */


class genRand_64{
public:
	static genRand_64* Instance(){if(_instance==NULL) _instance=new genRand_64();return _instance;}
	/* initializes mt[NN] with a seed */
	void init_genrand64(unsigned long long seed);

	/* initialize by an array with array-length */
	/* init_key is the array for initializing keys */
	/* key_length is its length */
	void init_by_array64(unsigned long long init_key[], 
				 unsigned long long key_length);

	/* generates a random number on [0, 2^64-1]-interval */
	unsigned long long genrand64_int64(void);


	/* generates a random number on [0, 2^63-1]-interval */
	long long genrand64_int63(void);

	/* generates a random number on [0,1]-real-interval */
	double genrand64_real1(void);

	/* generates a random number on [0,1)-real-interval */
	double genrand64_real2(void);

	/* generates a random number on (0,1)-real-interval */
	double genrand64_real3(void);
private:
	genRand_64()
	{
		init_by_array64(init_64, length_64);
	}
	static unsigned long long init_64[4], length_64;
	static genRand_64* _instance;
};



#endif