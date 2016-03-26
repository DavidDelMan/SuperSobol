
#ifndef MERSENNETWISTER_H
#define MERSENNETWISTER_H

#define NN 312
#define MM 156
#define MATRIX_A 0xB5026F5AA96619E9ULL
#define UM 0xFFFFFFFF80000000ULL /* Most significant 33 bits */
#define LM 0x7FFFFFFFULL /* Least significant 31 bits */

#include <cstdlib>
#include <iostream>
#include <ctime>  // time seed to init_genrand64()

class MersenneTwister {
 private:
  /* The array for the state vector */
  static unsigned long long mt[NN]; 
  /* mti==NN+1 means mt[NN] is not initialized */
  static int mti;
  void init_genrand64(unsigned long long seed);
  void init_by_array64(unsigned long long init_key[],
		       unsigned long long key_length);
  unsigned long long genrand64_int64(void);

 public:
  MersenneTwister();
  double genrand64_real3(void);
};
#endif

