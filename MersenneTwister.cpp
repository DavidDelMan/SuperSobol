#include "MersenneTwister.h"

/* init static member */
int MersenneTwister::mti = NN+1;
unsigned long long MersenneTwister::mt[NN] = {};

MersenneTwister::MersenneTwister() {
  unsigned long long init[4] = {0x12345ULL, 0x23456ULL, 0x34567ULL,
				0x45678ULL};
  unsigned long long length = 4;
  init_by_array64(init, length);
}

/* initializes mt[NN] with a seed */
void MersenneTwister::init_genrand64(unsigned long long seed)
{
  mt[0] = seed;
  for (mti=1; mti<NN; mti++) 
    mt[mti] =  (6364136223846793005ULL * (mt[mti-1] ^ (mt[mti-1] >> 62)) + mti);
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
void MersenneTwister::init_by_array64(unsigned long long init_key[],
		     unsigned long long key_length)
{
  unsigned long long i, j, k;
  init_genrand64(19650218ULL);
  i=1; j=0;
  k = (NN>key_length ? NN : key_length);
  for (; k; k--) {
    mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 62)) * 3935559000370003845ULL))
      + init_key[j] + j; /* non linear */
    i++; j++;
    if (i>=NN) { mt[0] = mt[NN-1]; i=1; }
    if (j>=key_length) j=0;
  }
  for (k=NN-1; k; k--) {
    mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 62)) * 2862933555777941757ULL))
      - i; /* non linear */
    i++;
    if (i>=NN) { mt[0] = mt[NN-1]; i=1; }
  }

  mt[0] = 1ULL << 63; /* MSB is 1; assuring non-zero initial array */ 
}

/* generates a random number on [0, 2^64-1]-interval */
/* Used in genrand64_real3() */
unsigned long long MersenneTwister::genrand64_int64(void)
{
  int i;
  unsigned long long x;
  static unsigned long long mag01[2]={0ULL, MATRIX_A};

  if (mti >= NN) { /* generate NN words at one time */

    /* if init_genrand64() has not been called, */

    if (mti == NN+1)
      /* time-based initial seed is used */
      //   init_genrand64(time(NULL));
      /* a default initial seed is used     */
      init_genrand64(5489ULL); 

    for (i=0;i<NN-MM;i++) {
      x = (mt[i]&UM)|(mt[i+1]&LM);
      mt[i] = mt[i+MM] ^ (x>>1) ^ mag01[(int)(x&1ULL)];
    }
    for (;i<NN-1;i++) {
      x = (mt[i]&UM)|(mt[i+1]&LM);
      mt[i] = mt[i+(MM-NN)] ^ (x>>1) ^ mag01[(int)(x&1ULL)];
    }
    x = (mt[NN-1]&UM)|(mt[0]&LM);
    mt[NN-1] = mt[MM-1] ^ (x>>1) ^ mag01[(int)(x&1ULL)];

    mti = 0;
  }
  
  x = mt[mti++];

  x ^= (x >> 29) & 0x5555555555555555ULL;
  x ^= (x << 17) & 0x71D67FFFEDA60000ULL;
  x ^= (x << 37) & 0xFFF7EEE000000000ULL;
  x ^= (x >> 43);

  return x;
}

/* generates a random number on (0,1)-real-interval */
double MersenneTwister::genrand64_real3(void)
{
  return ((genrand64_int64() >> 12) + 0.5) * (1.0/4503599627370496.0);
}
