
/* Class InverseTransformation implements the inverse transformation
 * method of drawing numbers from a distribution.  In computing Sobol'
 * indices, this class is used to draw model parameters from their
 * respective distributions.
 */

#ifndef INVERSETRANSFORMATION_H
#define INVERSETRANSFORMATION_H

#include <cmath>
#include <vector>
#include <algorithm>
#include "MersenneTwister.h"

typedef double Type;
class InverseTransformation
{
 private:
  MersenneTwister MT;
  /* values needed in BSM approximation of inverse normal CDF */
  static constexpr Type a0 = 2.50662823884;
  static constexpr Type a1 = -18.61500062529;
  static constexpr Type a2 = 41.39119773534;
  static constexpr Type a3 = -25.44106049637;

  static constexpr Type b0 = -8.4735109309;
  static constexpr Type b1 = 23.08336743743;
  static constexpr Type b2 = -21.06224101826;
  static constexpr Type b3 = 3.13082909833;

  static constexpr Type c0 = 0.3374754822726147;
  static constexpr Type c1 = 0.9761690190917186;
  static constexpr Type c2 = 0.1607979714918209;
  static constexpr Type c3 = 0.0276438810333863;
  static constexpr Type c4 = 0.0038405729373609;
  static constexpr Type c5 = 0.0003951896511919;
  static constexpr Type c6 = 0.0000321767881768;
  static constexpr Type c7 = 0.0000002888167364;
  static constexpr Type c8 = 0.0000003960315187;

 public:
  InverseTransformation();
  Type GenPareto(Type k, Type sigma, Type theta);
  Type Normal(Type u, Type mean, Type variance);
  Type Uniform(Type u, Type a, Type b);
  Type AndersonDarlingNormal(std::vector<Type> values, 
			     Type mean,
			     Type variance);
  Type NormCDF(Type x);
};
#endif
