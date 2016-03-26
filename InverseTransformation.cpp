#include "InverseTransformation.h"

/* Dflt ctor */
InverseTransformation::InverseTransformation() : MT() {}

/* Returns a sample from the Generalized Pareto distribution.  This
 * function is different from the others in that the pseudorandom
 * number is generated within the function, as opposed to passed in.
 *
 * Input:
 * k = shape parameter
 * sigma = scale parameter, > 0
 * theta = location parameter
 */
Type InverseTransformation::GenPareto(Type k, Type sigma, Type theta)
{
  Type u = MT.genrand64_real3();  // generate psuedorandom number
  Type result = theta + sigma/k * (std::pow(1.0-u,-k) - 1);

  if (result > 100)
    {
      std::cout << "big theta: \n";
      std::cout << "k = " << k << "\n";
      std::cout << "sigma = " << sigma << "\n";
      std::cout << "theta = " << theta << "\n";
      std::cout << "u = " << u << "\n";
      std::cout << "sigma/k = " << sigma/k << "\n";
      std::cout << "1.0-u = " << 1.0-u << "\n";
      std::cout << "pow(1.0-u,-k) = " << std::pow(1.0-u,-k) << "\n";
      std::cout << "pow(1.0-u,-k)-1 = " << std::pow(1.0-u,-k)-1 
		<< "\n";
      std::cout << "sigma/k * (pow(1.0-u,-k)-1) = " 
		<< sigma/k * (std::pow(1.0-u,-k)-1) << "\n";
      std::cout << "theta + sigma/k * (pow(1.0-u,-k)-1) = " 
		<< theta + sigma/k * (std::pow(1.0-u,-k)-1) << "\n";
      std::cout << "result = " << result << "\n\n";
    }

  if (result < theta)
    {
      std::cout << "small theta: \n";
      std::cout << "k = " << k << "\n";
      std::cout << "sigma = " << sigma << "\n";
      std::cout << "theta = " << theta << "\n";
      std::cout << "u = " << u << "\n";
      std::cout << "sigma/k = " << sigma/k << "\n";
      std::cout << "1.0-u = " << 1.0-u << "\n";
      std::cout << "pow(1.0-u,-k) = " << std::pow(1.0-u,-k) << "\n";
      std::cout << "pow(1.0-u,-k)-1 = " << std::pow(1.0-u,-k)-1 
		<< "\n";
      std::cout << "sigma/k * (pow(1.0-u,-k)-1) = " 
		<< sigma/k * (std::pow(1.0-u,-k)-1) << "\n";
      std::cout << "theta + sigma/k * (pow(1.0-u,-k)-1) = " 
		<< theta + sigma/k * (std::pow(1.0-u,-k)-1) << "\n";
      std::cout << "result = " << result << "\n\n";
    }

  return result;
}

/* Returns a normally-distributed pseudo-/quasi- random number.
 * Input:
 * 
 * u = Unif(0,1) random number
 * mean = mean of normal distro
 * variance = variance of normal distro
 */
Type InverseTransformation::Normal(Type u, Type mean, Type variance)
{
  Type y = u - 0.5;
  Type x;  // return value

  if (std::abs(y) < 0.42) {
    Type r = y*y;
    x = y*(((a3*r + a2)*r + a1)*r + a0) / 
      ((((b3*r + b2)*r + b1)*r + b0)*r + 1);
  }

  else {
    Type r = u;
    if (y > 0)
      r = 1 - u;
    r = log(-log(r));
    x = c0 + r*(c1 + r*(c2 + r*(c3 + r*(c4 + r*(c5 + r*(c6 + r*(c7 
							 + r*c8)))))));
    if (y < 0)
      x = -x;
  }

  return mean + sqrt(variance)*x;
}

/* Function Uniform transforms a Unif(0,1) random number to a
 * Unif(a,b) random number
 */
Type InverseTransformation::Uniform(Type u, Type a, Type b)
{
  return (b-a)*u + a;
}

/* Function AndersonDarlingNormal computes the Anderson Darling test
 * statistic for a standard normal distribution.  The vector "values"
 * is sorted in this function.  This function
 * normalizes the data to N(0,1) before computing the Anderson
 * Darling statistic for a standard normal, A^2.
 */
Type InverseTransformation::
AndersonDarlingNormal(std::vector<Type> values, 
		      Type mean,
		      Type variance)
{

  /* sort the values in increasing order */
  std::sort(values.begin(), values.end());

  /* number of random numbers of values vector */
  unsigned int N = values.size();

  /* normalize the data in the output member vector */
  for (int i = 0; i < N; ++i) {
    values[i] = (values[i] - mean) / sqrt(variance);
  }

  /* print sorted normalized data for debugging */
  // std::cout << "output in ADN after normalizing:\n";
  //   for (auto i : values)
  //     std::cout << i << " ";
  // std::cout << "\n";

  /* compute the standard normal cdf values of the normalized data */
  std::vector<Type> CDFValues;

  for (auto i : values) {
    Type val = NormCDF(i);  // NormCDF = standard normal CDF
    CDFValues.push_back(val);
  }

  /* print CDFValues of normalized data for debugging */
  // std::cout << "verifying CDFValues are sorted in ADN:\n";
  // for (auto i : CDFValues)
  //   std::cout << i << " ";
  // std::cout << "\n";

  /* compute sum in the Anderson-Darling statistic */
  Type sum = 0;
  //  Type sum2 = 0;
  for (int i = 1; i <= N; ++i) {
    sum += (2*i - 1) * ( log(CDFValues[i-1]) + log(1-CDFValues[N-i]) );
    // sum2 += (2*i-1)*log(CDFValues[i-1]) + 
    //   (2*N + 1 - 2*i)*log(1-CDFValues[i-1]);
  }

  Type ASquared =  -(1.0/N)*sum - N;
  //  Type ASquared2 = -(1.0/N)*sum2 - N;

  //  std::cout << "ASquared2: " << ASquared2 << "\n";

  return ASquared;
}

/* Approximates the standard normal cdf at x */
Type InverseTransformation::NormCDF(Type x)
{
  // constants
  Type a1 =  0.254829592;
  Type a2 = -0.284496736;
  Type a3 =  1.421413741;
  Type a4 = -1.453152027;
  Type a5 =  1.061405429;
  Type p  =  0.3275911;
 
  // Save the sign of x
  int sign = 1;
  if (x < 0)
    sign = -1;
  x = fabs(x)/sqrt(2.0);
 
  // A&S formula 7.1.26
  Type t = 1.0/(1.0 + p*x);
  Type y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);
 
  return 0.5*(1.0 + sign*y);
}
