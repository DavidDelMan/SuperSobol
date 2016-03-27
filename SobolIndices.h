/* Contains the functions necessary to compute Sobol sensitivity
 * indices for a real-valued function of a vector variable.
 * Using Owen's impementation. */

#ifndef SOBOLINDICES_H
#define SOBOLINDICES_H

#define _USE_MATH_DEFINES  /* M_PI for math constant pi */

#include <iostream>
#include <vector>
#include <set>
#include <fstream>
#include "Halton.h"
#include "MT64.h"
#include "InverseTransformation.h"

typedef double Type;


class SobolIndices
{
 private:
  Type (*model)(const std::vector<Type>&,
		const std::vector<Type>&);  /* model */
  int dim;  /* number of model parameters */
  unsigned int N_MC;  /* no. of MC runs to use */
  Type CoV;  /* coefficient of variation = std/mean */

  /* Sobol indices */
  Type lowerIndex, totalIndex, modelVariance, modelMean;
  std::vector<Type> x1, x2, arg1, arg2;  /* model args */
  std::vector<Type> constants;  /* model constants: K,r,... */
  std::set<int> indices;  /* index set to compute Sobol indices for */

  /* distribution params of model params */
  std::vector<std::vector<Type> > distroParams;
  halton *randomNumberGenerator;  /* halton (RASRAP) object */
  InverseTransformation *invTrans; /* inverse tarsnformation object */

 public:
  SobolIndices(Type (*model_)(const std::vector<Type>&,
			      const std::vector<Type>&),
	       const std::vector<Type> &constants_,
	       const std::set<int> &indices_,
	       const std::vector<std::vector<Type> >
	       &initialDistroParams_,
	       int dim_,
	       unsigned int N_MC_,
	       Type CoV_ = 1.0);
  void DisplayMembers();
  Type ComputeSensitivityIndices(const std::vector<Type> 
				 &uncertainties,
				 const std::set<int> &indices_
				 = std::set<int>());
  void AssignModelArguments(const std::set<int>& indices_);
  void TransformToModelDomain(const std::vector<Type> &uncertainties
			      = std::vector<Type>());
  std::vector<std::vector<Type> >
    PlotCoV(const std::vector<Type> &CoV_Vector, 
	    std::string &filename);

  void DisplayVector(const std::vector<Type>& vec);
  void DisplaySet(const std::set<int>& s);
  void DisplayVector(const std::vector<std::vector<Type> >& vec);
  Type GetLowerIndex() {return lowerIndex;}
  Type GetTotalIndex() {return totalIndex;}
  /* void SetDistroParams(const std::vector<std::vector<Type> >& */
  /* 		       distroParams_); */
  ~SobolIndices()
    {
      delete randomNumberGenerator;
      delete invTrans;
    }

};
#endif
