#include "SuperSobolIndices.h"

/* Ctor
 * Input:
 *
 * model_ = Type-valued function of a vector of Types.  First arg is
 *   vector of parameters to be drawn randomly.  Second arg is vector
 *   of fixed constants, like strike & interest rate.
 * constants_ = vector of constants for model, like strike price
 * indices_ = set of parameter indices to compute SIs for
 * paramUncertaintyDistroParams_ = "hyperparameters" giving the 
 *    parameters for the model parameters' uncertainties distros
 * dim_ = number of model parameters
 * N_MC_ = number of Monte Carlo runs to compute Sobol indices
 * N_Super_Sobol_ = number of Monte Carlo runs to compute Super Sobol
 *    indices
 */
SuperSobolIndices::
SuperSobolIndices(Type (*model_)(const std::vector<Type>&, 
				 const std::vector<Type>&),
		  const std::vector<Type> &constants_,
		  const std::set<int> &indices_,
		  const std::vector<std::vector<Type> >
		  &initialDistroParams_,
		  const std::vector<std::vector<Type> > 
		  &paramUncertaintyDistroParams_,
		  const unsigned int dim_,
		  const unsigned int N_MC_,
		  const unsigned int N_Super_Sobol_)
{
  // model = model_;
  // constants = constants_;
  indices = indices_;
  // initalDistroParams = initialDistroParams_;
  paramUncertaintyDistroParams = paramUncertaintyDistroParams_;
  dim = dim_;
  N_Super_Sobol = N_Super_Sobol_;

  // intialize Super Sobol indices
  lowerSuperIndex = 0;
  totalSuperIndex = 0;

  // allocate model argument vectors
  s1.resize(dim);
  s2.resize(dim);
  s_arg1.resize(dim);
  s_arg2.resize(dim);

  // construct halton (RASRAP) & InverseTransformation objects
  RNG = new halton();
  invTrans = new InverseTransformation();

  // construct SobolIndices object
  sobol = new SobolIndices(model_, constants_, indices, 
			   initialDistroParams_, dim, N_MC_);

  // init RNG: length of Halton vector, random start, random permute
  RNG->init(2*dim,true,true);
}

/* Displays member variables of the SobolIndices class */
void SuperSobolIndices::DisplayMembers()
{
  std::cout << "Members of SuperSobolIndices: \n\n";
  std::cout << "dim: " << dim << "\n";
  std::cout << "N_Super_Sobol: " << N_Super_Sobol << "\n";
  std::cout << "lowerSuperIndex: " << lowerSuperIndex << "\n";
  std::cout << "totalSuperIndex: " << totalSuperIndex << "\n";
  std::cout << "superModelVariance: " << superModelVariance << "\n";
  std::cout << "superModelMean: " << superModelMean << "\n";
  std::cout << "indices: \n";
  DisplaySet(indices);
  std::cout << "paramUncertaintyDistroParams: \n";
  DisplayVector(paramUncertaintyDistroParams);

  // std::cout << "randomNumbers: \n";
  // DisplayVector(randomNumbers);
  // std::cout << "x: \n";
  // DisplayVector(x);
  // std::cout << "arg1: \n";
  // DisplayVector(arg1);
  // std::cout << "arg2: \n";
  // DisplayVector(arg2);
  std::cout << "\n";
}


void SuperSobolIndices::
ComputeSuperSobolIndices()
{
  // MC accumulators
  Type f0_sum_super = 0, Dy_sum_super = 0, DT_sum_super = 0, 
    D_sum_super = 0;

  // model evaluations
  Type F, F2, F_model1, F_model2;

  for (unsigned int i = 0; i < N_Super_Sobol; ++i)
    {
      // std::cout << i << "\n";
      // generate 2*dim random numbers
      RNG->genHalton();

      // transform each random number to parameter uncertainty distro
      TransformToParamUncertaintyDomain();

      /* assign xformed RVs to proper model argument vectors, will now
       * have uncertainties for each parameter */
      AssignUncertaintyModelArguments();

      // compute Sobol index for given uncertainties
      F = sobol->ComputeSensitivityIndices(s1);
      F2 = sobol->ComputeSensitivityIndices(s2);
      F_model1 = sobol->ComputeSensitivityIndices(s_arg1);
      F_model2 = sobol->ComputeSensitivityIndices(s_arg2);

      // MC accumulations for Super Sobol indices
      f0_sum_super += F;
      D_sum_super += F*F;
      Dy_sum_super += F*(F_model1 - F2); 
      DT_sum_super += pow((F - F_model2), 2.0);
    }

  // compute Super Sobol indices
  superModelMean = f0_sum_super/N_Super_Sobol;
  superModelVariance = D_sum_super/N_Super_Sobol 
    - superModelMean*superModelMean;

  Type Dy_super = Dy_sum_super/N_Super_Sobol;
  Type DT_super = DT_sum_super/N_Super_Sobol;

  std::cout << "Dy_super = " << Dy_super << "\n";
  std::cout << "DT_super = " << DT_super << "\n";

  //  /* normalized */
  // lowerIndex = Dy_super/superModelVariance;
  // totalIndex = DT_super/(2.0*superModelVariance);

  /* non-normalized */
  lowerSuperIndex = Dy_super;
  totalSuperIndex = DT_super/2.0;
}

/* Fills the s_arg1 and s_arg2 member vectors that hold the 
 * uncertainties for the corresponding parameters according to the
 * parameter index for which we are computing Super Sobol indices for
 */
void SuperSobolIndices::
AssignUncertaintyModelArguments()
{
  for (int j = 0; j < dim; ++j)
    {
      // check if "j" is in index set to compute Super Sobol index for
      bool inIndexSet = indices.count(j+1);

      if (inIndexSet)
	{
	  s_arg1[j] = s1[j];
	  s_arg2[j] = s2[j];
	}
      else
	{
	  s_arg1[j] = s2[j];
	  s_arg2[j] = s1[j];
	}
    }
}

/* Retrieves the Unif(0,1) random numbers generates in 
 * ComputeSuperSobolIndices() and transforms them to the respective
 * distributions of the uncertainties, s1 and s2.
 */
void SuperSobolIndices::
TransformToParamUncertaintyDomain()
{
  // retrieve random numbers computed in ComputeSuperSobolIndices()
  for (int j = 0; j < dim; ++j)
    {
      s1[j] = RNG->get_rnd(j+1);
      s2[j] = RNG->get_rnd(j+1+dim);

      // left and right endpoints of uniform uncertainty distros
      Type a = paramUncertaintyDistroParams[j][0];
      Type b = paramUncertaintyDistroParams[j][1];

      // generate Unif(a,b) RVs from uncertainties
      s1[j] = invTrans->Uniform(s1[j],a,b);
      s2[j] = invTrans->Uniform(s2[j],a,b);
    }
}

// /* Changes the variance of each parameter by updating the second
//  * coordinate of member vector "distroParams[][]."
//  */
// void SuperSobolIndices::ChangeParameterUncertainty()
// {
//   for (int i = 0; i < dim; ++i)
//     {
//       distroParams[i][1] = s[i];
//     }

//   // vfy distroParams filled properly
//   DisplayVector(distroParams);
// }


void SuperSobolIndices::
DisplayVector(const std::vector<std::vector<Type> > &vec)
{
  for (const auto& i : vec)
    {
      for (const auto& j : i)
	{
	  std::cout << j << " ";
	}
      std::cout << "\n";
    }
  std::cout << "\n";
}

void SuperSobolIndices::DisplaySet(const std::set<int> &s)
{
  for (auto i : s)
    {
      std::cout << i << " ";
    }
  std::cout << "\n\n";
}
