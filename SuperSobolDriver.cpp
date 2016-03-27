#include "SuperSobolIndices.h"
#include <cmath>
#include <fstream>
#include <thread>  // std::this_thread::sleep_for
#include <chrono>  // std::chrono::seconds
#include <ctime>

/* Practice linear model, parameters.size() = 4 */
Type LinearModel(const std::vector<Type> &parameters,
		 const std::vector<Type> &constants)
{
  Type Y = 0;
  Type c = 1.0;
  for (int i = 0; i < 4; ++i)
    {
      Y += c*parameters[i];
    }
  return Y;
}

void DisplayVector(const std::vector<std::vector<Type> >& vec)
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

void DisplayVector(const std::vector<Type>& vec)
{
  for (auto& i : vec)
    {
      std::cout << i << " ";
    }
  std::cout << "\n\n";
}

void DisplaySet(const std::set<int>& s)
{
  for (auto i : s)
    {
      std::cout << i << " ";
    }
  std::cout << "\n\n";
}

int main(int argc, char** argv)
{
  /**** Simulation Parameters ****/
  int dim = 4;  // number of parameters in model
  std::vector<Type> constants = {};  // constant model parameters

  // /* verify constants vector is filled properly */
  // std::cout << "constants, in main: \n";
  // DisplayVector(constants);

  /* number of MC runs to use in Sobol indices approximation */
   int N_MC = 10000;
  //  int N_MC = atoi(argv[1]);

   /* number of MC runs to compute Super Sobol indices */
   int N_Super_Sobol = 10000;


  /* index set to compute Super Sobol index for */
  std::set<int> indices = {2};
  // std::cout << "indices: \n";
  // DisplaySet(indices);

  // /* coefficient of variation = std/mean for current parameter group */
  // Type CoV = 0.2;
  // std::vector<Type> CoV_Vector 
  //   = {0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35};

  // /* file name to plot indices to when using CoV routines */
  // std::string filename = "IndicesData.txt";

  /* specify INITIAL distribution parameters of model parameters */
  std::vector<std::vector<Type> > 
    initialDistroParams(dim, std::vector<Type>(2));

  initialDistroParams[0][0] = 0;
  initialDistroParams[0][1] = 1;
  initialDistroParams[1][0] = 0;
  initialDistroParams[1][1] = 4;
  initialDistroParams[2][0] = 0;
  initialDistroParams[2][1] = 9;
  initialDistroParams[3][0] = 0;
  initialDistroParams[3][1] = 16;

  /* Specify "hyperparameters" = distribution parameters of model
   * parameters' uncertainties.  For now, assuming the uncertainties
   * are Unif(alpha*sig_i, beta*sig_i) where
   *    sig_i = calibrated variance of parameter i,
   *    (0 < alpha < beta < 2) && (alpha + beta = 2).
   * The first constraint on alpha and beta ensure alpha acts to
   * decrease the variance and beta increases it.  The second
   * constraint is to ensure the mean of the Unif distro is sig_i.
*/
  Type alpha = 0.5;  // factor to decrease calibrated uncertainty by
  Type beta = 1.5;  // factor to increase calibrated uncertainty by
  std::vector<std::vector<Type> > 
    paramUncertaintyDistroParams(dim, std::vector<Type>(2));


  paramUncertaintyDistroParams[0][0] = alpha*initialDistroParams[0][1];
  paramUncertaintyDistroParams[0][1] = beta*initialDistroParams[0][1];
  paramUncertaintyDistroParams[1][0] = alpha*initialDistroParams[1][1];
  paramUncertaintyDistroParams[1][1] = beta*initialDistroParams[1][1];
  paramUncertaintyDistroParams[2][0] = alpha*initialDistroParams[2][1];
  paramUncertaintyDistroParams[2][1] = beta*initialDistroParams[2][1];
  paramUncertaintyDistroParams[3][0] = alpha*initialDistroParams[3][1];
  paramUncertaintyDistroParams[3][1] = beta*initialDistroParams[3][1];

  // verify paramUncertaintyDistroParams is correct
  std::cout << "paramUncertaintyDistroParams: \n";
  DisplayVector(paramUncertaintyDistroParams);


  /* construct a SobolIndices object */
  // SobolIndices sobol(Heston, constants, indices, initialDistroParams,
  // 		     dim, N_MC, CoV);
  // SobolIndices sobol(LinearModel,constants,indices, 
  // 		     initialDistroParams,dim,N_MC);

  // construct a SuperSobolIndices object
  SuperSobolIndices superSobol(LinearModel,constants,indices,
			       initialDistroParams,
			       paramUncertaintyDistroParams,dim,N_MC,
			       N_Super_Sobol);
  /* print member of SobolIndices object for verification */
  superSobol.DisplayMembers();

  clock_t tic = clock();

  // /* compute sensitivity indices */
  // std::cout << "computing sensitivity indices...\n\n";
  //  sobol.ComputeSensitivityIndices();
  // // std::vector<std::vector<Type> > results 
  // //   = sobol.PlotCoV(CoV_Vector, filename);

  // Compute Super Sobol indices
  std::cout << "Computing Super Sobol indices... \n\n";
  superSobol.ComputeSuperSobolIndices();

  std::cout << "...done.\n\n";

  Type toc = (Type)(clock() - tic) / CLOCKS_PER_SEC;
  std::cout << "total time: " << toc << "\n\n";

  /* display sensitivity indices */
  superSobol.DisplayMembers();

  // /* write to file */
  // std::ofstream File("sigma.txt", std::ios::app);
  // File << N_MC << " " << actual_N_MC << " " << sobol.GetLowerIndex() << " " << sobol.GetTotalIndex() << "\n";
  // File.close();
}
