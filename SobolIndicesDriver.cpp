#include "SobolIndices.h"
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
  Type c = 0.1;
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
  /* model dimension, i.e., number of parameters in model */
  int dim = 4;

  /* specify constant model parameters */
  std::vector<Type> constants = {};

  // /* verify constants vector is filled properly */
  // std::cout << "constants, in main: \n";
  // DisplayVector(constants);

  /* number of MC runs to use in Sobol indices approximation */
   int N_MC = 10000;
  //  int N_MC = atoi(argv[1]);

  /* index set to compute sensitivity indices for */
  std::set<int> indices = {1};
  // std::cout << "indices: \n";
  // DisplaySet(indices);

  // /* coefficient of variation = std/mean for current parameter group */
  // Type CoV = 0.2;
  // std::vector<Type> CoV_Vector 
  //   = {0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35};

  /* file name to plot indices to when using CoV routines */
  std::string filename = "IndicesData.txt";

  /* specify distribution parameters of model parameters */
  std::vector<std::vector<Type> > 
    distroParams(dim, std::vector<Type>(2));

  distroParams[0][0] = 0;
  distroParams[0][1] = 1;
  distroParams[1][0] = 0;
  distroParams[1][1] = 4;
  distroParams[2][0] = 0;
  distroParams[2][1] = 9;
  distroParams[3][0] = 0;
  distroParams[3][1] = 16;

  /* construct a SobolIndices object */
  // SobolIndices sobol(Heston, constants, indices, distroParams,
  // 		     dim, N_MC, CoV);
  SobolIndices sobol(LinearModel, constants, indices, distroParams,
  			 dim,N_MC);

  // /* print member of SobolIndices object for verification */
  // sobol.DisplayMembers();

  clock_t tic = clock();

  /* compute sensitivity indices */
  std::cout << "computing sensitivity indices...\n\n";
   sobol.ComputeSensitivityIndices();
  // std::vector<std::vector<Type> > results 
  //   = sobol.PlotCoV(CoV_Vector, filename);

  std::cout << "...done.\n\n";

  Type toc = (Type)(clock() - tic) / CLOCKS_PER_SEC;
  std::cout << "total time: " << toc << "\n\n";

  /* display sensitivity indices */
  sobol.DisplayMembers();

  // /* write to file */
  // std::ofstream File("sigma.txt", std::ios::app);
  // File << N_MC << " " << actual_N_MC << " " << sobol.GetLowerIndex() << " " << sobol.GetTotalIndex() << "\n";
  // File.close();
}
