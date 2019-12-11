/*! \file */

#include "abstract_solver.h"
#include "crank-nicholson_solver.h"
#include "dufort-frankel_solver.h"
#include "exact_solver.h"
#include "heat_diffusion_parameters.h"
#include "laasonen_simple_implicit_solver.h"
#include "richardson_solver.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

void resultToFile(std::string filename,
                  std::vector<std::vector<double> > *numericalSolution,
                  std::vector<std::vector<double> > *analyticalSolution,
                  double deltaX, double deltaT);
void addFirstStepRelsultToFile(
    std::fstream *outputFileStream,
    std::vector<std::vector<double> > *numericalSolution,
    std::vector<std::vector<double> > *analyticalSolution, double deltaX,
    double deltaT, std::string firstStepSolverName);
double uniform_norm(std::vector<std::vector<double> > *matrix);
double two_norm(std::vector<std::vector<double> > *matrix);

/**
 * @brief This programm aim to solve numericaly a one-space dimensional heat
 * conduction equation with several numerical schemes
 * The results should be in the ./Results folder
 *
 * First it will solve the problem for a given deltaX and deltaT grid with
 * different numerical schemes
 * -> Results/fullSolutionForSeveralSolvers for the results
 *
 * Then it will solve the same problem with the Laasonen implicit scheme for
 * several deltaT
 * -> Results/LaasonnenSeveralDeltat for the results
 *
 * Finally it will compute the first step for Richardson and Dufort-frankel
 * using different schemes
 * -> Results/FirstStepSolvers for the results
 *
 */
int main(int argc, const char **argv) {
  std::vector<std::vector<double> > numericalsolution, analyticalSolution;

  //==== Problem data =====
  double diffusivity = 93;         // cm²/hr
  double surfaceTemperature = 149; //°C
  double internalTemperature = 38; //°C
  double timeLimit = 0.5;          // hours
  double L = 31;                   // cm

  // ===== Grid size ======
  double deltaX = 0.05; // cm
  double deltaT = 0.01; // hours
  //====== Set up parameters =========

  HeatDiffusionParameters parameters = HeatDiffusionParameters();
  parameters.setDiffusivity(diffusivity);
  parameters.setSurfaceTemperature(surfaceTemperature);
  parameters.setInternalTemperature(internalTemperature);
  parameters.setWidth(L);
  parameters.setTimeLimit(timeLimit);

  // ====== Solving ========

  ExactSolver analyticalSolver = ExactSolver();
  LaasonenSolver laasonenSolver = LaasonenSolver();
  RichardsonSolver richardsonSolver = RichardsonSolver();
  CrankNicholsonSolver crankNicholsonSolver = CrankNicholsonSolver();
  DufortFrankelSolver dufortFrankelSolver = DufortFrankelSolver();

  /* SOLVE ENTIRELY WITH DIFFERENT SCHEMES */
  /*We use pointers to Abstract Solve because we cannot use a vector of
   * AbstractSolver as it is an abstract class*/
  std::vector<AbstractSolver *> solvers = {&laasonenSolver, &richardsonSolver,
                                           &crankNicholsonSolver,
                                           &dufortFrankelSolver};

  analyticalSolver.setParameters(parameters);
  analyticalSolver.solveRegularMeshes(deltaX, deltaT, &analyticalSolution);
  resultToFile("Results/fullSolutionForSeveralSolvers/Analytical.csv",
               &analyticalSolution, &analyticalSolution, deltaX, deltaT);

  for (auto &solverPtr : solvers) {
    (*solverPtr).setParameters(parameters);
    (*solverPtr).solveRegularMeshes(deltaX, deltaT, &numericalsolution);
    std::string filename = "Results/fullSolutionForSeveralSolvers/" +
                           (*solverPtr).getSchemeName() + ".csv";
    resultToFile(filename, &numericalsolution, &analyticalSolution, deltaX,
                 deltaT);
    numericalsolution.clear();
  }

  /* SOLVE WITH LAASONEN IMPLICIT SCHEME FOR DIFFERENT DELTA T */
  double laasonenDeltaTs[3] = {0.025, 0.05, 0.1};
  for (double laasonenDeltaT : laasonenDeltaTs) {
    std::string filename = "Results/LaasonnenSeveralDeltat/Laasonen Deltat = ";
    filename.append(std::to_string(laasonenDeltaT));
    filename.append(".csv");
    std::vector<std::vector<double> > laasonenDeltaTSolution,
        analyticalDeltaTSolution;
    laasonenSolver.solveRegularMeshes(deltaX, laasonenDeltaT,
                                      &laasonenDeltaTSolution);

    analyticalSolver.solveRegularMeshes(deltaX, laasonenDeltaT,
                                        &analyticalDeltaTSolution);
    resultToFile(filename, &laasonenDeltaTSolution, &analyticalDeltaTSolution,
                 deltaX, laasonenDeltaT);
  }

  /* COMPARE DIFFERENT SOLUTION FOR COMPUTING FIRST STEP FOR RICHARDSON AND
   * DUFORT FRANKEL THREE-LEVEL SCHEMES*/
  std::vector<AbstractSolver *> firstStepSolvers = {
      &analyticalSolver, &crankNicholsonSolver, &laasonenSolver};
  std::vector<AbstractSolver *> firstStepSolversWithRichardsonExtrapolation = {
      &crankNicholsonSolver, &laasonenSolver};
  std::vector<ExplicitSolver *> threeLevelsSolvers = {&dufortFrankelSolver,
                                                      &richardsonSolver};
  for (auto &threeLevelsSolver : threeLevelsSolvers) {
    parameters.setTimeLimit(5 * deltaT);
    threeLevelsSolver->setParameters(parameters);
    std::string filename = "Results/FirstStepSolvers/" +
                           threeLevelsSolver->getSchemeName() +
                           " with serveral first steps" + ".csv";
    std::fstream outputFileStream;
    outputFileStream.open(filename, std::fstream::out | std::fstream::trunc);
    outputFileStream << "deltaX:," << deltaX << ",deltaT:," << deltaT
                     << std::endl
                     << std::endl;

    for (auto &firstStepSolverPtr : firstStepSolvers) {
      threeLevelsSolver->setFirstStepSolver(firstStepSolverPtr, false);
      threeLevelsSolver->solveRegularMeshes(deltaX, deltaT, &numericalsolution);
      addFirstStepRelsultToFile(&outputFileStream, &numericalsolution,
                                &analyticalSolution, deltaX, deltaT,
                                firstStepSolverPtr->getSchemeName());
      numericalsolution.clear();
    }
    for (auto &firstStepSolverPtr :
         firstStepSolversWithRichardsonExtrapolation) {
      threeLevelsSolver->setFirstStepSolver(firstStepSolverPtr, true);
      threeLevelsSolver->solveRegularMeshes(deltaX, deltaT, &numericalsolution);
      addFirstStepRelsultToFile(
          &outputFileStream, &numericalsolution, &analyticalSolution, deltaX,
          deltaT, (firstStepSolverPtr->getSchemeName() + "with RE"));
      numericalsolution.clear();
    }
    outputFileStream.close();
  }
} /* End main*/

/*=========== NORMS ===============*/

/**
 * @brief Compute the two-norm / Euclidian norm of a matrix
 *
 * @param matrix : the matrix of which you want the norm
 * @return double : the two-norm of the matrix
 */
double two_norm(std::vector<std::vector<double> > *matrix) {
  int numberOfRows = (*matrix).size();
  if (numberOfRows == 0) {
    return (0);
  }
  int numberOfColumns = (*matrix)[0].size();
  if (numberOfColumns == 0) {
    return (0);
  }
  double norm = 0;
  for (std::vector<double> row : (*matrix)) {
    for (double elt : row) {
      norm += elt * elt;
    }
  }
  norm = sqrtf(norm);
  return norm;
}

/**
 * @brief Compute the uniform (maximum / infinity) norm of a matrix
 *
 * @param matrix : the matrix of which you want the norm
 * @return double : the uniform norm of the matrix
 */
double uniform_norm(std::vector<std::vector<double> > *matrix) {
  int numberOfRows = (*matrix).size();
  if (numberOfRows == 0) {
    return (0);
  }
  int numberOfColumns = (*matrix)[0].size();
  if (numberOfColumns == 0) {
    return (0);
  }
  double norm = 0;
  for (std::vector<double> row : (*matrix)) {
    for (double elt : row) {
      if (norm <= fabs(elt)) {
        norm = fabs(elt);
      }
    }
  }
  return norm;
}

/*=========== PRINT RESULT TO FILE ===============*/
/**
 * @brief Printing the numerical solution and its errors compared to the
 * analytical solution
 * - The file will start with the following information :
 * ```
 * deltaX:,<value of deltaX>,deltaT:,<value of deltaT>
 * ```
 *
 * - Then we get the uniform norm and the two norm for all the points computed :
 * ```
 * Errors measurements :
 * Uniform norm :,<value of the uniform norm>
 * Two norm :,<value of the two norm>
 * ```
 *
 * - Then all the points of the numerical solution will be printed :
 * ```
 * Numerical solution :
 * t\x, ........... <values of x_j> ...................
 * ...
 * ...
 * ...
 * <values of t^n>      T^n_j
 * ...
 * ...
 * ...
 * ```

 * - Finally the error for each points of the numerical solution and the two
 norm and uniform norm for each time step will be printed :
 * ```
 * Errors (analytical - numerical ):
 * t\x, ..... <values of x_j> ............  Uniform norm(t)  Two norm(t)
 * ...                                <uniform norm T(t0)> <two norm T(t0)>
 * ...                                              .             .
 * ...                                              .             .
 * ...                                              .             .
 * <values of t^n> Analytical(x_j,t^n) - T^n_j      .             .
 * ...                                              .             .
 * ...                                <uniform norm T(tn)> <two norm T(tn)>
 * ...                                              .             .
 * ```

 * @param filename : name of the output file, should end with .csv
 * @param numericalSolution : pointer to the numerical solution
 * @param analyticalSolution : pointer to the analytical solution. The
 analytical solution should be solved for the same grid as the numerical one
 * @param deltaX : space step size
 * @param deltaT : time step size
 */
void resultToFile(std::string filename,
                  std::vector<std::vector<double> > *numericalSolution,
                  std::vector<std::vector<double> > *analyticalSolution,
                  double deltaX, double deltaT) {

  std::fstream outputFileStream;
  outputFileStream.open(filename, std::fstream::out | std::fstream::trunc);
  outputFileStream << std::setprecision(16);
  int numberOfRows = (*numericalSolution).size();
  int numberOfColumns = (*numericalSolution)[0].size();
  std::vector<std::vector<double> > errorMatrix;
  std::vector<double> errorLine;
  for (int i = 0; i < numberOfRows; i++) {
    for (int j = 0; j < numberOfColumns; j++) {
      errorLine.push_back((*analyticalSolution)[i][j] -
                          (*numericalSolution)[i][j]);
    }
    errorMatrix.push_back(errorLine);
    errorLine.clear();
  }

  double uniformNorm = uniform_norm(&errorMatrix);
  double twoNorm = two_norm(&errorMatrix);

  outputFileStream << "deltaX:," << deltaX << ",deltaT:," << deltaT << std::endl
                   << std::endl;

  outputFileStream << "Errors measurements : " << std::endl;
  outputFileStream << "Uniform norm :," << uniformNorm << std::endl;
  outputFileStream << "Two norm :," << twoNorm << std::endl << std::endl;

  outputFileStream << "Numerical solution : " << std::endl;
  outputFileStream << "t\\x";
  for (int i = 0; i < numberOfColumns; i++) {
    outputFileStream << "," << i * deltaX;
  }
  outputFileStream << "\n";
  for (int i = 0; i < numberOfRows; i++) {
    outputFileStream << i * deltaT;
    for (int j = 0; j < numberOfColumns; j++) {
      outputFileStream << "," << (*numericalSolution)[i][j];
    }
    outputFileStream << std::endl;
  }

  outputFileStream << std::endl
                   << "Errors (analytical - numerical ): " << std::endl;
  outputFileStream << "t\\x";
  for (int i = 0; i < numberOfColumns; i++) {
    outputFileStream << "," << i * deltaX;
  }
  outputFileStream << ",,uniform_norm(t),two_norm(t)\n";
  for (int i = 0; i < numberOfRows; i++) {
    outputFileStream << i * deltaT;
    for (int j = 0; j < numberOfColumns; j++) {
      outputFileStream << "," << errorMatrix[i][j];
    }
    std::vector<std::vector<double> > line;
    line.push_back(errorMatrix[i]);
    outputFileStream << ",," << uniform_norm(&line) << "," << two_norm(&line)
                     << "\n";
  }
  outputFileStream.close();
}

/**
 * @brief Add the numerical solution and its errors compared to the
 * analytical solution in an opened file stream
 *
 * This function was meant to easily print in a file different numerical
 solution find with different first step solvers. It will add a section to an
 open file following this structure :
 *
 * - The new section will start with the following informations :
 * ```
 *  First Step Solver :, <name of the first step solver>
 *  Errors measurements :
 *  Uniform norm :,<value of the uniform norm>
 *  Two norm :,<value of the two norm>
 * ```
 * - Then all the points of the numerical solution will be printed :
 * ```
 * Numerical solution :
 * t\x, ........... <values of x_j> ...................
 * ...
 * ...
 * ...
 * <values of t^n>      T^n_j
 * ...
 * ...
 * ...
 * ```
 *
 * - Finally the error for each points of the numerical solution and the two
 norm and uniform norm for each time step will be printed :
 * ```
 * Errors (analytical - numerical ):
 * t\x, ..... <values of x_j> ............  Uniform norm(t)  Two norm(t)
 * ...                                <uniform norm T(t0)> <two norm T(t0)>
 * ...                                              .             .
 * ...                                              .             .
 * ...                                              .             .
 * <values of t^n> Analytical(x_j,t^n) - T^n_j      .             .
 * ...                                              .             .
 * ...                                <uniform norm T(tn)> <two norm T(tn)>
 * ...                                              .             .
 * ```
 *
 * @param outputFileStream : where the output data are sent
 * @param numericalSolution : pointer to the numerical solution
 * @param analyticalSolution : pointer to the analytical solution. The
 analytical solution should be solved for the same grid as the numerical one
 * @param deltaX : space step size
 * @param deltaT : time step size
 * @param firstStepSolverName : name of the solver used to compute the first
 step
 */

void addFirstStepRelsultToFile(
    std::fstream *outputFileStream,
    std::vector<std::vector<double> > *numericalSolution,
    std::vector<std::vector<double> > *analyticalSolution, double deltaX,
    double deltaT, std::string firstStepSolverName) {
  (*outputFileStream) << std::setprecision(16);
  int numberOfRows = (*numericalSolution).size();
  int numberOfColumns = (*numericalSolution)[0].size();
  std::vector<std::vector<double> > errorMatrix;
  std::vector<double> errorLine;
  for (int i = 0; i < numberOfRows; i++) {
    for (int j = 0; j < numberOfColumns; j++) {
      errorLine.push_back((*analyticalSolution)[i][j] -
                          (*numericalSolution)[i][j]);
    }
    errorMatrix.push_back(errorLine);
    errorLine.clear();
  }

  double uniformNorm = uniform_norm(&errorMatrix);
  double twoNorm = two_norm(&errorMatrix);
  (*outputFileStream) << std::endl
                      << " First Step Solver :," << firstStepSolverName
                      << std::endl;

  (*outputFileStream) << "Errors measurements : " << std::endl;
  (*outputFileStream) << "Uniform norm :," << uniformNorm << std::endl;
  (*outputFileStream) << "Two norm :," << twoNorm << std::endl << std::endl;

  (*outputFileStream) << "Numerical solution : " << std::endl;
  (*outputFileStream) << "t\\x";
  for (int i = 0; i < numberOfColumns; i++) {
    (*outputFileStream) << "," << i * deltaX;
  }
  (*outputFileStream) << "\n";
  for (int i = 0; i < numberOfRows; i++) {
    (*outputFileStream) << i * deltaT;
    for (int j = 0; j < numberOfColumns; j++) {
      (*outputFileStream) << "," << (*numericalSolution)[i][j];
    }
    (*outputFileStream) << std::endl;
  }

  (*outputFileStream) << std::endl
                      << "Errors (analytical - numerical ): " << std::endl;
  (*outputFileStream) << "t\\x";
  for (int i = 0; i < numberOfColumns; i++) {
    (*outputFileStream) << "," << i * deltaX;
  }
  (*outputFileStream) << ",,uniform_norm(t),two_norm(t)\n";
  for (int i = 0; i < numberOfRows; i++) {
    (*outputFileStream) << i * deltaT;
    for (int j = 0; j < numberOfColumns; j++) {
      (*outputFileStream) << "," << errorMatrix[i][j];
    }
    std::vector<std::vector<double> > line;
    line.push_back(errorMatrix[i]);
    (*outputFileStream) << ",," << uniform_norm(&line) << "," << two_norm(&line)
                        << "\n";
  }
}
