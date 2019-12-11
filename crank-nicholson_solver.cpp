#include "crank-nicholson_solver.h"
#include <iostream>

CrankNicholsonSolver::CrankNicholsonSolver() {
  schemeName = "Crank-Nicholson";
};

void CrankNicholsonSolver::initializeMatrixAForThomasAlgo() {
    // c=s/2 as defined in the report
  double c = parameters.getDiffusivity() * deltaT / (2 * deltaX * deltaX);
  int numberOfSpacePoints = (int)((parameters.getWidth()) / deltaX) + 1;
  matrixA = std::vector<std::vector<double> >(
      numberOfSpacePoints, std::vector<double>(numberOfSpacePoints, 0.));
  matrixA[0][0] = 1;
  matrixA[numberOfSpacePoints - 1][numberOfSpacePoints - 1] = 1;

  for (int spaceIndex = 1; spaceIndex < numberOfSpacePoints - 1; spaceIndex++) {
    matrixA[spaceIndex][spaceIndex] = 1;
    matrixA[spaceIndex][spaceIndex + 1] =
        -c / ((1 + (2 * c)) + c * matrixA[spaceIndex - 1][spaceIndex]);
  }
}

void CrankNicholsonSolver::initializeMatrixBForThomasAlgo(
    std::vector<double> *TpreviousTimeStep) {
// c=s/2 as defined in the report
  double c = parameters.getDiffusivity() * deltaT / (2 * deltaX * deltaX);
  int numberOfSpacePoints = (int)((parameters.getWidth()) / deltaX) + 1;
  matrixB = std::vector<double>(numberOfSpacePoints, -1);
  matrixB[0] = ((*TpreviousTimeStep)[0] / 1);
  for (int spaceIndex = 1; spaceIndex < numberOfSpacePoints - 1; spaceIndex++) {
    double di = (1 - (2 * c)) * (*TpreviousTimeStep)[spaceIndex] +
                c * (*TpreviousTimeStep)[spaceIndex + 1] +
                c * (*TpreviousTimeStep)[spaceIndex - 1];
    matrixB[spaceIndex] =
        ((di - (-c * matrixB[spaceIndex - 1])) /
         (1 + (2 * c) + c * matrixA[spaceIndex - 1][spaceIndex]));
  }
  double di = (*TpreviousTimeStep)[numberOfSpacePoints - 1];
  matrixB[numberOfSpacePoints - 1] =
      di / matrixA[numberOfSpacePoints - 1][numberOfSpacePoints - 1];
}
