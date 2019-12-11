#include "laasonen_simple_implicit_solver.h"
#include <iostream>

LaasonenSolver::LaasonenSolver() { schemeName = "Laasonen"; };

void LaasonenSolver::initializeMatrixAForThomasAlgo() {
  double s = parameters.getDiffusivity() * deltaT / (deltaX * deltaX);
  int numberOfSpacePoints = (int)(parameters.getWidth() / deltaX) + 1;
  matrixA = std::vector<std::vector<double> >(
      numberOfSpacePoints, std::vector<double>(numberOfSpacePoints, 0.));
  matrixA[0][0] = 1;
  matrixA[numberOfSpacePoints - 1][numberOfSpacePoints - 1] = 1;
  for (int spaceIndex = 1; spaceIndex < numberOfSpacePoints - 1; spaceIndex++) {
    matrixA[spaceIndex][spaceIndex] = 1;
    matrixA[spaceIndex][spaceIndex + 1] =
        -s / (1 + 2 * s + s * matrixA[spaceIndex - 1][spaceIndex]);
  }
}
void LaasonenSolver::initializeMatrixBForThomasAlgo(
    std::vector<double> *TpreviousTimeStep) {
  double s = parameters.getDiffusivity() * deltaT / (deltaX * deltaX);
  int numberOfSpacePoints = (int)((parameters.getWidth()) / deltaX) + 1;
  matrixB = std::vector<double>(numberOfSpacePoints, -1);
  matrixB[0] = ((*TpreviousTimeStep)[0] / 1);
  for (int spaceIndex = 1; spaceIndex < numberOfSpacePoints - 1; spaceIndex++) {
    double di = (*TpreviousTimeStep)[spaceIndex];
    matrixB[spaceIndex] =
        ((di - (-s * matrixB[spaceIndex - 1])) /
         (1 + (2 * s) + s * matrixA[spaceIndex - 1][spaceIndex]));
  }
  double di = (*TpreviousTimeStep)[numberOfSpacePoints - 1];
  matrixB[numberOfSpacePoints - 1] =
      di / matrixA[numberOfSpacePoints - 1][numberOfSpacePoints - 1];
}
