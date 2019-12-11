#include "implicit_solver.h"
#include <iostream>
#include <vector>

void ImplicitSolver::solveRegularMeshes(
    double pdeltaX, double pdeltaT,
    std::vector<std::vector<double> > *TSolutionPtr) {

  if (!parameters.checkInitialization()) {
    throw(std::invalid_argument(
        "Parameters have not yet been properly initialized"));
  }

  // variable used to store the solution at a step in time while we are still
  // computing it
  deltaX = pdeltaX;
  deltaT = pdeltaT;
  int numberOfSpacePoints = (int)(parameters.getWidth() / deltaX) + 1;
  std::vector<double> TSolutionAtOneTime(numberOfSpacePoints);

  // For t=0, we calculate the initial state of each point and store it
  std::fill(TSolutionAtOneTime.begin(), TSolutionAtOneTime.end(),
            parameters.getInternalTemperature());
  TSolutionAtOneTime[0] = parameters.getSurfaceTemperature();
  TSolutionAtOneTime[numberOfSpacePoints - 1] =
      parameters.getSurfaceTemperature();
  (*TSolutionPtr).push_back(TSolutionAtOneTime);
  TSolutionAtOneTime.clear();

  initializeMatrixAForThomasAlgo();

  // For the other time steps :
  for (int timeIndex = 1; (timeIndex * deltaT) <= parameters.getTimeStop();
       timeIndex++) {

    initializeMatrixBForThomasAlgo(&((*TSolutionPtr)[timeIndex - 1]));

    thomasAlgoSolve(&TSolutionAtOneTime);
    (*TSolutionPtr).push_back(TSolutionAtOneTime);
    TSolutionAtOneTime.clear();
  }
};

void ImplicitSolver::thomasAlgoSolve(std::vector<double> *TSolutionAtOneTime) {
  int sizeOfB = matrixB.size();
  (*TSolutionAtOneTime) = std::vector<double>(sizeOfB, -1);
  (*TSolutionAtOneTime)[sizeOfB - 1] = matrixB[sizeOfB - 1];
  for (int i = sizeOfB - 2; i >= 0; i--) {
    (*TSolutionAtOneTime)[i] =
        matrixB[i] - (*TSolutionAtOneTime)[i + 1] * (matrixA[i][i + 1]);
  }
}

ImplicitSolver::~ImplicitSolver(){};
