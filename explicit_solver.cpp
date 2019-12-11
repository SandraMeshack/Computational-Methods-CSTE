#include "explicit_solver.h"
#include "crank-nicholson_solver.h"
#include "exact_solver.h"
#include "laasonen_simple_implicit_solver.h"
#include <iostream>
#include <vector>

ExplicitSolver::ExplicitSolver() { firstStepSolver = &defaultFirstStepSolver; };

void ExplicitSolver::solveRegularMeshes(
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

  if (threeLevelScheme) {
    computeFirstStep(TSolutionPtr);
  } else {
    // For t=0, we calculate the initial state of each point and store it

    std::fill(TSolutionAtOneTime.begin(), TSolutionAtOneTime.end(),
              parameters.getInternalTemperature());
    TSolutionAtOneTime[0] = parameters.getSurfaceTemperature();
    TSolutionAtOneTime[numberOfSpacePoints - 1] =
        parameters.getSurfaceTemperature();
    (*TSolutionPtr).push_back(TSolutionAtOneTime);
  }
  // For the other time steps :
  for (int timeIndex = (threeLevelScheme ? 2 : 1);
       (timeIndex * deltaT) <= (parameters.getTimeStop()); timeIndex++) {
    TSolutionAtOneTime.clear();
    // using left boundary condition to get the first value
    TSolutionAtOneTime.push_back(parameters.getSurfaceTemperature());
    // using the explicit scheme implemented in nextStep to get the values
    // from the middle
    for (int spaceIndex = 1; spaceIndex < numberOfSpacePoints - 1;
         spaceIndex++) {
      TSolutionAtOneTime.push_back(
          nextStep(spaceIndex, timeIndex, TSolutionPtr));
    }
    // using right boundary condition to get the last value
    TSolutionAtOneTime.push_back(parameters.getSurfaceTemperature());

    (*TSolutionPtr).push_back(TSolutionAtOneTime);
  };

  // Add a two-step boolean member to the explicit class.
};

void ExplicitSolver::computeFirstStep(
    std::vector<std::vector<double> > *TSolutionPtr) {

  HeatDiffusionParameters firstStepParameters =
      HeatDiffusionParameters(parameters);
  firstStepParameters.setTimeLimit(1 * deltaT);
  (*firstStepSolver).setParameters(firstStepParameters);

  if (withRichardsonsExtrapolation) {
    // Computing with a Richardson Extrapolation Tc = 4/3 Tb - Ta/3. Working for
    // FTCS explicit, Laasonen, Cranck-Nicholson
    std::vector<std::vector<double> > TSolutionGridA, TSolutionGridB;
    std::vector<double> TSolutionOneStep;
    (*firstStepSolver).solveRegularMeshes(deltaX, deltaT, &TSolutionGridA);
    (*firstStepSolver)
        .solveRegularMeshes(deltaX / 2, deltaT / 4, &TSolutionGridB);

    // Getting initial state for t=0 from grid A:
    (*TSolutionPtr).push_back(TSolutionGridA[0]);

    // Getting richardson extrapolation for t=DeltaT;
    int numberOfSpacePoints = (int)((parameters).getWidth() / deltaX) + 1;
    for (int j = 0; j < numberOfSpacePoints; j++) {
      TSolutionOneStep.push_back(
          (4 * TSolutionGridB.back()[2 * j] - TSolutionGridA.back()[j]) / 3);
    }
    (*TSolutionPtr).push_back(TSolutionOneStep);
  }
  // Computing without Richarson Extrapolation
  else {
    (*firstStepSolver).solveRegularMeshes(deltaX, deltaT, TSolutionPtr);
  }
}

void ExplicitSolver::setFirstStepSolver(AbstractSolver *solver,
                                        bool useRichardsonsExtrapolation) {
  firstStepSolver = solver;
  withRichardsonsExtrapolation = useRichardsonsExtrapolation;
};

ExplicitSolver::~ExplicitSolver(){};
