#include "abstract_solver.h"
#include <exception>

void AbstractSolver::setParameters(HeatDiffusionParameters problemParameters) {
  if (problemParameters.checkInitialization()) {
    parameters = problemParameters;
  } else {
    throw(std::invalid_argument(
        "Parameters have not yet been properly initialized"));
  }
};

std::string AbstractSolver::getSchemeName() const { return (schemeName); };

AbstractSolver::~AbstractSolver(){};
