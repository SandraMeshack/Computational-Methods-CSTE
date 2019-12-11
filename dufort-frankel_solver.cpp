
#include "dufort-frankel_solver.h"
#include "laasonen_simple_implicit_solver.h"

DufortFrankelSolver::DufortFrankelSolver() {
  schemeName = "Dufort-Frankel";
  threeLevelScheme = true;
};

double DufortFrankelSolver::nextStep(
    int spaceStep, int timeStep,
    std::vector<std::vector<double> > *TSolution) const {

  // According to richardson scheme :
  double nextStep =
      ((*TSolution)[timeStep - 2][spaceStep] +
       2 * parameters.getDiffusivity() * (deltaT / (deltaX * deltaX)) *
           ((*TSolution)[timeStep - 1][spaceStep + 1] -
            (*TSolution)[timeStep - 2][spaceStep] +
            (*TSolution)[timeStep - 1][spaceStep - 1])) /
      (1 + 2 * parameters.getDiffusivity() * (deltaT / (deltaX * deltaX)));

  return (nextStep);
}
