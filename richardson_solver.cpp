#include "richardson_solver.h"

RichardsonSolver::RichardsonSolver() {
  schemeName = "Richardson";
  threeLevelScheme = true;
};

double
RichardsonSolver::nextStep(int spaceStep, int timeStep,
                           std::vector<std::vector<double> > *TSolution) const {
  // According to richardson scheme :
  double nextStep = (*TSolution)[timeStep - 2][spaceStep] +
                    2 * parameters.getDiffusivity() *
                        (deltaT / (deltaX * deltaX)) *
                        ((*TSolution)[timeStep - 1][spaceStep + 1] -
                         2 * (*TSolution)[timeStep - 1][spaceStep] +
                         (*TSolution)[timeStep - 1][spaceStep - 1]);

  return (nextStep);
}
