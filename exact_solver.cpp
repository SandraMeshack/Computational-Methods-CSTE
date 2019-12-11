#include "exact_solver.h"
#include "stdlib.h"

ExactSolver::ExactSolver() { schemeName = "Analytical"; };

double
ExactSolver::nextStep(int spaceStep, int timeStep,
                      std::vector<std::vector<double> > *TSolution) const {

  double nextStep = 0;
  double L = parameters.getWidth();
  for (int m = 61; m > 0; m -= 2) {
    nextStep += exp(-parameters.getDiffusivity() * (m * M_PI / L) *
                    (m * M_PI / L) * timeStep * deltaT) *
                2 / (m * M_PI) * sin(m * M_PI * spaceStep * deltaX / L);
  }
  nextStep = nextStep * 2 *
                 (parameters.getInternalTemperature() -
                  parameters.getSurfaceTemperature()) +
             parameters.getSurfaceTemperature();
  return (nextStep);
}
