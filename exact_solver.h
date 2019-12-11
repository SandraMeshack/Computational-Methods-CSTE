
#pragma once // Include guard
#include "explicit_solver.h"
#include <math.h>

class ExactSolver : public ExplicitSolver {
public:
  ExactSolver();

protected:
  double nextStep(int spaceStep, int timeStep,
                  std::vector<std::vector<double> > *TSolution) const override;
};
