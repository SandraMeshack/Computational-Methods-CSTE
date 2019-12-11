
#pragma once // Include guard
#include "explicit_solver.h"

class RichardsonSolver : public ExplicitSolver {
public:
  RichardsonSolver();

protected:
  double nextStep(int spaceStep, int timeStep,
                  std::vector<std::vector<double> > *TSolution) const override;
};
