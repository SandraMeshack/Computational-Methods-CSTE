#pragma once // Include guard
#include "explicit_solver.h"
class DufortFrankelSolver : public ExplicitSolver {
public:
  DufortFrankelSolver();

protected:
  double nextStep(int spaceStep, int timeStep,
                  std::vector<std::vector<double> > *TSolution) const override;
};
