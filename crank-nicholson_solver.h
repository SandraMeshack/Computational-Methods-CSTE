#pragma once // Include guard
#include "implicit_solver.h"
#include <string>

class CrankNicholsonSolver : public ImplicitSolver {
public:
  CrankNicholsonSolver();

protected:
  /**
   * Compute the A matrix of the linear system given by Crank-nicholson scheme
   * and apply transformations required for solving the system with Thomas
   * Algorithm
   *
   */
  void initializeMatrixAForThomasAlgo() override;
  /**
   * Compute the B matrix of the linear system given by Crank-nicholson scheme
   * and apply transformations required for solving the system with Thomas
   * Algorithm
   *
   */
  void initializeMatrixBForThomasAlgo(
      std::vector<double> *previousTimeStep) override;
};
