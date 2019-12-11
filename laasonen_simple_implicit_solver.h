#pragma once // Include guard
#include "implicit_solver.h"

class LaasonenSolver : public ImplicitSolver {
public:
  /**
   * @brief Construct a new Laasonen Solve object
   *
   */
  LaasonenSolver();

protected:
  /**
   * Compute the A matrix of the linear system given by Lassons scheme
   * and apply transformations required for solving the system with Thomas
   * Algorithm
   *
   */
  void initializeMatrixAForThomasAlgo() override;
  /**
   * Compute the B matrix of the linear system given by LAssonen scheme
   * and apply transformations required for solving the system with Thomas
   * Algorithm
   *
   */
  void initializeMatrixBForThomasAlgo(
      std::vector<double> *previousTimeStep) override;
};
