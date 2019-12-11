#pragma once // Include guard
#include "abstract_solver.h"

class ImplicitSolver : public AbstractSolver {
protected:
  /**
   * @brief The A matrix from AX=B equation of the Thomas Algorithm
   * A is a matrix with 1 it's main diagonal and other values in the diagonal
   * above it.
   * AX=B coresponding to the linear system that we are trying to solve
   * for each time step once the transformations required by the thomas
   * algorithm are applied to A and B
   */
  std::vector<std::vector<double> > matrixA;

  /**
   * @brief The B vector from AX=B equation of the Thomas Algorithm
   * AX=B coresponding to the linear system that we are trying to solve
   * for each time step once the transformations required by the thomas
   * algorithm are applied to A and B.
   */
  std::vector<double> matrixB;

  /**
   * @brief Initialized A matrix with transformations required by Thomas
   * Algorithm This matrix will be reuse for all time steps Specific to each
   * scheme
   */
  virtual void initializeMatrixAForThomasAlgo() = 0;

  /**
   * @brief Initialized B vector with transformations required by Thomas
   * Algorithm This vector needs to be recalculated for each new time step.
   * Specific to each scheme
   *
   * @param previousTimeStep  : Values of the temperature at the previous time
   * step
   */
  virtual void
  initializeMatrixBForThomasAlgo(std::vector<double> *previousTimeStep) = 0;

  /**
   * @brief Solve the linear system AX=B
   * with A = matrixA and B=matrixB following the Thomas Algorithm.
   * transformations to A and B matrix are applied before calling this
   * functions. They should be applied by initializeMatrix<X>ForThomasAlgo()
   * functions.
   *
   * @param TSolutionAtOneTime : where to store the solution given by the
   * algorithm
   */
  void thomasAlgoSolve(std::vector<double> *TSolutionAtOneTime);

public:
  /**
   * @brief Solve the issue if all the parameters have been defined
   *
   * Can throw exception if not all attributes have been properly initialized
   *
   * @param deltaX : size if the space step
   * @param deltaT : size of the time step
   *
   * @return std::vector<std::vector<double>>
   */
  void
  solveRegularMeshes(double deltaX, double deltaT,
                     std::vector<std::vector<double> > *TSolutionPtr) override;

  virtual ~ImplicitSolver();
};
