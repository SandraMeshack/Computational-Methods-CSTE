#pragma once // Include guard
#include "abstract_solver.h"
#include "laasonen_simple_implicit_solver.h"

class ExplicitSolver : public AbstractSolver {
  typedef void (ExplicitSolver::*ptrMethod)();

protected:
  /**
   * @brief Calculate the next point in space with the explicit scheme
   *
   * @param k position of the step in the array where the new value has to be
   * calculated
   */
  virtual double
  nextStep(int spaceStep, int timeStep,
           std::vector<std::vector<double> > *TSolution) const = 0;
  /**
   * @brief boolean to store if the explicit scheme is a three level scheme or
   * not
   *
   */
  bool threeLevelScheme = false;

  /**
   * @brief boolean to store if the first step should be computed with a
   * richardson extrapolation. Default solver and value are laasonen scheme
   * solver with richardson extrapolation
   *
   */
  bool withRichardsonsExtrapolation = true;

  /**
   * @brief Pointer to the solver that should be used for the first step
   *
   */
  AbstractSolver *firstStepSolver;

  /**
   * @brief Default first step solver for three-level scheme
   *
   */
  LaasonenSolver defaultFirstStepSolver;

  /**
   * @brief Compute the second level of data for three level scheme using
   * firstStepSolver and richardson's extrapolation according to the boolean
   * withRichardsonsExtrapolation
   *
   * @param TSolutionPtr
   */
  void computeFirstStep(std::vector<std::vector<double> > *TSolutionPtr);

public:
  /**
   * @brief Solve the issue if all the parameters have been defined
   *
   * Can throw exception if not all attributes have been properly
   * initialized
   *
   * @param deltaX : size if the space step
   * @param deltaT : size of the time step
   *
   * @return std::vector<std::vector<double>>
   */
  void
  solveRegularMeshes(double deltaX, double deltaT,
                     std::vector<std::vector<double> > *TSolutionPtr) override;

  /**
   * @brief Set the solver to compute the first step and if it should use a
   * richardson's extrapolation.
   *
   * @param solver : pointer to the solver that should be use to compute the
   * first step
   * @param useRichardsonsExtrapolation : whether or not use a richardson's
   * extrapolation
   */
  void setFirstStepSolver(AbstractSolver *solver,
                          bool useRichardsonsExtrapolation);

  /**
   * @brief Construct a new Explicit Solve object
   *
   */
  ExplicitSolver();

  /**
   * @brief Destroy the Explicit Solve object
   *
   */
  virtual ~ExplicitSolver();
};
