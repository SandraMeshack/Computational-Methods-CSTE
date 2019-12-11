#pragma once // Include guard
#include "heat_diffusion_parameters.h"
#include <functional>
#include <string>
#include <vector>
/**
 * @brief Abstract implementation of a solver
 * Aim to solve the unsteady one-space dimensional heat conduction equation in
 * Cartesian coordinates.
 *
 */

class AbstractSolver {
public:
  /**
   * @brief Set up parameters of the problem to solve
   *
   * @param problemParameters : parameters of the problem to solve
   */
  void setParameters(HeatDiffusionParameters problemParameters);

  /**
   * @brief Solve the issue if all the parameters have been defined
   *
   * Can throw exception if not all attributes have been properly initialized
   *
   * @param deltaX : size if the space step
   * @param deltaT : size of the time step
   * @param TSolutionPtr : pointer to a vector of vector of double to store the
   * solution
   *
   */
  virtual void
  solveRegularMeshes(double deltaX, double deltaT,
                     std::vector<std::vector<double> > *TSolutionPtr) = 0;

  /**
   * @brief Get the name of the scheme used by the solver
   *
   * @return std::string name of the scheme
   */
  std::string getSchemeName() const;

  /**
   * @brief Destroy the Abstract Solve object
   *
   */
  virtual ~AbstractSolver();

protected:
  /**
   * @brief parameters of the heat diffusion problem to solve
   *
   */
  HeatDiffusionParameters parameters;

  /**
   * @brief Space Step
   *
   */
  double deltaX;
  /**
   * @brief Time Step
   *
   */
  double deltaT;

  /**
   * @brief name of the scheme used by the solver
   *
   * We can't make a const out of it because we need to initialize it in the
   * derived classes
   */
  std::string schemeName;
};
