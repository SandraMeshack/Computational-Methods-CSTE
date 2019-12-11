#pragma once // Include guard
#include <functional>
#include <vector>

/**
 * @brief Parameters of 1D heat diffusion problem
 * Used to defined the problem that should be solved
 */

class HeatDiffusionParameters {

public:
  /**
   * @brief Construct a new Heat Diffusion Parameters object and initialized
   * internal boolean to check initialization to false
   *
   */
  HeatDiffusionParameters();

  /**
   * @brief Set up time limitation
   *
   * @param timeStop : time limit for stoping the simulation
   */
  void setTimeLimit(double timeStop);

  /**
   * @brief Set up size of the wall
   *
   * @param width : size of the wall
   */
  void setWidth(double width);

  /**
   * @brief Set up the diffusivity
   *
   * @param diffusivity : diffusivity of the material
   */
  void setDiffusivity(double diffusivity);

  /**
   * @brief Set up internal Temperature of the wall at t=0
   *
   * @param internalTemp : internal temperature of the wall at t=0
   */
  void setInternalTemperature(double internalTemp);

  /**
   * @brief Set up the temperature of both surfaces of the wall for any t
   *
   * @param surfaceTemp : temperature on both surfaces of the wall for any t
   */
  void setSurfaceTemperature(double surfaceTemp);

  /**
   * @brief Check if all parameters of the problem have been properly
   * initialized
   *
   */
  bool checkInitialization() const;

  /**
   * @brief Get the initial temperature inside the wall
   *
   */
  double getInternalTemperature() const;

  /**
   * @brief Get the temperature of both surfaces of the wall
   *
   * @return double
   */
  double getSurfaceTemperature() const;

  /**
   * @brief Get the time limit
   *
   * @return double
   */
  double getTimeStop() const;

  /**
   * @brief Get the width of the wall
   *
   * @return double
   */
  double getWidth() const;

  /**
   * @brief Get the Diffusivity object
   *
   * @return double
   */
  double getDiffusivity() const;

protected:
  /**
   * @brief Diffusivity of the material of the wall
   *
   */
  double diffusivity;
  /**
   * @brief width of the wall
   *
   */
  double width;

  /**
   * @brief time limit
   *
   */
  double timeStop;
  /**
   * @brief internal temperature of the wall at t=0
   *
   */
  double internalTemperature;
  /**
   * @brief temperature maintained on both surfaces of the wall
   *
   */
  double surfaceTemperature;

private:
  /**
   * @brief Check if time limit has been initialized
   *
   */
  bool timeLimitInitialized;

  /**
   * @brief Check if wall width has been initialized
   *
   */
  bool widthInitialized;

  /**
   * @brief Check if internal temperature has been initialized
   *
   */
  bool internalTemperatureInitialized;

  /**
   * @brief Check if surface temperature has been initialized
   *
   */
  bool surfaceTemperatureInitialized;

  /**
   * @brief Check if diffusivity has been initialized.
   *
   */
  bool diffusivityInitialized;
};
