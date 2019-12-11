#include "heat_diffusion_parameters.h"
#include <exception>

HeatDiffusionParameters::HeatDiffusionParameters() {
  timeLimitInitialized = false;
  internalTemperatureInitialized = false;
  surfaceTemperatureInitialized = false;
  diffusivityInitialized = false;
  widthInitialized = false;
};

void HeatDiffusionParameters::setTimeLimit(double newTimeStop) {
  if (0 < newTimeStop) {
    timeStop = newTimeStop;
    timeLimitInitialized = true;
  } else {
    throw(std::invalid_argument(
        "time limit is negative. As the resolution starts from t=0, time limit "
        "should be positive"));
  }
};

void HeatDiffusionParameters::setInternalTemperature(double Tin) {
  internalTemperature = Tin;
  internalTemperatureInitialized = true;
};

void HeatDiffusionParameters::setSurfaceTemperature(double Tsur) {
  surfaceTemperature = Tsur;
  surfaceTemperatureInitialized = true;
};

void HeatDiffusionParameters::setDiffusivity(double newDiffusivity) {
  if (newDiffusivity > 0) {
    diffusivity = newDiffusivity;
    diffusivityInitialized = true;
  } else {
    throw(std::invalid_argument("diffusivity should be positive"));
  }
};

void HeatDiffusionParameters::setWidth(double newWidth) {
  width = newWidth;
  widthInitialized = true;
}

bool HeatDiffusionParameters::checkInitialization() const {
  return (diffusivityInitialized && widthInitialized &&
          internalTemperatureInitialized && surfaceTemperatureInitialized &&
          timeLimitInitialized);
};

double HeatDiffusionParameters::getInternalTemperature() const {
  if (internalTemperatureInitialized) {
    return (internalTemperature);
  } else {
    throw(std::logic_error("internal Temperature not yet initialized"));
  }
};

double HeatDiffusionParameters::getSurfaceTemperature() const {
  if (surfaceTemperatureInitialized) {
    return (surfaceTemperature);
  } else {
    throw(std::logic_error("surface Temperature not yet initialized"));
  }
};

double HeatDiffusionParameters::getDiffusivity() const {
  if (diffusivityInitialized) {
    return (diffusivity);
  } else {
    throw(std::logic_error("difusivity not yet initialized"));
  }
};

double HeatDiffusionParameters::getTimeStop() const {
  if (timeLimitInitialized) {
    return (timeStop);
  } else {
    throw(std::logic_error("time limit not yet initialized"));
  }
};

double HeatDiffusionParameters::getWidth() const {
  if (widthInitialized) {
    return (width);
  } else {
    throw(std::logic_error("width not yet initialized"));
  }
};
