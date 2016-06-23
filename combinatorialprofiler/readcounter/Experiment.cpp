#include "Experiment.h"

Experiment::Experiment()
: ndsi(NDSIS::noNDSI)
{}

Experiment::Experiment(std::string n)
: name(std::move(n)), ndsi(NDSIS::noNDSI)
{}
