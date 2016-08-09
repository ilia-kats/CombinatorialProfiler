#include "Experiment.h"

Experiment::Experiment()
: dsi(DSIS::noDSI)
{}

Experiment::Experiment(std::string n)
: name(std::move(n)), dsi(DSIS::noDSI)
{}
