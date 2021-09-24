#ifndef REACTION_H
#define REACTION_H
#include "cantera/zerodim.h"
#include "cantera/thermo/IdealGasPhase.h"
//#include "cantera/numerics.h"

#include "global_settings.h"
#include <time.h>
#include "example_utils.h"
using namespace Cantera;

void chemicalinitiation(shared_ptr<Cantera::Solution>* sol, shared_ptr<Cantera::ThermoPhase>* gas,Reactor*	 r,ReactorNet* sim, int cell);
#endif
