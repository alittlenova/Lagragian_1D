#ifndef INITIALCONDITIONFIX_H
#define INITIALCONDITIONFIX_H
#include "cantera/zerodim.h"
#include "cantera/thermo/IdealGasPhase.h"
//#include "cantera/numerics.h"
#include "global_settings.h"

#include <time.h>
#include "example_utils.h"
#include "struct.h"
#include "ReimannSolver.h"
#include "initialcondition.h"
#include "reaction.h"
using namespace Cantera;
void place_shock_condition(cell* CS, int place_lcell, cell* lcell,cell* rcell,int nc,shared_ptr<Cantera::Solution>* sol, shared_ptr<Cantera::ThermoPhase>* gas);
#endif
