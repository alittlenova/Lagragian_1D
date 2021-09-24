#ifndef EXCHANGING_H
#define EXCHANGING_H
#include "cantera/zerodim.h"
#include "cantera/thermo/IdealGasPhase.h"
//#include "cantera/numerics.h"
#include "global_settings.h"
#include "initialconditionfix.h"

#include <time.h>
#include "example_utils.h"
#include "struct.h"
using namespace Cantera;
void exchange_to_hydro_part(cell *CS, int nc, shared_ptr<Cantera::ThermoPhase>* gass, Cantera::Reactor* reac);
void exchange_to_chemical_part(cell *CS, int nc, shared_ptr<Cantera::ThermoPhase>* gass, Cantera::Reactor* reac);
void Try_temperature( shared_ptr<Cantera::ThermoPhase>* gas,double intE);
#endif
