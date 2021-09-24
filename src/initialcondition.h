#ifndef INITIALCONDITION_H
#define INITIALCONDITION_H
#include "cantera/zerodim.h"
//#include "cantera/numerics.h"
#include "global_settings.h"
//#include "initialconditionfix.h"

#include <time.h>
#include "example_utils.h"
#include "struct.h"
#include "dataread.h"
#include "reaction.h"
using namespace Cantera;
void initial_to_chemical_part(cell *CS, int nc,shared_ptr<Cantera::Solution>* sols, shared_ptr<Cantera::ThermoPhase>* gass,Cantera::Reactor* reac);
void initial_to_znd(cell *CS, int nc,shared_ptr<Cantera::Solution>* sols, shared_ptr<Cantera::ThermoPhase>* gass,Cantera::Reactor** reac);
#endif
