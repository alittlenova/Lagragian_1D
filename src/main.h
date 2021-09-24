#ifndef MAIN_H
#define MAIN_H
#include "cantera/zerodim.h"
#include "cantera/thermo/IdealGasPhase.h"
//#include "cantera/numerics.h"
#include "global_settings.h"
#include "initialcondition.h"
#include <time.h>
#include <sys/time.h>
#include "example_utils.h"
#include "ReimannSolver.h"
#include "initialconditionfix.h"
#include "exchanging.h"
#include "reaction.h"
#include "struct.h"

#include "dataread.h"
int kinetics2(int np, void* p);
#endif
