#ifndef REIMANNSOLVER_H
#define REIMANNSOLVER_H
#include "cantera/zerodim.h"
//#include "cantera/numerics.h"

#include "global_settings.h"
#include <time.h>
#include "example_utils.h"
#include "exchanging.h"
#include "struct.h"
#include "math.h"
#include "ApproximateSolver_lagrangian.h"
using namespace Cantera;
double initialiter(double cl,double cr,double al,double ar,double pl,double pr,double ul,double ur,double gammal, double gammar);
double iteration(double cl,double cr,double al,double ar,double pl,double pr,double ul,double ur,double gammal, double gammar);
void fluidcons(double *rl,double *rr,double *ul,double *ur,double *pl,double *pr,double *f0,double *f1,double *f2,double gammal,double gammar);
double maxdt(cell *CS,int number);
double march(cell *CS,join *JS,cell *GS,cell *CS2,join *JS2,cell *GS2,int nc, double time,double endtime,shared_ptr<Cantera::ThermoPhase>* gass,Cantera::Reactor* reac,shared_ptr<Cantera::ThermoPhase>* gass2,Cantera::Reactor* reac2,Cantera::ReactorNet* sim,Cantera::ReactorNet* sim2);
double marchtest(cell *CS,join *JS,cell *GS,cell *CS2,join *JS2,cell *GS2,cell *left,int nc);
cell* initiate(int nc, double *rini, double *uini, double *pini, double *gaini);
double shockjumpp(double up,double ur,double pr,double rr,double gammar);
double expansionjumpp(double up,double ur,double pr,double rr,double gammar);
double shockjumpr(double up,double ur,double pr,double rr,double gammar);
void plot(cell *CS,join *JS,int nc,shared_ptr<Cantera::ThermoPhase>* GS,int plotnum,double tm);
void copyreactor(shared_ptr<Cantera::ThermoPhase>* gass,Cantera::Reactor* reac,shared_ptr<Cantera::ThermoPhase>* gass2,Cantera::Reactor* reac2,double tm,int nc);
void copycell(cell *CS1,cell *CS2);
void copyjoin(join *JS1,join *JS2);
#endif
