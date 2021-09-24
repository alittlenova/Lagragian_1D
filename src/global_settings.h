#ifndef GLOBALSETTINGS_INITIAL
#define GLOBALSETTINGS_INITIAL
#include "struct.h"
#include "math.h"
#include "cantera/zerodim.h"
//#include "cantera/numerics.h"
void globalsettings_initial(void);
double funcp(double time);
extern int nctest;
extern double phil;
extern double phir;
extern double xl;
extern double gau;
extern int ncell;
extern double final_time;
extern int total_image;
#endif
