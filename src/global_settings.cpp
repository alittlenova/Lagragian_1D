#include "global_settings.h"
int		nctest;
double	phil;
double	phir;
double	xl;
double	gau;
int		ncell;
double	final_time;
int total_image;
//Validation Parameters
/*void globalsettings_initial(void)
{
	nctest=50;
	phil=0;
	phir=1;
	xl=0.0;
	gau=1.4;
	ncell=4000;
	endtime=0.0002;
}
double funcp(double time)
{
	double mean=1500;
	double frequency=0.1;
	double Amplitude=10;
	return mean+Amplitude*sin(frequency*time);
}
*/
/*void globalsettings_initial(void) //H2 up=1500m/s reaction double region
{
	nctest=50;
	phil=0;
	phir=1e-2;
	xl=0.0;
	gau=1.4;
	ncell=1600;
	endtime=5.0e-6;
}*/
/*void globalsettings_initial(void)	//Test 3500
{
	nctest=50;
	phil=0;
	phir=2.0e-3;
	xl=0.0;
	gau=1.4;
	ncell=200;
	final_time=5.0e-6;
	total_image=1000;
}*/

/*void globalsettings_initial(void)	//Strong Ignition 1100K H2O2
{
	nctest=50;
	phil=0;
	phir=4.0e-3;
	xl=0.0;
	gau=1.4;
	ncell=2000;
	final_time=8.0e-5;
	total_image=1000;
}*/
/*void globalsettings_initial(void)	//C2H4 No fluctuation
{
	nctest=50;
	phil=0;
	phir=3.5e-1;
	xl=0.0;
	gau=1.4;
	ncell=2000;
	final_time=2.0e-3;
	total_image=2000;
}*/

/*void globalsettings_initial(void)	//C2H4 No fluctuation Zoom in
{
	nctest=50;
	phil=0;
	phir=5.0e-3;
	xl=0.0;
	gau=1.4;
	ncell=200;
	final_time=3.2e-3;
	total_image=2000;
}*/
/*void globalsettings_initial(void)	//C2H4 cascade case
{
	nctest=50;
	phil=0;
	phir=1.5e-2;
	xl=0.0;
	gau=1.4;
	ncell=5000;
	final_time=1.0e-4;
	total_image=2000;
}*/
void globalsettings_initial(void)	//Test Case 1 
{
	nctest=50;
	phil=0;
	phir=2.0;
	xl=0.0;
	gau=1.4;
	ncell=2000;
	final_time=0.5;
	total_image=2000;
}
/*void globalsettings_initial(void)
{
	nctest=50;
	phil=0;
	phir=1e-2;
	xl=0.0;
	gau=1.4;
	ncell=1600;
	endtime=5.0e-6;
}*/
double funcp(double time)      //My previous system settings
{
//	double mean=1800;
//	double mean=1436.35;//950K H2
//	double mean=1626.345;//1100K H2
//	double mean=1259.0552;//C2H4 Mt=4.5 Saif
//	double mean=2582.9;
	double mean=0;
//	double frequency=453.5*1000; //45.35khz
	double frequency=200*1000;
	double Amplitude=0.20*mean;
//	double Amplitude=0;
	return mean+Amplitude*sin(2*M_PI*frequency*time);
//	return mean;
}
/*double funcp(double time)
{
	double mean=1258.71;
	double frequency=80000;
	double Amplitude=0.05*mean;
	*/
