#include "main.h"
using namespace Cantera;
using std::cout;
using std::endl;
int kinetics()
{
    cout << "Constant-pressure ignition of a "
         << "hydrogen/oxygen/nitrogen"
         " mixture \nbeginning at T = 300 K and P = 1 atm." << endl;


    // create an ideal gas mixture that corresponds to GRI-Mech
    // 3.0
	// set up gas models
	shared_ptr<Cantera::Solution> sol[ncell];
	shared_ptr<Cantera::ThermoPhase> gas[ncell];
	
	// set up image gas models
	shared_ptr<Cantera::Solution> sol2[ncell];
	shared_ptr<Cantera::ThermoPhase> gas2[ncell];
	// create reactors
    Reactor r[ncell];
	// create image reactors
    Reactor r2[ncell];
    
    // create a container object to run the simulation
    // and add the reactor to it
    ReactorNet *sim=new ReactorNet[ncell];
    // and add the image reactor to it
    ReactorNet *sim2=new ReactorNet[ncell];
    
	chemicalinitiation(sol,gas,r,sim,ncell);
	chemicalinitiation(sol2,gas2,r2,sim2,ncell);
	//Initial Hydro part
	cell *CS,*CS2;
	cell *GC,*GC2;
	join *JS,*JS2;
	
	cell *TCS;

	double rini[ncell];
	double uini[ncell];
	double pini[ncell];
	double gaini[ncell];
    double tm,tc,tp;
    double dt;
    
	for(int i=0;i<ncell;i++)
	{
		if(i<1000||i>1020)
		{
			//rini[i]=shockjumpr(funcp(0),0,gas[i]->pressure(),gas[i]->density(),gas[i]->cp_mole()/gas[i]->cv_mole());
			//uini[i]=funcp(0);
			//pini[i]=shockjumpp(funcp(0),0,gas[i]->pressure(),gas[i]->density(),gas[i]->cp_mole()/gas[i]->cv_mole());
			//rini[i]=2*gas[i]->density();
			//uini[i]=0;
			//pini[i]=3*gas[i]->pressure();
			//gaini[i]=gas[i]->cp_mole()/gas[i]->cv_mole();
			
			rini[i]=1;
			uini[i]=0;
			pini[i]=1;
			gaini[i]=1.4;
		}
		else
		{
			//rini[i]=gas[i]->density();
			//uini[i]=0;
			//pini[i]=gas[i]->pressure();
			//gaini[i]=gas[i]->cp_mole()/gas[i]->cv_mole();
			
			rini[i]=1;
			uini[i]=0;
			pini[i]=2;
			gaini[i]=1.4;
		}
	}
	
	CS=initiate(ncell,rini,uini,pini,gaini);
	GC=CS->lc;
	JS=CS->lj;
	
	CS2=initiate(ncell,rini,uini,pini,gaini);
	GC2=CS2->lc;
	JS2=CS2->lj;
	
    // main loop
    clock_t t0 = clock();        // save start time
    system("rm ./output/*.csv");
    int pic=0;	
    int picnumber=0;
    tm=0;
    
    
    GC->ga=CS->ga;
    GC->prim[0]=CS->prim[0];
    GC->prim[1]=funcp(0);
	GC->prim[2]=shockjumpp(GC->prim[1],CS->prim[1],CS->prim[2],1/(CS->prim[0]),CS->ga);
	
	GC2->ga=CS2->ga;
	GC2->prim[0]=CS2->prim[0];
	GC2->prim[1]=funcp(0);
	GC2->prim[2]=shockjumpp(GC2->prim[1],CS2->prim[1],CS2->prim[2],1/(CS2->prim[0]),CS2->ga);
	
	//initial_to_chemical_part(CS,ncell,sol,gas,r);
	for(int i=0;i<ncell;i++)
	{
		(CS+i)->react=1;
	}
//	initial_to_znd(CS,ncell,sol,gas,r);
//	exchange_to_hydro_part(CS,ncell,gas,r);
//	The function below is the initial condition fix, normally it is not running.
//	place_shock_condition(CS,0, GC,CS,ncell,sol,gas);  
//	exchange_to_chemical_part(CS,ncell,gas,r);
	for(int i=0;i<ncell;i++)
	{
		sim[i].setInitialTime(tm);
	}
//	initial_to_znd(CS,ncell,gas,r);
//	place_shock_condition(CS, 63, CS+63,CS+64,ncell,gas); 
	{
		//initial_to_chemical_part(CS2,ncell,sol2,gas2,r2);
//		initial_to_znd(CS,ncell,sol,gas,r);
		//exchange_to_hydro_part(CS2,ncell,gas2,r2);
	//	The function below is the initial condition fix, normally it is not running.
//		place_shock_condition(CS2,0, GC2,CS2,ncell,sol2,gas2);  
//		exchange_to_chemical_part(CS2,ncell,gas2,r2);
	//	initial_to_znd(CS,ncell,gas,r);
	//	place_shock_condition(CS, 63, CS+63,CS+64,ncell,gas); 
		for(int i=0;i<ncell;i++)
		{
			sim2[i].setInitialTime(tm);
		}
	}

//	std::cout<<(CS)->rj<<" "<<(CS+1)->lj<<" "<<std::endl;
	plot(CS,JS,ncell,gas,picnumber,0);
	
	//one step time marching. The Hydrodynamic solver and chemical solver are using two different memory sets
	//Hydro solver would only store the thermodynamic data.
	//Chemical solver would store both thermodynamic data and chemical component data.
//	double final_time=endtime;



	for(int ploti=1;ploti<=total_image;ploti++)
	{
		
		double endtime=final_time*ploti/total_image;
		//double endtime=final_time;
		while(tm<endtime)
		{
			//copy thermaldynamic properties from chemical solver's memory to hydrodynamics solver's memory
			//exchange_to_hydro_part(CS,ncell,gas,r);
			//hydrodynamic step activate (the marching time would be computed at this part)
			
			dt=march(CS,JS,GC,CS2,JS2,GC2,ncell,tm,endtime,gas,r,gas2,r2,sim,sim2);
			//dt=endtime;
			cout<<"dt="<<dt<<" tm="<<tm<< "   "<< tm*100/final_time<<" %"<<endl;
			//copy thermaldynamic properties from hydrodynamics solver's memory to chemical solver's memory
			
			//exchange_to_chemical_part(CS,ncell,gas,r);
			//std::cout<<"dt = "<<dt<<" T="<<r[2]->temperature()<<std::endl;
			//Initialize the chemical reactor, otherwise reactor would refuse the new thermodynamic properties update from last commend
			#pragma omp parallel for
			for(int i=0;i<ncell;i++)
			{
				sim[i].setInitialTime(tc);
			}

			//update the current time
			tm = tm+dt;
			//chemical step activate
			if(pic%1==0)
			{
				#pragma omp parallel for
				for(int i=0;i<ncell;i++)
				{
					if((CS+i)->react==1)
					{
						//sim[i].advance(tm);
					}
 				}
				tc=tm;
			}
			//Below this part is the plot stuff, not related to result
			//If you only wanna check the computation process, this part is not important
			

			
			//std::cout<<" dt = "<<dt<<" T="<<r[2]->temperature()<<std::endl;
			//if(pic%4==0)

			cout<<"dt="<<dt<<" tm="<<tm<<endl;
			
		
		}//end of a step
		pic++;
		picnumber=pic/1;
		plot(CS,JS,ncell,gas,picnumber,tm);
		tp=tm;	
    }
    clock_t t1 = clock();        // save end time
    
    return 0;
}


int main()
{
	globalsettings_initial();
    try {
		struct timeval t1,t2;
		double timeuse;	
		gettimeofday(&t1,NULL);
        int retn = kinetics();
        gettimeofday(&t2,NULL);
        timeuse = (t2.tv_sec - t1.tv_sec)*1000.0 + (t2.tv_usec - t1.tv_usec)/1000.0;
        int day,hour,min,second,timesec;
        {
			timesec=t2.tv_sec - t1.tv_sec;
			day=timesec/86400;
			timesec=timesec%86400;
			hour=timesec/3600;
			timesec=timesec%3600;
			min=timesec/60;
			second=timesec%60;
		}
		if(day!=0)
		{
			std::cout << day << "day(s)  ";
		}
		if(hour!=0)
		{
			std::cout << hour << "hour(s)  ";
		}
		if(min!=0)
		{
			std::cout << min << "min(s)  ";
		}
		if(second!=0)
		{
			std::cout << second << "second(s)  ";
		}

        //std::cout << timeuse << std::endl;
        appdelete();
        return retn;
    }
    // handle exceptions thrown by Cantera
    catch (CanteraError& err) {
        std::cout << err.what() << std::endl;
        cout << " terminating... " << endl;
        appdelete();
        return -1;
    }
}


