#include "initialconditionfix.h"
void place_shock_condition(cell* CS, int place_lcell, cell* left,cell* right,int nc,shared_ptr<Cantera::Solution>* sol, shared_ptr<Cantera::ThermoPhase>* gas)
{
	cell* TCS;
	cell* TCS2;
	cell* TGC;
	cell* TGC2;
	cell* TCE;
	join* TJS;
	join* TJS2;
	bool righttoward=true;
	
	if(left->prim[2]<right->prim[2])
	{
		left=CS+place_lcell+1;
		right=CS+place_lcell;
		righttoward=false;
		if(place_lcell<nc-7)
		{
			std::cout<<left->prim[2]<<" "<<right->prim[2]<<std::endl;
			return;
		}
	}
	
	std::cout<<"Step1"<<std::endl;
	// create an ideal gas mixture that corresponds to GRI-Mech
    // 3.0
	// set up gas models
	shared_ptr<Cantera::Solution> solt[nctest];
	shared_ptr<Cantera::ThermoPhase> gast[nctest];
	
	// create reactors
    Cantera::Reactor rt[nctest];

    // create a container object to run the simulation
    // and add the reactor to it
    Cantera::ReactorNet simt[nctest];
    
    chemicalinitiation(solt,gast,rt,simt,nctest);
    
    	
	double trini[nctest];
	double tuini[nctest];
	double tpini[nctest];
	double tgaini[nctest];
	std::cout<<"Step2"<<std::endl;
	for(int i=0;i<nctest;i++)
	{
		trini[i]=1.0/right->prim[0];
		tuini[i]=right->prim[1];
		tpini[i]=right->prim[2];
		tgaini[i]=right->ga;
	}
	TCS=initiate(nctest,trini,tuini,tpini,tgaini);
	TCS2=initiate(nctest,trini,tuini,tpini,tgaini);
	TGC=TCS->lc;
	TGC2=TCS2->lc;
	TGC->cellcopy(left);
	TGC2->cellcopy(left);

	TJS=TCS->lj;
	TJS2=TCS2->lj;
	TCE=TCS+nctest-3;
	    
	double dt,tm,dfpressure,valve,pvalve;
	double dfdensity=left->prim[2]-right->prim[2];
	valve=0.00001;
	pvalve=valve*dfdensity+right->prim[2];
	initial_to_chemical_part(TCS,nctest,solt,gast,rt);
	for(int i=0;i<nctest;i++)
	{
		simt[i].setInitialTime(tm);
	}
	std::cout<<"Step3"<<std::endl;
	while(TCE->prim[2]<=pvalve)
	{
		
		exchange_to_hydro_part(TCS,nctest,gast,rt);
		//dt=marchtest(TCS,TJS,TGC,TCS2,TJS2,TGC2,left,nctest,tm,endtime,gas,r,gas2,r2,sim,sim2);
		dt=marchtest(TCS,TJS,TGC,TCS2,TJS2,TGC2,left,nctest);

		//exchange thermaldynamic properties
        exchange_to_chemical_part(TCS,nctest,gast,rt);
		for(int i=0;i<nctest;i++)
		{
			simt[i].setInitialTime(tm);
		}
        tm = tm+dt;

        for(int i=0;i<nctest;i++)
		{
			//simt[i].advance(tm);
		}		

        std::cout<<"TEST t="<<tm<<" "<<TCE->prim[2]<<" "<<pvalve<<std::endl;
	}
	cell* TCP;
	cell* CP;
	shared_ptr<Cantera::ThermoPhase>* TgasP;
	shared_ptr<Cantera::ThermoPhase>* gasP;
	if(righttoward==true)
	{
		for(int i=0;i<9;i++)
		{
			gasP=gas+place_lcell+i;
			TgasP=gast+nctest-3-7+i;
			CP=CS+place_lcell+i;
			TCP=TCS+nctest-3-7+i;
			CP->cell::cellcopy(TCP);
//			*gasP = *TgasP;
		}
	}
	else
	{
		for(int i=0;i<9;i++)
		{
			gasP=gas+place_lcell+i;
			TgasP=gast+nctest-3-7+i;
			CP=CS+place_lcell+1-i;
			TCP=TCS+nctest-3-7+i;
			CP->cell::cellcopy(TCP);
//			*gasP = *TgasP;
		}
	}
	for(int i=0;i<nctest;i++)
	{
		TCE=TCS+i;
		if(i==0)
		{
			std::cout<<"GS "<<TCE->lc->ga<<" "<<1.0/TCE->lc->prim[0]<<" "<<TCE->lc->prim[1]<<" "<<TCE->lc->prim[2]<<" "<<TCE->lc->prim[3]<<std::endl;
		}
		std::cout<<"CS "<<TCE->ga<<" "<<1.0/TCE->prim[0]<<" "<<TCE->prim[1]<<" "<<TCE->prim[2]<<" "<<TCE->prim[3]<<std::endl;
	}
}
