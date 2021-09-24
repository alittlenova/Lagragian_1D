#include "exchanging.h"
#include <iostream>
//internal energy e Exchange
/*void exchange_to_hydro_part(cell *CS, int nc,shared_ptr<Cantera::ThermoPhase>** gass,Cantera::Reactor** reac)
{
	//CS is the first cell in hydro solver
	//gass is the first cell in chemical solver 
	cell *CP;
	Cantera::IdealGasMix* gas;
	Cantera::Reactor* react;
	for(int i=0;i<nc;i++)
	{
		CP=CS+i;
		gas=gass[i];
		react=reac[i];
//		gas->setDensity(react->density());
//		gas->setTemperature(react->temperature());
		CP->prim[2]=gas->pressure(); 
		CP->prim[4]=gas->intEnergy_mole()/gas->meanMolecularWeight() ;
		CP->prim[3]=CP->prim[4]+(0.5*CP->prim[1]*CP->prim[1]);
		CP->ga=gas->cp_mole()/gas->cv_mole();
	}
}*/
/*
void exchange_to_hydro_part(cell *CS, int nc,shared_ptr<Cantera::ThermoPhase>* gass,Cantera::Reactor* reac)
{
	//CS is the first cell in hydro solver
	//gass is the first cell in chemical solver 
	double T;
	cell *CP;
	shared_ptr<Cantera::ThermoPhase>* gas;
	Cantera::Reactor* react;
	for(int i=0;i<nc;i++)
	{
		CP=CS+i;
		gas=&gass[i];
		react=&reac[i];
//		gas->setDensity(react->density());
//		gas->setTemperature(react->temperature());
		T=(*gas)->temperature();
//		T=CP->prim[2]*(*gas)->meanMolecularWeight()*CP->prim[0]/8314.4621;
		CP->prim[2]=T*8314.4621/((*gas)->meanMolecularWeight()*CP->prim[0]);
//		CP->prim[4]=(*gas)->intEnergy_mole()/(*gas)->meanMolecularWeight() ;
		CP->prim[4]=(CP->prim[2])*(CP->prim[0])/((CP->ga)-1);
		CP->prim[3]=CP->prim[4]+(0.5*CP->prim[1]*CP->prim[1]);
		CP->ga=(*gas)->cp_mole()/(*gas)->cv_mole();
		if((CP->ga)<1)
		{
			std::cout<<"Gamma Error  Gamma="<<CP->ga<<"  on Cell "<<i<<std::endl;
			std::cout<<"P="<<gas[0]->pressure()<<"  T="<<gas[0]->temperature()<<std::endl;
			std::cout<<"cp="<<gas[0]->cp_mole()<<"  cv="<<gas[0]->cv_mole()<<std::endl;
		}
	}
}*/
//temperature T exchange

void exchange_to_hydro_part(cell *CS, int nc,shared_ptr<Cantera::ThermoPhase>* gass,Cantera::Reactor* reac)
{
	//CS is the first cell in hydro solver
	//gass is the first cell in chemical solver 
	cell *CP;
	shared_ptr<Cantera::ThermoPhase>* gas;
	Cantera::Reactor* react;
	double E0,E1,T,dedt,De,intE;
	for(int i=0;i<nc;i++)
	{
		if((CS+i)->react==1)
		{
			CP=CS+i;
			gas=&gass[i];
			react=&reac[i];
	//		gas->setDensity(react->density());
	//		gas->setTemperature(react->temperature());
			CP->prim[2]=(*gas)->pressure();
			
	//		std::cout<<(*gas)->pressure()<<std::endl;
	//		CP->prim[4]=(CP->prim[2])*(CP->prim[0])/((CP->ga)-1);
			
			CP->ga=(*gas)->cp_mole()/(*gas)->cv_mole();
			T=(*gas)->temperature();
			if(T<999 ||T>1001)
			{
				CP->prim[5]=(*gas)->intEnergy_mole()/(*gas)->meanMolecularWeight() ;
			}
			else
			{
				(*gas)->setTemperature(999);
				E0=(*gas)->intEnergy_mole();
				(*gas)->setTemperature(1001);
				E1=(*gas)->intEnergy_mole();
				dedt=(E1-E0)/2.0;
				intE=((T-999)*dedt)+E0;
				CP->prim[5]=intE/(*gas)->meanMolecularWeight() ;
				(*gas)->setTemperature(T);
			}
			
			CP->prim[4]=(1/((CP->ga)-1))*(*gas)->pressure()/(*gas)->density();
			CP->prim[3]=CP->prim[4]+(0.5*CP->prim[1]*CP->prim[1]);

			CP->prim[6]=CP->prim[4]-CP->prim[5];
			//std::cout<<"prim4="<<(*gas)->intEnergy_mole()<<"  on Cell "<<i<<std::endl;
			if((CP->ga)<1)
			{
				//std::cout<<"Gamma Error  Gamma="<<CP->ga<<"  on Cell "<<i<<std::endl;
				//std::cout<<"P="<<gas[0]->pressure()<<"  T="<<gas[0]->temperature()<<std::endl;
				//std::cout<<"cp="<<gas[0]->cp_mole()<<"  cv="<<gas[0]->cv_mole()<<std::endl;
			}
		}
	}
}

//internal energy e exchange
/*void exchange_to_chemical_part(cell *CS, int nc,Cantera::IdealGasMix** gass,Cantera::Reactor** reac)
{
	//CS is the first cell in hydro solver
	//gass is the first cell in chemical solver
	cell *CP;
	Cantera::IdealGasMix* gas;
	Cantera::Reactor* react;
	for(int i=0;i<nc;i++)
	{
		CP=CS+i;
		gas=gass[i];
		react=reac[i];
		if(CP->prim[0]<0)
		{
			std::cout<<i<<std::endl;
		}
		gas->setDensity(1/CP->prim[0]);
		CP->prim[4]=CP->prim[3]-(0.5*CP->prim[1]*CP->prim[1]);
		Try_temperature(gas, CP->prim[4]* gas->meanMolecularWeight());
		react->setThermoMgr(*gass[i]);
	}
}*/
//temperature T exchange
void exchange_to_chemical_part(cell *CS, int nc,shared_ptr<Cantera::ThermoPhase>* gass,Cantera::Reactor* reac)
{
	//CS is the first cell in hydro solver
	//gass is the first cell in chemical solver
	const double cutoff=4980;
	cell *CP;
	double T;
	shared_ptr<Cantera::ThermoPhase>* gas;
	Cantera::Reactor* react;
	for(int i=0;i<nc;i++)
	{
		if((CS+i)->react==1)
		{
			CP=CS+i;
			gas=&gass[i];
			react=&reac[i];
			if(CP->prim[0]<0)
			{
				std::cout<<i<<std::endl;
			}
			//std::cout<<"Rho="<<1/CP->prim[0]<<std::endl;
			gas[0]->setDensity(1/CP->prim[0]);
			CP->prim[4]=CP->prim[3]-(0.5*CP->prim[1]*CP->prim[1]);
			CP->prim[5]=CP->prim[4]-CP->prim[6];
			//T=CP->prim[2]*(*gas)->meanMolecularWeight()*CP->prim[0]/8314.4621;
			//std::cout<<"CP prim 4 ="<<CP->prim[4]<<std::endl;
			//(*gas)->setTemperature(T);
			if(gas[0]->temperature()>cutoff)
			{
				(CS+i)->react=0;
			}
			Try_temperature(gas, CP->prim[5]* (*gas)->meanMolecularWeight());
			react->setThermoMgr(*gas[0]);
		}
	}
}
void Try_temperature(shared_ptr<Cantera::ThermoPhase>* gas,double intE)
{
	//Use Newton's method to compute the right temperature for certain state with certain internal energy
	double E0,E1,T,dedt,De;
	E0=gas[0]->intEnergy_mole();
	T=gas[0]->temperature();
	gas[0]->setTemperature(T+0.001);
	E1=gas[0]->intEnergy_mole();
	dedt=(E1-E0)/0.001;
	De=intE-E0;
	const double TTrans=1000;
	const double FitRange=2;
	
	while(De<-0.0001||De>0.0001)
	{
		if(T<TTrans-FitRange ||T>TTrans+FitRange)
		{
			if(T+(De/dedt)<0.01)
			{
				gas[0]->setTemperature(0.01);
			}
			else
			{
				gas[0]->setTemperature(T+(De/dedt));
			}
			E0=gas[0]->intEnergy_mole();
			T=gas[0]->temperature();
			gas[0]->setTemperature(T+(0.001*T));
			E1=gas[0]->intEnergy_mole();
			dedt=(E1-E0)/(0.001*T);
			De=intE-E0;
		}
		else
		{
			gas[0]->setTemperature(TTrans-FitRange);
			E0=gas[0]->intEnergy_mole();
			gas[0]->setTemperature(TTrans+FitRange);
			E1=gas[0]->intEnergy_mole();
			dedt=(E1-E0)/(2.0*FitRange);
			De=intE-E0;
			gas[0]->setTemperature(TTrans-FitRange+(De/dedt));
			T=gas[0]->temperature();
			De=0;
		}
	}
	gas[0]->setTemperature(T+(De/dedt));
}
