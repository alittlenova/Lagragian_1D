#include "initialcondition.h"
void initial_to_chemical_part(cell *CS, int nc,shared_ptr<Cantera::Solution>* sols,shared_ptr<Cantera::ThermoPhase>* gass,Cantera::Reactor* reac)
{
	cell *CP;
	shared_ptr<Cantera::ThermoPhase>* gas;
	Cantera::Reactor* react;
	//for(int i=0;i<nc;i++)
	int i=0;
	{
		CP=CS+i;
		gas=&gass[i]; 
		react=reac;
		(*gas)->setDensity(1/CP->prim[0]);
		(*gas)->setTemperature(CP->prim[2] * ((*gas)->meanMolecularWeight()) * CP->prim[0]/ (Cantera::GasConstant));
		react->insert(sols[i]);
	}
}
void initial_to_znd(cell *CS, int nc,shared_ptr<Cantera::Solution>* sols,shared_ptr<Cantera::ThermoPhase>* gass,Cantera::Reactor** reac)
{
	double **hydrodata;
	double **chemdata;
	int hydrom,chemm,hydron,chemn,target;
	int preprocess=1;
	int skip1=50;//tail skip
	int skip2=10;//Head skip
	double left_x,normaled_speed;
	double shock_speed=2242.3;
	cell* CP;
	cell* CP1;
	shared_ptr<Cantera::ThermoPhase>* gas;
	Cantera::Reactor* react;
	hydrodata=dataread("./H2.plt", &hydrom, &hydron);
	chemdata=dataread("./H2_sp.plt", &chemm, &chemn);
	doublereal chemdata_pure[chemn-2];
	std::cout<<chemm<<" "<<hydrom<<std::endl;
	const int d=2;
	for(int i=0;i<(hydrom-skip1-skip2)/d;i++)
	{
		CP=CS+i;
		gas=&gass[i];
		react=reac[i];
		if(i<preprocess)
		{
			target=hydrom-skip1;
			left_x=hydrodata[target][1];
			normaled_speed= shock_speed-hydrodata[target][7];
			
			for(int j=2;j<chemn;j++)
			{
				chemdata_pure[j-2]=chemdata[target][j];
			}
			(*gas)->setMassFractions(chemdata_pure);
			(*gas)->setDensity(hydrodata[target][5]);
			(*gas)->setTemperature(hydrodata[target][4]);
			
			CP->prim[0]=1/hydrodata[target][5];
			CP->prim[1]=normaled_speed;
			CP->prim[2]=hydrodata[target][3];
			CP->ga=hydrodata[target][9];
			CP->lj->x=left_x-hydrodata[target][1];
			CP->dx=hydrodata[target][1]-hydrodata[target-d][1];
			CP->x=CP->lj->x+(CP->dx/2.0);
			CP->rj->x=CP->x+(CP->dx/2.0);
			CP->dphi=(CP->dx)/(CP->prim[0]);
		}
		else
		{
			target=hydrom-skip1-(d*i);		
			normaled_speed= shock_speed-hydrodata[target][7];
			
			for(int j=2;j<chemn;j++)
			{
				chemdata_pure[j-2]=chemdata[target][j];
			}
			(*gas)->setMassFractions(chemdata_pure);
			(*gas)->setDensity(hydrodata[target][5]);
			(*gas)->setTemperature(hydrodata[target][4]);
			
			CP->prim[0]=1/hydrodata[target][5];
			CP->prim[1]=normaled_speed;
			CP->prim[2]=hydrodata[target][3];
			CP->ga=hydrodata[target][9];
			CP->x=left_x-hydrodata[target][1];
			CP->rj->x=CP->x+(0.5*(hydrodata[target][1]-hydrodata[target-d][1]));
			CP->dx=(CP->rj->x)-(CP->lj->x);
			CP->dphi=(CP->dx)/(CP->prim[0]);
			CP->phi=(CP->lc->phi)+0.5*(CP->lc->dphi)+0.5*(CP->dphi);
		}
		react->insert(sols[i]);
	}
	CP1=CP;
	for(int i=(hydrom-skip1-skip2)/d;i<nc;i++)
	{
		CP=CS+i;
		gas=&gass[i];
		react=reac[i];
		CP->dphi=CP1->dphi;
		CP->phi=(CP->lc->phi)+0.5*(CP->lc->dphi)+0.5*(CP->dphi);
		CP->dx=(CP->dphi)*(CP->prim[0]);
		CP->x=(CP->lj->x)+0.5*(CP->dx);
		CP->rj->x=(CP->x)+0.5*(CP->dx);
	}
	std::cout<<(hydrom-skip1)/d<<std::endl;
	
}
