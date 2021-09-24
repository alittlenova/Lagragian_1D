/////////////////////////////////////////////////////////////
//
//  zero-dimensional kinetics example program
//
//  copyright California Institute of Technology 2002
//
/////////////////////////////////////////////////////////////

#include "reaction.h"

using namespace Cantera;
using std::cout;
using std::endl;
void chemicalinitiation(shared_ptr<Cantera::Solution>* sol, shared_ptr<Cantera::ThermoPhase>* gas,Reactor* r,ReactorNet* sim,int cell)
{
	for(int i=0;i<cell;i++)
    {
		//set the gas class
		//sol[i] = newSolution("h2o2.yaml", "ohmech", "None");
		//sol[i] = newSolution("h2o2mech_Lietal_2003.yaml", "gas", "None");
		
		
		//gas[i] = new IdealGasMix("sandiego2014_WT.yaml", "gas");
		sol[i] = newSolution("sandiego2014_WT_Short_Parameter_V0.yaml", "gas", "None");
		gas[i] = sol[i]->thermo();
	}
	
	for(int i=0;i<cell;i++)
    {
		//set the reactor
		//r[i] = new Reactor;
	}
	
	for(int i=0;i<cell;i++)
    {
		// set the state
		//gas[i]->setState_TPX(1100, OneAtm, "H2:2.0, O2:1.0");
		gas[i]->setState_TPX(300, 6200, "C2H4:1.0, O2:3"); //For C2H4
		//gas[i]->setState_TPX(300, 0.07305*OneAtm, "H2:2.0, O2:1.0"); //For H2O2 950K 1 atm
		//gas[i]->setState_TPX(300, 0.05818*OneAtm, "H2:2.0, O2:1.0"); //For H2O2 1100K 1 atm
		//gas[i]->setState_TPX(1100, OneAtm, "H2:2.0, O2:1.0"); //For H2O2 1100K 1 atm
		//gas[i]->setState_TPX(293, OneAtm, "N2:3.76");
		// insert gas
		r[i].insert(sol[i]);
		//the commend below separate the energy computation for sundial integrator. 
		//Because it is a private component, we cannot change it outside.
		//but we could turn it off.
		r[0].setEnergy(0);
		//r[i]->setEnergy(0);
	}		
	
	for(int i=0;i<cell;i++)
    {
		sim[i].addReactor(r[i]);
	}
	
}

