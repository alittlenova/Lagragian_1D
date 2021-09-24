#include<fstream>
#include<string>
#include<math.h>
#include<iostream>
#include<stdlib.h>
#include "ReimannSolver.h"
 

double initialiter(double cl,double cr,double al,double ar,double pl,double pr,double ul,double ur,double gammal, double gammar)
{
	//find iteration initial state
	double sigma,z,ulbar,urbar,ustar0;
	if(pl>=pr)
	{
		sigma=gammal;
	}
	else
	{
		sigma=gammar;
	}
	z=((gammal-1.0)/(gammar-1.0))*(ar/al)* pow( (pl/pr), ((sigma-1.0)/(2.0*sigma)) );
	ulbar=ul+al*(2.0/(gammal-1.0));
	urbar=ur-ar*(2.0/(gammar-1.0));
	ustar0=(ulbar*z+urbar)/(1.0+z);
	return ustar0;
}
double iteration(double cl,double cr,double al,double ar,double pl,double pr,double ul,double ur,double gammal, double gammar)
{
	//This function Will solve the reimann problem and will only return ustar
	double ustar,plstar,prstar,plstarprime,prstarprime,wl,wr,alstar,arstar,errors,trytime;
	ustar=initialiter(cl,cr,al,ar,pl,pr,ul,ur,gammal,gammar);
	if(ustar<ul)
	{
		wl=((gammal+1.0)/4.0)*((ustar-ul)/al)-pow((((gammal+1.0)/4.0)*((gammal+1.0)/4.0)*((ustar-ul)/al)*((ustar-ul)/al)+1.0),0.5);
		plstar=pl+(cl*(ustar-ul)*wl);
		plstarprime=(2.0*cl*pow(wl,3.0))/(wl*wl+1.0);
	}
	else
	{
		alstar=al-(gammal-1.0)*(ustar-ul)/2.0;
		plstar=pl*pow((alstar/al),2.0*gammal/(gammal-1.0));
		plstarprime=-gammal*plstar/alstar;
	}
	if(ustar>=ur)
	{
		wr=((gammar+1.0)/4.0)*((ustar-ur)/ar)+pow((((gammar+1.0)/4.0)*((gammar+1.0)/4.0)*((ustar-ur)/ar)*((ustar-ur)/ar)+1.0),0.5);
		prstar=pr+(cr*(ustar-ur)*wr);
		prstarprime=(2.0*cr*pow(wr,3.0))/(wr*wr+1.0);
	}
	else
	{
		arstar=ar+(gammar-1.0)*(ustar-ur)/2.0;
		prstar=pr*pow((arstar/ar),2.0*gammar/(gammar-1.0));
		prstarprime=gammar*prstar/arstar;
	}
	errors=plstar-prstar;
	trytime=0;
	
	while(((errors>0.001)||(errors<-0.001)) && (trytime<50))
	{
		ustar=ustar-( (plstar-prstar)/(plstarprime-prstarprime) );
		
	if(ustar<ul)
	{
		wl=((gammal+1.0)/4.0)*((ustar-ul)/al)-pow((((gammal+1.0)/4.0)*((gammal+1.0)/4.0)*((ustar-ul)/al)*((ustar-ul)/al)+1.0),0.5);
		plstar=pl+(cl*(ustar-ul)*wl);
		plstarprime=(2.0*cl*pow(wl,3.0))/(wl*wl+1.0);
	}
	else
	{
		alstar=al-(gammal-1.0)*(ustar-ul)/2.0;
		plstar=pl*pow((alstar/al),(2.0*gammal/(gammal-1.0)));
		plstarprime=-gammal*plstar/alstar;
	}
	if(ustar>ur)
	{
		wr=((gammar+1.0)/4.0)*((ustar-ur)/ar)+pow((((gammar+1.0)/4.0)*((gammar+1.0)/4.0)*((ustar-ur)/ar)*((ustar-ur)/ar)+1.0),0.5);
		prstar=pr+(cr*(ustar-ur)*wr);
		prstarprime=(2.0*cr*pow(wr,3.0))/(wr*wr+1.0);
	}
	else
	{
		arstar=ar+(gammar-1.0)*(ustar-ur)/2.0;
		prstar=pr*pow((arstar/ar),2.0*gammar/(gammar-1.0));
		prstarprime=gammar*prstar/arstar;
	}
		errors=prstar-plstar;
		trytime++;
	}
	if(trytime==50)
	{
		ustar=initialiter(cl,cr,al,ar,pl,pr,ul,ur,gammal,gammar);
	}
	return ustar;
}
void fluidcons(double *rl,double *rr,double *ul,double *ur,double *pl,double *pr,double *f0,double *f1,double *f2,double gammal,double gammar)
{
	//use iteration method to find the flux on the cell interface and only return f1 f2 f3(the three premitive variables rho p u by sequence) 
	//This function only deal the two local cells but not global
	//the result of this function gives the flux direcly
	double ustar,astar,pstar,rstar,al,ar,cl,cr,lr,wl,wr;
	al=pow(gammal*(*pl)/(*rl), 0.5);
	ar=pow(gammar*(*pr)/(*rr), 0.5);
	cl=gammal*(*pl)/al;
	cr=gammar*(*pr)/ar;
	ustar=iteration(cl,cr,al,ar,*pl,*pr,*ul,*ur,gammal,gammar);
	if(ustar>*ul)
	{
		astar=al-((gammal-1)/2)*(ustar-*ul);
		pstar=*pl*pow( (astar/al), ((2*gammal)/(gammal-1)) );
		rstar=gammal*pstar/(astar*astar);
	}
	else if(ustar<=*ur)
	{
		astar=ar+((gammar-1)/2)*(ustar-*ur);
		pstar=*pr*pow( (astar/ar), ((2*gammar)/(gammar-1)) );
		rstar=gammar*pstar/(astar*astar);
	}
	else if(*ul>*ur)
	{
		if(2*ustar-*ul-*ur<0)
		{
			wr=((gammar+1.0)/4.0)*((ustar-*ur)/ar)+pow((((gammar+1.0)/4.0)*((gammar+1.0)/4.0)*((ustar-*ur)/ar)*((ustar-*ur)/ar)+1.0),0.5);
			pstar=*pr+(cr*(ustar-*ur)*wr);
		}
		else
		{
			wl=((gammal+1.0)/4.0)*((ustar-*ul)/al)-pow((((gammal+1.0)/4.0)*((gammal+1.0)/4.0)*((ustar-*ul)/al)*((ustar-*ul)/al)+1.0),0.5);
			pstar=*pl+(cl*(ustar-*ul)*wl);
		}
	}
	else
	{
		if(2*ustar-*ul-*ur<0)
		{
			astar=al-((gammal-1)/2)*(ustar-*ul);
			pstar=*pl*pow( (astar/al), ((2*gammal)/(gammal-1)) );
		}
		else
		{
			astar=ar+((gammar-1)/2)*(ustar-*ur);
			pstar=*pr*pow( (astar/ar), ((2*gammar)/(gammar-1)) );
		}
	}
/*	*f0=rstar;
	*f1=ustar;
	*f2=pstar*ustar;*/
	*f0=-ustar; 
	*f1=pstar;
	*f2=pstar*ustar;
}
 

double maxdt(cell *CS,int number)
{
	cell *CI;
	double local,global,global1;
	global=10000;
	for(int i=0;i<number;i++)
	{
		CI=CS+i;
		local=(CI->dx)/pow(((CI->ga)*CI->prim[2]*std::max(CI->prim[0],std::max(CI->lc->prim[0],CI->rc->prim[0]))),0.5);
		if(local<global)
		{
			global=local;
			//std::cout<<"global="<<global<<" i="<<i<<" CI->prim[0]="<<CI->prim[0]<<" CI->lc->prim[0]="<<CI->lc->prim[0]<<" CI->rc->prim[0]="<<CI->rc->prim[0]<<std::endl;
			//std::cout<<"fl="<< CI->lj->f[0]<<"  fr="<< CI->rj->f[0]<<std::endl;
		}
//		global1=(CI->dphi)/global;
/*		if(local>-1)
		{
		}
		else
		{
			std::cout<<CI->phi<<std::endl;
		}*/
	}
	return global;
}
double march(cell *CS,join *JS,cell *GS,cell *CS2,join *JS2,cell *GS2,int nc, double time,double endtime,shared_ptr<Cantera::ThermoPhase>* gass,Cantera::Reactor* reac,shared_ptr<Cantera::ThermoPhase>* gass2,Cantera::Reactor* reac2,Cantera::ReactorNet* sim,Cantera::ReactorNet* sim2)
{
	double maxt;
	double up=funcp(time);
	
	cell *CP,*CP1;
	join *JP;
	CP=GS;
	
	CP=GS+1;
	CP1=CS+nc-1;
	CP->cellcopy(CP1);

	double rl,rr,ul,ur,pl,pr,f0,f1,f2,gammal,gammar;
	for(int c=1;c<=nc;c++)
	{
		JP=JS+c;
		rl=1/(JP->lc->prim[0]);
		rr=1/(JP->rc->prim[0]);
		ul=JP->lc->prim[1];
		ur=JP->rc->prim[1];
		pl=JP->lc->prim[2];
		pr=JP->rc->prim[2];
		gammal=JP->lc->ga;
		gammar=JP->rc->ga;
		//rl=1.0/vl;
		//rr=1.0/vr;
		//fluidcons(&rl,&rr,&ul,&ur,&pl,&pr,&f0,&f1,&f2,gammal,gammar);
		HLLE(&(JP->lc->prim[0]),&(JP->lc->prim[1]),&(JP->lc->prim[3]),&(JP->lc->prim[2]),&(JP->rc->prim[0]),&(JP->rc->prim[1]),&(JP->rc->prim[3]),&(JP->rc->prim[2]),&f0,&f1,&f2,&gammal,&gammar);
		JP->f[0]=f0;
		JP->f[1]=f1;
		JP->f[2]=f2;
	}
	maxt=0.69*maxdt(CS,nc);
	
	double plghost=shockjumpp(up,CS->prim[1],CS->prim[2],1.0/CS->prim[0],CS->ga);
	double rlghost=shockjumpr(up,CS->prim[1],CS->prim[2],1.0/CS->prim[0],CS->ga);
	double alghost=pow(plghost*(CS->ga)/rlghost, 0.5);
	double clghost=alghost*rlghost;
	//std::cout<<"alghost="<<alghost<<"  rlghost="<<rlghost<<"  clghost="<<clghost<<" CS->dphi="<<CS->dphi<<std::endl;
	double maxt1=0.69*(CS->dphi)/clghost;
	
	if(maxt>maxt1)
	{
		maxt=maxt1;
	}
	if(time+maxt>endtime)
	{
		maxt=endtime-time;
	}
	for(int c=1;c<=nc;c++)
	{
		JP=JS+c;
		JP->x=(JP->x)-(JP->f[0])*maxt;
	}
	for(int c=0;c<nc;c++)
	{
		CP=CS+c;
		
		if(c==0)
		{
			CP->lj->f[0]=(-up);
			if(up>=CP->prim[1])
			{
				CP->lj->f[1]=shockjumpp(up,CP->prim[1],CP->prim[2],1/(CP->prim[0]),CP->ga);
			}
			else
			{
				CP->lj->f[1]=expansionjumpp(up,CP->prim[1],CP->prim[2],1/(CP->prim[0]),CP->ga);
			}
			CP->prim[0]=CP->prim[0]+((CP->lj->f[0])-(CP->rj->f[0]))*maxt/(CP->dphi);
			CP->prim[1]=CP->prim[1]+((CP->lj->f[1])-(CP->rj->f[1]))*maxt/(CP->dphi);
			CP->prim[3]=CP->prim[3]+(-(CP->rj->f[2])-((CP->lj->f[0])*(CP->lj->f[1])))*maxt/(CP->dphi);
			CP->prim[4]=CP->prim[3]-(0.5*CP->prim[1]*CP->prim[1]);
			CP->prim[2]=CP->prim[4]*((CP->ga)-1.0)/(CP->prim[0]);
			CP->dx=(CP->dphi)*(CP->prim[0]);
			CP->lj->x=CP->lj->x-(CP->lj->f[0])*maxt;
			CP->x=(CP->lj->x)+0.5*(CP->dx);
			if(CP->prim[0]>0)
			{}
			else
			{
				std::cout<<"Error at Marching Cell U"<<std::endl;
				std::cout<<"cell="<<c<<" v="<<CP->prim[0]<<" fl="<<CP->lj->f[0]<<" fr="<<CP->rj->f[0]<<std::endl;
				std::cout<<"cell="<<c<<" u="<<CP->prim[1]<<" fl="<<CP->lj->f[1]<<" fr="<<CP->rj->f[1]<<std::endl;
			}
		}
		else
		{
			CP->prim[0]=CP->prim[0]+((CP->lj->f[0])-(CP->rj->f[0]))*maxt/(CP->dphi);
			CP->prim[1]=CP->prim[1]+((CP->lj->f[1])-(CP->rj->f[1]))*maxt/(CP->dphi);
			CP->prim[3]=CP->prim[3]+((CP->lj->f[2])-(CP->rj->f[2]))*maxt/(CP->dphi);
			CP->prim[4]=CP->prim[3]-(0.5*CP->prim[1]*CP->prim[1]);
			CP->prim[2]=CP->prim[4]*((CP->ga)-1.0)/(CP->prim[0]);
			CP->dx=(CP->dphi)*(CP->prim[0]);
			CP->x=(CP->lj->x)+0.5*(CP->dx);
			if(CP->prim[0]>0)
			{}
			else
			{
				std::cout<<"Error at Marching Cell U"<<std::endl;
				std::cout<<"cell="<<c<<" v="<<CP->prim[0]<<" fl="<<CP->lj->f[0]<<" fr="<<CP->rj->f[0]<<std::endl;
				std::cout<<"cell="<<c<<" u="<<CP->prim[1]<<" fl="<<CP->lj->f[1]<<" fr="<<CP->rj->f[1]<<std::endl;
				std::cout<<"cell="<<c<<" u="<<CP->prim[3]<<" fl="<<CP->lj->f[2]<<" fr="<<CP->rj->f[2]<<std::endl;
			}
		}
		
	}
	
	return maxt;
}
double marchtest(cell *CS,join *JS,cell *GS,cell *CS2,join *JS2,cell *GS2,cell *left,int nc)
{
	double maxt;
	double up=funcp(0);
	cell *CP,*CP1;
	join *JP;
	CP=GS;
	
	for(int i=0;i<3;i++)
	{
		CP->prim[i]=left->prim[i];
	}
	
	
	
	CP=GS+1;
	CP1=CS+nc-1;
	CP->cellcopy(CP1);
	double rl,rr,ul,ur,pl,pr,f0,f1,f2,gammal,gammar;
	for(int c=1;c<=nc;c++)
	{
		JP=JS+c;
		rl=1/(JP->lc->prim[0]);
		rr=1/(JP->rc->prim[0]);
		ul=JP->lc->prim[1];
		ur=JP->rc->prim[1];
		pl=JP->lc->prim[2];
		pr=JP->rc->prim[2];
		gammal=JP->lc->ga;
		gammar=JP->rc->ga;
		//fluidcons(&rl,&rr,&ul,&ur,&pl,&pr,&f0,&f1,&f2,gammal,gammar);
		HLLE(&(JP->lc->prim[0]),&(JP->lc->prim[1]),&(JP->lc->prim[3]),&(JP->lc->prim[2]),&(JP->rc->prim[0]),&(JP->rc->prim[1]),&(JP->rc->prim[3]),&(JP->rc->prim[2]),&f0,&f1,&f2,&gammal,&gammar);
		JP->f[0]=f0;
		JP->f[1]=f1;
		JP->f[2]=f2;
	}
	maxt=0.69*maxdt(CS,nc);

	double plghost=shockjumpp(up,CS->prim[1],CS->prim[2],1.0/CS->prim[0],CS->ga);
	double rlghost=shockjumpr(up,CS->prim[1],CS->prim[2],1.0/CS->prim[0],CS->ga);
	double alghost=pow(plghost*(CS->ga)/rlghost, 0.5);
	double clghost=alghost*rlghost;
	double maxt1=0.69*(CS->dphi)/clghost;
	if(maxt>maxt1)
	{
		maxt=maxt1;
	}
	for(int c=1;c<=nc;c++)
	{
		JP=JS+c;
		JP->x=(JP->x)-(JP->f[0])*maxt;
	}
	for(int c=0;c<nc;c++)
	{
		CP=CS+c;
		
		if(c==0)
		{
			CP->lj->f[0]=(-up);
			if(up>=CP->prim[1])
			{
				CP->lj->f[1]=shockjumpp(up,CP->prim[1],CP->prim[2],1/(CP->prim[0]),CP->ga);
			}
			else
			{
				CP->lj->f[1]=expansionjumpp(up,CP->prim[1],CP->prim[2],1/(CP->prim[0]),CP->ga);
			}
			CP->prim[0]=CP->prim[0]+((CP->lj->f[0])-(CP->rj->f[0]))*maxt/(CP->dphi);
			CP->prim[1]=CP->prim[1]+((CP->lj->f[1])-(CP->rj->f[1]))*maxt/(CP->dphi);
			CP->prim[3]=CP->prim[3]+(-(CP->rj->f[2])-((CP->lj->f[0])*(CP->lj->f[1])))*maxt/(CP->dphi);
			CP->prim[4]=CP->prim[3]-(0.5*CP->prim[1]*CP->prim[1]);
			CP->prim[2]=CP->prim[4]*((CP->ga)-1.0)/(CP->prim[0]);
			CP->dx=(CP->dphi)*(CP->prim[0]);
			CP->lj->x=CP->lj->x-(CP->lj->f[0])*maxt;
			CP->x=(CP->lj->x)+0.5*(CP->dx);
		}
		else
		{
			CP->prim[0]=CP->prim[0]+((CP->lj->f[0])-(CP->rj->f[0]))*maxt/(CP->dphi);
			CP->prim[1]=CP->prim[1]+((CP->lj->f[1])-(CP->rj->f[1]))*maxt/(CP->dphi);
			CP->prim[3]=CP->prim[3]+((CP->lj->f[2])-(CP->rj->f[2]))*maxt/(CP->dphi);
			CP->prim[4]=CP->prim[3]-(0.5*CP->prim[1]*CP->prim[1]);
			CP->prim[2]=CP->prim[4]*((CP->ga)-1.0)/(CP->prim[0]);
			CP->dx=(CP->dphi)*(CP->prim[0]);
			CP->x=(CP->lj->x)+0.5*(CP->dx);
		}
		
	}
	return maxt;
}
cell* initiate(int nc, double *rini, double *uini, double *pini, double *gaini)
{
	cell *CA = new cell[nc];
	join *JA = new join[nc+1];
	cell *GC = new cell[2];
	cell *C;
	join *J;
	for(int i=0;i<=nc;i++)
	{
		
		if(i==0)
		{
			C=CA+i;
			J=JA+i;
			C->lj=J;
			C->rj=J+1;
			C->lc=GC;
			C->rc=C+1;
			J->lc=GC;
			J->rc=C;
			J->rj=J+1;
		}
		else if(i==(nc-1))
		{
			C=CA+i;
			J=JA+i;
			C->lj=J;
			C->rj=J+1;
			C->lc=C-1;
			C->rc=GC+1;
			J->lc=C-1;
			J->rc=C;
			J->lj=J-1;
			J->rj=J+1;
		}
		else if(i==nc)
		{
			C=CA+i-1;
			J=JA+i;
			J->lc=C;
			J->rc=GC+1;
			J->lj=J-1;
		}
		else
		{
			C=CA+i;
			J=JA+i;
			C->lj=J;
			C->rj=J+1;
			C->lc=C-1;
			C->rc=C+1;
			J->lc=C-1;
			J->rc=C;
			J->lj=J-1;
			J->rj=J+1;
		}
	}
	double dp=(phir-phil)/nc;
	CA->x=0;

	for(int i=0;i<nc;i++)
	{
		C=CA+i;
		C->dphi=dp;
		C->phi=phil+i*dp;
			C->prim[0]=1/rini[i];
			C->prim[1]=uini[i];
			C->prim[2]=pini[i];
			C->ga=gaini[i];
		C->dx=(C->dphi)*(C->prim[0]);
		C->prim[3]=(((C->prim[2])*(C->prim[0]))/((C->ga)-1.0))+((C->prim[1])*(C->prim[1])*0.5);
	}
	for(int i=0;i<nc;i++)
	{
		C=CA+i;
		if(i==0)
		{
			C->x=xl+0.5*(C->dx);
		}
		else
		{
			C->x=(C->lc->x)+(0.5*(C->lc->dx))+(0.5*(C->dx));	
		}
	}
	for(int i=0;i<=nc;i++)
	{
		J=JA+i;
		if(i==0)
		{
			J->x=(J->rc->x)-(0.5*(J->rc->dx));
		}
		else
		{
			J->x=(J->lc->x)+(0.5*(J->lc->dx));
		}
	}
	return CA;
}
void plot(cell *CS,join *JS,int nc,shared_ptr<Cantera::ThermoPhase>* GS,int plotnum,double tm)
{
	cell *CI;
	join *JI;
	shared_ptr<Cantera::ThermoPhase>* GI;
	char t[5];
	char str[80];
	std::ofstream output1("./outputdata.csv");
	sprintf(t, "%d", plotnum);
	for(int i=0;i<nc;i++)
	{
		GI=GS+i;
		CI=CS+i;
		JI=JS+i;
		output1<<CI->x<<","<<(1/(CI->prim[0]))<<","<<CI->prim[1]<<","<<CI->prim[2]<<","<<CI->prim[3]<<","<<(*GI)->temperature()<<","<<CI->ga<<","<<CI->phi;
		//1 xlocation, 2 density, 3 velocity, 4 pressure, 5 internalE per mass, 6 Temperature, 7 Gamma, 8 Phi
		for(int j=0;j<((*GI)->nSpecies());j++)
		{
			output1<<","<<(*GI)->moleFraction(j);
		}
		output1<<","<<tm;
		//output1<<","<<CI->lj->f[0]<<","<<CI->lj->f[1]<<","<<CI->lj->f[2];
		output1<<std::endl;
	}
	output1.close();
	strcpy (str,"mv ./outputdata.csv ./output/");
	strcat (str,t);
	strcat (str,".csv");
		system(str);
}
double shockjumpp(double up,double ur,double pr,double rr,double gammar)
{
	double cr=pow(gammar*pr/rr, 0.5);
	double ut=(up-ur)/cr;
	double s1=(gammar+1.0)*ut;
	double M=(s1+pow(((s1*s1)+16),0.5))/4.0;
	double pratio=(2.0*gammar*(M*M-1.0))/(gammar+1);
	return (pratio*pr)+pr;
}
double expansionjumpp(double up,double ur,double pr,double rr,double gammar)
{
	double cr=pow(gammar*pr/rr, 0.5);
	double ud=up-ur;
	double pratio=pow(1+((gammar-1)*ud/(2*cr)),2*gammar/(gammar-1));
	return pratio*pr;
}
double shockjumpr(double up,double ur,double pr,double rr,double gammar)
{
	double cr=pow(gammar*pr/rr, 0.5);
	double ut=(up-ur)/cr;
	double s1=(gammar+1.0)*ut;
	double M=(s1+pow(((s1*s1)+16),0.5))/4.0;
	double rratio=(gammar+1.0)*M*M/((gammar-1)*M*M+2);
	return rr*rratio;
}
void copyreactor(shared_ptr<Cantera::ThermoPhase>* gass,Cantera::Reactor* reac,shared_ptr<Cantera::ThermoPhase>* gass2,Cantera::Reactor* reac2,double tm,int nc)
{
	shared_ptr<Cantera::ThermoPhase>* gas;
	Cantera::Reactor* react;
	shared_ptr<Cantera::ThermoPhase>* gas2;
	Cantera::Reactor* react2;
	double cahe[200];
	for(int i=0;i<nc;i++)
	{
		gas=&gass[i];
		react=&reac[i];
		gas2=&gass2[i];
		react2=&reac2[i];
		
		(*gas)-> getMoleFractions(cahe);
		(*gas2)-> setMoleFractions(cahe);
		
		(*gas2)->setDensity((*gas)->density());
		(*gas2)->setTemperature((*gas)->temperature());
		react2->setThermoMgr(*gas2[0]);
	}
}
void copycell(cell *CS1,cell *CS2)
{
	CS2->prim[0]=CS1->prim[0];
	CS2->prim[1]=CS1->prim[1];
	CS2->prim[2]=CS1->prim[2];
	CS2->prim[3]=CS1->prim[3];
	CS2->prim[4]=CS1->prim[4];
}
void copyjoin(join *JS1,join *JS2)
{
	JS2->f[0]=JS1->f[0];
	JS2->f[1]=JS1->f[1];
	JS2->f[2]=JS1->f[2];
	JS2->f[3]=JS1->f[3];
	JS2->x=JS1->x;
}
