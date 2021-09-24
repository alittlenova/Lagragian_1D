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
	
	cell *CP,*CP1,*CPX;
	join *JP,*JPX;
	CP=GS;
	
	CP=GS+1;
	CP1=CS+nc-1;
	CP->cellcopy(CP1);
	CP=GS2+1;
	CP1=CS2+nc-1;
	CP->cellcopy(CP1);

	double rl,rr,ul,ur,pl,pr,f0,f1,f2,gammal,gammar,vl,vr,el,er;
	double fslope,bslope,slope;
	double f0v[5000];
	double f1v[5000];
	double f2v[5000];
	//{//second order
	#pragma omp parallel for
		for(int c=1;c<=nc;c++)
		{
			JP=JS+c;
			if(c<2||c>nc-1)
			{
				//1st order variables on boundary
				ul=JP->lc->prim[1];
				ur=JP->rc->prim[1];
				pl=JP->lc->prim[2];
				pr=JP->rc->prim[2];
				gammal=JP->lc->ga;
				gammar=JP->rc->ga;
				vl=JP->lc->prim[0];
				vr=JP->rc->prim[0];
				el=JP->lc->prim[3];
				er=JP->rc->prim[3];
				rl=1.0/vl;
				rr=1.0/vr;
			}
			else
			{
				//Reconstruct 2nd order variables
				
				//left v
				fslope=((JP->rc->prim[0])-(JP->lc->prim[0]))/(JP->lc->dphi);
				bslope=((JP->lc->prim[0])-(JP->lc->lc->prim[0]))/(JP->lc->dphi);
				slope=std::max(fslope*bslope/((fslope*fslope)+(bslope*bslope)+0.0000000001),0.0);
				slope=slope*(fslope+bslope);
				vl=JP->lc->prim[0]+(0.5*(JP->lc->dphi)*slope);
				
				//right v
				fslope=((JP->rc->rc->prim[0])-(JP->rc->prim[0]))/(JP->rc->dphi);
				bslope=((JP->rc->prim[0])-(JP->lc->prim[0]))/(JP->rc->dphi);
				slope=std::max(fslope*bslope/((fslope*fslope)+(bslope*bslope)+0.0000000001),0.0);
				slope=slope*(fslope+bslope);
				vr=JP->rc->prim[0]-(0.5*(JP->rc->dphi)*slope);
				
				/*
				//left rho
				fslope=((1.0/(JP->rc->prim[0]))-(1.0/(JP->lc->prim[0])))/(JP->lc->dphi);
				bslope=((1.0/(JP->lc->prim[0]))-(1.0/(JP->lc->lc->prim[0])))/(JP->lc->dphi);
				slope=std::max(fslope*bslope/((fslope*fslope)+(bslope*bslope)+0.0000000001),0.0);
				slope=slope*(fslope+bslope);
				rl=(1.0/(JP->lc->prim[0]))+(0.5*(JP->lc->dphi)*slope);
				
				//right rho
				fslope=((1.0/(JP->rc->rc->prim[0]))-(1.0/(JP->rc->prim[0])))/(JP->rc->dphi);
				bslope=((1.0/(JP->rc->prim[0]))-(1.0/(JP->lc->prim[0])))/(JP->rc->dphi);
				slope=std::max(fslope*bslope/((fslope*fslope)+(bslope*bslope)+0.0000000001),0.0);
				slope=slope*(fslope+bslope);
				rr=(1.0/(JP->rc->prim[0]))-(0.5*(JP->rc->dphi)*slope);
				*/
				
				//left u
				fslope=((JP->rc->prim[1])-(JP->lc->prim[1]))/(JP->lc->dphi);
				bslope=((JP->lc->prim[1])-(JP->lc->lc->prim[1]))/(JP->lc->dphi);
				slope=std::max(fslope*bslope/((fslope*fslope)+(bslope*bslope)+0.0000000001),0.0);
				slope=slope*(fslope+bslope);
				ul=JP->lc->prim[1]+(0.5*(JP->lc->dphi)*slope);
				//right u
				fslope=((JP->rc->rc->prim[1])-(JP->rc->prim[1]))/(JP->rc->dphi);
				bslope=((JP->rc->prim[1])-(JP->lc->prim[1]))/(JP->rc->dphi);
				slope=std::max(fslope*bslope/((fslope*fslope)+(bslope*bslope)+0.0000000001),0.0);
				slope=slope*(fslope+bslope);
				ur=JP->rc->prim[1]-(0.5*(JP->rc->dphi)*slope);
				
				
				//left e
				fslope=((JP->rc->prim[3])-(JP->lc->prim[3]))/(JP->lc->dphi);
				bslope=((JP->lc->prim[3])-(JP->lc->lc->prim[3]))/(JP->lc->dphi);
				slope=std::max(fslope*bslope/((fslope*fslope)+(bslope*bslope)+0.0000000001),0.0);
				slope=slope*(fslope+bslope);
				el=JP->lc->prim[3]+(0.5*(JP->lc->dphi)*slope);
				//right e
				fslope=((JP->rc->rc->prim[3])-(JP->rc->prim[3]))/(JP->rc->dphi);
				bslope=((JP->rc->prim[3])-(JP->lc->prim[3]))/(JP->rc->dphi);
				slope=std::max(fslope*bslope/((fslope*fslope)+(bslope*bslope)+0.0000000001),0.0);
				slope=slope*(fslope+bslope);
				er=JP->rc->prim[3]-(0.5*(JP->rc->dphi)*slope);
				/*
				//left p
				fslope=((JP->rc->prim[2])-(JP->lc->prim[2]))/(JP->lc->dphi);
				bslope=((JP->lc->prim[2])-(JP->lc->lc->prim[2]))/(JP->lc->dphi);
				slope=std::max(fslope*bslope/((fslope*fslope)+(bslope*bslope)+0.0000000001),0.0);
				slope=slope*(fslope+bslope);
				pl=JP->lc->prim[2]+(0.5*(JP->lc->dphi)*slope);
				//right p
				fslope=((JP->rc->rc->prim[2])-(JP->rc->prim[2]))/(JP->rc->dphi);
				bslope=((JP->rc->prim[2])-(JP->lc->prim[2]))/(JP->rc->dphi);
				slope=std::max(fslope*bslope/((fslope*fslope)+(bslope*bslope)+0.0000000001),0.0);
				slope=slope*(fslope+bslope);
				pr=JP->rc->prim[2]-(0.5*(JP->rc->dphi)*slope);
				*/
				//gamma
				gammal=JP->lc->ga;
				gammar=JP->rc->ga;
				
			}
			//flux solver
			//vl=1.0/rl;
			//vr=1.0/rr;
			rl=1.0/vl;
			rr=1.0/vr;
			pl=(el-(0.5*ul*ul))*rl*(gammal-1);
			pr=(er-(0.5*ur*ur))*rr*(gammar-1);
			//el=(pl/(rl*(gammal-1)))+(0.5*ul*ul);
			//er=(pr/(rr*(gammal-1)))+(0.5*ur*ur);
			
			//fluidcons(&rl,&rr,&ul,&ur,&pl,&pr,&f0,&f1,&f2,gammal,gammar);
			HLLE(&vl,&ul,&el,&pl,&vr,&ur,&er,&pr,&f0,&f1,&f2,&gammal,&gammar);
			
			//save 1st time result for 2nd order time march
			f0v[c-1]=f0;
			f1v[c-1]=f1;
			f2v[c-1]=f2;
			JP->f[0]=f0;
			JP->f[1]=f1;
			JP->f[2]=f2;
			
		}
		//std::cout<<"A"<<std::endl;
		for(int c=0;c<=nc;c++)
		{
			JP=JS+c;
			JPX=JS2+c;
			copyjoin(JP,JPX);
		}
		maxt=0.69*maxdt(CS,nc);
		//std::cout<<"B"<<std::endl;
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
		#pragma omp parallel for
		for(int c=0;c<nc;c++)
		{
			CP=CS+c;
			CPX=CS2+c;
			copycell(CP,CPX);
		}
		//std::cout<<"C"<<std::endl;
		//First time march CS2 update
		#pragma omp parallel for
		for(int c=0;c<nc;c++)
		{
			CP=CS2+c;
			
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

			}
			else
			{
				CP->prim[0]=CP->prim[0]+((CP->lj->f[0])-(CP->rj->f[0]))*maxt/(CP->dphi);
				CP->prim[1]=CP->prim[1]+((CP->lj->f[1])-(CP->rj->f[1]))*maxt/(CP->dphi);
				CP->prim[3]=CP->prim[3]+((CP->lj->f[2])-(CP->rj->f[2]))*maxt/(CP->dphi);
				CP->prim[4]=CP->prim[3]-(0.5*CP->prim[1]*CP->prim[1]);
				CP->prim[2]=CP->prim[4]*((CP->ga)-1.0)/(CP->prim[0]);
				CP->dx=(CP->dphi)*(CP->prim[0]);
			}
		}
		//std::cout<<"D"<<std::endl;
		//First time Chemical
		//copy gas1 to gas2
		copyreactor(gass,reac,gass2,reac2,time,nc);
		//set up the first hydro march
		exchange_to_chemical_part(CS2,nc,gass2,reac2);
		for(int i=0;i<nc;i++)
		{
			sim2[i].setInitialTime(time);
		}
		#pragma omp parallel for
		for(int i=0;i<ncell;i++)
		{
			//sim2[i].advance(time+maxt);
		}
		exchange_to_hydro_part(CS2,nc,gass2,reac2);
		//}
		
		//std::cout<<"E"<<std::endl;
		//second order Next time
		//{
		#pragma omp parallel for
		for(int c=1;c<=nc;c++)
		{
			JP=JS2+c;
			if(c<2||c>nc-1)
			{
				//1st order variables on boundary
				ul=JP->lc->prim[1];
				ur=JP->rc->prim[1];
				pl=JP->lc->prim[2];
				pr=JP->rc->prim[2];
				gammal=JP->lc->ga;
				gammar=JP->rc->ga;
				vl=JP->lc->prim[0];
				vr=JP->rc->prim[0];
				el=JP->lc->prim[3];
				er=JP->rc->prim[3];
				rl=1.0/vl;
				rr=1.0/vr;
			}
			else
			{
				//Reconstruct 2nd order variables
				
				//left v
				fslope=((JP->rc->prim[0])-(JP->lc->prim[0]))/(JP->lc->dphi);
				bslope=((JP->lc->prim[0])-(JP->lc->lc->prim[0]))/(JP->lc->dphi);
				slope=std::max(fslope*bslope/((fslope*fslope)+(bslope*bslope)+0.0000000001),0.0);
				slope=slope*(fslope+bslope);
				vl=JP->lc->prim[0]+(0.5*(JP->lc->dphi)*slope);
				//right v
				fslope=((JP->rc->rc->prim[0])-(JP->rc->prim[0]))/(JP->rc->dphi);
				bslope=((JP->rc->prim[0])-(JP->lc->prim[0]))/(JP->rc->dphi);
				slope=std::max(fslope*bslope/((fslope*fslope)+(bslope*bslope)+0.0000000001),0.0);
				slope=slope*(fslope+bslope);
				vr=JP->rc->prim[0]-(0.5*(JP->rc->dphi)*slope);
				
				/*
				//left rho
				fslope=((1.0/(JP->rc->prim[0]))-(1.0/(JP->lc->prim[0])))/(JP->lc->dphi);
				bslope=((1.0/(JP->lc->prim[0]))-(1.0/(JP->lc->lc->prim[0])))/(JP->lc->dphi);
				slope=std::max(fslope*bslope/((fslope*fslope)+(bslope*bslope)+0.0000000001),0.0);
				slope=slope*(fslope+bslope);
				rl=(1.0/(JP->lc->prim[0]))+(0.5*(JP->lc->dphi)*slope);
				
				//right rho
				fslope=((1.0/(JP->rc->rc->prim[0]))-(1.0/(JP->rc->prim[0])))/(JP->rc->dphi);
				bslope=((1.0/(JP->rc->prim[0]))-(1.0/(JP->lc->prim[0])))/(JP->rc->dphi);
				slope=std::max(fslope*bslope/((fslope*fslope)+(bslope*bslope)+0.0000000001),0.0);
				slope=slope*(fslope+bslope);
				rr=(1.0/(JP->rc->prim[0]))-(0.5*(JP->rc->dphi)*slope);
				*/
				//left u
				fslope=((JP->rc->prim[1])-(JP->lc->prim[1]))/(JP->lc->dphi);
				bslope=((JP->lc->prim[1])-(JP->lc->lc->prim[1]))/(JP->lc->dphi);
				slope=std::max(fslope*bslope/((fslope*fslope)+(bslope*bslope)+0.0000000001),0.0);
				slope=slope*(fslope+bslope);
				ul=JP->lc->prim[1]+(0.5*(JP->lc->dphi)*slope);
				//right u
				fslope=((JP->rc->rc->prim[1])-(JP->rc->prim[1]))/(JP->rc->dphi);
				bslope=((JP->rc->prim[1])-(JP->lc->prim[1]))/(JP->rc->dphi);
				slope=std::max(fslope*bslope/((fslope*fslope)+(bslope*bslope)+0.0000000001),0.0);
				slope=slope*(fslope+bslope);
				ur=JP->rc->prim[1]-(0.5*(JP->rc->dphi)*slope);
				
				
				//left e
				fslope=((JP->rc->prim[3])-(JP->lc->prim[3]))/(JP->lc->dphi);
				bslope=((JP->lc->prim[3])-(JP->lc->lc->prim[3]))/(JP->lc->dphi);
				slope=std::max(fslope*bslope/((fslope*fslope)+(bslope*bslope)+0.0000000001),0.0);
				slope=slope*(fslope+bslope);
				el=JP->lc->prim[3]+(0.5*(JP->lc->dphi)*slope);
				//right e
				fslope=((JP->rc->rc->prim[3])-(JP->rc->prim[3]))/(JP->rc->dphi);
				bslope=((JP->rc->prim[3])-(JP->lc->prim[3]))/(JP->rc->dphi);
				slope=std::max(fslope*bslope/((fslope*fslope)+(bslope*bslope)+0.0000000001),0.0);
				slope=slope*(fslope+bslope);
				er=JP->rc->prim[3]-(0.5*(JP->rc->dphi)*slope);
				
				/*
				//left p
				fslope=((JP->rc->prim[2])-(JP->lc->prim[2]))/(JP->lc->dphi);
				bslope=((JP->lc->prim[2])-(JP->lc->lc->prim[2]))/(JP->lc->dphi);
				slope=std::max(fslope*bslope/((fslope*fslope)+(bslope*bslope)+0.0000000001),0.0);
				slope=slope*(fslope+bslope);
				pl=JP->lc->prim[2]+(0.5*(JP->lc->dphi)*slope);
 
				//right p
				fslope=((JP->rc->rc->prim[2])-(JP->rc->prim[2]))/(JP->rc->dphi);
				bslope=((JP->rc->prim[2])-(JP->lc->prim[2]))/(JP->rc->dphi);
				slope=std::max(fslope*bslope/((fslope*fslope)+(bslope*bslope)+0.0000000001),0.0);
				slope=slope*(fslope+bslope);
				pr=JP->rc->prim[2]-(0.5*(JP->rc->dphi)*slope);
				*/
				//gamma
				gammal=JP->lc->ga;
				gammar=JP->rc->ga;
			}
			//flux solver
			
			//vl=1.0/rl;
			//vr=1.0/rr;
			
			rl=1.0/vl;
			rr=1.0/vr;
			pl=(el-(0.5*ul*ul))*rl*(gammal-1);
			pr=(er-(0.5*ur*ur))*rr*(gammar-1);
			//el=(pl/(rl*(gammal-1)))+(0.5*ul*ul);
			//er=(pr/(rr*(gammal-1)))+(0.5*ur*ur);
			//fluidcons(&rl,&rr,&ul,&ur,&pl,&pr,&f0,&f1,&f2,gammal,gammar);
			HLLE(&vl,&ul,&el,&pl,&vr,&ur,&er,&pr,&f0,&f1,&f2,&gammal,&gammar);

			JP->f[0]=(0.5*f0)+(0.5*f0v[c-1]);
			JP->f[1]=(0.5*f1)+(0.5*f1v[c-1]);
			JP->f[2]=(0.5*f2)+(0.5*f2v[c-1]);
			
		}
		//std::cout<<"F"<<std::endl;
		#pragma omp parallel for
		//JS2->x=(JS2->x)+up*maxt;
		for(int c=1;c<=nc;c++)
		{
			JP=JS2+c;
			JP->x=(JP->x)-(JP->f[0])*maxt;
		}
		
		
		#pragma omp parallel for
		for(int c=0;c<=nc;c++)
		{
			JP=JS+c;
			JPX=JS2+c;
			copyjoin(JPX,JP);
		}
		//std::cout<<"G"<<std::endl;
/*		for(int c=0;c<nc;c++)
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
				CP->prim[0]=CP->prim[0]-((CP->lj->f[0])-f0v[c])*maxt/(CP->dphi);
				CP->prim[1]=CP->prim[1]-((CP->lj->f[1])-f1v[c])*maxt/(CP->dphi);
				CP->prim[3]=CP->prim[3]-(-(f2v[c])-((CP->lj->f[0])*(CP->lj->f[1])))*maxt/(CP->dphi);
				CP->prim[4]=CP->prim[3]-(0.5*CP->prim[1]*CP->prim[1]);
				CP->prim[2]=CP->prim[4]*((CP->ga)-1.0)/(CP->prim[0]);
				CP->dx=(CP->dphi)*(CP->prim[0]);
			}
			else
			{
				CP->prim[0]=CP->prim[0]-(f0v[c-1]-f0v[c])*maxt/(CP->dphi);
				CP->prim[1]=CP->prim[1]-(f1v[c-1]-f1v[c])*maxt/(CP->dphi);
				CP->prim[3]=CP->prim[3]-(f2v[c-1]-f2v[c])*maxt/(CP->dphi);
				CP->prim[4]=CP->prim[3]-(0.5*CP->prim[1]*CP->prim[1]);
				CP->prim[2]=CP->prim[4]*((CP->ga)-1.0)/(CP->prim[0]);
				CP->dx=(CP->dphi)*(CP->prim[0]);
				CP->x=(CP->lj->x)+0.5*(CP->dx);
			}
			
		}*/
		#pragma omp parallel for
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
		//std::cout<<"H"<<std::endl;
		
	//}
	return maxt;
}
double marchtest(cell *CS,join *JS,cell *GS,cell *CS2,join *JS2,cell *GS2,cell *left,int nc)
{
	double maxt;
	double up=funcp(0);
	
	cell *CP,*CP1,*CPX;
	join *JP,*JPX;
	CP=GS;
	
	CP=GS+1;
	CP1=CS+nc-1;
	CP->cellcopy(CP1);
	CP=GS2+1;
	CP1=CS2+nc-1;
	CP->cellcopy(CP1);

	double rl,rr,ul,ur,pl,pr,f0,f1,f2,gammal,gammar,vl,vr,el,er;
	double fslope,bslope,slope;
	double f0v[5000];
	double f1v[5000];
	double f2v[5000];
	//{//second order
	#pragma omp parallel for
		for(int c=1;c<=nc;c++)
		{
			JP=JS+c;
			if(c<2||c>nc-1)
			{
				//1st order variables on boundary
				ul=JP->lc->prim[1];
				ur=JP->rc->prim[1];
				pl=JP->lc->prim[2];
				pr=JP->rc->prim[2];
				gammal=JP->lc->ga;
				gammar=JP->rc->ga;
				vl=JP->lc->prim[0];
				vr=JP->rc->prim[0];
				el=JP->lc->prim[3];
				er=JP->rc->prim[3];
				rl=1.0/vl;
				rr=1.0/vr;
			}
			else
			{
				//Reconstruct 2nd order variables
				
				//left v
				fslope=((JP->rc->prim[0])-(JP->lc->prim[0]))/(JP->lc->dphi);
				bslope=((JP->lc->prim[0])-(JP->lc->lc->prim[0]))/(JP->lc->dphi);
				slope=std::max(fslope*bslope/((fslope*fslope)+(bslope*bslope)+0.0000000001),0.0);
				slope=slope*(fslope+bslope);
				vl=JP->lc->prim[0]+(0.5*(JP->lc->dphi)*slope);
				
				//right v
				fslope=((JP->rc->rc->prim[0])-(JP->rc->prim[0]))/(JP->rc->dphi);
				bslope=((JP->rc->prim[0])-(JP->lc->prim[0]))/(JP->rc->dphi);
				slope=std::max(fslope*bslope/((fslope*fslope)+(bslope*bslope)+0.0000000001),0.0);
				slope=slope*(fslope+bslope);
				vr=JP->rc->prim[0]-(0.5*(JP->rc->dphi)*slope);
				
				/*
				//left rho
				fslope=((1.0/(JP->rc->prim[0]))-(1.0/(JP->lc->prim[0])))/(JP->lc->dphi);
				bslope=((1.0/(JP->lc->prim[0]))-(1.0/(JP->lc->lc->prim[0])))/(JP->lc->dphi);
				slope=std::max(fslope*bslope/((fslope*fslope)+(bslope*bslope)+0.0000000001),0.0);
				slope=slope*(fslope+bslope);
				rl=(1.0/(JP->lc->prim[0]))+(0.5*(JP->lc->dphi)*slope);
				
				//right rho
				fslope=((1.0/(JP->rc->rc->prim[0]))-(1.0/(JP->rc->prim[0])))/(JP->rc->dphi);
				bslope=((1.0/(JP->rc->prim[0]))-(1.0/(JP->lc->prim[0])))/(JP->rc->dphi);
				slope=std::max(fslope*bslope/((fslope*fslope)+(bslope*bslope)+0.0000000001),0.0);
				slope=slope*(fslope+bslope);
				rr=(1.0/(JP->rc->prim[0]))-(0.5*(JP->rc->dphi)*slope);
				*/
				
				//left u
				fslope=((JP->rc->prim[1])-(JP->lc->prim[1]))/(JP->lc->dphi);
				bslope=((JP->lc->prim[1])-(JP->lc->lc->prim[1]))/(JP->lc->dphi);
				slope=std::max(fslope*bslope/((fslope*fslope)+(bslope*bslope)+0.0000000001),0.0);
				slope=slope*(fslope+bslope);
				ul=JP->lc->prim[1]+(0.5*(JP->lc->dphi)*slope);
				//right u
				fslope=((JP->rc->rc->prim[1])-(JP->rc->prim[1]))/(JP->rc->dphi);
				bslope=((JP->rc->prim[1])-(JP->lc->prim[1]))/(JP->rc->dphi);
				slope=std::max(fslope*bslope/((fslope*fslope)+(bslope*bslope)+0.0000000001),0.0);
				slope=slope*(fslope+bslope);
				ur=JP->rc->prim[1]-(0.5*(JP->rc->dphi)*slope);
				
				
				//left e
				fslope=((JP->rc->prim[3])-(JP->lc->prim[3]))/(JP->lc->dphi);
				bslope=((JP->lc->prim[3])-(JP->lc->lc->prim[3]))/(JP->lc->dphi);
				slope=std::max(fslope*bslope/((fslope*fslope)+(bslope*bslope)+0.0000000001),0.0);
				slope=slope*(fslope+bslope);
				el=JP->lc->prim[3]+(0.5*(JP->lc->dphi)*slope);
				//right e
				fslope=((JP->rc->rc->prim[3])-(JP->rc->prim[3]))/(JP->rc->dphi);
				bslope=((JP->rc->prim[3])-(JP->lc->prim[3]))/(JP->rc->dphi);
				slope=std::max(fslope*bslope/((fslope*fslope)+(bslope*bslope)+0.0000000001),0.0);
				slope=slope*(fslope+bslope);
				er=JP->rc->prim[3]-(0.5*(JP->rc->dphi)*slope);
				/*
				//left p
				fslope=((JP->rc->prim[2])-(JP->lc->prim[2]))/(JP->lc->dphi);
				bslope=((JP->lc->prim[2])-(JP->lc->lc->prim[2]))/(JP->lc->dphi);
				slope=std::max(fslope*bslope/((fslope*fslope)+(bslope*bslope)+0.0000000001),0.0);
				slope=slope*(fslope+bslope);
				pl=JP->lc->prim[2]+(0.5*(JP->lc->dphi)*slope);
				//right p
				fslope=((JP->rc->rc->prim[2])-(JP->rc->prim[2]))/(JP->rc->dphi);
				bslope=((JP->rc->prim[2])-(JP->lc->prim[2]))/(JP->rc->dphi);
				slope=std::max(fslope*bslope/((fslope*fslope)+(bslope*bslope)+0.0000000001),0.0);
				slope=slope*(fslope+bslope);
				pr=JP->rc->prim[2]-(0.5*(JP->rc->dphi)*slope);
				*/
				//gamma
				gammal=JP->lc->ga;
				gammar=JP->rc->ga;
				
			}
			//flux solver
			//vl=1.0/rl;
			//vr=1.0/rr;
			rl=1.0/vl;
			rr=1.0/vr;
			pl=(el-(0.5*ul*ul))*rl*(gammal-1);
			pr=(er-(0.5*ur*ur))*rr*(gammar-1);
			//el=(pl/(rl*(gammal-1)))+(0.5*ul*ul);
			//er=(pr/(rr*(gammal-1)))+(0.5*ur*ur);
			
			//fluidcons(&rl,&rr,&ul,&ur,&pl,&pr,&f0,&f1,&f2,gammal,gammar);
			HLLE(&vl,&ul,&el,&pl,&vr,&ur,&er,&pr,&f0,&f1,&f2,&gammal,&gammar);
			
			//save 1st time result for 2nd order time march
			f0v[c-1]=f0;
			f1v[c-1]=f1;
			f2v[c-1]=f2;
			JP->f[0]=f0;
			JP->f[1]=f1;
			JP->f[2]=f2;
			
		}
		//std::cout<<"A"<<std::endl;
		for(int c=0;c<=nc;c++)
		{
			JP=JS+c;
			JPX=JS2+c;
			copyjoin(JP,JPX);
		}
		maxt=0.69*maxdt(CS,nc);
		//std::cout<<"B"<<std::endl;
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
		
		#pragma omp parallel for
		for(int c=0;c<nc;c++)
		{
			CP=CS+c;
			CPX=CS2+c;
			copycell(CP,CPX);
		}
		//std::cout<<"C"<<std::endl;
		//First time march CS2 update
		#pragma omp parallel for
		for(int c=0;c<nc;c++)
		{
			CP=CS2+c;
			
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

			}
			else
			{
				CP->prim[0]=CP->prim[0]+((CP->lj->f[0])-(CP->rj->f[0]))*maxt/(CP->dphi);
				CP->prim[1]=CP->prim[1]+((CP->lj->f[1])-(CP->rj->f[1]))*maxt/(CP->dphi);
				CP->prim[3]=CP->prim[3]+((CP->lj->f[2])-(CP->rj->f[2]))*maxt/(CP->dphi);
				CP->prim[4]=CP->prim[3]-(0.5*CP->prim[1]*CP->prim[1]);
				CP->prim[2]=CP->prim[4]*((CP->ga)-1.0)/(CP->prim[0]);
				CP->dx=(CP->dphi)*(CP->prim[0]);
			}
		}
		//std::cout<<"D"<<std::endl;
		//First time Chemical
		//copy gas1 to gas2
		//copyreactor(gass,reac,gass2,reac2,time,nc);
		//set up the first hydro march
		//exchange_to_chemical_part(CS2,nc,gass2,reac2);
		for(int i=0;i<nc;i++)
		{
			//sim2[i].setInitialTime(time);
		}
		#pragma omp parallel for
		for(int i=0;i<ncell;i++)
		{
			//sim2[i].advance(time+maxt);
		}
		//exchange_to_hydro_part(CS2,nc,gass2,reac2);
		//}
		
		//std::cout<<"E"<<std::endl;
		//second order Next time
		//{
		#pragma omp parallel for
		for(int c=1;c<=nc;c++)
		{
			JP=JS2+c;
			if(c<2||c>nc-1)
			{
				//1st order variables on boundary
				ul=JP->lc->prim[1];
				ur=JP->rc->prim[1];
				pl=JP->lc->prim[2];
				pr=JP->rc->prim[2];
				gammal=JP->lc->ga;
				gammar=JP->rc->ga;
				vl=JP->lc->prim[0];
				vr=JP->rc->prim[0];
				el=JP->lc->prim[3];
				er=JP->rc->prim[3];
				rl=1.0/vl;
				rr=1.0/vr;
			}
			else
			{
				//Reconstruct 2nd order variables
				
				//left v
				fslope=((JP->rc->prim[0])-(JP->lc->prim[0]))/(JP->lc->dphi);
				bslope=((JP->lc->prim[0])-(JP->lc->lc->prim[0]))/(JP->lc->dphi);
				slope=std::max(fslope*bslope/((fslope*fslope)+(bslope*bslope)+0.0000000001),0.0);
				slope=slope*(fslope+bslope);
				vl=JP->lc->prim[0]+(0.5*(JP->lc->dphi)*slope);
				//right v
				fslope=((JP->rc->rc->prim[0])-(JP->rc->prim[0]))/(JP->rc->dphi);
				bslope=((JP->rc->prim[0])-(JP->lc->prim[0]))/(JP->rc->dphi);
				slope=std::max(fslope*bslope/((fslope*fslope)+(bslope*bslope)+0.0000000001),0.0);
				slope=slope*(fslope+bslope);
				vr=JP->rc->prim[0]-(0.5*(JP->rc->dphi)*slope);
				
				/*
				//left rho
				fslope=((1.0/(JP->rc->prim[0]))-(1.0/(JP->lc->prim[0])))/(JP->lc->dphi);
				bslope=((1.0/(JP->lc->prim[0]))-(1.0/(JP->lc->lc->prim[0])))/(JP->lc->dphi);
				slope=std::max(fslope*bslope/((fslope*fslope)+(bslope*bslope)+0.0000000001),0.0);
				slope=slope*(fslope+bslope);
				rl=(1.0/(JP->lc->prim[0]))+(0.5*(JP->lc->dphi)*slope);
				
				//right rho
				fslope=((1.0/(JP->rc->rc->prim[0]))-(1.0/(JP->rc->prim[0])))/(JP->rc->dphi);
				bslope=((1.0/(JP->rc->prim[0]))-(1.0/(JP->lc->prim[0])))/(JP->rc->dphi);
				slope=std::max(fslope*bslope/((fslope*fslope)+(bslope*bslope)+0.0000000001),0.0);
				slope=slope*(fslope+bslope);
				rr=(1.0/(JP->rc->prim[0]))-(0.5*(JP->rc->dphi)*slope);
				*/
				//left u
				fslope=((JP->rc->prim[1])-(JP->lc->prim[1]))/(JP->lc->dphi);
				bslope=((JP->lc->prim[1])-(JP->lc->lc->prim[1]))/(JP->lc->dphi);
				slope=std::max(fslope*bslope/((fslope*fslope)+(bslope*bslope)+0.0000000001),0.0);
				slope=slope*(fslope+bslope);
				ul=JP->lc->prim[1]+(0.5*(JP->lc->dphi)*slope);
				//right u
				fslope=((JP->rc->rc->prim[1])-(JP->rc->prim[1]))/(JP->rc->dphi);
				bslope=((JP->rc->prim[1])-(JP->lc->prim[1]))/(JP->rc->dphi);
				slope=std::max(fslope*bslope/((fslope*fslope)+(bslope*bslope)+0.0000000001),0.0);
				slope=slope*(fslope+bslope);
				ur=JP->rc->prim[1]-(0.5*(JP->rc->dphi)*slope);
				
				
				//left e
				fslope=((JP->rc->prim[3])-(JP->lc->prim[3]))/(JP->lc->dphi);
				bslope=((JP->lc->prim[3])-(JP->lc->lc->prim[3]))/(JP->lc->dphi);
				slope=std::max(fslope*bslope/((fslope*fslope)+(bslope*bslope)+0.0000000001),0.0);
				slope=slope*(fslope+bslope);
				el=JP->lc->prim[3]+(0.5*(JP->lc->dphi)*slope);
				//right e
				fslope=((JP->rc->rc->prim[3])-(JP->rc->prim[3]))/(JP->rc->dphi);
				bslope=((JP->rc->prim[3])-(JP->lc->prim[3]))/(JP->rc->dphi);
				slope=std::max(fslope*bslope/((fslope*fslope)+(bslope*bslope)+0.0000000001),0.0);
				slope=slope*(fslope+bslope);
				er=JP->rc->prim[3]-(0.5*(JP->rc->dphi)*slope);
				
				/*
				//left p
				fslope=((JP->rc->prim[2])-(JP->lc->prim[2]))/(JP->lc->dphi);
				bslope=((JP->lc->prim[2])-(JP->lc->lc->prim[2]))/(JP->lc->dphi);
				slope=std::max(fslope*bslope/((fslope*fslope)+(bslope*bslope)+0.0000000001),0.0);
				slope=slope*(fslope+bslope);
				pl=JP->lc->prim[2]+(0.5*(JP->lc->dphi)*slope);
 
				//right p
				fslope=((JP->rc->rc->prim[2])-(JP->rc->prim[2]))/(JP->rc->dphi);
				bslope=((JP->rc->prim[2])-(JP->lc->prim[2]))/(JP->rc->dphi);
				slope=std::max(fslope*bslope/((fslope*fslope)+(bslope*bslope)+0.0000000001),0.0);
				slope=slope*(fslope+bslope);
				pr=JP->rc->prim[2]-(0.5*(JP->rc->dphi)*slope);
				*/
				//gamma
				gammal=JP->lc->ga;
				gammar=JP->rc->ga;
			}
			//flux solver
			
			//vl=1.0/rl;
			//vr=1.0/rr;
			
			rl=1.0/vl;
			rr=1.0/vr;
			pl=(el-(0.5*ul*ul))*rl*(gammal-1);
			pr=(er-(0.5*ur*ur))*rr*(gammar-1);
			//el=(pl/(rl*(gammal-1)))+(0.5*ul*ul);
			//er=(pr/(rr*(gammal-1)))+(0.5*ur*ur);
			//fluidcons(&rl,&rr,&ul,&ur,&pl,&pr,&f0,&f1,&f2,gammal,gammar);
			HLLE(&vl,&ul,&el,&pl,&vr,&ur,&er,&pr,&f0,&f1,&f2,&gammal,&gammar);

			JP->f[0]=(0.5*f0)+(0.5*f0v[c-1]);
			JP->f[1]=(0.5*f1)+(0.5*f1v[c-1]);
			JP->f[2]=(0.5*f2)+(0.5*f2v[c-1]);
			
		}
		//std::cout<<"F"<<std::endl;
		#pragma omp parallel for
		for(int c=1;c<=nc;c++)
		{
			JP=JS2+c;
			JP->x=(JP->x)-(JP->f[0])*maxt;
		}
		
		
		#pragma omp parallel for
		for(int c=0;c<=nc;c++)
		{
			JP=JS+c;
			JPX=JS2+c;
			copyjoin(JPX,JP);
		}
		//std::cout<<"G"<<std::endl;
/*		for(int c=0;c<nc;c++)
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
				CP->prim[0]=CP->prim[0]-((CP->lj->f[0])-f0v[c])*maxt/(CP->dphi);
				CP->prim[1]=CP->prim[1]-((CP->lj->f[1])-f1v[c])*maxt/(CP->dphi);
				CP->prim[3]=CP->prim[3]-(-(f2v[c])-((CP->lj->f[0])*(CP->lj->f[1])))*maxt/(CP->dphi);
				CP->prim[4]=CP->prim[3]-(0.5*CP->prim[1]*CP->prim[1]);
				CP->prim[2]=CP->prim[4]*((CP->ga)-1.0)/(CP->prim[0]);
				CP->dx=(CP->dphi)*(CP->prim[0]);
			}
			else
			{
				CP->prim[0]=CP->prim[0]-(f0v[c-1]-f0v[c])*maxt/(CP->dphi);
				CP->prim[1]=CP->prim[1]-(f1v[c-1]-f1v[c])*maxt/(CP->dphi);
				CP->prim[3]=CP->prim[3]-(f2v[c-1]-f2v[c])*maxt/(CP->dphi);
				CP->prim[4]=CP->prim[3]-(0.5*CP->prim[1]*CP->prim[1]);
				CP->prim[2]=CP->prim[4]*((CP->ga)-1.0)/(CP->prim[0]);
				CP->dx=(CP->dphi)*(CP->prim[0]);
				CP->x=(CP->lj->x)+0.5*(CP->dx);
			}
			
		}*/
		#pragma omp parallel for
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
		//std::cout<<"H"<<std::endl;
		
	//}
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
		C->prim[4]=C->prim[3]-(0.5*C->prim[1]*C->prim[1]);
		C->prim[5]=C->prim[3]-(0.5*C->prim[1]*C->prim[1]);
		C->prim[6]=0;
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
		output1<<CI->x<<","<<(1/(CI->prim[0]))<<","<<CI->prim[1]<<","<<CI->prim[2]<<","<<CI->prim[3]<<","<<(*GI)->temperature()<<","<<CI->ga<<","<<CI->phi<<","<<pow((1/(CI->prim[0]))*CI->ga*CI->prim[2],0.5);
		//1 xlocation, 2 density, 3 velocity, 4 pressure, 5 internalE per mass, 6 Temperature, 7 Gamma, 8 Phi,9 c
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
	CS2->x=CS1->x;
}
void copyjoin(join *JS1,join *JS2)
{
	JS2->f[0]=JS1->f[0];
	JS2->f[1]=JS1->f[1];
	JS2->f[2]=JS1->f[2];
	JS2->f[3]=JS1->f[3];
	JS2->x=JS1->x;
}

