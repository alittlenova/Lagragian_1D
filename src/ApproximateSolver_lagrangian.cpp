#include<fstream>
#include<string>
#include<math.h>
#include<iostream>
#include<stdlib.h>
#include "ApproximateSolver_lagrangian.h"
void primitivecompute(double *u1,double *u2,double *u3,double *rho,double *u,double *h,double *p)
{
	//You could call it step 1
	const double gamma=1.4;
	*rho=*u1;
	*u=(*u2)/(*u1);
	*p=(gamma-1)*((*u3)-(*u1)*(*u)*(*u)/2);
	*h=((*u3)+(*p))/(*rho);
}
void primitivecomputeHLLE(double *u1,double *u2,double *u3,double *p,double *rho,double *u,double *a,double *h,double *gamma)
{
	//You could call it step 1
	*rho=(1.0/(*u1));
	*u=*u2;
	*h=((*gamma)*(*p)/((*gamma)-1)+0.5*(*u)*(*u)*(*rho))/(*rho);
	//*a=pow((*p)*(*gamma)/(*rho),0.5);
	*a=pow(((*gamma)*(*p)*(*rho)),0.5);
}
void roeaveragecompute(double *rhol,double *rhor,double *ul,double *ur,double *hl,double *hr,double *pl,double *pr,double *gal,double *gar,double *rhorl,double *url,double *arl,double *garl)
{
	//You could call it step 2
	static double srhol,srhor;
	static double hrl;
	srhol=pow(*rhol,0.5);
	srhor=pow(*rhor,0.5);
	*rhorl=srhol*srhor;
	*url=(srhol*(*ul)+srhor*(*ur))/(srhol+srhor);
	//prl=(srhol*(*pl)+srhor*(*pr))/(srhol+srhor);
	hrl=(srhol*(*hl)+srhor*(*hr))/(srhol+srhor);
	*garl=(srhol*(*gal)+srhor*(*gar))/(srhol+srhor);
	//*arl=(srhol*(*al)+srhor*(*ar))/(srhol+srhor);
	//*arl=pow(((*garl)*(prl)*(*rhorl)),0.5);
	*arl=pow((hrl*(*rhorl)-0.5*(*url)*(*url)*(*rhorl))*((*garl)-1)*(*rhorl),0.5);
}
void wavespeedcompute(double *rhorl,double *arl,double *lambda1,double *lambda2,double *lambda3)
{
	//You could call it step 3
	*lambda1=0;
	*lambda2=(*arl);
	*lambda3=-(*arl);
}
/*void wavespeedcomputieHLL(double *lambda1,double *lambda2,double *lambda3,double vl,double vr,double ul,double ur)
{
	//You could call it step 3
	static double eta,eta0;
	eta0=1;
	eta=eta0*(std::max(0,ul-ur)/(vl+vr));
	
}*/
void wavespeedcomputieHLLE(double *lambda1,double *lambda2,double *lambda3,double *al,double *ar,double *rl,double *rr)
{
	if(*lambda3>(-*al))
	{
		*lambda3=(-*al);
	}
	if(*lambda2<(*ar))
	{
		*lambda2=(*ar);
	}
}
void wavestrengthcompute(double *rhol,double *rhor,double *pl,double *pr,double *ul,double *ur,double *arl,double *rhorl,double *v1,double *v2,double *v3)
{
	//You could call it step 4
	*v1=((*rhor)-(*rhol))-((*pr)-(*pl))/((*arl)*(*arl));
	*v2=((*ur)-(*ul))+((*pr)-(*pl))/((*rhorl)*(*arl));
	*v3=((*ur)-(*ul))-((*pr)-(*pl))/((*rhorl)*(*arl));
}
void fluxcompute(double *rho,double *u,double *p,double *f1,double *f2,double *f3)
{
	*f1=-(*u);
	*f2=*p;
	*f3=(*p)*(*u);
}
void characteristiccompute(double *r1,double *r2,double *r3,double *rhorl,double *url,double *hrl,double *arl)
{
	//You could call it step 5
	*r1=1;
	*(r1+1)=*url;
	*(r1+2)=(*url)*(*url)*0.5;
	
	*r2=(*rhorl)*0.5/(*arl);
	*(r2+1)=(*r2)*((*url)+(*arl));
	*(r2+2)=(*r2)*((*hrl)+(*arl)*(*url));
	
	*r3=-(*rhorl)*0.5/(*arl);
	*(r3+1)=(*r3)*((*url)-(*arl));
	*(r3+2)=(*r3)*((*hrl)-(*arl)*(*url));
}
void computesurfaceHLLE(double *f1l,double *f2l,double *f3l,double *f1r,double *f2r,double *f3r,double *u1l,double *u2l,double *u3l,double *u1r,double *u2r,double *u3r,double *lambda2,double *lambda3,double *f1m,double *f2m,double *f3m)
{
	//You can call it step 5*
	/*if((*lambda3)>0)
	{
		*f1m=*f1l;
		*f2m=*f2l;
		*f3m=*f3l;
	}
	else if((*lambda2)<0)
	{
		*f1m=*f1r;
		*f2m=*f2r;
		*f3m=*f3r;		
	}
	else*/
	{
		*f1m=(((*lambda2)*(*f1l)-(*lambda3)*(*f1r))/((*lambda2)-(*lambda3)))+(((*u1r)-(*u1l))*(((*lambda2)*(*lambda3))/((*lambda2)-(*lambda3))));
		*f2m=(((*lambda2)*(*f2l)-(*lambda3)*(*f2r))/((*lambda2)-(*lambda3)))+(((*u2r)-(*u2l))*(((*lambda2)*(*lambda3))/((*lambda2)-(*lambda3))));
		*f3m=(((*lambda2)*(*f3l)-(*lambda3)*(*f3r))/((*lambda2)-(*lambda3)))+(((*u3r)-(*u3l))*(((*lambda2)*(*lambda3))/((*lambda2)-(*lambda3))));
	}
	
}
void Thetacomputation(double *theta1,double *theta2,double *theta3,double *hl,double *hr,double *ul,double *ur)
{
	//Entropy Fix Parameter
	static double al,ar,n1l,n2l,n3l,n1r,n2r,n3r;
	const double gamma=1.4;
	al=pow((gamma-1)*((*hl)-(*ul)*(*ul)*0.5),0.5);
	ar=pow((gamma-1)*((*hr)-(*ur)*(*ur)*0.5),0.5);
	n1l=*ul;
	n1r=*ur;
	n2l=*ul+al;
	n2r=*ur+ar;
	n3l=*ul-al;
	n3r=*ur-ar;
	
	*theta1=std::max(0.0,(2*(n1r-n1l)));
	*theta2=std::max(0.0,(2*(n2r-n2l)));
	*theta3=std::max(0.0,(2*(n3r-n3l)));
}
void Deltacomputation(double *theta1,double *theta2,double *theta3,double *delta1,double *delta2,double *delta3,double *lambda1,double *lambda2,double *lambda3)
{
	//using for entropy fix
	static double n1,n2,n3;
	n1=abs(*lambda1);
	n2=abs(*lambda2);
	n3=abs(*lambda3);
	if(n1<*theta1)
	{
		*delta1=0.5*((n1*n1)/(*theta1)+(*theta1))-n1;
	}
	else
	{
		*delta1=0;
	}
	if(n2<*theta1)
	{
		*delta2=0.5*((n2*n2)/(*theta2)+(*theta2))-n2;
	}
	else
	{
		*delta2=0;
	}
	if(n3<*theta1)
	{
		*delta3=0.5*((n3*n3)/(*theta3)+(*theta3))-n3;
	}
	else
	{
		*delta3=0;
	}
}
void computeall(double *u1l,double *u2l,double *u3l,double *u1r,double *u2r,double *u3r,double *lambda1,double *lambda2,double *lambda3,double *r2,double *r3,double *v2,double *v3,double *u1m,double *u2m,double *u3m)
{
	//You could call it step 6, I do not use this part for any computation
	if((*lambda3)>0)
	{
		*u1m=*u1l;
		*u2m=*u2l;
		*u3m=*u3l;
	}
	else if((*lambda1)>0)
	{
		*u1m=*u1l+(*r3)*(*v3);
		*u2m=*u2l+(*(r3+1))*(*v3);
		*u3m=*u3l+(*(r3+2))*(*v3);
	}
	else if((*lambda2)>0)
	{
		*u1m=*u1r-(*r2)*(*v2);
		*u2m=*u2r-(*(r2+1))*(*v2);
		*u3m=*u3r-(*(r2+2))*(*v2);
	}
	else
	{
		*u1m=*u1r;
		*u2m=*u2r;
		*u3m=*u3r;
	}
}
void computesurface(double *rhol,double *ul,double *hl,double *rhor,double *ur,double *hr,double *lambda1,double *lambda2,double *lambda3,double *r1start,double *r2start,double *r3start,double *v1,double *v2,double *v3,double *f1m,double *f2m,double *f3m,double *pl,double *pr)
{
	*f1m=(*rhol)*(*ul)+(*rhor)*(*ur)-(*r1start)*abs(*lambda1)*(*v1)-(*r2start)*abs(*lambda2)*(*v2)-(*r3start)*abs(*lambda3)*(*v3);
	*f2m=(*rhol)*(*ul)*(*ul)+(*pl)+(*rhor)*(*ur)*(*ur)+(*pr)-(*(r1start+1))*abs(*lambda1)*(*v1)-(*(r2start+1))*abs(*lambda2)*(*v2)-(*(r3start+1))*abs(*lambda3)*(*v3);
	*f3m=(*rhol)*(*ul)*(*hl)+(*rhor)*(*ur)*(*hr)-(*(r1start+2))*abs(*lambda1)*(*v1)-(*(r2start+2))*abs(*lambda2)*(*v2)-(*(r3start+2))*abs(*lambda3)*(*v3);
	*f1m=*f1m*0.5;
	*f2m=*f2m*0.5;
	*f3m=*f3m*0.5;
}
void wavespeedbalence(double *n1,double *delta)
{
	//Entropy Fix
	static double n12,delta12;
	n12=((*n1)*(*n1));
	delta12=((*delta)*(*delta));
	if(n12<delta12)
	{
		(*n1)=(n12+delta12)/(2*(*delta));
	}
	else
	{
		(*n1)=abs(*n1);
	}
}
void wavespeedtotoalbacence(double *n1,double *n2,double *n3,double *delta1,double *delta2,double *delta3)
{
	//Entropy fix for three waves
	wavespeedbalence(n1,delta1);
	wavespeedbalence(n2,delta2);
	wavespeedbalence(n3,delta3);
}


void HLLE(double *u1l,double *u2l,double *u3l,double *pl,double *u1r,double *u2r,double *u3r,double *pr,double *f1m,double *f2m,double *f3m,double *gal,double *gar)
{
	static double rhol,rhor,ul,ur,al,ar,hr,hl; //should be computed at step 1
	static double rhorl,url,arl,garl,hrl; //should be computed at step 2
	static double lambda1,lambda2,lambda3;//should be computed at step 3
	static double f1l,f2l,f3l,f1r,f2r,f3r;//should be computed at step 4
	primitivecomputeHLLE(u1l,u2l,u3l,pl,&rhol,&ul,&al,&hl,gal);//step 1
	primitivecomputeHLLE(u1r,u2r,u3r,pr,&rhor,&ur,&ar,&hr,gar);//step 1
	roeaveragecompute(&rhol,&rhor,&ul,&ur,&hl,&hr,pl,pr,gal,gar,&rhorl,&url,&arl,&garl);//step 2
	wavespeedcompute(&rhorl,&arl,&lambda1,&lambda2,&lambda3);//step 3

	fluxcompute(&rhol,&ul,pl,&f1l,&f2l,&f3l);//step 4*left
	fluxcompute(&rhor,&ur,pr,&f1r,&f2r,&f3r);//step 4*right
	//std::cout<<f1l<<" "<<f2l<<" "<<f3l<<std::endl;
	wavespeedcomputieHLLE(&lambda1,&lambda2,&lambda3,&al,&ar,&rhol,&rhor);//Wave Speed Caliberation
	computesurfaceHLLE(&f1l,&f2l,&f3l,&f1r,&f2r,&f3r,u1l,u2l,u3l,u1r,u2r,u3r,&lambda2,&lambda3,f1m,f2m,f3m);//step 5*
	
/*	if(f1l!=f1r)
	{
		std::cout<<std::endl;
		std::cout<<f2l<<" "<<*f2m<<" "<<f2r<<std::endl;
		std::cout<<std::endl;
		std::cout<<lambda1<<" "<<lambda2<<" "<<lambda3<<std::endl;
		std::cout<<std::endl;
		std::cout<<((lambda2)*(f2l)-(lambda3)*(f2r))/((lambda2)-(lambda3))<<std::endl;
		std::cout<<std::endl;
		std::cout<<((*u2r)-(*u2l))*(((lambda2)*(lambda3))/((lambda2)-(lambda3)))<<std::endl;
		std::cout<<std::endl;
	}*/
	/*if((*u3l>*u3r)&&(-(*f1m)<*u2l)&&(-(*f1m)<*u2r))
	{
		std::cout<<"Error"<<std::endl;
		std::cout<<"P="<<*u3l<<" "<<" "<<*u3r<<std::endl;
		std::cout<<"f2="<<f3l<<" "<<(*f3m)<<" "<<f3r<<std::endl;
		std::cout<<"lambda="<<lambda2<<" "<<arl<<" "<<lambda3<<std::endl;
		std::cout<<"part_f= "<<(((lambda2)*(f3l)-(lambda3)*(f3r))/((lambda2)-(lambda3)))<<std::endl;
		std::cout<<"part_u= "<<(*u3r)-(*u3l)<<" "<<(((*u3r)-(*u3l))*(((lambda2)*(lambda3))/((lambda2)-(lambda3))))<<std::endl;
	}*/
	//std::cout<<*f1m<<" "<<*f2m<<" "<<*f3m<<std::endl;
}

