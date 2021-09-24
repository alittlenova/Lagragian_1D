#ifndef APPROXIMATESOLVERLAG_H
#define APPROXIMATESOLVERLAG_H
void primitivecompute(double *u1,double *u2,double *u3,double *rho,double *u,double *h,double *p);
void primitivecomputeHLLE(double *u1,double *u2,double *u3,double *p,double *rho,double *u,double *a,double *h,double *gamma);
void roeaveragecompute(double *rhol,double *rhor,double *ul,double *ur,double *hl,double *hr,double *pl,double *pr,double *gar,double *gal,double *rhorl,double *url,double *arl,double *garl);
void wavespeedcompute(double *rhorl,double *arl,double *lambda1,double *lambda2,double *lambda3);
void wavestrengthcompute(double *rhol,double *rhor,double *pl,double *pr,double *ul,double *ur,double *arl,double *rhorl,double *v1,double *v2,double *v3);
void fluxcompute(double *rho,double *u,double *p,double *f1,double *f2,double *f3);
void characteristiccompute(double *r1,double *r2,double *r3,double *rhorl,double *url,double *hrl,double *arl);
void computesurfaceHLLE(double *f1l,double *f2l,double *f3l,double *f1r,double *f2r,double *f3r,double *u1l,double *u2l,double *u3l,double *u1r,double *u2r,double *u3r,double *lambda2,double *lambda3,double *f1m,double *f2m,double *f3m);
void Thetacomputation(double *theta1,double *theta2,double *theta3,double *hl,double *hr,double *ul,double *ur);
void Deltacomputation(double *theta1,double *theta2,double *theta3,double *delta1,double *delta2,double *delta3,double *lambda1,double *lambda2,double *lambda3);
void computeall(double *u1l,double *u2l,double *u3l,double *u1r,double *u2r,double *u3r,double *lambda1,double *lambda2,double *lambda3,double *r2,double *r3,double *v2,double *v3,double *u1m,double *u2m,double *u3m);
void computesurface(double *rhol,double *ul,double *hl,double *rhor,double *ur,double *hr,double *lambda1,double *lambda2,double *lambda3,double *r1start,double *r2start,double *r3start,double *v1,double *v2,double *v3,double *f1m,double *f2m,double *f3m,double *pl,double *pr);
void wavespeedbalence(double *n1,double *delta);
void wavespeedtotoalbacence(double *n1,double *n2,double *n3,double *delta1,double *delta2,double *delta3);
void wavespeedcomputieHLLE(double *lambda1,double *lambda2,double *lambda3,double *al,double *ar,double rl,double rr);
void RoeNonFix(double *u1l,double *u2l,double *u3l,double *u1r,double *u2r,double *u3r,double *f1m,double *f2m,double *f3m);
void RoeFix(double *u1l,double *u2l,double *u3l,double *u1r,double *u2r,double *u3r,double *f1m,double *f2m,double *f3m);
void HLLE(double *u1l,double *u2l,double *u3l,double *pl,double *u1r,double *u2r,double *u3r,double *pr,double *f1m,double *f2m,double *f3m,double *gal,double *gar);
#endif
