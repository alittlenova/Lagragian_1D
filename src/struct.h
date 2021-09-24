#ifndef STRUCT_H
#define STRUCT_H

class join{
public:
	double f[4];
	struct cell *lc,*rc;
	struct join *lj,*rj;
	double x;
};
class cell 
{
public:
	double prim[7];//v,u,p,e+ke,e
	double ga;
	double x;
	double dx;
	double phi;
	double dphi;
	struct cell *lc,*rc;
	struct join *lj,*rj;
	void cellcopy(cell *cell1)
	{
		for(int i=0;i<6;i++)
		{
			prim[i]=cell1->prim[i];
		}
		ga=cell1->ga;
	}
	int react;
};

#endif
