#include "dataread.h"

double** dataread(const char filename[], int* row, int* col)
{
    std::ifstream file;
    file.open(filename);

    // You forgot to read the Col Header out.
    // Just read them into strings.
    int m,n;
    if(file.good())
    {
		double data1;
		int i=0;
		int j=0;
		int find_column=0;
	
		while (!file.eof())
		{	  
			if(find_column==0&&file.peek()==10)
			{
				n=j;
				find_column=1;
			}
			i++;
			file >> data1;
			if(find_column==0)
			{
				j++;
			}
		}
		m=(i-1)/j;
	}
	file.close();
    // Don't do manual memory management.
    // A vector of vectors gives you a 2D array.
    // The initialize(r) is slightly more complex but given the example above
    // I think you should be able to see the outer vector is initialized with
    // n copies of a vector with m elements. 
    //std::vector<std::vector<double> >   data(n,std::vector<double>(m));
	file.open(filename);

	double **data=new double*[m];
	for(int i=0;i<m;i++)
	{
		double *dat=new double[n];
		data[i]=dat;
	}

	for(int i=0;i<m;i++)
	{
		for(int j=0;j<n;j++)
		{
			file >> data[i][j];
		}
	}
	*row=m;
	*col=n;
    file.close();
    return data;
}
/*int main()
{
	int a,b;
	char filename[]="./C2H4_sp.plt";
	double **matrix=dataread(filename, &a, &b);
	std::cout<<matrix[1][0]<<std::endl;
	return 0;
}*/
