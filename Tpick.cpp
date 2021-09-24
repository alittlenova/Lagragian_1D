#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <cmath>
using namespace std;

int main()
{
	typedef vector<vector<double> > vvs;
	int row, col, k, n,leftphimc,rightphimc,leftphipc,rightphipc;
	double rho,pressure,gamma,phim,phip,leftphim,rightphim,leftphip,rightphip,cminus,cplus;
	int axis=0;//0 for space cut and 1 for time cut
	double pickphi=2e-3;
	double pickt=3e-5;
	int ipick_phi,ipick_t;
	char inputfile[50];
	char outputfile1[50]="./globaloutput/output_Tpick.csv";
	ofstream output1(outputfile1);
	n=sprintf (inputfile, "./globaloutput/output_uprT.csv");
	ifstream infile(inputfile);
	infile>>noskipws;//不忽略空白且把每行最后那个'n'也读进来
	//判断行数和列数,文件中每列数据以空格隔开
	
	row = col = k = 0;

	char chr;
	while(infile>>chr)
	{
		//cout<<chr<<endl;
		switch(chr)
		{
		 case '\n':  //判断读入字符是否为换行符
			 row++;  //是换行符则行数＋1
			 k=1;
			break;
		 case ',':  //判断读入字符是否为空格
			if(k==0)
			{
				col++;  //则列数＋1
			}
			break;
		 default:;
		}
	 }
	infile.close();  //关闭文本文件
	
	col++;  //读文件得到的行列数均加1才是真正的行列数(注意:文件最后没有空白行)
	
	vvs data; // 定义二维字符串数组
	data.resize(row);
		
	for(k = 0; k < row; k ++)
		data[k].resize(col);
 
	ifstream infile2(inputfile);

	int trow, tcol;
	trow = tcol = 0;

	//string str;
	string str;
	const char *chr1;
	while(trow < row)//循环比较，以行数作为限制
	{
		while(tcol < col-1)
		{
			getline(infile2,str,','); //以空格为分隔对象，进行数据读入
			chr1 = str.c_str();
			data[trow][tcol] = atof(chr1);  //将空格前的数据保存到相应的字符串数组中
			tcol++;      //保存好数据后列数＋1
		}
		getline(infile2,str,'\n');  //将回车符作为对象进行数据输入，以判断一行是否结束
		chr1 = str.c_str();
		data[trow][tcol] = atof(chr1);   //将回车符前面的数据读入到相应的字符串数组中
		tcol = 0;      //列数赋值为1，准备下一行的判断
		trow ++;      //行数＋1
	}
	//cout<<data[0][0]<<" "<<data[0][col-10]<<endl;
	infile2.close();
	//find the dimension size
	int find=0;
	double max_phi,max_t;
	int max_phidim,max_tdim;
	for(int i=0;i<row;i++)
	{
		if(find==0 && data[i][0]!=data[0][0])
		{
			find=1;
			max_phi=data[i-1][1];
			max_phidim=i;
		}
	}
	max_tdim=(row/max_phidim)-1;
	max_t=data[max_phidim*(max_tdim-1)-1][0];
	
	std::cout<<"matrix dimension "<<row<<" "<<col<<std::endl;
	std::cout<<"Phi dimension "<<max_phi<<" "<<max_phidim<<std::endl;
	std::cout<<"t dimension "<<max_t<<" "<<max_tdim<<std::endl;
	if(axis==0)
	{
		find=0;
		for(int i=0;i<max_phidim;i++)
		{
			if(find==0 && data[i][1]>pickphi)
			{
				find=1;
				ipick_phi=i;
			}
		}
		std::cout<<"phi pick "<<ipick_phi<<" "<<data[ipick_phi][1]<<std::endl;
	}
	else
	{
	}
	output1.close();	
	return 0;
}
