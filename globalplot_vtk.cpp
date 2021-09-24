#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <cmath>
#include <cstdlib>
using namespace std;

int main()
{
	typedef vector<vector<double> > vvs;
	double tm;
//	const double up=1500;
	const int celln=5000;
	int picnum=2000;
	const double dphi=1.5e-2/celln;
	const double dt=1.0e-4/picnum;
	const int characterdensity=10;
	char inputfile[50];
	char outputfile1[50]="./globaloutput/meshgrid_vtk.vtk";
	char outputfile2[50]="./globaloutput/datapoints_vtk.vtk";

	int row, col, k, n,leftphimc,rightphimc,leftphipc,rightphipc;
	double rho,pressure,gamma,phim,phip,leftphim,rightphim,leftphip,rightphip,cminus,cplus;
	ofstream output1(outputfile1);
	ofstream output2(outputfile2);

	int characterspace=picnum/(characterdensity+1);
	int oldcharacters=0;
	
		output1<<"# vtk DataFile Version 2.0"<<endl;
		output1<<"An example vtk file for Wentian"<<endl;
		output1<<"ASCII"<<endl;
		output1<<"DATASET STRUCTURED_GRID"<<endl;
		output1<<"DIMENSIONS "<<celln<<" "<<picnum<<" "<<"1"<<endl;
		output1<<"POINTS "<<celln*picnum<<" double"<<endl;
		
		
		output2<<"POINT_DATA "<<celln*picnum<<endl;
		output2<<"SCALARS Temperature[K] double 1"<<endl;
		output2<<"LOOKUP_TABLE default"<<endl;
	for(int i=0;i<picnum;i++)
	{
		n=sprintf (inputfile, "./output/%d.csv", i);
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
		double characterdataminus[1000];
		double characterdataplus[1000];
		
		data.resize(row);

		
		for(k = 0; k < row; k ++)
			data[k].resize(col);
	 
		ifstream infile2(inputfile);

		int trow, tcol;
		trow = tcol = 0;
		int characters=0;

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
		//数据分割处理
		//(1) 算时间
		tm=data[0][col-1];
		//(2)分别载入密度，速度，压强以及温度数据
		for(int j=0;j<row;j++)
		{
			output1<<data[j][0]<<" "<<tm<<" "<<"0"<<endl;
			output2<<data[j][5]<<endl;
		}
		cout<<"Progress in "<<i<<"/"<<picnum<<endl;

		
		
	}
	output1<<endl;
	output1.close();
	output2.close();
	system("cat ./globaloutput/meshgrid_vtk.vtk ./globaloutput/datapoints_vtk.vtk > ./globaloutput/output.vtk");
	return 0;
}
