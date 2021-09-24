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
	double tm;
//	const double up=1500;
	const double dphi=1.5e-2/5000;
	int picnum=2000;
	const double dt=1.0e-4/picnum;
	const int characterdensity=10;
	char inputfile[50];
	char outputfile1[50]="./globaloutput/output_uprT.csv";
	char outputfile2[50]="./globaloutput/incidentshock.csv";
	char outputfile3[50]="./globaloutput/mcharacters.csv";
	char outputfile4[50]="./globaloutput/pcharacters.csv";
	int row, col, k, n,leftphimc,rightphimc,leftphipc,rightphipc;
	double rho,pressure,gamma,phim,phip,leftphim,rightphim,leftphip,rightphip,cminus,cplus;
	ofstream output1(outputfile1);
	ofstream output2(outputfile2);
	ofstream output3(outputfile3);
	ofstream output4(outputfile4);
	int characterspace=picnum/(characterdensity+1);
	int oldcharacters=0;
	
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
			output1<<tm<<","<<j*dphi<<","<<data[j][0]<<","<<data[j][1]<<","<<data[j][2]<<","<<data[j][3]<<","<<data[j][5]<<","<<data[j][6]<<","<<data[j][8]<<endl;
		}
		output1<<endl;
		cout<<"Progress in "<<i<<"/"<<picnum<<endl;
		//Calculate Shock position
		
		int track1=row-4;
		int shockfind=0;
		
		while(track1>4 && shockfind==0)
		{
			if(data[track1-3][1]>data[track1+3][1]*2)
			{
				shockfind=1;
				output2<<tm<<","<<track1*dphi<<endl;
			}
			track1=track1-1;
		}
		
		if(shockfind==0)
		{
			track1=row-4;
		}
		//Calculate Character lines
		if(i>0)
		{
			
			characters=i/characterspace;
			cout<<"Progress on "<<characters<<"/"<<characterdensity<<endl;
			
			
			//New characters
			if(characters>oldcharacters)
			{
				characterdataminus[characters-1] = track1*dphi;
				characterdataplus[characters-1] = 0;
			}
				
			//Update old characters
			if(characters>0)
			{

				for(int characterlinetrack=0; characterlinetrack<oldcharacters; characterlinetrack++)
				{
					
					phim=characterdataminus[characterlinetrack];
					phip=characterdataplus[characterlinetrack];

					//C-
					if(phim>dphi)
					{
						for(int trackm=0;trackm<row-4;trackm++)
						{
							if((trackm*dphi<=phim) && ((trackm+1)*dphi>phim))
							{
								leftphimc=trackm;
								rightphimc=trackm+1;
								leftphim=trackm*dphi;
								rightphim=(trackm+1)*dphi;

							}
							
						}
						//mean variable
						rho=data[leftphimc][1]+((data[rightphimc][1]-data[leftphimc][1])*( (phim-leftphim)/(rightphim-leftphim) ));
						pressure=data[leftphimc][3]+((data[rightphimc][3]-data[leftphimc][3])*( (phim-leftphim)/(rightphim-leftphim) ));
						gamma=data[leftphimc][6]+((data[rightphimc][6]-data[leftphimc][6])*( (phim-leftphim)/(rightphim-leftphim) ));
						cminus=pow(rho*pressure*gamma,0.5);
						

						characterdataminus[characterlinetrack]=characterdataminus[characterlinetrack]-cminus*dt;
						cout<<"phim="<<cminus<<endl;
						cout<<"phim="<<dt<<endl;
					}
					else
					{
						characterdataminus[characterlinetrack]=0.0/0.0;
					}
					
					//C+
					
					if(phip<(track1-3)*dphi)
					{
						for(int trackp=0;trackp<row-4;trackp++)
						{
							if((trackp*dphi<=phip) && ((trackp+1)*dphi>phip))
							{
								leftphipc=trackp;
								rightphipc=trackp+1;
								leftphip=trackp*dphi;
								rightphip=(trackp+1)*dphi;

							}
							
						}
						//mean variable
						rho=data[leftphipc][1]+((data[rightphipc][1]-data[leftphipc][1])*( (phip-leftphip)/(rightphip-leftphip) ));
						pressure=data[leftphipc][3]+((data[rightphipc][3]-data[leftphipc][3])*( (phip-leftphip)/(rightphip-leftphip) ));
						gamma=data[leftphipc][6]+((data[rightphipc][6]-data[leftphipc][6])*( (phip-leftphip)/(rightphip-leftphip) ));
						cplus=pow(rho*pressure*gamma,0.5);
						

						characterdataplus[characterlinetrack]=characterdataplus[characterlinetrack]+cplus*dt;
					}
					else
					{
						characterdataplus[characterlinetrack]=0.0/0.0;
					}
					
					
				}


				
				output3<<tm;
				output4<<tm;
				
				for(int characterlinetrack=0; characterlinetrack<characters; characterlinetrack++)
				{				
					output3<<","<<characterdataminus[characterlinetrack];
					output4<<","<<characterdataplus[characterlinetrack];
				}
				output3<<endl;
				output4<<endl;
			}
		}

		oldcharacters=characters;
	}
	output1.close();
	output2.close();
	output3.close();
	output4.close();
	return 0;
}
