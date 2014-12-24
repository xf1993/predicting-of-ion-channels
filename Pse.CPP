//---------伪氨基酸程序特征提取部分（by夏凤）--------
#include<iostream>
#include<string>
#include<cstring>
#include<fstream>
#include<iomanip>
#include<cmath>
#include<vector>
#include<cstdlib>

using namespace std;
using std::vector;
using std::string;

#define k 2						//-----------理化因子个数--------
#define G 15					//-----------阶数--------
#define num1 302


int main()
{
	char AA[20]={'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'};
	static double H0[k][20],H[k][20];
	//---------从文件中导入理化参数-----------
	ifstream data1("./Dataset/H0.txt");
	int i=0,j=0;
	for(i=0;i<20;i++)
		for(j=0;j<k;j++)
		{
			data1>>H0[j][i];
			//cout<<H0[j][i]<<endl;getchar();
		}
	data1.close();
	//-----------------理化因子进行标准化------------
	double h0[k]={0},s2[k]={0};
	for(j=0;j<k;j++)
	{
		for(i=0;i<20;i++) {h0[j]+=H0[j][i];}
		h0[j]/=20;
	}
	for(j=0;j<k;j++)
	{
		for(i=0;i<20;i++)
		{s2[j]+=((H0[j][i]-h0[j])*(H0[j][i]-h0[j]));}
		s2[j]/=20;
	}
	for(i=0;i<20;i++)
		for(j=0;j<k;j++)
		{H[j][i]=(H0[j][i]-h0[j])/sqrt(s2[j]);}
	ofstream data2("./Dataset/H.txt");
	for(i=0;i<20;i++)
	{
		for(j=0;j<k;j++)
		{data2<<H[j][i]<<'\t';}
		data2<<endl;
	}
	data2.close();
	//----------------读取文件中的蛋白质序列------------
	int len1[num1/2]={0};
	string No1[num1/2],protien_1[num1/2];
	ifstream intext1("./Dataset/456.txt");
	if(!intext1)
	{
		cout<<"can't open 456.txt."<<endl;
		exit(1);
	}


	ofstream train1("./Dataset/train1.txt");
	int p=0,q=0;
	for(i=0,p=0,q=0;i<num1;i++)
	{
		if(i%2==0){getline(intext1,No1[p]);train1<<No1[p]<<endl;p++;continue;}
		else
		{
			getline(intext1,protien_1[q]);
			len1[q]=protien_1[q].length();
			train1<<protien_1[q]<<endl;
			q++;
		}
	}
	intext1.close();
	train1.close();


	//--------------构造容器计算关联值（1~G阶）-------------
	vector<double> seq[k],CF;
	double cf=0,f=0;
	int g;
	ofstream F1("./Result/f1.txt");
	F1<<fixed<<"这是1-"<<G<<"阶数据："<<'\n';
	for(q=0;q<num1/2;q++)
	{
		for(i=0;i<len1[q];i++)
			for(j=0;j<20;j++)
				if(protien_1[q][i]==AA[j])
				{
					for(p=0;p<k;p++) seq[p].push_back(H[p][j]);
					break;
				}
		for(g=1;g<=G;g++)
		{
			for(i=0;i<len1[q]-g;i++)
			{
				for(p=0;p<k;p++) cf+=((seq[p][i+g]-seq[p][i])*(seq[p][i+g]-seq[p][i]));
				CF.push_back(cf/k);
				f+=CF[i];
				cf=0;
			}
			f/=CF.size();
			F1<<f<<'\t';
			f=0;
			CF.clear();
		}
		F1<<endl;
		for(p=0;p<k;p++) seq[p].clear();
	}
	F1.close();
	return 0;
}