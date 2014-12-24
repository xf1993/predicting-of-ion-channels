//--------伪氨基酸程序分类预测部分（by夏凤）--------
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

#define row 404
#define PC 20
#define PseAAC 15
#define w1 0.05
#define w2 0.96

int main()
{
	int p,q;
	double f[row][PC],g[row][PseAAC],
			F[row]={0},G[row]={0};
	ifstream data1("./Dataset/AAC.txt");
	if(!data1)
	{cout<<"can't open the file AAC.txt."<<endl;
	exit(1);}
	ifstream data2("./Dataset/PseAAC.txt");
	if(!data2)
	{cout<<"can't open the file PseAAC.txt."<<endl;
	exit(1);}
	for(p=0;p<row;p++)
	{
		for(q=0;q<PC;q++)
		{
			data1>>f[p][q];
			F[p]+=f[p][q];
		}
		//F[p]/=PC;
		for(q=0;q<PseAAC;q++)
		{
			data2>>g[p][q];
		}cout<<g[p][q-1]<<endl;getchar();
	}
	data1.close();
	data2.close();
	double W;
	vector<double> X;
	for(int t=1;t<=10;t++)  //-----------10jier-------
	{
		for(double w=w1;w<=w2;w=w+0.1)
		{
			ofstream out("./out");
			cout<<t<<"-jier,w="<<w<<":\n";
			for(p=0;p<row;p++)
			{
				for(q=0;q<t;q++)
				{
					G[p]+=g[p][q];
				}
				//G[p]/=t;
				W=1.0+w*G[p];									//---W=F[p]+w*G[p];
				for(q=0;q<PC;q++) X.push_back(f[p][q]/W);
				for(q=0;q<t;q++) X.push_back(w*g[p][q]/W);     //----------------g[p][q]/W
				if(p<178) out<<fixed<<"1 ";
				else out<<fixed<<"2 ";
				for(q=0;q<X.size();q++)
				{
					out<<q+1<<':'<<X[q]<<' ';
				}
				out<<endl;
				X.clear();
				G[p]=0;
			}
			out.close();
			//------------------调用SVM----------------
			cout<<"svm-scale -l 0 -u 1 -s ./Result/rule ./out>./out.scale"<<endl;
			system("svm-scale -l 0 -u 1 -s ./Result/rule ./out>./out.scale");
			cout<<"tools/grid.py -v 404 ./out.scale"<<endl;
			system("tools/grid.py -v 404 ./out.scale");			//----------------mark,折数修改
			ofstream grid("./No.txt",ios_base::app);
			grid<<t<<"-jier,w="<<w<<endl;
			grid.close();
			cout<<endl;
		}
		cout<<"end of "<<t<<"-jier.press enter to continue";//getchar();
	}
	return 0;
}