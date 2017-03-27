#include <iostream>
#include <fstream>
#include "Rls.h"
using namespace std;

int main(int argc, char **argv)
{
	ifstream fi("dataTest.csv");

	vector<double> inp;
	double v;
	while(fi>>v)
		inp.push_back(v);


	int N=600;
	Rls rls(N+1,1e14,1);

	int k0=0, k1=N;

	double *S_k= new double[N+1];

	while(k1<inp.size())
	{
		double y_k=inp[k1];
		for(int i=0;i<N;i++)
		{
			S_k[i]=inp[k1-i-1];
		}
		S_k[N]=1;

		rls.Step(S_k,y_k);

		cout<<k0<<endl;

		k0++;
		k1++;
	}

	delete[] S_k;


  for(int i=0;i<rls.getN();i++)
		cout<<rls.getX()[i]<<endl;



	return 0;
}
