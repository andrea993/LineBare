#include "Rls.h"
#include <assert.h>
#include <iostream>

void Rls::Step(const double S[], double yd)
{
	double yest=0;
	for(int i=0;i<N;i++)
		yest+=S[i]*X[i];

	double e=yd-yest;

	double K[N];
	KalmanGain(K,S);

	for(int i=0;i<N;i++)
		X[i]+=K[i]*e;


	MatrixMul(K,S);


}


void Rls::MatrixMul(const double K[],const double S[])
{
	double tmp[N];
	double V_jk;
	for (int i=0; i<N; i++)
	{

		for(int k=0; k<N; k++)
			tmp[k]=cov[k][i]/lambda;

		for (int j=0; j<N; j++)
		{
			cov[j][i]=0;
			for (int k=0; k<N; k++)
			{
				V_jk = j!=k ? (-K[j]*S[k]) : (1-K[j]*S[k]);
				cov[j][i]+=V_jk*tmp[k];
			}

		}
	}
}

void Rls::KalmanGain(double K[],const double S[])
{
	double den=lambda;

	for(int i=0;i<N;i++)
	{
		double a=0;
		for(int j=0;j<N;j++)
			a+=S[j]*cov[j][i];

		den+=S[i]*a;
	}

	for(int i=0;i<N;i++)
	{
		double num=0;
		for(int j=0;j<N;j++)
			num+=S[j]*cov[i][j];


		K[i]=num/den;
	}
}
