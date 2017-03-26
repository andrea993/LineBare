#include "Rls.h"
#include <assert.h>
#include <iostream>

void Rls::Step(const double S[], double yd)
{
	double yest=0;
	for(unsigned i=0;i<N;i++)
		yest+=S[i]*X[i];

	double e=yd-yest;

	double *K = new double[N];
	KalmanGain(K,S);

	for(unsigned i=0;i<N;i++)
		X[i]+=K[i]*e;


	MatrixMul(K,S);

	delete[] K;

}


void Rls::MatrixMul(const double K[],const double S[])
{
	double *tmp = new double[N];
	double V_jk;
	for (unsigned i=0; i<N; i++)
	{

		for(unsigned k=0; k<N; k++)
			tmp[k]=cov[k][i]/lambda;

		for (unsigned j=0; j<N; j++)
		{
			cov[j][i]=0;
			for (unsigned k=0; k<N; k++)
			{
				V_jk = j!=k ? (-K[j]*S[k]) : (1-K[j]*S[k]);
				cov[j][i]+=V_jk*tmp[k];
			}

		}
	}
	delete[] tmp;
}

void Rls::KalmanGain(double K[],const double S[])
{
	double den=lambda;

	for(unsigned i=0;i<N;i++)
	{
		double a=0;
		for(unsigned j=0;j<N;j++)
			a+=S[j]*cov[j][i];

		den+=S[i]*a;
	}

	for(unsigned i=0;i<N;i++)
	{
		double num=0;
		for(unsigned j=0;j<N;j++)
			num+=S[j]*cov[i][j];


		K[i]=num/den;
	}
}
