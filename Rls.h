#ifndef RLS_H_
#define RLS_H_

#include <vector>
#include <algorithm>

using namespace std;

class Rls
{
public:
	Rls(): N(0),lambda(1),X(nullptr),cov(nullptr) {}

	Rls(unsigned N,double sigma2,double lambda): X(nullptr),cov(nullptr)
	{
		Reset(N,sigma2,lambda);
	}

	virtual ~Rls()
	{
		freeAll();
	}

	void Reset(unsigned N, double sigma2,double lambda)
	{
		freeAll();
		this->N=N;
		this->lambda=lambda;

		X=new double[N];
		fill(X,X+N,0);

		cov=new double*[N];
		for (unsigned i=0; i<N; i++)
		{
			cov[i]=new double[N];
			fill(cov[i],cov[i]+N,0);
			cov[i][i]=sigma2;
		}
	}

	void Step(const double S[], double yd);

	double* getX() const { return X; }
	unsigned getN() const { return N; }
	double getLambda() const { return lambda; }





private:

	unsigned N;
	double lambda;
	double *X;
	double **cov;


	inline void KalmanGain(double K[], const double S[]);
	inline void MatrixMul(const double K[],const double S[]);

	void freeAll()
	{
		if(X)
			delete[] X;

		if(cov)
		{
			for(unsigned i=0;i<N;i++)
				if(cov[i])
					delete[] cov[i];

			delete[] cov;
		}
	}



};

#endif /* RLS_H_ */
