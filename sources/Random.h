#ifndef RANDOM_H
#define RANDOM_H

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>

#define MT_LEN       624

#define MT_IA           397
#define MT_IB           (MT_LEN - MT_IA)
#define UPPER_MASK      0x80000000
#define LOWER_MASK      0x7FFFFFFF
#define MATRIX_A        0x9908B0DF
#define TWIST(b,i,j)    ((b)[i] & UPPER_MASK) | ((b)[j] & LOWER_MASK)
#define MAGIC(s)        (((s)&1)*MATRIX_A)


using namespace std;

const double gammacoefs[] = {0.9999999999995183,676.5203681218835,-1259.139216722289,771.3234287757674,-176.6150291498386,12.50734324009056,-0.1385710331296526,0.9934937113930748e-05,0.1659470187408462e-06};
static const double Logroot2pi =0.918938533204673;

class Random {
 
	public:

	Random();
  	Random(int seed = -1);

	static void InitRandom(int seed = -1);

	static int GetSeed();

	static double Uniform();
	static double Gamma(double alpha,double beta);
  	static double sNormal(void);
  	static double sExpo(void);
  	static double sGamma(double);
  	static double sGammanew(double);

 	static int Choose(int);
	static int FiniteDiscrete(int n, const double* probarray);
  	static void DrawFromUrn(int*, int n, int N);
	static int DrawFromDiscreteDistribution(double* p, int n);

	static double logGamma(double a);
	static double Bessel_I(double a,  double b);

	static int rrindex(int i, int j, int N);

	private:

	static int Seed;
	static int mt_index;
	static unsigned long mt_buffer[MT_LEN];

  };

#endif // RANDOM_H
