#include <iomanip>
#include <iostream>
#include <fstream>
#include<sstream>

#include <cstdlib>

#include <cmath>
//#include "phylobayes.h"

using namespace std;

const double Pi = 3.14155927;
const double defaultCI=95;

enum	Switch		{No = 0, Yes = 1};

class Correlation{
 private:

  string chainName;
  int burnin;
  long int nbsample;
  int nbparameter;
  double **parameters;
  double **sortparameters;
  double **covparam;
  double **covnorm;
  double *variance;
  double *meanparam;
  double *weight;
  double *effectiveSize;
  double CI;
  double **CIbuffer;
  string *parameterName;

  Switch *isConstant;
  Switch isSort;

  void createParameterBuffer();
  void createCovarianceBuffer();
  void createWeight();
  void computeMean();
  void init();


  void sortParameters();
  void reccurquicksort(double *to, double *buf, int deb, int fin);
  void quicksort(double *from, double *to, int s);
  void getCI(double *distrib,int s, double ci, double *inf, double* sup);
 public:
  Correlation(double ci=-1);
  ~Correlation();
  void getParameters(string chainName,int start,int stop);
  void computeCovariance();
  void computeWeight();
  void computeEffectiveSize();
  void outputResult();
  void computeCI(double ci=-1);
  int getNbParameter();

  string getParameterName(int n);
  double getVariance(int n);
  double getEffectiveSize(int n);
  double getMean(int n);
  double getNbSample();
  double getInfCI(int n);
  double getSupCI(int n);

};
//#endif
