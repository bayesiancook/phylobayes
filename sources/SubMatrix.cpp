#include "phylo.h"
#include "linalg.h"


// ---------------------------------------------------------------------------
//		 SubMatrix()
// ---------------------------------------------------------------------------

void SubMatrix::BasicAllocation()	{

	mArray = new double[4 * NstateNstate + 2 * Nstate];

	Q = new double*[Nstate];
	for (int i=0; i<Nstate; i++)	{
		Q[i] = mArray + i*Nstate;
	}

	a = new double*[Nstate];
	for (int i=0; i<Nstate; i++)	{
		a[i] = mArray + (Nstate + i) * Nstate;
	}

	u = new double*[Nstate];
	for (int i=0; i<Nstate; i++)	{
		u[i] = mArray + (2*Nstate + i) * Nstate;
	}

	invu = new double*[Nstate];
	for (int i=0; i<Nstate; i++)	{
		invu[i] = mArray + (3 * Nstate + i) * Nstate;
	}

	v = mArray + 4 * NstateNstate;
	vi = mArray + 4 * NstateNstate + Nstate;

}

SubMatrix::SubMatrix(MCParameters* inParam, int inNstate, double* rr, double* stat)	{

	mStationary0 = 0;
	mParam = inParam;
	if (mParam)	{
		Normalise = mParam->Normalise;
	}
	else	{
		Normalise = No;
	}
	Nstate = inNstate;
	NstateNstate = Nstate * Nstate;

	if (stat)	{
		mStationary = stat;
		StatAlloc = 0;
	}
	else	{
		mStationary = new double[Nstate];
		for (int i=0; i<Nstate; i++)	{
			mStationary[i] = 1.0 / Nstate;
		}
		StatAlloc = 1;
	}
	mRelativeRate = rr;
	RelRateAlloc = 0;
	MatrixStationary = mStationary;
	MatrixStatAlloc = 0;

	mPow = 0;
	xi = 0;
	invprob = 0;
	rate = 0;
	HeteroMode = Homo;
	BasicAllocation();
}

SubMatrix::SubMatrix(MCParameters* inParam, int inNstate, double* rr, double* stat0, double* stat, double* f, double a)	{

	if (! inParam)	{
		cerr << "error in subMatrix creation\n";
		exit(1);
		// Normalise = No;
	}
	mParam = inParam;
	Normalise = mParam->Normalise;
	Nstate = inNstate;
	NstateNstate = Nstate * Nstate;
	GWf = f;

	if (mParam->MutMode == 0)	{
		mStationary0 = 0;
		if (stat)	{
			mStationary = stat;
			StatAlloc = 0;
			Stat0Alloc = 0;
		}
		else	{
			mStationary = new double[Nstate];
			for (int i=0; i<Nstate; i++)	{
				mStationary[i] = 1.0 / Nstate;
			}
			StatAlloc = 1;
			Stat0Alloc = 0;
		}
		mRelativeRate = rr;
		RelRateAlloc = 0;
	}
	else if (mParam->MutMode == 1)	{
		kappa = 0;
		mStationary = stat;
		mStationary0 = stat0;
		StatAlloc = 0;
		Stat0Alloc = 0;
		mRelativeRate = rr;
		RelRateAlloc = 0;
	}
	else if (mParam->MutMode == 2)	{
		mStationary = stat;
		kappa = rr;
		mRR_ACGT = 0;
		acgt = stat0;
		mStationary0 = new double[Nstate];
		Stat0Alloc = 1;
		StatAlloc = 0;
		mRelativeRate = new double[Nstate*(Nstate-1)/2];
		RelRateAlloc = 1;
		MatrixStationary = new double[Nstate];
		MatrixStatAlloc = 1;
	}
	else if (mParam->MutMode == 3)	{
		mStationary = stat;
		mRR_ACGT = rr;
		kappa = 0;
		acgt = stat0;
		mStationary0 = new double[Nstate];
		Stat0Alloc = 1;
		StatAlloc = 0;
		mRelativeRate = new double[Nstate*(Nstate-1)/2];
		RelRateAlloc = 1;
		MatrixStationary = new double[Nstate];
		MatrixStatAlloc = 1;
	}

	if (mParam->SelMode)	{
		MatrixStationary = new double[Nstate];
		MatrixStatAlloc = 1;
	}
	else	{
		MatrixStationary = mStationary;
		MatrixStatAlloc = 0;
	}

	mPow = 0;
	xi = 0;
	invprob = 0;
	rate = 0;
	HeteroMode = Homo;
	BasicAllocation();
}

SubMatrix::SubMatrix(MCParameters* inParam, int inNstate, double* inrate, double* ininvprob, double* inxi, double* inStationary)	{

	mParam = inParam;
	mStationary0 = 0;
	Normalise = No;
	HeteroMode = Covarion;
	NormaliseCov = No;
	if (mParam)	{
		NormaliseCov = mParam->NormaliseCov;
	}

	Nstate = 2 * inNstate;
	NstateNstate = Nstate * Nstate;

	StatAlloc = 0;
	Stat0Alloc = 0;
	mStationary = inStationary;

	MatrixStatAlloc = 1;
	MatrixStationary = new double[Nstate];

	RelRateAlloc = 0;
	mRelativeRate = 0;

	xi = inxi;
	invprob = ininvprob;
	rate = inrate;

	mPow = 0;
	BasicAllocation();
}

SubMatrix::SubMatrix(MCParameters* inParam, int inNstate, double* inrate, double* ininvprob, double* inxi, double* rr, double* stat)	{

	mParam = inParam;
	mStationary0 = 0;
	Normalise = No;
	HeteroMode = Covarion;

	Nstate = 2 * inNstate;
	NstateNstate = Nstate * Nstate;

	StatAlloc = 0;
	Stat0Alloc = 0;
	mStationary = stat;
	
	MatrixStatAlloc = 1;
	MatrixStationary = new double[Nstate];

	RelRateAlloc = 0;
	mRelativeRate = rr;

	xi = inxi;
	invprob = ininvprob;
	rate = inrate;

	mPow = 0;
	BasicAllocation();
}

void SubMatrix::CreatePowers()	{

	if (mParam && mParam->UniSub)	{
		mPow = new double**[mParam->UniSubNmax];
		for (int n=0; n<mParam->UniSubNmax; n++)	{
			mPow[n] = 0;
		}
	}
	npow = 0;
}

void SubMatrix::CreatePowers(int n)	{

	mPow[n] = new double*[Nstate];
	for (int i=0; i<Nstate; i++)	{
		mPow[n][i] = new double[Nstate];
		for (int j=0; j<Nstate; j++)	{
			mPow[n][i][j] = 0;
		}
	}
}

void SubMatrix::DeletePowers()	{

	if (mPow)	{
		for (int n=0; n<mParam->UniSubNmax; n++)	{
			if (mPow[n])	{
				for (int i=0; i<Nstate; i++)	{
					delete[] mPow[n][i];
				}
				delete [] mPow[n];
			}
		}
		delete[] mPow;
		mPow = 0;
	}
}

// ---------------------------------------------------------------------------
//		 void Copy(SubMatrix* from)
// ---------------------------------------------------------------------------

void	SubMatrix::Copy( SubMatrix* from)	{

	if (HeteroMode == Covarion)	{
		for (int i=0; i<Nstate/2; i++)	{
			mStationary[i] = from->mStationary[i];
		}
		for (int i=0; i<Nstate; i++)	{
			MatrixStationary[i] = from->MatrixStationary[i];
		}
		if (mRelativeRate)	{
			for (int i=0; i<Nstate/2-1; i++)	{
				for (int j=i+1; j<Nstate/2; j++)	{
					RelativeRate(i,j) = from->RelativeRate(i,j);
				}
			}
		}
	}
	else	{
		for (int i=0; i<Nstate; i++)	{
			mStationary[i] = from->mStationary[i];
		}
		if (mRelativeRate)	{
			for (int i=0; i<Nstate-1; i++)	{
				for (int j=i+1; j<Nstate; j++)	{
					RelativeRate(i,j) = from->RelativeRate(i,j);
				}
			}
		}
		if (mStationary0)	{
			for (int i=0; i<Nstate; i++)	{
				mStationary0[i] = from->mStationary0[i];
				MatrixStationary[i] = from->MatrixStationary[i];
			}
		}
	}

	for (int i=0; i< 4 * NstateNstate + 2 * Nstate ; i++)	{
		mArray[i] = from->mArray[i];
	}
	UniMu = from->UniMu;
}

// ---------------------------------------------------------------------------
//		 void ScalarMul()
// ---------------------------------------------------------------------------

void	SubMatrix::ScalarMul(double e)	{

	for (int i=0; i< Nstate ; i++)	{
		v[i] *= e;
		vi[i] *= e;
	}
	UniMu *= e;
}


// ---------------------------------------------------------------------------
//		 ~SubMatrix()
// ---------------------------------------------------------------------------

SubMatrix::~SubMatrix()	{

	if (RelRateAlloc)	{
		delete[] mRelativeRate;
	}
	if (StatAlloc)	{
		delete[] mStationary;
	}
	if (Stat0Alloc)	{
		delete[] mStationary0;
	}
	if (MatrixStatAlloc)	{
		delete[] MatrixStationary;
	}

	delete[] mArray;
	delete[] Q;
	delete[] a;
	delete[] u;
	delete[] invu;

	DeletePowers();
}

// ---------------------------------------------------------------------------
//		 void ComputeExponential(double range, double** expo)
// ---------------------------------------------------------------------------

void	SubMatrix::ComputeExponential(double range, double** expo)	{

	// ComputeArray();
	
	double** temp = new double*[Nstate];
	for (int i=0; i<Nstate; i++)	{
		temp[i] = new double[Nstate];
		for (int j=0; j<Nstate; j++)	{
			expo[i][j] = 0;
		}
	}

	for (int i=0; i<Nstate; i++)	{
		expo[i][i] = exp(range * v[i]);
	}

	for (int i=0; i<Nstate; i++)	{
		for (int j=0; j<Nstate; j++)	{
			double t = 0;
			for (int k=0; k<Nstate; k++)	{
				t += expo[i][k] * invu[k][j];
			}
			temp[i][j] = t;
		}
	}

	for (int i=0; i<Nstate; i++)	{
		for (int j=0; j<Nstate; j++)	{
			double t = 0;
			for (int k=0; k<Nstate; k++)	{
				t += u[i][k] *temp[k][j];
			}
			expo[i][j] = t;
		}
	}

	// deleting temp array
	for (int i=0; i<Nstate; i++)	{
		delete[] temp[i];
	}
	delete[] temp;

}


// ---------------------------------------------------------------------------
//		 double& operator()(int ,int)
// ---------------------------------------------------------------------------

double&	SubMatrix::operator()(int i, int j)	{
	return mArray[Nstate * j + i];
}
double&	SubMatrix::operator()(int i, int j)	const {
	return mArray[Nstate * j + i];
}

// ---------------------------------------------------------------------------
//		 double& Stationary(int i)
// ---------------------------------------------------------------------------

double& SubMatrix::Stationary(int i)	{
	
	return MatrixStationary[i];
}


// ---------------------------------------------------------------------------
//		 double& RelativeRate(int i, int j)
// ---------------------------------------------------------------------------

 double& SubMatrix::RelativeRate(int i, int j)	{

	if (! mRelativeRate)	{
		cerr << "error : cannot call rel rate under cov poisson\n";
		exit(1);
	}
	if (HeteroMode == Covarion)	{
	return (i<j) ?
		mRelativeRate[(Nstate - i - 1) * i / 2 + j - i - 1] :
		mRelativeRate[(Nstate - j - 1) * j / 2 + i - j - 1] ;
	}
	return (i<j) ?
		mRelativeRate[(2 * Nstate - i - 1) * i / 2 + j - i - 1] :
		mRelativeRate[(2 * Nstate - j - 1) * j / 2 + i - j - 1] ;
}

// ---------------------------------------------------------------------------
//		 void ComputeArray()
// ---------------------------------------------------------------------------


int	SubMatrix::ComputeArray()	{

	if (HeteroMode == Covarion)	{
		ComputeArrayHetero();
	}
	else if (mParam->MutMode)	{
		ComputeArrayHomoMutSel();
	}
	else	{
		ComputeArrayHomo();
	}
	int tmp = Diagonalise();
	return tmp;
}

int SubMatrix::ComputeArrayHomo()	{

	double norm = 1;
	if (Normalise)	{
		norm = ComputeRate();
	}

	if (mRelativeRate)	{
		for (int i=0; i<Nstate; i++)	{
			for (int j=0; j<Nstate; j++)	{
				if (i==j)	{
					double temp = 0;
					for (int k=0; k<Nstate; k++)	{
						if (k != i)	{
							temp -= RelativeRate(i,k) * Stationary(k);
						}
					}
					operator()(i,j) = temp;
				}
				else	{
					operator()(i,j) = RelativeRate(i,j) * Stationary(j);
				}
			}
		}
	}
	else	{
		for (int i=0; i<Nstate; i++)	{
			for (int j=0; j<Nstate; j++)	{
				if (i==j)	{
					double temp = 0;
					for (int k=0; k<Nstate; k++)	{
						if (k != i)	{
							temp -= Stationary(k);
						}
					}
					operator()(i,j) = temp;
				}
				else	{
					operator()(i,j) = Stationary(j);
				}
			}
		}
	}
	if (Normalise)	{
		for (int i=0; i<Nstate; i++)	{
			for (int j=0; j<Nstate; j++)	{
				operator()(i,j) /= norm;
			}
		}
	}

	return 0;
}

void SubMatrix::ComputeGWStationary()	{

	if (mParam->SelMode != 1)	{
		cerr << "error in submatrix compute gw stat\n";
		exit(1);
	}
	double total = 0;
	for (int k=0; k<Nstate; k++)	{
		MatrixStationary[k] = mStationary0[k] * mStationary[k];
		total += MatrixStationary[k];
	}
	for (int k=0; k<Nstate; k++)	{
		MatrixStationary[k] /= total;
	}
}

int SubMatrix::IsTransition(int a, int b)	{

	if ((a == 0) && (b == 2))	{
		return 1;
	}
	if ((a == 2) && (b == 0))	{
		return 1;
	}
	if ((a == 1) && (b == 3))	{
		return 1;
	}
	if ((a == 3) && (b == 1))	{
		return 1;
	}
	return 0;
}

int SubMatrix::position(int i, int j)	{

	// if two codons differ at position 0, 1 or 2	

	// one of them is stop codon
	if (mParam->CodonCode[i] == -1)	{
		return -1;
	}
	if (mParam->CodonCode[j] == -1)	{
		return -1;
	}
	
	// identical
	if ((codonpos[0][i] == codonpos[0][j]) && (codonpos[1][i] == codonpos[1][j]) && (codonpos[2][i] == codonpos[2][j]))	{
		return -1;
	}
	if (codonpos[0][i] != codonpos[0][j])	{
		if ((codonpos[1][i] == codonpos[1][j]) && (codonpos[2][i] == codonpos[2][j]))	{
			return 0;
		}
		else	{
			return 3;
		}
	}
	if (codonpos[1][i] != codonpos[1][j])	{
		if ((codonpos[0][i] == codonpos[0][j]) && (codonpos[2][i] == codonpos[2][j]))	{
			return 1;
		}
		else	{
			return 3;
		}
	}
	if (codonpos[2][i] != codonpos[2][j])	{
		if ((codonpos[0][i] == codonpos[0][j]) && (codonpos[1][i] == codonpos[1][j]))	{
			return 2;
		}
		else	{
			return 3;
		}
	}
	return 3;
}

double SubMatrix::RR_ACGT(int i, int j)	{

	return (i<j) ?
		mRR_ACGT[(2 * Nnuc - i - 1) * i / 2 + j - i - 1] :
		mRR_ACGT[(2 * Nnuc - j - 1) * j / 2 + i - j - 1] ;

}

void SubMatrix::ComputeCodonProjection()	{

	double codonmatrix[Ncodon][Ncodon];
	for (int i=0; i<Ncodon; i++)	{
		for (int j=0; j<Ncodon; j++)	{
			double rate = 0;
			int pos = position(i,j);
			if ((pos != -1) && (pos != 3))	{
				int a = codonpos[pos][i];
				int b = codonpos[pos][j];
				if (a == b)	{
					cerr << "error in SubMatrix::ComputeCodonProjection\n";
					exit(1);
				}
				if (mParam->MutMode == 2)	{
					if (IsTransition(a,b))	{
						rate = *kappa;
					}
					else	{
						rate = 1;
					}
				}
				else	{
					rate = RR_ACGT(a,b);
				}
				rate *= acgt[b];
				codonmatrix[i][j] = rate;
			}
			else	{
				codonmatrix[i][j] = 0;
			}
		}
	}
	double total = 0;
	double codonstat[Ncodon];
	for (int i=0; i<Ncodon; i++)	{
		if (mParam->CodonCode[i] == -1)	{
			codonstat[i] = 0;
		}
		else	{ 
			codonstat[i] = acgt[codonpos[0][i]] * acgt[codonpos[1][i]] * acgt[codonpos[2][i]];
		}
		total += codonstat[i];
	}
	for (int i=0; i<Ncodon; i++)	{
		codonstat[i] /= total;
	}
	total = 0;
	for (int i=0; i<Nstate; i++)	{
		mStationary0[i] = 0;
	}
	for (int i=0; i<Ncodon; i++)	{
		mStationary0[mParam->CodonCode[i]] += codonstat[i];
	}
	double totalrr = 0;
	for (int l=0; l<Nstate; l++)	{
		for (int m=l+1; m<Nstate; m++)	{
			double tot = 0;
			for (int i=0; i<Ncodon; i++)	{
				if (mParam->CodonCode[i] == l)	{
					for (int j=0; j<Ncodon; j++)	{
						if (mParam->CodonCode[j] == m)	{
							tot += codonstat[i] * codonmatrix[i][j];
						}
					}
				}
			}
			tot /= mStationary0[l] * mStationary0[m];
			RelativeRate(l,m) = tot;
			totalrr += tot;
		}
	}
	/*
	for (int l=0; l<Nstate; l++)	{
		for (int m=l+1; m<Nstate; m++)	{
			cout << RelativeRate(l,m) / totalrr << '\n';
		}
	}
	exit(1);
	*/
}

double SubMatrix::SelectionFactor(int i, int j)	{

	if (mParam->SelMode == 1)	{
		return exp(*GWf * log(Stationary(j)) - (1-*GWf) * log(Stationary(i)));
	}
	double temp = mStationary[j]*mStationary0[i]/mStationary[i]/mStationary0[j];
	if (fabs(temp-1) < 1e-8)	{
		return 1;
	}
	return temp * log(temp) /(temp -1);
}

int SubMatrix::ComputeArrayHomoMutSel()	{

	double norm = 1;
	if (Normalise)	{
		norm = ComputeRate();
	}

	if ((mParam->MutMode == 2) || (mParam->MutMode == 3))	{
		ComputeCodonProjection();
	}
	if (mParam->SelMode == 1)	{
		ComputeGWStationary();
	}
	for (int i=0; i<Nstate; i++)	{
		for (int j=0; j<Nstate; j++)	{
			if (i==j)	{
				double temp = 0;
				for (int k=0; k<Nstate; k++)	{
					if (k != i)	{
						temp -= RelativeRate(i,k) * mStationary0[k] * SelectionFactor(i,k);
					}
				}
				operator()(i,j) = temp;
			}
			else	{
				operator()(i,j) = RelativeRate(i,j) * mStationary0[j] * SelectionFactor(i,j);
			}
		}
	}

	if (Normalise)	{
		for (int i=0; i<Nstate; i++)	{
			for (int j=0; j<Nstate; j++)	{
				operator()(i,j) /= norm;
			}
		}
	}
	return 0;
}

double SubMatrix::GetOnRate()	{

	if (NormaliseCov)	{
		return *rate / (1 - *invprob);
	}
	return *rate;
}

int SubMatrix::ComputeArrayHetero()	{
	
	double r = GetOnRate();

	for (int i=0; i<Nstate; i++)	{
		for (int j=0; j<Nstate; j++)	{
			operator()(i,j) = 0;
		}
	}

	int N = Nstate / 2;
	
	for (int i=0; i<N; i++)	{
		MatrixStationary[i] = *invprob * mStationary[i];
		MatrixStationary[N + i] = (1 - *invprob) * mStationary[i];
	}

	if (mRelativeRate)	{
		double norm = 1;
		if (Normalise)	{
			norm = ComputeRate();
		}

		for (int i=0; i<N; i++)	{
			for (int j=0; j<N; j++)	{
				if (i==j)	{
					double temp = 0;
					for (int k=0; k<N; k++)	{
						if (k != i)	{
							temp -= r * RelativeRate(i,k) * mStationary[k];
						}
					}
					operator()(i+N,j+N) = temp / norm;
				}
				else	{
					operator()(i+N,j+N) = r * RelativeRate(i,j) * mStationary[j] / norm;
				}
			}
		}
	}
	else	{
		for (int i=0; i<N; i++)	{
			for (int j=0; j<N; j++)	{
				operator()(i+N,j+N) = r * mStationary[j];
				if (i==j)	{
					operator()(i+N,j+N) -= r;
				}
			}
		}
	}	

	double s10 = *xi * *invprob;
	double s01 = *xi * (1 - *invprob);
	for (int i=0; i<N; i++)	{
		operator()(i,i) = -s01;
		operator()(i,i+N) = s01;
		operator()(i+N, i) = s10;
		operator()(i+N, i+N) -= s10 ;
	}

	return 0;
}

double SubMatrix::Power(int n, int i, int j)	{

	if (!n)	{
		return (i == j);
	}
	if (n > mParam->UniSubNmax)	{
		return Stationary(j);
	}
	if (n > npow)	{
		ComputePowers(n);
	}
	return mPow[n-1][i][j];
}

// ---------------------------------------------------------------------------
//		 ComputePowers()
// ---------------------------------------------------------------------------

void SubMatrix::ComputePowers()	{

	UniMu = 0;
	if (HeteroMode == Covarion)	{
		double r = GetOnRate();
		if (mRelativeRate)	{
			int N = Nstate / 2;
			for (int i=0; i<N; i++)	{
				double total = 0;
				for (int k=0; k<N; k++)	{
					if (k!= i)	{
						total += r * RelativeRate(i,k) * mStationary[k];
					}
				}
				if (UniMu < total)	{
					UniMu = total;
				}
			}
			UniMu += *xi;
		}
		else	{
			UniMu = r + *xi;
		}
	}
	else	{
		for (int i=0; i<Nstate; i++)	{
			if (UniMu < fabs(operator()(i,i)))	{
				UniMu = fabs(operator()(i,i));
			}
		}
	}

	CreatePowers(0);
	for (int i=0; i<Nstate; i++)	{
		mPow[0][i][i] = 1;
	}
	for (int i=0; i<Nstate; i++)	{
		for (int j=0; j<Nstate; j++)	{
			mPow[0][i][j] += operator()(i,j) / UniMu;
			if (mPow[0][i][j] < 0)	{
				cerr << "error in SubMatrix::ComputePowers: negative prob : ";
				cerr << i << '\t' << j << '\t' << mPow[0][i][j] << '\n';
				cerr << "Nstate : " << Nstate << '\n';
				exit(1);
			}
		}
	}
	npow = 1;
}

void SubMatrix::ComputePowers(int N)	{

	for (int n=npow; n<N; n++)	{
		CreatePowers(n);
		for (int i=0; i<Nstate; i++)	{
			for (int j=0; j<Nstate; j++)	{
				double& t = mPow[n][i][j];
				for (int k=0; k<Nstate; k++)	{
					t += mPow[n-1][i][k] * mPow[0][k][j];
				}
			}
		}
	}
	npow = N;
}

// ---------------------------------------------------------------------------
//		 ComputeRate()
// ---------------------------------------------------------------------------

double SubMatrix::ComputeRate()	{

	int N = Nstate;
	if (HeteroMode == Covarion)	{
		N /= 2;
	}

	double norm = 0;
	for (int i=0; i<N-1; i++)	{
		for (int j=i+1; j<N; j++)	{
			norm += MatrixStationary[i] * operator()(i,j);
		}
	}
	return 2 * norm;
}

// ---------------------------------------------------------------------------
//		 Diagonalise()
// ---------------------------------------------------------------------------

int SubMatrix::Diagonalise()	{

	bool failed = true;
	double * w = new double[Nstate];
	int* iw = new int[Nstate];

	if ((mParam->HeteroMode == Covarion) && (! mRelativeRate))	{

		double q = 1 - *invprob;
		double p = *invprob;
		double x = *xi;
		double r = GetOnRate();

		int N = Nstate / 2;
		double vald[2][2];
		double vecd[2][2][2];

		vald[0][0] = 0;
		vald[0][1] = -x;
		vecd[0][0][0] = -x*p;
		vecd[0][1][0] = -x*q;
		vecd[0][0][1] = x*q;
		vecd[0][1][1] = -x*q;

		double sd = sqrt((r-x)*(r-x) + 4*x*p*r);
		vald[1][0] = (-r-x+sd)/2;
		vald[1][1] = (-r-x-sd)/2;
		vecd[1][0][0] = (-r+x-sd)/2 -x*p;
		vecd[1][1][0] = -x*q;
		vecd[1][0][1] = (-r+x+sd)/2 -x*p;
		vecd[1][1][1] = -x*q;

		double valm[N];
		double vecm[N][N];
		for (int i=0; i<N; i++)	{
			valm[i] = -1;
			for (int j=0; j<N; j++)	{
				vecm[i][j] = 0;
			}
		}
		valm[0] = 0;
		for (int i=0; i<N; i++)	{
			vecm[i][0] = mStationary[i];
		}
		for (int i=1; i<N; i++)	{
			vecm[0][i] = -1;
			vecm[i][i] = 1;
		}

		for (int j=0; j<N; j++)	{
			int m = 0;
			if (valm[j] == -1)	{
				m = 1;
			}
			for (int i=0; i<2; i++)	{
				v[i*N+j] = vald[m][i];
				vi[i*N+j] = 0;
				for (int k=0; k<2; k++)	{
					for (int l=0; l<N; l++)	{
						u[k*N+l][i*N+j] = vecd[m][k][i] * vecm[l][j];
					}
				}	
			}
		}
				
		LinAlg::Gauss(u,Nstate,invu);
	}
	else if (!mRelativeRate)	{
		cerr << "compute array for poisson matrices\n";
		cerr << "deprecated\n";
		for (int i=0; i<Nstate; i++)	{
			v[i] = -1;
			vi[i] = 0;
		}
		v[0] = 0;
		for (int i=0; i<Nstate; i++)	{
			for (int j=0; j<Nstate; j++)	{
				u[i][j] = 0;
			}
		}
		for (int i=0; i<Nstate; i++)	{
			u[i][0] = mStationary[i];
		}
		for (int i=1; i<Nstate; i++)	{
			u[0][i] = -1;
			u[i][i] = 1;
		}
		// copy u into a
		for (int i=0; i<Nstate; i++)	{
			for (int j=0; j<Nstate; j++)	{
				a[j][i] = u[j][i];
			}
		}
		// invert a into invu
		// InvertMatrix(a, Nstate, w, iw, invu);
	}
	else	{	

		// transpose-copy Q into a :
		// (stupid convention chosen early on...)
		for (int i=0; i<Nstate; i++)	{
			for (int j=0; j<Nstate; j++)	{
				a[j][i] = Q[i][j];
			}
		}
			
		// diagonalise a into v and u
		int nmax = 1000;
		double epsilon = 1e-20;
		int n = LinAlg::DiagonalizeRateMatrix(a,MatrixStationary,Nstate,v,u,invu,nmax,epsilon);
		failed = (n == nmax);
	
		if (failed)	{
			cerr << "stat:\n";
			for (int i=0; i<Nstate; i++)	{
				cerr << MatrixStationary[i] << '\n';
			}
			cerr << "\n";
			for (int i = 0; i<Nstate; i++)	{
				for (int j=0; j<Nstate; j++)	{
					cerr << operator()(i,j) << '\t';
				}
				cerr << '\n';
			}
			cerr << '\n';
			exit(1);
		}

		// transpose eigenvector matrices
		// (stupid convention chosen early on, see above...)
		for (int i=0; i<Nstate; i++)	{
			for (int j=0; j<Nstate; j++)	{
				double tmp = u[i][j];
				u[i][j] = invu[j][i];
				invu[j][i] = tmp;
			}
		}

	}

	delete[] w;
	delete[] iw;
	return failed;
}


// ---------------------------------------------------------------------------
//		 CheckDiag()
// ---------------------------------------------------------------------------

void SubMatrix::CheckDiag()	{

	int width = 10;

	cout << "checking matrix diagonalisation\n\n";

	cout << "initial matrix\n\n";

	cout << "stationaries\n";
	for (int i=0; i<Nstate; i++)	{
		// cout << setw(width) << ((double) ((int) (Stationary(i)*1000))) / 1000 << '\n';
		cout << Stationary(i) << '\n';
	}
	cout << '\n';

	if (mRelativeRate)	{
		cout << "relative rates\n";
		for (int i=0; i<Nstate; i++)	{
			for (int j = i + 1; j<Nstate; j++)	{
				cout << setw(width) << ((double) ((int) (RelativeRate(i,j)*1000))) / 1000 ;
			}
			cout << '\n';
		}
	}

	cout << "initial matrix : \n";

	for (int i=0; i<Nstate; i++)	{
		for (int j=0; j<Nstate; j++)	{
			cout << setw(width) << ((double) ((int) (Q[j][i]*1000))) / 1000;
		}
		cout << '\n';
	}

	cout << "eigenvalues (real and imaginary part, the latter should be zero)\n";
	for (int i=0; i<Nstate; i++)	{
		cout << v[i] << '\t' << vi[i] << '\n';
		// cout 	<< setw(width) << ((double) ((int) (v[i]*1000))) / 1000
		// 		<< setw(width) << ((double) ((int) (vi[i]*1000))) / 1000 << '\n';
	}
	cout << '\n';

	cout << "eigenvector matrix\n";
	for (int i=0; i<Nstate; i++)	{
		for (int j = 0; j<Nstate; j++)	{
			cout << setw(width) << ((double) ((int) (u[j][i]*1000))) / 1000;
		}
		cout << '\n';
	}

	cout << "eigenvector matrix\n";
	for (int i=0; i<Nstate; i++)	{
		for (int j = 0; j<Nstate; j++)	{
			cout << setw(width) << ((double) ((int) (invu[j][i]*1000))) / 1000;
		}
		cout << '\n';
	}

	cout << '\n';
	
	cout << "check inverse\n";

	for (int i=0; i<Nstate; i++)	{
		for (int j = 0; j<Nstate; j++)	{
			double temp = 0;
			for (int k=0; k<Nstate; k++)	{
				temp += u[i][k] * invu[k][j];
			}
			// cout << setw(width) << ((double) ((int) (temp*1000))) / 1000;
			if (fabs(temp) < 1e-4)	{
				temp = 0;
			}
			cout << temp << '\t';
		}
		cout << '\n';
	}


	cout << "comput P-1 * Q * P \n\n";

	// now, let b = Q

	double* barray = new double[NstateNstate];
	double** b = new double*[Nstate];
	for (int i=0; i<Nstate; i++)	{
		b[i] = barray + i * Nstate;
		for (int j=0; j<Nstate; j++)	{
			b[i][j] = Q[i][j];
		}
	}

	// make a new matrix c
	// and let c = P-1 * b

	double* carray = new double[NstateNstate];
	double ** c = new double*[Nstate];

	for (int i=0; i<Nstate; i++)	{
		c[i] = carray + i * Nstate;
		for (int j=0; j<Nstate; j++)	{
			double temp = 0;
			for (int k=0; k<Nstate; k++)	{
				temp += invu[i][k] * b[k][j];
			}
			c[i][j] = temp;
		}
	}

	// comput c * P --> should be the diagonal matrix

	for (int i=0; i<Nstate; i++)	{
		for (int j=0; j<Nstate; j++)	{
			double temp = 0;
			for (int k = 0; k<Nstate; k++)	{
				temp += c[i][k] * u[k][j];
			}
			// cout << setw(width) << ((double) ((int) (temp*1000))) / 1000;
			if (fabs(temp) < 1e-4)	{
				temp = 0;
			}
			cout << temp << '\t';
		}
		cout << '\n';
	}

	cout << '\n';
}


