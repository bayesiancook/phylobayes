
class MCParameters;

class SubMatrix		{

	public:

			// this will create a standard time-reversible Nstate x Nstate substitution matrix
			// whose instant rates are defined by the set of relative rates inRR (dim = Nstate(Nstate-1)/2)
			// and the set of stationary probabilities inStationary (dim = Nstate, sum to 1)
			SubMatrix(MCParameters* inParam, int inNstate, double* inRR, double* inStationary);
			SubMatrix(MCParameters* inParam, int inNstate, double* inRR, double* inStationary0, double* inStationary, double *f, double a);

			// this will create a covarion matrix, based on
			// a substitution matrix defined by Nstate, inRR and inStationary
			// and a covarion process of parameters xi0 and invprob
			// you also need to specify the rate of the substitution process when in the on state
			// through the inrate pointer
			// if inParam->NormaliseCov,
			// then, the rate of the substitution process will actually be (*inrate) / (1 - *invprob)
			SubMatrix(MCParameters* inParam, int inNstate, double* inrate, double* ininvprob, double* inxi, double* inRR, double* inStationary);
			
			// same thing as above, except that no relative rates are specified
			// in that case, they are all equal to 1
			SubMatrix(MCParameters* inParam, int inNstate, double* inrate, double* ininvprob, double* inxi, double* inStationary);
			~SubMatrix();

	void		BasicAllocation();

	void		Copy(SubMatrix* from);
	void		CopyStationaries(double* from);

	int 		ComputeArray();
	int 		ComputeArrayHomo();
	int 		ComputeArrayHomoMutSel();
	int 		ComputeArrayHetero();

	void		ComputeGWStationary();
	void		ComputeCodonProjection();

	int 		Diagonalise();
	double		ComputeRate();

	void		CreatePowers();
	void		CreatePowers(int n);
	void		DeletePowers();
	int		npow;

	void		ComputePowers();
	void		ComputePowers(int n);
	double		Power(int n, int i, int j);

	void		ComputeExponential(double range, double** expo);
	double**	GetGenerator()	{return Q;}

	double& 	Stationary(int i);
	double& 	RelativeRate(int i, int j);

	double* 	GetStationaries()	{ return MatrixStationary;}

	double 	EigenValue(int i) { return v[i];}

	double*	GetRelRatePtr() { return mRelativeRate;}
	void	ScalarMul(double e);

	// access to matrix Q colums and rows
	// inline double& operator()(int i, int j);			//
	// inline double& operator()(int i, int j) const ;
	double& operator()(int i, int j);			//
	double& operator()(int i, int j) const ;


	void	CheckDiag();

	double SelectionFactor(int i, int j);
	//-------------
	// data fields
	//-------------

	double* mArray;

	// matrices as lists of rows
	// Q : the infinitesimal generator matrix
	// a : the copy of it -> to feed linalg ( as a list of rows)
	// v : eigenvalues
	// vi : imaginary part
	// u : the matrix of eigen vectors
	// invu : the inverse of u
	// wi : work vector
	// w : work vector
	//
	// total : 4 Nstate*Nstate + 4 * Nstate

	double ** Q;
	double ** a;
	double ** u;
	double ** invu;
	double * v;
	double * vi;

	double* MatrixStationary;
	double*	MatrixRelativeRate;		

	double* mStationary;
	double* mStationary0;
	double*	mRelativeRate;	

	int RelRateAlloc;
	int StatAlloc;
	int Stat0Alloc;
	int MatrixStatAlloc;

	int Nstate;
	int NstateNstate;
	MCParameters*	mParam;

	double UniMu;
	double*** mPow;

	Switch Normalise;
	Switch NormaliseCov;
	double GetOnRate();

	HeterotachyMode HeteroMode;
	double* xi;
	double* invprob;
	double* rate;
	double* GWf;
	double* kappa;
	double* acgt;
	double* mRR_ACGT;
	double RR_ACGT(int i, int j);
	int IsTransition(int i, int j);
	int position(int a, int b);

};
