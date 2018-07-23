

class MCParameters;
class PhyloBayes;
class TaxaParameters;

class Sample	{


	public :

			Sample(string inSampleName);
			Sample(string inChainName, int burnin, int every, int until = -1, string Path = "");
			~Sample();

	MCParameters*	GetParameters()	{return mParam;}
	int		GetSize()	{return mSize;}
	int		GetChainNumber()	{return Nchain;}

	PhyloBayes*	GetNextPB();
	void		ToFile(string Name = "");

	void		ReadCATGTRPoint();
	void		MeanRR();
	void		Dating(int ps = 0, int verbose = 0, double alpha = 0.05);
	void		ReadCont(double meanlogT, int ps = 0);
	void		ReadCont2(int ps = 0);
	void		ReadPopEff(int ps = 0);
	void		ReadBranchFreqs();

	void		ReadSiteLogLikelihood();
	void 		ReadSummedLogLikelihood();
	void		ReadFiniteTimeEntropy(double timemax, int N);
	// private :

	void		OpenChainStream();
	void		Reset();

	void		ComputeStatePostProb(double*** ret = 0, int singlefreq = 0);

	void		ReadRoot(int verbose);
	void		Coaff(int type, int step = 10);

	void		GrepModes();
	void		TrueLength(int nrep, int sitebranch);

	void		Check(PhyloBayes* truepb, int verbose, int topo, int ncat, int* truehisto, int* tothisto, double& alpha, double& gamma, double& length, double& statent, double& rrent);

	int		ReadBench();
	int		Read(int rates = 0, int modes = 0, int sitestat = 0, double cutoff = 0.5, int ncat = 20, int ps = 0);
	int 		ClusterModes(int SizeThreshold=1, double ClusteringRange=0.01, int ps = 0);
	void		ReadModes();
	void		ReadNHCov();
	void		ReadNHVar();
	void		ReadContCharacters();

	void		SiteCount(int nrep = 1);
	void		SampleSub(int nrep, int prred, int priormode, int priorratemode, int priorrootstate, int* sitemask);
	void		SampleSub(double*** statepostprob = 0, int* sitemask = 0);
	void		SampleSubWS(int priormode, int priorratemode, int priorrootstate);

	void		ClockRate(int nrep = 1);
	void		ClockRateAutoCorrel(int nrep = 1);
	void		PosteriorPredictive(int nrep = 1, int priormode=1, int priorratemode=1, int priorrootstate=1);
	void 		Homoplasy(int nrep, int priormode, int priorratemode, int priorrootstate);
	void 		HomoplasyStochastic(int nrep, int ncat, int priormode, int priorratemode, int priorrootstate, int redraw, double cutoff);
	void 		GenePostPred(string genepartition, int nrep, int ncat, int priormode, int priorratemode, int priorrootstate, int redraw);
	void 		postpredLRT(int nrep, int ncat, int priormode, int priorratemode, int priorrootstate);
	void 		DiscrepancyTest(int nrep, int ncat, int priormode, int priorratemode, int priorrootstate);
	void 		postpredDoubleCount(int nrep, int ncat, int priormode, int priorratemode, int priorrootstate);
	void 		Diversity(int nrep, int ncat, int priormode=1, int priorratemode=1, int priorrootstate=1);
	void 		DifferentialDiversity(string taxonmask, int nrep, int ncat, int priormode=1, int priorratemode=1, int priorrootstate=1);
	void		HomogeneityTest(int nrep = 1, double threshold = -1, int priormode = 1, int priorratemode = 1, int priorrootstate = 1, int skew = 0);
	double		HomogeneityMeasure(int** data, double* dev, double& mean);
	double		CodonHomogeneityMeasure(int** data, double* dev, double& mean);
	double		SerineSkew(int** data, double* dev);
	double		LeucineSkew(int** data, double* dev);
	double		ArginineSkew(int** data, double* dev);

	double		CV(string testfile);
	void		Cov();
	void		CovWoOutgroup();

	void		ConstantSiteCorrection(double rho);

	int		mSize;
	MCParameters*	mParam;
	int		Nchain;

	int		Burnin;
	int		Every;

	int		extract;
	ifstream*	Sample_is;
	ifstream*	Chain_is;

	string 		SampleName;
	string		ChainName;

	PhyloBayes*	currentPB;
	int		skip;

	int		currentChain;
	int		currentIndex;

	int 		eof;


}
;

