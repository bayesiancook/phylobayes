
class PhyloBayes;
class Tree;

class MCParameters {

	public:

				MCParameters();
				~MCParameters();

	PhyloBayes*		GetCurrentState()	{return currentState[0];}
	PhyloBayes*		GetNextState()		{return nextState[0];}

	PhyloBayes*		GetCurrentState(int i)	{return currentState[i];}
	PhyloBayes*		GetNextState(int i)		{return nextState[i];}

	void			SetCurrentState(int i, PhyloBayes* inPB)	{ currentState[i] = inPB;}
	void			SetNextState(int i, PhyloBayes* inPB)	{ nextState[i] = inPB;}

	void			SetInitTree(Tree& tree);

	void			Reset();				// called before the first start
	void			Update(); 				// called before any start

	void			InitFromFile(string filename);
	void			ReadMoveType(istream& is);
	void			ReadSuperMoveType(istream& is);
	void			ReadTopoMoveType(istream& is);
	void			ReadCalibration(istream& is);
	void			ReadCalibration(string filename);
	void			ReadPartition(istream& is);
	void			ReadPartition(string filename);
	int			ReadContDataFromFile(string datafilename);
	int			ReadDataFromFile(string datafilename, int forceinterleaved = 0);
	int			ReadNexus();
	int			ReadSpecial();
	int			TestPhylipSequential();
	int			TestPhylip(int repeattaxa);
	void			ReadPhylipSequential();
	void			ReadPhylip(int repeattaxa);
	void			MissingData();
	double			PropInvariant();
	double			MeanDiversity(int* histo = 0, double* tmp = 0);
	double			MeanDifferentialDiversity(int* taxonmask);
	double			SiteDiversity(int* site);
	void			WriteDataToFile(string datafilename, int N1=0, int N2=0);
	void			WriteDataToFile(string datafilename, int * Mask);
	void			WriteNexus(string datafilename);
	void			WriteDataToFileSiteMask(string datafilename, int * Mask);
	void			SwapData(MCParameters* from);
	void			RegisterWithData();
	void			ComputeZipArrays();
	void			EstimateEmpiricalFrequencies();
	void			EliminateConstantPositions();
	void			EliminateUnknownColumns();
	void			GetSequenceProfile(ostream& os, double pseudocount);
	int			ConstantColumns();

	void 			InitMove();
	void 			InitMoveBF();
	void			PushMove(const MoveType movetype, const int Nrep, const double delta, const int N);
	void 			DesactivateLengthMoves();

	void			Move();
	void			Move(int chain);
	void			TopoMove(int chain, int topo);
	void			SuperMove(int chain, int topo, int super);
	void			Move(int chain, int topo, int super, int move);
	void			Swap(int chain);

	friend ostream& operator<<( ostream& os , MCParameters& param);
	friend istream& operator>>( istream& is , MCParameters& param);

	int			Version;
	int			SubVersion;

	// -----------------------------
	// Model Parameters
	// ------------------------

	// data matrix

	string			DataFileSpec;				// nexus file containing the aligned sequences
	string			ContDataFileSpec;				// nexus file containing the aligned sequences
	int			ContNchar;

	Chrono			ch1;
	Chrono			ch2;

	int 			Ntaxa;
	int 			Nnode;
	int			Ngene;
	int*			Gene;
	int 			Nsite;
	int*			GeneSize;
	int*			GeneFirstSite;
	int**			GeneTaxa;
	int			mGetGeneTaxaFromData;
	void			GetGeneTaxaFromData();

	string*			SpeciesNames;

	int**	 		Data;						// Data[i][j] : taxon i, site j
	double**	 	ContData;						// Data[i][j] : taxon i, site j
	int**		 	ContMissingData;						// Data[i][j] : taxon i, site j
	int 			Nstate;
	char*			Alphabet;
	char*			AlphabetSet;
	int			NAlphabetSet;

	int**	 		RawData;						// Data[i][j] : taxon i, site j
	int 			RawNstate;
	const char*		RawAlphabet;

	RecodingMode		Recoding; // 0: no recoding 1 : dayhoff5 2: Dayhoff4: 2: hp 3: custom
	string			RecodingFile;
	void			RecodeData();
	void			LoadRecoding();
	int*			RecodingTable;
	char*			RecodingAlphabet;
	int			RecodingNstate;
	
	Boolean**		Orbit;
	int*			OrbitSize;
	double			MeanOrbitSize;
	int*			ZipSize;
	int**			Indices;
	int**			ZipIndices;
	int**			ZipData;
	int*			ConstantState;

	double*			EmpiricalFreq;
	int**			SiteEmpiricalCount;

	Switch			DeleteConstant;
	Switch			Normalise;					// 1 : self-substitutions excluded while measuring lengths
	Switch			NormaliseCov;	
										// 0 : self-substitutions included
	Switch			SavePartialLogLikelihoods;
	Switch			ModeFastCompute;
	Switch			RefFastCompute;
	Switch			RefPoisson;
	Switch			ModePoisson;
	Prior			LengthPrior;					// 0 : flat
	double			LengthMin;
	double			LengthMax;

	Prior			RatePrior;					// 0 : flat (only possible value as for now)
	Prior			GammaPrior;
	double 			RateMin;
	double			RateMax;
	double			GammaMin;
	double			GammaMax;
	double			XiMin;
	double			XiMax;

	double 			RRAlphaMin;
	double 			RRAlphaMax;
	
	Prior			ModeStatPrior;					// idem
	Prior			ModeStatAlphaPrior;
	double			StatMin;
	double			StatMax;
	double			StatAlphaMin;
	double			StatAlphaMax;

	Prior			AlphaPrior;					// 0 : flat
	double			AlphaMin;
	double			AlphaMax;
	double			BetaMin;
	double			BetaMax;

	Prior			RateAlphaPrior;					// 0 : flat
	double			RateAlphaMin;
	double			RateAlphaMax;

	double			TooSmall;
	double			TooLarge;
	double			InfProb;

	// ----------------------
	// chain parameters
	// ----------------------

	// Moves

	int 			TopoMoveTypeNumber;
	int**			TopoNIterationArray;
	int**			TopoNArray;
	TopoMoveType*		TopoMoveTypeArray;

	int* 			SuperMoveTypeNumber;
	int***			SuperNIterationArray;
										//  	Global, Local, Rates, SwitchMode...
	int** 			MoveTypeNumber;
	MoveType***		MoveTypeArray;					// type of nth. move (MoveType can take one of the values :
										//  	Global, Local, Rates, SwitchMode...
	int****			NIterationArray;				// mean number of calls to nth. move per cycle

	double****		deltaArray;					// real tuning parameter for the nth. move
										// 		(the meaning of this parameter is MoveType dependent)
	int****			NArray;						// integer tuning parameter for the nth. move
										//		(idem)
	// Saves and Stops
	int 			SaveEvery;					// Number of calls defining one chain's cycle
										// after each cycle, the chain saves its current state
										// and appends this state to the sample
	int			StopAfter;					// the chain will stop after having recorded
										// StopAfter points (burn in included)
										// if StopAfter == -1 	: the chain never stops
	int			HowManySaved;
	double			SwapFreq;

	// thermodynamic integration :

	ModelSwitchMode		MSMode;

	int			BurnIn;
	int 			TopoBurnin;

	double			RASModelSwitch;
	double			SUBModelSwitch;
	double			HeteroModelSwitch;
	double			SeparateModelSwitch;

	
	double			InitBeta;
	double			FinalBeta;
	double			BetaStep;
	int			BurnInDone;

	double*			Beta;
	double			Beta0;
	int			Nchain;


	// Monte Carlo statistics

	// the following counters keep track of

	int			RateInfCount;					// # times maxrate/minrate exceeded authorised upper limit
	int			StatInfCount;					// # times maxstat/minstat exceeded authorised upper limit
	int			LengthInfCount;					// # times max branch length overflowed
	int			LogProbInfCount;				// # times a probability was found too small for its log to be computed
	int 			SubOverflowCount;	// # same as StatInfCount, but specific for switch mode move

	int****			NCallArray;
	double****		TimeArray;
	double****		SuccessArray;

	int*			TrialSwapArray;
	int*			AcceptedSwapArray;

	PhyloBayes**		currentState;					// the current point in hypothesis space
	PhyloBayes**		nextState;						// its back up copy

	// ----------------------------------------
	// accesssory parameters
	// used for compaction of the state representation (in the case of catfish)
	// ------------------------------------------------



	int			DataOK;
	double 			SpeedFactor;

	int*			SwapPermut;

	int			Parallel;
	int			SwapMode;
	int			StorePartialLikelihoods;
	int			SaveAllChains;

	HeterotachyMode		HeteroMode;

	ClockModelType		ClockModel;
	FlexClockModelType	FlexClockModel;

	TimePriorType		TimePrior;

	int			ImproperLowerBound;
	double			LowerP;
	double			LowerC;

	Prior			ModePrior;
	Prior			RateModePrior;

	string			Path;
	void			SetPath(string inPath){Path = inPath;}
	string			GetPath(string inPath){return Path;}

	Switch			Qmode;
	Switch			GeneGammaMode;
	Switch			GeneRateMode;
	Switch			GeneBLMode;
	Switch			GeneStationaryMode;
	Switch			GeneBLMultiplier;
	
	int			GammaNcat; // if 0 : continuous gamma

	Switch			ActivateSumOverModes;
	Switch			ActivateSumOverRateModes;
	Switch			SumOverRateModes;
	Switch			SumOverModes;
	Switch			SimpleSampling;
	
	int			NmodeMax;
	int			NRateModeMax;

	int 			Rao;
	int			StatNrep;
	double			StatEpsilon;
	int			StatExpectNrep;
	int			StatNcycle;
	int 			SaveStat;
	int 			SaveAll;

	int			FixGamma;
	int			FixMeanLength;
	int			FixStatCenter;
	int			FixStat;
	int			FixRR;
	int 			FixLength;
	int			FixRate;
	int			FixNmode;
	int			FixTopo;
	int			FixPconst;
	void 			Init(RRMode rr, int empfreq, int ncat, int statcenter, RASMode ras, int discrate, int fixmeanlength, double meanlength = 0.1);
	void			InitClock(int clock, int genemode);

	void			MoveAll(int conjugate = 1);
	void			MoveStat(int conjugate = 1);
	void			MoveRate(int conjugate = 1);
	void			MoveSiteVariables(int conjugate = 1);


	void			ActivateDiscreteGamma(int discrate);

	int 			ReducedCounts;

	Switch			NormalApprox;
	int			NormalConcMode; // 0 : not defined ; 1: Use one matrix ; 2 : use gene specific matrices

	Switch			ActivateClock;
	double			FlexClockModelSwitch;
	double			BLModelSwitch;
	double			ClockPriorModelSwitch;
	double			AutoCorrelModelSwitch;

	double			Zeta;
	double			ArcTh(double zeta);
	double			ArcdH(double zeta);
	double			ArcH(double zeta);
	double			ArcW(double zeta);
	double			ArcZ();
	int			ArcThermo;
	double 			alpha1;
	double 			alpha2;
	

	Tree*			TemplateTree;

	double**		InvCov;
	double			LogDetCov;
	double			LogLMax;
	double*			MaxBL;
	Tree*			ConcTree;
	int*			MapBL;
	int			Nbranch;

	double***		GeneInvCov;
	double*			GeneLogDetCov;
	double*			GeneLogLMax;
	double**		GeneMaxBL;
	Tree**			GeneTree;
	int**			GeneMapBL;

	int*			GeneNtaxa;
	int*			GeneNbranch;
	string**		GeneSpeciesNames;
	string*			GeneNames;

	int			RootDegree;
	// return root degree
	int			BaseReadCov(string datafile, string*& speciesnames, Tree*& tree, double*& maxbl, double**& incov, double& logdetcov, double& loglmax, int& Ntaxa, int& Nbranch);

	void			ReadCov(string datafile);
	void			ReadSepCov(string datafile);

	string			ConcCovFileName;
	string			SepCovFileName;
	int			ConcCov;
	int			SepCov;
	void			ReadOutGroup(string filespec);
	int			ReadCatConstraints(string filename);
	int**			StatFlag;

	double			RhoMin;
	double			RhoMax;
	double			SigmaMin;
	double			SigmaMax;
	double			ThetaMin;
	double			ThetaMax;
	double			BesselMax;
	double			Mumax;

	Chrono			movetime;

	int			NCalib;
	int*			isCalibrated;
	double*			CalibLower;
	double*			CalibUpper;
	string*			CalibTaxon1;
	string*			CalibTaxon2;
	int*			CalibIndex;

	int			SoftBounds;
	double			Softa;

	Switch			OptimizeNormalLogSampling;
	Switch			SeparateRhoPrior;

	string*			OutGroup;
	int			OutNtaxa;

	void			RegisterCov();
	void			RegisterCovDichotomous();
	void			RegisterCovTrichotomous();
	
	int			ZipSub;
	int			Ncat;
	double**		statfix;
	double**		statrr;
	double*			weight;
	void			ReadMixture(string filename);
	void			ReadMatFix(string filename);
	void			ReadStatFix(string filename);
	void			WriteMixture(ostream& os);

	int			Nrho;
	Switch			DCM;

	Switch			UniSub;
	Switch			UniSubOnly;
	int			UniSubNmax;
	int			ObservedUniSubNmax;
	double			logBF;
	double			logBFmin;
	double			logBFmax;
	double			logBFrev;
	double			logBFminrev;
	double			logBFmaxrev;
	int			BFcount;
	Prior			ScalePrior;
	double			MeanScale;	
	double			VarScale;	

	int			LikelihoodOffset;
	int			FixedWeights;
	int			FixedLengths;
	int			FixedAlpha;
	int			ActivateEM;
	int			RandomInitLength;

	int 			WithPartial;

	Switch			ExternalCov;

	int			Conjugate;

	int			NH; // 0 : homogeneous 1: full NH 2: BNH 3 : FNH
	int			NHNcatMax;
	int			NHPrior; // 0: dirichlet  1 : diagonal gaussian  2 : multivariate gaussian

	int			EmpiricalDP;
	int			IncrementalDP;
	int			IncrementalRateDP;
	int			Ngen;
	int 			NgenPerCycle();
	int 			ZipGTR;
	int 			ZipGTRDP;
	int 			ZipPrior;
	int 			NZipModeMax;

	int			Tempering;
	double			PriorOcc1;
	double			PriorOcc0;

	int			MBL;
	int			NMixBLMax;
	int			IncrementalMixBLDP;
	int 			MBLPrior;
	
	int 			LengthGammaPrior;

	double			FixRateAlpha;

	void			CompositionalCovariation();
	int			NHClampFirst;
	int			NHWishart;
	int			NHHalfSum;
	int			NHAutoRegressive;
	int			NHLengthSwitch;
	double*			ContVar;
	void			ComputeContVar();
	int 			KeepTimes;
	Prior			RRPrior;
	int			SelMode;
	int			MutMode;
	double			GWf;
	int			StartDP;
	double 			meannpow;

	int 			PopEff;
	int 			NHPop;

	double*			CustomRR;
	void			ReadCustomRR(string filename);

	GeneticCodeType		GeneticCode;
	int*			CodonCode;
	double 			Chi;
	double			Chi2;

	double			AlphaChi;
	double			BetaChi;
	double			AlphaChi2;
	double			BetaChi2;

	Switch			ResampleStateAtRoot;
	Switch			RedundantPath;
	
	int			FixRoot;
	int 			midpointrooting;
}
;
