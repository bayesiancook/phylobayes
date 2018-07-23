
#include "phylo.h"
#include "linalg.h"


//---------------------------------------------------------------------------------
//		 MCParameters()
//---------------------------------------------------------------------------------

MCParameters::MCParameters()	{

	Version = 4;
	SubVersion = 1;

	// only for ancestral sequences
	// under homogeneous models
	midpointrooting = 0;

	FixRateAlpha = 0;
	NormalConcMode = 0;
	RandomInitLength = 1;	
	OutGroup = 0;

	StatFlag = 0;
	FixRoot = 0;

	OutNtaxa = 0;
	FixTopo = 0;

	ResampleStateAtRoot = Yes;

	RedundantPath = Yes;

	Chi = -1;
	Chi2 = -1;

	AlphaChi = 1;
	BetaChi = 1;
	AlphaChi2 = 1;
	BetaChi2 = 1;

	ImproperLowerBound = 1;
	LowerC = 1;
	LowerP = 0.1;
	
	CustomRR = 0;
	GeneticCode = Universal;
	CodonCode = 0;

	RRPrior = Exponential;
	// RRPrior = GammaDistributed;

	MutMode = 0;
	// 0 : regular GTR
	// >0 : mutation selection
	// 1 : GTR + selection
	// 2 : Codon HKY
	// 3 : Codon GTR Nuc

	SelMode = 0;
	// if MutMode != 0
	// 0 : HB Halpern Bruno
	// 1 : GW Goldman Whelan

	GWf = 1;
	PopEff = 0;
	// 1 : root constrained to equal 1
	// 2 : iid gamma including root
	// 3 : with exchangeability move 

	NHPop = 0;
	// 1 : both NH and pop

	KeepTimes = 0;
	MBL = 0;
	NMixBLMax = 0;
	MBLPrior = 1;

	LengthGammaPrior = 1;
	
	NHClampFirst = 1;
	NHWishart = 1;
	NHHalfSum = 1;
	NHAutoRegressive = 1;
	NHLengthSwitch = 1;

	PriorOcc1 = 1;
	PriorOcc0 = 1;
	Tempering = 0;

	ZipGTR = 0;
	// = 1 : modes can have 0 probability for certain amino-acids
	// = 2 : each site has a 0/1 mask
	// = 3 : minimal mask
	// = 4 : projection, with unknown state

	ZipGTRDP = 1;
	// = 1 : masks are drawn from a dirichlet process
	
	ZipPrior = 0;
	// = 0 : uniform prior on the number of 1 in mask
	// = 1 : geometric
	// = 2 : mixture of independent binomial
	
	EmpiricalDP = 0;
	IncrementalDP = 1;
	IncrementalRateDP = 1;
	StartDP = 0;
	IncrementalMixBLDP = 1;
	ch1.Reset();
	ch2.Reset();
	
	NH = 0;
	// 1 : fixed number
	// 2 : bnh : branch specific modulators
	// 3 : fnh : foster (mixture model of branch modulators)
	// 4 : anh : autocorrelated
	// 5 : cnh : ?
	// 6 : dnh : ?

	NHNcatMax = 20;
	NHPrior = 0;

	Conjugate = 1;
	WithPartial = 0;
	DataFileSpec = "None";
	ContDataFileSpec = "None";
	ContData = 0;
	ContMissingData = 0;

	ActivateEM = 0;
	LikelihoodOffset = 1;
	FixedWeights = 0;
	FixedLengths = 0;
	FixedAlpha = 0;

	MeanScale = 1;
	VarScale = 1;
	ScalePrior = Flat;

	Zeta = 0;
	ArcThermo = 0;
	alpha1 = 1000;
	alpha2 = 1000;

	DCM = Yes;
	Ncat = 0;
	statfix = 0;

	UniSubOnly = No;
	UniSub = Yes;
	UniSubNmax = 5000;
	ObservedUniSubNmax = 0;

	RawNstate = 0;
	RawData = 0;
	RawAlphabet = 0;
	Recoding = NoRecoding;
	RecodingFile = "";
	RecodingTable = 0;
	RecodingAlphabet = 0;
	Alphabet = 0;
	
	ZipSub = 1;
	SeparateRhoPrior = No;
	OptimizeNormalLogSampling = Yes;
	Ngene = 0;
	mGetGeneTaxaFromData = 1;
	GeneTaxa = 0;
	NCalib = 0;
	CalibTaxon1 = 0;
	CalibTaxon2 = 0;
	CalibIndex = 0;
	CalibUpper = 0;
	CalibLower = 0;
	isCalibrated = 0;

	SoftBounds = 0;
	Softa = 0.025;

	ReducedCounts = 1;
	Path = "";

	movetime.Reset();

	ActivateClock = No;
	NormalApprox = No;
	ConcCov = 0;
	SepCov = 0;
	
	FlexClockModelSwitch = 1;
	BLModelSwitch = 1;
	ClockModel = CIRrigid;
	FlexClockModel = WhiteNoise;
	ClockPriorModelSwitch = 1;
	AutoCorrelModelSwitch = 1;
	Nrho = 0;
	
	TimePrior = Unif;

	ConcTree = 0;
	MaxBL = 0;
	InvCov = 0;
	MapBL = 0;
	
	RootDegree = 0;

	GeneTree = 0;
	GeneInvCov = 0;
	GeneMaxBL = 0;
	GeneMapBL = 0;
	GeneNbranch = 0;
	GeneNtaxa = 0;	
	GeneSpeciesNames = 0;
	GeneFirstSite = 0;
	GeneSize = 0;
	Gene = 0;

	SaveAll = 1;
	SaveStat = 1;
	GammaNcat = 0;

	StatNrep = 200;
	StatExpectNrep = 1000;
	StatNcycle = 10;
	StatEpsilon = 5000;
	Rao = 0;

	Qmode = No;
	GeneBLMode = No;
	GeneBLMultiplier = No;
	GeneGammaMode = No;
	GeneStationaryMode = No;
	GeneRateMode = No;
	NmodeMax = -1;
	NRateModeMax = 32;
	
	ActivateSumOverRateModes = No;
	SumOverRateModes = No;
	ActivateSumOverModes = No;
	SumOverModes = No;

	SimpleSampling = No;
	SwapPermut = 0;

	Parallel = 1;
	SwapMode = 1;
	SaveAllChains = 1;

	// phylogeny
	LengthMin = DefaultLengthMin;
	LengthMax = DefaultLengthMax;
	// LengthPrior = Exponential;
	LengthPrior = GammaDistributed;

	// rates
	RatePrior = DefaultRatePrior;
	GammaPrior = Cauchy;
	RateMin = DefaultRateMin;
	RateMax = DefaultRateMax;

	GammaMin = DefaultGammaMin;
	GammaMax = DefaultGammaMax;
	RRAlphaMin = 0.1;
	RRAlphaMax = 100;
	XiMin = DefaultXiMin;
	XiMax = DefaultXiMax;

	RhoMin = 0;
	RhoMax = 100;
	SigmaMin = 0;
	SigmaMax = 100;
	ThetaMin = 0.05;
	ThetaMax = 200;
	BesselMax = 500;
	Mumax = 500;

	ModePrior = Flat;
	RateModePrior = DirichletProcess;
	// substitution model
	ModeStatPrior = DefaultModeStatPrior;
	ModeStatAlphaPrior = Exponential;
	// StatMin = 0;
	StatMin = DefaultStatMin;
	StatMax = DefaultStatMax;
	StatAlphaMin = DefaultStatAlphaMin;
	StatAlphaMax = DefaultStatAlphaMax;

	// Dirichlet process
	AlphaPrior = DefaultAlphaPrior;
	// AlphaMin = DefaultAlphaMin;
	AlphaMin = 0;
	// AlphaMax = DefaultAlphaMax;
	AlphaMax = 100;

	// general
	TooSmall = DefaultTooSmall;
	TooLarge = DefaultTooLarge;
	InfProb = DefaultInfProb;
	DeleteConstant = DefaultDeleteConstant;
	// Normalise = DefaultNormalise;
	Normalise = No;
	SavePartialLogLikelihoods = No;
	
	ModeFastCompute = Yes;
	RefFastCompute = No;
	ModePoisson = Yes;
	RefPoisson = No;

	EmpiricalFreq = 0;

	SpeciesNames = 0;
	Data = 0;
	ContData = 0;
	ContNchar = 0;
	Orbit = 0;
	OrbitSize = 0;
	ZipSize = 0;
	Indices = 0;
	ZipData = 0;

	DataOK = 0;

	TopoMoveTypeNumber = 0;
	TopoNIterationArray = new int*[MaxTopoMoveTypeNumber];
	TopoMoveTypeArray = new TopoMoveType[MaxTopoMoveTypeNumber];
	TopoNArray = new int*[MaxTopoMoveTypeNumber];
	SuperMoveTypeNumber = new int[MaxTopoMoveTypeNumber];
	for (int i=0; i< MaxTopoMoveTypeNumber; i++)	{
		SuperMoveTypeNumber[i] = 0;
		TopoNIterationArray[i] = new int[MaxChainNumber];
		TopoNArray[i] = new int[MaxChainNumber];
		for (int j=0; j<MaxChainNumber; j++)	{
			TopoNIterationArray[i][j] = 0;
			TopoNArray[i][j] = 0;
		}
	}
	SuperNIterationArray = new int**[MaxTopoMoveTypeNumber];
	MoveTypeNumber = new int*[MaxTopoMoveTypeNumber];
	for (int i=0; i< MaxTopoMoveTypeNumber; i++)	{
		MoveTypeNumber[i] = new int[MaxSuperMoveTypeNumber];
		for (int j=0; j < MaxSuperMoveTypeNumber; j++)	{
			MoveTypeNumber[i][j] = 0;
		}
		SuperNIterationArray[i] = new int*[MaxSuperMoveTypeNumber];
		for (int j=0; j<MaxSuperMoveTypeNumber; j++)	{
			SuperNIterationArray[i][j] = new int[MaxChainNumber];
		}
	}

	MoveTypeArray = new MoveType**[MaxTopoMoveTypeNumber];
	NIterationArray = new int***[MaxTopoMoveTypeNumber];
	deltaArray = new double***[MaxTopoMoveTypeNumber];
	NArray = new int***[MaxTopoMoveTypeNumber];
	SuccessArray = new double***[MaxTopoMoveTypeNumber];
	TimeArray = new double***[MaxTopoMoveTypeNumber];
	NCallArray = new int***[MaxTopoMoveTypeNumber];
	for (int i=0; i < MaxTopoMoveTypeNumber ; i++)	{
		MoveTypeArray[i] = new MoveType*[MaxSuperMoveTypeNumber];
		NIterationArray[i] = new int**[MaxSuperMoveTypeNumber];
		deltaArray[i] = new double**[MaxSuperMoveTypeNumber];
		NArray[i] = new int**[MaxSuperMoveTypeNumber];
		SuccessArray[i] = new double**[MaxSuperMoveTypeNumber];
		TimeArray[i] = new double**[MaxSuperMoveTypeNumber];
		NCallArray[i] = new int**[MaxSuperMoveTypeNumber];
		for (int j=0; j<MaxSuperMoveTypeNumber; j++)	{
			MoveTypeArray[i][j] = new MoveType[MaxMoveTypeNumber];
			NIterationArray[i][j] = new int*[MaxMoveTypeNumber];
			NArray[i][j] = new int*[MaxMoveTypeNumber];
			deltaArray[i][j] = new double*[MaxMoveTypeNumber];
			SuccessArray[i][j] = new double*[MaxMoveTypeNumber];
			TimeArray[i][j] = new double*[MaxMoveTypeNumber];
			NCallArray[i][j] = new int*[MaxMoveTypeNumber];
			for (int k=0; k<MaxMoveTypeNumber; k++)	{
				NIterationArray[i][j][k] = new int[MaxChainNumber];
				NArray[i][j][k] = new int[MaxChainNumber];
				deltaArray[i][j][k] = new double[MaxChainNumber];
				SuccessArray[i][j][k] = new double[MaxChainNumber];
				TimeArray[i][j][k] = new double[MaxChainNumber];
				NCallArray[i][j][k] = new int[MaxChainNumber];
			}
		}
	}

	TrialSwapArray = new int[MaxChainNumber];
	AcceptedSwapArray = new int[MaxChainNumber];

	StopAfter = -1;
	BurnIn = 0;
	InitBeta = 0;
	FinalBeta = 1;
	BetaStep = 0;
	BurnInDone = 0;

	SaveEvery = 100;
	SwapFreq = 0;

	RASModelSwitch = 1;
	SUBModelSwitch = 1;
	HeteroModelSwitch = 0;
	SeparateModelSwitch = 0;
	MSMode = DefaultModelSwitchMode;
	
	HeteroMode = Homo;
	ExternalCov = Yes;
	NormaliseCov = Yes;

	Beta = new double[MaxChainNumber];
	Beta0 = 1;
	for (int l=0; l<MaxChainNumber; l++)	{
		Beta[l] = 1;
	}
	Nchain = 1;

	currentState = 0;
	nextState = 0;
	SiteEmpiricalCount = 0;
}


//---------------------------------------------------------------------------------
//		 ~MCParameters()
//---------------------------------------------------------------------------------


MCParameters::~MCParameters()	{

	if (currentState)	{
		for (int chain=0; chain<Nchain; chain++)	{
			delete currentState[chain];
			delete nextState[chain];
		}
		delete[] currentState;
		delete[] nextState;
	}
		
	delete[] EmpiricalFreq;
	if (SiteEmpiricalCount)	{
		for (int i=0; i<Nsite; i++)	{
			delete[] SiteEmpiricalCount[i];
		}
		delete[] SiteEmpiricalCount;
	}

	delete[] SpeciesNames;

	if (RawData)	{
		for (int i=0; i<Ntaxa; i++)	{
			delete[] RawData[i];
		}
		delete[] RawData;
		delete[] RawAlphabet;
	}
	delete[] Alphabet;
	delete[] AlphabetSet;


	if (ZipData)	{

		for (int i=0; i<Ntaxa; i++)		{
			delete[] ZipData[i];
		}
		delete[] ZipData;
	}

	delete[] SwapPermut;

	if (Data)	{
		for (int i=0; i<Ntaxa; i++)	{
			delete[] Data[i];
		}
		delete[] Data;
	}
	if (ContData)	{
		for (int i=0; i<Ntaxa; i++)	{
			delete[] ContData[i];
		}
		delete[] ContData;
	}
	if (ContMissingData)	{
		for (int i=0; i<Ntaxa; i++)	{
			delete[] ContMissingData[i];
		}
		delete[] ContMissingData;
	}

	delete[] OrbitSize;
	delete[] ZipSize;
	if (Orbit)	{
		for (int i=0; i<Nsite; i++)		{
			delete[] Orbit[i];
		}
		delete[] Orbit;
	}

	if (Indices)	{
		for (int i=0; i<Nsite; i++)		{
			delete[] Indices[i];
		}
		delete[] Indices;
	}
	
	for (int i=0; i<MaxTopoMoveTypeNumber; i++)	{
		for (int j=0; j<MaxSuperMoveTypeNumber; j++)	{
			for (int k=0; k<MaxMoveTypeNumber; k++)	{
				delete[] NIterationArray[i][j][k];
				delete[] deltaArray[i][j][k];
				delete[] SuccessArray[i][j][k];
				delete[] NCallArray[i][j][k];
				delete[] TimeArray[i][j][k];
				delete[] NArray[i][j][k];
			}
			delete[] NIterationArray[i][j];
			delete[] deltaArray[i][j];
			delete[] SuccessArray[i][j];
			delete[] NCallArray[i][j];
			delete[] TimeArray[i][j];
			delete[] NArray[i][j];

			delete[] SuperNIterationArray[i][j];
			delete[] MoveTypeArray[i][j];
		}
		delete[] NIterationArray[i];
		delete[] deltaArray[i];
		delete[] SuccessArray[i];
		delete[] NCallArray[i];
		delete[] TimeArray[i];
		delete[] NArray[i];

		delete[] SuperNIterationArray[i];
		delete[] MoveTypeArray[i];

		delete[] TopoNIterationArray[i];
		delete[] TopoNArray[i];
		delete[] MoveTypeNumber[i];
	}

	delete[] NIterationArray;
	delete[] deltaArray;
	delete[] SuccessArray;
	delete[] NCallArray;
	delete[] TimeArray;
	delete[] NArray;

	delete[] SuperNIterationArray;
	delete[] MoveTypeArray;

	delete[] TopoNIterationArray;
	delete[] TopoNArray;
	delete[] MoveTypeNumber;

	delete[] TopoMoveTypeArray;

	delete[] TrialSwapArray;
	delete[] AcceptedSwapArray;

	delete[] Beta;

	delete[] Gene;
	delete[] GeneFirstSite;
	delete[] GeneSize;

	delete[] RecodingAlphabet;
	delete[] RecodingTable;

	if (Ncat)	{
		for (int i=0; i<Ncat; i++)	{
			delete[] statfix[i];
		}
		delete[] statfix;
		delete[] weight;
		if (StatFlag)	{
			for (int i=0; i<Ncat; i++)	{
				delete[] StatFlag[i];
			}
			delete[] StatFlag;
		}
	}

	if (verbose)	{
		cerr << "mParam delete ok\n";
		cerr.flush();
	}
}


//---------------------------------------------------------------------------------
//		 Reset()
//---------------------------------------------------------------------------------

void	MCParameters::Reset()	{

	for (int k=0; k<TopoMoveTypeNumber; k++)	{
		for (int j=0; j<SuperMoveTypeNumber[k]; j++)	{
			for (int i=0; i<MoveTypeNumber[k][j]; i++)	{
				for (int l=0; l<MaxChainNumber; l++)	{
					SuccessArray[k][j][i][l] = 0;
					NCallArray[k][j][i][l] = 0;
					TimeArray[k][j][i][l] = 0;
				}
			}
		}
	}

	for (int chain = 0; chain < Nchain ; chain++)	{
		TrialSwapArray[chain] = 0;
		AcceptedSwapArray[chain] = 0;
	}
	
	SwapPermut = new int[Nchain];
	for (int i=0; i<Nchain; i++)	{
		SwapPermut[i] = i;
	}

	StatInfCount = 0;
	RateInfCount = 0;
	LengthInfCount = 0;
	LogProbInfCount = 0;
	SubOverflowCount = 0;

	HowManySaved = 0;
	Ngen = 0;

	Zeta = InitBeta;
	double initbeta = ArcTh(Zeta);
	
	if (MSMode == Thermo)	{
		Beta0 = initbeta;
	}
	else if (MSMode == SA)	{
		Beta0 = initbeta;
	}
	else if (MSMode == RAS)	{
		RASModelSwitch = initbeta;
	}
	else if (MSMode == BL)	{
		BLModelSwitch = initbeta;
	}
	else if (MSMode == FlexClock)	{
		FlexClockModelSwitch = initbeta;
	}
	else if (MSMode == SUB)	{
		SUBModelSwitch = initbeta;
	}
	else if (MSMode == Separate)	{
		SeparateModelSwitch = initbeta;
	}
	else if (MSMode == Hetero)	{
		HeteroModelSwitch = initbeta;
	}
	else if (MSMode == ClockPrior)	{
		ClockPriorModelSwitch = initbeta;
	}
	else if (MSMode == Autocorrel)	{
		AutoCorrelModelSwitch = initbeta;
	}
	else if (MSMode == LengthMixGamma)	{
		GetCurrentState()->LengthGamma = InitBeta;
		GetNextState()->LengthGamma = InitBeta;
	}
	else if (MSMode == RateGamma)	{
		GetCurrentState()->gamma = InitBeta;
		GetNextState()->gamma = InitBeta;
	}
	else if (MSMode == CATAlpha)	{
		GetCurrentState()->alpha = InitBeta;
		GetNextState()->alpha = InitBeta;
		/*
		if (AlphaMin == 0)	{
			GetCurrentState()->alpha = AlphaMax;
			GetNextState()->alpha = AlphaMax;
		}
		else	{
			double alpha = exp((1 - InitBeta) * log(AlphaMin) + InitBeta * log(AlphaMax)) ;
			cerr << "init alpha : " << alpha << '\n';
			GetCurrentState()->alpha = alpha;
			GetNextState()->alpha = alpha;
		}
		*/
	}
	logBF = 0;
	logBFmin = 0;
	logBFmax = 0;
	logBFrev = 0;
	logBFminrev = 0;
	logBFmaxrev = 0;
	BFcount = 0;

	if (ActivateSumOverModes)	{
		SumOverModes = Yes;
	}
	if (ActivateSumOverRateModes)	{
		SumOverRateModes = Yes;
	}
	
	if (NormalApprox == Yes)	{
		RegisterCov();
		// create phylobayes here
		currentState = new PhyloBayes*[Nchain];
		nextState = new PhyloBayes*[Nchain];
		for (int chain=0; chain<Nchain; chain++)	{
			currentState[chain] = new PhyloBayes(this);
			nextState[chain] = new PhyloBayes(this);
			currentState[chain]->SetClockTree();
		}
	}
	Update();
}

//---------------------------------------------------------------------------------
//		 Update()
//---------------------------------------------------------------------------------

void	MCParameters::Update()	{

	if (! CodonCode)	{
		CodonCode = new int[Ncodon];
		if (GeneticCode == Universal)	{
			for (int k=0; k<Ncodon; k++)	{
				CodonCode[k] = UniCodonCode[k];
			}
		}
		else if (GeneticCode == MtMam)	{
			for (int k=0; k<Ncodon; k++)	{
				CodonCode[k] = MtMamCodonCode[k];
			}
		}
	}

	if (MSMode == None) {
		if (	(!ActivateSumOverModes) && (!ActivateSumOverRateModes) && 
			((RASModelSwitch == 0) || (RASModelSwitch == 1)) &&
			((SeparateModelSwitch == 0) || (SeparateModelSwitch == 1)) &&
			((SUBModelSwitch == 0) || (SUBModelSwitch == 1))	)	{
			
				SimpleSampling = Yes;
		}
		else	{
				SimpleSampling = No;
		}
	}

	if (MSMode == None)	{
		if (Nchain > 1)	{
			for (int chain=0; chain<Nchain; chain++)	{
				currentState[chain]->SetBeta(Beta[chain]);
				nextState[chain]->SetBeta(Beta[chain]);
				nextState[chain]->Clone(currentState[chain]);
				currentState[chain]->Update();
				nextState[chain]->Update();
			}
		}
		else	{
			currentState[0]->SetBeta(Beta0);
			nextState[0]->SetBeta(Beta0);
			nextState[0]->Clone(currentState[0]);
			currentState[0]->Update();
			nextState[0]->Update();
		}
	}

	else	{
		if (Nchain > 1)	{
			cerr << "error : thermo only with one chain\n";
			exit(1);
		}
		currentState[0]->SetBeta(Beta0);
		nextState[0]->SetBeta(Beta0);
		nextState[0]->Clone(currentState[0]);
		currentState[0]->Update();
		nextState[0]->Update();
	}
}

void MCParameters::SetInitTree(Tree& tree)	{

	for (int chain=0; chain<Nchain; chain++)	{
		GetCurrentState(chain)->SetTree(tree);
		// GetCurrentState(chain)->IsInitTree = Yes;
		// GetCurrentState(chain)->InitTree = Fixed;
	}
}


double MCParameters::ArcTh(double zeta)	{

	if (zeta == 0)	{
		return 0;
	}
	if (zeta == 1)	{
		return 1;
	}
	if (ArcThermo)	{
		double x = zeta - 0.5;
		return exp(alpha1 * x) / (exp(alpha1 * x) + exp(-alpha2 * x));
	}
	else	{
		return zeta;
	}
}

double MCParameters::ArcZ()	{
	if (alpha1 != alpha2)	{
		cerr << "error in ArcZ: alpha1 should be equal to alpha2\n";
		exit(1);
	}
	double a = alpha1;
	if (ArcThermo)	{
		return 1.0/4/a/a * (2*a + 2*exp(a) - 2*exp(-a) + exp(2*a) - exp(-2*a));
	}
	return 1;
}

double MCParameters::ArcH(double zeta)	{
	double a = alpha1;
	double x = zeta-0.5;
	if (ArcThermo)	{
		if ((zeta == 0)	|| (zeta == 1)) return 0;
		return 1.0/2/a * (1 + 2*exp(-2*a*x) + exp(-4*a*x));
	}
	return 1;
}

double MCParameters::ArcdH(double zeta)	{
	double a = alpha1;
	double x = zeta-0.5;
	if (ArcThermo)	{
		return -2 * (exp(-2*a*x) + exp(-4*a*x));
	}
	return 0;
}

double MCParameters::ArcW(double zeta)	{
	return -2 * ArcH(zeta) * ArcdH(zeta);
}

//---------------------------------------------------------------------------------
//		 Move()
//---------------------------------------------------------------------------------

// returns 1 if Move has been accepted

void	MCParameters::Move()	{

	// first, move the chains from right to left
	for (int chain=0; chain < Nchain; chain++)	{
		Move(chain);
	}

	// second, swap the chains 
	if ((Nchain > 1) && SwapFreq) 	{
		for (int chain=Nchain-1; chain>0; chain--)	{
			Swap(chain);
		}
		for (int chain=1; chain<Nchain; chain++)	{
			Swap(chain);
		}
	}
}


//---------------------------------------------------------------------------------
//		 NgenPerCycle()
//---------------------------------------------------------------------------------

int	MCParameters::NgenPerCycle()	{
	int n = 0;
	int chain = 0;
	for (int topo=0; topo<TopoMoveTypeNumber; topo++)	{
		for (int j=0; j<TopoNIterationArray[topo][chain]; j++)	{
			for (int super=0; super<SuperMoveTypeNumber[topo]; super++)	{
				for (int j=0; j<SuperNIterationArray[topo][super][chain]; j++)	{
					for (int move=0; move<MoveTypeNumber[topo][super]; move++)	{
						n+= NIterationArray[topo][super][move][chain];
					}
				}
			}
		}
	}
	return n;
}

//---------------------------------------------------------------------------------
//		 Move(int chain)
//---------------------------------------------------------------------------------

void	MCParameters::Move(int chain)	{


	for (int topo=0; topo<TopoMoveTypeNumber; topo++)	{
		for (int j=0; j<TopoNIterationArray[topo][chain]; j++)	{
			TopoMove(chain, topo);
		}
	}
}


//---------------------------------------------------------------------------------
//		 TopoMove(int chain, int topo)
//---------------------------------------------------------------------------------

void	MCParameters::TopoMove(int chain, int topo)	{

	for (int super=0; super<SuperMoveTypeNumber[topo]; super++)	{
		for (int j=0; j<SuperNIterationArray[topo][super][chain]; j++)	{
			SuperMove(chain, topo, super);
		}
	}
}

//---------------------------------------------------------------------------------
//		 SuperMove(int chain, int topo, int super)
//---------------------------------------------------------------------------------

void	MCParameters::SuperMove(int chain, int topo, int super)	{

	for (int move=0; move<MoveTypeNumber[topo][super]; move++)	{
		for (int j=0; j<NIterationArray[topo][super][move][chain]; j++)	{
			Move(chain, topo, super, move);
		}
	}
}

//---------------------------------------------------------------------------------
//		 Move(int chain, int topo, int super, int move)
//---------------------------------------------------------------------------------


void	MCParameters::Move(int chain, int topo, int super, int move)	{

	int check = 0;
	MoveType movetype = MoveTypeArray[topo][super][move];
	if (check)	{
		cerr << MoveTypeName[movetype] << '\n';
		cerr.flush();
	}
	movetime.Reset();
	movetime.Start();

	double success = GetNextState(chain)->Move(
						MoveTypeArray[topo][super][move],
						deltaArray[topo][super][move][chain],
						NArray[topo][super][move][chain],
						GetCurrentState(chain));
	movetime.Stop();
	double duration = movetime.GetTime();
	NCallArray[topo][super][move][chain]++;
	SuccessArray[topo][super][move][chain]+= success;
	TimeArray[topo][super][move][chain]+= duration;
	
	if ((Nchain > 1) && SwapFreq) 	{
		for (int chain=Nchain-1; chain>0; chain--)	{
			Swap(chain);
		}
		for (int chain=1; chain<Nchain; chain++)	{
			Swap(chain);
		}
	}
	// check
	if (check)	{
		cerr << "move ok, check\n";
		cerr.flush();

		GetCurrentState()->checkbl();
		if ((movetype != GammaIntegral) && (movetype != ResampleStatMove) && (movetype != ResampleLengthMove) && (movetype != ResampleRateMove) && (movetype != ResampleRelativeRateMove) && (movetype != ResampleSubMove)) {
		double logpriorn1 = GetNextState()->mLogPrior;
		double logsampn1 = GetNextState()->mLogSampling;
		GetNextState()->Update();
		double logpriorn2 = GetNextState()->mLogPrior;
		double logsampn2 = GetNextState()->mLogSampling;
		double logpriorc1 = GetCurrentState()->mLogPrior;
		double logsampc1 = GetCurrentState()->mLogSampling;
		GetCurrentState()->Update();
		double logpriorc2 = GetCurrentState()->mLogPrior;
		double logsampc2 = GetCurrentState()->mLogSampling;
		if (fabs(logpriorn2 - logpriorn1) > 1e-6)	{
			cerr << "logprior next : " << logpriorn1 << '\t' << logpriorn2 << '\n';
		}
		if (fabs(logsampn2 - logsampn1) > 1e-6)	{
			cerr << "logsamp next : " << logsampn1 << '\t' << logsampn2 << '\n';
		}
		if (fabs(logpriorc2 - logpriorc1) > 1e-6)	{
			cerr << "logprior current: " << logpriorc1 << '\t' << logpriorc2 << '\n';
		}
		if (fabs(logsampc2 - logsampc1) > 1e-6)	{
			cerr << "logsamp current: " << logsampc1 << '\t' << logsampc2 << '\n';
		}
		if (fabs(logpriorc2 - logpriorn2) > 1e-6)	{
			cerr << "logprior clone: " << logpriorc2 << '\t' << logpriorn2 << '\n';
		}
		if (fabs(logsampc2 - logsampn2) > 1e-6)	{
			cerr << "logsamp clone: " << logsampc2 << '\t' << logsampn2 << '\n';
		}

		}
	}	
}

void	MCParameters::Swap(int chain1)	{

		int chain2 = chain1 - 1;
		double beta1 = GetCurrentState(chain1)->Beta;
		double beta2 = GetCurrentState(chain2)->Beta;

		double logRatio = 0;

		double logsamp1 = GetCurrentState(chain1)->logSampling();
		double logsamp2 = GetCurrentState(chain2)->logSampling();

		logRatio = - (logsamp1 - logsamp2)* (beta1 - beta2);

		int Accepted = (-log(Random::Uniform()) > logRatio);

		if (Accepted)	{

			PhyloBayes* temp = GetCurrentState(chain1);
			SetCurrentState(chain1,GetCurrentState(chain2));
			SetCurrentState(chain2,temp);

			temp = GetNextState(chain1);
			SetNextState(chain1, GetNextState(chain2));
			SetNextState(chain2,temp);

			GetCurrentState(chain1)->SetBeta(beta1);
			GetNextState(chain1)->SetBeta(beta1);

			GetCurrentState(chain2)->SetBeta(beta2);
			GetNextState(chain2)->SetBeta(beta2);

			AcceptedSwapArray[chain1]++;

			int tmp = SwapPermut[chain1];
			SwapPermut[chain1] = SwapPermut[chain2];
			SwapPermut[chain2] = tmp;
		}

		TrialSwapArray[chain1]++;
	}

//---------------------------------------------------------------------------------
//		 Init()
//---------------------------------------------------------------------------------


void MCParameters::InitClock(int clock, int genemode)	{

	if (clock != 0)	{
		ActivateClock = Yes;
	}
	BLModelSwitch = 0;
	FlexClockModelSwitch = 1;
	ClockModel = KT;
	FlexClockModel = WhiteNoise;
	if (clock == 1)	{
		ClockModel = Strict;
	}
	else if (clock == 2)	{
		ClockModel = KT;
	}
	else if (clock == 3)	{
		ClockModel = CIRrigid;
	}
	else if (clock == 4)	{
		ClockModel = Strict;
		BLModelSwitch = 1;
		FlexClockModel = WhiteNoise;
	}
	else if (clock == 5)	{
		ClockModel = Strict;
		BLModelSwitch = 1;
		FlexClockModel = UGam;
	}
	else if (clock == 6)	{
		ClockModel = AutoRegressive;
	}
	else if (clock == 7)	{
		ClockModel = Strict;
		BLModelSwitch = 1;
		FlexClockModel = CIRflex;
	}
	else if (clock == 8)	{
		ClockModel = UndebiasedKT;
	}

	else	{
		cerr << "error in MCParamaters::Init: do not recognise clock model\n";
		exit(1);
	}

	if (genemode == 0)	{
		SeparateModelSwitch = 0;
	}
	else if (genemode == 1)	{
		SeparateModelSwitch = 1;
		GeneBLMultiplier = No;
		GeneBLMode = Yes;
		SeparateRhoPrior = No;
	}
	else if (genemode == 2)	{
		SeparateModelSwitch = 1;
		GeneBLMultiplier = Yes;
		GeneBLMode = Yes;
		SeparateRhoPrior = No;
		BLModelSwitch = 1;
		ClockModel = CIRrigid;
	}
}


void MCParameters::Init(RRMode rr, int empfreq, int ncat, int statcenter, RASMode ras, int discrate, int fixmeanlength, double meanlength)	{
			
	RegisterWithData();

	NZipModeMax = Nsite +10 ;

	if ((NH == 2) || (NH >= 4))	{
		NHNcatMax = Nnode;
	}

	FixNmode = 1;
	// FixNRateMode = 1;
	int catfix = 0;
	ActivateSumOverModes = No;
	ActivateSumOverRateModes = No;
	if (statcenter)	{	
		ModeStatPrior = MultiGamma;
	}
	else	{
		ModeStatPrior = Flat;
	}
	FixStatCenter = 1;

	RefFastCompute = Yes;
	if (rr == poisson)	{
		ModeFastCompute = Yes;
		ModePoisson = Yes;
	}
	else	{
		if (ZipGTR)	{
			ModeFastCompute = Yes;
		}
		else	{
			ModeFastCompute = No;
		}
		ModePoisson = No;
	}
 	if ((ncat == -3) || (ncat == -2)) {
		NmodeMax = Ncat;
		/*
		if (Ncat < 100)	{
			ActivateSumOverModes = Yes;
			SumOverModes = Yes;
		}
		*/
		if (ActivateEM)	{
			ActivateSumOverModes = Yes;
			SumOverModes = Yes;
		}
		if (Tempering)	{
			ActivateSumOverModes = Yes;
			SumOverModes = Yes;
		}
		if (Ncat == 0)	{
			cerr << "error in MCParamaters::Init: Ncat = " << Ncat << '\n';
			exit(1);
		}
		if (ncat == -3)	{
			FixedWeights = 1;
		}
		ncat = Ncat;
		catfix = 1;
	}
	else if (ncat == -1)	{ // max
		NmodeMax = Nsite;
		ncat = Nsite;
		FixStatCenter = 0;
	}
	else if (ncat == 0)	{ // cat
		ModePrior = DirichletProcess;
		NmodeMax = Nsite + 5;
		if (ModePoisson && (HeteroMode == Homo))	{
			ncat = Nsite;
		}
		else	{
			ncat = 1;
		}
		FixNmode = 0;
		// if (ModePoisson)	{
			FixStatCenter = 0;
		// }
		if (ModeFastCompute && (! ModePoisson))	{
			EmpiricalDP = 0;
			EmpiricalDP = 1;
			IncrementalDP = 0;
		}
		else	{
			EmpiricalDP = 1;
		}
		if (Beta0 == 0)	{
			EmpiricalDP = 0;
			IncrementalDP = 0;
		}
	}
	else	{
		NmodeMax = ncat;
		if ((ncat > 1) && (ncat < 100)) 	{
			ActivateSumOverModes = Yes;
			SumOverModes = Yes;
			ModeStatPrior = Flat;
			FixStatCenter = 1;
		}
		else	{
			ActivateSumOverModes = No;
			SumOverModes = No;
			FixStatCenter = 0;
		}
	}
	NRateModeMax = 1;
	if (ModeStatPrior == Flat)	{
		FixStatCenter = 1;
	}


	RASModelSwitch = 1;
	RatePrior = GammaInv;
	if (discrate)	{
		GammaNcat = discrate;
		ActivateSumOverRateModes = Yes;
		SumOverRateModes = Yes;
		RASModelSwitch = 0;
		NRateModeMax = discrate;
	}
	if (ras == dprate)	{
		RASModelSwitch = 0;
		NRateModeMax = Nsite + 5;
	}
	if (ras == uni)	{
		RASModelSwitch = 0;
		NRateModeMax = 1;
		GammaNcat = 1;
	}
		
	if (! DataOK)	{
		cerr << "should specify data file before init\n";
		exit(1);
	}
	if (IncrementalDP)	{
		IncrementalDP = Nchain;
	}

	if (MBL)	{
		NMixBLMax = Nsite + 10;
		// NMixBLMax = Nsite/10;
	}
	currentState = new PhyloBayes*[Nchain];
	nextState = new PhyloBayes*[Nchain];
	for (int chain=0; chain<Nchain; chain++)	{
		currentState[chain] = new PhyloBayes(this);
		nextState[chain] = new PhyloBayes(this);
	
		PhyloBayes* pb = currentState[chain];

		pb->InitModeStat = Empirical;
		if (rr == poisson)	{
			pb->InitModeRR = Poisson;
		}
		else if (rr == wag)	{
			pb->InitModeRR = WAG;
			pb->InitModeStat = WAG;
		}
		else if (rr == waghssp)	{
			pb->InitModeRR = WAGHSSP;
			pb->InitModeStat = WAGHSSP;
		}
		else if (rr == lg)	{
			pb->InitModeRR = LG;
			pb->InitModeStat = LG;
		}
		else if (rr == mtrev)	{
			pb->InitModeRR = mtREV;
			pb->InitModeStat = mtREV;
		}
		else if (rr == mtzoa)	{
			pb->InitModeRR = mtZOA;
			pb->InitModeStat = mtZOA;
		}
		else if (rr == mtart)	{
			pb->InitModeRR = mtART;
			pb->InitModeStat = mtART;
		}
		else if (rr == jtt)	{
			pb->InitModeRR = JTT;
			pb->InitModeStat = JTT;
		}
		else if (rr == gtr)	{
			pb->InitModeRR = WAG;
		}
		else if (rr == custom)	{
			pb->InitModeRR = Cstom;
		}
		if (empfreq == 2)	{
			pb->InitModeStat = SiteEmpirical;
		}
		if (empfreq == 1)	{
			pb->InitModeStat = Empirical;
		}
		if (catfix)	{
			pb->InitModeStat = CatFix;
			/*
			ActivateSumOverModes = Yes;
			SumOverModes = Yes;
			*/
		}
		if ((ras != uni) && (RASModelSwitch == 1))	{
			pb->InitRate = Random;
		}
		else	{
			pb->InitRate = Uniform;
		}

		if (ActivateEM)	{
			if (pb->InitGamma == Random)	{
				pb->gamma = Random::sExpo();
			}
		}
		else if (FixRateAlpha)	{
			pb->gamma = FixRateAlpha;
		}
		else if (MSMode == RateGamma)	{
			pb->gamma = InitBeta;
		}
		else	{
			pb->gamma = 1;
		}
		pb->IsInitGamma = Yes;

		pb->MeanLength = meanlength;
		pb->IsInitMeanLength = Yes;
		pb->ModeStatAlpha = 20;

		pb->Nmode = ncat;
		pb->pconst = 0;

		if (ncat)	{
			pb->IsInitNmode = Yes;
		}
		if (ActivateEM)	{
			if (RandomInitLength)	{
				pb->IsInitBranchLength = No;
			}
			else	{
				pb->IsInitBranchLength = Yes;
			}
		}
		else	{
			pb->IsInitBranchLength = No;
		}

		pb->IsInitModeStatAlpha = Yes;
				
		pb->Initialise();
	}

	// determine what is fixed
	FixStat = 0;
	FixGamma = 0;
	if (FixRateAlpha)	{
		FixGamma = 1;
	}
	FixMeanLength = fixmeanlength;
	FixRR = 0;
	FixLength = 0;
	FixRate = 0;
	FixTopo = 0;
	FixPconst = 1;

	if (ras == uni)	{
		FixGamma = 1;
		FixRate = 1;
	}

	if (MSMode == RateGamma)	{
		FixGamma = 1;
	}

	if (ras == invgam)	{
		FixPconst = 0;
		cerr << "pconst will move\n";
	}

	if (discrate)	{
		FixRate = 1;
	}

	if ((catfix && (Qmode)) || (rr != gtr))	{
		FixRR = 1;
	}
	
	if (empfreq || catfix)	{
		FixStat = 1;
	}
}

void ::MCParameters::ActivateDiscreteGamma(int discrate)	{

	GammaNcat = discrate;
	ActivateSumOverRateModes = Yes;
	RASModelSwitch = 0;
	NRateModeMax = discrate;
	
	GetCurrentState()->ActivateDiscreteGamma();
	GetNextState()->ActivateDiscreteGamma();
	SumOverRateModes = Yes;
	Update();
}


	
void MCParameters::MoveSiteVariables(int conjugate)	{

	if ((! FixStat) && (NmodeMax > 1))	{
		MoveStat(conjugate);
	}
	if (! FixRate)	{
		MoveRate(conjugate);
	}
}

	
void MCParameters::MoveRate(int conjugate)	{

	if (SumOverRateModes == Yes)	{
		cerr << "error in MCParam::MoveRate: not with sum over rate modes\n";
		exit(1);
	}
	if (conjugate)	{
		GetNextState()->ResampleSub();
		GetNextState()->ResampleRate(GetCurrentState());
	}
	else	{
		GetNextState()->Move(Rate, 1, 1, GetCurrentState());
		GetNextState()->Move(Rate, 1, 5, GetCurrentState());
	}
}

void MCParameters::MoveStat(int conjugate)	{
	
	int Nrep = 10;
	int N = NmodeMax;
	if (N>10)	{
		N = 10;
	}
	if (ActivateSumOverModes)	{
		GetNextState()->Move(UpdateSumModeMove, 1, 10, GetCurrentState());
		GetNextState()->Move(SumModeStationary, 1, 10, GetCurrentState());
		GetNextState()->Move(SumModeStationary, 0.1, 10, GetCurrentState());
		for (int rep=0; rep<Nrep; rep++)	{
			GetNextState()->Move(SumModeWeight, 100, N, GetCurrentState());
			GetNextState()->Move(SumModeWeight, 1000, N, GetCurrentState());
			GetNextState()->Move(SumModeWeight, 10000, N, GetCurrentState());
		}
		if (conjugate)	{
			GetNextState()->Move(ResampleModeAffMove, 1, 1, GetCurrentState());
			GetNextState()->Move(ResampleSubMove, 1, 1, GetCurrentState());
			GetNextState()->Move(SwitchModeIntegral, 5, 1, GetCurrentState());
			GetNextState()->Move(ResampleStatMove, 1, 1, GetCurrentState());
		}
	}
	else	{
		if (conjugate)	{
			GetNextState()->ResampleSub();
			GetNextState()->ResampleStat(1,10,GetCurrentState());
		}
		else	{
			for (int rep=0; rep<Nrep; rep++)	{
				GetNextState()->Move(ModeStationary, 1, 2, GetCurrentState());
				// GetNextState()->Move(SumModeStationary, 0.1, 2, GetCurrentState());
			}
		}
	}
}


void MCParameters::InitMoveBF()	{

	TopoMoveTypeArray[TopoMoveTypeNumber] = FixedTopo;
	TopoNIterationArray[TopoMoveTypeNumber][0] = 1;
	int topo = TopoMoveTypeNumber;
	int super = SuperMoveTypeNumber[topo];
	SuperNIterationArray[topo][super][0] = 100;

	PushMove(muMove,1,1,1);
	//PushMove(muRho,1,1,1);
	PushMove(sigmaMove,1,1,1);
	if ((ClockModel == CIRrigid) || (FlexClockModel == CIRflex)) 	{
		PushMove(thetaMove,1,1,1);
	}
	PushMove(muMove,1,0.1,1);
	//PushMove(muRho,1,0.1,1);
	PushMove(sigmaMove,1,0.1,1);
	if ((ClockModel == CIRrigid) || (FlexClockModel == CIRflex)) 	{
		PushMove(thetaMove,1,1,1);
	}
	PushMove(addRhoMove,1,1,1);
	PushMove(addRhoMove,1,0.1,1);
	PushMove(OneBranchLength,1,10,1);
	PushMove(OneBranchLength,1,1,1);
	PushMove(OneBranchLength,1,0.1,1);
	PushMove(oneTimeNodeMove,1,1,1);
	PushMove(oneTimeNodeMove,1,0.1,1);
	// PushMove(oneTimeNodeMove,1,0.01,1);
	if (Chi == -1)	{
		if (TimePrior == BD)	{
			PushMove(chi,100,10,1);
			PushMove(chi,100,1,1);
			PushMove(chi,100,0.1,1);
			PushMove(chi2,100,10,1);
			PushMove(chi2,100,1,1);
			PushMove(chi2,100,0.1,1);
		}
		else if (TimePrior == Dirich)	{
			PushMove(chi,100,1,1);
			PushMove(chi,100,0.1,1);
		}
	}
	SuperMoveTypeNumber[topo]++;
	TopoMoveTypeNumber++;
}


void MCParameters::InitMove()	{

	TopoMoveTypeArray[TopoMoveTypeNumber] = FixedTopo;
	for (int l=0; l<Nchain; l++)	{
		TopoNIterationArray[TopoMoveTypeNumber][l] = 1;
	}
	int topo = TopoMoveTypeNumber;
	int super = SuperMoveTypeNumber[topo];

	if (! NormalApprox)	{
		for (int l=0; l<Nchain; l++)	{
			SuperNIterationArray[topo][super][l] = 1;
		}

		int fixnratemode = 1;
		if (NRateModeMax == Nsite + 5)	{
			fixnratemode = 0;
		}
		int meanlength = (! FixMeanLength) && (! ActivateClock);
		int statcenter = ! FixStatCenter;
		int rr = ! FixRR;
		int length = ! FixLength;
		int stat = ! FixStat;
		int rate = ! FixRate; 
		int Niter = 10;
		int freetopo = (! FixTopo) && (! ActivateClock);
		int catfix = (FixStat && GetCurrentState()->Nmode==Ncat);
		/*
		int gamma = ! FixGamma;
		if (gamma)	{
			PushMove(Gamma,10,1,1);
			PushMove(Gamma,10,0.1,1);
		}
		*/

		/*
		if (covarion)	{
			PushMove(Xi,1,1,1);
			PushMove(Xi,1,0.1,1);
			PushMove(InvProb,1,1,1);
			PushMove(InvProb,1,0.1,1);
		}
		*/

		if (ActivateSumOverModes && (! Tempering))	{
			// PushMove(UpdateSumModeMove,1,1,1);
			PushMove(ResampleModeAffMove,1,1,1);
		}
		if (meanlength)	{
			PushMove(MeanLengthM,Niter,1,1);
			PushMove(MeanLengthM,Niter,0.1,1);
			if (LengthPrior == GammaDistributed)	{
				PushMove(VarLengthM,Niter,1,1);
				PushMove(VarLengthM,Niter,0.1,1);
			}
		}
		if (statcenter)	{
			if (ModeStatPrior == Flat)	{
				cerr << "error in MCParam::InitMove: moving statcenter, but modestatprior flat\n";
				exit(1);
			}
			// PushMove(ModeStatCenterM,2,10000,5);
			// PushMove(ModeStatCenterM,2,100000,5);
			// PushMove(ModeStatCenterM,2,1000000,5);
			PushMove(ModeStatCenterM,20,1,2);
			PushMove(ModeStatCenterM,20,0.1,2);
			PushMove(ModeStatCenterM,20,0.01,2);
			PushMove(ModeStatAlphaM,20,1,1);
			PushMove(ModeStatAlphaM,20,0.1,1);
			PushMove(ModeStatAlphaM,20,0.01,1);
		}
		if ((! FixNmode) || (GetCurrentState()->Nmode > 1))	{
			if (MSMode != CATAlpha)	{
				PushMove(Alpha,Niter,10,1);
				PushMove(Alpha,Niter,1,1);
				PushMove(Alpha,Niter,0.1,1);
			}
		}
		if (! fixnratemode)	{
			PushMove(RateAlphaM,Niter,10,1);
			PushMove(RateAlphaM,Niter,1,1);
			PushMove(RateAlphaM,Niter,0.1,1);
		}
		if (NH && NHPrior)	{
			if (! NHWishart)	{
				PushMove(NHVarM,10,1,1);
				PushMove(NHVarM,10,0.1,1);
			}
			else	{
				PushMove(NHContVarM,10,1,1);
				PushMove(NHContVarM,10,0.1,1);
			}
			if (catfix && NHAutoRegressive)	{
				PushMove(NHStatCenterM,10,1,1);
				PushMove(NHStatCenterM,10,0.1,1);
				PushMove(NHStatAlphaM,10,1,1);
				PushMove(NHStatAlphaM,10,0.1,1);
			}
			else if (catfix && (NH<4))	{
				PushMove(NHStatCenterM,10,1,1);
				PushMove(NHStatCenterM,10,0.1,1);
			}
				
			// PushMove(NHVarStatAlphaM,10,1,1);
			// PushMove(NHVarStatAlphaM,10,0.1,1);
			if ((NHPrior == 2) && (!NHWishart))	{
				PushMove(NHCovM,2,1,10);
				PushMove(NHCovM,2,0.1,10);
				PushMove(NHCovM,10,1,2);
				PushMove(NHCovM,10,0.1,2);
			}
		}
		if (! freetopo)	{
			if ((! ActivateClock) && length && WithPartial && ((NH>=4) || (! Conjugate) || NHLengthSwitch))	{
				PushMove(EnterPartialM,1,1,1);
				PushMove(OneBranchLengthPartial,1,1,1);
				PushMove(OneBranchLengthPartial,1,0.1,1);
				PushMove(LeavePartialM,1,1,1);
			}
		}
		if (freetopo)	{
			if (WithPartial)	{
				double epsilon = 1;
				if (MBL == 2)	{
					epsilon = 0;
				}
				if ((MBL != 2) && length && ((NH>=4) || (! Conjugate) || NHLengthSwitch))	{
					PushMove(EnterPartialM,1,1,1);
					PushMove(OneBranchLengthPartial,1,1,1);
					PushMove(OneBranchLengthPartial,1,0.1,1);
				}
				PushMove(SPRLocal,1,epsilon/10,10);
				PushMove(SPRLocal,1,epsilon/10,0);
				PushMove(GibbsSPR,1,0,0);

				// ???
				if (NH>=4)	{
					PushMove(SPRLocal,1,0.3,10);
					PushMove(SPRLocal,1,0.3,0);
					PushMove(SPRLocal,1,1,10);
					PushMove(SPRLocal,1,1,0);
				}

				if (NH){
					if ((NH != 2) && (NH < 4)) {
						PushMove(NHSPRLocalPartial,1,0,10);
						// PushMove(NHSPRPartial,1,0,3);
						// ???
						PushMove(GibbsSPR,1,0,0);
					}
				}
				PushMove(LeavePartialM,1,1,1);

				if (! FixRoot)	{
					/*
					if (NH)	{
						if ((NH != 2) && (NH < 4)) {
							PushMove(NHRoot,2,0,0);
						}
						else	{
							PushMove(Root,2,0,0);
						}
					}
					*/
					if (PopEff && (!FixRoot))	{
						PushMove(Root,2,0,0);
					}
				}
			}
			else	{
				PushMove(SPR,2,0,0);
				PushMove(SPR,2,0,3);
				PushMove(SPR,2,0,6);
				PushMove(SPR,2,0,9);
				PushMove(SPR,2,0,-1);
				if (PopEff && (!FixRoot))	{
					PushMove(Root,2,0,0);
				}
				if (NH)	{
					if ((NH != 2) && (NH < 4)) {
						// PushMove(NHRoot,1,0,0);
						PushMove(NHSPR,1,0,-1);
						// PushMove(Root,2,0,0);
					}
					else	{
						//  PushMove(Root,2,0,0);
					}
				}

				PushMove(TBR,2,1,1);
				PushMove(Global,2,0.1,1);
				PushMove(Global,2,0.01,1);
				PushMove(Local,2,1,1);
				PushMove(NodeSliding,2,10,1);

				PushMove(ResampleSubMove,1,1,1);
				PushMove(ResampleLengthMove,1,1,1);
			}
		}
		if (ActivateSumOverRateModes && (! Tempering))	{
			PushMove(UpdateSumRateModeMove,1,1,1);
			PushMove(ResampleRateModeAffMove,1,1,1);
		}
		int nrep = 1;
		if (Qmode && rr)	{
			nrep = 2;
		}
		for (int rep=0; rep<nrep; rep++)	{	
		if (rr)	{
			if (Conjugate && (MutMode != 2) && (MutMode != 3))	{
				PushMove(ResampleSubMove,1,1,1);
				if (Qmode && (RRPrior == GammaDistributed))	{
					PushMove(RefRelativeRate,10,1,1);
					PushMove(RefRelativeRate,10,0.1,1);
					PushMove(RefRelativeRate,10,0.01,1);
					PushMove(RRAlphaM,10,1,1);
					PushMove(RRAlphaM,10,0.1,1);
					PushMove(RRAlphaM,10,0.01,1);
				}
				else	{
					if (ZipGTR == 4)	{
						PushMove(ResampleRelativeRateMove,1,0.5,5);
						PushMove(ResampleRelativeRateMove,1,0.3,10);
						PushMove(ResampleRelativeRateMove,1,0.1,10);
					}
					else	{
						PushMove(ResampleRelativeRateMove,1,1,1);
					}
					PushMove(LengthRelRate,1,1,1);
					PushMove(LengthRelRate,1,0.1,1);
					PushMove(LengthRelRate,1,0.01,1);
					if ((!Qmode) && (RRPrior == GammaDistributed))	{
						PushMove(RRAlphaM,10,1,1);
						PushMove(RRAlphaM,10,0.1,1);
						PushMove(RRAlphaM,10,0.01,1);
					}
				}
				PushMove(UpdateModeMove,1,1,1);
			}
			/*
			else if (MutMode == 3)	{
				PushMove(ResampleSubMove,1,0.1,1);
				for (int i=0; i<5; i++)	{
					PushMove(ACGT,10,0.1,1);
					PushMove(ACGT,10,0.04,1);
					PushMove(ACGT,10,0.01,1);
					PushMove(ACGT,10,0.1,2);
					PushMove(ACGT,10,0.04,2);
					PushMove(ACGT,10,0.01,2);
					PushMove(RRACGT,10,0.1,1);
					PushMove(RRACGT,10,0.01,1);
					PushMove(LengthRelRate,10,1,1);
					PushMove(LengthRelRate,10,0.1,1);
					PushMove(LengthRelRate,10,0.01,1);
					PushMove(UpdateModeMove,1,1,1);
				}
			}
			*/
		}
		if (length && WithPartial && (!ActivateClock) && ((!NH) || ((NH<4) && (!NHLengthSwitch))))	{
			if (Conjugate)	{
				PushMove(ResampleSubMove,1,1,1);
				PushMove(ResampleLengthMove,1,1,1);
			}
		}
		if (length && (!WithPartial) && ((!NH) || ((NH<4) && (!NHLengthSwitch))))	{
			if (MBL == 2)	{
				PushMove(AllBranchLength,3,1,1);
				PushMove(AllBranchLength,3,0.1,1);
			}
			else	{
				if (Conjugate)	{
					PushMove(ResampleSubMove,1,1,1);
					PushMove(ResampleLengthMove,1,1,1);
				}
				else	{
					PushMove(OneBranchLength,10,1,1);
					PushMove(OneBranchLength,10,0.1,1);
				}
			}
		}
		if (NH)	{
			if (Conjugate)	{
				for (int i=0; i<1; i++)	{
				PushMove(ResampleSubMove,1,1,1);

				if (PopEff)	{
					PushMove(PopEffSizeM,1,1,10);
					PushMove(PopEffSizeM,1,0.1,10);
					PushMove(PopEffSizeM,1,0.01,10);
					PushMove(PopAlphaM,10,10,1);
					PushMove(PopBetaM,10,10,1);
					PushMove(PopAlphaM,10,1,1);
					PushMove(PopBetaM,10,1,1);
					PushMove(PopAlphaM,10,0.1,1);
					PushMove(PopBetaM,10,0.1,1);
					if (PopEff == 3)	{
						// PushMove(PopEffModeStatM,10,1,1);
						PushMove(PopEffModeStatM,10,0.1,1);
						PushMove(PopEffModeStatM,10,0.01,1);
					}
				}
				if (!PopEff || NHPop)	{
					PushMove(NHBranchStatM,1,0.2,10);
					PushMove(NHBranchStatM,1,0.05,10);
					PushMove(NHBranchStatM,1,0.2,5);
					PushMove(NHBranchStatM,1,0.5,2);
					PushMove(NHBranchStatM,1,1,1);
					PushMove(NHBranchStatM,1,0.02,10);
					PushMove(NHBranchStatM,1,0.005,10);
					PushMove(NHBranchStatM,1,0.02,5);
					PushMove(NHBranchStatM,1,0.05,2);
					PushMove(NHBranchStatM,1,0.1,1);
				}
				/*
				PushMove(NHContBranchStat,1,10,1);
				PushMove(NHContBranchStat,1,1,1);
				PushMove(NHContBranchStat,1,0.1,1);
				*/

				/*
				PushMove(NHBranchStatTranslate,10,1,1);
				PushMove(NHBranchStatTranslate,10,0.1,1);
				PushMove(NHBranchStatTranslate,10,0.01,1);
				PushMove(NHBranchStatTranslate,10,0.001,1);
				*/
				// PushMove(NHBranchStatM,1,0.4,1);
				if (! catfix)	{
					if (ZipGTR >= 2)	{
						PushMove(ResampleStatMove,1,1,2);
						PushMove(ResampleStatMove,1,0.1,2);
					}
					else	{
						PushMove(NHModeStatM,2,1,1);
						PushMove(NHModeStatM,2,0.1,3);
						PushMove(NHModeStatM,2,1,1);
						PushMove(NHModeStatM,2,0.1,3);
						PushMove(NHModeStatM,2,1,1);
					}
				}
				// PushMove(NHBranchStatM,1,0.2,2);

				/*
				// PushMove(SwitchModeIntegral,1,1,5);
					PushMove(NHStatCenterM,10,1,2);
					PushMove(NHStatAlphaM,10,1,1);
					PushMove(NHStatCenterM,10,0.1,2);
					PushMove(NHStatCenterM,10,0.1,4);
					PushMove(NHStatCenterM,10,0.1,10);
					PushMove(NHStatAlphaM,10,0.1,1);
				PushMove(ResampleSubMove,1,1,1);
				PushMove(NHBranchStatM,5,0.2,10);
				if (ZipGTR >= 2)	{
					PushMove(ResampleStatMove,1,1,2);
					PushMove(ResampleStatMove,1,0.1,2);
				}
				else	{
					PushMove(NHModeStatM,2,1,1);
				}
				// PushMove(NHBranchStatM,1,0.1,3);
				// PushMove(SwitchModeIntegral,1,1,5);
					PushMove(NHStatCenterM,10,1,2);
					PushMove(NHStatAlphaM,10,1,1);
					PushMove(NHStatCenterM,10,0.1,2);
					PushMove(NHStatCenterM,10,0.1,4);
					PushMove(NHStatCenterM,10,0.1,10);
					PushMove(NHStatAlphaM,10,0.1,1);
				PushMove(ResampleSubMove,1,1,1);
				PushMove(NHBranchStatM,5,1,5);
				if (ZipGTR >= 2)	{
					PushMove(ResampleStatMove,1,1,2);
					PushMove(ResampleStatMove,1,0.1,2);
				}
				else	{
					PushMove(NHModeStatM,2,1,1);
				}
				// PushMove(NHBranchStatM,1,0.1,4);
				// PushMove(SwitchModeIntegral,1,1,5);
					PushMove(NHStatCenterM,10,1,2);
					PushMove(NHStatAlphaM,10,1,1);
					PushMove(NHStatCenterM,10,0.1,2);
					PushMove(NHStatCenterM,10,0.1,4);
					PushMove(NHStatCenterM,10,0.1,10);
					PushMove(NHStatAlphaM,10,0.1,1);
				PushMove(ResampleSubMove,1,1,1);
				PushMove(NHBranchStatM,5,1,2);
				if (ZipGTR >= 2)	{
					PushMove(ResampleStatMove,1,1,2);
					PushMove(ResampleStatMove,1,0.1,2);
				}
				else	{
					PushMove(NHModeStatM,2,1,1);
				}
				// PushMove(NHBranchStatM,1,0.1,5);
				// PushMove(SwitchModeIntegral,1,1,5);
					PushMove(NHStatCenterM,10,1,2);
					PushMove(NHStatAlphaM,10,1,1);
					PushMove(NHStatCenterM,10,0.1,2);
					PushMove(NHStatCenterM,10,0.1,4);
					PushMove(NHStatCenterM,10,0.1,10);
					PushMove(NHStatAlphaM,10,0.1,1);
				PushMove(ResampleSubMove,1,1,1);
				PushMove(NHBranchStatM,5,1,1);
				if (ZipGTR >= 2)	{
					PushMove(ResampleStatMove,1,1,2);
					PushMove(ResampleStatMove,1,0.1,2);
				}
				else	{
					PushMove(NHModeStatM,2,1,1);
				}
				*/
				if ((NH != 2) && (NH < 4)) {
					// PushMove(NHAllocationM,5,1,1);
					PushMove(NHAllocationM,5,1,0);
					PushMove(NHAllocationM,5,1,0);
				}
				if (! FixNmode)	{
					PushMove(SwitchModeIntegral,1,1,5);
					if (ZipGTR >= 2)	{
						PushMove(ResampleStatMove,1,1,2);
						PushMove(ResampleStatMove,1,0.1,2);
						if (ZipGTR < 3)	{
							PushMove(ResampleMaskMove,1,1,10);
							if (ZipGTRDP)	{
								PushMove(ZipSwitch,1,1,10);
								PushMove(ZipAlphaM,10,1,1);
								PushMove(ZipAlphaM,10,0.1,1);
							}
							else	{
								cerr << "zipgtrdp 0 \n";
								exit(1);
							}
						}
					}
				}
				}
			}
			else	{
				cerr << "error : nh moves only under data augmentation mcmc\n";
				exit(1);
			}
			if (NH == 1)	{
				PushMove(NHpswitchM,10,1,1);
			}
			if (NH < 4)	{
				if (! NHPrior)	{
				// if ((! NHPrior) && (! PopEff))	{
					PushMove(NHStatCenterM,10,1,2);
				}
			}
			// if (! ((NH < 4) && (NHPrior)))	{ 
			if (! NHPrior)	{
			// if ((! NHPrior) && (! PopEff))	{
				PushMove(NHStatAlphaM,10,1,1);
			}
			if (NH == 1)	{
				PushMove(NHpswitchM,10,0.1,1);
			}
			if (NH < 4)	{
				if (! NHPrior)	{
				// if ((! NHPrior) && (! PopEff))	{
					PushMove(NHStatCenterM,10,0.1,2);
				}
			}
			// if (! ((NH < 4) && (NHPrior)))	{ 
			if (! NHPrior)	{
			// if ((! NHPrior) && (! PopEff))	{
				PushMove(NHStatAlphaM,10,0.1,1);
			}
			PushMove(LeaveNHMove,1,1,1);
		}
		else if (stat)	{
			if (Conjugate)	{
				if (ModeFastCompute && (! ModePoisson))	{
					if (ZipGTR == 1)	{
						PushMove(ResampleSubMove,1,1,1);
						if (! FixNmode)	{
							PushMove(SwitchModeIntegral,1,1,5);
							// PushMove(ZipProb,1,1,1);
						}
						PushMove(ResampleStatMove,10,1,50);
						PushMove(ResampleStatMove,10,0.1,50);
					}
					else if (ZipGTR == 4)	{
						PushMove(ResampleSubMove,1,1,1);
						if (! FixNmode)	{
							PushMove(SwitchModeIntegral,1,1,5);
							// PushMove(ZipProb,1,1,1);
						}
						PushMove(ResampleStatMove,10,1,1);
						PushMove(ResampleStatMove,10,0.1,5);
						//PushMove(ResampleStatMove,5,1,1);
						//PushMove(ResampleStatMove,5,0.1,5);
					}
					else if (ZipGTR >= 2)	{
						for (int rep=0; rep<5; rep++)	{
							PushMove(ResampleSubMove,1,1,1);
							if (! FixNmode)	{
								PushMove(SwitchModeIntegral,1,1,5);
							}
							if (ZipGTR < 3)	{
								PushMove(ResampleMaskMove,1,1,10);
								if (ZipGTRDP)	{
									PushMove(ZipSwitch,1,1,10);
									PushMove(ZipAlphaM,10,1,1);
									PushMove(ZipAlphaM,10,0.1,1);
								}
								else	{
									cerr << "zipgtrdp 0 \n";
									exit(1);
								}
							}
							PushMove(ResampleStatMove,1,1,5);
							PushMove(ResampleStatMove,1,0.1,5);

							if (ZipPrior == 2)	{
								PushMove(ZipSwitch,1,1,1);
								PushMove(ZipAlphaM,10,1,1);
								PushMove(ZipAlphaM,10,0.1,1);
								PushMove(ZipProb,10,1,1);
								PushMove(ZipProb,10,0.1,1);
							}
						}
					}
					else	{
						cerr << "error in initmove: unknown zipgtr value\n";
						exit(1);
					}
				}
				else	{
					PushMove(ResampleSubMove,1,1,1);
					if (FixNmode && (!ActivateSumOverModes) && (!catfix))	{
						PushMove(SwitchModeIntegral,1,1,5);
					}
					if (! FixNmode)	{
						if (MutMode > 1)	{
							PushMove(SwitchModeIntegral,1,1,5);
						}
						else	{
							PushMove(SwitchModeIntegral,1,1,5);
						}
					}
					PushMove(ResampleStatMove,1,1,1);
					if (StatFlag)	{
						PushMove(ResampleStatMove,1,0.1,1);
					}
				}
			}
			else	{
				PushMove(UpdateLogProbMove,1,1,1);
				PushMove(ModeStationary,5,1,2);
				PushMove(ModeStationary,5,0.1,2);
				if (! FixNmode)	{
					PushMove(SwitchMode,2,1,5);
				}
			}
		}
		}
		if (rate)	{
			PushMove(ResampleSubMove,1,1,1);
			if (! fixnratemode)	{
				PushMove(RateSwitchModeIntegral,1,1,1);
			}
			PushMove(ResampleRateMove,1,1,1);
		}
		if (GeneBLMultiplier)	{
			PushMove(ResampleSubMove,1,1,1);
			PushMove(ResampleGeneBLMulMove,1,1,1);
		}
		if (MBL == 1)	{
			PushMove(ResampleSubMove,1,1,1);
			PushMove(MixBLSwitchModeIntegral,1,1,5);
			PushMove(ResampleMixBLMove,1,1,1);
			PushMove(ResampleSubMove,1,1,1);
			PushMove(MixBLSwitchModeIntegral,1,1,5);
			PushMove(ResampleMixBLMove,1,1,1);
			PushMove(ResampleSubMove,1,1,1);
			PushMove(MixBLSwitchModeIntegral,1,1,5);
			PushMove(ResampleMixBLMove,1,1,1);
			PushMove(ResampleSubMove,1,1,1);
			PushMove(MixBLSwitchModeIntegral,1,1,5);
			PushMove(ResampleMixBLMove,1,1,1);
			PushMove(ResampleSubMove,1,1,1);
			PushMove(MixBLSwitchModeIntegral,1,1,5);
			PushMove(ResampleMixBLMove,1,1,1);
		}
		if (MBL == 2)	{
			PushMove(ResampleSubMove,1,1,1);
			// for (int i=0; i<5; i++)	{
			PushMove(MixBLSwitchModeIntegral,10,1,5);
			PushMove(ResampleMixBLMove,10,1,1);
			PushMove(ResampleMixBLMove,10,0.1,1);
			PushMove(ResampleMixBLMove,10,0.1,5);
			PushMove(ResampleMixBLMove,10,0.01,5);
			PushMove(ResampleMixBLMove,10,0.01,10);
			PushMove(ResampleMixBLMove,10,0.001,10);
			// }
			/*
			PushMove(ResampleSubMove,1,1,1);
			PushMove(MixBLSwitchModeIntegral,1,1,5);
			PushMove(ResampleMixBLMove,1,1,1);
			PushMove(ResampleSubMove,1,1,1);
			PushMove(MixBLSwitchModeIntegral,1,1,5);
			PushMove(ResampleMixBLMove,1,1,1);
			PushMove(ResampleSubMove,1,1,1);
			PushMove(MixBLSwitchModeIntegral,1,1,5);
			PushMove(ResampleMixBLMove,1,1,1);
			PushMove(ResampleSubMove,1,1,1);
			PushMove(MixBLSwitchModeIntegral,1,1,5);
			PushMove(ResampleMixBLMove,1,1,1);
			*/
		}
		if (rate || GeneBLMultiplier || MBL)	{
			PushMove(UpdateLogProbMove,1,1,1);
		}
		if (GeneBLMultiplier)	{
			if (MSMode != LengthMixGamma)	{
				PushMove(LengthGammaM,10,1,1);
				PushMove(LengthGammaM,10,0.1,1);
			}
			PushMove(LengthGeneBL,10,0.1,1);
			PushMove(LengthGeneBL,10,1,1);
		}
		if (MBL)	{
			if (MSMode != LengthMixGamma)	{
				PushMove(LengthGammaM,10,1,1);
				PushMove(LengthGammaM,10,0.1,1);
			}
			if (MBL == 1)	{
				PushMove(LengthMixBL,10,0.1,1);
				PushMove(LengthMixBL,10,1,1);
			}
			PushMove(MixBLAlphaM,10,10,1);
			PushMove(MixBLAlphaM,10,1,1);
			PushMove(MixBLAlphaM,10,0.1,1);
		}
		if (catfix && (! Tempering))	{
			if (NH)	{
				PushMove(UpdateLogProbMove,1,1,0);
			}
			PushMove(SwitchMode,1,1,0);
		}
		SuperMoveTypeNumber[topo]++;
		super = SuperMoveTypeNumber[topo];
	}

	if (ActivateClock)	{
		if (! NormalApprox)	{
			for (int l=0; l<Nchain; l++)	{
				SuperNIterationArray[topo][super][l] = 1;
			}
			PushMove(ResampleSubMove,1,1,1);
			SuperMoveTypeNumber[topo]++;
			super = SuperMoveTypeNumber[topo];
		}

		if (((SUBModelSwitch == 1) && (ModeFastCompute)) || ((SUBModelSwitch == 0) && (RefFastCompute)))	{
			for (int l=0; l<Nchain; l++)	{
				SuperNIterationArray[topo][super][l] = 50;
			}
		}
		else	{
			for (int l=0; l<Nchain; l++)	{
				SuperNIterationArray[topo][super][l] = 500;
			}
		}

		PushMove(muMove,1,1,1);
		PushMove(muRho,1,1,1);

		PushMove(sigmaMove,1,1,1);
		if ((ClockModel == CIRrigid) || (FlexClockModel == CIRflex)) 	{
			PushMove(thetaMove,1,1,1);
		}
		if ((SeparateModelSwitch == 1) && (GeneBLMultiplier))	{
			PushMove(mulSigmaMove,1,1,1);
			PushMove(mulThetaMove,1,1,1);
		}
		PushMove(muMove,1,0.1,1);
		PushMove(muRho,1,0.1,1);

		PushMove(sigmaMove,1,0.1,1);
		if ((ClockModel == CIRrigid) || (FlexClockModel == CIRflex)) 	{
			PushMove(thetaMove,1,1,1);
		}
		if ((SeparateModelSwitch == 1) && (GeneBLMultiplier))	{
			PushMove(mulSigmaMove,1,0.1,1);
			PushMove(mulThetaMove,1,0.1,1);
		}
		PushMove(muMove,1,0.01,1);
		PushMove(muRho,1,0.01,1);

		PushMove(sigmaMove,1,0.01,1);
		if ((ClockModel == CIRrigid) || (FlexClockModel == CIRflex)) 	{
			PushMove(thetaMove,1,1,1);
		}
		if ((SeparateModelSwitch == 1) && (GeneBLMultiplier))	{
			PushMove(mulSigmaMove,1,0.01,1);
			PushMove(mulThetaMove,1,0.01,1);
		}
		PushMove(addRhoMove,1,1,1);
		PushMove(addRhoMove,1,0.1,1);
		PushMove(addRhoMove,1,0.01,1);
		if (ContNchar)	{
			PushMove(addRhoMove,1,10,1);
			PushMove(ContRhoCenterM,1,1,1);
			PushMove(ContRhoCenterM,1,0.1,1);
		}
		if ((BLModelSwitch == 1) || (NormalApprox))	{
			PushMove(OneBranchLength,10,10,1);
			PushMove(OneBranchLength,10,1,1);
			PushMove(OneBranchLength,10,0.1,1);
		}
		if (! KeepTimes)	{
			PushMove(oneTimeNodeMove,1,1,1);
			PushMove(oneTimeNodeMove,1,0.1,1);
			// PushMove(oneTimeNodeMove,1,0.01,1);
		}
		if (Chi == -1)	{
			if (TimePrior == BD)	{
				PushMove(chi,10,10,1);
				PushMove(chi,10,1,1);
				PushMove(chi,10,0.1,1);
				PushMove(chi2,10,10,1);
				PushMove(chi2,10,1,1);
				PushMove(chi2,10,0.1,1);
			}
			else if (TimePrior == Dirich)	{
				PushMove(chi,10,1,1);
				PushMove(chi,10,0.1,1);
			}
		}
		/*
		if (Ngene)	{
			PushMove(AddGeneRho,1,1,1);
			PushMove(AddGeneRho,1,0.1,1);
			if (BLModelSwitch == 1)	{
				PushMove(OneGeneBranch,1,10,1);
				PushMove(OneGeneBranch,1,1,1);
				PushMove(OneGeneBranch,1,0.1,1);
			}
		}
		*/
		if (NCalib)	{
			/*
			if (TimePrior == BD)	{
				PushMove(AbsBDScaleM,3,1000,1);
				PushMove(AbsBDScaleM,3,100,1);
				PushMove(AbsBDScaleM,3,10,1);
				PushMove(AbsBDScaleM,3,1,1);
			}
			*/
			PushMove(ScaleM,10,10,1);
			PushMove(ScaleM,10,1,1);
			PushMove(ScaleM,10,0.1,1);
			PushMove(ScaleM,10,0.01,1);
		}
		SuperMoveTypeNumber[topo]++;
		super = SuperMoveTypeNumber[topo];
	}

	if (! NormalApprox)	{
		for (int l=0; l<Nchain; l++)	{
			SuperNIterationArray[topo][super][l] = 1;
		}
		PushMove(UpdateLogProbMove,1,1,1);
		if (MutMode == 2)	{
			for (int i=0; i<5; i++)	{
				PushMove(KappaM,1,0.1,1);
				PushMove(KappaM,1,0.01,1);
				PushMove(ACGT,1,0.1,1);
				PushMove(ACGT,1,0.04,1);
				PushMove(ACGT,1,0.01,1);
				PushMove(ACGT,1,0.1,2);
				PushMove(ACGT,1,0.04,2);
				PushMove(ACGT,1,0.01,2);
			}
		}
		else if (MutMode == 3)	{
			/*
			PushMove(RRACGT,5,0.1,1);
			PushMove(RRACGT,5,0.01,1);
			PushMove(LengthRelRate,10,1,1);
			PushMove(LengthRelRate,10,0.1,1);
			PushMove(LengthRelRate,10,0.01,1);
			PushMove(UpdateModeMove,1,1,1);
			PushMove(ACGT,5,0.1,1);
			PushMove(ACGT,5,0.04,1);
			PushMove(ACGT,5,0.01,1);
			PushMove(ACGT,5,0.1,2);
			PushMove(ACGT,5,0.04,2);
			PushMove(ACGT,5,0.01,2);
			*/
			for (int i=0; i<5; i++)	{
				PushMove(RRACGT,1,0.1,1);
				PushMove(RRACGT,1,0.01,1);
				PushMove(LengthRelRate,2,1,1);
				PushMove(LengthRelRate,2,0.1,1);
				PushMove(LengthRelRate,2,0.01,1);
				PushMove(UpdateModeMove,1,1,1);
				PushMove(ACGT,1,0.1,1);
				PushMove(ACGT,1,0.04,1);
				PushMove(ACGT,1,0.01,1);
				PushMove(ACGT,1,0.1,2);
				PushMove(ACGT,1,0.04,2);
				PushMove(ACGT,1,0.01,2);
			}
		}
		/*
		if (MutMode == 2)	{
			PushMove(KappaM,15,0.1,1);
			PushMove(KappaM,15,0.01,1);
			PushMove(ACGT,15,0.1,1);
			PushMove(ACGT,15,0.04,1);
			PushMove(ACGT,15,0.01,1);
			PushMove(ACGT,15,0.1,2);
			PushMove(ACGT,15,0.04,2);
			PushMove(ACGT,15,0.01,2);
		}
		else if (MutMode == 3)	{
			PushMove(RRACGT,15,0.1,1);
			PushMove(RRACGT,15,0.01,1);
			PushMove(LengthRelRate,1,1,1);
			PushMove(LengthRelRate,1,0.1,1);
			PushMove(LengthRelRate,1,0.01,1);
			PushMove(UpdateModeMove,1,1,1);
			PushMove(ACGT,15,0.1,1);
			PushMove(ACGT,15,0.04,1);
			PushMove(ACGT,15,0.01,1);
			PushMove(ACGT,15,0.1,2);
			PushMove(ACGT,15,0.04,2);
			PushMove(ACGT,15,0.01,2);
		}
		else if (MutMode == 1)	{
			PushMove(RefStationaryM,30,1,1);
			PushMove(RefStationaryM,30,0.1,1);
			PushMove(RefStationaryM,30,0.1,3);
			PushMove(RefStationaryM,30,0.1,3);
			PushMove(RefStationaryM,30,0.01,3);
			PushMove(RefStationaryM,30,0.1,10);
			PushMove(RefStationaryM,30,0.01,10);
			PushMove(UpdateModeMove,1,1,1);
		}
		*/
		/*
		if ((MutMode == 1) && (GWf==1))	{
			PushMove(GWF,3,1,1);
			PushMove(GWF,3,0.1,1);
			PushMove(GWF,3,0.01,1);
		}
		*/
		int gamma = ! FixGamma;
		if (gamma && (HeteroMode == Homo))	{
			if (GammaNcat && (! Tempering))	{
				PushMove(ResampleSubMove,1,1,1);
				PushMove(GammaIntegral,10,1,1);
				PushMove(GammaIntegral,10,0.1,1);
				PushMove(GammaIntegral,10,0.01,1);
				if (! FixPconst)	{
					PushMove(Pconst,1,1,1);
				}
				
			}
			else	{
				if (Tempering)	{
					PushMove(Gamma,4,1,1);
					PushMove(Gamma,4,0.1,1);
					PushMove(Gamma,4,0.01,1);
				}
				else	{
					PushMove(Gamma,10,1,1);
					PushMove(Gamma,10,0.1,1);
					PushMove(Gamma,10,0.01,1);
				}
				if (NRateModeMax == Nsite + 5)	{
					PushMove(LengthRate,10,0.1,1);
					PushMove(LengthRate,10,1,1);
				}
			}
		}
		if (HeteroMode == Covarion)	{
			if (ExternalCov)	{
				PushMove(ResampleSubMove,1,1,1);
				PushMove(GammaIntegral,10,1,1);
				PushMove(GammaIntegral,10,0.1,1);
				PushMove(GammaIntegral,10,0.01,1);
				PushMove(UpdateLogProbMove,1,1,1);
			}
			else	{
				PushMove(UpdateLogProbMove,1,1,1);
				PushMove(Gamma,3,1,1);
				PushMove(Gamma,3,0.1,1);
				PushMove(Gamma,3,0.01,1);
			}
			PushMove(Xi,3,1,1);
			PushMove(InvProb,3,1,1);
			PushMove(Xi,3,0.1,1);
			PushMove(InvProb,3,0.1,1);
			PushMove(Xi,3,0.01,1);
			PushMove(InvProb,3,0.01,1);
		}
		if (ActivateSumOverRateModes && (! Tempering))	{
			PushMove(UpdateSumRateModeMove,1,1,1);
			if (! FixPconst)	{
				PushMove(GammaPconst,1,1,1);
				PushMove(GammaPconst,1,0.1,1);
			}
		}
		if (ActivateSumOverModes && (! Tempering))	{
			PushMove(ResampleModeWeightMove,1,1,1);
			// PushMove(UpdateSumModeMove,1,1,1);
			// PushMove(SumModeStationary,1,1,2);
			// PushMove(SumModeStationary,1,0.1,2);
		}
		if (!ActivateSumOverRateModes && ! ActivateSumOverModes)	{
			PushMove(UpdateLogProbMove,1,1,1);
		}
		SuperMoveTypeNumber[topo]++;
	}
	TopoMoveTypeNumber++;
}


void MCParameters::PushMove(const MoveType movetype, const int Nrep, const double delta, const int N)	{

	int topo = TopoMoveTypeNumber;
	int super = SuperMoveTypeNumber[topo];
	int move = MoveTypeNumber[topo][super];
	
	MoveTypeArray[topo][super][move] = movetype;
	for (int l=0; l<Nchain; l++)	{
		NIterationArray[topo][super][move][l] = Nrep;
	}
	for (int l=0; l<Nchain; l++)	{
		deltaArray[topo][super][move][l] = delta;
	}
	for (int l=0; l<Nchain; l++)	{
		NArray[topo][super][move][l] = N;
	}
	
	MoveTypeNumber[topo][super]++;

}

void MCParameters::DesactivateLengthMoves()	{
	for (int topo=0; topo<TopoMoveTypeNumber; topo++)	{
		for (int super=0; super<SuperMoveTypeNumber[topo]; super++)	{
			for (int move=0; move<MoveTypeNumber[topo][super]; move++)	{
				
				MoveType movetype = MoveTypeArray[topo][super][move];
				if (movetype == ResampleLengthMove)	{
					NIterationArray[topo][super][move] = 0;
					if (! move)	{
						cerr << "error in mc::desactivate\n";
						exit(1);
					}
					NIterationArray[topo][super][move-1] = 0;
				}
			}
		}
	}
}

void MCParameters::MoveAll(int conjugate)	{

	int gamma = ! FixGamma;
	int meanlength = ! FixMeanLength;
	int statcenter = ! FixStatCenter;
	int rr = ! FixRR;
	int length = ! FixLength;
	int stat = ! FixStat;
	int rate = ! FixRate; 

	PhyloBayes* copy = currentState[0];

	int Niter = 10;
	double gammasuccess = 0;
	double meanlengthsuccess = 0;
	if (gamma)	{
		if (! copy)	{
			cerr << "error in MCParam::Move()\n";
			exit(1);
		}
		for (int i=0; i<Niter; i++)	{
			gammasuccess += GetNextState()->Move(Gamma, 1, 1, copy);
			gammasuccess += GetNextState()->Move(Gamma, 0.1, 1, copy);
		}
	}

	if (meanlength)	{
		if (! copy)	{
			cerr << "error in MCParam::Move()\n";
			exit(1);
		}
		for (int i=0; i<Niter; i++){
			meanlengthsuccess += GetNextState()->Move(MeanLengthM, 1, 1, copy);
			meanlengthsuccess += GetNextState()->Move(MeanLengthM, 0.1, 1, copy);
		}
	}
	if (statcenter)	{
		if (! copy)	{
			cerr << "error in MCParam::Move()\n";
			exit(1);
		}
		for (int i=0; i<Niter; i++){
			GetNextState()->Move(ModeStatCenterM, 10000, 5, copy);
			GetNextState()->Move(ModeStatAlphaM, 1, 1, copy);
			GetNextState()->Move(ModeStatCenterM, 100000, 5, copy);
			GetNextState()->Move(ModeStatAlphaM, 0.1, 1, copy);
			GetNextState()->Move(ModeStatCenterM, 1000000, 5, copy);
			GetNextState()->Move(ModeStatAlphaM, 0.01, 1, copy);
		}

	}
	if (rr)	{
		if (conjugate)	{
			GetNextState()->ResampleSub();
			GetNextState()->ResampleRelativeRate(copy);
		}
		else	{
			GetNextState()->Move(ModeRelativeRate, 1, 1, copy);
		}
		for (int i=0; i<Niter; i++)	{
			GetNextState()->Move(LengthRelRate, 1, 1, copy);
			GetNextState()->Move(LengthRelRate, 0.1, 1, copy);
			GetNextState()->Move(RateRelRate, 1, 1, copy);
			GetNextState()->Move(RateRelRate, 0.1, 1, copy);
			GetNextState()->Move(LengthRate, 1, 1, copy);
			GetNextState()->Move(LengthRate, 0.1, 1, copy);
		}
		
	}
	if (length)	{
		if (conjugate)	{
			GetNextState()->ResampleSub();
			GetNextState()->ResampleLength(copy);
		}
		else	{
			GetNextState()->Move(AllBranchLength, 1, 1, copy);
			GetNextState()->Move(AllBranchLength, 0.1, 1, copy);
		}
		
	}
	if (rate)	{
		MoveRate();
	}
	if (stat)	{
		MoveStat();
	}
	Update();
}


//---------------------------------------------------------------------------------
//		 InitFromFile(string filename)
//---------------------------------------------------------------------------------

void	MCParameters::InitFromFile(string filename)	{

	ifstream is(filename.c_str());
	string temp = "";

	try	{

		while (! is.eof())	{

			is >> temp;
			cerr << temp << '\n';

			if (temp == "ModePoisson")	{
				is >> temp;
				if (temp == "Yes")	{
					ModePoisson = Yes;
				}
				else if (temp == "No")	{
					ModePoisson = No;
				}
				else	{
					cerr << "error after ModePoisson\n";
					cerr.flush();
					throw;
				}
			}

			else if (temp == "RefPoisson")	{
				is >> temp;
				if (temp == "Yes")	{
					RefPoisson = Yes;
				}
				else if (temp == "No")	{
					RefPoisson = No;
				}
				else	{
					cerr << "error after RefPoisson\n";
					cerr.flush();
					throw;
				}
			}

			else if (temp == "ModeFastCompute")	{
				is >> temp;
				if (temp == "Yes")	{
					ModeFastCompute = Yes;
				}
				else if (temp == "No")	{
					ModeFastCompute = No;
				}
				else	{
					cerr << "error after ModeFastCompute\n";
					cerr.flush();
					throw;
				}
			}

			else if (temp == "RefFastCompute")	{
				is >> temp;
				if (temp == "Yes")	{
					RefFastCompute = Yes;
				}
				else if (temp == "No")	{
					RefFastCompute = No;
				}
				else	{
					cerr << "error after RefFastCompute\n";
					throw;
				}
			}

			else if (temp == "Calibrations")	{
				ReadCalibration(is);
			}

			else if (temp == "MCMCMC")	{
				is >> Nchain;
				double step;
				if (Nchain * step > 1)	{
					cerr << "error after MCMCMC\n";
					exit(1);
				}
				is >> step;
				for (int l=0; l<Nchain; l++)	{
					Beta[l] = 1.0 - l *step;
					// Beta[l] = 1.0 / (1 + l* step);
				}
			}

			else if (temp == "Parallel")	{
				is >> temp;
				if (temp == "Yes")	{
					Parallel = 1;
				}
				else if (temp == "No")	{
					Parallel = 0;
				}
				else	{
					cerr << "error : " << temp << " after Parallel\n";
					exit(1);
				}
			}

			else if (temp == "SwapMode")	{
				is >> temp;
				if (temp == "Wave")	{
					SwapMode = 1;
				}
				else if (temp == "Point")	{
					SwapMode = 0;
				}
				else	{
					cerr << "error after SwapMode\n";
					exit(1);
				}
			}


			else if (temp == "Qmode")	{
				string temp;
				is >> temp;
				if (temp == "Yes")	{
					Qmode = Yes;
				}
				else if (temp == "No")	{
					Qmode = No;
				}
				else	{
					cerr << "error after Qmode\n";
					exit(1);
				}
			}
			
			else if (temp == "SeparateRhoPrior")	{
				is >> temp;
				if (temp == "Yes")	{
					SeparateRhoPrior = Yes;
				}
				else if (temp == "No")	{
					SeparateRhoPrior = No;
				}
				else	{
					cerr << "error after SeparateRhoPrior\n";
					exit(1);
				}
			}

			else if (temp == "NmodeMax")	{
				int temp;
				is >> temp;
				NmodeMax = temp;
			}
			
			else if (temp == "NRateModeMax")	{
				int temp;
				is >> temp;
				NRateModeMax = temp;
			}
			
			else if (temp == "SumOverRateModes")	{
				string tmp;
				is >> tmp;
				if (tmp == "Yes")	{
					ActivateSumOverRateModes = Yes;
				}
				else if (temp == "No")	{
					ActivateSumOverRateModes = No;
				}
				else {
					cerr << "error after SumOverRateModes\n";
					exit(1);
				}
			}
			else if (temp == "DiscreteGamma")	{
				is >> GammaNcat;
				ActivateSumOverRateModes = Yes;
				RASModelSwitch = 0;
			}	
			
			else if (temp == "DCM")	{
				is >> temp;
				if (temp == "Yes")	{
					DCM = Yes;
				}
				else if (temp == "No")	{
					DCM = No;
				}
				else 	{
					cerr << "error after DCM\n";
				}
			}

			else if (temp == "SumOverModes")	{
				string tmp;
				is >> tmp;
				if (tmp == "Yes")	{
					ActivateSumOverModes = Yes;
				}
				else if (tmp == "No")	{
					ActivateSumOverModes = No;
				}
				else {
					cerr << "error after SumOverModes\n";
					exit(1);
				}
			}

			else if (temp == "GeneBLMode")	{
				string tmp;
				is >> tmp;
				if (tmp == "Yes")	{
					GeneBLMode = Yes;
				}
				else if (tmp == "No")	{
					GeneBLMode = No;
				}
				else if (tmp == "Multipliers")	{
					GeneBLMode = Yes;
					GeneBLMultiplier = Yes;
				}
				else {
					cerr << "error after GeneBLMode\n";
					exit(1);
				}
			}

			else if (temp == "GeneGammaMode")	{
				string tmp;
				is >> tmp;
				if (tmp == "Yes")	{
					GeneGammaMode = Yes;
				}
				else if (tmp == "No")	{
					GeneGammaMode = No;
				}
				else {
					cerr << "error after GeneGammaMode\n";
					exit(1);
				}
			}

			else if (temp == "GeneStationaryMode")	{
				string tmp;
				is >> tmp;
				if (tmp == "Yes")	{
					GeneStationaryMode = Yes;
				}
				else if (tmp == "No")	{
					GeneStationaryMode = No;
				}
				else {
					cerr << "error after GeneStationaryMode\n";
					exit(1);
				}
			}

			else if (temp == "GeneRateMode")	{
				string tmp;
				is >> tmp;
				if (tmp == "Yes")	{
					GeneRateMode = Yes;
				}
				else if (tmp == "No")	{
					GeneRateMode = No;
				}
				else {
					cerr << "error after GeneRateMode\n";
					exit(1);
				}
			}

			else if (temp == "BypassSub")	{
				Rao = 0;
			}
			
			else if (temp == "SaveAllChains")	{
				is >> temp;
				if (temp == "Yes")	{
					SaveAllChains = 1;
				}
				else if (temp == "No")	{
					SaveAllChains = 0;
				}
				else	{
					cerr << "error : " << temp << " after SaveAllChains\n";
					exit(1);
				}
			}

			else if (temp == "QuasiStatic")	{
				is >> InitBeta;
				is >> FinalBeta;
				is >> BetaStep;
				is >> BurnIn;
			}

			else if (temp == "Burnin")	{
				is >> BurnIn;
			}

			else if (temp == "Beta")	{
				is >> Beta0;
			}

			else if (temp == "SwapFreq")	{
				is >> SwapFreq;
				if ((SwapFreq < 0) || (SwapFreq > 1))	{
					cerr << "error : SwapFreq = : " << SwapFreq << '\n';
					exit(1);
				}
			}

			else if (temp == "MSMode")	{
				is >> temp;
				if (temp == "None")	{
					MSMode = None;
				}
				else if (temp == "FlexClock")	{
					MSMode = FlexClock;
				}
				else if (temp == "BL")	{
					MSMode = BL;
				}
				else if (temp == "AutoCorrel")	{
					MSMode = Autocorrel;
				}
				else if (temp == "ClockPrior")	{
					MSMode = ClockPrior;
				}
				else if (temp == "Thermo")	{
					MSMode = Thermo;
				}
				else if (temp == "SA")	{
					MSMode = SA;
				}
				else if (temp == "RAS")	{
					MSMode = RAS;
				}
				else if (temp == "SUB")	{
					MSMode = SUB;
				}
				else if (temp == "Separate")	{
					MSMode = Separate;
				}
				else if (temp == "Hetero")	{
					MSMode = Hetero;
				}
				else if (temp == "CATAlpha")	{
					MSMode = CATAlpha;
				}
				else	{
					cerr << "error after MSMode\n";
					throw;
				}
			}
	
			else if (temp == "ArcThermo")	{
				ArcThermo = 1;
				is >> alpha1 >> alpha2;
			}

			else if (temp == "Cov")	{
				is >> temp;
				ReadCov(temp);
			}

			else if (temp == "SepCov")	{
				is >> temp;
				ReadSepCov(temp);
			}

			else if (temp == "OutGroup")	{
				is >> temp;
				ReadOutGroup(temp);
			}

			else if (temp == "NormalConcMode")	{
				is >> temp;
				if (temp == "FromConc")	{
					NormalConcMode = 1;
				}
				else if (temp == "FromSep")	{
					NormalConcMode = 2;
				}
			}

			else if (temp == "DataFile")	{
				is >> temp;
				ReadDataFromFile(temp);
   			}

			else if (temp == "Recoding")	{
				is >> temp;
				if (temp == "NoRecoding")	{
					Recoding = NoRecoding;
				}
				else if (temp == "dayhoff6")	{
					Recoding = Dayhoff6;
				}
				else if (temp == "dayhoff4")	{
					Recoding = Dayhoff4;
				}
				else if (temp == "HP")	{
					Recoding = HP;
				}
				else if (temp == "Custom")	{
					Recoding = Custom;
					is >> RecodingFile;
				}
				LoadRecoding();
			}

			else if (temp == "SaveEvery")	{
				is >> SaveEvery;
			}

			else if (temp == "StopAfter")	{
				is >> StopAfter;
			}

			else if (temp == "TopoMoveType")	{
				ReadTopoMoveType(is);
			}
			
			else if (temp == "SavePartialLogLikelihoods")	{
				is >> temp;
				if (temp == "Yes")	{
					SavePartialLogLikelihoods = Yes;
				}
				else if (temp == "No")	{
					SavePartialLogLikelihoods = No;
				}
				else	{
					cerr << "error after SavePartialLogLikelihoods\n";
					throw;
				}
			}

			else if (temp == "EliminateConstantPositions")	{
				is >> temp;
				if (temp == "Yes")	{
					DeleteConstant = Yes;
				}
				else if (temp == "No")	{
					DeleteConstant = No;
				}
				else	{
					cerr << "error after EliminateConstantPositions\n";
					throw;
				}
			}

			else if (temp == "Normalise")	{
				is >> temp;
				if (temp == "Yes")	{
					Normalise = Yes;
				}
				else if (temp == "No")	{
					Normalise = No;
				}
				else	{
					cerr << "error after Normalise\n";
					throw;
				}
			}

			else if (temp == "SaveStat")	{
				is >> temp;
				if (temp == "Yes")	{
					SaveStat = 1;
				}
				else if (temp == "No")	{
					SaveStat = 0;
				}
				else	{
					cerr << "error after SaveStat\n";
					throw;
				}
			}
	
			else if (temp == "AlphaPrior")	{
				is >> temp;
				if (temp == "Flat")	{
					AlphaPrior = Flat;
				}
				else if (temp == "Exponential")	{
					AlphaPrior = Exponential;
				}
				else	{
					cerr << "error after AlphaPrior\n";
					throw;
				}
			}

			else if (temp == "ModePrior")	{
				is >> temp;
				if (temp == "Flat")	{
					ModePrior = Flat;
				}
				else if (temp == "DirichletProcess")	{
					ModePrior = DirichletProcess;
				}
				else	{
					cerr << "error after ModePrior : " << temp << '\n';
					exit(1);
				}
			}

			else if (temp == "RateModePrior")	{
				is >> temp;
				if (temp == "Flat")	{
					RateModePrior = Flat;
				}
				else if (temp == "DirichletProcess")	{
					RateModePrior = DirichletProcess;
				}
				else	{
					cerr << "error after ModePrior : " << temp << '\n';
					exit(1);
				}
			}

			else if (temp == "LengthPrior")	{
				is >> temp;
				if (temp == "Flat")	{
					LengthPrior = Flat;
				}
				else if (temp == "Exponential")	{
					LengthPrior = Exponential;
				}
				else if (temp == "PowerLaw")	{
					LengthPrior = PowerLaw;
				}
				else if (temp == "Gamma")	{
					LengthPrior = GammaDistributed;
				}
				else	{
					cerr << "error in choice of Length prior\n";
					throw;
				}
			}

			else if (temp == "RatePrior")	{
				is >> temp;
				if (temp == "GammaInv")	{
					RatePrior = GammaInv;
				}
				else if (temp == "Dirichlet")	{
					RatePrior = Dirichlet;
				}
				else if (temp == "Flat")	{
					RatePrior = Flat;
					is >> RateMin >> RateMax; 
				}
				else if (temp == "Dirac")	{
					RatePrior = Dirac;
				}
				else	{
					cerr << "error in choice of Rate prior\n";
					throw;
				}
			}

			else if (temp == "ModeStatPrior")	{
				is >> temp;
				if (temp == "Flat")	{
					ModeStatPrior = Flat;
				}
				else if (temp == "MultiGamma")	{
					ModeStatPrior = MultiGamma;
				}
				else	{
					cerr << "error in choive of Mode Stat Prior\n";
					throw;
				}
			}

			else if (temp == "TimePrior")	{
				is >> temp;
				if (temp == "Unif")	{
					TimePrior = Unif;
				}
				else if (temp == "BD")	{
					TimePrior = BD;
				}
				else if (temp == "Dir")	{
					TimePrior = Dirich;
				}
				else	{
					cerr << "error after TimePrior\n";
					exit(1);
				}
			}

			else if (temp == "RASModelSwitch")	{
				is >> RASModelSwitch;
			}

			else if (temp == "SUBModelSwitch")	{
				is >> SUBModelSwitch;
			}

			else if (temp == "SeparateModelSwitch")	{
				is >> SeparateModelSwitch;
			}

			else if (temp == "FlexClockModelSwitch")	{
				is >> FlexClockModelSwitch;
			}

			else if (temp == "BLModelSwitch")	{
				is >> BLModelSwitch;
			}

			else if (temp == "ClockPriorModelSwitch")	{
				is >> ClockPriorModelSwitch;
			}

			else if (temp == "AutoCorrelModelSwitch")	{
				is >> AutoCorrelModelSwitch;
			}

			else if (temp == "ActivateClock")	{
				is >> temp;
				if (temp == "Yes")	{
					ActivateClock = Yes;
				}
				else if (temp == "No")	{
					ActivateClock = No;
				}
				else	{
					cerr << "error after activate clock\n";
					exit(1);
				}
			}

			else if (temp == "FlexClockModel")	{
				ActivateClock = Yes;
				is >> temp;
				if (temp == "WhiteNoise")	{
					FlexClockModel = WhiteNoise;
				}
				else if (temp == "UGam")	{
					FlexClockModel = UGam;
				}
				else if (temp == "AGam")	{
					FlexClockModel = AGam;
				}
				else	{
					cerr << "error after ClockModel\n";
					exit(1);
				}
			}

			else if (temp == "ClockModel")	{
				ActivateClock = Yes;
				is >> temp;
				if (temp == "CIR")	{
					ClockModel = CIRrigid;
				}
				else if (temp == "KT")	{
					ClockModel = KT;
				}
				else if (temp == "UndebiasedKT")	{
					ClockModel = UndebiasedKT;
				}
				else if (temp == "Strict")	{
					ClockModel = Strict;
				}
				else if (temp == "DCIR")	{
					ClockModel = CIRrigid;
					is >> Nrho;
				}
				else	{
					cerr << "error after ClockModel\n";
					exit(1);
				}
			}

			else if (temp == "HeteroModelSwitch")	{
				is >> HeteroModelSwitch;
			}

			else if (temp == "HeteroMode")	{
				is >> temp;
				if (temp == "Homo")	{
					HeteroMode = Homo;
				}
				else if (temp == "Covarion")	{
					HeteroMode = Covarion;
				}
				else	{
					cerr << "error in choice of Heterotachy mode\n";
					throw;
				}
			}

			else if (temp == "XiMax")	{
				is >> XiMax;
			}

			else if (temp == "XiMin")	{
				is >> XiMin;
			}

			else if (temp == "AlphaMax")	{
				is >> AlphaMax;
			}

			else if (temp == "AlphaMin")	{
				is >> AlphaMin;
			}

			else if (temp == "StatMax")	{
				is >> StatMax;
			}

			else if (temp == "StatMin")	{
				is >> StatMin;
			}

			else if (temp == "LengthMin")	{
				is >> LengthMin;
			}

			else if (temp == "LengthMax")	{
				is >> LengthMax;
			}

			else if (temp == "Partition")	{
				ReadPartition(is);
			}

			else if (temp == "InitState")	{
				if (! DataOK)	{
					cerr << "should specify data file before initstate\n";
					throw;
				}
				if (! Nchain)	{
					cerr << "should specify number of chains before initstate\n";
					throw;
				}
				RegisterWithData();
				if (NmodeMax == -1)	{
					NmodeMax =  Nsite + 20;
				}
				if (NRateModeMax == -1)	{
					NRateModeMax = Nsite + 20;
				}
				currentState = new PhyloBayes*[Nchain];
				nextState = new PhyloBayes*[Nchain];
				for (int chain = 0; chain < Nchain; chain ++)	{
					currentState[chain] = new PhyloBayes(this);
					nextState[chain] = new PhyloBayes(this);
				}
				if (verbose)	{
					cerr << "init phylobayes from stream\n";
					cerr.flush();
				}
				currentState[0]->InitFromStream(is);
				if (verbose)	{
					cerr << "cloning state to coupled chains\n";
					cerr.flush();
				}
				for (int chain = 1; chain < Nchain; chain ++)	{
					currentState[chain]->Clone(currentState[0]);
					cerr << "MAKE RANDOM TREES FOR ALL CHAINS\n";
					currentState[chain]->MakeRandomTree();
					currentState[chain]->DrawLengths();
				}
			}

			else	{
				if (temp != ""){
					cerr << "error : cannot read \"" << temp << "\"\n";
					throw;
				}
			}
			char c;
			while ( ((c = is.peek()) != EOF) && ( (c == ' ') || (c == '\n') || (c == '\t')) )	{
				is.get();
			}
		}
	}
	catch(...)	{
		cerr << "InitFile processing failed\n";
		exit(1);
	}
}


//---------------------------------------------------------------------------------
//		 LoadRecoding()
//---------------------------------------------------------------------------------

void MCParameters::LoadRecoding()	{

	if (Recoding == NoRecoding)	{
	}
	else if (Recoding == Dayhoff6)	{
		RecodingNstate = 6;
		RecodingAlphabet = new char[6];
		for (int i=0; i<6; i++)	{
			RecodingAlphabet[i] = Dayhoff6Alphabet[i];
		}
		RecodingTable = new int[20];
		for (int i=0; i<20; i++)	{
			RecodingTable[i] = Dayhoff6Table[i];
		}
	}
	else if (Recoding == Dayhoff4)	{
		RecodingNstate = 4;
		RecodingAlphabet = new char[4];
		for (int i=0; i<4; i++)	{
			RecodingAlphabet[i] = Dayhoff4Alphabet[i];
		}
		RecodingTable = new int[20];
		for (int i=0; i<20; i++)	{
			RecodingTable[i] = Dayhoff4Table[i];
		}
	}
	else if (Recoding == HP)	{
		RecodingNstate = 2;
		RecodingAlphabet = new char[2];
		for (int i=0; i<2; i++)	{
			RecodingAlphabet[i] = HPAlphabet[i];
		}
		RecodingTable = new int[20];
		for (int i=0; i<20; i++)	{
			RecodingTable[i] = HPTable[i];
		}
	}
	else if (Recoding == Custom)	{
		if (RecodingFile == "")	{
			cerr << "error: cannot find recoding file\n";
		       	exit(1);
		}
		ifstream is(RecodingFile.c_str());
		RecodingAlphabet = new char[RecodingNstate];
		RecodingTable = new int[Nstate];
		char a[Nstate];
		char b[Nstate];
		for (int i=0; i<Nstate; i++)	{
			is >> a[i] >> b[i];
		}
		// check that all entries correspond to bona fide states
		for (int i=0; i<Nstate; i++)	{
			char c = a[i];
			int k = 0;
			while ((k<Nstate) && (c!=Alphabet[k])) k++;
			if (k==Nstate)	{
				cerr << "error in recoding table : " << c << "is not a recognised state\n";
				exit(1);
			}
		}
		// check that all states are present
		for (int i=0; i<Nstate; i++)	{
			char c = Alphabet[i];
			int k = 0;
			while ((k<Nstate) && (c!=a[k])) k++;
			if (k==Nstate)	{
				cerr << "error in recoding table : did not found state  " << c << '\n';
				exit(1);
			}
		}
		// compute the size of the final alphabet
		int i = 0;
		while ((i<Nstate) && ((b[i]=='-') || (b[i]=='?') || (b[i]=='X')))	{
			i++;
		}
		if (i==Nstate)	{
			cerr << "?? : recoding alphabet of size 0\n";
			exit(1);
		}
		int size = 1;
		int imin = i;
		i++;
		while (i<Nstate)	{
			char c = b[i];
			if ((c!='-') && (c!='?') && (c!='X'))	{
				int k = 0;
				while ((k<i) && (c!=b[k])) k++;
				if (k==i) size++;
			}
			i++;
		}
		// compute final alphabet
		RecodingNstate = size;
		RecodingAlphabet = new char[RecodingNstate];
		RecodingAlphabet[0] = b[imin];
		size = 1;
		for (int i=imin+1; i<Nstate; i++)	{
			char c = b[i];
			if ((c!='-') && (c!='?') && (c!='X'))	{
				int k = 0;
				while ((k<i) && (c!=b[k])) k++;
				if (k==i)	{
					RecodingAlphabet[size] = c;
					size++;
				}
			}
		}
		/*
		for (int i=0; i<RecodingNstate; i++)	{
			cerr << RecodingAlphabet[i] << ' ';
		}
		cerr << '\n';
		*/
		// load translation table
		for (int i=0; i<Nstate; i++)	{
			char c = Alphabet[i];
			int k = 0;
			while ((k<Nstate) && (c!=a[k])) k++;
			if (k==Nstate)	{
				cerr << "error in recoding table : did not found state " << c << '\n';
				exit(1);
			}
			char d = b[k];
			if ((d!='-') && (d!='?') && (d!='X'))	{
				int l = 0;
				while ((l<RecodingNstate) && (d!=RecodingAlphabet[l])) l++;
				if (l==RecodingNstate)	{
					cerr << "error when loading translation table\n";
					cerr << "does not recognise recoded state : " << l << '\n';
					exit(1);
				}
				RecodingTable[k] = l;
			}
			else	{
				RecodingTable[k] = -1;
			}
		}
	}
}

//---------------------------------------------------------------------------------
//		 RecodeData()
//---------------------------------------------------------------------------------

void MCParameters::RecodeData()	{

	RawData = new int*[Ntaxa];	
	for (int j=0; j<Ntaxa; j++)	{
		RawData[j] = new int[Nsite];
		for (int i=0; i<Nsite; i++)	{
			if (Data[j][i] == unknown)	{
				RawData[j][i] = unknown;
			}
			else	{
				RawData[j][i] = RecodingTable[Data[j][i]];
			}
		}
	}
	RawNstate = Nstate;
	RawAlphabet = Alphabet;
	/*
	int* tmp = 0;
	for (int j=0; j<Ntaxa; j++)	{
		tmp = RawData[j];
		RawData[j] = Data[j];
		Data[j] = tmp;
	}
	*/
	int** tmp = RawData;
	RawData = Data;
	Data = tmp;
	Nstate = RecodingNstate;
	Alphabet = RecodingAlphabet;
	cerr << "data recoded. new alphabet size : " << Nstate << '\n';
	//EliminateUnknownColumns();
}


//---------------------------------------------------------------------------------
//		 ReadCustomRR(string filename)
//---------------------------------------------------------------------------------

void MCParameters::ReadCustomRR(string filename)	{
	
	int Nrr = Nstate * (Nstate - 1)/2;
	CustomRR = new double[Nrr];
	for (int i=0; i<Nrr; i++)	{
		CustomRR[i] = 0;
	}
	double* rr = CustomRR;

	if ((filename == "CGR10") || (filename == "cgr10"))	{
		double total = 0;
		for (int i=0; i<Nrr; i++)	{
			rr[i]= CG10RR[i];
			total += rr[i];
		}
		for (int i=0; i<Nrr; i++)	{
			rr[i] /= total / Nrr;
		}
	}
	else if ((filename == "CGR20") || (filename == "cgr20"))	{
		double total = 0;
		for (int i=0; i<Nrr; i++)	{
			rr[i]= CG20RR[i];
			total += rr[i];
		}
		for (int i=0; i<Nrr; i++)	{
			rr[i] /= total / Nrr;
		}
	}
	else if ((filename == "CGR30") || (filename == "cgr30"))	{
		double total = 0;
		for (int i=0; i<Nrr; i++)	{
			rr[i]= CG30RR[i];
			total += rr[i];
		}
		for (int i=0; i<Nrr; i++)	{
			rr[i] /= total / Nrr;
		}
	}
	else if ((filename == "CGR40") || (filename == "cgr40"))	{
		double total = 0;
		for (int i=0; i<Nrr; i++)	{
			rr[i]= CG40RR[i];
			total += rr[i];
		}
		for (int i=0; i<Nrr; i++)	{
			rr[i] /= total / Nrr;
		}
	}
	else if ((filename == "CGR50") || (filename == "cgr50"))	{
		double total = 0;
		for (int i=0; i<Nrr; i++)	{
			rr[i]= CG50RR[i];
			total += rr[i];
		}
		for (int i=0; i<Nrr; i++)	{
			rr[i] /= total / Nrr;
		}
	}
	else if ((filename == "CGR60") || (filename == "cgr60"))	{
		double total = 0;
		for (int i=0; i<Nrr; i++)	{
			rr[i]= CG60RR[i];
			total += rr[i];
		}
		for (int i=0; i<Nrr; i++)	{
			rr[i] /= total / Nrr;
		}
	}
	else	{

	ifstream is(filename.c_str());
	if (!is)	{
		cerr << "error when reading custom rr: cannot open file\n";
		exit(1);
	}
	// read alphabet
	int permut[Nstate];
	for (int k=0; k<Nstate; k++)	{
		char c;
		is >> c;
		int l=0;
		while ((l<Nstate) && (c != Alphabet[l])) l++;
		if (l == Nstate)	{
			cerr << "error when reading empirical mixture in " << filename << ": does not recognise letter " << c << '\n';
			exit(1); 
		}
		permut[k] = l;
	}
	for (int k=0; k<Nstate-1; k++)	{
		for (int l=k+1; l<Nstate; l++)	{
			is >> CustomRR[Random::rrindex(permut[k],permut[l],Nstate)];
		}
	}
	}
	/*
	for (int k=0; k<Nrr; k++)	{
		if (! is.eof())	{
			char a,b;
			is >> a >> b;
			int i = 0;
			while ((i<Nstate) && (a != Alphabet[i])) i++;
			if (i == Nstate)	{
				cerr << "error when reading custom rr: does not recognise state " << a << '\n';
				exit(1);
			}
			int j = 0;
			while ((j<Nstate) && (b != Alphabet[j])) j++;
			if (j == Nstate)	{
				cerr << "error when reading custom rr: does not recognise state " << j << '\n';
				exit(1);
			}
			if (i == j)	{
				cerr << "error : not relative rate should be specified between identical amino-acids\n";
			}
			int l = (2 * Nstate - i - 1) * i / 2 + j - i - 1;
			if (i > j)	{
				l = (2 * Nstate - j - 1) * i / 2 + i - j - 1;
			}
			is >> CustomRR[l];
		}
	}
	*/
}

	
//---------------------------------------------------------------------------------
//		 ReadMixture(string filename)
//---------------------------------------------------------------------------------

void MCParameters::ReadMixture(string filename)	{

	if (Qmode)	{
		ReadMatFix(filename);
	}
	else	{
		ReadStatFix(filename);
	}
}
	
int MCParameters::ReadCatConstraints(string filename)	{

	ifstream is(filename.c_str());
	is >> Ncat;
	StatFlag = new int*[Ncat];
	cerr << Nstate << '\n';
	for (int i=0; i<Ncat; i++)	{
		StatFlag[i] = new int[Nstate];
		for (int j=0; j<Nstate; j++)	{
			StatFlag[i][j] = 0;
		}
		string tmp;
		is >> tmp;
		for (unsigned int j=0; j<tmp.length(); j++)	{
			int k = 0;
			while ((k<Nstate) && (tmp[j] != Alphabet[k])) k++;
			if (k == Nstate)	{
				cerr << "error in read cat constraints: " << tmp[j] << '\n';
				exit(1);
			}
			StatFlag[i][k] = 1;
		}
	}
	return Ncat;
}

//---------------------------------------------------------------------------------
//		 ReadStatFix(istream& is)
//---------------------------------------------------------------------------------

void MCParameters::ReadStatFix(string filename)	{
	if (filename == "ef")	{
		Ncat = EFN;
		statfix = new double*[Ncat];
		weight = new double[Ncat];
		for (int i=0; i<Ncat; i++)	{
			statfix[i] = new double[Nstate];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] = EFStatFix[i][k];
				if (statfix[i][k]<StatMin)	{
					statfix[i][k] = StatMin;
				}
				total += statfix[i][k];
			}
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] /= total;
			}
			weight[i] = 1.0 / Ncat;
		}
	}
	else if ((filename == "WLSR5") || (filename == "wlsr5"))	{
		Ncat = WLSR5N;
		statfix = new double*[Ncat];
		weight = new double[Ncat];
		for (int i=0; i<Ncat-1; i++)	{
			statfix[i] = new double[Nstate];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] = WLSR5StatFix[i][k];
				if (statfix[i][k]<StatMin)	{
					statfix[i][k] = StatMin;
				}
				total += statfix[i][k];
			}
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] /= total;
			}
		}
		statfix[Ncat-1] = new double[Nstate];	
		for (int k=0; k<Nstate; k++)	{
			statfix[Ncat-1][k] = -1;
		}
		for (int i=0; i<Ncat; i++)	{
			weight[i] = 1.0 / WLSR5N;
		}
	}
	else if ((filename == "c10") || (filename == "C10"))	{
		Ncat = C10N;
		statfix = new double*[Ncat];
		weight = new double[Ncat];
		for (int i=0; i<Ncat; i++)	{
			statfix[i] = new double[Nstate];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] = C10StatFix[i][k];
				if (statfix[i][k]<StatMin)	{
					statfix[i][k] = StatMin;
				}
				total += statfix[i][k];
			}
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] /= total;
			}
			weight[i] = C10StatWeight[i];
		}
	}
	else if ((filename == "c20") || (filename == "C20"))	{
		Ncat = C20N;
		statfix = new double*[Ncat];
		weight = new double[Ncat];
		for (int i=0; i<Ncat; i++)	{
			statfix[i] = new double[Nstate];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] = C20StatFix[i][k];
				if (statfix[i][k]<StatMin)	{
					statfix[i][k] = StatMin;
				}
				total += statfix[i][k];
			}
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] /= total;
			}
			weight[i] = C20StatWeight[i];
		}
	}
	else if ((filename == "c30") || (filename == "C30"))	{
		Ncat = C30N;
		statfix = new double*[Ncat];
		weight = new double[Ncat];
		for (int i=0; i<Ncat; i++)	{
			statfix[i] = new double[Nstate];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] = C30StatFix[i][k];
				if (statfix[i][k]<StatMin)	{
					statfix[i][k] = StatMin;
				}
				total += statfix[i][k];
			}
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] /= total;
			}
			weight[i] = C30StatWeight[i];
		}
	}
	else if ((filename == "c30") || (filename == "C30"))	{
		Ncat = C30N;
		statfix = new double*[Ncat];
		weight = new double[Ncat];
		for (int i=0; i<Ncat; i++)	{
			statfix[i] = new double[Nstate];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] = C30StatFix[i][k];
				if (statfix[i][k]<StatMin)	{
					statfix[i][k] = StatMin;
				}
				total += statfix[i][k];
			}
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] /= total;
			}
			weight[i] = C30StatWeight[i];
		}
	}
	else if ((filename == "c40") || (filename == "C40"))	{
		Ncat = C40N;
		statfix = new double*[Ncat];
		weight = new double[Ncat];
		for (int i=0; i<Ncat; i++)	{
			statfix[i] = new double[Nstate];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] = C40StatFix[i][k];
				if (statfix[i][k]<StatMin)	{
					statfix[i][k] = StatMin;
				}
				total += statfix[i][k];
			}
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] /= total;
			}
			weight[i] = C40StatWeight[i];
		}
	}
	else if ((filename == "c50") || (filename == "C50"))	{
		Ncat = C50N;
		statfix = new double*[Ncat];
		weight = new double[Ncat];
		for (int i=0; i<Ncat; i++)	{
			statfix[i] = new double[Nstate];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] = C50StatFix[i][k];
				if (statfix[i][k]<StatMin)	{
					statfix[i][k] = StatMin;
				}
				total += statfix[i][k];
			}
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] /= total;
			}
			weight[i] = C50StatWeight[i];
		}
	}
	else if ((filename == "c60") || (filename == "C60"))	{
		Ncat = C60N;
		statfix = new double*[Ncat];
		weight = new double[Ncat];
		for (int i=0; i<Ncat; i++)	{
			statfix[i] = new double[Nstate];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] = C60StatFix[i][k];
				if (statfix[i][k]<StatMin)	{
					statfix[i][k] = StatMin;
				}
				total += statfix[i][k];
			}
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] /= total;
			}
			weight[i] = C60StatWeight[i];
		}
	}
	else if ((filename == "cgr10") || (filename == "CGR10") || (filename == "CG10") || (filename == "cg10"))	{
		Ncat = 10;
		statfix = new double*[Ncat];
		weight = new double[Ncat];
		for (int i=0; i<Ncat; i++)	{
			statfix[i] = new double[Nstate];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] = CG10StatFix[i][k];
				if (statfix[i][k]<StatMin)	{
					statfix[i][k] = StatMin;
				}
				total += statfix[i][k];
			}
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] /= total;
			}
			weight[i] = CG10StatWeight[i];
		}
	}
	else if ((filename == "cgr20") || (filename == "CGR20") || (filename == "CG20") || (filename == "cg20"))	{
		Ncat = 20;
		statfix = new double*[Ncat];
		weight = new double[Ncat];
		for (int i=0; i<Ncat; i++)	{
			statfix[i] = new double[Nstate];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] = CG20StatFix[i][k];
				if (statfix[i][k]<StatMin)	{
					statfix[i][k] = StatMin;
				}
				total += statfix[i][k];
			}
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] /= total;
			}
			weight[i] = CG20StatWeight[i];
		}
	}
	else if ((filename == "cgr30") || (filename == "CGR30") || (filename == "CG30") || (filename == "cg30"))	{
		Ncat = 30;
		statfix = new double*[Ncat];
		weight = new double[Ncat];
		for (int i=0; i<Ncat; i++)	{
			statfix[i] = new double[Nstate];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] = CG30StatFix[i][k];
				if (statfix[i][k]<StatMin)	{
					statfix[i][k] = StatMin;
				}
				total += statfix[i][k];
			}
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] /= total;
			}
			weight[i] = CG30StatWeight[i];
		}
	}
	else if ((filename == "cgr40") || (filename == "CGR40") || (filename == "CG40") || (filename == "cg40"))	{
		Ncat = 40;
		statfix = new double*[Ncat];
		weight = new double[Ncat];
		for (int i=0; i<Ncat; i++)	{
			statfix[i] = new double[Nstate];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] = CG40StatFix[i][k];
				if (statfix[i][k]<StatMin)	{
					statfix[i][k] = StatMin;
				}
				total += statfix[i][k];
			}
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] /= total;
			}
			weight[i] = CG40StatWeight[i];
		}
	}
	else if ((filename == "cgr50") || (filename == "CGR50") || (filename == "CG50") || (filename == "cg50"))	{
		Ncat = 50;
		statfix = new double*[Ncat];
		weight = new double[Ncat];
		for (int i=0; i<Ncat; i++)	{
			statfix[i] = new double[Nstate];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] = CG50StatFix[i][k];
				if (statfix[i][k]<StatMin)	{
					statfix[i][k] = StatMin;
				}
				total += statfix[i][k];
			}
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] /= total;
			}
			weight[i] = CG50StatWeight[i];
		}
	}
	else if ((filename == "cgr60") || (filename == "CGR60") || (filename == "CG60") || (filename == "cg60"))	{
		Ncat = 60;
		statfix = new double*[Ncat];
		weight = new double[Ncat];
		for (int i=0; i<Ncat; i++)	{
			statfix[i] = new double[Nstate];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] = CG60StatFix[i][k];
				if (statfix[i][k]<StatMin)	{
					statfix[i][k] = StatMin;
				}
				total += statfix[i][k];
			}
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] /= total;
			}
			weight[i] = CG60StatWeight[i];
		}
	}
	else	{
		ifstream is(filename.c_str());
		if (!is)	{
			cerr << "error : cannot find file : " << filename << '\n';
			exit(1);
		}
		// read alphabet
		int permut[Nstate];
		for (int k=0; k<Nstate; k++)	{
			char c;
			is >> c;
			int l=0;
			while ((l<Nstate) && (c != Alphabet[l])) l++;
			if (l == Nstate)	{
				cerr << "error when reading empirical mixture in " << filename << ": does not recognise letter " << c << '\n';
				cerr << "file should be formatted as follows\n";
				cerr << "list of amino acids (e.g. A C D E ...)\n";
				cerr << "Ncat\n";
				cerr << "weight_1 (1 number) profile_1 (20 numbers)\n";
				cerr << "weight_2 (1 number) profile_2 (20 numbers)\n";
				cerr << "...\n";
				cerr << "weight_Ncat (1 number) profile_Ncat (20 numbers)\n";
				exit(1); 
			}
			permut[k] = l;
		}
		is >> Ncat;
		statfix = new double*[Ncat];
		weight = new double[Ncat];
		for (int i=0; i<Ncat; i++)	{
			statfix[i] = new double[Nstate];
			is >> weight[i];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				is >> statfix[i][permut[k]];
				if (statfix[i][permut[k]]<StatMin)	{
					statfix[i][permut[k]] = StatMin;
				}
				total += statfix[i][permut[k]];
			}
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] /= total;
			}
		}
		double total = 0;
		for (int i=0; i<Ncat; i++)	{
			total += weight[i];
		}
		for (int i=0; i<Ncat; i++)	{
			weight[i] /= total;
		}
	}
}


//---------------------------------------------------------------------------------
//		 ReadMatFix(string filename)
//---------------------------------------------------------------------------------

void MCParameters::ReadMatFix(string filename)	{

	int Nrr = Nstate * (Nstate-1) / 2;
	if ((filename == "UL2") || (filename == "ul2"))	{
		int Nrr = Nstate * (Nstate-1) / 2;
		Ncat = UL2N;
		statfix = new double*[Ncat];
		weight = new double[Ncat];
		statrr = new double*[Ncat];
		for (int i=0; i<Ncat; i++)	{
			statfix[i] = new double[Nstate];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] = UL2StatFix[i][k];
				if (statfix[i][k]<StatMin)	{
					statfix[i][k] = StatMin;
				}
				total += statfix[i][k];
			}
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] /= total;
			}
			weight[i] = 1.0 / Ncat;
			statrr[i] = new double[Nrr];
			for (int k=0; k<Nrr; k++)	{
				statrr[i][k] = UL2RRFix[i][k];
			}
		}
	}
	else if ((filename == "UL3") || (filename == "ul3"))	{
		Ncat = UL3N;
		statfix = new double*[Ncat];
		weight = new double[Ncat];
		statrr = new double*[Ncat];
		for (int i=0; i<Ncat; i++)	{
			statfix[i] = new double[Nstate];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] = UL3StatFix[i][k];
				if (statfix[i][k]<StatMin)	{
					statfix[i][k] = StatMin;
				}
				total += statfix[i][k];
			}
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] /= total;
			}
			weight[i] = 1.0 / Ncat;
			statrr[i] = new double[Nrr];
			for (int k=0; k<Nrr; k++)	{
				statrr[i][k] = UL3RRFix[i][k];
			}
		}
	}
	else	{
		ifstream is(filename.c_str());
		if (!is)	{
			cerr << "error : cannot find file : " << filename << '\n';
			exit(1);
		}
		// read alphabet
		int permut[Nstate];
		for (int k=0; k<Nstate; k++)	{
			char c;
			is >> c;
			int l=0;
			while ((l<Nstate) && (c != Alphabet[l])) l++;
			if (l == Nstate)	{
				cerr << "error when reading empirical mixture in " << filename << ": does not recognise letter " << c << '\n';
				exit(1); 
			}
			permut[k] = l;
		}
		is >> Ncat;
		statfix = new double*[Ncat];
		weight = new double[Ncat];
		statrr = new double*[Ncat];
		for (int i=0; i<Ncat; i++)	{
			statfix[i] = new double[Nstate];
			is >> weight[i];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				is >> statfix[i][permut[k]];
				if (statfix[i][permut[k]]<StatMin)	{
					statfix[i][permut[k]] = StatMin;
				}
				total += statfix[i][permut[k]];
			}
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] /= total;
			}

			statrr[i] = new double[Nrr];
			for (int k=0; k<Nstate-1; k++)	{
				for (int l=k+1; l<Nstate; l++)	{
					is >> statrr[i][Random::rrindex(permut[k],permut[l],Nstate)];
				}
			}
		}
		double total = 0;
		for (int i=0; i<Ncat; i++)	{
			total += weight[i];
		}
		for (int i=0; i<Ncat; i++)	{
			weight[i] /= total;
		}
	}
}


//---------------------------------------------------------------------------------
//		 WriteMixture(ostream& os)
//---------------------------------------------------------------------------------

void MCParameters::WriteMixture(ostream& os)	{

	for (int k=0; k<Nstate; k++)	{
		os << Alphabet[Nstate] << ' ';
	}
	os << '\n';
	os << "Ncat: " << Ncat << '\n';
	os << '\n';

	for (int i=0; i<Ncat; i++)	{
		os << "weight: " << weight[i] << '\n';
		os << '\n';
		os << "stationaries :\n";
		for (int k=0; k<Nstate; k++)	{
			os << statfix[i][k] << '\t';
		}
		os << '\n';
		if (Qmode)	{
			os << "stationaries :\n";
			int rri=0;
			for (int k=0; k<Nstate-1; k++)	{
				for (int l=k+1; l<Nstate; l++)	{
					os << statrr[i][rri] << '\t';
					rri++;
				}
				os << '\n';
			}
		}
		os << '\n';
		os << '\n';
	}
}


//---------------------------------------------------------------------------------
//		 ReadPartition(istream& is)
//---------------------------------------------------------------------------------

void MCParameters::ReadPartition(string filename)	{
	ifstream is(filename.c_str());
	if (!is)	{
		cerr << "error when reading partition: cannot open file\n";
		exit(1);
	}
	ReadPartition(is);
}

void MCParameters::ReadPartition(istream& is)	{

	cerr << "read partition\n";

	if (! DataOK)	{
		cerr << "should specify data file before partition\n";
		throw;
	}
	is >> Ngene;
	delete[] GeneSize;
	delete[] GeneFirstSite;
	delete[] Gene;
	GeneSize = new int[Ngene];	
	GeneFirstSite = new int[Ngene];
	Gene = new int[Nsite];
	int total = 0;
	for (int i=0; i<Ngene; i++)	{
		GeneFirstSite[i] = total;
		int temp;
		is >> temp;
		GeneSize[i] = temp;
		for (int j=0; j<temp; j++)	{
			Gene[j+total] = i;
		}
		total += temp;
	}
}

void MCParameters::ReadCalibration(string filename)	{
	ifstream is(filename.c_str());
	if (! is)	{
		cerr << "error when reading calibrations\n";
		exit(1);
	}
	ReadCalibration(is);
}

void MCParameters::ReadCalibration(istream& is)	{

	is >> NCalib;
	if (NCalib <= 0)	{
		cerr << "error in calibrations: incorrect number of calibrations : " << NCalib << '\n';
		exit(1);
	}	
	CalibIndex = new int[NCalib];
	CalibUpper = new double[NCalib];
	CalibLower = new double[NCalib];
	CalibTaxon1 = new string[NCalib];
	CalibTaxon2 = new string[NCalib];
	for (int i=0; i<NCalib; i++)	{
		is >> CalibTaxon1[i] >> CalibTaxon2[i];
		is >> CalibUpper[i] >> CalibLower[i];
	}
}

//---------------------------------------------------------------------------------
//		 ReadTopoMoveType(istream& is)
//---------------------------------------------------------------------------------

void	MCParameters::ReadTopoMoveType(istream& is)	{

	string temp;
	is >> temp;

	if (temp == "FixedTopo")	{
		TopoMoveTypeArray[TopoMoveTypeNumber] = FixedTopo;
		if (Parallel)	{
			int tmp;
			is >> tmp;
			for (int l=0; l<Nchain; l++)	{
				TopoNIterationArray[TopoMoveTypeNumber][l] = tmp;
			}
		}
		else	{
			for (int l=0; l<Nchain; l++)	{
				is >> TopoNIterationArray[TopoMoveTypeNumber][l];
			}
		}
	}

	else	{
		cerr << "error : only Fixed Topo Move has been implemented thus far\n";
		exit(1);
	}

	is >> temp;
	while (temp == "SuperMoveType")	{
		ReadSuperMoveType(is);
		is >> temp;
	}

	if (temp != "End")	{
		cerr << "error in TopoMoveType\n";
		cerr << temp << '\n';
		exit(1);
	}
	TopoMoveTypeNumber++;

}



//---------------------------------------------------------------------------------
//		 ReadSuperMoveType(istream& is)
//---------------------------------------------------------------------------------

void	MCParameters::ReadSuperMoveType(istream& is)	{

	int topo = TopoMoveTypeNumber;
	int super = SuperMoveTypeNumber[topo];


	if (Parallel)	{
		int tmp;
		is >> tmp;
		for (int l=0; l<Nchain; l++)	{
			SuperNIterationArray[topo][super][l] = tmp;
		}
	}
	else	{
		for (int l=0; l<Nchain; l++)	{
			is >> SuperNIterationArray[topo][super][l];
		}
	}

	string temp;
	is >> temp;

	while (temp == "MoveType")	{
		ReadMoveType(is);
		is >> temp;
	}

	if (temp != "End")	{
		cerr << "error in SuperMoveType\n";
		cerr << temp << '\n';
		exit(1);
	}
	SuperMoveTypeNumber[topo]++;

}


//---------------------------------------------------------------------------------
//		 ReadMoveType(istream& is)
//---------------------------------------------------------------------------------

void	MCParameters::ReadMoveType(istream& is)	{

	string temp;
	is >> temp;

	int topo = TopoMoveTypeNumber;
	int super = SuperMoveTypeNumber[topo];
	int move = MoveTypeNumber[topo][super];
	
	int Identified = 0;

	for (int i=0; i<ConstMoveTypeNumber; i++)	{
		if (temp == MoveTypeName[i])	{
			Identified = 1;
			MoveTypeArray[topo][super][move] = (MoveType) i;
			if (Parallel)	{
				int tmp;
				is >> tmp;
				for (int l=0; l<Nchain; l++)	{
					NIterationArray[topo][super][move][l] = tmp;
				}
				if (1)	{
					double tmp;
					is >> tmp;
					for (int l=0; l<Nchain; l++)	{
						deltaArray[topo][super][move][l] = tmp;
					}
				}
				if (1)	{
					int tmp;
					is >> tmp;
					for (int l=0; l<Nchain; l++)	{
						NArray[topo][super][move][l] = tmp;
					}
				}
			}
			else	{
				for (int l=0; l<Nchain; l++)	{
					is >> NIterationArray[topo][super][move][l];
				}
				if (1)	{
					for (int l=0; l<Nchain; l++)	{
						is >> deltaArray[topo][super][move][l];
					}
				}
				if (1)	{
					for (int l=0; l<Nchain; l++)	{
						is >> NArray[topo][super][move][l];
					}
				}
			}
		}
	}

	if (! Identified)	{
		cerr << "error : cannot recognise " << temp << '\n';
		exit(1);
	}
	
	MoveTypeNumber[topo][super]++;
}


// ---------------------------------------------------------------------------
//		 MissingData
// ---------------------------------------------------------------------------

void MCParameters::MissingData()	{

	for (int i=0; i<Ntaxa; i++)	{
		int n = 0;
		for (int j=0; j<Nsite; j++)	{
			if (Data[i][j] == unknown)	{
				n++;
			}
		}
		cout << SpeciesNames[i] << '\t' << 100 * ((double) n) / Nsite << '\n';
	}
}


// ---------------------------------------------------------------------------
//		 WriteDataToFile() (with mask)
// ---------------------------------------------------------------------------

void

MCParameters::WriteDataToFile (string filename, int * SpeciesMask)	{


	int N = 0;
	for (int i=0; i<Ntaxa; i++)	{
		N += SpeciesMask[i];
	}

	ofstream os(filename.c_str());

	os << N << '\t' << Nsite << '\n';

	unsigned int maxsize = 0;
	for (int i=0; i<Ntaxa; i++)	{
		if (SpeciesMask[i])	{
			if (maxsize < SpeciesNames[i].length())	{
				maxsize = SpeciesNames[i].length();
			}
		}
	}

	maxsize += 5;
	for (int i=0; i<Ntaxa; i++)	{
		if (SpeciesMask[i])	{
			os << SpeciesNames[i];
			for (unsigned int j=0; j<maxsize-SpeciesNames[i].length(); j++)	{
				os << ' ';
			}

			for (int j=0; j<Nsite; j++)	{
				
				if (Data[i][j] == unknown)	{
					os << '-';
				}
				else	{
					os << Alphabet[Data[i][j]];
				}

			}
			os << '\n';
		}
	}
}


void

MCParameters::WriteNexus (string filename)	{


	ofstream os(filename.c_str());

	os << "#NEXUS\n";
	os << "begin data;\n";
	os << "dimensions ntax=" << Ntaxa << " nchar=" << Nsite << ";\n";
	os << "format datatype=protein gap=?;\n";
	os << "Matrix\n";

	unsigned int maxsize = 0;
	for (int i=0; i<Ntaxa; i++)	{
		if (maxsize < SpeciesNames[i].length())	{
			maxsize = SpeciesNames[i].length();
		}
	}

	maxsize += 5;
	for (int i=0; i<Ntaxa; i++)	{
		os << SpeciesNames[i];
		for (unsigned int j=0; j<maxsize-SpeciesNames[i].length(); j++)	{
			os << ' ';
		}

		for (int j=0; j<Nsite; j++)	{
			
			if (Data[i][j] == unknown)	{
				os << '-';
			}
			else	{
				os << Alphabet[Data[i][j]];
			}

		}
		os << '\n';
	}
	os << ";\n";
	os << "end;\n";
}


// ---------------------------------------------------------------------------
//		 WriteDataToFileSiteMask() (with mask)
// ---------------------------------------------------------------------------

void

MCParameters::WriteDataToFileSiteMask(string filename, int * SiteMask)	{


	int N = 0;
	for (int i=0; i<Nsite; i++)	{
		N += SiteMask[i];
	}

	ofstream os(filename.c_str());

	os << Ntaxa << '\t' << N << '\n';
	cerr << Ntaxa << '\t' << N << '\n';

	unsigned int maxsize = 0;
	for (int i=0; i<Ntaxa; i++)	{
		if (maxsize < SpeciesNames[i].length())	{
			maxsize = SpeciesNames[i].length();
		}
	}

	maxsize += 5;
	for (int i=0; i<Ntaxa; i++)	{
		os << SpeciesNames[i];
		for (unsigned int j=0; j<maxsize-SpeciesNames[i].length(); j++)	{
			os << ' ';
		}

		for (int j=0; j<Nsite; j++)	{
			if (SiteMask[j])	{
				if (Data[i][j] == unknown)	{
					os << '-';
				}
				else	{
					os << Alphabet[Data[i][j]];
				}
			}
		}
		os << '\n';
	}
}


// ---------------------------------------------------------------------------
//		 SwapData()
// ---------------------------------------------------------------------------

void MCParameters::SwapData(MCParameters* from)	{

	if (Ntaxa != from->Ntaxa)	{
		cerr << "error in swap data: Ntaxa should be the same\n";
		exit(1);
	}

	int tmp = Nsite;
	Nsite = from->Nsite;
	from->Nsite = tmp;

	int** bkd = Data;
	Data = from->Data;
	from->Data = bkd;
	// if (ModeFastCompute)	{

	if (1)	{
		bkd = ZipData;
		ZipData = from->ZipData;
		from->ZipData = bkd;
		
		int* tmp = OrbitSize;
		OrbitSize = from->OrbitSize;
		from->OrbitSize = tmp;
		tmp = ZipSize;
		ZipSize = from->ZipSize;
		from->ZipSize = tmp;

		Boolean** bko = Orbit;
		Orbit = from->Orbit;
		from->Orbit = bko;

		bkd = Indices;
		Indices = from->Indices;
		from->Indices = bkd;
	}

	GetCurrentState()->SetData(-1);
	GetNextState()->SetData(-1);
}
	

// ---------------------------------------------------------------------------
//		 PropInvariant()
// ---------------------------------------------------------------------------

double MCParameters::PropInvariant()	{

	int obs[Nstate];
	double total  = 0;
	for (int i=0; i<Nsite; i++)	{
		for (int k=0; k<Nstate; k++)	{
			obs[k] = 0;
		}
		for (int j=0; j<Ntaxa; j++)	{
			if (Data[j][i] != unknown)	{
				obs[Data[j][i]] = 1;
			}
		}
		int orbitsize = 0;
		for (int k=0; k<Nstate; k++)	{
			orbitsize += obs[k];
		}
		if (orbitsize == 1)	{
			total++;
		}
	}
	return total / Nsite;
}
		

// ---------------------------------------------------------------------------
//		 MeanDiversity()
// ---------------------------------------------------------------------------

double MCParameters::MeanDifferentialDiversity(int* taxonmask )	{

	double total  = 0;
	int obs[Nstate];
	for (int i=0; i<Nsite; i++)	{
		int tot = 0;
		for (int j=0; j<Nstate; j++)	{
			obs[j] = 0;
		}
		for (int j=0; j<Ntaxa; j++)	{
			if (taxonmask[j])	{
				if (Data[j][i] != unknown)	{
					obs[Data[j][i]] = 1;
				}
			}
		}
		for (int j=0; j<Ntaxa; j++)	{
			if (!taxonmask[j])	{
				if (Data[j][i] != unknown)	{
					obs[Data[j][i]] = 1 - obs[Data[j][i]];
				}
			}
		}
		for (int j=0; j<Nstate; j++)	{
			tot += obs[j];
		}
		total += tot;
	}
	return total / Nsite;
}
	

double
MCParameters::MeanDiversity (int* histo, double* tmp)	{

	double total  = 0;
	int obs[Nstate];
	if (histo)	{
		for (int j=0; j<Nstate; j++)	{
			histo[j] = 0;
		}
	}
	for (int i=0; i<Nsite; i++)	{
		int tot = 0;
		for (int j=0; j<Nstate; j++)	{
			obs[j] = 0;
		}
		for (int j=0; j<Ntaxa; j++)	{
			if (Data[j][i] != unknown)	{
				obs[Data[j][i]] = 1;
			}
		}
		for (int j=0; j<Nstate; j++)	{
			tot += obs[j];
		}
		total += tot;
		if (histo)	{
			histo[tot] ++;
		}
		if (tmp)	{
			tmp[i] = tot;
		}
	}
	return total / Nsite;
}



// ---------------------------------------------------------------------------
//		 SiteDiversity()
// ---------------------------------------------------------------------------

double
MCParameters::SiteDiversity (int* site)	{

	double total  = 0;
	int obs[Nstate];
	for (int i=0; i<Nsite; i++)	{
		int tot = 0;
		for (int j=0; j<Nstate; j++)	{
			obs[j] = 0;
		}
		for (int j=0; j<Ntaxa; j++)	{
			if (Data[j][i] != unknown)	{
				obs[Data[j][i]] = 1;
			}
		}
		for (int j=0; j<Nstate; j++)	{
			tot += obs[j];
		}
		site[i] = tot;
		total += tot;
	}
	return total / Nsite;
}

// ---------------------------------------------------------------------------
//		 ReadContDataFromFile()
// ---------------------------------------------------------------------------

int MCParameters::ReadContDataFromFile (string filespec)	{
	
	if (! DataOK)	{
		cerr << "error: should read continuous data after data matrix\n";
		exit(1);
	}
	ContDataFileSpec = filespec;
	ifstream is(filespec.c_str());
	if (!is)	{
		cerr << "error when reading cont data: cannot open file\n";
		exit(1);
	}
	int ntaxa;
	is >> ntaxa;
	/*
	if (ntaxa > Ntaxa)	{
		cerr << "error in ReadContData: non matching number of taxa\n";
		exit(1);
	}
	*/
	is >> ContNchar;
	ContMissingData = new int*[Ntaxa];
	for (int j=0; j<Ntaxa; j++)	{
		ContMissingData[j] = new int[ContNchar];
		for (int i=0; i<ContNchar; i++)	{
			ContMissingData[j][i] = 1;
		}
	}
	ContData = new double*[Ntaxa];
	for (int j=0; j<Ntaxa; j++)	{
		ContData[j] = new double[ContNchar];
	}
	int totalfound = 0;
	for (int j=0; j<ntaxa; j++)	{
		string temp;
		is >> temp;
		int k = 0;
		while ((k<Ntaxa) && (SpeciesNames[k] != temp)) k++;
		if (k < Ntaxa)	{
			for (int i=0; i<ContNchar; i++)	{
				/*
				is >> temp;
				if (IsFloat(temp))	{
					ContData[k][i] = log(Double(temp));
					// ContData[k][i] = Double(temp);
					ContMissingData[k][i] = 0;
					totalfound++;
				}
				else	{
					// missing
				}
				*/
				double tmp;
				is >> tmp;
				ContData[k][i] = log(tmp);
				ContMissingData[k][i] = 0;
			}
			totalfound++;	
		}
		else	{
			double tmp;
			for (int i=0; i<ContNchar; i++)	{
				is >> tmp;
			}
		}	
	}
	// cerr << "Total found : " << totalfound << '\n';
	ComputeContVar();
	return 1;
}

void MCParameters::ComputeContVar()	{

	int ContNstate = Nstate + ContNchar;
	ContVar = new double[ContNstate];
	double ContMean[ContNstate];
	for (int k=0; k<ContNstate; k++)	{
		ContMean[k] = 0;
		ContVar[k] = 0;
	}
	double freq[Ntaxa][Nstate];
	for (int j=0; j<Ntaxa; j++)	{
		for (int k=0; k<Nstate; k++)	{
			freq[j][k] = 0;
		}
	}
	for (int j=0; j<Ntaxa; j++)	{
		int n = 0;
		for (int i=0; i<Nsite; i++)	{
			if (Data[j][i] != unknown)	{
				n++;
				freq[j][Data[j][i]]++;
			}
		}
		for (int k=0; k<Nstate; k++)	{
			freq[j][k] /= n;
			ContMean[k] += freq[j][k]; 
			ContVar[k] += freq[j][k] * freq[j][k]; 
		}
	}
	for (int k=0; k<Nstate; k++)	{
		ContMean[k] /= Ntaxa;
		ContVar[k] /= Ntaxa;
		ContVar[k] -= ContMean[k] * ContMean[k];
	}
	for (int k=0; k<ContNchar; k++)	{
		int n = 0;
		for (int j=0; j<Ntaxa; j++)	{
			if (! ContMissingData[j][k])	{
				n++;
				ContMean[Nstate+k] += ContData[j][k];
				ContVar[Nstate+k] += ContData[j][k] * ContData[j][k];
			}
		}
		ContMean[Nstate+k] /= n;
		ContVar[Nstate+k] /= n;
		ContVar[Nstate+k] -= ContMean[Nstate+k] * ContMean[Nstate +k];
	}
	for (int k=0; k<ContNstate; k++)	{
		cerr << k << '\t' << ContMean[k] << '\t' << ContVar[k] << '\n';
	}
	cerr << '\n';
}

// ---------------------------------------------------------------------------
//		 WriteDataToFile()
// ---------------------------------------------------------------------------

void

MCParameters::WriteDataToFile (string filename, int N1, int N2)	{

	if ((!N1) && (!N2))	{
		N1 = 0;
		N2 = Nsite-1;
	}
	int N= N2 - N1;
	if ( (N1<0) || (N1>Nsite) || (N2<0) || (N2>Nsite) || (N<=0) )	{
		cerr << "error in MCParameters::WriteDataToFile : bad boundaries (" << N1 << ";" << N2 << ")\n";
		exit(1);
	}

	ofstream os(filename.c_str());

	os << Ntaxa << '\t' << N+1 << '\n';

	unsigned int maxsize = 0;
	for (int i=0; i<Ntaxa; i++)	{
		if (maxsize < SpeciesNames[i].length())	{
			maxsize = SpeciesNames[i].length();
		}
	}

	maxsize += 5;
	for (int i=0; i<Ntaxa; i++)	{
		os << SpeciesNames[i];
		for (unsigned int j=0; j<maxsize-SpeciesNames[i].length(); j++)	{
			os << ' ';
		}

		for (int j=N1; j<=N2; j++)	{
			if (Data[i][j] == unknown)	{
				os << '-';
			}
			else	{
				os << Alphabet[Data[i][j]];
			}
		}
		os << '\n';
	}

}


// ---------------------------------------------------------------------------
//		 GetGeneTaxaFromData()
// ---------------------------------------------------------------------------

void MCParameters::GetGeneTaxaFromData()	{

	GeneTaxa = new int*[Ngene];
	for (int gene=0; gene<Ngene; gene++)	{
		GeneTaxa[gene] = new int[Ntaxa];
		for (int j=0; j<Ntaxa; j++)	{
			GeneTaxa[gene][j] = 1;
		}
	}
	if (NormalApprox)	{
		cerr << "WARNING : no get gene taxa under normal approx\n";
	}
	else	{
		for (int gene=0; gene<Ngene; gene++)	{
			for (int j=0; j<Ntaxa; j++)	{
				int nothing = 1;
				for (int i=0; i<GeneSize[gene]; i++)	{
					nothing &= (Data[j][GeneFirstSite[gene] + i] == unknown);
				}
				GeneTaxa[gene][j] = 1 - nothing;
			}
		}
	}
}


// ---------------------------------------------------------------------------
//		 ReadOutGroup
// ---------------------------------------------------------------------------

void MCParameters::ReadOutGroup(string filespec)	{

	ifstream is(filespec.c_str());
	if (!is)	{
		cerr << "error when reading outgroup: cannot open file\n";
		exit(1);
	}
	is >> OutNtaxa;
	OutGroup = new string[OutNtaxa];
	for (int i=0; i<OutNtaxa; i++)	{
		is >> OutGroup[i];
	}
}


// ---------------------------------------------------------------------------
//		 ReadCov
// ---------------------------------------------------------------------------

void MCParameters::ReadCov(string filespec)	{

	ConcCovFileName = filespec;
	NormalApprox = Yes;
	RootDegree = BaseReadCov(filespec, SpeciesNames, ConcTree, MaxBL, InvCov, LogDetCov, LogLMax, Ntaxa, Nbranch);
	ConcCov = 1;
	DataOK = 1;
}


// ---------------------------------------------------------------------------
//		 ReadSepCov
// ---------------------------------------------------------------------------

void MCParameters::ReadSepCov(string filespec)	{

	SepCovFileName = filespec;
	NormalApprox = Yes;
	ifstream is((Path + filespec).c_str());
	if (!is)	{
		cerr << "error when reading covariance file: cannot open file\n";
		exit(1);
	}
	is >> Ngene;
	GeneNames = new string[Ngene];
	for (int gene=0; gene<Ngene; gene++)	{
		is >> GeneNames[gene];
	}
	cerr << "number of partitions : " << Ngene << '\n';
	cerr.flush();
	GeneMaxBL = new double*[Ngene];
	GeneInvCov = new double**[Ngene];
	GeneLogDetCov = new double[Ngene];
	GeneLogLMax = new double[Ngene];
	GeneNtaxa = new int[Ngene];
	GeneNbranch = new int[Ngene];
	GeneSpeciesNames = new string*[Ngene];
	GeneTree = new Tree*[Ngene];

	RootDegree = 0;
	for (int gene=0; gene<Ngene; gene++)	{
		int temp = BaseReadCov(GeneNames[gene], GeneSpeciesNames[gene], GeneTree[gene], GeneMaxBL[gene], GeneInvCov[gene], GeneLogDetCov[gene], GeneLogLMax[gene], GeneNtaxa[gene], GeneNbranch[gene]);
		if (!RootDegree)	{
			RootDegree = temp;
		}
		else {
			if (RootDegree != temp)	{
				cerr << "error in ReadSepCov: gene trees should be either all dichotomous, or all trichotomous, at the root\n";
				exit(1);
			}
		}
	}

	SepCov = 1;
	DataOK = 1;
}

// ---------------------------------------------------------------------------
//		 BaseReadCov
// ---------------------------------------------------------------------------

int MCParameters::BaseReadCov(string filespec, string*& speciesnames, Tree*& tree, double*& maxbl, double**& invcov, double& logdetcov, double& loglmax, int& ntaxa, int& nbranch)	{

	ifstream is((Path + filespec).c_str());
	if (! is)	{
		cerr << "error: non existing file " << filespec << '\n';
		cerr << '\n';
		exit(1);
	}	
	tree = new Tree;
	tree->ReadFromStream(is,1);
	ntaxa = tree->GetSize();
	int nnode = 2*ntaxa-1;
	int noderange = 0;
	int rootdegree = tree->GetRoot()->GetDegree();
	if (rootdegree == 2)	{
		nbranch = 2 * ntaxa - 2;
		noderange = ntaxa - 1;
	}
	else	{
		nbranch = 2 * ntaxa - 3;
		noderange = ntaxa - 3;
	}
	speciesnames = new string[ntaxa];

	int range = nbranch;
	invcov = new double*[range];
	for (int i=0; i<range; i++)	{
		invcov[i] = new double[range];
	}
	maxbl = new double[range];

	GoPastNextWord(is, "follows:");
	for (int i=0; i<ntaxa; i++)	{
		int index;
		string name;
		is >> name >> index;
		speciesnames[index] = name;
	}

	GoPastNextWord(is, "follows:");
	int** nodes = new int*[nnode];
	for (int i=0; i<nnode; i++)	{
		nodes[i] = new int[nnode];
		for (int j=0; j<nnode; j++)	{
			nodes[i][j] = -1;
		}
	}
	for (int i=0; i<noderange; i++)	{
		int i1, i2, i3;
		is >> i1 >> i2 >> i3;
		if ((i1 == -1) || (i2 == -1))	{
			cerr << "error when reading cov file\n";
			cerr << i1 << '\t' << i2 << '\t' << i3 << '\n';
			exit(1);
		}
		nodes[i1][i2] = i3;
		nodes[i2][i1] = i3;
	}
	if (rootdegree == 3)	{
		int tmp;
		is >> tmp >> tmp >> tmp >> tmp;
	}

	string temp;
	double** Cov = new double*[range];
	for (int i=0; i<range; i++)	{
		Cov[i] = new double[range];
	}	
	
	tree->RegisterWithDivTime(speciesnames,nodes,maxbl,ntaxa);

	string prior;
	is >> prior;
	if (prior == "prior")	{
		Beta0 = 0;
	}
	else	{
		GoPastNextWord(is, "follows:");
		for (int i=0; i<range; i++)	{
			for (int j=0; j<range; j++)	{
				is >> Cov[i][j];
			}
		}

		double total = LinAlg::Gauss(Cov, range, invcov);

		if (rootdegree == 2)	{
			logdetcov = total;
			is >> loglmax;
			loglmax *= -1;
		}
		else	{
			logdetcov = total;
		}
	}

	// deleting
	for (int i=0; i<range; i++)	{
		delete[] Cov[i];
	}
	delete[] Cov;
	for (int i=0; i<nnode; i++)	{
		delete[] nodes[i];
	}
	delete[] nodes;
	return rootdegree;
}


// ---------------------------------------------------------------------------
//		 RegisterCov
// ---------------------------------------------------------------------------

void MCParameters::RegisterCov()	{

	if (RootDegree == 2)	{
		RegisterCovDichotomous();
	}
	else	{
		RegisterCovTrichotomous();
	}

	// checks
	if (! ConcCov)	{
		NormalConcMode = 2;
		if ((SeparateModelSwitch != 1)	|| (MSMode == Separate)) {
			// cerr << "assuming concatenate log sampling is computed from gene specific bl and cov matrices\n";
		}
	}
	else	{
		NormalConcMode = 1;
	}
	if (SepCov)	{
	}
	else	{
		if ((SeparateModelSwitch != 0)	|| (MSMode == Separate)) {
			cerr << "error : cannot use separate model under normal approx if no gene specific bl and cov matrices are specified\n";
			exit(1);
		}
	}

	double conclaplace = 0;
	double seplaplace = 0;
	if (ConcCov)	{
		conclaplace = LogLMax + 0.5 * LogDetCov;
		// cerr << "Laplace concatenate   : " << conclaplace << '\n';
	}
	if (SepCov)	{
		for (int gene=0; gene<Ngene; gene++)	{
			seplaplace += GeneLogLMax[gene] + 0.5 * GeneLogDetCov[gene];
			// cerr << gene << '\t' << GeneLogLMax[gene] << '\t' << GeneLogDetCov[gene] << '\n';
		}
		// cerr << "Laplace separate      : " << seplaplace << '\n';
	}
	if (ConcCov && SepCov)	{
		cerr << "ln BF in favor of sep : " << -seplaplace + conclaplace << '\n';
	}
}

// ---------------------------------------------------------------------------
//		 RegisterCovTrichotomous
// ---------------------------------------------------------------------------

void MCParameters::RegisterCovTrichotomous()	{

	cerr << ConcCov << '\n';
	if (ConcCov)	{
		// simple copy
		TemplateTree = new Tree(ConcTree);
	}
	else	{
		if (! SepCov)	{
			cerr << "error : cannot work under normal approx as long as no cov matrix is specified\n";
			exit(1);
		}
		// take first gene with largest set of taxa: 
		// this defines the reference tree
		// (it is assumed that this ref tree encompasses all other trees)

		Nsite = 0;
		Ntaxa = 0;
		int refgene = -1;
		string refname = "";

		for (int gene=0; gene<Ngene; gene++)	{
			if (Ntaxa < GeneNtaxa[gene])	{
				Ntaxa = GeneNtaxa[gene];
				refgene = gene;
				refname = GeneNames[gene];
			}
		}
		Nbranch = 2 * Ntaxa - 3;
		ifstream refgene_is((Path + refname).c_str());
		TemplateTree = new Tree(GeneTree[refgene]);
		SpeciesNames = new string[Ntaxa];
		TemplateTree->GetRoot()->GetSpeciesNames(SpeciesNames);
	}
	Nnode = 2 * Ntaxa - 1;

	TemplateTree->Dichotomise(1); // "1" means: arrange labels at the root
	if (TemplateTree->mRoot->label != 2 * Ntaxa - 2)	{
		cerr << "error in RegisterCovTri : root label does not match label number\n";
		exit(1);
	}

	if (ConcCov)	{
		MapBL = new int[Nnode];
		for (int i=0; i<Nnode-2; i++)	{
			MapBL[i] = i;
		}

		int i1 = TemplateTree->mRoot->down->label;
		int i2 = TemplateTree->mRoot->down->next->label;
		if (i1 == Nnode-2)	{
			MapBL[Nnode-2] = i2;
		}
		else if (i2 == Nnode-2)	{
			MapBL[Nnode-2] = i1;
		}
		else	{
			cerr << "error in MCParameters::RegisterCovTri\n";
			cerr << "labels at the root: " << i1 << ' ' << i2 << '\n';
			exit(1);
		}
		MapBL[Nnode-1] = -1;
	}
	if (SepCov)	{

		// make bipartition list
		int** bplist = new int*[Nnode-1];
		for (int i=0; i<Nnode-1; i++)	{
			bplist[i] = new int[Ntaxa];
			for (int j=0; j<Ntaxa; j++)	{
				bplist[i][j] = 0;
			}
		}
		TemplateTree->mRoot->GetRootedBPList(bplist, Nnode-1, Ntaxa);
		for (int i=0; i<Nnode-1; i++)	{
			if (bplist[i][0] == 1)	{
				for (int j=0; j<Ntaxa; j++)	{
					bplist[i][j] = 1 - bplist[i][j];
				}
			}
		}

		// for each tree
		// make bipartition list associated with branch indices
		// for each bipartition : find the corresponding bipartition in reference tree
		// --> GeneMapBL[refbip] = local index

		GeneMapBL = new int*[Ngene];
		for (int gene=0; gene<Ngene; gene++)	{
			// make a map of taxa
			int* taxonmap = new int[GeneNtaxa[gene]];
			for (int i=0; i<GeneNtaxa[gene]; i++)	{
				int j=0;
				while((j<Ntaxa) && (GeneSpeciesNames[gene][i] != SpeciesNames[j]))	{
					j++;
				}
				if (j == Ntaxa)	{
					cerr << "error in taxon mapping for gene " << gene << '\n';
					cerr << GeneSpeciesNames[gene][i] << '\n';
					exit(1);
				}
				taxonmap[i] = j;
			}
			
			GeneMapBL[gene] = new int[Nnode];	

			int** genebplist = new int*[2*GeneNtaxa[gene]-3];
			for (int i=0; i<2*GeneNtaxa[gene]-3; i++)	{
				genebplist[i] = new int[GeneNtaxa[gene]];
				for (int j=0; j<GeneNtaxa[gene]; j++)	{
					genebplist[i][j] = 0;
				}
			}
			GeneTree[gene]->mRoot->GetRootedBPList(genebplist, 2*GeneNtaxa[gene]-3, GeneNtaxa[gene]);
			for (int i=0; i<2*GeneNtaxa[gene]-3; i++)	{
				if (genebplist[i][taxonmap[0]] == 1)	{
					for (int j=0; j<GeneNtaxa[gene]; j++)	{
						genebplist[i][j] = 1 - genebplist[i][j];
					}
				}
			}

			// root
			GeneMapBL[gene][Nnode-1] = -1;

			for (int i=0; i<Nnode-1; i++)	{
				int j=0;
				int found = 0;
				while ((j<2*GeneNtaxa[gene]-3) && (!found))	{

					int compatible = 1;
					int k = 0;
					while ((k<GeneNtaxa[gene]) && compatible)	{
						compatible &= (genebplist[j][k] == bplist[i][taxonmap[k]]);
						k++;
					}
					if (compatible)	{
						found = 1;
					}
					else	{
						j++;
					}
				}
				if (!found)	{
					GeneMapBL[gene][i] = -1;
				}
				else	{
					GeneMapBL[gene][i] = j;
				}
			}

			for (int i=0; i<2*GeneNtaxa[gene]-3; i++)	{
				delete[] genebplist[i];
			}
			delete[] genebplist;
			delete[] taxonmap;
		}
	}
	DataOK = 1;
}


// ---------------------------------------------------------------------------
//		 RegisterCovDichotomous
// ---------------------------------------------------------------------------

void MCParameters::RegisterCovDichotomous()	{

	if (ConcCov)	{
		// simple copy
		TemplateTree = new Tree(ConcTree);
	}
	else	{
		if (! SepCov)	{
			cerr << "error : cannot work under normal approx as long as no cov matrix is specified\n";
			exit(1);
		}
		// take first gene with largest set of taxa: 
		// this defines the reference tree
		// (it is assumed that this ref tree encompasses all other trees)

		Nsite = 0;
		Ntaxa = 0;
		int refgene = -1;
		string refname = "";

		for (int gene=0; gene<Ngene; gene++)	{
			if (Ntaxa < GeneNtaxa[gene])	{
				Ntaxa = GeneNtaxa[gene];
				refgene = gene;
				refname = GeneNames[gene];
			}
		}
		Nbranch = 2 * Ntaxa - 2;

		ifstream refgene_is((Path + refname).c_str());
		TemplateTree = new Tree(GeneTree[refgene]);
		SpeciesNames = new string[Ntaxa];
		TemplateTree->GetRoot()->GetSpeciesNames(SpeciesNames);

	}
	Nnode = 2 * Ntaxa - 1;

	if (ConcCov)	{
		MapBL = new int[Nnode];
		for (int i=0; i<Nnode-1; i++)	{
			MapBL[i] = i;
		}
		MapBL[Nnode-1] = -1;
	}

	if (SepCov)	{

		// make bipartition list
		int** bplist = new int*[Nnode-1];
		for (int i=0; i<Nnode-1; i++)	{
			bplist[i] = new int[Ntaxa];
			for (int j=0; j<Ntaxa; j++)	{
				bplist[i][j] = 0;
			}
		}
		TemplateTree->mRoot->GetRootedBPList(bplist, Nnode-1, Ntaxa);

		// for each tree
		// make bipartition list associated with branch indices
		// for each bipartition : find the corresponding bipartition in reference tree
		// --> GeneMapBL[refbip] = local index

		GeneMapBL = new int*[Ngene];
		for (int gene=0; gene<Ngene; gene++)	{
			// make a map of taxa
			int* taxonmap = new int[GeneNtaxa[gene]];
			for (int i=0; i<GeneNtaxa[gene]; i++)	{
				int j=0;
				while((j<Ntaxa) && (GeneSpeciesNames[gene][i] != SpeciesNames[j]))	{
					j++;
				}
				if (j == Ntaxa)	{
					cerr << "error in taxon mapping for gene " << gene << '\n';
					cerr << GeneSpeciesNames[gene][i] << '\n';
					exit(1);
				}
				taxonmap[i] = j;
			}
			
			GeneMapBL[gene] = new int[Nnode];	

			int** genebplist = new int*[2*GeneNtaxa[gene]-2];
			for (int i=0; i<2*GeneNtaxa[gene]-2; i++)	{
				genebplist[i] = new int[GeneNtaxa[gene]];
				for (int j=0; j<GeneNtaxa[gene]; j++)	{
					genebplist[i][j] = 0;
				}
			}
			GeneTree[gene]->mRoot->GetRootedBPList(genebplist, 2*GeneNtaxa[gene]-2, GeneNtaxa[gene]);

			// root
			GeneMapBL[gene][Nnode-1] = -1;

			for (int i=0; i<Nnode-1; i++)	{
				int j=0;
				int found = 0;
				while ((j<2*GeneNtaxa[gene]-2) && (!found))	{
					int compatible = 1;
					int k = 0;
					while ((k<GeneNtaxa[gene]) && compatible)	{
						compatible &= (genebplist[j][k] == bplist[i][taxonmap[k]]);
						k++;
					}
					if (compatible)	{
						found = 1;
					}
					else	{
						j++;
					}
				}
				if (!found)	{
					GeneMapBL[gene][i] = -1;
				}
				else	{
					GeneMapBL[gene][i] = j;
				}
			}

			for (int i=0; i<2*GeneNtaxa[gene]-2; i++)	{
				delete[] genebplist[i];
			}
			delete[] genebplist;
			delete[] taxonmap;
		}
	}
	DataOK = 1;
}


// ---------------------------------------------------------------------------
//		 ReadDataFromFile()
// ---------------------------------------------------------------------------

int MCParameters::ReadDataFromFile (string filespec, int forceinterleaved)	{
	
	DataFileSpec = filespec;
	NormalApprox = No;
	string tmp;
	ifstream is((Path + filespec).c_str());
	if (!is)	{
		cerr << "error : cannot find data file " << filespec << '\n';
		cerr << "\n";
		exit(1);
	}
	is >> tmp;
	try{
		if (tmp == "#SPECIALALPHABET")	{
			ReadSpecial();
			return 1;
		}
		else if (tmp == "#NEXUS")	{
			ReadNexus();
			return 1;
		}
		else	{
			if (! forceinterleaved)	{
				// cerr << "test sequential phylip\n";
				int returnvalue = TestPhylipSequential();
				if (returnvalue)	{
					// cerr << "assuming phylip sequential format\n";
					ReadPhylipSequential();
					return 1;
				}
			}
			// cerr << "test interleaved phylip\n";
			int returnvalue = TestPhylip(1);
			if (returnvalue)	{
				// cerr << "assuming phylip interleaved format\n";
				ReadPhylip(1);
				return 1;
			}
			TestPhylip(0);
			// cerr << "assuming phylip interleaved format\n";
			ReadPhylip(0);
			return 1;
		}
	}
	catch(...)	{
		exit(1);
		return 0;
	}
	return 1;
}

int MCParameters::ReadSpecial()	{

	int returnvalue = 0;
	try	{

		ifstream theStream((Path + DataFileSpec).c_str());

		string tmp;
		theStream >> tmp;
		theStream >> Ntaxa;
		theStream >> Nsite;
		theStream >> tmp;
		cerr << tmp << '\n';
		Nstate = tmp.length();
		
		Alphabet = new char[Nstate];
		NAlphabetSet = Nstate+5;
		Alphabet = new char[Nstate];
		AlphabetSet = new char[NAlphabetSet];
		cerr << "alphabet size : " << Nstate << '\n';
		cerr << "alphabet : ";
		for (int i=0; i<Nstate; i++)	{
			Alphabet[i] = tmp[i];
			AlphabetSet[i] = tmp[i];
			cerr << Alphabet[i] << ' ';
		}
		cerr << '\n';
		returnvalue = 4;

		AlphabetSet[Nstate] = '?';
		AlphabetSet[Nstate+1] = '-';
		AlphabetSet[Nstate+2] = '*';
		AlphabetSet[Nstate+3] = 'X';
		AlphabetSet[Nstate+4] = 'x';

		Nnode = 2*Ntaxa -1;

		if (Data)	{
			for (int i=0; i<Ntaxa; i++)	{
				delete Data[i];
			}
			delete[] Data;
		}
		Data = new int*[Ntaxa];
		for (int i=0; i<Ntaxa; i++)	{
			Data[i] = new int[Nsite];
		}

		delete[] SpeciesNames;
		SpeciesNames = new string[Ntaxa];

		int ntaxa = 0;
		string temp;
		while ((!theStream.eof()) && (ntaxa<Ntaxa))	{
			theStream >> temp;
			SpeciesNames[ntaxa] = temp;
			int nsite = 0;

			char c;
			do	{
				c = theStream.get();
				if ((!theStream.eof()) && (c != ' ') && (c != '\n') && (c != '\t') && (c!=13))	{
					if (c == '(')	{
						Data[ntaxa][nsite] = unknown;
						while (c != ')')	{
							theStream >> c;
						}
					}
					else if (c == '{')	{
						Data[ntaxa][nsite] = unknown;
						while (c != '}')	{
							theStream >> c;
						}
					}
					else	{
						int p =0;
						while ((p < NAlphabetSet) && (c != AlphabetSet[p])) p++;
						if (p == NAlphabetSet)	{
							cout << "error: does not recognise character. taxon " << ntaxa << '\t' << SpeciesNames[ntaxa] << "  site  " << nsite << '\t' << c << '\n';
							exit(1);
						}
						if (p >= Nstate)	{
							Data[ntaxa][nsite] = unknown;
						}
						else	{
							for (int l=0; l<Nstate; l++)		{
								if (c == Alphabet[l])	{
									Data[ntaxa][nsite] = l;
								}
							}
						}
					}
					nsite++;
				}
			}
			while ((!theStream.eof()) && (nsite < Nsite));
			ntaxa++;
		}
		DataOK = 1;
	}
	catch(...)	{
		cerr << "error while reading data file\n";
		return 0;
	}
	return returnvalue;
}

int MCParameters::ReadNexus()	{

	int returnvalue = 0;
	try	{

		if (verbose)	{
			cerr << "read data from file : " << DataFileSpec << "\n";
			cerr.flush();
		}

		ifstream theStream((Path + DataFileSpec).c_str());

		GoPastNextWord(theStream, "dimensions");
		GoPastNext(theStream, '=');
		theStream >> Ntaxa;
		GoPastNext(theStream, '=');
		theStream >> Nsite;
		GoPastNextWord(theStream, "format");
		string tmp;
		theStream >> tmp;
		if (tmp == "symbols")	{
			theStream >> tmp >> tmp;
			cerr << tmp << '\n';
			if (tmp[0] != '\"')	{
				cerr << "error: does not recognize symbol list : " << tmp << '\n';
				exit(1);
			}
			unsigned int l = tmp.length();
			if (tmp[l-2] == '\"')	{
				Nstate = l-3;
			}
			else if (tmp[l-1] == '\"')	{
				Nstate = l-2;
			}
			else	{
				cerr << "error: does not recognize symbol list : " << tmp << '\n';
				exit(1);
			}
			
			Alphabet = new char[Nstate];
			NAlphabetSet = Nstate;
			Alphabet = new char[Nstate];
			AlphabetSet = new char[Nstate];
			cerr << "alphabet size : " << Nstate << '\n';
			cerr << "alphabet : ";
			for (int i=0; i<Nstate; i++)	{
				Alphabet[i] = tmp[i+1];
				AlphabetSet[i] = tmp[i+1];
				cerr << Alphabet[i] << ' ';
			}
			cerr << '\n';
			returnvalue = 4;
		}
		else	{
		// GoPastNextWord(theStream, "datatype");
		GoPastNext(theStream, '=');
		string type;
		theStream >> type;
		if (EquivalentStrings(type,"protein"))	{
			Nstate = 20;
			Alphabet = new char[Nstate];
			for (int i=0; i<Nstate; i++)	{
				Alphabet[i] = AminoAcids[i];
			}
			NAlphabetSet = AAN;
			AlphabetSet = new char[NAlphabetSet];
			for (int i=0; i<NAlphabetSet; i++)	{
				AlphabetSet[i] = AAset[i];
			}
			// cerr << "protein sequences\n";
			returnvalue = 3;
		}
		else if (EquivalentStrings(type,"dna"))	{
			Nstate = 4;
			Alphabet = new char[Nstate];
			for (int i=0; i<Nstate; i++)	{
				Alphabet[i] = DNAletters[i];
			}
			NAlphabetSet = DNAN;
			AlphabetSet = new char[NAlphabetSet];
			for (int i=0; i<NAlphabetSet; i++)	{
				AlphabetSet[i] = DNAset[i];
			}
			cerr << "dna sequences\n";
			returnvalue = 1;
		}
		else if (EquivalentStrings(type,"rna"))	{
			Nstate = 4;
			Alphabet = new char[Nstate];
			for (int i=0; i<Nstate; i++)	{
				Alphabet[i] = RNAletters[i];
			}
			NAlphabetSet = RNAN;
			AlphabetSet = new char[NAlphabetSet];
			for (int i=0; i<NAlphabetSet; i++)	{
				AlphabetSet[i] = RNAset[i];
			}
			cerr << "rna sequences\n";
			returnvalue = 2;
		}
		else	{
			cerr << "error in nexus: cannot recognise data type\n";
			cerr << type << "\n";
			exit(1);
		}
		}
		
		if (verbose)	{
			cerr << "Ntaxa	:" << Ntaxa << '\n';
			cerr << "Nsite	:" << Nsite << '\n';
			cerr.flush();
		}

		Nnode = 2*Ntaxa -1;

		if (Data)	{
			for (int i=0; i<Ntaxa; i++)	{
				delete Data[i];
			}
			delete[] Data;
		}
		Data = new int*[Ntaxa];
		for (int i=0; i<Ntaxa; i++)	{
			Data[i] = new int[Nsite];
		}

		delete[] SpeciesNames;
		SpeciesNames = new string[Ntaxa];

		GoPastNextWord(theStream, "Matrix");

		int l = 0;
		while (l<Nsite)	{
			cerr.flush();
			int m = 0;
			for (int i=0; i<Ntaxa; i++)	{
				
				string temp;
				theStream >> temp;
				while (temp == "[")	{
					unsigned char c;
					c = 'i';
					while (c != ']') c = theStream.get();
					theStream >> temp;
				}

				if (!l)	{
					SpeciesNames[i] = temp;
				}
				else	{
					if (temp != SpeciesNames[i])	{
						cerr << "error when reading tree base: " << temp << '\t' << SpeciesNames[i] << '\n';
						exit(1);
					}
				}

				unsigned char c;
				int k = l;
				do	{
					c = theStream.get();
					if (c == '[')	{
						while (c != ']') c = theStream.get();
						c = theStream.get();
					}
					if ((c != ' ') && (c != '\t') && (c != '\n') && (c != 13))	{
						if (c == '(')	{
							Data[i][k] = unknown;
							while (c != ')')	{
								theStream >> c;
							}
						}
						else if (c == '{')	{
							Data[i][k] = unknown;
							while (c != '}')	{
								theStream >> c;
							}
						}
						else	{
							if (returnvalue == 4)	{
								if ((c == '-') || (c == '?') || (c == '*'))	{
									Data[i][k] = unknown;
								}
								else	{
									int p =0;
									while ((p < NAlphabetSet) && (c != AlphabetSet[p])) p++;
									if (p == NAlphabetSet)	{
										cout << "error: does not recognise character. taxon " << i << '\t' << SpeciesNames[i] << "  site  " << k << '\t' << c << '\n';
										exit(1);
									}
									for (int l=0; l<Nstate; l++)		{
										if ((c == Alphabet[l]) || (c == Alphabet[l]+32))	{
											Data[i][k] = l;
										}
									}
								}
							}
							else	{
								int p =0;
								while ((p < NAlphabetSet) && (c != AlphabetSet[p])) p++;
								if (p == NAlphabetSet)	{
									cout << "error: does not recognise character. taxon " << i << '\t' << SpeciesNames[i] << "  site  " << k << '\t' << c << '\n';
									exit(1);
								}
								if (p >= 2*Nstate)	{
									Data[i][k] = unknown;
								}
								else	{
									for (int l=0; l<Nstate; l++)		{
										if ((c == Alphabet[l]) || (c == Alphabet[l]+32))	{
											Data[i][k] = l;
										}
									}
								}
							}
						}
						k++;
					}
				}
				while ((!theStream.eof()) && (c != '\n') && (c != 13));
				if (theStream.eof())	{
					if (i < Ntaxa-1)	{
						cerr << "error : found " << i << " taxa instead of " << Ntaxa << " in datafile\n";
						exit(1);
					}
				}
				if (!m)	{
					m = k;
				}
				else	{
					if (m != k)	{
						cerr << "error when reading nexus : " << m << '\t' << k << '\n';
						cerr << "taxa : " << i << '\t' << SpeciesNames[i] << '\n';
						if (m > k)	{
							while (k != m)	{
								Data[i][k] = unknown;
								k++;
							}
						}
					}
				}
			}
			l= m;
		}

		/*
		EliminateUnknownColumns();
		if (DeleteConstant)	{
			EliminateConstantPositions();
		}
		*/
		DataOK = 1;
	}
	catch(...)	{
		cerr << "error while reading data file\n";
		return 0;
	}
	return returnvalue;
}


// ---------------------------------------------------------------------------
//		 EliminateUnnownColumns
// ---------------------------------------------------------------------------


void MCParameters::EliminateUnknownColumns()	{
		
		// cerr << "eliminate unknown columns\n";
	
		// delete all the sites that are '-' for every species

		int i=0;
		int j=0;
		int Eliminated = 0;
		while (i<Nsite)	{
			int k=0;
			Boolean test = true;
			while (test && 	k<Ntaxa)	{
				test &= (Data[k][i] == unknown);
				k++;
			}
			if (! test)	{
				for (int k=0; k<Ntaxa; k++)	{
					Data[k][j] = Data[k][i];
				}
				j++;
			}
			else	{
				GeneSize[Gene[i]]--;
				Eliminated ++;
			}
			i++;
		}
		Nsite -= Eliminated;
		if (Eliminated)	{
			cerr << Eliminated << " columns completely undetermined (full of gaps, or unknown): eliminated\n";
		}
}

// ---------------------------------------------------------------------------
//		 ReadPhylip()
// ---------------------------------------------------------------------------

int MCParameters::TestPhylipSequential ()	{

	try	{

		// cerr << "beware: phylip data sets only for amino acids for the moment\n";
		// cerr.flush();
	
		ifstream theStream((Path + DataFileSpec).c_str());

		string temp;
		theStream >> temp;
		if (!IsInt(temp))	{
			cerr << "error when reading data\n";
			cerr << "data should be formatted as follows:\n";
			cerr << "#taxa #sites\n";
			cerr << "name1 seq1.....\n";
			cerr << "name2 seq2.....\n";
			cerr << "...\n";
			cerr << '\n';
			exit(1);
		}
		Ntaxa = Int(temp);
		theStream >> temp;
		if (!IsInt(temp))	{
			cerr << "error when reading data\n";
			cerr << "data should be formatted as follows:\n";
			cerr << "#taxa #sites\n";
			cerr << "name1 seq1.....\n";
			cerr << "name2 seq2.....\n";
			cerr << "...\n";
			cerr << '\n';
			exit(1);
		}
		Nsite = Int(temp);
		if (verbose)	{
			cerr << "Ntaxa	:" << Ntaxa << '\n';
			cerr << "Nsite	:" << Nsite << '\n';
			cerr.flush();
		}

		delete[] SpeciesNames;
		SpeciesNames = new string[Ntaxa];

		int AAcomp = 1;
		int DNAcomp = 1;
		int RNAcomp = 1;

		int ntaxa = 0;
		while ((!theStream.eof()) && (ntaxa<Ntaxa))	{
			theStream >> temp;
			SpeciesNames[ntaxa] = temp;
			// cerr << temp << '\n';
			int nsite = 0;

			char c = ' ';
			do	{
				c = theStream.get();
				if ((!theStream.eof()) && (c != ' ') && (c != '\n') && (c!='\t') && (c != 13))	{
					if (c == '(')	{
						while (c != ')')	{
							theStream >> c;
						}
					}
					else if (c == '{')	{
						while (c != '}')	{
							theStream >> c;
						}
					}
					else	{
						int p =0;
						if (DNAcomp)	{
							while ((p < DNAN) && (c != DNAset[p])) p++;
							if (p == DNAN) {
								// cerr << c << '\t';
								DNAcomp = 0;
							}
						}
						p =0;
						if (RNAcomp)	{
							while ((p < RNAN) && (c != RNAset[p])) p++;
							if (p == RNAN) RNAcomp = 0;
						}
						p =0;
						if (AAcomp)	{
							while ((p < AAN) && (c != AAset[p])) p++;
							if (p == AAN)	{
								AAcomp = 0;
							}
						}
					}
					nsite++;
				}
			}
			while ((!theStream.eof()) && (nsite < Nsite));
			if (theStream.eof())	{
				if (nsite < Nsite)	{
					// cerr << "taxon : " << ntaxa << " : " << SpeciesNames[ntaxa]  << "contain only " << nsite << " sites\n"; 
					return 0;
				}
			}
			ntaxa++;
		}
		if (theStream.eof())	{
			if (ntaxa < Ntaxa)	{
				// cerr << "found only " << ntaxa << " taxa\n";
				return 0;
			}
		}
		if (DNAcomp)	{
			Nstate = 4;
			Alphabet = new char[Nstate];
			for (int i=0; i<Nstate; i++)	{
				Alphabet[i] = DNAletters[i];
			}
			NAlphabetSet = DNAN;
			AlphabetSet = new char[NAlphabetSet];
			for (int i=0; i<NAlphabetSet; i++)	{
				AlphabetSet[i] = DNAset[i];
			}
			cerr << "dna sequences\n";
		}
		else if (RNAcomp)	{
			Nstate = 4;
			Alphabet = new char[Nstate];
			for (int i=0; i<Nstate; i++)	{
				Alphabet[i] = RNAletters[i];
			}
			NAlphabetSet = RNAN;
			AlphabetSet = new char[NAlphabetSet];
			for (int i=0; i<NAlphabetSet; i++)	{
				AlphabetSet[i] = RNAset[i];
			}
			cerr << "rna sequences\n";
		}
		else if (AAcomp)	{
			Nstate = 20;
			Alphabet = new char[Nstate];
			for (int i=0; i<Nstate; i++)	{
				Alphabet[i] = AminoAcids[i];
			}
			NAlphabetSet = AAN;
			AlphabetSet = new char[NAlphabetSet];
			for (int i=0; i<NAlphabetSet; i++)	{
				AlphabetSet[i] = AAset[i];
			}
			// cerr << "protein sequences\n";
		}
		else	{
			// cerr << "format not recognised\n";
			return 0;
		}
	}
	catch(...)	{
		return 0;
	}
	return 1;
}


void MCParameters::ReadPhylipSequential ()	{

	try	{

		// cerr << "beware: phylip data sets only for amino acids for the moment\n";
		// cerr.flush();
	
		ifstream theStream((Path + DataFileSpec).c_str());

		string temp;
		theStream >> temp;
		if (!IsInt(temp))	{
			cerr << "error when reading data\n";
			cerr << "data should be formatted as follows:\n";
			cerr << "#taxa #sites\n";
			cerr << "name1 seq1.....\n";
			cerr << "name2 seq2.....\n";
			cerr << "...\n";
			cerr << '\n';
			exit(1);
		}
		Ntaxa = Int(temp);
		theStream >> temp;
		if (!IsInt(temp))	{
			cerr << "error when reading data\n";
			cerr << "data should be formatted as follows:\n";
			cerr << "#taxa #sites\n";
			cerr << "name1 seq1.....\n";
			cerr << "name2 seq2.....\n";
			cerr << "...\n";
			cerr << '\n';
			exit(1);
		}
		Nsite = Int(temp);
		if (verbose)	{
			cerr << "Ntaxa	:" << Ntaxa << '\n';
			cerr << "Nsite	:" << Nsite << '\n';
			cerr.flush();
		}

		Nnode = 2*Ntaxa -1;
				
		if (Data)	{
			for (int i=0; i<Ntaxa; i++)	{
				delete Data[i];
			}
			delete[] Data;
		}
		Data = new int*[Ntaxa];
		for (int i=0; i<Ntaxa; i++)	{
			Data[i] = new int[Nsite];
		}

		delete[] SpeciesNames;
		SpeciesNames = new string[Ntaxa];

		int ntaxa = 0;
		while ((!theStream.eof()) && (ntaxa<Ntaxa))	{
			theStream >> temp;
			SpeciesNames[ntaxa] = temp;
			int nsite = 0;

			char c;
			do	{
				c = theStream.get();
				if ((!theStream.eof()) && (c != ' ') && (c != '\n') && (c != '\t') && (c!=13))	{
					if (c == '(')	{
						Data[ntaxa][nsite] = unknown;
						while (c != ')')	{
							theStream >> c;
						}
					}
					else if (c == '{')	{
						Data[ntaxa][nsite] = unknown;
						while (c != '}')	{
							theStream >> c;
						}
					}
					else	{
						int p =0;
						while ((p < NAlphabetSet) && (c != AlphabetSet[p])) p++;
						if (p == NAlphabetSet)	{
							cout << "error: does not recognise character. taxon " << ntaxa << '\t' << SpeciesNames[ntaxa] << "  site  " << nsite << '\t' << c << '\n';
							exit(1);
						}
						if (p >= 2*Nstate)	{
							Data[ntaxa][nsite] = unknown;
						}
						else	{
							for (int l=0; l<Nstate; l++)		{
								if ((c == Alphabet[l]) || (c == Alphabet[l]+32))	{
									Data[ntaxa][nsite] = l;
								}
							}
						}
					}
					nsite++;
				}
			}
			while ((!theStream.eof()) && (nsite < Nsite));
			ntaxa++;
		}

		// delete all the sites that are '-' for every species
		//EliminateUnknownColumns();

		DataOK = 1;
	}
	catch(...)	{
		cerr << "error while reading data file\n";
		exit(1);
	}
}

int MCParameters::TestPhylip (int repeattaxa)	{

	int returnvalue = 0;
	try	{

		// cerr << "beware: phylip data sets only for amino acids for the moment\n";
		// cerr.flush();
	
		ifstream theStream((Path + DataFileSpec).c_str());

		string temp;
		theStream >> temp;
		if (!IsInt(temp))	{
			cerr << "error when reading data\n";
			cerr << "data should be formatted as follows:\n";
			cerr << "#taxa #sites\n";
			cerr << "name1 seq1.....\n";
			cerr << "name2 seq2.....\n";
			cerr << "...\n";
			cerr << '\n';
			exit(1);
		}
		Ntaxa = Int(temp);
		theStream >> temp;
		if (!IsInt(temp))	{
			cerr << "error when reading data\n";
			cerr << "data should be formatted as follows:\n";
			cerr << "#taxa #sites\n";
			cerr << "name1 seq1.....\n";
			cerr << "name2 seq2.....\n";
			cerr << "...\n";
			cerr << '\n';
			exit(1);
		}
		Nsite = Int(temp);
		if (verbose)	{
			cerr << "Ntaxa	:" << Ntaxa << '\n';
			cerr << "Nsite	:" << Nsite << '\n';
			cerr.flush();
		}

		delete[] SpeciesNames;
		SpeciesNames = new string[Ntaxa];

		int AAcomp = 1;
		int DNAcomp = 1;
		int RNAcomp = 1;

		int l = 0;
		int block = 0;
		while (l<Nsite)	{
			block++;
			int m = 0;
			for (int i=0; i<Ntaxa; i++)	{
				
				if ((!l) || repeattaxa)	{
					string temp;
					theStream >> temp;
					if (!l)	{
						SpeciesNames[i] = temp;
					}
					else	{
						if (temp != SpeciesNames[i])	{
							return 0;
						}
					}
				}

				unsigned char c;
				int k = l;
				do	{
					c = theStream.get();
					if ((!theStream.eof()) && (c != ' ') && (c != '\n') && (c!='\t') && (c != 13))	{
						if (c == '(')	{
							while (c != ')')	{
								theStream >> c;
							}
						}
						else if (c == '{')	{
							while (c != '}')	{
								theStream >> c;
							}
						}
						else	{
							int p =0;
							if (DNAcomp)	{
								while ((p < DNAN) && (c != DNAset[p])) p++;
								if (p == DNAN) DNAcomp = 0;
							}
							p =0;
							if (RNAcomp)	{
								while ((p < RNAN) && (c != RNAset[p])) p++;
								if (p == RNAN) RNAcomp = 0;
							}
							p =0;
							if (AAcomp)	{
								while ((p < AAN) && (c != AAset[p])) p++;
								if (p == AAN)	{
									AAcomp = 0;
								}
							}
						}
						k++;
					}
				}
				while ((!theStream.eof()) && (c != '\n') && (c != 13) && (c != 10));
				if (theStream.eof())	{
					if (i < Ntaxa-1)	{
						cerr << "error : found " << i << " taxa instead of " << Ntaxa << " in datafile\n";
						exit(1);
					}
				}
				c = theStream.peek();
				while ((!theStream.eof()) && ((c == '\n') || (c == 13)))	{
					c = theStream.get();
					c = theStream.peek();
				}
				if (!m)	{
					m = k;
				}
				else	{
					if (m != k)	{
						cerr << "in test phylip\n";
						cerr << "error when reading data non matching number of sequences in block number " << block << " for taxon " << i+1 << " " << SpeciesNames[i] << '\n';
						cerr << "taxa : " << i << '\t' << SpeciesNames[i] << '\n';
						cerr << "read " << k << " instead of " << m << "characters\n";
						exit(1);
					}
				}
			}
			l= m;
		}
		if (l<Nsite)	{
			cerr << "error : reached end of stream \n";
			cerr << "data should be formatted as follows:\n";
			cerr << "#taxa #sites\n";
			cerr << "name1 seq1.....\n";
			cerr << "name2 seq2.....\n";
			cerr << "...\n";
			cerr << '\n';
			exit(1);
		}
		if (DNAcomp)	{
			Nstate = 4;
			Alphabet = new char[Nstate];
			for (int i=0; i<Nstate; i++)	{
				Alphabet[i] = DNAletters[i];
			}
			NAlphabetSet = DNAN;
			AlphabetSet = new char[NAlphabetSet];
			for (int i=0; i<NAlphabetSet; i++)	{
				AlphabetSet[i] = DNAset[i];
			}
			cerr << "dna sequences\n";
			returnvalue = 1;
		}
		else if (RNAcomp)	{
			Nstate = 4;
			Alphabet = new char[Nstate];
			for (int i=0; i<Nstate; i++)	{
				Alphabet[i] = RNAletters[i];
			}
			NAlphabetSet = RNAN;
			AlphabetSet = new char[NAlphabetSet];
			for (int i=0; i<NAlphabetSet; i++)	{
				AlphabetSet[i] = RNAset[i];
			}
			cerr << "rna sequences\n";
			returnvalue = 1;
		}
		else if (AAcomp)	{
			Nstate = 20;
			Alphabet = new char[Nstate];
			for (int i=0; i<Nstate; i++)	{
				Alphabet[i] = AminoAcids[i];
			}
			NAlphabetSet = AAN;
			AlphabetSet = new char[NAlphabetSet];
			for (int i=0; i<NAlphabetSet; i++)	{
				AlphabetSet[i] = AAset[i];
			}
			// cerr << "protein sequences\n";
			returnvalue = 1;
		}
		else	{
			cerr << "error : cannot recognise data type\n";
			exit(1);
		}
	}
	catch(...)	{
		cerr << "error while reading data file\n";
		return 0;
	}
	return returnvalue;
}

void
MCParameters::ReadPhylip (int repeattaxa)	{

	try	{

		// cerr << "beware: phylip data sets only for amino acids for the moment\n";
		// cerr.flush();
	
		ifstream theStream((Path + DataFileSpec).c_str());

		string temp;
		theStream >> temp;
		if (!IsInt(temp))	{
			cerr << "error when reading data\n";
			cerr << "data should be formatted as follows:\n";
			cerr << "#taxa #sites\n";
			cerr << "name1 seq1.....\n";
			cerr << "name2 seq2.....\n";
			cerr << "...\n";
			cerr << '\n';
			exit(1);
		}
		Ntaxa = Int(temp);
		theStream >> temp;
		if (!IsInt(temp))	{
			cerr << "error when reading data\n";
			cerr << "data should be formatted as follows:\n";
			cerr << "#taxa #sites\n";
			cerr << "name1 seq1.....\n";
			cerr << "name2 seq2.....\n";
			cerr << "...\n";
			cerr << '\n';
			exit(1);
		}
		Nsite = Int(temp);
		if (verbose)	{
			cerr << "Ntaxa	:" << Ntaxa << '\n';
			cerr << "Nsite	:" << Nsite << '\n';
			cerr.flush();
		}

		Nnode = 2*Ntaxa -1;
				
		if (Data)	{
			for (int i=0; i<Ntaxa; i++)	{
				delete Data[i];
			}
			delete[] Data;
		}
		Data = new int*[Ntaxa];
		for (int i=0; i<Ntaxa; i++)	{
			Data[i] = new int[Nsite];
		}

		delete[] SpeciesNames;
		SpeciesNames = new string[Ntaxa];

		int l = 0;
		int block = 0;
		while (l<Nsite)	{
			block++;
			cerr.flush();
			int m = 0;
			for (int i=0; i<Ntaxa; i++)	{
				
				if ((!l) || repeattaxa)	{
					string temp;
					theStream >> temp;
					if (!l)	{
						SpeciesNames[i] = temp;
					}
					else	{
						if (temp != SpeciesNames[i])	{
							cerr << "error when reading data: read " << temp << " instead of " << SpeciesNames[i] << '\n';
							exit(1);
						}
					}
				}

				unsigned char c;
				int k = l;
				do	{

					c = theStream.get();
					if ((!theStream.eof()) && (c != ' ') && (c != '\n') && (c != '\t') && (c!=13))	{
						if (c == '(')	{
							Data[i][k] = unknown;
							while (c != ')')	{
								theStream >> c;
							}
						}
						else if (c == '{')	{
							Data[i][k] = unknown;
							while (c != '}')	{
								theStream >> c;
							}
						}
						else	{
							int p =0;
							while ((p < NAlphabetSet) && (c != AlphabetSet[p])) p++;
							if (p == NAlphabetSet)	{
								cout << "error: does not recognise character. taxon " << i << '\t' << SpeciesNames[i] << "  site  " << k << '\t' << c << '\n';
								exit(1);
							}
							if (p >= 2*Nstate)	{
								Data[i][k] = unknown;
							}
							else	{
								for (int l=0; l<Nstate; l++)		{
									if ((c == Alphabet[l]) || (c == Alphabet[l]+32))	{
										Data[i][k] = l;
									}
								}
							}
						}
						k++;
					}
				}
				while ((!theStream.eof()) && (c != '\n') && (c != 13));
				if (theStream.eof())	{
					if (i < Ntaxa-1)	{
						cerr << "error : found " << i << " taxa instead of " << Ntaxa << " in datafile\n";
						exit(1);
					}
				}
				c = theStream.peek();
				while ((!theStream.eof()) && ((c == '\n') || (c == 13)))	{
					c = theStream.get();
					c = theStream.peek();
				}
				
				if (!m)	{
					m = k;
				}
				else	{
					if (m != k)	{
						cerr << "error when reading data non matching number of sequences in block number " << block << " for taxon " << i << " " << SpeciesNames[i] << '\n';
						cerr << "taxa : " << i << '\t' << SpeciesNames[i] << '\n';
						cerr << "read " << k << " instead of " << m << "characters\n";
						exit(1);
					}
				}
			}
			l= m;
		}
		if (l<Nsite)	{
			cerr << "error : reached end of stream \n";
			cerr << "data should be formatted as follows:\n";
			cerr << "#taxa #sites\n";
			cerr << "name1 seq1.....\n";
			cerr << "name2 seq2.....\n";
			cerr << "...\n";
			cerr << '\n';
			exit(1);
		}

		// delete all the sites that are '-' for every species
		//EliminateUnknownColumns();

		DataOK = 1;
	}
	catch(...)	{
		cerr << "error while reading data file\n";
	}
}

// ---------------------------------------------------------------------------
//		 EstimateEmpiricalFrequencies()
// ---------------------------------------------------------------------------

void MCParameters::EstimateEmpiricalFrequencies()	{

	// cerr << "estimate emp freq\n";
	for (int k=0; k<Nstate; k++)	{
		EmpiricalFreq[k] = 1;
	}
	for (int i=0; i<Ntaxa; i++)	{
		for (int j=0; j<Nsite; j++)	{
			if (Data[i][j] != unknown)	{
				if ( (Data[i][j] < 0 ) || (Data[i][j] >= Nstate) )	{
					cerr << "error in data matrix : " << Data[i][j] << '\n';
					exit(1);
				}
				EmpiricalFreq[Data[i][j]] += 1.0;
			}
		}
	}
	double total = 0;
	for (int k = 0; k<Nstate; k++)	{
		total += EmpiricalFreq[k];
	}
	for (int k=0; k<Nstate; k++)	{
		EmpiricalFreq[k] /= total;
	}
}

// ---------------------------------------------------------------------------
//		 ConstantColumns()
// ---------------------------------------------------------------------------

int MCParameters::ConstantColumns()	{

	int n=0;
	for (int i=0; i<Nsite; i++)	{
		int j = 0;
		int test = 1;
		while ( (j<Ntaxa) && (Data[j][i] == unknown))	{
			j++;
		}
		if (j != Ntaxa)	{
			int k = j+1;
			while (k < Ntaxa)	{
				if (Data[k][i] != unknown)	{
					if (Data[j][i] != Data[k][i])	{
						test = 0;
					}
				}
				k++;
			}
		}
		if (test)	{
			n++;
		}
	}
	return n;
}


// ---------------------------------------------------------------------------
//		 EliminateConstantPositions()
// ---------------------------------------------------------------------------

void MCParameters::EliminateConstantPositions()	{

		cerr << "eliminate constant positions\n";
	int i=0;
	int j=0;
	int Eliminated = 0;
	while (i<Nsite)	{
		int k = 0;
		while ((k<Ntaxa) && (Data[k][i] == unknown)) k++;
		if (k<Ntaxa)	{
			int a = Data[k][i];
			k++;
			while ((k<Ntaxa) && ((Data[k][i] == unknown) || (Data[k][i] == a))) k++;
			if (k==Ntaxa)	{
				Eliminated ++;
				GeneSize[Gene[i]]--;
			}
			else	{
				for (int k=0; k<Ntaxa; k++)	{
					Data[k][j] = Data[k][i];
				}
				j++;
			}
		}
		i++;
	}

	Nsite -= Eliminated;
	cout << "number of positions eliminated : " << Eliminated << '\n';
	WriteDataToFile("woconst.puz");
}

// ---------------------------------------------------------------------------
//		 GetSequenceProfile(ostream& os)
// ---------------------------------------------------------------------------

void MCParameters::GetSequenceProfile(ostream& os, double pseudocount)	{

	os << Nsite << '\n';
	double* Stat = new double[Nstate];
	for (int i=0; i<Nsite; i++)	{

		for (int j=0; j<Nstate; j++)	{
			Stat[j] = pseudocount ;
		}

		for (int j=0; j<Ntaxa; j++)	{
			if (Data[j][i] == unknown)	{
				/*
				for (int k=0; k<Nstate; k++)	{
					Stat[k] += 1.0 / Nstate;
				}
				*/
			}
			else	{
				Stat[Data[j][i]] += 1.00;
			}
		}

		for (int j=0; j<Nstate; j++)	{
			os << Stat[j] << '\t';
		}
		os << '\n';

	}
	delete[] Stat;
}


// ---------------------------------------------------------------------------
//		 RegisterWithData()
// ---------------------------------------------------------------------------

void MCParameters::RegisterWithData()	{

	// test that no taxon name appears twice
      	for (int i=0; i<Ntaxa; i++)     {
                for (int j=i+1; j<Ntaxa; j++)   {
                        if (SpeciesNames[i] == SpeciesNames[j]) {
                                cerr << "error: taxa " << i << " and " << j << " in datafile have same name : " << SpeciesNames[i] << "\n";
                                exit(1);
                        }
                }
        }

	if (Ngene == 0)	{
		Ngene = 1;
		Gene = new int[Nsite];
		for (int i=0; i<Nsite; i++)	{
			Gene[i] =0;
		}
		GeneSize = new int[Ngene];
		GeneSize[0] = Nsite;
		GeneFirstSite = new int[Ngene];
		GeneFirstSite[0] = 0;
	}

	if (Recoding)	{
		RecodeData();
	}
	EliminateUnknownColumns();
	int total = 0;
	for (int i=0; i<Ngene; i++)	{
		GeneFirstSite[i] = total;
		int temp = GeneSize[i];
		for (int j=0; j<temp; j++)	{
			Gene[j+total] = i;
		}
		total += temp;
	}
	if (total != Nsite)	{
		cerr << "total number of sites does not match\n";
		exit(1);
	}
	if (DeleteConstant)	{
		EliminateConstantPositions();
		int total = 0;
		for (int i=0; i<Ngene; i++)	{
			GeneFirstSite[i] = total;
			int temp = GeneSize[i];
			for (int j=0; j<temp; j++)	{
				Gene[j+total] = i;
			}
			total += temp;
		}
		if (total != Nsite)	{
			cerr << "total number of sites does not match\n";
			exit(1);
		}
	}
	ComputeZipArrays();
	EstimateEmpiricalFrequencies();
	int nconst = 0;
	for (int i=0; i<Nsite; i++)	{
		if (OrbitSize[i] == 1)	{
			nconst ++;
		}
	}

	/*
	for (int j=0; j<Ntaxa; j++)	{
		int missing = 0;
		for (int i=0; i<Nsite; i++)	{
			if (Data[j][i] == unknown)	{
				missing++;
			}
		}
		cout << SpeciesNames[j] << '\t' << ((double) missing) / Nsite * 100 << '\n';
	}
	*/

	// cerr << "empirical proportion of invariable sites : " << ((double) nconst) / Nsite << '\n';

}

// ---------------------------------------------------------------------------
//		 ComputeZipArrays()
// ---------------------------------------------------------------------------

void MCParameters::ComputeZipArrays()	{

	EmpiricalFreq = new double[Nstate];
	for (int i=0; i<Nstate; i++)	{
		EmpiricalFreq[i] = 1.0 / Nstate;
	}
	SiteEmpiricalCount = new int*[Nsite];
	for (int i=0; i<Nsite; i++)	{
		SiteEmpiricalCount[i] = new int[Nstate];
		for (int j=0; j<Nstate; j++)	{
			SiteEmpiricalCount[i][j] = 0;
		}
		for (int j=0; j<Ntaxa; j++)	{
			if (Data[j][i] != unknown)	{
				SiteEmpiricalCount[i][Data[j][i]]++;
			}
		}
	}

	Indices = new int*[Nsite];
	ZipIndices = new int*[Nsite];
	for (int i=0; i<Nsite; i++)	{
		Indices[i] = new int[Nstate];
		ZipIndices[i] = new int[Nstate];
	}

	ZipSize = new int[Nsite];

	OrbitSize = new int[Nsite];
	Orbit = new int *[Nsite];
	for (int i=0; i<Nsite; i++)	{
		Orbit[i] = new int[Nstate];
	}

	ZipData = new int*[Ntaxa];
	for (int i=0; i<Ntaxa; i++)	{
		ZipData[i] = new int[Nsite];
	}

	for (int i=0; i<Nsite; i++)	{

		OrbitSize[i] = 0;
		for (int k=0 ; k< Nstate; k++)	{
			Orbit[i][k]= false;
		}

		for (int j=0; j<Ntaxa; j++)	{
			int d = Data[j][i];
			if (d != unknown)	{
				if (! Orbit[i][d])	{
					Orbit[i][d] = true;
					Indices[i][OrbitSize[i]] = d;
					OrbitSize[i] ++;
				}
			}
		}

		// sort Indices[i]
		for (int j=0; j<OrbitSize[i]; j++)	{
			for (int k=OrbitSize[i]-1; k>j; k--)	{
				if (Indices[i][j] > Indices[i][k])	{
					int tmp = Indices[i][j];
					Indices[i][j] = Indices[i][k];
					Indices[i][k] = tmp;
				}
			}
		}

		// reverse translation table
		for (int j=0; j<OrbitSize[i]; j++)	{
			ZipIndices[i][Indices[i][j]] = j;
		}
		for (int j=0; j<Nstate; j++)	{
			if (! Orbit[i][j])	{
				ZipIndices[i][j] = OrbitSize[i];
			}
		}	

		if (OrbitSize[i] == 0)	{
			cerr << "PhyloParameters::RegisterWithData : missing column";
			// exit(1);
		}

		if (OrbitSize[i] < Nstate)	{
			ZipSize[i] = OrbitSize[i] + 1;
		}
		else	{
			ZipSize[i] = OrbitSize[i];
		}

		for (int j=0; j<Ntaxa; j++)	{
			ZipData[j][i] = -2;
			int d = (Data)[j][i];
			if (d == unknown)	{
				ZipData[j][i] = unknown;
			}
			else	{
				for (int k=0; k<OrbitSize[i]; k++)	{
					if (Indices[i][k] == d)	{
						ZipData[j][i] = k;
					}
				}

				// Zip identity:
				// Indices[i][ZipData[j][i]] == mParam->Data[j][i] , for every i and j
				// within their respective range

			}

			// here, may be check that ZipData != -2

			if (ZipData[j][i] == -2)	{
				cerr << "PhyloParameters::RegisterWithData : error in zip data making\n";
			}
		}
		
		int observed[Nstate];
		int orbitsize= OrbitSize[i];
		for (int k=0; k<Nstate; k++)	{
			observed[k] = 0;
		}
		for (int k=0; k<orbitsize; k++)	{
			observed[Indices[i][k]] = 1;
		}
		for (int k=0; k<Nstate; k++)	{
			if (! observed[k])	{
				Indices[i][orbitsize++] = k;
			}
		}
		if (orbitsize != Nstate)	{
			cerr << "error in MCParameters::RegisterWithData\n";
			cerr << "site : " << i << '\n';
			/*
			for (int l=0; l<Nsite; l++)	{
				cerr << ":" << Alphabet[Data[l][i]] << ":";
			}
			*/
			cerr << '\n';
			exit(1);
		}
	}

	double temp = 0;
	for (int i=0; i<Nsite; i++)	{
		temp += ((double) ZipSize[i] * ZipSize[i]) / Nstate / Nstate;
	}
	SpeedFactor =  temp / Nsite ;
	// cerr << "speed factor : " << SpeedFactor << '\n';
	int nconst = 0;
	for (int i=0; i<Nsite; i++)	{
		if (OrbitSize[i] == 1)	{
			nconst++;
		}
	}
	// cerr << "number of constant columns: " << nconst << " (" << ((double) nconst) / Nsite * 100 << ")" << '\n';
	MeanOrbitSize = 0;
	for (int i=0; i<Nsite; i++)	{
		MeanOrbitSize += OrbitSize[i];
	}
	MeanOrbitSize /= Nsite;

}

void MCParameters::CompositionalCovariation()	{

	for (int j=0; j<Ntaxa; j++)	{
		int missing = 0;
		for (int i=0; i<Nsite; i++)	{
			if (Data[j][i] == unknown)	{
				missing++;
			}
		}
		cout << SpeciesNames[j] << '\t' << ((double) missing) / Nsite * 100 << '\n';
	}

	double taxfreq[Ntaxa][Nstate];
	double freq[Nstate];
	double cov[Nstate][Nstate];
	for (int j=0; j<Ntaxa; j++)	{
		for (int k=0; k<Nstate; k++)	{
			taxfreq[j][k] = 0;
		}
	}
	for (int i=0; i<Nstate; i++)	{
		freq[i] = 0;
	}
	for (int i=0; i<Nstate; i++)	{
		for (int j=0; j<Nstate; j++)	{
			cov[i][j] = 0;
		}
	}
	int count = 0;
	for (int j=0; j<Ntaxa; j++)	{
		int taxcount = 0;
		for (int i=0; i<Nsite; i++)	{
			if (Data[j][i] != unknown)	{
				taxcount++;
				taxfreq[j][Data[j][i]]++;
				count++;
				freq[Data[j][i]]++;
			}
		}
		for (int k=0; k<Nstate; k++)	{
			taxfreq[j][k] /= taxcount;
		}
	}
	for (int k=0; k<Nstate; k++)	{
		freq[k] /= count;
	}
	
	for (int i=0; i<Nstate; i++)	{
		for (int j=i; j<Nstate; j++)	{
			for (int k=0; k<Ntaxa; k++)	{
				cov[i][j] += (taxfreq[k][i] - freq[i]) * (taxfreq[k][j] - freq[j]);
			}
			cov[i][j] /= Ntaxa;
		}
	}
	/*
	for (int i=0; i<Nstate; i++)	{
		for (int j=i; j<Nstate; j++)	{
			cout << i << '\t' << j << '\t' << Alphabet[i] << '\t' << Alphabet[j] << '\t' << cov[i][j] << '\t' << cov[i][j] / sqrt(cov[i][i]*cov[j][j]) << '\n';
		}
	}
	*/
	cout << '\n';
	cout << '\t';
	for (int i=0; i<Nstate; i++)	{
		cout << Alphabet[i] << '\t';
	}
	cout << '\n';
	cout << '\n';
	for (int i=0; i<Nstate; i++)	{
		cout << Alphabet[i] << '\t';
		for (int j=0; j<i; j++)	{
			cout << '\t';
		}
		for (int j=i; j<Nstate; j++)	{
			cout << (int) (1000000 *  cov[i][j]) << '\t';
			// cout << (int) (100 *  cov[i][j] / sqrt(cov[i][i]*cov[j][j])) << '\t';
		}
		cout << Alphabet[i] << '\n';
	}

}



// ---------------------------------------------------------------------------------
//		 Streams
// ---------------------------------------------------------------------------------


ostream& operator<<( ostream& os , MCParameters& param)	{

	// data file
	os << param.Version << '\n';
	os << param.SubVersion << '\n';
	os << param.DataFileSpec << '\n';
	os << param.ContDataFileSpec << '\n';

	// model specifications

	os << param.NormalApprox << '\n';
	os << param.NormalConcMode << '\n';
	os << param.ClockModel << '\n';
	os << param.WithPartial << '\n';
	os << param.FlexClockModel << '\n';
	os << param.ActivateClock << '\n';
	os << param.SeparateRhoPrior << '\n';

	os << param.NH << '\n';
	os << param.DeleteConstant << '\n';
	os << param.Normalise << '\n';
	os << param.NormaliseCov << '\n';
	os << param.ModeFastCompute << '\n';

	os << param.ZipGTR << '\n';
	os << param.ZipPrior << '\n';
	os << param.ZipGTRDP << '\n';
	os << param.NZipModeMax << '\n';
	os << param.ModePoisson << '\n';
	if (param.PopEff)	{
		os << 10 << '\n';
	}
	else	{
		os << param.NHPrior << '\n';
	}

	// specific to version 3.3
	os << param.PopEff << '\n';
	os << param.NHPop << '\n';
	os << param.FixRoot << '\n';

	os << param.MBL << '\n';
	os << param.NMixBLMax << '\n';
	os << param.IncrementalMixBLDP << '\n';

	os << param.RefFastCompute << '\n';
	os << param.Parallel << '\n';
	os << param.NHNcatMax << '\n';
	os << param.ExternalCov << '\n';
	os << param.Ngen << '\n';
	os << param.SaveStat << '\n';
	os << param.SaveAll << '\n';
	os << param.NmodeMax << '\n';
	os << param.NRateModeMax << '\n';
	os << param.GammaNcat << '\n';
	os << param.ActivateSumOverModes << '\n';
	os << param.SumOverModes << '\n';
	os << param.ActivateSumOverRateModes << '\n';
	os << param.SumOverRateModes << '\n';
	os << param.Qmode << '\n';
	os << param.GeneBLMode << '\n';
	os << param.GeneBLMultiplier << '\n';
	os << param.GeneRateMode << '\n';
	os << param.GeneGammaMode << '\n';
	os << param.GeneStationaryMode << '\n';

	os << param.LengthPrior << '\n';
	os << param.LengthMin << '\n';
	os << param.LengthMax << '\n';

	os << param.RatePrior << '\n';
	os << param.RateMin << '\n';
	os << param.RateMax << '\n';
	os << param.GammaMin << '\n';
	os << param.GammaMax << '\n';
	os << param.XiMin << '\n';
	os << param.XiMax << '\n';
	os << param.GammaPrior << '\n';

	os << param.ModePrior << '\n';
	os << param.ModeStatPrior << '\n';
	os << param.StatMin << '\n';
	os << param.StatMax << '\n';
	os << param.StatAlphaMin << '\n';
	os << param.StatAlphaMax << '\n';

	os << param.AlphaPrior << '\n';
	os << param.AlphaMin << '\n';
	os << param.AlphaMax << '\n';

	os << param.TooSmall << '\n';
	os << param.TooLarge << '\n';
	os << param.InfProb << '\n';
	os << param.OutNtaxa << '\n';
	for (int i=0; i<param.OutNtaxa; i++)	{
		os << param.OutGroup[i] << '\n';
	}
	if (param.NormalApprox == Yes)	{
	       	os << param.ConcCov << '\n';
		if (param.ConcCov)	{
			os << param.ConcCovFileName << '\n';
		}
		os << param.SepCov << '\n';
		if (param.SepCov)	{
			os << param.SepCovFileName << '\n';
		}
	}
	else	{
		os << param.Recoding << '\n';
		if (param.Recoding)	{
			os << param.RawNstate << '\n';
			os << param.Nstate << '\n';
			for (int i=0; i<param.RawNstate; i++)	{
				os << param.RecodingTable[i] << '\t';
			}
			os << '\n';
			for (int i=0; i<param.Nstate; i++)	{
				os << param.RecodingAlphabet[i] << '\t';
			}
			os << '\n';
		}

		os << param.Ngene << '\n';
		for (int i=0; i<param.Ngene; i++)	{
			os << param.GeneSize[i] << '\t';
			os << param.GeneFirstSite[i] << '\n';
		}
		os << '\n';
		for (int i=0; i<param.Nsite; i++)	{
			os << param.Gene[i] << '\t' ;
		}
		os << '\n';

	}
		
	os << param.NCalib << '\n';
	if (param.NCalib)	{
		for (int i=0; i<param.NCalib; i++)	{
			os << param.CalibTaxon1[i] << '\t' << param.CalibTaxon2[i] << '\t' << param.CalibIndex[i] << '\t' << param.CalibUpper[i] << '\t' << param.CalibLower[i] << '\n';
		}
	}

	os << param.SaveEvery << '\n';
	os << param.StopAfter << '\n';
	os << param.HowManySaved << '\n';
	os << param.SwapFreq << '\n';

	os << param.BurnIn << '\n';
	os << param.InitBeta << '\n';
	os << param.FinalBeta << '\n';
	os << param.BetaStep << '\n';
	os << param.BurnInDone << '\n';
	os << param.Beta0 << '\n';
	os << param.Nchain << '\n';
	for (int l=0; l<param.Nchain; l++)	{
		os << param.Beta[l] << '\n';
	}

	for (int i=0; i<param.Nchain; i++)	{
		os << param.SwapPermut[i] << '\n';
	}

	os << param.MSMode << '\n';
	os << param.RASModelSwitch << '\n';
	os << param.SUBModelSwitch << '\n';
	os << param.SeparateModelSwitch << '\n';
	os << param.HeteroModelSwitch << '\n';
	os << param.FlexClockModelSwitch << '\n';
	os << param.BLModelSwitch << '\n';

	os << param.HeteroMode << '\n';

	// chain parameters

	os << param.TopoMoveTypeNumber << '\n';
	for (int k=0; k<param.TopoMoveTypeNumber; k++)	{
		os << param.TopoMoveTypeArray[k] << '\n';
		os << param.SuperMoveTypeNumber[k] << '\n';
		for (int l=0; l<param.Nchain; l++)	{
			os << param.TopoNIterationArray[k][l] << '\n';
			os << param.TopoNArray[k][l] << '\n';
		}
		for (int i=0; i<param.SuperMoveTypeNumber[k]; i++)	{
			os << param.MoveTypeNumber[k][i] << '\n';
			for (int l=0; l<param.Nchain; l++)	{
				os << param.SuperNIterationArray[k][i][l] << '\n';
			}
			for (int j=0; j<param.MoveTypeNumber[k][i]; j++)	{
				os << param.MoveTypeArray[k][i][j] << '\n';
				for (int l=0; l<param.Nchain; l++)	{
					os << param.deltaArray[k][i][j][l] << '\n';
					os << param.NIterationArray[k][i][j][l] << '\n';
					os << param.NArray[k][i][j][l] << '\n';
				}
			}
		}
	}

	// chain statistics

	os << param.RateInfCount << '\n';
	os << param.StatInfCount << '\n';
	os << param.LengthInfCount << '\n';
	os << param.LogProbInfCount << '\n';
	os << param.SubOverflowCount << '\n';

	for (int k=0; k<param.TopoMoveTypeNumber; k++)	{
		for (int i=0; i<param.SuperMoveTypeNumber[k]; i++)	{
			for (int l=0; l<param.MoveTypeNumber[k][i]; l++)	{
				for (int m=0; m<param.Nchain; m++)	{
					os << param.NCallArray[k][i][l][m]<< '\t';
					os << param.SuccessArray[k][i][l][m]<< '\t';
					os << param.TimeArray[k][i][l][m] << '\n';
				}
			}
		}
	}
	for (int m=0; m<param.Nchain; m++)	{
		os << param.TrialSwapArray[m] << '\t';
		os << param.AcceptedSwapArray[m] << '\n';
	}

	os << param.FixGamma << '\n';
	os << param.FixMeanLength << '\n';
	os << param.FixStatCenter << '\n';
	os << param.FixStat << '\n';
	os << param.FixRR << '\n';
	os << param.FixLength << '\n';
	os << param.FixRate << '\n';
	os << param.FixNmode << '\n';
	// os << param.FixNRateMode << '\n';
	os << param.FixTopo << '\n';
	os << param.FixPconst << '\n';

	os << param.TimePrior << '\n';
	os << param.ClockPriorModelSwitch << '\n';
	os << param.AutoCorrelModelSwitch << '\n';
	os << param.ArcThermo << '\n';
	os << param.alpha1 << '\n';
	os << param.alpha2 << '\n';
	os << param.Zeta << '\n';
	os << param.logBF << '\n';
	os << param.logBFmin << '\n';
	os << param.logBFmax << '\n';
	os << param.logBFrev << '\n';
	os << param.logBFminrev << '\n';
	os << param.logBFmaxrev << '\n';
	os << param.BFcount << '\n';
	os << param.ScalePrior << '\n';
	os << param.MeanScale << "\n";
	os << param.MutMode << '\n';
	os << param.SelMode << '\n';

	os << param.SoftBounds << '\n';
	os << param.Softa << '\n';
	os << param.VarScale << '\n';

	os << param.ImproperLowerBound << '\n';
	os << param.LowerC << '\n';
	os << param.LowerP << '\n';

	os << param.AlphaChi << '\n';
	os << param.BetaChi << '\n';
	os << param.AlphaChi2 << '\n';
	os << param.BetaChi2 << '\n';

	os << param.ZipSub << '\n';

	os << param.LengthGammaPrior << '\n';

	os << 1 << '\n';
	os << param.ModeStatAlphaPrior << '\n';

	os << 0 << '\n';

	int tmp = param.SaveAll;
	param.SaveAll = 1;
	for (int chain = 0; chain < param.Nchain; chain++)	{
		os << *(param.currentState[chain]);
	}
	param.SaveAll = tmp;

	return os;
}


istream& operator>>( istream& is , MCParameters& param)	{

	is >> param.Version;
	is >> param.SubVersion;
	if (param.Version < 2)	{
		cerr << "error: version 1 is not supported by current version of phylobayes\n";
		cerr << "sorry\n";
		exit(1);
	}
	if ((param.Version == 3) && (param.SubVersion == 0))	{
		// data file
		is >> param.DataFileSpec;
		is >> param.ContDataFileSpec;


		// model specifications
		is >> param.NormalApprox;
		is >> param.NormalConcMode;
		is >> param.ClockModel;
		is >> param.WithPartial;
		is >> param.FlexClockModel;
		is >> param.ActivateClock;
		is >> param.SeparateRhoPrior;

		is >> param.NH;
		is >> param.DeleteConstant;
		is >> param.Normalise;
		is >> param.NormaliseCov;
		is >> param.ModeFastCompute;

		is >> param.ZipGTR;
		is >> param.ZipPrior;
		is >> param.ZipGTRDP;
		is >> param.NZipModeMax;
		is >> param.ModePoisson;
		is >> param.NHPrior;
		if (param.NHPrior == 10)	{
			param.NHPrior = 0;
			param.PopEff = 1;
		}
		is >> param.MBL;
		is >> param.NMixBLMax;
		is >> param.IncrementalMixBLDP;

		is >> param.RefFastCompute;
		is >> param.Parallel;
		is >> param.NHNcatMax;
		is >> param.ExternalCov;
		is >> param.Ngen;
		is >> param.SaveStat;
		is >> param.SaveAll;

		is >> param.NmodeMax;
		is >> param.NRateModeMax;
		is >> param.GammaNcat;
		is >> param.ActivateSumOverModes;
		is >> param.SumOverModes;
		is >> param.ActivateSumOverRateModes;
		is >> param.SumOverRateModes;
		is >> param.Qmode;
		is >> param.GeneBLMode;
		is >> param.GeneBLMultiplier;
		is >> param.GeneRateMode;
		is >> param.GeneGammaMode;
		is >> param.GeneStationaryMode;

		is >> param.LengthPrior;
		is >> param.LengthMin;
		is >> param.LengthMax;

		is >> param.RatePrior;
		is >> param.RateMin;
		is >> param.RateMax;
		is >> param.GammaMin;
		is >> param.GammaMax;
		is >> param.XiMin;
		is >> param.XiMax;

		is >> param.ModePrior;
		is >> param.ModeStatPrior;
		is >> param.StatMin;
		is >> param.StatMax;
		is >> param.StatAlphaMin;
		is >> param.StatAlphaMax;

		is >> param.AlphaPrior;
		is >> param.AlphaMin;
		is >> param.AlphaMax;

		is >> param.TooSmall;
		is >> param.TooLarge;
		is >> param.InfProb;
		is >> param.OutNtaxa;
		if (param.OutNtaxa)	{
			param.OutGroup = new string[param.OutNtaxa];
			for (int i=0; i<param.OutNtaxa; i++)	{
				is >> param.OutGroup[i];
			}
		}
		if (param.NormalApprox == Yes)	{
			int tmp;
			is >> tmp;
			if (tmp)	{
				string temp;
				is >> temp;
				param.ReadCov(temp);
			}
			is >> tmp;
			if (tmp)	{
				string temp;
				is >> temp;
				param.ReadSepCov(temp);
			}
			param.RegisterCov();
		}
		else	{
			param.ReadDataFromFile(param.DataFileSpec);

			is >> param.Recoding;
			if (param.Recoding)	{
				is >> param.Nstate;
				is >> param.RecodingNstate;
				param.RecodingTable = new int[param.Nstate];
				for (int i=0; i<param.Nstate; i++)	{
					is >> param.RecodingTable[i];
				}
				param.RecodingAlphabet = new char[param.RecodingNstate];
				for (int i=0; i<param.RecodingNstate; i++)	{
					is >> param.RecodingAlphabet[i];
				}
			}

			param.RegisterWithData();

			is >> param.Ngene;
			delete[] param.Gene;
			delete[] param.GeneSize;
			delete[] param.GeneFirstSite;
			param.GeneSize = new int[param.Ngene];
			param.GeneFirstSite = new int[param.Ngene];
			param.Gene = new int[param.Nsite];

			for (int i=0; i<param.Ngene; i++)	{
				is >> param.GeneSize[i];
				is >> param.GeneFirstSite[i];
			}
			for (int i=0; i<param.Nsite; i++)	{
				is >> param.Gene[i];
			}
			if (param.ContDataFileSpec != "None")	{
				param.ReadContDataFromFile(param.ContDataFileSpec);
			}	
		}

		is >> param.NCalib;
		if (param.NCalib)	{
			param.isCalibrated = new int[param.Nnode];
			for (int j=0; j<param.Nnode; j++)	{
				param.isCalibrated[j] = 0;
			}
			param.CalibTaxon1 = new string[param.NCalib];
			param.CalibTaxon2 = new string[param.NCalib];
			param.CalibIndex = new int[param.NCalib];
			param.CalibUpper = new double[param.NCalib];
			param.CalibLower = new double[param.NCalib];
			for (int i=0; i<param.NCalib; i++)	{
				is >> param.CalibTaxon1[i] >> param.CalibTaxon2[i] >> param.CalibIndex[i] >> param.CalibUpper[i] >> param.CalibLower[i];
				param.isCalibrated[param.CalibIndex[i]] = 1;
			}
		}

		is >> param.SaveEvery;
		is >> param.StopAfter;
		is >> param.HowManySaved;
		is >> param.SwapFreq;

		is >> param.BurnIn;
		is >> param.InitBeta;
		is >> param.FinalBeta;
		is >> param.BetaStep;
		is >> param.BurnInDone;
		is >> param.Beta0;
		is >> param.Nchain;

		for (int l=0; l<param.Nchain; l++)	{
			is >> param.Beta[l];
		}

		param.SwapPermut = new int[param.Nchain];
		for (int i=0; i<param.Nchain; i++)	{
			is >> param.SwapPermut[i];
		}

		is >> param.MSMode;
		is >> param.RASModelSwitch;
		is >> param.SUBModelSwitch;
		is >> param.SeparateModelSwitch;
		is >> param.HeteroModelSwitch;
		is >> param.FlexClockModelSwitch;
		is >> param.BLModelSwitch;

		is >> param.HeteroMode;

		// chain parameters

		is >> param.TopoMoveTypeNumber;
		for (int k=0; k<param.TopoMoveTypeNumber; k++)	{
			is >> param.TopoMoveTypeArray[k];
			is >> param.SuperMoveTypeNumber[k];
			for (int l=0; l<param.Nchain; l++)	{
				is >> param.TopoNIterationArray[k][l];
				is >> param.TopoNArray[k][l];
			}
			for (int i=0; i<param.SuperMoveTypeNumber[k]; i++)	{
				is >> param.MoveTypeNumber[k][i];
				for (int l=0; l<param.Nchain; l++)	{
					is >> param.SuperNIterationArray[k][i][l];
				}
				for (int j=0; j<param.MoveTypeNumber[k][i]; j++)	{
					is >> param.MoveTypeArray[k][i][j];
					for (int l=0; l<param.Nchain; l++)	{
						is >> param.deltaArray[k][i][j][l];
						is >> param.NIterationArray[k][i][j][l];
						is >> param.NArray[k][i][j][l];
					}
				}
			}
		}

		// chain statistics

		is >> param.RateInfCount;
		is >> param.StatInfCount;
		is >> param.LengthInfCount;
		is >> param.LogProbInfCount;
		is >> param.SubOverflowCount;


		for (int k=0; k<param.TopoMoveTypeNumber; k++)	{
			for (int i=0; i<param.SuperMoveTypeNumber[k]; i++)	{
				for (int l=0; l<param.MoveTypeNumber[k][i]; l++)	{
					for (int m=0; m<param.Nchain; m++)	{
						is >> param.NCallArray[k][i][l][m];
						is >> param.SuccessArray[k][i][l][m];
						is >> param.TimeArray[k][i][l][m];
					}
				}
			}
		}
		for (int m=0; m<param.Nchain; m++)	{
			is >> param.TrialSwapArray[m];
			is >> param.AcceptedSwapArray[m];
		}


		is >> param.FixGamma;
		is >> param.FixMeanLength;
		is >> param.FixStatCenter;
		is >> param.FixStat;
		is >> param.FixRR;
		is >> param.FixLength;
		is >> param.FixRate;
		is >> param.FixNmode;
		is >> param.FixTopo;
		is >> param.FixPconst;
		is >> param.TimePrior;
		is >> param.ClockPriorModelSwitch;
		is >> param.AutoCorrelModelSwitch;
		is >> param.ArcThermo;
		is >> param.alpha1;
		is >> param.alpha2;
		is >> param.Zeta;
		is >> param.logBF;
		is >> param.logBFmin;
		is >> param.logBFmax;
		is >> param.logBFrev;
		is >> param.logBFminrev;
		is >> param.logBFmaxrev;
		is >> param.BFcount;
		is >> param.ScalePrior;
		is >> param.MeanScale;
		is >> param.MutMode;
		is >> param.SelMode;

		int cont;
		is >> cont;
		if (cont)	{
			cerr << "error in MCParameters::icstream\n";
			exit(1);
		}

		// current state
		int tmp = param.SaveAll;
		param.SaveAll = 1;
		param.currentState = new PhyloBayes*[param.Nchain];
		param.nextState = new PhyloBayes*[param.Nchain];
		for (int chain = 0; chain < param.Nchain; chain++)	{
			param.currentState[chain] = new PhyloBayes(&param);
			param.nextState[chain] = new PhyloBayes(&param);
			is >> *param.currentState[chain];
		}
		param.SaveAll = tmp;

		return is;
	}
	else if ((param.Version == 3) && (param.SubVersion == 1))	{

		// data file
		is >> param.DataFileSpec;
		is >> param.ContDataFileSpec;


		// model specifications
		is >> param.NormalApprox;
		is >> param.NormalConcMode;
		is >> param.ClockModel;
		is >> param.WithPartial;
		is >> param.FlexClockModel;
		is >> param.ActivateClock;
		is >> param.SeparateRhoPrior;

		is >> param.NH;
		is >> param.DeleteConstant;
		is >> param.Normalise;
		is >> param.NormaliseCov;
		is >> param.ModeFastCompute;

		is >> param.ZipGTR;
		is >> param.ZipPrior;
		is >> param.ZipGTRDP;
		is >> param.NZipModeMax;
		is >> param.ModePoisson;
		is >> param.NHPrior;
		if (param.NHPrior == 10)	{
			param.NHPrior = 0;
			param.PopEff = 1;
		}
		is >> param.MBL;
		is >> param.NMixBLMax;
		is >> param.IncrementalMixBLDP;

		is >> param.RefFastCompute;
		is >> param.Parallel;
		is >> param.NHNcatMax;
		is >> param.ExternalCov;
		is >> param.Ngen;
		is >> param.SaveStat;
		is >> param.SaveAll;

		is >> param.NmodeMax;
		is >> param.NRateModeMax;
		is >> param.GammaNcat;
		is >> param.ActivateSumOverModes;
		is >> param.SumOverModes;
		is >> param.ActivateSumOverRateModes;
		is >> param.SumOverRateModes;
		is >> param.Qmode;
		is >> param.GeneBLMode;
		is >> param.GeneBLMultiplier;
		is >> param.GeneRateMode;
		is >> param.GeneGammaMode;
		is >> param.GeneStationaryMode;

		is >> param.LengthPrior;
		is >> param.LengthMin;
		is >> param.LengthMax;

		is >> param.RatePrior;
		is >> param.RateMin;
		is >> param.RateMax;
		is >> param.GammaMin;
		is >> param.GammaMax;
		is >> param.XiMin;
		is >> param.XiMax;

		is >> param.ModePrior;
		is >> param.ModeStatPrior;
		is >> param.StatMin;
		is >> param.StatMax;
		is >> param.StatAlphaMin;
		is >> param.StatAlphaMax;

		is >> param.AlphaPrior;
		is >> param.AlphaMin;
		is >> param.AlphaMax;

		is >> param.TooSmall;
		is >> param.TooLarge;
		is >> param.InfProb;
		is >> param.OutNtaxa;
		if (param.OutNtaxa)	{
			param.OutGroup = new string[param.OutNtaxa];
			for (int i=0; i<param.OutNtaxa; i++)	{
				is >> param.OutGroup[i];
			}
		}
		if (param.NormalApprox == Yes)	{
			int tmp;
			is >> tmp;
			if (tmp)	{
				string temp;
				is >> temp;
				param.ReadCov(temp);
			}
			is >> tmp;
			if (tmp)	{
				string temp;
				is >> temp;
				param.ReadSepCov(temp);
			}
			param.RegisterCov();
		}
		else	{
			param.ReadDataFromFile(param.DataFileSpec);

			is >> param.Recoding;
			if (param.Recoding)	{
				is >> param.Nstate;
				is >> param.RecodingNstate;
				param.RecodingTable = new int[param.Nstate];
				for (int i=0; i<param.Nstate; i++)	{
					is >> param.RecodingTable[i];
				}
				param.RecodingAlphabet = new char[param.RecodingNstate];
				for (int i=0; i<param.RecodingNstate; i++)	{
					is >> param.RecodingAlphabet[i];
				}
			}

			param.RegisterWithData();

			is >> param.Ngene;
			delete[] param.Gene;
			delete[] param.GeneSize;
			delete[] param.GeneFirstSite;
			param.GeneSize = new int[param.Ngene];
			param.GeneFirstSite = new int[param.Ngene];
			param.Gene = new int[param.Nsite];

			for (int i=0; i<param.Ngene; i++)	{
				is >> param.GeneSize[i];
				is >> param.GeneFirstSite[i];
			}
			for (int i=0; i<param.Nsite; i++)	{
				is >> param.Gene[i];
			}
			if (param.ContDataFileSpec != "None")	{
				param.ReadContDataFromFile(param.ContDataFileSpec);
			}	
		}

		is >> param.NCalib;
		if (param.NCalib)	{
			param.isCalibrated = new int[param.Nnode];
			for (int j=0; j<param.Nnode; j++)	{
				param.isCalibrated[j] = 0;
			}
			param.CalibTaxon1 = new string[param.NCalib];
			param.CalibTaxon2 = new string[param.NCalib];
			param.CalibIndex = new int[param.NCalib];
			param.CalibUpper = new double[param.NCalib];
			param.CalibLower = new double[param.NCalib];
			for (int i=0; i<param.NCalib; i++)	{
				is >> param.CalibTaxon1[i] >> param.CalibTaxon2[i] >> param.CalibIndex[i] >> param.CalibUpper[i] >> param.CalibLower[i];
				/*
				if (param.CalibIndex[i] <0)	{
					cerr << "negative index\n";
					exit(1);
				}
				if (param.CalibIndex[i] >= param.Nnode)	{
					cerr << "index overflow\n";
					exit(1);
				}
				cerr << param.CalibIndex[i] << '\n';
				cerr << "ok\n";
				*/
				param.isCalibrated[param.CalibIndex[i]] = 1;
			}
		}

		is >> param.SaveEvery;
		is >> param.StopAfter;
		is >> param.HowManySaved;
		is >> param.SwapFreq;

		is >> param.BurnIn;
		is >> param.InitBeta;
		is >> param.FinalBeta;
		is >> param.BetaStep;
		is >> param.BurnInDone;
		is >> param.Beta0;
		is >> param.Nchain;

		for (int l=0; l<param.Nchain; l++)	{
			is >> param.Beta[l];
		}

		param.SwapPermut = new int[param.Nchain];
		for (int i=0; i<param.Nchain; i++)	{
			is >> param.SwapPermut[i];
		}

		is >> param.MSMode;
		is >> param.RASModelSwitch;
		is >> param.SUBModelSwitch;
		is >> param.SeparateModelSwitch;
		is >> param.HeteroModelSwitch;
		is >> param.FlexClockModelSwitch;
		is >> param.BLModelSwitch;

		is >> param.HeteroMode;

		// chain parameters

		is >> param.TopoMoveTypeNumber;
		for (int k=0; k<param.TopoMoveTypeNumber; k++)	{
			is >> param.TopoMoveTypeArray[k];
			is >> param.SuperMoveTypeNumber[k];
			for (int l=0; l<param.Nchain; l++)	{
				is >> param.TopoNIterationArray[k][l];
				is >> param.TopoNArray[k][l];
			}
			for (int i=0; i<param.SuperMoveTypeNumber[k]; i++)	{
				is >> param.MoveTypeNumber[k][i];
				for (int l=0; l<param.Nchain; l++)	{
					is >> param.SuperNIterationArray[k][i][l];
				}
				for (int j=0; j<param.MoveTypeNumber[k][i]; j++)	{
					is >> param.MoveTypeArray[k][i][j];
					for (int l=0; l<param.Nchain; l++)	{
						is >> param.deltaArray[k][i][j][l];
						is >> param.NIterationArray[k][i][j][l];
						is >> param.NArray[k][i][j][l];
					}
				}
			}
		}

		// chain statistics

		is >> param.RateInfCount;
		is >> param.StatInfCount;
		is >> param.LengthInfCount;
		is >> param.LogProbInfCount;
		is >> param.SubOverflowCount;


		for (int k=0; k<param.TopoMoveTypeNumber; k++)	{
			for (int i=0; i<param.SuperMoveTypeNumber[k]; i++)	{
				for (int l=0; l<param.MoveTypeNumber[k][i]; l++)	{
					for (int m=0; m<param.Nchain; m++)	{
						is >> param.NCallArray[k][i][l][m];
						is >> param.SuccessArray[k][i][l][m];
						is >> param.TimeArray[k][i][l][m];
					}
				}
			}
		}
		for (int m=0; m<param.Nchain; m++)	{
			is >> param.TrialSwapArray[m];
			is >> param.AcceptedSwapArray[m];
		}


		is >> param.FixGamma;
		is >> param.FixMeanLength;
		is >> param.FixStatCenter;
		is >> param.FixStat;
		is >> param.FixRR;
		is >> param.FixLength;
		is >> param.FixRate;
		is >> param.FixNmode;
		is >> param.FixTopo;
		is >> param.FixPconst;
		is >> param.TimePrior;
		is >> param.ClockPriorModelSwitch;
		is >> param.AutoCorrelModelSwitch;
		is >> param.ArcThermo;
		is >> param.alpha1;
		is >> param.alpha2;
		is >> param.Zeta;
		is >> param.logBF;
		is >> param.logBFmin;
		is >> param.logBFmax;
		is >> param.logBFrev;
		is >> param.logBFminrev;
		is >> param.logBFmaxrev;
		is >> param.BFcount;
		is >> param.ScalePrior;
		is >> param.MeanScale;
		is >> param.MutMode;
		is >> param.SelMode;

		int cont;
		is >> cont;
		if (cont)	{
			is >> param.SoftBounds;
			is >> param.Softa;
			is >> param.VarScale;
			int cont;
			is >> cont;
			if (cont)	{
				is >> param.ImproperLowerBound;
				is >>  param.LowerC;
				is >>  param.LowerP;

				is >> param.AlphaChi;
				is >> param.BetaChi;
				is >> param.AlphaChi2;
				is >> param.BetaChi2;

				int cont;
				is >> cont;
				if (cont)	{
					is >> param.ZipSub;
					is >> cont;
					if (cont)	{
						cerr << "error in MCParameters::icstream\n";
						exit(1);
					}
				}
			}
		}

		// current state
		int tmp = param.SaveAll;
		param.SaveAll = 1;
		param.currentState = new PhyloBayes*[param.Nchain];
		param.nextState = new PhyloBayes*[param.Nchain];
		for (int chain = 0; chain < param.Nchain; chain++)	{
			param.currentState[chain] = new PhyloBayes(&param);
			param.nextState[chain] = new PhyloBayes(&param);
			is >> *param.currentState[chain];
		}
		param.SaveAll = tmp;
	}
	else if (param.Version >= 4)	{

		// data file
		is >> param.DataFileSpec;
		is >> param.ContDataFileSpec;


		// model specifications
		is >> param.NormalApprox;
		is >> param.NormalConcMode;
		is >> param.ClockModel;
		is >> param.WithPartial;
		is >> param.FlexClockModel;
		is >> param.ActivateClock;
		is >> param.SeparateRhoPrior;

		is >> param.NH;
		is >> param.DeleteConstant;
		is >> param.Normalise;
		is >> param.NormaliseCov;
		is >> param.ModeFastCompute;

		is >> param.ZipGTR;
		is >> param.ZipPrior;
		is >> param.ZipGTRDP;
		is >> param.NZipModeMax;
		is >> param.ModePoisson;
		is >> param.NHPrior;
		if (param.NHPrior == 10)	{
			param.NHPrior = 0;
			param.PopEff = 1;
		}

		// specific to version 3.3
		is >> param.PopEff;
		is >> param.NHPop;
		is >> param.FixRoot;

		is >> param.MBL;
		is >> param.NMixBLMax;
		is >> param.IncrementalMixBLDP;

		is >> param.RefFastCompute;
		is >> param.Parallel;
		is >> param.NHNcatMax;
		is >> param.ExternalCov;
		is >> param.Ngen;
		is >> param.SaveStat;
		is >> param.SaveAll;

		is >> param.NmodeMax;
		is >> param.NRateModeMax;
		is >> param.GammaNcat;
		is >> param.ActivateSumOverModes;
		is >> param.SumOverModes;
		is >> param.ActivateSumOverRateModes;
		is >> param.SumOverRateModes;
		is >> param.Qmode;
		is >> param.GeneBLMode;
		is >> param.GeneBLMultiplier;
		is >> param.GeneRateMode;
		is >> param.GeneGammaMode;
		is >> param.GeneStationaryMode;

		is >> param.LengthPrior;
		is >> param.LengthMin;
		is >> param.LengthMax;

		is >> param.RatePrior;
		is >> param.RateMin;
		is >> param.RateMax;
		is >> param.GammaMin;
		is >> param.GammaMax;
		is >> param.XiMin;
		is >> param.XiMax;
		is >> param.GammaPrior;

		is >> param.ModePrior;
		is >> param.ModeStatPrior;
		is >> param.StatMin;
		is >> param.StatMax;
		is >> param.StatAlphaMin;
		is >> param.StatAlphaMax;

		is >> param.AlphaPrior;
		is >> param.AlphaMin;
		is >> param.AlphaMax;

		is >> param.TooSmall;
		is >> param.TooLarge;
		is >> param.InfProb;
		is >> param.OutNtaxa;
		if (param.OutNtaxa)	{
			param.OutGroup = new string[param.OutNtaxa];
			for (int i=0; i<param.OutNtaxa; i++)	{
				is >> param.OutGroup[i];
			}
		}
		if (param.NormalApprox == Yes)	{
			int tmp;
			is >> tmp;
			if (tmp)	{
				string temp;
				is >> temp;
				param.ReadCov(temp);
			}
			is >> tmp;
			if (tmp)	{
				string temp;
				is >> temp;
				param.ReadSepCov(temp);
			}
			param.RegisterCov();
		}
		else	{
			param.ReadDataFromFile(param.DataFileSpec);

			is >> param.Recoding;
			if (param.Recoding)	{
				is >> param.Nstate;
				is >> param.RecodingNstate;
				param.RecodingTable = new int[param.Nstate];
				for (int i=0; i<param.Nstate; i++)	{
					is >> param.RecodingTable[i];
				}
				param.RecodingAlphabet = new char[param.RecodingNstate];
				for (int i=0; i<param.RecodingNstate; i++)	{
					is >> param.RecodingAlphabet[i];
				}
			}

			param.RegisterWithData();

			is >> param.Ngene;
			delete[] param.Gene;
			delete[] param.GeneSize;
			delete[] param.GeneFirstSite;
			param.GeneSize = new int[param.Ngene];
			param.GeneFirstSite = new int[param.Ngene];
			param.Gene = new int[param.Nsite];

			for (int i=0; i<param.Ngene; i++)	{
				is >> param.GeneSize[i];
				is >> param.GeneFirstSite[i];
			}
			for (int i=0; i<param.Nsite; i++)	{
				is >> param.Gene[i];
			}
			if (param.ContDataFileSpec != "None")	{
				param.ReadContDataFromFile(param.ContDataFileSpec);
			}	
		}

		is >> param.NCalib;
		if (param.NCalib)	{
			param.isCalibrated = new int[param.Nnode];
			for (int j=0; j<param.Nnode; j++)	{
				param.isCalibrated[j] = 0;
			}
			param.CalibTaxon1 = new string[param.NCalib];
			param.CalibTaxon2 = new string[param.NCalib];
			param.CalibIndex = new int[param.NCalib];
			param.CalibUpper = new double[param.NCalib];
			param.CalibLower = new double[param.NCalib];
			for (int i=0; i<param.NCalib; i++)	{
				is >> param.CalibTaxon1[i] >> param.CalibTaxon2[i] >> param.CalibIndex[i] >> param.CalibUpper[i] >> param.CalibLower[i];
				/*
				if (param.CalibIndex[i] <0)	{
					cerr << "negative index\n";
					exit(1);
				}
				if (param.CalibIndex[i] >= param.Nnode)	{
					cerr << "index overflow\n";
					exit(1);
				}
				cerr << param.CalibIndex[i] << '\n';
				cerr << "ok\n";
				*/
				param.isCalibrated[param.CalibIndex[i]] = 1;
			}
		}

		is >> param.SaveEvery;
		is >> param.StopAfter;
		is >> param.HowManySaved;
		is >> param.SwapFreq;

		is >> param.BurnIn;
		is >> param.InitBeta;
		is >> param.FinalBeta;
		is >> param.BetaStep;
		is >> param.BurnInDone;
		is >> param.Beta0;
		is >> param.Nchain;

		for (int l=0; l<param.Nchain; l++)	{
			is >> param.Beta[l];
		}

		param.SwapPermut = new int[param.Nchain];
		for (int i=0; i<param.Nchain; i++)	{
			is >> param.SwapPermut[i];
		}

		is >> param.MSMode;
		is >> param.RASModelSwitch;
		is >> param.SUBModelSwitch;
		is >> param.SeparateModelSwitch;
		is >> param.HeteroModelSwitch;
		is >> param.FlexClockModelSwitch;
		is >> param.BLModelSwitch;

		is >> param.HeteroMode;

		// chain parameters

		is >> param.TopoMoveTypeNumber;
		for (int k=0; k<param.TopoMoveTypeNumber; k++)	{
			is >> param.TopoMoveTypeArray[k];
			is >> param.SuperMoveTypeNumber[k];
			for (int l=0; l<param.Nchain; l++)	{
				is >> param.TopoNIterationArray[k][l];
				is >> param.TopoNArray[k][l];
			}
			for (int i=0; i<param.SuperMoveTypeNumber[k]; i++)	{
				is >> param.MoveTypeNumber[k][i];
				for (int l=0; l<param.Nchain; l++)	{
					is >> param.SuperNIterationArray[k][i][l];
				}
				for (int j=0; j<param.MoveTypeNumber[k][i]; j++)	{
					is >> param.MoveTypeArray[k][i][j];
					for (int l=0; l<param.Nchain; l++)	{
						is >> param.deltaArray[k][i][j][l];
						is >> param.NIterationArray[k][i][j][l];
						is >> param.NArray[k][i][j][l];
					}
				}
			}
		}

		// chain statistics

		is >> param.RateInfCount;
		is >> param.StatInfCount;
		is >> param.LengthInfCount;
		is >> param.LogProbInfCount;
		is >> param.SubOverflowCount;


		for (int k=0; k<param.TopoMoveTypeNumber; k++)	{
			for (int i=0; i<param.SuperMoveTypeNumber[k]; i++)	{
				for (int l=0; l<param.MoveTypeNumber[k][i]; l++)	{
					for (int m=0; m<param.Nchain; m++)	{
						is >> param.NCallArray[k][i][l][m];
						is >> param.SuccessArray[k][i][l][m];
						is >> param.TimeArray[k][i][l][m];
					}
				}
			}
		}
		for (int m=0; m<param.Nchain; m++)	{
			is >> param.TrialSwapArray[m];
			is >> param.AcceptedSwapArray[m];
		}


		is >> param.FixGamma;
		is >> param.FixMeanLength;
		is >> param.FixStatCenter;
		is >> param.FixStat;
		is >> param.FixRR;
		is >> param.FixLength;
		is >> param.FixRate;
		is >> param.FixNmode;
		is >> param.FixTopo;
		is >> param.FixPconst;
		is >> param.TimePrior;
		is >> param.ClockPriorModelSwitch;
		is >> param.AutoCorrelModelSwitch;
		is >> param.ArcThermo;
		is >> param.alpha1;
		is >> param.alpha2;
		is >> param.Zeta;
		is >> param.logBF;
		is >> param.logBFmin;
		is >> param.logBFmax;
		is >> param.logBFrev;
		is >> param.logBFminrev;
		is >> param.logBFmaxrev;
		is >> param.BFcount;
		is >> param.ScalePrior;
		is >> param.MeanScale;
		is >> param.MutMode;
		is >> param.SelMode;

		is >> param.SoftBounds;
		is >> param.Softa;
		is >> param.VarScale;
		is >> param.ImproperLowerBound;
		is >>  param.LowerC;
		is >>  param.LowerP;

		is >> param.AlphaChi;
		is >> param.BetaChi;
		is >> param.AlphaChi2;
		is >> param.BetaChi2;

		is >> param.ZipSub;
		is >> param.LengthGammaPrior;
		
		int cont;

		is >> cont;
		if (cont)	{
			is >> param.ModeStatAlphaPrior;
			// cerr << param.ModeStatAlphaPrior << '\n';
			is >> cont;
			if (cont)	{
				cerr << "error in MCParameters::icstream\n";
				exit(1);
			}
		}

		// current state
		int tmp = param.SaveAll;
		param.SaveAll = 1;
		param.currentState = new PhyloBayes*[param.Nchain];
		param.nextState = new PhyloBayes*[param.Nchain];
		for (int chain = 0; chain < param.Nchain; chain++)	{
			param.currentState[chain] = new PhyloBayes(&param);
			param.nextState[chain] = new PhyloBayes(&param);
			is >> *param.currentState[chain];
		}
		param.SaveAll = tmp;
	}
	else if ((param.Version == 3) && (param.SubVersion == 3))	{

		// data file
		is >> param.DataFileSpec;
		is >> param.ContDataFileSpec;


		// model specifications
		is >> param.NormalApprox;
		is >> param.NormalConcMode;
		is >> param.ClockModel;
		is >> param.WithPartial;
		is >> param.FlexClockModel;
		is >> param.ActivateClock;
		is >> param.SeparateRhoPrior;

		is >> param.NH;
		is >> param.DeleteConstant;
		is >> param.Normalise;
		is >> param.NormaliseCov;
		is >> param.ModeFastCompute;

		is >> param.ZipGTR;
		is >> param.ZipPrior;
		is >> param.ZipGTRDP;
		is >> param.NZipModeMax;
		is >> param.ModePoisson;
		is >> param.NHPrior;
		if (param.NHPrior == 10)	{
			param.NHPrior = 0;
			param.PopEff = 1;
		}

		// specific to version 3.3
		is >> param.PopEff;
		is >> param.NHPop;
		is >> param.FixRoot;

		is >> param.MBL;
		is >> param.NMixBLMax;
		is >> param.IncrementalMixBLDP;

		is >> param.RefFastCompute;
		is >> param.Parallel;
		is >> param.NHNcatMax;
		is >> param.ExternalCov;
		is >> param.Ngen;
		is >> param.SaveStat;
		is >> param.SaveAll;

		is >> param.NmodeMax;
		is >> param.NRateModeMax;
		is >> param.GammaNcat;
		is >> param.ActivateSumOverModes;
		is >> param.SumOverModes;
		is >> param.ActivateSumOverRateModes;
		is >> param.SumOverRateModes;
		is >> param.Qmode;
		is >> param.GeneBLMode;
		is >> param.GeneBLMultiplier;
		is >> param.GeneRateMode;
		is >> param.GeneGammaMode;
		is >> param.GeneStationaryMode;

		is >> param.LengthPrior;
		is >> param.LengthMin;
		is >> param.LengthMax;

		is >> param.RatePrior;
		is >> param.RateMin;
		is >> param.RateMax;
		is >> param.GammaMin;
		is >> param.GammaMax;
		is >> param.XiMin;
		is >> param.XiMax;

		is >> param.ModePrior;
		is >> param.ModeStatPrior;
		is >> param.StatMin;
		is >> param.StatMax;
		is >> param.StatAlphaMin;
		is >> param.StatAlphaMax;

		is >> param.AlphaPrior;
		is >> param.AlphaMin;
		is >> param.AlphaMax;

		is >> param.TooSmall;
		is >> param.TooLarge;
		is >> param.InfProb;
		is >> param.OutNtaxa;
		if (param.OutNtaxa)	{
			param.OutGroup = new string[param.OutNtaxa];
			for (int i=0; i<param.OutNtaxa; i++)	{
				is >> param.OutGroup[i];
			}
		}
		if (param.NormalApprox == Yes)	{
			int tmp;
			is >> tmp;
			if (tmp)	{
				string temp;
				is >> temp;
				param.ReadCov(temp);
			}
			is >> tmp;
			if (tmp)	{
				string temp;
				is >> temp;
				param.ReadSepCov(temp);
			}
			param.RegisterCov();
		}
		else	{
			param.ReadDataFromFile(param.DataFileSpec);

			is >> param.Recoding;
			if (param.Recoding)	{
				is >> param.Nstate;
				is >> param.RecodingNstate;
				param.RecodingTable = new int[param.Nstate];
				for (int i=0; i<param.Nstate; i++)	{
					is >> param.RecodingTable[i];
				}
				param.RecodingAlphabet = new char[param.RecodingNstate];
				for (int i=0; i<param.RecodingNstate; i++)	{
					is >> param.RecodingAlphabet[i];
				}
			}

			param.RegisterWithData();

			is >> param.Ngene;
			delete[] param.Gene;
			delete[] param.GeneSize;
			delete[] param.GeneFirstSite;
			param.GeneSize = new int[param.Ngene];
			param.GeneFirstSite = new int[param.Ngene];
			param.Gene = new int[param.Nsite];

			for (int i=0; i<param.Ngene; i++)	{
				is >> param.GeneSize[i];
				is >> param.GeneFirstSite[i];
			}
			for (int i=0; i<param.Nsite; i++)	{
				is >> param.Gene[i];
			}
			if (param.ContDataFileSpec != "None")	{
				param.ReadContDataFromFile(param.ContDataFileSpec);
			}	
		}

		is >> param.NCalib;
		if (param.NCalib)	{
			param.isCalibrated = new int[param.Nnode];
			for (int j=0; j<param.Nnode; j++)	{
				param.isCalibrated[j] = 0;
			}
			param.CalibTaxon1 = new string[param.NCalib];
			param.CalibTaxon2 = new string[param.NCalib];
			param.CalibIndex = new int[param.NCalib];
			param.CalibUpper = new double[param.NCalib];
			param.CalibLower = new double[param.NCalib];
			for (int i=0; i<param.NCalib; i++)	{
				is >> param.CalibTaxon1[i] >> param.CalibTaxon2[i] >> param.CalibIndex[i] >> param.CalibUpper[i] >> param.CalibLower[i];
				/*
				if (param.CalibIndex[i] <0)	{
					cerr << "negative index\n";
					exit(1);
				}
				if (param.CalibIndex[i] >= param.Nnode)	{
					cerr << "index overflow\n";
					exit(1);
				}
				cerr << param.CalibIndex[i] << '\n';
				cerr << "ok\n";
				*/
				param.isCalibrated[param.CalibIndex[i]] = 1;
			}
		}

		is >> param.SaveEvery;
		is >> param.StopAfter;
		is >> param.HowManySaved;
		is >> param.SwapFreq;

		is >> param.BurnIn;
		is >> param.InitBeta;
		is >> param.FinalBeta;
		is >> param.BetaStep;
		is >> param.BurnInDone;
		is >> param.Beta0;
		is >> param.Nchain;

		for (int l=0; l<param.Nchain; l++)	{
			is >> param.Beta[l];
		}

		param.SwapPermut = new int[param.Nchain];
		for (int i=0; i<param.Nchain; i++)	{
			is >> param.SwapPermut[i];
		}

		is >> param.MSMode;
		is >> param.RASModelSwitch;
		is >> param.SUBModelSwitch;
		is >> param.SeparateModelSwitch;
		is >> param.HeteroModelSwitch;
		is >> param.FlexClockModelSwitch;
		is >> param.BLModelSwitch;

		is >> param.HeteroMode;

		// chain parameters

		is >> param.TopoMoveTypeNumber;
		for (int k=0; k<param.TopoMoveTypeNumber; k++)	{
			is >> param.TopoMoveTypeArray[k];
			is >> param.SuperMoveTypeNumber[k];
			for (int l=0; l<param.Nchain; l++)	{
				is >> param.TopoNIterationArray[k][l];
				is >> param.TopoNArray[k][l];
			}
			for (int i=0; i<param.SuperMoveTypeNumber[k]; i++)	{
				is >> param.MoveTypeNumber[k][i];
				for (int l=0; l<param.Nchain; l++)	{
					is >> param.SuperNIterationArray[k][i][l];
				}
				for (int j=0; j<param.MoveTypeNumber[k][i]; j++)	{
					is >> param.MoveTypeArray[k][i][j];
					for (int l=0; l<param.Nchain; l++)	{
						is >> param.deltaArray[k][i][j][l];
						is >> param.NIterationArray[k][i][j][l];
						is >> param.NArray[k][i][j][l];
					}
				}
			}
		}

		// chain statistics

		is >> param.RateInfCount;
		is >> param.StatInfCount;
		is >> param.LengthInfCount;
		is >> param.LogProbInfCount;
		is >> param.SubOverflowCount;


		for (int k=0; k<param.TopoMoveTypeNumber; k++)	{
			for (int i=0; i<param.SuperMoveTypeNumber[k]; i++)	{
				for (int l=0; l<param.MoveTypeNumber[k][i]; l++)	{
					for (int m=0; m<param.Nchain; m++)	{
						is >> param.NCallArray[k][i][l][m];
						is >> param.SuccessArray[k][i][l][m];
						is >> param.TimeArray[k][i][l][m];
					}
				}
			}
		}
		for (int m=0; m<param.Nchain; m++)	{
			is >> param.TrialSwapArray[m];
			is >> param.AcceptedSwapArray[m];
		}


		is >> param.FixGamma;
		is >> param.FixMeanLength;
		is >> param.FixStatCenter;
		is >> param.FixStat;
		is >> param.FixRR;
		is >> param.FixLength;
		is >> param.FixRate;
		is >> param.FixNmode;
		is >> param.FixTopo;
		is >> param.FixPconst;
		is >> param.TimePrior;
		is >> param.ClockPriorModelSwitch;
		is >> param.AutoCorrelModelSwitch;
		is >> param.ArcThermo;
		is >> param.alpha1;
		is >> param.alpha2;
		is >> param.Zeta;
		is >> param.logBF;
		is >> param.logBFmin;
		is >> param.logBFmax;
		is >> param.logBFrev;
		is >> param.logBFminrev;
		is >> param.logBFmaxrev;
		is >> param.BFcount;
		is >> param.ScalePrior;
		is >> param.MeanScale;
		is >> param.MutMode;
		is >> param.SelMode;

		is >> param.SoftBounds;
		is >> param.Softa;
		is >> param.VarScale;
		is >> param.ImproperLowerBound;
		is >>  param.LowerC;
		is >>  param.LowerP;

		is >> param.AlphaChi;
		is >> param.BetaChi;
		is >> param.AlphaChi2;
		is >> param.BetaChi2;

		is >> param.ZipSub;
		is >> param.LengthGammaPrior;
		
		int cont;

		is >> cont;
		if (cont)	{
			is >> param.ModeStatAlphaPrior;
			// cerr << param.ModeStatAlphaPrior << '\n';
			is >> cont;
			if (cont)	{
				cerr << "error in MCParameters::icstream\n";
				exit(1);
			}
		}

		// current state
		int tmp = param.SaveAll;
		param.SaveAll = 1;
		param.currentState = new PhyloBayes*[param.Nchain];
		param.nextState = new PhyloBayes*[param.Nchain];
		for (int chain = 0; chain < param.Nchain; chain++)	{
			param.currentState[chain] = new PhyloBayes(&param);
			param.nextState[chain] = new PhyloBayes(&param);
			is >> *param.currentState[chain];
		}
		param.SaveAll = tmp;
	}
	else	{
		if (param.Version == 2)	{

			// data file
			is >> param.DataFileSpec;
			is >> param.NormalApprox;
			is >> param.NormalConcMode;
			is >> param.ClockModel;
			is >> param.WithPartial;
			is >> param.FlexClockModel;
			is >> param.ActivateClock;
			is >> param.SeparateRhoPrior;

			// model specifications

			is >> param.NH;
			if ((param.Version <= 2) && (param.SubVersion < 3))	{
				param.NH = 0;
			}
			is >> param.DeleteConstant;
			is >> param.Normalise;
			is >> param.NormaliseCov;
			is >> param.ModeFastCompute;
			param.ModePoisson = param.ModeFastCompute;
			is >> param.RefFastCompute;
			is >> param.Parallel;
			is >> param.NHNcatMax;
			is >> param.ExternalCov;
			is >> param.Ngen;
			is >> param.SaveStat;
			is >> param.SaveAll;

			is >> param.NmodeMax;
			is >> param.NRateModeMax;
			is >> param.GammaNcat;
			is >> param.ActivateSumOverModes;
			is >> param.SumOverModes;
			is >> param.ActivateSumOverRateModes;
			is >> param.SumOverRateModes;
			is >> param.Qmode;
			is >> param.GeneBLMode;
			is >> param.GeneBLMultiplier;
			is >> param.GeneRateMode;
			is >> param.GeneGammaMode;
			is >> param.GeneStationaryMode;

			is >> param.LengthPrior;
			is >> param.LengthMin;
			is >> param.LengthMax;

			is >> param.RatePrior;
			is >> param.RateMin;
			is >> param.RateMax;
			is >> param.GammaMin;
			is >> param.GammaMax;
			is >> param.XiMin;
			is >> param.XiMax;

			is >> param.ModePrior;
			is >> param.ModeStatPrior;
			is >> param.StatMin;
			is >> param.StatMax;
			is >> param.StatAlphaMin;
			is >> param.StatAlphaMax;

			is >> param.AlphaPrior;
			is >> param.AlphaMin;
			is >> param.AlphaMax;

			is >> param.TooSmall;
			is >> param.TooLarge;
			is >> param.InfProb;
			is >> param.OutNtaxa;
			if (param.OutNtaxa)	{
				param.OutGroup = new string[param.OutNtaxa];
				for (int i=0; i<param.OutNtaxa; i++)	{
					is >> param.OutGroup[i];
				}
			}
			if (param.NormalApprox == Yes)	{
				int tmp;
				is >> tmp;
				if (tmp)	{
					string temp;
					is >> temp;
					param.ReadCov(temp);
				}
				is >> tmp;
				if (tmp)	{
					string temp;
					is >> temp;
					param.ReadSepCov(temp);
				}
				param.RegisterCov();
			}
			else	{
				param.ReadDataFromFile(param.DataFileSpec);

				is >> param.Recoding;
				if (param.Recoding)	{
					is >> param.Nstate;
					is >> param.RecodingNstate;
					param.RecodingTable = new int[param.Nstate];
					for (int i=0; i<param.Nstate; i++)	{
						is >> param.RecodingTable[i];
					}
					param.RecodingAlphabet = new char[param.RecodingNstate];
					for (int i=0; i<param.RecodingNstate; i++)	{
						is >> param.RecodingAlphabet[i];
					}
				}

				param.RegisterWithData();

				is >> param.Ngene;
				delete[] param.Gene;
				delete[] param.GeneSize;
				delete[] param.GeneFirstSite;
				param.GeneSize = new int[param.Ngene];
				param.GeneFirstSite = new int[param.Ngene];
				param.Gene = new int[param.Nsite];

				for (int i=0; i<param.Ngene; i++)	{
					is >> param.GeneSize[i];
					is >> param.GeneFirstSite[i];
				}
				for (int i=0; i<param.Nsite; i++)	{
					is >> param.Gene[i];
				}

			}

			is >> param.NCalib;
			if (param.NCalib)	{
				param.isCalibrated = new int[param.Nnode];
				for (int j=0; j<param.Nnode; j++)	{
					param.isCalibrated[j] = 0;
				}
				param.CalibTaxon1 = new string[param.NCalib];
				param.CalibTaxon2 = new string[param.NCalib];
				param.CalibIndex = new int[param.NCalib];
				param.CalibUpper = new double[param.NCalib];
				param.CalibLower = new double[param.NCalib];
				for (int i=0; i<param.NCalib; i++)	{
					is >> param.CalibTaxon1[i] >> param.CalibTaxon2[i] >> param.CalibIndex[i] >> param.CalibUpper[i] >> param.CalibLower[i];
					param.isCalibrated[param.CalibIndex[i]] = 1;
				}
			}

			is >> param.SaveEvery;
			is >> param.StopAfter;
			is >> param.HowManySaved;
			is >> param.SwapFreq;

			is >> param.BurnIn;
			is >> param.InitBeta;
			is >> param.FinalBeta;
			is >> param.BetaStep;
			is >> param.BurnInDone;
			is >> param.Beta0;
			is >> param.Nchain;

			for (int l=0; l<param.Nchain; l++)	{
				is >> param.Beta[l];
			}

			param.SwapPermut = new int[param.Nchain];
			for (int i=0; i<param.Nchain; i++)	{
				is >> param.SwapPermut[i];
			}

			param.currentState = new PhyloBayes*[param.Nchain];
			param.nextState = new PhyloBayes*[param.Nchain];

			is >> param.MSMode;
			is >> param.RASModelSwitch;
			is >> param.SUBModelSwitch;
			is >> param.SeparateModelSwitch;
			is >> param.HeteroModelSwitch;
			is >> param.FlexClockModelSwitch;
			is >> param.BLModelSwitch;

			is >> param.HeteroMode;

			// chain parameters

			is >> param.TopoMoveTypeNumber;
			for (int k=0; k<param.TopoMoveTypeNumber; k++)	{
				is >> param.TopoMoveTypeArray[k];
				is >> param.SuperMoveTypeNumber[k];
				for (int l=0; l<param.Nchain; l++)	{
					is >> param.TopoNIterationArray[k][l];
					is >> param.TopoNArray[k][l];
				}
				for (int i=0; i<param.SuperMoveTypeNumber[k]; i++)	{
					is >> param.MoveTypeNumber[k][i];
					for (int l=0; l<param.Nchain; l++)	{
						is >> param.SuperNIterationArray[k][i][l];
					}
					for (int j=0; j<param.MoveTypeNumber[k][i]; j++)	{
						is >> param.MoveTypeArray[k][i][j];
						for (int l=0; l<param.Nchain; l++)	{
							is >> param.deltaArray[k][i][j][l];
							is >> param.NIterationArray[k][i][j][l];
							is >> param.NArray[k][i][j][l];
						}
					}
				}
			}

			// chain statistics

			is >> param.RateInfCount;
			is >> param.StatInfCount;
			is >> param.LengthInfCount;
			is >> param.LogProbInfCount;
			is >> param.SubOverflowCount;


			for (int k=0; k<param.TopoMoveTypeNumber; k++)	{
				for (int i=0; i<param.SuperMoveTypeNumber[k]; i++)	{
					for (int l=0; l<param.MoveTypeNumber[k][i]; l++)	{
						for (int m=0; m<param.Nchain; m++)	{
							is >> param.NCallArray[k][i][l][m];
							is >> param.SuccessArray[k][i][l][m];
							is >> param.TimeArray[k][i][l][m];
						}
					}
				}
			}
			for (int m=0; m<param.Nchain; m++)	{
				is >> param.TrialSwapArray[m];
				is >> param.AcceptedSwapArray[m];
			}

			int tmp = param.SaveAll;
			param.SaveAll = 1;
			for (int chain = 0; chain < param.Nchain; chain++)	{
				param.currentState[chain] = new PhyloBayes(&param);
				param.nextState[chain] = new PhyloBayes(&param);
				is >> *param.currentState[chain];
			}
			param.SaveAll = tmp;

			is >> param.FixGamma;
			is >> param.FixMeanLength;
			is >> param.FixStatCenter;
			is >> param.FixStat;
			is >> param.FixRR;
			is >> param.FixLength;
			is >> param.FixRate;
			is >> param.FixNmode;
			is >> param.FixTopo;
			is >> param.FixPconst;

			int cont;
			is >> cont;
			if (cont)	{
				is >> param.TimePrior;
				is >> param.ClockPriorModelSwitch;
				is >> param.AutoCorrelModelSwitch;
				is >> param.ArcThermo;
				is >> param.alpha1;
				is >> param.alpha2;
				is >> param.Zeta;
				is >> param.logBF;
				is >> param.logBFmin;
				is >> param.logBFmax;
				is >> param.logBFrev;
				is >> param.logBFminrev;
				is >> param.logBFmaxrev;
				is >> param.BFcount;
				is >> param.ScalePrior;
				is >> param.MeanScale;

				int cont;
				is >> cont;
				if (cont)	{
					cerr << "error in MCParameters::icstream\n";
					exit(1);
				}
			}
		}
		else	{

			// data file
			is >> param.DataFileSpec;
			is >> param.ContDataFileSpec;
			is >> param.NormalApprox;
			is >> param.NormalConcMode;
			is >> param.ClockModel;
			is >> param.WithPartial;
			is >> param.FlexClockModel;
			is >> param.ActivateClock;
			is >> param.SeparateRhoPrior;

			// model specifications

			is >> param.NH;
			if ((param.Version <= 2) && (param.SubVersion < 3))	{
				param.NH = 0;
			}
			is >> param.DeleteConstant;
			is >> param.Normalise;
			is >> param.NormaliseCov;
			is >> param.ModeFastCompute;
			if (param.Version >=3)	{
				is >> param.ZipGTR;
				is >> param.ZipPrior;
				is >> param.ZipGTRDP;
				is >> param.NZipModeMax;
				is >> param.ModePoisson;
				is >> param.NHPrior;
				is >> param.MBL;
				is >> param.NMixBLMax;
				is >> param.IncrementalMixBLDP;
			}
			else	{
				param.ModePoisson = param.ModeFastCompute;
			}
			is >> param.RefFastCompute;
			is >> param.Parallel;
			is >> param.NHNcatMax;
			is >> param.ExternalCov;
			is >> param.Ngen;
			is >> param.SaveStat;
			is >> param.SaveAll;

			is >> param.NmodeMax;
			is >> param.NRateModeMax;
			is >> param.GammaNcat;
			is >> param.ActivateSumOverModes;
			is >> param.SumOverModes;
			is >> param.ActivateSumOverRateModes;
			is >> param.SumOverRateModes;
			is >> param.Qmode;
			is >> param.GeneBLMode;
			is >> param.GeneBLMultiplier;
			is >> param.GeneRateMode;
			is >> param.GeneGammaMode;
			is >> param.GeneStationaryMode;

			is >> param.LengthPrior;
			is >> param.LengthMin;
			is >> param.LengthMax;

			is >> param.RatePrior;
			is >> param.RateMin;
			is >> param.RateMax;
			is >> param.GammaMin;
			is >> param.GammaMax;
			is >> param.XiMin;
			is >> param.XiMax;

			is >> param.ModePrior;
			is >> param.ModeStatPrior;
			is >> param.StatMin;
			is >> param.StatMax;
			is >> param.StatAlphaMin;
			is >> param.StatAlphaMax;

			is >> param.AlphaPrior;
			is >> param.AlphaMin;
			is >> param.AlphaMax;

			is >> param.TooSmall;
			is >> param.TooLarge;
			is >> param.InfProb;
			is >> param.OutNtaxa;
			if (param.OutNtaxa)	{
				param.OutGroup = new string[param.OutNtaxa];
				for (int i=0; i<param.OutNtaxa; i++)	{
					is >> param.OutGroup[i];
				}
			}
			if (param.NormalApprox == Yes)	{
				int tmp;
				is >> tmp;
				if (tmp)	{
					string temp;
					is >> temp;
					param.ReadCov(temp);
				}
				is >> tmp;
				if (tmp)	{
					string temp;
					is >> temp;
					param.ReadSepCov(temp);
				}
				param.RegisterCov();
			}
			else	{
				param.ReadDataFromFile(param.DataFileSpec);

				is >> param.Recoding;
				if (param.Recoding)	{
					is >> param.Nstate;
					is >> param.RecodingNstate;
					param.RecodingTable = new int[param.Nstate];
					for (int i=0; i<param.Nstate; i++)	{
						is >> param.RecodingTable[i];
					}
					param.RecodingAlphabet = new char[param.RecodingNstate];
					for (int i=0; i<param.RecodingNstate; i++)	{
						is >> param.RecodingAlphabet[i];
					}
				}

				param.RegisterWithData();

				is >> param.Ngene;
				delete[] param.Gene;
				delete[] param.GeneSize;
				delete[] param.GeneFirstSite;
				param.GeneSize = new int[param.Ngene];
				param.GeneFirstSite = new int[param.Ngene];
				param.Gene = new int[param.Nsite];

				for (int i=0; i<param.Ngene; i++)	{
					is >> param.GeneSize[i];
					is >> param.GeneFirstSite[i];
				}
				for (int i=0; i<param.Nsite; i++)	{
					is >> param.Gene[i];
				}
				if (param.ContDataFileSpec != "None")	{
					param.ReadContDataFromFile(param.ContDataFileSpec);
				}	
			}

			is >> param.NCalib;
			if (param.NCalib)	{
				param.isCalibrated = new int[param.Nnode];
				for (int j=0; j<param.Nnode; j++)	{
					param.isCalibrated[j] = 0;
				}
				param.CalibTaxon1 = new string[param.NCalib];
				param.CalibTaxon2 = new string[param.NCalib];
				param.CalibIndex = new int[param.NCalib];
				param.CalibUpper = new double[param.NCalib];
				param.CalibLower = new double[param.NCalib];
				for (int i=0; i<param.NCalib; i++)	{
					is >> param.CalibTaxon1[i] >> param.CalibTaxon2[i] >> param.CalibIndex[i] >> param.CalibUpper[i] >> param.CalibLower[i];
					param.isCalibrated[param.CalibIndex[i]] = 1;
				}
			}

			is >> param.SaveEvery;
			is >> param.StopAfter;
			is >> param.HowManySaved;
			is >> param.SwapFreq;

			is >> param.BurnIn;
			is >> param.InitBeta;
			is >> param.FinalBeta;
			is >> param.BetaStep;
			is >> param.BurnInDone;
			is >> param.Beta0;
			is >> param.Nchain;

			for (int l=0; l<param.Nchain; l++)	{
				is >> param.Beta[l];
			}

			param.SwapPermut = new int[param.Nchain];
			for (int i=0; i<param.Nchain; i++)	{
				is >> param.SwapPermut[i];
			}

			param.currentState = new PhyloBayes*[param.Nchain];
			param.nextState = new PhyloBayes*[param.Nchain];

			is >> param.MSMode;
			is >> param.RASModelSwitch;
			is >> param.SUBModelSwitch;
			is >> param.SeparateModelSwitch;
			is >> param.HeteroModelSwitch;
			is >> param.FlexClockModelSwitch;
			is >> param.BLModelSwitch;

			is >> param.HeteroMode;

			// chain parameters

			is >> param.TopoMoveTypeNumber;
			for (int k=0; k<param.TopoMoveTypeNumber; k++)	{
				is >> param.TopoMoveTypeArray[k];
				is >> param.SuperMoveTypeNumber[k];
				for (int l=0; l<param.Nchain; l++)	{
					is >> param.TopoNIterationArray[k][l];
					is >> param.TopoNArray[k][l];
				}
				for (int i=0; i<param.SuperMoveTypeNumber[k]; i++)	{
					is >> param.MoveTypeNumber[k][i];
					for (int l=0; l<param.Nchain; l++)	{
						is >> param.SuperNIterationArray[k][i][l];
					}
					for (int j=0; j<param.MoveTypeNumber[k][i]; j++)	{
						is >> param.MoveTypeArray[k][i][j];
						for (int l=0; l<param.Nchain; l++)	{
							is >> param.deltaArray[k][i][j][l];
							is >> param.NIterationArray[k][i][j][l];
							is >> param.NArray[k][i][j][l];
						}
					}
				}
			}

			// chain statistics

			is >> param.RateInfCount;
			is >> param.StatInfCount;
			is >> param.LengthInfCount;
			is >> param.LogProbInfCount;
			is >> param.SubOverflowCount;


			for (int k=0; k<param.TopoMoveTypeNumber; k++)	{
				for (int i=0; i<param.SuperMoveTypeNumber[k]; i++)	{
					for (int l=0; l<param.MoveTypeNumber[k][i]; l++)	{
						for (int m=0; m<param.Nchain; m++)	{
							is >> param.NCallArray[k][i][l][m];
							is >> param.SuccessArray[k][i][l][m];
							is >> param.TimeArray[k][i][l][m];
						}
					}
				}
			}
			for (int m=0; m<param.Nchain; m++)	{
				is >> param.TrialSwapArray[m];
				is >> param.AcceptedSwapArray[m];
			}

			int tmp = param.SaveAll;
			param.SaveAll = 1;
			for (int chain = 0; chain < param.Nchain; chain++)	{
				param.currentState[chain] = new PhyloBayes(&param);
				param.nextState[chain] = new PhyloBayes(&param);
				is >> *param.currentState[chain];
			}
			param.SaveAll = tmp;

			is >> param.FixGamma;
			is >> param.FixMeanLength;
			is >> param.FixStatCenter;
			is >> param.FixStat;
			is >> param.FixRR;
			is >> param.FixLength;
			is >> param.FixRate;
			is >> param.FixNmode;
			// is >> param.FixNRateMode;
			is >> param.FixTopo;
			is >> param.FixPconst;

			int cont;
			is >> cont;
			if (cont)	{
				is >> param.TimePrior;
				is >> param.ClockPriorModelSwitch;
				is >> param.AutoCorrelModelSwitch;
				is >> param.ArcThermo;
				is >> param.alpha1;
				is >> param.alpha2;
				is >> param.Zeta;
				is >> param.logBF;
				is >> param.logBFmin;
				is >> param.logBFmax;
				is >> param.logBFrev;
				is >> param.logBFminrev;
				is >> param.logBFmaxrev;
				is >> param.BFcount;
				is >> param.ScalePrior;
				is >> param.MeanScale;

				int cont;
				is >> cont;
				if (cont)	{
					double tmp;
					is >> tmp;
					int cont;
					is >> cont;
					if (cont)	{
						cerr << "error in MCParameters::icstream\n";
						exit(1);
					}
				}
			}
		}
	}

	return is;
}

