#include "phylo.h"


// ***************************************************************************
// ***************************************************************************
//		 Constructing 
// ***************************************************************************
// ***************************************************************************


// ---------------------------------------------------------------------------
//		 PhyloBayes( const MCParameters* inParam)
// ---------------------------------------------------------------------------

PhyloBayes::PhyloBayes(MCParameters* inParam) 	{


	Version = 4;
	SubVersion = 1;

	StatePostProb = 0;
	NHNcat = 0;
	Integral = 0;
	
	mParam = inParam;

	RecomputeMatrix = 1;
	GWf = mParam->GWf;

	One = 1;
	Zero = 0;

	MissingFlag = 0;
	Partialpdown = 0;

	NDir = 0;
	DCM = 0;
	Scale = 1000;
	VarScale = 1;


	Beta = 1;
	Nstate = mParam->Nstate;
	Nrr = (Nstate * (Nstate-1)) / 2;
	ContNstate = Nstate + mParam->ContNchar;
	ContNrr = (ContNstate * (ContNstate-1)) / 2;

	Alphabet = mParam->Alphabet;

	InitTree = Random;
	InitBranchLength = Random;
	InitInternalLength = DefaultInitInternalLength;
	InitLeafLength = DefaultInitLeafLength;

	InitRate = Random;
	InitGamma = Random;
	InitPconst = Random;
	InitXi = Random;

	InitRefRR = LG;
	InitRefStat = Uniform;
	InitRefStatCenter = Random;
	InitRefStatAlpha = Random;

	InitModeRR = Random;
	InitModeStat = Random;
	InitModeStatCenter = Empirical;
	InitModeStatAlpha = Fixed;

	InitModeAffiliation = Random;
	InitNmode = Random;
	InitAlpha = Random;
	InitNRateMode = Random;
	InitRateAlpha = Random;

	IsInitTree = No;
	IsInitBranchLength = No;
	IsInitMeanLength = No;
	IsInitNH = No;

	IsInitRate = No;
	IsInitGamma = No;
	IsInitXi = No;
	pconst = 0;
	IsInitPconst = Yes;

	IsInitRefRR = No;
	IsInitRefStat = No;
	IsInitRefStatCenter = No;
	IsInitRefStatAlpha = No;

	IsInitModeWeight = No;
	IsInitModeRR = No;
	IsInitModeStat = No;
	IsInitModeStatCenter = No;
	IsInitModeStatAlpha = No;

	IsInitModeAffiliation = No;
	IsInitNmode = No;
	IsInitAlpha = No;
	IsInitNRateMode = No;
	IsInitRateAlpha = No;

	if (mParam->ExternalCov)	{
		Ncov = 1;
	}
	else	{
		Ncov = mParam->NRateModeMax;
	}
	//
	//--------------
	// allocations
	//--------------
	MixBL = 0;
	MixBLAlpha = 1;
	
	ZipGTRData = 0;

	BKData = 0;
	BKZipData = 0;
	TrueState = 0;
	CreateTree();
	CreateGene();
	if (mParam->ActivateClock)	{
		CreateClock();
	}
	if (mParam->NormalApprox == No)	{
		CreateHetero();
		CreateRate();
		CreateRateMode();
		CreateRef();
		CreateMutSelSub();
		CreateMode();
		CreateNH();
		SubCreated = 0;
		CreateSub();
		CreateBaseFields();
	}

	SetData(-1);
	CreateLogProbs();		
	NHNsub = 0;
	TempZipStationary = 0;
	
}

// ---------------------------------------------------------------------------
//		 CreateStatePostProb()
// ---------------------------------------------------------------------------

void PhyloBayes::CreateStatePostProb()	{

	if (! StatePostProb)	{
		StatePostProb = new double**[mParam->Nnode];
		for (int i=0; i<mParam->Nnode; i++)	{
			StatePostProb[i] = new double*[mParam->Nsite];
			for (int j=0; j<mParam->Nsite; j++)	{
				StatePostProb[i][j] = new double[Nstate];
			}
		}
	}
}

// ---------------------------------------------------------------------------
//		 DeleteStatePostProb()
// ---------------------------------------------------------------------------

void PhyloBayes::DeleteStatePostProb()	{

	if (StatePostProb)	{
		for (int i=0; i<mParam->Nnode; i++)	{
			for (int j=0; j<mParam->Nsite; j++)	{
				delete[] StatePostProb[i][j];
			}
			delete[] StatePostProb[i];
		}
		delete[] StatePostProb;
		StatePostProb = 0;
	}
}

string PhyloBayes::GetLeftNodeName(int j)	{

	if (tree[j].isLeaf())	{
		return mParam->SpeciesNames[j];
	}
	return GetLeftNodeName(tree[j].left->label);
}

string PhyloBayes::GetRightNodeName(int j)	{
	if (tree[j].isLeaf())	{
		return mParam->SpeciesNames[j];
	}
	return GetRightNodeName(tree[j].right->label);

}

// ---------------------------------------------------------------------------
//		 CreateTree()
// ---------------------------------------------------------------------------

void PhyloBayes::CreateTree()	{

	if (verbose)	{
		cerr << "create tree\n";
	}

	// tree structure
	tree = new AmphiNode[mParam->Nnode];
	for (int i=0; i<mParam->Nnode;i++)	{
		tree[i].label = i;
	}
	BL = new double[mParam->Nnode];
	RigidBL = new double[mParam->Nnode];
	for (int j=0; j<mParam->Nnode; j++)	{
		BL[j] = 0;
		RigidBL[j] = 0;
	}

	if (mParam->MBL)	{
		MixBL = new double*[mParam->NMixBLMax];
		for (int i=0; i<mParam->NMixBLMax; i++)	{
			MixBL[i] = new double[mParam->Nnode];
		}
		MixBLSiteNumber = new int[mParam->NMixBLMax];
		MixBLMode = new int[mParam->Nsite];
	}

	MeanLength = 0.1;
	VarLength = 1;
	LengthGamma = 1;
}	


// ---------------------------------------------------------------------------
//		 CreateClock()
// ---------------------------------------------------------------------------

void PhyloBayes::CreateClock()	{
	
	Rho = new double[mParam->Nnode];
	if (mParam->ClockModel == AutoRegressive)	{
		for (int i=0; i<mParam->Nnode;i++)	{
			Rho[i] = 0;
		}
	}
	else	{
		for (int i=0; i<mParam->Nnode;i++)	{
			Rho[i] = 1;
		}
	}
	if (mParam->ContNchar)	{
		ContRho = new double*[mParam->ContNchar+1];
		ContRho[0] = Rho;
		for (int i=0; i<mParam->Nnode;i++)	{
			Rho[i] = 0;
		}
		for (int k=0; k<mParam->ContNchar; k++)	{
			ContRho[k+1] = new double[mParam->Nnode];
			for (int i=0; i<mParam->Nnode;i++)	{
				if ((i<mParam->Ntaxa) && (! mParam->ContMissingData[i][k]))	{
					ContRho[k+1][i] = mParam->ContData[i][k];
				}
				else	{
					ContRho[k+1][i] = 0;
				}
			}
		}
		cerr << "OK\n";
		ContRhoCenter = new double[mParam->ContNchar+1];
		for (int k=0; k<mParam->ContNchar+1; k++)	{
			ContRhoCenter[k] = 0;
		}
		ContVar = new double[mParam->ContNchar+1];
		ContCov = new double*[mParam->ContNchar+1];
		for (int k=0; k<mParam->ContNchar+1; k++)	{
			ContCov[k] = new double[mParam->ContNchar+1];
		}
	}
	if (mParam->Ngene)	{
		GeneRho = new double*[mParam->Ngene];
		for (int k=0; k<mParam->Ngene; k++)	{
			GeneRho[k] = new double[mParam->Nnode];
			for (int i=0; i<mParam->Nnode;i++)	{
				GeneRho[k][i] = 1;
			}
		}
		GeneTheta = new double[mParam->Ngene];
		GeneSigma = new double[mParam->Ngene];
		GeneMu = new double[mParam->Ngene];
		for (int k=0; k<mParam->Ngene; k++)	{
			GeneTheta[k] = 1;
			GeneSigma[k] = 1;
			GeneMu[k] = 1;
		}
	}
	
	MulSigma = 1;
	MulTheta = 1;
	MeanSigma = 1;
	MeanTheta = 1;
	VarSigma = 1;
	VarTheta = 1;
	Sigma = 1;
	Theta = 1;
	Mu = 0.001;
	Chi3 = 0.01;
	Chi2 = 2;
	if (mParam->Chi == -1)	{
		Chi = 1;
		Chi2 = 30;
	}
	else	{
		Chi = mParam->Chi;
		Chi2 = mParam->Chi2;
	}
	if (mParam->Nrho)	{
		Drho = new double[mParam->Nrho];
		UpperRhoLogL = new double*[mParam->Nnode];
		LowerRhoLogL = new double*[mParam->Nnode];
		RhoLogNorm = new double*[mParam->Nnode];
		for (int j=0; j<mParam->Nnode; j++)	{
			UpperRhoLogL[j] = new double[mParam->Nrho];
			LowerRhoLogL[j] = new double[mParam->Nrho];
			RhoLogNorm[j] = new double[mParam->Nrho];
		}
	}

	if (mParam->TimePrior == BD)	{
		mBDLogGG = new double[mParam->Nnode];
		for (int j=0; j<mParam->Nnode; j++)	{
			mBDLogGG[j] = 0;
		}
	}
}


double PhyloBayes::GetEfflogRho()	{

	double Scale0 = 1000;
	return log(mParam->Ntaxa) - Chi * Scale / Scale0;
}

double PhyloBayes::GetEffLambda()	{

	return exp(Chi2 - GetEfflogRho());
}

double PhyloBayes::GetEffMu()	{

	return GetEffLambda() - Chi;
}
 
// ---------------------------------------------------------------------------
//		 CreateRate()
// ---------------------------------------------------------------------------

void PhyloBayes::CreateRate()	{
	
	rate = new double[mParam->Nsite];
	unirate = new double[mParam->Nsite];

	for (int i=0; i<mParam->Nsite; i++)	{
		unirate[i] = 1.0;
		rate[i] = 1.0;
	}
}


// ---------------------------------------------------------------------------
//		 CreateGene()
// ---------------------------------------------------------------------------

void PhyloBayes::CreateGene()	{
	
	if (mParam->Ngene)	{
		GeneBranch = new int*[mParam->Ngene];
		for (int i=0; i<mParam->Ngene; i++)	{
			GeneBranch[i] = new int[mParam->Nnode];
		}
		GeneBL = new double*[mParam->Ngene];
		GeneRigidBL = new double*[mParam->Ngene];
		for (int i=0; i<mParam->Ngene; i++)	{
			GeneBL[i] = new double[mParam->Nnode];
			GeneRigidBL[i] = new double[mParam->Nnode];
			for (int j=0; j<mParam->Nnode; j++)	{
				GeneBL[i][j] = 0;
				GeneRigidBL[i][j] = 0;
			}
		}

		if (mParam->NormalApprox == No)	{
			GeneGamma = new double[mParam->Ngene];
			GeneRate = new double[mParam->Ngene];
			GeneStationary = new double*[mParam->Ngene];
			for (int i=0; i<mParam->Ngene; i++)	{
				GeneStationary[i] = new double[Nstate];
			}
			mGeneMatrixArray = new SubMatrix*[mParam->Ngene];
			GeneRateFactor = new double[mParam->Ngene];
			GeneZipStationary = new double*[mParam->Nsite];
			for (int i=0; i<mParam->Ngene; i++)	{
				mGeneMatrixArray [i] = 0;
				GeneRateFactor[i] = 1;
			}
			for (int i=0; i<mParam->Nsite; i++)	{
				GeneZipStationary[i] = 0;
			}

			if (! mParam->RefFastCompute)	{
				for (int i=0; i<mParam->Ngene; i++)	{
					mGeneMatrixArray[i] = new SubMatrix(mParam,Nstate,RefRR,GeneStationary[i]);
				}
			}
			else	{
				for (int i=0; i<mParam->Nsite; i++)	{
					GeneZipStationary[i] = new double[Nstate];
				}
			}
		}
	}
}



// ---------------------------------------------------------------------------
//		 CreateHetero()
// ---------------------------------------------------------------------------

void PhyloBayes::CreateHetero()	{

	xi0 = 1;
	invprob = 0.5;
}


// ---------------------------------------------------------------------------
//		 CreateTrueSub()
// ---------------------------------------------------------------------------

void PhyloBayes::CreateTrueSub()	{

	BranchSiteTotalTrueSub = new int*[mParam->Nnode];
	for (int i=0; i<mParam->Nnode; i++)	{
		BranchSiteTotalTrueSub[i] = new int[mParam->Nsite];
	}
	TotalTrueSub = new int[mParam->Nsite];
	TrueSub = new int*[mParam->Nsite];
	TimeSpent = new double*[mParam->Nsite];
	for (int i=0; i<mParam->Nsite; i++)	{
		TrueSub[i] = new int[mParam->Nstate];
		TimeSpent[i] = new double[mParam->Nstate];
	}
	BranchTotalTrueSub = new int[mParam->Nnode];
	if (mParam->HeteroMode == Covarion)	{
		BranchSiteCovSwitch = new int*[mParam->Nnode];
		BranchSiteTimeOn = new double*[mParam->Nnode];
		BranchSiteTimeOff = new double*[mParam->Nnode];
		for (int j=0; j<mParam->Nnode; j++)	{
			BranchSiteCovSwitch[j] = new int[mParam->Nsite];
			BranchSiteTimeOn[j] = new double[mParam->Nsite];
			BranchSiteTimeOff[j] = new double[mParam->Nsite];
		}
	}
}

// ---------------------------------------------------------------------------
//		 WriteSiteCount(string filename)
// ---------------------------------------------------------------------------

void PhyloBayes::WriteSiteCount(string filename)	{

	cerr << filename << '\n';
	ofstream os(filename.c_str());
	os << mParam->Nsite << '\t' << Nstate << '\n';
	for (int i=0; i<mParam->Nsite; i++)	{
		os << 1 << '\t';
		for (int j=0; j<Nstate; j++)	{
			os << '\t' << Nsub[i][j];
		}
		os << '\n';
	}
	os.close();
}


// ---------------------------------------------------------------------------
//		 CreateSub()
// ---------------------------------------------------------------------------

void PhyloBayes::CreateSub()	{

	if (verbose)	{
		cerr << "create sub\n";
	}

	BranchTotalTrueSub = 0;
	TotalTrueSub = 0;
	TrueSub = 0;
	TimeSpent = 0;
	BranchSiteTotalTrueSub = 0;

	State = new int*[mParam->Nnode];
	BranchSiteTotalSub = new int*[mParam->Nnode];
	for (int i=0; i<mParam->Nnode; i++)	{
		State[i] = new int[mParam->Nsite];
		BranchSiteTotalSub[i] = new int[mParam->Nsite];
	}
	if (mParam->ModeFastCompute)	{
		TrueState = new int*[mParam->Nnode];
		for (int i=0; i<mParam->Nnode; i++)	{
			TrueState[i] = new int[mParam->Nsite];
		}
	}
	if (mParam->HeteroMode == Covarion)	{
		// Those statistics are used for data-augmentation based moves
		// but these are disactivated anyway
		NSiteOnSub = new int[mParam->Nsite];
		NSiteOffSub = new int[mParam->Nsite];
		SiteXiEffLength = new double[mParam->Nsite];
	}

	Nsub = new int*[mParam->Nsite];
	for (int i=0; i<mParam->Nsite; i++)	{
		Nsub[i] = new int[Nstate];
	}

	TotalSub = new int[mParam->Nsite];
	BranchTotalSub = new int[mParam->Nnode];
	BranchEffLength = new double[mParam->Nnode];
	if (mParam->MBL)	{
		MixBLBranchEffLength = new double*[mParam->NMixBLMax];
		MixBLBranchTotalSub = new int*[mParam->NMixBLMax];
		for (int i=0; i<mParam->NMixBLMax; i++)	{
			MixBLBranchEffLength[i] = new double[mParam->Nnode];
			MixBLBranchTotalSub[i] = new int[mParam->Nnode];
		}
	}

	RefNsub = new int[Nstate];

	StatBeta = new double[Nstate];
	RelRateBeta = new double[Nrr];
	RelRateNsub = new int[Nrr];

	SiteStatBeta = new double*[mParam->Nsite];
	SiteRelRateBeta = new double*[mParam->Nsite];
	SiteRelRateNsub = new int*[mParam->Nsite];
	for (int i=0; i<mParam->Nsite; i++)	{
		SiteStatBeta[i] = new double[Nstate];
		SiteRelRateBeta[i] = new double[Nrr];
		SiteRelRateNsub[i] = new int[Nrr];
	}

	SiteRRFactor = 0;
	if (mParam->Qmode)	{
		SiteRRFactor = new double**[mParam->Nsite];
		for (int i=0; i<mParam->Nsite; i++)	{
			SiteRRFactor[i] = new double*[Nstate];
			for (int j=0; j<Nstate; j++)	{
				SiteRRFactor[i][j] = new double[Nstate];
			}
		}
	}

	CreateGeneSub();
		
	SubCreated = 1;

	NDoubleSub = 0;
	NTripleSub = 0;
	countmultiplesub = 0;
}

void PhyloBayes::CreateGeneSub()	{

	GeneBranchTotalSub = new int*[mParam->Ngene];	
	GeneBranchEffLength = new double*[mParam->Ngene];
	GeneTotalSub = new int[mParam->Ngene];
	GeneNsub = new int*[mParam->Ngene];
	GeneRateEffLength = new double[mParam->Ngene];
	for (int gene=0; gene<mParam->Ngene; gene++)	{
		GeneNsub[gene] = new int[Nstate];
		GeneBranchTotalSub[gene] = new int[mParam->Nnode];
		GeneBranchEffLength[gene] = new double[mParam->Nnode];
	}

	GeneStatBeta = new double*[mParam->Ngene];
	for (int gene=0; gene<mParam->Ngene; gene++)	{
		GeneStatBeta[gene] = new double[Nstate];
	}
}

void PhyloBayes::CreateMultipleSub()	{

	NDoubleSub = new int**[mParam->Nsite];
	for (int i=0; i<mParam->Nsite; i++)	{
		NDoubleSub[i] = new int*[Nstate];
		for (int k=0; k<Nstate; k++)	{
			NDoubleSub[i][k] = new int[Nstate];
		}
	}

	NTripleSub = new int***[mParam->Nsite];
	for (int i=0; i<mParam->Nsite; i++)	{
		NTripleSub[i] = new int**[Nstate];
		for (int k=0; k<Nstate; k++)	{
			NTripleSub[i][k] = new int*[Nstate];
			for (int l=0; l<Nstate; l++)	{
				NTripleSub[i][k][l] = new int[Nstate];
			}
		}
	}
}
	
void PhyloBayes::DeleteMultipleSub()	{

	for (int i=0; i<mParam->Nsite; i++)	{
		for (int k=0; k<Nstate; k++)	{
			delete[] NDoubleSub[i][k];
		}
		delete[] NDoubleSub[i];
	}
	delete[] NDoubleSub;

	for (int i=0; i<mParam->Nsite; i++)	{
		for (int k=0; k<Nstate; k++)	{
			for (int l=0; l<Nstate; l++)	{
				delete[] NTripleSub[i][k][l];
			}
			delete[] NTripleSub[i][k];
		}
		delete[] NTripleSub[i];
	}
	delete[] NTripleSub;
}
	

void PhyloBayes::ResetMultipleSub()	{

	for (int i=0; i<mParam->Nsite; i++)	{
		for (int k=0; k<Nstate; k++)	{
			for (int l=0; l<Nstate; l++)	{
				NDoubleSub[i][k][l] = 0;
			}
		}
	}

	for (int i=0; i<mParam->Nsite; i++)	{
		for (int k=0; k<Nstate; k++)	{
			for (int l=0; l<Nstate; l++)	{
				for (int m=0; m<Nstate; m++)	{
					NTripleSub[i][k][l][m] = 0;
				}
			}
		}
	}
}
	
// ---------------------------------------------------------------------------
//		 CreateBaseFields()
// ---------------------------------------------------------------------------

void PhyloBayes::CreateBaseFields()	{

	if (verbose)	{
		cerr << "create base fields\n";
	}
	BaseRate = new double*[mParam->Nsite];
	BaseGeneRate = new double*[mParam->Nsite];
	BaseRateFactor = new double*[mParam->Nsite];
	BaseRefRateFactor = new double*[mParam->Nsite];
	BaseStationary = new double*[mParam->Nsite];
	BaseRefStationary = new double*[mParam->Nsite];
	BaseBL = new double*[mParam->Nsite];
	BaseBLMul = new double*[mParam->Nsite];
	
	BaseRefMatrix = new SubMatrix*[mParam->Nsite];
	BaseMatrix = new SubMatrix*[mParam->Nsite];
	
	UniBLMul = new double[mParam->Nnode];
	EffLength = new double*[mParam->Nnode];
	for (int j=0; j<mParam->Nnode; j++)	{
		EffLength[j] = new double[mParam->Nsite];
		UniBLMul[j] = 1;
	}
	Pconst0 = 0;
	UniGeneRate = 1;
	// ??? 
	BasePoisson = No;

	Basetbl = new double*[mParam->Nnode];
	for (int j=0; j<mParam->Nnode; j++)	{
		Basetbl[j] = 0;
	}
	Localtbl = new double**[mParam->Nnode];
	for (int j=0; j<mParam->Nnode; j++)	{
		Localtbl[j] = 0;
	}
	DCMCenter = new int[mParam->Nnode];
}


// ---------------------------------------------------------------------------
//		 CreateMutSelSub()
// ---------------------------------------------------------------------------

void PhyloBayes::CreateMutSelSub()	{

	if ((mParam->ZipGTR == 4)||(mParam->MutMode)){
		GWNsub = new double*[mParam->Nsite];
		for (int i=0; i<mParam->Nsite; i++)	{
			GWNsub[i] = new double[Nstate];
		}
		GWRootNsub = new double*[mParam->Nsite];
		for (int i=0; i<mParam->Nsite; i++)	{
			GWRootNsub[i] = new double[Nstate];
		}
		GWStatBeta = new double**[mParam->Nsite];
		for (int i=0; i<mParam->Nsite; i++)	{
			GWStatBeta[i] = new double*[mParam->Nstate];
			for (int j=0; j<Nstate; j++)	{
				GWStatBeta[i][j] = new double[mParam->Nstate];
			}
		}

		NucStatBeta = new double**[mParam->Nsite];
		for (int i=0; i<mParam->Nsite; i++)	{
			NucStatBeta[i] = new double*[mParam->Nstate];
			for (int j=0; j<Nstate; j++)	{
				NucStatBeta[i][j] = new double[mParam->Nstate];
			}
		}

		GWDoubleNsub = new int**[mParam->Nsite];
		for (int i=0; i<mParam->Nsite; i++)	{
			GWDoubleNsub[i] = new int*[Nstate];
			for (int j=0; j<Nstate; j++)	{
				GWDoubleNsub[i][j] = new int[Nstate];
			}
		}

		ModeGWDoubleNsub = new int**[mParam->NmodeMax];
		for (int i=0; i<mParam->NmodeMax; i++)	{
			ModeGWDoubleNsub[i] = new int*[Nstate];
			for (int j=0; j<Nstate; j++)	{
				ModeGWDoubleNsub[i][j] = new int[Nstate];
			}
		}

		ModeGWNsub = new double*[mParam->NmodeMax];
		for (int i=0; i<mParam->NmodeMax; i++)	{
			ModeGWNsub[i] = new double[Nstate];
		}
		ModeGWRootNsub = new double*[mParam->NmodeMax];
		for (int i=0; i<mParam->NmodeMax; i++)	{
			ModeGWRootNsub[i] = new double[Nstate];
		}
		ModeGWStatBeta = new double**[mParam->NmodeMax];
		for (int i=0; i<mParam->NmodeMax; i++)	{
			ModeGWStatBeta[i] = new double*[mParam->Nstate];
			for (int j=0; j<Nstate; j++)	{
				ModeGWStatBeta[i][j] = new double[mParam->Nstate];
			}
		}
		ModeNucStatBeta = new double**[mParam->NmodeMax];
		for (int i=0; i<mParam->NmodeMax; i++)	{
			ModeNucStatBeta[i] = new double*[mParam->Nstate];
			for (int j=0; j<Nstate; j++)	{
				ModeNucStatBeta[i][j] = new double[mParam->Nstate];
			}
		}
	}
}

// ---------------------------------------------------------------------------
//		 CreateMode()
// ---------------------------------------------------------------------------

void PhyloBayes::CreateMode()	{

	SiteNumber = new int[mParam->NmodeMax];
	for (int i=0; i<mParam->NmodeMax ; i++)	{
		SiteNumber[i] = 0;
	}
	Mode = new int[mParam->Nsite];

	Stationary = new double*[mParam->NmodeMax];
	for (int i=0; i<mParam->NmodeMax; i++)	{
		Stationary[i] = new double[Nstate];
	}
	RR = 0;
	RRAlpha = 1;

	Kappa = 1;
	RefACGT = 0;
	RefRRACGT = 0;
	if (mParam->MutMode == 2)	{
		RefACGT = new double[Nnuc];
		for (int i=0; i<Nnuc; i++)	{
			RefACGT[i] = 1.0 / Nnuc;
		}
	}
	if (mParam->MutMode == 3)	{
		RefACGT = new double[Nnuc];
		for (int i=0; i<Nnuc; i++)	{
			RefACGT[i] = 1.0 / Nnuc;
		}
		RefRRACGT = new double[Nnuc * (Nnuc-1) / 2];
		for (int i=0; i<Nnuc * (Nnuc-1) / 2; i++)	{
			RefRRACGT[i] = 1.0;
		}
	}

	ModeWeight = new double[mParam->NmodeMax];
	for (int i=0; i<mParam->NmodeMax; i++)	{
		ModeWeight[i] = 1;
	}
	ModeStatCenter = new double[Nstate];
	ModeRR = new double[Nrr];

	mMatrixArray = new SubMatrix*[mParam->NmodeMax];
	ModeRateFactor = new double[mParam->NmodeMax];
	for (int i=0; i<mParam->NmodeMax; i++)	{
		mMatrixArray[i] = 0;
		ModeRateFactor[i] = 1;
	}
	mZipMatrixArray = new SubMatrix*[mParam->Nsite];
	for (int i=0; i<mParam->Nsite; i++)	{
		mZipMatrixArray[i] = 0;
	}

	// Covarion : creating the arrays of matrices

	mCovMatrixArray = new SubMatrix**[mParam->NmodeMax];
	for (int i=0; i<mParam->NmodeMax; i++)	{
		mCovMatrixArray[i] = new SubMatrix*[Ncov];
		for (int k=0; k<Ncov; k++)	{
			mCovMatrixArray[i][k] = 0;
		}
	}
	mCovZipMatrixArray= new SubMatrix**[mParam->Nsite];
	for (int i=0; i<mParam->Nsite; i++)	{
		mCovZipMatrixArray[i] = new SubMatrix*[Ncov];
		for (int k=0; k<Ncov; k++)	{
			mCovZipMatrixArray[i][k] = 0;
		}
	}	

	if (! mParam->ModeFastCompute)	{
		ZipStationary = 0;
		if (mParam->ModePoisson)	{
			cerr << "poisson computation: only under fast computation\n";
			exit(1);
		}
		if (mParam->Qmode == Yes)	{
			RR = new double*[mParam->NmodeMax];
			for (int i=0; i<mParam->NmodeMax; i++)	{
				RR[i] = new double[Nrr];
			}
		}
		if (mParam->HeteroMode == Covarion)	{
			for (int i=0; i<mParam->NmodeMax; i++)	{
				double* rr = 0;
				if (mParam->Qmode)	{
					rr = RR[i];
				}
				else	{
					rr = ModeRR;
				}


				if (mParam->ExternalCov)	{
					mCovMatrixArray[i][0] = new SubMatrix(mParam,Nstate,&One,&invprob,&xi0,rr,Stationary[i]);
				}
				else	{
					for (int k=0; k<Ncov; k++)	{
						mCovMatrixArray[i][k] = new SubMatrix(mParam,Nstate,&ModeRate[k],&invprob,&xi0,rr,Stationary[i]);
					}
				}
			}
		}
		else	{
			for (int i=0; i<mParam->NmodeMax; i++)	{
				double* rr = 0;
				if (mParam->Qmode)	{
					rr = RR[i];
				}
				else	{
					rr = ModeRR;
				}
				int a = 1;
				if (mParam->MutMode == 1)	{
					mMatrixArray[i] = new SubMatrix(mParam,Nstate,rr,RefStationary,Stationary[i],&GWf,a);
				}
				else if (mParam->MutMode == 2)	{
					mMatrixArray[i] = new SubMatrix(mParam,Nstate,&Kappa,RefACGT,Stationary[i],&GWf,a);
				}
				else if (mParam->MutMode == 3)	{
					mMatrixArray[i] = new SubMatrix(mParam,Nstate,RefRRACGT,RefACGT,Stationary[i],&GWf,a);
				}
				else	{
					mMatrixArray[i] = new SubMatrix(mParam,Nstate,rr,Stationary[i]);
				}
			}
		}
	}
	else	{
		if (mParam->ModePoisson)	{
			ZipStationary = new double*[mParam->Nsite];
			for (int i=0; i<mParam->Nsite; i++)	{
				ZipStationary[i] = new double[Nstate];
			}
			if (mParam->HeteroMode == Covarion)	{
				for (int i=0; i<mParam->Nsite; i++)	{
					if (mParam->ExternalCov)	{
						mCovZipMatrixArray[i][0] = new SubMatrix(mParam,mParam->ZipSize[i],&One,&invprob,&xi0,ZipStationary[i]);
					}
					else	{
						for (int k=0; k<Ncov; k++)	{
							mCovZipMatrixArray[i][k] = new SubMatrix(mParam,mParam->ZipSize[i],&ModeRate[k],&invprob,&xi0,ZipStationary[i]);
						}
					}
				}
			}
		}
		else	{
			if (mParam->Qmode == Yes)	{
				RR = new double*[mParam->NmodeMax];
				for (int i=0; i<mParam->NmodeMax; i++)	{
					RR[i] = new double[Nrr];
				}
			}
			if (mParam->ZipGTR == 4)	{
				ZipStationary = new double*[mParam->Nsite];
				for (int i=0; i<mParam->Nsite; i++)	{
					ZipStationary[i] = new double[mParam->ZipSize[i]];
				}
				SiteZipRR = new double*[mParam->Nsite];
				for (int i=0; i<mParam->Nsite; i++)	{
					int nrr = mParam->ZipSize[i] * (mParam->ZipSize[i] - 1) / 2;
					SiteZipRR[i] = new double[nrr];
				}
				for (int i=0; i<mParam->Nsite; i++)	{
					mZipMatrixArray[i] = new SubMatrix(mParam,mParam->ZipSize[i],SiteZipRR[i],ZipStationary[i]);
				}
			}
			else	{
				ModeZipProb = new double[Nstate];
				ZipGTRData = new int*[mParam->Ntaxa];
				for (int i=0; i<mParam->Ntaxa; i++)	{
					ZipGTRData[i] = new int[mParam->Nsite];
				}
				ZipOcc = new int*[mParam->NZipModeMax];
				for (int i=0; i<mParam->NZipModeMax; i++)	{
					ZipOcc[i] = new int[Nstate];
				}
				GlobalZipOcc = new int[Nstate];
				ZipA0 = new double[Nstate];
				ZipA1 = new double[Nstate];
				ZipMode = new int[mParam->Nsite];
				ZipSiteNumber = new int[mParam->NZipModeMax];

				for (int i=0; i<mParam->NmodeMax; i++)	{
					double* rr = 0;
					if (mParam->Qmode)	{
						rr = RR[i];
					}
					else	{
						rr = ModeRR;
					}
					mMatrixArray[i] = new SubMatrix(mParam,Nstate,rr,Stationary[i]);
				}

				if (mParam->ZipGTR >= 2)	{
					SiteZipSize = new int[mParam->Nsite];
					SiteAASupport = new int*[mParam->Nsite];
					SiteAAZip = new int*[mParam->Nsite];
					SiteZipAA = new int*[mParam->Nsite];	
					ZipStationary = new double*[mParam->Nsite];
					SiteZipRR = new double*[mParam->Nsite];
					for (int i=0; i<mParam->Nsite; i++)	{
						SiteAASupport[i] = new int[Nstate];
						SiteAAZip[i] = new int[Nstate];
						SiteZipAA[i] = new int[Nstate];
						ZipStationary[i] = new double[Nstate];
						SiteZipRR[i] = new double[Nrr];
					}
					ZipMode = new int[mParam->Nsite];

					if (mParam->ZipGTRDP)	{
						ModeAASupport = new int*[mParam->NZipModeMax];
						for (int i=0; i<mParam->NZipModeMax; i++)	{
							ModeAASupport[i] = new int[Nstate];
						}
						ModeZipSize = new int[mParam->NZipModeMax];
					}
				}
				else	{
					ModeZipSize = new int[mParam->NmodeMax];
					ModeAAZip = new int*[mParam->NmodeMax];
					ModeZipAA = new int*[mParam->NmodeMax];
					for (int i=0; i<mParam->NmodeMax; i++)	{
						ModeAAZip[i] = new int[Nstate];
						ModeZipAA[i] = new int[Nstate];
					}
					ModeAASupport = new int*[mParam->NmodeMax];
					for (int i=0; i<mParam->NmodeMax; i++)	{
						ModeAASupport[i] = new int[Nstate];
					}
					ModeZipStationary = new double*[mParam->NmodeMax];
					ModeZipRR = new double*[mParam->NmodeMax];
					for (int i=0; i<mParam->NmodeMax; i++)	{
						ModeZipStationary[i] = new double[Nstate];
						ModeZipRR[i] = new double[Nrr];
					}
					ZipStationary = 0;
				}
			}
		}
	}
	ModeGrandTotal = new int[mParam->NmodeMax];
	ModeTotal = new int*[mParam->NmodeMax];
	for (int i=0; i<mParam->NmodeMax; i++)	{
		ModeTotal[i] = new int[Nstate];
	}
	ModeStatBeta = new double*[mParam->NmodeMax];
	for (int i=0; i<mParam->NmodeMax; i++)	{
		ModeStatBeta[i] = new double[Nstate];
	}
	ModeRelRateBeta = 0;
	ModeRelRateNsub = 0;
	if (mParam->Qmode)	{
		ModeRelRateBeta = new double*[mParam->NmodeMax];
		ModeRelRateNsub = new int*[mParam->NmodeMax];
		for (int i=0; i<mParam->NmodeMax; i++)	{
			ModeRelRateBeta[i] = new double[Nrr];
			ModeRelRateNsub[i] = new int[Nrr];
		}
	}
}

void PhyloBayes::RegisterTempZipStationary()	{

	if (! TempZipStationary)	{
		TempZipStationary = new double**[mParam->Nsite];
		for (int i=0; i<mParam->Nsite; i++)	{
			TempZipStationary[i] = new double*[Nmode];
			for (int j=0; j<Nmode; j++)	{
				TempZipStationary[i][j] = new double[mParam->ZipSize[i]];
			}
		}	
	}
	for (int i=0; i<mParam->Nsite; i++)	{
		for (int mode=0; mode<Nmode; mode++)	{
			for (int j=0; j<mParam->OrbitSize[i]; j++)	{
				TempZipStationary[i][mode][j] =  Stationary[mode][mParam->Indices[i][j]];
			}
			if (mParam->ZipSize[i]>mParam->OrbitSize[i])	{
				double total1 = 0;
				for (int j=0; j<Nstate; j++)	{
					if (! mParam->Orbit[i][j])	{
						double temp = Stationary[mode][j];
						total1 += temp;
					}
				}
				TempZipStationary[i][mode][mParam->OrbitSize[i]] = total1;
			}
		}
	}
}

// ---------------------------------------------------------------------------
//		UpdateZipOccupationNumber()
// ---------------------------------------------------------------------------

void PhyloBayes::UpdateZipOccupationNumber()	{

	if (mParam->TimePrior == BD)	{
		RefreshBDLogGG();
	}

	if (mParam->ZipGTR == 4)	{
		cerr << "error in update zip occ number\n";
		exit(1);
	}
	if (mParam->ZipGTR >= 2)	{
		for (int mode=0; mode<NZipMode; mode++)	{
			UpdateZipOccupationNumber(mode);
		}
		if (mParam->ZipGTRDP)	{
			if (mParam->ZipPrior == 1)	{
				for (int k=0; k<Nstate; k++)	{
					GlobalZipOcc[k] = 0;
				}
				for (int i=0; i<NZipMode; i++)	{
					for (int k=0; k<Nstate; k++)	{
						GlobalZipOcc[k] += ModeAASupport[i][k];
					}
				}
			}
			for (int zipmode=0; zipmode<NZipMode; zipmode++)	{
				ModeZipSize[zipmode] = 0;
				for (int k=0; k<Nstate; k++)	{
					ModeZipSize[zipmode] += ModeAASupport[zipmode][k];
				}
			}
		}
	}
	else if ((mParam->ZipGTR == 1) && (mParam->ZipPrior == 1))	{
		for (int k=0; k<Nstate; k++)	{
			GlobalZipOcc[k] = 0;
		}
		for (int i=0; i<Nmode; i++)	{
			for (int k=0; k<Nstate; k++)	{
				GlobalZipOcc[k] += ModeAASupport[i][k];
			}
		}
	}
}

void PhyloBayes::UpdateZipOccupationNumber(int mode)	{
	for (int k=0; k<Nstate; k++)	{
		ZipOcc[mode][k] = 0;
	}
	ZipSiteNumber[mode] = 0;
	for (int i=0; i<mParam->Nsite; i++)	{
		if (ZipMode[i] == mode)	{
			ZipSiteNumber[mode]++;
			for (int k=0; k<Nstate; k++)	{
				ZipOcc[mode][k] += SiteAASupport[i][k];
			}
		}
	}
}
	
// ---------------------------------------------------------------------------
//		 RecreateSiteMatrices()
// ---------------------------------------------------------------------------


void PhyloBayes::RecreateSiteMatrices()	{
	for (int i=0; i<mParam->Nsite; i++)	{
		RecreateSiteMatrix(i);
	}
}

void PhyloBayes::RecreateSiteMatrix(int site)	{
	if (mParam->NH)	{
		for (int n=0; n<NHNcat; n++)	{
			delete mNHMatrixArray[site][n];
			mNHMatrixArray[site][n] = new SubMatrix(mParam,SiteZipSize[site],SiteZipRR[site],NHZipStationary[site][n]);
		}
	}
	else	{
		delete mZipMatrixArray[site];
		mZipMatrixArray[site] = new SubMatrix(mParam,SiteZipSize[site],SiteZipRR[site],ZipStationary[site]);
	}
}

void PhyloBayes::SiteZip()	{
	for (int i=0; i<mParam->Nsite; i++)	{
		SiteZip(i);
	}
}

void PhyloBayes::NHSiteZip(int site, int n)	{

	double total = 0;
	int mode = Mode[site];
	for (int i=0; i<SiteZipSize[site]; i++)	{
		NHZipStationary[site][n][i] = NHStationary[mode][n][SiteZipAA[site][i]];
		total += NHZipStationary[site][n][i];
	}
	for (int i=0; i<SiteZipSize[site]; i++)	{
		NHZipStationary[site][n][i] /= total;
	}
}


void PhyloBayes::SiteZip(int site)	{
	// given information given by ModeAASupport, make the zipping

	// develop recoding
	if (mParam->NH)	{
		int zipsize = 0;
		int mode = Mode[site];
		for (int i=0; i<Nstate; i++)	{
			 if (SiteAASupport[site][i])	{
				for (int n=0; n<NHNcat; n++)	{
					NHZipStationary[site][n][zipsize] = NHStationary[mode][n][i];
				}
				SiteZipAA[site][zipsize] = i;
				SiteAAZip[site][i] = zipsize;
				zipsize++;
			 }
		}
		SiteZipSize[site] = zipsize;

		// make zip stat
		for (int n=0; n<NHNcat; n++)	{
			double total = 0;
			for (int i=0; i<zipsize; i++)	{
				total += NHZipStationary[site][n][i];
			}
			for (int i=0; i<zipsize; i++)	{
				NHZipStationary[site][n][i] /= total;
			}
		}
	}
	else	{
		int zipsize = 0;
		double total = 0;
		int mode = Mode[site];
		for (int i=0; i<Nstate; i++)	{
			 if (SiteAASupport[site][i])	{
				ZipStationary[site][zipsize] = Stationary[mode][i];
				SiteZipAA[site][zipsize] = i;
				SiteAAZip[site][i] = zipsize;
				total += Stationary[mode][i];
				zipsize++;
			 }
		}
		SiteZipSize[site] = zipsize;

		// make zip stat
		for (int i=0; i<zipsize; i++)	{
			ZipStationary[site][i] /= total;
		}
	}
		
	// make zip rel rates
	int k = 0;
	double* rr = 0;
	if (mParam->Qmode)	{
		rr = RR[Mode[site]];
	}
	else	{
		rr = ModeRR;
	}
	
	for (int i=0; i<SiteZipSize[site]; i++)	{
		for (int j=i+1; j<SiteZipSize[site]; j++)	{
			int aa = SiteZipAA[site][i];
			int bb = SiteZipAA[site][j];
			SiteZipRR[site][k] = rr[(2 * Nstate - aa - 1) * aa / 2 + bb - aa - 1];
			k++;
		}
	}
}	

// ---------------------------------------------------------------------------
//		 RecreateModeMatrices()
// ---------------------------------------------------------------------------


void PhyloBayes::RecreateModeMatrices()	{
	for (int i=0; i<Nmode; i++)	{
		RecreateModeMatrix(i);
	}
	UpdateBaseFields();
}

void PhyloBayes::RecreateModeMatrix(int mode)	{
	if (mParam->NH)	{
		for (int n=0; n<NHNcat; n++)	{
			delete mNHMatrixArray[mode][n];
			mNHMatrixArray[mode][n] = new SubMatrix(mParam,ModeZipSize[mode],ModeZipRR[mode],NHModeZipStationary[mode][n]);
		}
	}
	else	{
		delete mMatrixArray[mode];
		mMatrixArray[mode] = new SubMatrix(mParam,ModeZipSize[mode],ModeZipRR[mode],ModeZipStationary[mode]);
	}
}

void PhyloBayes::ModeZip()	{
	for (int i=0; i<Nmode; i++)	{
		ModeZip(i);
	}
}

void PhyloBayes::NHModeZip(int mode, int n)	{

	double total = 0;
	for (int i=0; i<ModeZipSize[mode]; i++)	{
		NHModeZipStationary[mode][n][i] = NHStationary[mode][n][ModeZipAA[mode][i]];
		total += NHModeZipStationary[mode][n][i];
	}
	for (int i=0; i<ModeZipSize[mode]; i++)	{
		NHModeZipStationary[mode][n][i] /= total;
	}
}

void PhyloBayes::ModeZip(int mode)	{
	// given information given by ModeAASupport, make the zipping

	// develop recoding
	if (mParam->NH)	{
		int zipsize = 0;
		for (int i=0; i<Nstate; i++)	{
			if (ModeAASupport[mode][i])	{
				for (int n=0; n<NHNcat; n++)	{
					NHModeZipStationary[mode][n][zipsize] = NHStationary[mode][n][i];
				}
				ModeZipAA[mode][zipsize] = i;
				ModeAAZip[mode][i] = zipsize;
				zipsize++;
			}
			else	{
				Stationary[mode][i] = 0;
				for (int n=0; n<NHNcat; n++)	{
					NHStationary[mode][n][i] = 0;
				}
			}
		}
		ModeZipSize[mode] = zipsize;

		for (int n=0; n<NHNcat; n++)	{
			double total  = 0;
			for (int i=0; i<zipsize; i++)	{
				total += NHModeZipStationary[mode][n][i];
			}
			for (int i=0; i<zipsize; i++)	{
				NHModeZipStationary[mode][n][i] /= total;
			}
			total = 0;
			for (int i=0; i<Nstate; i++)	{
				total += NHStationary[mode][n][i];
			}
			for (int i=0; i<Nstate; i++)	{
				NHStationary[mode][n][i] /= total;
			}
		}
		double total = 0;
		for (int k=0; k<Nstate; k++)	{
			total += Stationary[mode][k];
		}
		for (int k=0; k<Nstate; k++)	{
			Stationary[mode][k] /= total;
		}
	}	
	else	{	
		int zipsize = 0;
		double total = 0;
		for (int i=0; i<Nstate; i++)	{
			if (ModeAASupport[mode][i])	{
				ModeZipStationary[mode][zipsize] = Stationary[mode][i];
				ModeZipAA[mode][zipsize] = i;
				ModeAAZip[mode][i] = zipsize;
				total += Stationary[mode][i];
				zipsize++;
			}
			else	{
				Stationary[mode][i] = 0;
			}
		}
		ModeZipSize[mode] = zipsize;

		// make zip stat
		for (int i=0; i<zipsize; i++)	{
			ModeZipStationary[mode][i] /= total;
		}
		for (int k=0; k<Nstate; k++)	{
			Stationary[mode][k] /= total;
		}
	}
	
	// make zip rel rates
	int k = 0;
	double* rr = 0;
	if (mParam->Qmode)	{
		rr = RR[mode];
	}
	else	{
		rr = ModeRR;
	}
	
	for (int i=0; i<ModeZipSize[mode]; i++)	{
		for (int j=i+1; j<ModeZipSize[mode]; j++)	{
			int aa = ModeZipAA[mode][i];
			int bb = ModeZipAA[mode][j];
			ModeZipRR[mode][k] = rr[(2 * Nstate - aa - 1) * aa / 2 + bb - aa - 1];
			k++;
		}
	}
}	

int PhyloBayes::ZipGTRCompatible(int site, int mode)	{

	int k = 0;
	int cont = 1;
	while (cont && (k<Nstate))	{
		if (!ModeAASupport[mode][k])	{
			cont = ! (mParam->Orbit[site][k]);
		}
		k++;
	}
	return cont;
}

int PhyloBayes::ZipConstrained(int zipmode, int k)	{
	int site=0;
	while ((site<mParam->Nsite) && (!Nsub[site][k])) site++;
	return (site<mParam->Nsite);
}

// ---------------------------------------------------------------------------
//		 CreateRateMode()
// ---------------------------------------------------------------------------

void PhyloBayes::CreateRateMode()	{

	if (verbose)	{
		cerr << "create rate mode\n";
	}
	// rate mode related structures
	RateSiteNumber = new int[mParam->NRateModeMax];
	RateMode = new int[mParam->Nsite];
	ModeRate = new double[mParam->NRateModeMax];
	for (int i=0; i<mParam->NRateModeMax; i++)	{
		RateSiteNumber[i] = 0;
		ModeRate[i] =1;
	}
	RateModeWeight = new double[mParam->NRateModeMax];
	for (int i=0; i<mParam->NRateModeMax; i++)	{
		RateModeWeight[i] = 1;
	}
	SiteEffLength = new double[mParam->Nsite];
	RateModeEffLength = new double[mParam->NRateModeMax];
	RateAlpha = 1;
	RateModeTotal = new int[mParam->NRateModeMax];

	NRateMode = 1;
	for (int i=0; i<mParam->Nsite; i++)	{
		RateMode[i] = 0;
	}

	ConstantStatus = new int[mParam->Nsite];
	for (int i=0; i<mParam->Nsite; i++)	{
		ConstantStatus[i] = 0;
	}
}

// ---------------------------------------------------------------------------
//		 CreateNH()
// ---------------------------------------------------------------------------

void PhyloBayes::CreateNH()	{

	if (mParam->NH)	{
	if (mParam->HeteroMode == Covarion)	{
		cerr << "cov nh not yet implemented\n";
		exit(1);
	}

	NHNcat = mParam->NHNcatMax;
	NHcat = new int[mParam->Nnode];
	for (int j=0; j<mParam->Nnode; j++)	{
		NHcat[j] = 0;
	}

	if (mParam->PopEff)	{
		PopEffSize = new double[mParam->NHNcatMax];
		for (int j=0; j<mParam->NHNcatMax; j++)	{
			PopEffSize[j] = 1;
		}
		PopAlpha = 1;
		PopBeta = 1;
	}
	NHStatDistorter = new double*[mParam->NHNcatMax];
	for (int n=0; n<mParam->NHNcatMax; n++)	{
		NHStatDistorter[n] = new double[Nstate];
	}
	logNHStatDistorter = new double*[mParam->NHNcatMax];
	for (int n=0; n<mParam->NHNcatMax; n++)	{
		logNHStatDistorter[n] = new double[ContNstate];
	}
	
	NHWeight = new double[mParam->NHNcatMax];
	for (int n=0; n<mParam->NHNcatMax; n++)	{
		NHWeight[n] = 1;
	}
	NHStatCenter = new double[ContNstate];
	for (int k=0; k<ContNstate; k++)	{
		NHStatCenter[k] = 1.0/Nstate;
	}
	NHStatAlpha = Nstate;
	mNHMatrixArray = 0;
	NHZipStationary = 0;
	
	NHStationary = new double**[mParam->NmodeMax];
	for (int i=0; i<mParam->NmodeMax; i++)	{
		NHStationary[i] = new double*[mParam->NHNcatMax];
		for (int n=0; n<mParam->NHNcatMax; n++)	{
			NHStationary[i][n] = new double[Nstate];
		}
	}

	NHLogStationary = new double**[mParam->NmodeMax];
	for (int i=0; i<mParam->NmodeMax; i++)	{
		NHLogStationary[i] = new double*[mParam->NHNcatMax];
		for (int n=0; n<mParam->NHNcatMax; n++)	{
			NHLogStationary[i][n] = new double[Nstate];
		}
	}
	if (mParam->NHPrior)	{
		NHVar = new double[ContNstate];
		NHContVar = new double[mParam->ContNchar + 1];
		for (int k=0; k<mParam->ContNchar+1; k++)	{
			NHContVar[k] = 1;
		}
		if (mParam->NHPrior == 2)	{
			NHCovIndex = new double[ContNrr];
			NHCov = new double*[ContNstate];
			NHInvCov = new double*[ContNstate];
			for (int i=0; i<ContNstate; i++)	{
				NHCov[i] = new double[ContNstate];
				NHInvCov[i] = new double[ContNstate];
			}
		}
	}

	if (mParam->ModeFastCompute)	{
		if (mParam->ModePoisson)	{
			NHZipStationary = new double**[mParam->Nsite];
			for (int i=0; i<mParam->Nsite; i++)	{
				NHZipStationary[i] = new double*[mParam->NHNcatMax];
				for (int n=0; n<mParam->NHNcatMax; n++)	{
					// NHZipStationary[i][n] = new double[mParam->ZipSize[i]];
					NHZipStationary[i][n] = new double[Nstate];
				}
			}
		}
		else	{
			if (mParam->ZipGTR >= 2)	{
				NHZipStationary = new double**[mParam->Nsite];
				for (int i=0; i<mParam->Nsite; i++)	{
					NHZipStationary[i] = new double*[mParam->NHNcatMax];
					for (int n=0; n<mParam->NHNcatMax; n++)	{
						NHZipStationary[i][n] = new double[Nstate];
					}
				}
				
				mNHMatrixArray = new SubMatrix**[mParam->Nsite];
				for (int i=0; i<mParam->Nsite; i++)	{
					double* rr = 0;
					if (mParam->Qmode)	{
						rr = RR[i];
					}
					else	{
						rr = ModeRR;
					}
					mNHMatrixArray[i] = new SubMatrix*[mParam->NHNcatMax];
					for (int n=0; n<mParam->NHNcatMax; n++)	{
						// mNHMatrixArray[i][n] = 0;
						mNHMatrixArray[i][n] = new SubMatrix(mParam,2,rr,NHZipStationary[i][n]);
					}
				}
			}
			else if (mParam->ZipGTR == 1)	{
				NHModeZipStationary = new double**[mParam->NmodeMax];
				for (int i=0; i<mParam->NmodeMax; i++)	{
					NHModeZipStationary[i] = new double*[mParam->NHNcatMax];
					for (int n=0; n<mParam->NHNcatMax; n++)	{
						NHModeZipStationary[i][n] = new double[Nstate];
					}
				}
				
				mNHMatrixArray = new SubMatrix**[mParam->NmodeMax];
				for (int i=0; i<mParam->NmodeMax; i++)	{
					double* rr = 0;
					if (mParam->Qmode)	{
						rr = RR[i];
					}
					else	{
						rr = ModeRR;
					}
					mNHMatrixArray[i] = new SubMatrix*[mParam->NHNcatMax];
					for (int n=0; n<mParam->NHNcatMax; n++)	{
						mNHMatrixArray[i][n] = new SubMatrix(mParam,Nstate,rr,NHModeZipStationary[i][n]);
					}
				}
			}
		}
	}
	else	{
		mNHMatrixArray = new SubMatrix**[mParam->NmodeMax];
		for (int i=0; i<mParam->NmodeMax; i++)	{
			double* rr = 0;
			if (mParam->Qmode)	{
				rr = RR[i];
			}
			else	{
				rr = ModeRR;
			}
			mNHMatrixArray[i] = new SubMatrix*[mParam->NHNcatMax];
			for (int n=0; n<mParam->NHNcatMax; n++)	{
				mNHMatrixArray[i][n] = new SubMatrix(mParam,Nstate,rr,NHStationary[i][n]);
			}
		}
	}
	}

}	

void PhyloBayes::CreateNHTotal()	{

	NHTotal = new int**[Nmode];
	for (int i=0; i<Nmode; i++)	{
		NHTotal[i] = new int*[NHNcat];
		for (int n=0; n<NHNcat; n++)	{
			NHTotal[i][n] = new int[Nstate];
		}
	}
	if (! mParam->ModePoisson)	{
		NHModeStatBeta = new double**[Nmode];
		for (int i=0; i<Nmode; i++)	{
			NHModeStatBeta[i] = new double*[NHNcat];
			for (int n=0; n<NHNcat; n++)	{
				NHModeStatBeta[i][n] = new double[Nstate];
			}
		}
	}
}

void PhyloBayes::DeleteNHTotal()	{

	for (int i=0; i<Nmode; i++)	{
		for (int n=0; n<NHNcat; n++)	{
			delete[] NHTotal[i][n];
		}
		delete[] NHTotal[i];
	}
	delete[] NHTotal;
	NHTotal = 0;

	if (! mParam->ModePoisson)	{
		for (int i=0; i<Nmode; i++)	{
			for (int n=0; n<NHNcat; n++)	{
				delete[] NHModeStatBeta[i][n];
			}
			delete[] NHModeStatBeta[i];
		}
		delete[] NHModeStatBeta;
		NHModeStatBeta = 0;
	}
}

void PhyloBayes::CreateNHSub()	{

	if (! NHNsub)	{
		NHNsub = new int**[mParam->Nnode];
		for (int j=0; j<mParam->Nnode; j++)	{
			NHNsub[j] = new int*[mParam->Nsite];
			for (int i=0; i<mParam->Nsite; i++)	{
				NHNsub[j][i] = new int[Nstate];
			}
		}
		NHStatBeta = new double**[mParam->Nnode];
		for (int j=0; j<mParam->Nnode; j++)	{
			NHStatBeta[j] = new double*[mParam->Nsite];
			for (int i=0; i<mParam->Nsite; i++)	{
				NHStatBeta[j][i] = new double[Nstate];
			}
		}
	}
}

void PhyloBayes::DeleteNHSub()	{

	if (NHNsub)	{
		for (int j=0; j<mParam->Nnode; j++)	{
			for (int i=0; i<mParam->Nsite; i++)	{
				delete[] NHNsub[j][i];
			}
			delete[] NHNsub[j];
		}
		delete[] NHNsub;
		NHNsub = 0;
	}
	if (NHStatBeta)	{
		for (int j=0; j<mParam->Nnode; j++)	{
			for (int i=0; i<mParam->Nsite; i++)	{
				delete[] NHStatBeta[j][i];
			}
			delete[] NHStatBeta[j];
		}
		delete[] NHStatBeta;
		NHStatBeta = 0;
	}
}

void PhyloBayes::CreateNHBranch()	{

	NHBranchNsub = new int**[mParam->Nnode];
	for (int j=0; j<mParam->Nnode; j++)	{
		NHBranchNsub[j] = new int*[Nmode];
		for (int mode=0; mode<Nmode; mode++)	{
			NHBranchNsub[j][mode] = new int[Nstate];
			for (int k=0; k<Nstate; k++)	{
				NHBranchNsub[j][mode][k] = 0;
			}
		}
		for (int i=0; i<mParam->Nsite; i++)	{
			int mode = Mode[i];
			for (int k=0; k<Nstate; k++)	{
				NHBranchNsub[j][mode][k] += NHNsub[j][i][k];
			}
		}
	}

	NHBranchStatBeta = 0;
	if (!mParam->ModePoisson)	{
		NHBranchStatBeta = new double**[mParam->Nnode];
		for (int j=0; j<mParam->Nnode; j++)	{
			NHBranchStatBeta[j] = new double*[Nmode];
			for (int mode=0; mode<Nmode; mode++)	{
				NHBranchStatBeta[j][mode] = new double[Nstate];
				for (int k=0; k<Nstate; k++)	{
					NHBranchStatBeta[j][mode][k] = 0;
				}
			}
			for (int i=0; i<mParam->Nsite; i++)	{
				int mode = Mode[i];
				for (int k=0; k<Nstate; k++)	{
					NHBranchStatBeta[j][mode][k] += NHStatBeta[j][i][k];
				}
			}
		}
	}
}

void PhyloBayes::DeleteNHBranch()	{

	for (int j=0; j<mParam->Nnode; j++)	{
		for (int mode=0; mode<Nmode; mode++)	{
			delete[] NHBranchNsub[j][mode];
		}
		delete[] NHBranchNsub[j];
	}
	delete[] NHBranchNsub;
	NHBranchNsub = 0;

	if (!mParam->ModePoisson)	{
		for (int j=0; j<mParam->Nnode; j++)	{
			for (int mode=0; mode<Nmode; mode++)	{
				delete[] NHBranchStatBeta[j][mode];
			}
			delete[] NHBranchStatBeta[j];
		}
		delete[] NHBranchStatBeta;
	}
	NHBranchStatBeta = 0;
}


void PhyloBayes::CreateNHSite()	{

	NHSiteNsub = new int**[mParam->Nsite];
	for (int i=0; i<mParam->Nsite; i++)	{
		NHSiteNsub[i] = new int*[NHNcat];
		for (int n=0; n<NHNcat; n++)	{
			NHSiteNsub[i][n] = new int[Nstate];
			for (int k=0; k<Nstate; k++)	{
				NHSiteNsub[i][n][k] = 0;
			}
		}
		for (int j=0; j<mParam->Nnode; j++)	{
			int n = NHcat[j];
			for (int k=0; k<Nstate; k++)	{
				NHSiteNsub[i][n][k] += NHNsub[j][i][k];
			}
		}
	}

	if (!mParam->ModePoisson)	{
		NHSiteStatBeta = new double**[mParam->Nsite];
		for (int i=0; i<mParam->Nsite; i++)	{
			NHSiteStatBeta[i] = new double*[NHNcat];
			for (int n=0; n<NHNcat; n++)	{
				NHSiteStatBeta[i][n] = new double[Nstate];
				for (int k=0; k<Nstate; k++)	{
					NHSiteStatBeta[i][n][k] = 0;
				}
			}
			for (int j=0; j<mParam->Nnode; j++)	{
				int n = NHcat[j];
				for (int k=0; k<Nstate; k++)	{
					NHSiteStatBeta[i][n][k] += NHStatBeta[j][i][k];
				}
			}
		}
	}
}

void PhyloBayes::DeleteNHSite()	{

	for (int i=0; i<mParam->Nsite; i++)	{
		for (int n=0; n<NHNcat; n++)	{
			delete[] NHSiteNsub[i][n];
		}
		delete[] NHSiteNsub[i];
	}
	delete[] NHSiteNsub;
	NHSiteNsub = 0;

	if (!mParam->ModePoisson)	{
		for (int i=0; i<mParam->Nsite; i++)	{
			for (int n=0; n<NHNcat; n++)	{
				delete[] NHSiteStatBeta[i][n];
			}
			delete[] NHSiteStatBeta[i];
		}
		delete[] NHSiteStatBeta;
		NHSiteStatBeta = 0;
	}
}

// ---------------------------------------------------------------------------
//		 DeleteNH()
// ---------------------------------------------------------------------------

void PhyloBayes::DeleteNH()	{

	if (mParam->NH)	{
	if (mParam->HeteroMode == Covarion)	{
		cerr << "cov nh not yet implemented\n";
		exit(1);
	}

	if (mParam->PopEff)	{
		delete[] PopEffSize;
	}

	delete[] NHcat;
	for (int n=0; n<mParam->NHNcatMax; n++)	{
		delete[] NHStatDistorter[n];
	}
	for (int n=0; n<mParam->NHNcatMax; n++)	{
		delete[] logNHStatDistorter[n];
	}
	delete[] logNHStatDistorter;
	
	delete[] NHStatDistorter;
	delete[] NHWeight;
	delete[] NHStatCenter;
	
	for (int i=0; i<mParam->NmodeMax; i++)	{
		for (int n=0; n<mParam->NHNcatMax; n++)	{
			delete[] NHStationary[i][n];
		}
		delete[] NHStationary[i];
	}
	delete[] NHStationary;

	for (int i=0; i<mParam->NmodeMax; i++)	{
		for (int n=0; n<mParam->NHNcatMax; n++)	{
			delete[] NHLogStationary[i][n];
		}
		delete[] NHLogStationary[i];
	}
	delete[] NHLogStationary;

	if (mParam->NHPrior)	{
		delete[] NHVar;
		delete[] NHContVar;
		if (mParam->NHPrior == 2)	{
			for (int i=0; i<ContNstate; i++)	{
				delete[] NHCov[i];
				delete[] NHInvCov[i];
			}
			delete[] NHCov;
			delete[] NHInvCov;
			delete[] NHCovIndex;
		}
	}
	if (mParam->ModeFastCompute)	{
		if (mParam->ModePoisson)	{
			for (int i=0; i<mParam->Nsite; i++)	{
				for (int n=0; n<mParam->NHNcatMax; n++)	{
					delete[] NHZipStationary[i][n];
				}
				delete[] NHZipStationary[i];
			}
			delete[] NHZipStationary;
		}
		else	{
			if (mParam->ZipGTR >= 2)	{
				for (int i=0; i<mParam->Nsite; i++)	{
					for (int n=0; n<mParam->NHNcatMax; n++)	{
						delete[] NHZipStationary[i][n];
					}
					delete[] NHZipStationary[i];
				}
				delete[] NHZipStationary;
				for (int i=0; i<mParam->Nsite; i++)	{
					for (int n=0; n<mParam->NHNcatMax; n++)	{
						delete mNHMatrixArray[i][n];
					}
					delete[] mNHMatrixArray[i];
				}
				delete[] mNHMatrixArray;
			}
			else	{
				for (int i=0; i<mParam->NmodeMax; i++)	{
					for (int n=0; n<mParam->NHNcatMax; n++)	{
						delete[] NHModeZipStationary[i][n];
					}
					delete[] NHModeZipStationary[i];
				}
				delete[] NHModeZipStationary;
				for (int i=0; i<mParam->NmodeMax; i++)	{
					for (int n=0; n<mParam->NHNcatMax; n++)	{
						delete mNHMatrixArray[i][n];
					}
					delete[] mNHMatrixArray[i];
				}
				delete[] mNHMatrixArray;
			}
		}
	}
	else	{
		for (int i=0; i<mParam->NmodeMax; i++)	{
			for (int n=0; n<mParam->NHNcatMax; n++)	{
				delete[] mNHMatrixArray[i][n];
			}
			delete[] mNHMatrixArray[i];
		}
		delete[] mNHMatrixArray;
	}
	}
}	

// ---------------------------------------------------------------------------
//		 CloneNH(PhyloBayes* from)
// ---------------------------------------------------------------------------

void PhyloBayes::CloneNH(PhyloBayes* from)	{

	if (mParam->NH)	{
		/*
		if (mParam->PopEff)	{
			for (int j=0; j<mParam->Nnode; j++)	{
				PopEffSize[j] = from->PopEffSize[j];
			}
			PopAlpha = from->PopAlpha;
		}
		*/
		PopAlpha = from->PopAlpha;
		PopBeta = from->PopBeta;
		NHpswitch = from->NHpswitch;
		NHNcat = from->NHNcat;
		for (int j=0; j<mParam->Nnode; j++)	{
			NHcat[j] = from->NHcat[j];
		}

		NHStatAlpha = from->NHStatAlpha;
		for (int k=0; k<ContNstate; k++)	{
			NHStatCenter[k] = from->NHStatCenter[k];
		}

		if (mParam->NHPrior)	{
			for (int k=0; k<ContNstate; k++)	{
				NHVar[k] = from->NHVar[k];
			}
			for (int k=0; k<mParam->ContNchar+1; k++)	{
				NHContVar[k] = from->NHContVar[k];
			}
			if (mParam->NHPrior == 2)	{
				for (int k=0; k<ContNrr; k++)	{
					NHCovIndex[k] = from->NHCovIndex[k];
				}
				for (int k=0; k<ContNstate; k++)	{
					for (int l=0; l<ContNstate; l++)	{
						NHCov[k][l] = from->NHCov[k][l];
						NHInvCov[k][l] = from->NHInvCov[k][l];
					}
				}
				NHInvertible = from->NHInvertible;
				NHLogDet = from->NHLogDet;
			}
		}

		for (int n=0; n<NHNcat; n++)	{
			CloneNH(from,n);
		}
	}
}

void PhyloBayes::ClonePopEff(PhyloBayes* from)	{
	for (int n=0; n<NHNcat; n++)	{
		PopEffSize[n] = from->PopEffSize[n];
	}
}

void PhyloBayes::CloneNH(PhyloBayes* from, int n)	{
	if (mParam->PopEff)	{
		PopEffSize[n] = from->PopEffSize[n];
	}
	NHWeight[n] = from->NHWeight[n];
	for (int k=0; k<Nstate; k++)	{
		NHStatDistorter[n][k] = from->NHStatDistorter[n][k];
	}
	if (mParam->NHPrior)	{
		for (int k=0; k<ContNstate; k++)	{
			logNHStatDistorter[n][k] = from->logNHStatDistorter[n][k];
		}
	}
	for (int mode=0; mode<Nmode; mode++)	{
		CloneNH(from,mode,n);
	}
}

void PhyloBayes::CloneNH(PhyloBayes* from, int mode, int n)	{

	CloneNHStat(from,mode,n);
	if (mParam->ModePoisson)	{
		for (int i=0; i<mParam->Nsite; i++)	{
			if (Mode[i] == mode)	{
				CloneZipNH(from,i,n);
			}
		}
	}
	else	{
		if (mParam->ModeFastCompute)	{
			if (mParam->ZipGTR >= 2)	{
				for (int i=0; i<mParam->Nsite; i++)	{
					if (Mode[i] == mode)	{
						mNHMatrixArray[i][n]->Copy(from->mNHMatrixArray[i][n]);
					}
				}
			}
			else	{
				mNHMatrixArray[mode][n]->Copy(from->mNHMatrixArray[mode][n]);
			}
		}
		else	{
			mNHMatrixArray[mode][n]->Copy(from->mNHMatrixArray[mode][n]);
		}
	}
}
	
void PhyloBayes::CloneNHStat(PhyloBayes* from, int mode, int n)	{
	for (int k=0; k<Nstate; k++)	{
		NHStationary[mode][n][k] = from->NHStationary[mode][n][k];
		NHLogStationary[mode][n][k] = from->NHLogStationary[mode][n][k];
	}
}

void PhyloBayes::CloneZipNH(PhyloBayes* from, int site, int n)	{
	for (int k=0; k<mParam->ZipSize[site]; k++)	{
		NHZipStationary[site][n][k] = from->NHZipStationary[site][n][k];
	}
}

// ---------------------------------------------------------------------------
//		 UpdateNH()
// ---------------------------------------------------------------------------

void PhyloBayes::ComputeNHStat()	{
	for (int mode=0; mode<Nmode; mode++)	{
		for (int n=0; n<NHNcat; n++)	{
			ComputeNHStat(mode,n);
		}
	}
}

void PhyloBayes::ComputeNHStat(int n)	{
	for (int mode=0; mode<Nmode; mode++)	{
		ComputeNHStat(mode,n);
	}
}

void PhyloBayes::ComputeNHModeStat(int mode)	{
	for (int n=0; n<NHNcat; n++)	{
		ComputeNHStat(mode,n);
	}
}

void PhyloBayes::UpdateNH()	{

	if (mParam->NH)	{
		if (mParam->NHPrior)	{
			ComputeNHStatDistorter();
			if (mParam->NHPrior == 2)	{
				ComputeNHCov();
			}
		}
		for (int mode=0; mode<Nmode; mode++)	{
			for (int n=0; n<NHNcat; n++)	{
				UpdateNH(mode,n);
			}
		}
	}
}

void PhyloBayes::ComputeNHStatDistorter()	{

	for (int n=0; n<NHNcat; n++)	{
		ComputeNHStatDistorter(n);
	}
}

void PhyloBayes::ComputeNHStatDistorter(int n)	{

	double total = 0;
	if ((mParam->NH >=4) && (mParam->NHHalfSum))	{
		int j = 0;
		while ((j<mParam->Nnode) && (NHcat[j] != n))	{
			j++;
		}
		if (j == mParam->Nnode)	{
			cerr << "error in compute nh stat distorter\n";
			exit(1);
		}
		if (tree[j].isRoot())	{
			for (int k=0; k<Nstate; k++)	{
				NHStatDistorter[n][k] = exp(logNHStatDistorter[n][k]);
				total += NHStatDistorter[n][k];
			}
		}
		else	{
			int nup = NHcat[tree[j].up->label];
			for (int k=0; k<Nstate; k++)	{
				NHStatDistorter[n][k] = exp(0.5 * (logNHStatDistorter[n][k] + logNHStatDistorter[nup][k]));
				total += NHStatDistorter[n][k];
			}
		}
	}
	else	{
		for (int k=0; k<Nstate; k++)	{
			NHStatDistorter[n][k] = exp(logNHStatDistorter[n][k]);
			total += NHStatDistorter[n][k];
		}
	}
	for (int k=0; k<Nstate; k++)	{
		NHStatDistorter[n][k] /= total;
	}
}

void PhyloBayes::ComputeNHCov()	{

	cerr << "compute NH Cov deprecated\n";
	exit(1);
}



void PhyloBayes::UpdateNH(int n)	{
	for (int mode=0; mode<Nmode; mode++)	{
		UpdateNH(mode,n);
	}
}

void PhyloBayes::UpdateModeNH(int mode)	{
	for (int n=0; n<NHNcat; n++)	{
		UpdateNH(mode,n);
	}
}

void PhyloBayes::UpdateNH(int mode, int n)	{

	ComputeNHStat(mode,n);

	if (mParam->ModeFastCompute)	{
		if (mParam->ModePoisson)	{
			for (int i=0; i<mParam->Nsite; i++)	{
				if (Mode[i] == mode)	{
					UpdateZipNH(i,n);
				}
			}
		}
		else	{
			if (mParam->ZipGTR == 1)	{
				NHModeZip(mode,n);
				mNHMatrixArray[mode][n]->ComputeArray();
				
			}
			else if (mParam->ZipGTR >= 2)	{
				for (int i=0; i<mParam->Nsite; i++)	{
					if (Mode[i] == mode)	{
						NHSiteZip(i,n);
						mNHMatrixArray[i][n]->ComputeArray();
					}
				}
			}
		}
	}
	else	{
		mNHMatrixArray[mode][n]->ComputeArray();
	}
}
		
void PhyloBayes::ComputeNHStat(int mode, int n)	{

	// compute eq frequencies
	double total = 0;
	if (mParam->PopEff && ! mParam->NH)	{
		for (int k=0; k<Nstate; k++)	{
			NHStationary[mode][n][k] = exp(PopEffSize[n] * log(Stationary[mode][k]));
			total += NHStationary[mode][n][k];
		}
		for (int k=0; k<Nstate; k++)	{
			NHStationary[mode][n][k] /= total;
			NHLogStationary[mode][n][k] = log(NHStationary[mode][n][k]);
		}
	}
	else	{
		if (mParam->PopEff)	{
			for (int k=0; k<Nstate; k++)	{
				NHStationary[mode][n][k] = exp(PopEffSize[n] * log(Stationary[mode][k])) * NHStatDistorter[n][k];
				total += NHStationary[mode][n][k];
			}
		}
		else	{
			for (int k=0; k<Nstate; k++)	{
				NHStationary[mode][n][k] = Stationary[mode][k] * NHStatDistorter[n][k];
				total += NHStationary[mode][n][k];
			}
		}
		for (int k=0; k<Nstate; k++)	{
			NHStationary[mode][n][k] /= total;
			NHLogStationary[mode][n][k] = log(NHStationary[mode][n][k]);
		}
	}
}

// ---------------------------------------------------------------------------
//		 CreateRef()
// ---------------------------------------------------------------------------

void PhyloBayes::CreateRef()	{

	if (verbose)	{
		cerr << "create ref\n";
	}
	RefRR = new double[Nrr];
	RefStationary = new double[Nstate];
	RefZipStationary = 0;
	mSubMatrix = 0;
	RefRateFactor = 1;
	RefZipStationary = new double*[mParam->Nsite];
	for (int i=0; i<mParam->Nsite; i++)	{
		RefZipStationary[i] = 0;
	}

	if (! mParam->RefFastCompute)	{
		// allocating slow computation related structures (matrices)
		mSubMatrix = new SubMatrix(mParam,Nstate,RefRR,RefStationary);
	}
	else	{

		for (int i=0; i<mParam->Nsite; i++)	{
			RefZipStationary[i] = new double[Nstate];
		}
	}
}


// ---------------------------------------------------------------------------
//		 CreatePartial()
// ---------------------------------------------------------------------------

void PhyloBayes::CreatePartial()	{

	if (mParam->Tempering)	{
		NPartialMode = Nmode;
	}
	else	{
		NPartialMode = 1;
	}
	if (mParam->SumOverRateModes)	{
		NPartialRateMode = NRateMode;
	}
	else	{
		NPartialRateMode = 1;
	}
	Partialpdown = new double****[mParam->Nnode];
	Partialqup = new double****[mParam->Nnode];
	for (int i=0; i<mParam->Nnode; i++)	{
		Partialpdown[i] = new double***[mParam->Nsite];
		Partialqup[i] = new double***[mParam->Nsite];
		for (int j=0; j<mParam->Nsite; j++)	{
			Partialpdown[i][j] = new double**[NPartialMode];
			Partialqup[i][j] = new double**[NPartialMode];
			for (int k=0; k<NPartialMode; k++)	{
				Partialpdown[i][j][k] = new double*[NPartialRateMode];
				Partialqup[i][j][k] = new double*[NPartialRateMode];
				for (int l=0; l<NPartialRateMode; l++)	{
					if (mParam->ModeFastCompute)	{
						if (mParam->ModePoisson)	{
							if (mParam->HeteroMode == Covarion) {
								Partialpdown[i][j][k][l] = new double[2*mParam->ZipSize[j] + 1];
								Partialqup[i][j][k][l] = new double[2*mParam->ZipSize[j] + 1];
							}
							else	{
								Partialpdown[i][j][k][l] = new double[mParam->ZipSize[j] + 1];
								Partialqup[i][j][k][l] = new double[mParam->ZipSize[j] + 1];
							}
						}
						else	{
							int nstate = 0;
							if (mParam->ZipGTR == 4)	{
								nstate = mParam->ZipSize[j];
							}
							else if (mParam->ZipGTR >= 2)	{
								nstate = SiteZipSize[j];
							}
							else	{
								nstate = ModeZipSize[Mode[j]];
							}
							Partialpdown[i][j][k][l] = new double[nstate + 1];
							Partialqup[i][j][k][l] = new double[nstate + 1];
						}
					}
					else	{
						if (mParam->HeteroMode == Covarion) {
							Partialpdown[i][j][k][l] = new double[2*Nstate + 1];
							Partialqup[i][j][k][l] = new double[2*Nstate + 1];
						}
						else	{
							Partialpdown[i][j][k][l] = new double[Nstate + 1];
							Partialqup[i][j][k][l] = new double[Nstate + 1];
						}
					}
				}
			}
		}
	}

	Partialqdownleft = new double***[mParam->Nsite];
	Partialqdownright = new double***[mParam->Nsite];
	Partialpup = new double***[mParam->Nsite];
	for (int j=0; j<mParam->Nsite; j++)	{
		Partialqdownleft[j] = new double**[NPartialMode];
		Partialqdownright[j] = new double**[NPartialMode];
		Partialpup[j] = new double**[NPartialMode];
		for (int k=0; k<NPartialMode; k++)	{
			Partialqdownleft[j][k] = new double*[NPartialRateMode];
			Partialqdownright[j][k] = new double*[NPartialRateMode];
			Partialpup[j][k] = new double*[NPartialRateMode];
			for (int l=0; l<NPartialRateMode; l++)	{
				if (mParam->ModeFastCompute)	{
					if (mParam->ModePoisson)	{
						if (mParam->HeteroMode == Covarion) {
							Partialqdownleft[j][k][l] = new double[2*mParam->ZipSize[j] + 1];
							Partialqdownright[j][k][l] = new double[2*mParam->ZipSize[j] + 1];
							Partialpup[j][k][l] = new double[2*mParam->ZipSize[j] + 1];
						}
						else	{
							Partialqdownleft[j][k][l] = new double[mParam->ZipSize[j] + 1];
							Partialqdownright[j][k][l] = new double[mParam->ZipSize[j] + 1];
							Partialpup[j][k][l] = new double[mParam->ZipSize[j] + 1];
						}
					}
					else	{
						int nstate = 0;
						if (mParam->ZipGTR == 4)	{
							nstate = mParam->ZipSize[j];
						}
						else if (mParam->ZipGTR >= 2)	{
							nstate = SiteZipSize[j];
						}
						else	{
							nstate = ModeZipSize[Mode[j]];
						}
						Partialqdownleft[j][k][l] = new double[nstate + 1];
						Partialqdownright[j][k][l] = new double[nstate + 1];
						Partialpup[j][k][l] = new double[nstate + 1];
					}
				}
				else	{
					if (mParam->HeteroMode == Covarion) {
						Partialqdownleft[j][k][l] = new double[2*Nstate + 1];
						Partialqdownright[j][k][l] = new double[2*Nstate + 1];
						Partialpup[j][k][l] = new double[2*Nstate + 1];
					}
					else	{
						Partialqdownleft[j][k][l] = new double[Nstate + 1];
						Partialqdownright[j][k][l] = new double[Nstate + 1];
						Partialpup[j][k][l] = new double[Nstate + 1];
					}
				}
			}
		}
	}

	PartialUpdateFlag = new int[mParam->Nnode];

	PartialBaseRate = new double***[mParam->Nsite];
	PartialBaseStationary = new double***[mParam->Nsite];
	PartialBaseMatrix = new SubMatrix***[mParam->Nsite];
	for (int i=0; i<mParam->Nsite; i++)	{
		PartialBaseRate[i] = new double**[NPartialMode];
		PartialBaseStationary[i] = new double**[NPartialMode];
		PartialBaseMatrix[i] = new SubMatrix**[NPartialMode];
		for (int k=0; k<NPartialMode; k++)	{
			PartialBaseRate[i][k] = new double*[NPartialRateMode];
			PartialBaseStationary[i][k] = new double*[NPartialRateMode];
			PartialBaseMatrix[i][k] = new SubMatrix*[NPartialRateMode];
		}
	}
	ResetPartialUpdateFlags();
}


// ---------------------------------------------------------------------------
//		 DeletePartial()
// ---------------------------------------------------------------------------

void PhyloBayes::DeletePartial()	{

	for (int i=0; i<mParam->Nnode; i++)	{
		for (int j=0; j<mParam->Nsite; j++)	{
			for (int k=0; k<NPartialMode; k++)	{
				for (int l=0; l<NPartialRateMode; l++)	{
					delete[] Partialpdown[i][j][k][l];
					delete[] Partialqup[i][j][k][l];
				}
				delete[] Partialpdown[i][j][k];
				delete[] Partialqup[i][j][k];
			}
			delete[] Partialpdown[i][j];
			delete[] Partialqup[i][j];
		}
		delete[] Partialpdown[i];
		delete[] Partialqup[i];
	}
	delete[] Partialpdown;
	Partialpdown = 0;
	delete[] Partialqup;

	for (int j=0; j<mParam->Nsite; j++)	{
		for (int k=0; k<NPartialMode; k++)	{
			for (int l=0; l<NPartialRateMode; l++)	{
				delete[] Partialqdownleft[j][k][l];
				delete[] Partialqdownright[j][k][l];
				delete[] Partialpup[j][k][l];
			}
			delete[] Partialqdownleft[j][k];
			delete[] Partialqdownright[j][k];
			delete[] Partialpup[j][k];
		}
		delete[] Partialqdownleft[j];
		delete[] Partialqdownright[j];
		delete[] Partialpup[j];
	}
	delete[] Partialqdownleft;
	delete[] Partialqdownright;
	delete[] Partialpup;

	delete[] PartialUpdateFlag;

	for (int i=0; i<mParam->Nsite; i++)	{
		for (int k=0; k<NPartialMode; k++)	{
			delete[] PartialBaseRate[i][k];
			delete[] PartialBaseStationary[i][k];
			delete[] PartialBaseMatrix[i][k];
		}
		delete[] PartialBaseRate[i];
		delete[] PartialBaseStationary[i];
		delete[] PartialBaseMatrix[i];
	}
	delete[] PartialBaseRate;
	delete[] PartialBaseStationary;
	delete[] PartialBaseMatrix;
}


// ---------------------------------------------------------------------------
//		 ResetPartialUpdateFlags()
// ---------------------------------------------------------------------------

void PhyloBayes::ResetPartialUpdateFlags()	{

	for (int j=0; j<mParam->Nnode; j++)	{
		PartialUpdateFlag[j] = 0;
	}
}


// ---------------------------------------------------------------------------
//		 UpdatePartialBaseFields()
// ---------------------------------------------------------------------------

void PhyloBayes::UpdatePartialBaseFields()	{

	for (int i=0; i<mParam->Nsite; i++)	{
		if (mParam->Tempering)	{
			RegisterTempZipStationary();
			for (int k=0; k<NRateMode; k++)	{
				for (int l=0; l<Nmode; l++)	{
					PartialBaseRate[i][l][k] = &ModeRate[k];
					if (mParam->ModeFastCompute)	{
						if (mParam->ModePoisson)	{
							PartialBaseMatrix[i][l][k] = 0;
							PartialBaseStationary[i][l][k] = TempZipStationary[i][l];
						}
						else	{
							cerr << "fast partial gtr not implemented\n";
							exit(1);
						}
					}
					else	{
						PartialBaseMatrix[i][l][k] = mMatrixArray[l];
						PartialBaseStationary[i][l][k] = PartialBaseMatrix[i][l][k]->GetStationaries();
					}
				}
			}
		}
		else	{
			
		for (int k=0; k<NPartialRateMode; k++)	{
			if (mParam->HeteroMode == Covarion)	{
				if (mParam->ModeFastCompute)	{
					if (mParam->ExternalCov)	{
						PartialBaseRate[i][0][k] = &ModeRate[k];
						PartialBaseMatrix[i][0][k] = mCovZipMatrixArray[i][0];
					}
					else	{
						PartialBaseRate[i][0][k] = &One;
						PartialBaseMatrix[i][0][k] = mCovZipMatrixArray[i][k];
					}
				}
				else	{
					if (mParam->ExternalCov)	{
						PartialBaseRate[i][0][k] = &ModeRate[k];
						PartialBaseMatrix[i][0][k] = mCovMatrixArray[Mode[i]][0];
					}
					else	{
						PartialBaseRate[i][0][k] = &One;
						PartialBaseMatrix[i][0][k] = mCovMatrixArray[Mode[i]][k];
					}
				}
				PartialBaseStationary[i][0][k] = PartialBaseMatrix[i][0][k]->GetStationaries();
			}
			else	{
				if (mParam->SumOverRateModes)	{
					PartialBaseRate[i][0][k] = &ModeRate[k];
				}
				else	{
					PartialBaseRate[i][0][k] = &unirate[i];
				}
				if (mParam->ModeFastCompute)	{
					if (mParam->ModePoisson)	{
						PartialBaseMatrix[i][0][k] = 0;
						PartialBaseStationary[i][0][k] = ZipStationary[i];
					}
					else	{
						if (mParam->ZipGTR >= 2)	{
							if (mParam->NH)	{
								PartialBaseMatrix[i][0][k] = mNHMatrixArray[i][NHcat[root->label]];
							}
							else	{
								PartialBaseMatrix[i][0][k] = mZipMatrixArray[i];
							}
						}
						else	{
							if (mParam->NH)	{
								PartialBaseMatrix[i][0][k] = mNHMatrixArray[Mode[i]][NHcat[root->label]];
							}
							else	{
								PartialBaseMatrix[i][0][k] = mMatrixArray[Mode[i]];
							}
						}
						PartialBaseStationary[i][0][k] = PartialBaseMatrix[i][0][k]->GetStationaries();
					}
				}
				else	{
					if (mParam->NH)	{
						PartialBaseMatrix[i][0][k] = mNHMatrixArray[Mode[i]][NHcat[root->label]];
					}
					else	{
						PartialBaseMatrix[i][0][k] = mMatrixArray[Mode[i]];
					}
					PartialBaseStationary[i][0][k] = PartialBaseMatrix[i][0][k]->GetStationaries();
				}
			}
		}
		}
	}
}

// ---------------------------------------------------------------------------
//		 CreateLogProbs()
// ---------------------------------------------------------------------------

void PhyloBayes::CreateLogProbs()	{

	mLogSampling = -1;
	mLogPrior = -1;
	mLogPosterior = -1;

	if (verbose)	{
		cerr << "create log probs \n";
	}
	
	if (mParam->NormalApprox)	{
		if (mParam->Ngene)	{
			mFlexGeneLogSampling = new double[mParam->Ngene];
			mRigidGeneLogSampling = new double[mParam->Ngene];
			mGeneLogSampling = new double[mParam->Ngene];
		}
	}
	else	{
		// structures related to the computation of the log probs
		mSiteLogSampling = new double[mParam->Nsite];
		mSiteLogSampling0 = new double[mParam->Nsite];
		for (int i=0; i<mParam->Nsite; i++)	{
			mSiteLogSampling[i] = 0;
			mSiteLogSampling0[i] = 0;
		}
		mModeSiteLogSampling = 0;
		if (mParam->ActivateSumOverModes)	{
			mModeSiteLogSampling = new double*[mParam->NmodeMax];
			for (int i=0; i<mParam->NmodeMax; i++)	{
				mModeSiteLogSampling[i] = new double[mParam->Nsite];
				for (int j=0; j<mParam->Nsite; j++)	{
					mModeSiteLogSampling[i][j] = 0;
				}
			}
		}
		mRateModeSiteLogSampling = 0;
		if (mParam->ActivateSumOverRateModes)	{
			mRateModeSiteLogSampling = new double*[mParam->NRateModeMax];
			for (int i=0; i<mParam->NRateModeMax; i++)	{
				mRateModeSiteLogSampling[i] = new double[mParam->Nsite];
				for (int j=0; j<mParam->Nsite; j++)	{
					mRateModeSiteLogSampling[i][j] = 0;
				}
			}
		}

		tbl = new double*[mParam->Nnode];
		tbloffset = new double[mParam->Nnode];
		for (int i=0; i<mParam->Nnode; i++)	{
			tbloffset[i] = 0;
		}

		for (int i=0; i<mParam->Nnode; i++)	{
			if (mParam->HeteroMode == Covarion)	{
				tbl[i] = new double[2*Nstate];
			}
			else	{
				tbl[i] = new double[Nstate];
			}
		}
		if (mParam->HeteroMode == Covarion)	{
			aux = new double[2*Nstate];
			taux = new double[2*Nstate];
			taux2 = new double[2*Nstate];
		}
		else	{
			aux = new double[Nstate];
			taux = new double[Nstate];
			taux2 = new double[Nstate];
		}

		if (mParam->DCM)	{
			DCMtbl = new double****[mParam->Nnode];  // [node][mode][ratemode][site][state]
			for (int j=0; j<mParam->Nnode; j++)	{
				DCMtbl[j] = 0;
			}
			DCMFlag = new int[mParam->Nnode];
			for (int j=0; j<mParam->Nnode; j++)	{
				DCMFlag[j] = 0;
			}
		}
		if (mParam->ZipSub == 2)	{
			MissingFlag = new int*[mParam->Nnode];
			for (int i=0; i<mParam->Nnode; i++)	{
				MissingFlag[i] = new int[mParam->Nsite];
			}
		}
	}
}

double** PhyloBayes::CreateLocaltbl()	{

	double** t= new double*[mParam->Nsite];
	for (int i=0; i<mParam->Nsite; i++)	{
		int nstate = Nstate;
		if (mParam->ModeFastCompute)	{
			nstate = mParam->ZipSize[i];
		}	
		t[i] = new double[nstate];
	}
	return t;
}

void PhyloBayes::DeleteLocaltbl(double** t)	{

	for (int i=0; i<mParam->Nsite; i++)	{
		delete[] t[i];
	}
	delete[] t;
	t = 0;
}


void PhyloBayes::CreateNodeDCMtbl(int j)	{

	int nmode = 1;
	if (mParam->ActivateSumOverModes)	{
		nmode = Nmode;
	}
	int nratemode = 1;
	if (mParam->ActivateSumOverRateModes)	{
		nratemode = NRateMode;
	}
	DCMtbl[j] = new double***[nmode];
	for (int mode=0; mode<nmode; mode++)	{
		DCMtbl[j][mode] = new double**[nratemode];
		for (int ratemode=0; ratemode<nratemode; ratemode++)	{
			DCMtbl[j][mode][ratemode] = new double*[mParam->Nsite];
			for (int i=0; i<mParam->Nsite; i++)	{
				int nstate = Nstate;
				if (mParam->ModeFastCompute)	{
					nstate = mParam->ZipSize[i];
				}	
				DCMtbl[j][mode][ratemode][i] = new double[nstate];
			}
		}
	}
	DCMFlag[j] = 1;
}

void PhyloBayes::DeleteNodeDCMtbl(int j)	{

	if (DCMFlag[j])	{
		int nmode = 1;
		if (mParam->ActivateSumOverModes)	{
			nmode = Nmode;
		}
		int nratemode = 1;
		if (mParam->ActivateSumOverRateModes)	{
			nratemode = NRateMode;
		}
		for (int mode=0; mode<nmode; mode++)	{
			for (int ratemode=0; ratemode<nratemode; ratemode++)	{
				for (int i=0; i<mParam->Nsite; i++)	{
					delete[] DCMtbl[j][mode][ratemode][i];
				}
				delete[] DCMtbl[j][mode][ratemode];
			}
			delete[] DCMtbl[j][mode];
		}
		delete[] DCMtbl[j];
		DCMtbl[j] = 0;
		DCMFlag[j] = 0;
	}
}

void PhyloBayes::CreateDCM(double p)	{

	int size = 5;
	do	{
		for (int j=0; j<mParam->Nnode; j++)	{
			DCMCenter[j] = 0;
			DCMFlag[j] = 0;
		}
		DCMCenterSize = 1;
		DCMCenter[root->label] = 1;
		PropagateDCM(root->left->label,p);
		PropagateDCM(root->right->label,p);
		size = DCMCenterSize;
		if (size == 5){
			DeleteDCM();
		}
	} while (size == 5);
	if (mParam->DCM)	{
		DCM = 2;
		UpdateLogProbs();
		DCM = 1;
	}
}

void PhyloBayes::PropagateDCM(int j, double p)	{
	
	if (! tree[j].isLeaf())	{
		if ((tree[j].up->isRoot()) || (Random::Uniform() < p))	{
			PropagateDCM(tree[j].left->label,p);
			PropagateDCM(tree[j].right->label,p);
		}
		else	{
			if (mParam->DCM)	{
				CreateNodeDCMtbl(j);
			}
		}
	}
	DCMCenter[j] = 1;
	DCMCenterSize++;
}

void PhyloBayes::SwitchOffDCM(int j, int* dcm, double& length)	{
	if (dcm[j])	{
		dcm[j] = 0;
		length -= tree[j].branchLength;
	}
	if (! tree[j].isLeaf())	{
		SwitchOffDCM(tree[j].left->label, dcm, length);
		SwitchOffDCM(tree[j].right->label, dcm, length);
	}
}

void PhyloBayes::SwitchOnDCM(int j, int* dcm)	{
	dcm[j] = 1;
	if (! tree[j].isLeaf())	{
		SwitchOnDCM(tree[j].left->label, dcm);
		SwitchOnDCM(tree[j].right->label, dcm);
	}
}

void PhyloBayes::DeleteDCM()	{

	for (int j=0; j<mParam->Nnode; j++)	{
		DeleteNodeDCMtbl(j);
		DCMCenter[j] = 0;
	}
	DCMCenterSize = 0;
	DCM = 0;
}


void PhyloBayes::ActivateDiscreteGamma(int ncat)	{

	mParam->NRateModeMax = ncat;
	DeleteRateMode();
	CreateRateMode();

	if (mRateModeSiteLogSampling)	{
		cerr << "error in activate discrete gamma: array of rate-specific likelihoods already exists\n";
		exit(1);
	}

	// cerr << "number of discrete categories: " << mParam->NRateModeMax << '\n';
	mRateModeSiteLogSampling = new double*[mParam->NRateModeMax];
	for (int i=0; i<mParam->NRateModeMax; i++)	{
		mRateModeSiteLogSampling[i] = new double[mParam->Nsite];
	}
	NRateMode = mParam->NRateModeMax;
	mParam->RASModelSwitch = 0;
}




// ***************************************************************************
// ***************************************************************************
//		 Destructing
// ***************************************************************************
// ***************************************************************************



// ---------------------------------------------------------------------------
//		 ~PhyloBayes()
// ---------------------------------------------------------------------------

PhyloBayes::~PhyloBayes () 	{

	if (mParam->ActivateClock)	{
		DeleteClock();
	}
	DeleteTree();
	DeleteGene();
	if (Partialpdown)	{
		DeletePartial();
	}
	if (mParam->NormalApprox == No)	{
		DeleteRef();
		DeleteNH();
	 	DeleteMutSelSub();
		DeleteRateMode();
		DeleteRate();
		if (mParam->HeteroMode == Covarion)	{
			DeleteHetero();
		}
		if (SubCreated)	{
			DeleteSub();
		}
		DeleteLogProbs();		
		DeleteBaseFields();
	 	DeleteMode();
	}
	DeleteStatePostProb();
}

// ---------------------------------------------------------------------------
//		 DeleteClock()
// ---------------------------------------------------------------------------

void PhyloBayes::DeleteClock()	{

	if (mParam->TimePrior == BD)	{
		delete[] mBDLogGG;
	}

	delete[] Rho;
	if (mParam->ContNchar)	{
		for (int k=0; k<mParam->ContNchar; k++)	{
			delete[] ContRho[k+1];
		}
		delete[] ContRho;
		delete[] ContRhoCenter;
		for (int k=0; k<mParam->ContNchar+1; k++)	{
			delete[] ContCov[k];
		}
		delete[] ContCov;
		delete[] ContVar;
	}
	if (mParam->Ngene)	{
		for (int k=0; k<mParam->Ngene; k++)	{
			delete[] GeneRho[k];
		}
		delete[] GeneRho;
		delete[] GeneTheta;
		delete[] GeneSigma;
		delete[] GeneMu;
	}

	if (mParam->Nrho)	{
		delete[] Drho;
		for (int j=0; j<mParam->Nnode; j++)	{
			delete[] UpperRhoLogL[j];
			delete[] LowerRhoLogL[j];
			delete[] RhoLogNorm[j];
		}
		delete[] UpperRhoLogL;
		delete[] LowerRhoLogL;
		delete[] RhoLogNorm;
	}
	if (NDir)	{
		// delete dirichlet prior structures
		for (int i=0; i<NDir; i++)	{
			delete[] DirMap[i];
		}
		delete[] DirMap;
		delete[] DirSize;
		delete[] DirInvMap;
	}
}

// ---------------------------------------------------------------------------
//		 DeleteTree()
// ---------------------------------------------------------------------------

void PhyloBayes::DeleteTree()	{
	delete[] tree;
	delete[] BL;
	delete[] RigidBL;
	if (mParam->MBL)	{
		for (int i=0; i<mParam->NMixBLMax; i++)	{
			delete[] MixBL[i];
		}
		delete[] MixBL;
		delete[] MixBLSiteNumber;
		delete[] MixBLMode;
	}
}

// ---------------------------------------------------------------------------
//		 DeleteRate()
// ---------------------------------------------------------------------------

void PhyloBayes::DeleteRate()	{

	delete[] rate;
	delete[] unirate;
}
	
// ---------------------------------------------------------------------------
//		 DeleteHetero()
// ---------------------------------------------------------------------------

void PhyloBayes::DeleteHetero()	{

}

// ---------------------------------------------------------------------------
//		 DeleteGene()
// ---------------------------------------------------------------------------

void PhyloBayes::DeleteGene()	{

	if (mParam->Ngene)	{
	for (int i=0; i<mParam->Ngene; i++)	{
		delete[] GeneBL[i];
		delete[] GeneRigidBL[i];
		delete[] GeneBranch[i];
	}
	delete[] GeneBL;
	delete[] GeneRigidBL;
	delete[] GeneBranch;
	if (mParam->NormalApprox == No)	{
		delete[] GeneGamma;
		delete[] GeneRate;
		if (GeneStationary)	{
			for (int i=0; i<mParam->Ngene; i++)	{
				delete[] GeneStationary[i];
			}
			delete[] GeneStationary;
		}
		delete[] GeneRateFactor;
		if (mGeneMatrixArray)	{

			for (int i=0; i<mParam->Ngene; i++)	{
				delete mGeneMatrixArray[i];
			}
			delete[] mGeneMatrixArray;
		}
		if (GeneZipStationary)	{

			for (int i=0; i<mParam->Nsite; i++)		{
				delete[] GeneZipStationary[i];
			}
			delete[] GeneZipStationary;
		}
	}
	}
}


// ---------------------------------------------------------------------------
//		 DeleteMode()
// ---------------------------------------------------------------------------

void PhyloBayes::DeleteMode()	{

	if (mParam->ModeFastCompute && (! mParam->ModePoisson))	{
		if (mParam->ZipGTR == 4)	{
		}
		else if (mParam->ZipGTR >= 2)	{
			delete[] SiteZipSize;
			
			for (int i=0; i<mParam->NmodeMax; i++)	{
				delete[] SiteAAZip[i];
				delete[] SiteZipAA[i];
			}
			delete[] SiteAAZip;
			delete[] SiteZipAA;
			
			for (int i=0; i<mParam->NmodeMax; i++)	{
				delete[] SiteAASupport[i];
			}
			delete[] SiteAASupport;
			
			for (int i=0; i<mParam->NmodeMax; i++)	{
				delete[] ZipStationary[i];
				delete[] SiteZipRR[i];
			}
			delete[] ZipStationary;
			delete[] SiteZipRR;

		}
		else	{
			delete[] ModeZipSize;
			
			for (int i=0; i<mParam->NmodeMax; i++)	{
				delete[] ModeAAZip[i];
				delete[] ModeZipAA[i];
			}
			delete[] ModeAAZip;
			delete[] ModeZipAA;
			
			for (int i=0; i<mParam->NmodeMax; i++)	{
				delete[] ModeAASupport[i];
			}
			delete[] ModeAASupport;
			
			for (int i=0; i<mParam->NmodeMax; i++)	{
				delete[] ModeZipStationary[i];
				delete[] ModeZipRR[i];
			}
			delete[] ModeZipStationary;
			delete[] ModeZipRR;
		}
		delete[] ModeZipProb;
		for (int i=0; i<mParam->Ntaxa; i++)	{
			delete[] ZipGTRData[i];
		}
		delete[] ZipGTRData;
		for (int i=0; i<mParam->NZipModeMax; i++)	{
			delete[] ZipOcc[i];
		}
		delete[] ZipOcc;
		delete[] GlobalZipOcc;
	}

	if (mMatrixArray)	{
		for (int i=0; i<mParam->NmodeMax; i++)	{
			delete mMatrixArray[i];
		}
		delete[] mMatrixArray;
	}
	if (mCovMatrixArray)	{
		for (int i=0; i<mParam->NmodeMax; i++)	{
			for (int k=0; k<Ncov; k++)	{
				delete mCovMatrixArray[i][k];
			}
			delete[] mCovMatrixArray[i];
		}
		delete[] mCovMatrixArray;
	}
	if (mCovZipMatrixArray)	{
		for (int i=0; i<mParam->Nsite; i++)	{
			for (int k=0; k<Ncov; k++)	{
				delete mCovZipMatrixArray[i][k];
			}
			delete[] mCovZipMatrixArray[i];
		}
		delete[] mCovZipMatrixArray;
	}
	if (mZipMatrixArray)	{
		for (int i=0; i<mParam->Nsite; i++)	{
			delete mZipMatrixArray[i];
		}
		delete[] mZipMatrixArray;
	}

	if (Stationary)	{
		for (int i=0; i<mParam->NmodeMax ; i++)		{
			delete[] Stationary[i];
		}
		delete[] Stationary;
	}

	delete[] ModeRateFactor;
	if (ZipStationary)	{
		for (int i=0; i<mParam->Nsite; i++)		{
			delete[] ZipStationary[i];
		}
		delete[] ZipStationary;
	}

	delete[] Mode;
	delete[] ModeWeight;

	delete[] SiteNumber;
	if (RR)	{
		for (int i=0; i<mParam->Ngene; i++)	{
			delete[] RR[i];
		}
		delete[] RR;
	}

	delete[] ModeRR;
	delete[] ModeStatCenter;

	delete[] ModeGrandTotal;
	for (int i=0; i<mParam->NmodeMax; i++)	{
		delete[] ModeTotal[i];
		delete[] ModeStatBeta[i];
		if (mParam->Qmode)	{
			delete[] ModeRelRateBeta[i];
			delete[] ModeRelRateNsub[i];
		}
	}
	delete[] ModeTotal;
	delete[] ModeStatBeta;
	delete[] ModeRelRateBeta;
	delete[] ModeRelRateNsub;
}

void PhyloBayes::DeleteMutSelSub()	{

	if ((mParam->ZipGTR == 4) || (mParam->MutMode))	{
		for (int i=0; i<mParam->NmodeMax; i++)	{
			for (int j=0; j<Nstate; j++)	{
				delete[] ModeGWDoubleNsub[i][j];
			}
			delete[] ModeGWDoubleNsub[i];
			delete[] ModeGWNsub[i];
		}
		delete[] ModeGWNsub;
		delete[] ModeGWDoubleNsub;

		for (int i=0; i<mParam->NmodeMax; i++)	{
			delete[] ModeGWRootNsub[i];
		}
		delete[] ModeGWRootNsub;
		for (int i=0; i<mParam->NmodeMax; i++)	{
			for (int j=0; j<Nstate; j++)	{
				delete[] ModeGWStatBeta[i][j];
			}
			delete[] ModeGWStatBeta[i];
		}
		delete[] ModeGWStatBeta;
		for (int i=0; i<mParam->NmodeMax; i++)	{
			for (int j=0; j<Nstate; j++)	{
				delete[] ModeNucStatBeta[i][j];
			}
			delete[] ModeNucStatBeta[i];
		}
		delete[] ModeNucStatBeta;

		for (int i=0; i<mParam->Nsite; i++)	{
			for (int j=0; j<Nstate; j++)	{
				delete[] GWDoubleNsub[i][j];
			}
			delete[] GWDoubleNsub[i];
			delete[] GWNsub[i];
		}
		delete[] GWDoubleNsub;
		delete[] GWNsub;

		for (int i=0; i<mParam->Nsite; i++)	{
			delete[] GWRootNsub[i];
		}
		delete[] GWRootNsub;
		for (int i=0; i<mParam->Nsite; i++)	{
			for (int j=0; j<Nstate; j++)	{
				delete[] GWStatBeta[i][j];
			}
			delete[] GWStatBeta[i];
		}
		delete[] GWStatBeta;
		for (int i=0; i<mParam->Nsite; i++)	{
			for (int j=0; j<Nstate; j++)	{
				delete[] NucStatBeta[i][j];
			}
			delete[] NucStatBeta[i];
		}
		delete[] NucStatBeta;
	}
}

// ---------------------------------------------------------------------------
//		 DeleteRateMode()
// ---------------------------------------------------------------------------

void PhyloBayes::DeleteRateMode()	{

	delete[] ModeRate;
	delete[] RateMode;
	delete[] RateSiteNumber;
	delete[] RateModeWeight;
	delete[] RateModeTotal;
	delete[] SiteEffLength;
	delete[] RateModeEffLength;
}

// ---------------------------------------------------------------------------
//		 DeleteLogProbs()
// ---------------------------------------------------------------------------

void PhyloBayes::DeleteLogProbs()	{

	if (mParam->NormalApprox)	{
		delete[] mFlexGeneLogSampling;
		delete[] mRigidGeneLogSampling;
		delete[] mGeneLogSampling;
	}
	else	{
		if (tbl)	{
			for (int i=0; i<mParam->Nnode; i++)	{
				delete[] tbl[i];
			}
			delete[] tbl;
		}
		delete[] aux;
		delete[] taux;
		delete[] taux2;
		delete[] mSiteLogSampling;
		delete[] mSiteLogSampling0;

		if (mParam->ActivateSumOverModes)	{
			for (int i=0; i<mParam->NmodeMax; i++)	{
				delete[] mModeSiteLogSampling[i];
			}
			delete[] mModeSiteLogSampling;
		}
		if (mParam->ActivateSumOverRateModes)	{
			for (int i=0; i<mParam->NRateModeMax; i++)	{
				delete[] mRateModeSiteLogSampling[i];
			}
			delete[] mRateModeSiteLogSampling;
		}
		if (mParam->ZipSub == 2)	{
			for (int i=0; i<mParam->Nnode; i++)	{
				delete[] MissingFlag[i];
			}
			delete[] MissingFlag;
			MissingFlag = 0;
		}
	}
}

	
// ---------------------------------------------------------------------------
//		 DeleteRef()
// ---------------------------------------------------------------------------

void PhyloBayes::DeleteRef()	{

	if (RefZipStationary)	{

		for (int i=0; i<mParam->Nsite; i++)		{
			delete[] RefZipStationary[i];
		}
		delete[] RefZipStationary;
	}

	delete[] RefRR;
	delete[] RefStationary;
	delete mSubMatrix;
}

	
// ---------------------------------------------------------------------------
//		 DeleteSub()
// ---------------------------------------------------------------------------

void PhyloBayes::DeleteSub()	{
	
	for (int i=0; i<mParam->Nnode; i++)	{
		delete[] State[i];
		delete[] BranchSiteTotalSub[i];
	}
	delete[] State;
	delete[] BranchSiteTotalSub;

	if (mParam->ModeFastCompute)	{
		for (int i=0; i<mParam->Nnode; i++)	{
			delete[] TrueState[i];
		}
		delete[] TrueState;
	}

	for (int i=0; i<mParam->Nsite; i++)	{
		delete[] Nsub[i];
	}
	
	if (mParam->HeteroMode == Covarion)	{
		delete[] NSiteOnSub;
		delete[] NSiteOffSub;
		delete[] SiteXiEffLength;
	}

	delete[] TotalSub;
	delete[] BranchTotalSub;
	delete[] BranchEffLength;
	delete[] GeneTotalSub;
	delete[] GeneRateEffLength;
	for (int gene=0; gene<mParam->Ngene; gene++)	{
		delete[] GeneNsub[gene];
		delete[] GeneBranchTotalSub[gene];
		delete[] GeneBranchEffLength[gene];
	
		delete[] GeneStatBeta[gene];
	}
	delete[] GeneNsub;
	delete[] GeneBranchTotalSub;
	delete[] GeneBranchEffLength;
	delete[] GeneStatBeta;

	if (mParam->MBL)	{
		for (int i=0; i<mParam->NMixBLMax; i++)	{
			delete[] MixBLBranchEffLength[i];
			delete[] MixBLBranchTotalSub[i];
		}
		delete[] MixBLBranchEffLength;
		delete[] MixBLBranchTotalSub;
	}
	
	delete[] RefNsub;

	delete[] StatBeta;
	delete[] RelRateNsub;
	delete[] RelRateBeta;

	for (int i=0; i<mParam->Nsite; i++)	{	
		delete[] SiteStatBeta[i];
		delete[] SiteRelRateBeta[i];
		delete[] SiteRelRateNsub[i];
	}
	delete[] SiteStatBeta;
	delete[] SiteRelRateBeta;
	delete[] SiteRelRateNsub;
	if (mParam->Qmode)	{
		for (int i=0; i<mParam->Nsite; i++)	{	
			for (int j=0; j<Nstate; j++)	{
				delete[] SiteRRFactor[i][j];
			}
			delete[] SiteRRFactor[i];
		}
		delete[] SiteRRFactor;
	}
}

// ---------------------------------------------------------------------------
//		 DeleteBaseFields()
// ---------------------------------------------------------------------------

void PhyloBayes::DeleteBaseFields()	{

	delete[] BaseRate;
	delete[] BaseGeneRate;
	delete[] BaseRateFactor;
	delete[] BaseRefRateFactor;
	delete[] BaseStationary;
	delete[] BaseRefStationary;
	delete[] BaseMatrix;
	delete[] BaseRefMatrix;

	delete[] BaseBL;
	delete[] BaseBLMul;
	delete[] UniBLMul;

	if (EffLength)	{
		for (int j=0; j<mParam->Nnode; j++)	{
			delete[] EffLength[j];	
		}
		delete[] EffLength;
	}
}



// ***************************************************************************
// ***************************************************************************
//		 Cloning
// ***************************************************************************
// ***************************************************************************


// ---------------------------------------------------------------------------
//		 Clone(PhyloBayes* from)
// --------------------------------------------------------------------------

void PhyloBayes::Clone(PhyloBayes* from)	{

	// check whether not cloning from self
	if (this != from)	{

		CloneTree(from);	
		if (mParam->ActivateClock)	{
			CloneClock(from);
		}
		CloneGene(from);

		if (mParam->NormalApprox == No)	{
			CloneRate(from);
			CloneRateMode(from);	
			CloneNH(from);
			CloneMode(from);
			if (mParam->HeteroMode != Homo)	{
				CloneHetero(from);
			}
		}
		CloneLogProbs(from);
	}
}


// ---------------------------------------------------------------------------
//		 CloneTree(PhyloBayes* from)
// ---------------------------------------------------------------------------

void PhyloBayes::CloneTree(PhyloBayes* from)	{
	
	root = cloneTree(tree, from->tree, from->root);
	for (int j=0; j<mParam->Nnode; j++)	{
		BL[j] = from->BL[j];
		RigidBL[j] = from->RigidBL[j];
	}
	if (mParam->GeneBLMultiplier)	{
		for (int gene=0; gene<mParam->Ngene; gene++)	{
			for (int j=0; j<mParam->Nnode; j++)	{
				GeneBL[gene][j] = from->GeneBL[gene][j];
			}
		}
	}
	CloneMBL(from);
	MeanLength = from->MeanLength;
	VarLength = from->VarLength;
	LengthGamma = from->LengthGamma;
}

// ---------------------------------------------------------------------------
//		 CloneMBL()
// ---------------------------------------------------------------------------

void PhyloBayes::CloneMBL(PhyloBayes* from)	{
	if (mParam->MBL)	{
		MixBLAlpha = from->MixBLAlpha;
		NMixBL = from->NMixBL;
		for (int i=0; i<NMixBL; i++)	{
			for (int j=0; j<mParam->Nnode; j++)	{
				MixBL[i][j] = from->MixBL[i][j];
			}
			MixBLSiteNumber[i] = from->MixBLSiteNumber[i];
		}
		for (int i=0; i<mParam->Nsite; i++)	{
			MixBLMode[i] = from->MixBLMode[i];
		}
	}
}

// ---------------------------------------------------------------------------
//		 cloneTree()
// ---------------------------------------------------------------------------

AmphiNode*
PhyloBayes::cloneTree(AmphiNode* inTree, AmphiNode* fromTree, AmphiNode* fromRoot)	{

	AmphiNode* p = inTree;
	AmphiNode* fromp = fromTree;

	for (int i=0; i< mParam->Nnode; i++)	{

		p->label = i;
		p->branchLength = fromp->branchLength;
		p->up = (fromp->up) ? &inTree[fromp->up->label] : 0;
		p->left = (fromp->left) ? &inTree[fromp->left->label] : 0;
		p->right = (fromp->right) ? &inTree[fromp->right->label] : 0;
		p++;
		fromp++;
	}

	return &inTree[fromRoot->label];
}


// ---------------------------------------------------------------------------
//		 CloneRate
// ---------------------------------------------------------------------------

void
PhyloBayes::CloneRate(PhyloBayes* from)	{

	gamma = from->gamma;
	pconst = from->pconst;
	for (int i=0; i< mParam->Nsite; i++)	{
		rate[i] = from->rate[i];
		unirate[i] = from->unirate[i];
		ConstantStatus[i] = from->ConstantStatus[i];
	}
}

// ---------------------------------------------------------------------------
//		 CloneMode
// ---------------------------------------------------------------------------

void
PhyloBayes::CloneMode(PhyloBayes* from)	{

	CloneRef(from);
	GWf = from->GWf;
	RRAlpha = from->RRAlpha;
	if (mParam->MutMode == 2)	{
		Kappa = from->Kappa;
		for (int i=0; i<Nnuc; i++)	{
			RefACGT[i] = from->RefACGT[i];
		}
	}
	else if (mParam->MutMode == 3)	{
		for (int i=0; i<Nnuc; i++)	{
			RefACGT[i] = from->RefACGT[i];
		}
		for (int i=0; i<Nnuc * (Nnuc-1) / 2; i++)	{
			RefRRACGT[i] = from->RefRRACGT[i];
		}
	}
	Nmode = from->Nmode;
	alpha = from->alpha;
	for (int i=0; i< mParam->Nsite; i++)	{
		Mode[i]= from->Mode[i];
	}
	for (int j=0; j<Nrr; j++)	{
		ModeRR[j] = from->ModeRR[j];
	}
	for (int j=0; j<Nstate; j++)	{
		ModeStatCenter[j] = from->ModeStatCenter[j];
	}
	if( (mParam->ModeFastCompute) && (! mParam->ModePoisson))	{
		if (mParam->ZipGTR == 4)	{
		}
		else	{
			NZipMode = from->NZipMode;
			ZipAlpha = from->ZipAlpha;
			for (int j=0; j<Nstate; j++)	{
				ModeZipProb[j] = from->ModeZipProb[j];
				ZipA0[j] = from->ZipA0[j];
				ZipA1[j] = from->ZipA1[j];
			}
			for (int i=0; i<mParam->Nsite; i++)	{
				ZipMode[i] = from->ZipMode[i];
			}
			ZipRate = from->ZipRate;
			
			if (mParam->ZipGTRDP)	{
				for (int i=0; i<NZipMode; i++)	{
					ModeZipSize[i] = from->ModeZipSize[i];
					for (int k=0; k<Nstate; k++)	{
						ModeAASupport[i][k] = from->ModeAASupport[i][k];
					}
				}
			}
		}
	}
	ModeStatAlpha = from->ModeStatAlpha;

	for (int i=0; i<Nmode; i++)	{
		CloneMode(from,i);
	}

}

void PhyloBayes::WriteQMM(ostream& os, int mode) {
    os << ((double) SiteNumber[mode]) / mParam->Nsite << '\t';
    for (int j=0; j<Nrr; j++)   {
        os << RR[mode][j] << '\t';
    }
    for (int j=0; j<Nstate; j++)    {
        os << Stationary[mode][j] << '\t';
    }
    os << '\n';
}

void
PhyloBayes::CloneMode(PhyloBayes* from, int mode)	{

	for (int j=0; j<Nstate; j++)	{
		Stationary[mode][j] = from->Stationary[mode][j];
	}
	ModeWeight[mode] = from->ModeWeight[mode];
	SiteNumber[mode] = from->SiteNumber[mode];

	if (! mParam->ModeFastCompute)	{
		if (mParam->Qmode)	{
			for (int j=0; j<Nrr; j++)	{
				RR[mode][j] = from->RR[mode][j];
			}
		}
		if (mParam->HeteroMode == Covarion)	{
			// Copy all the covarion matrices
			for (int k=0; k<Ncov; k++)	{
				mCovMatrixArray[mode][k]->Copy(from->mCovMatrixArray[mode][k]);
			}
		}
		else	{
			if (mParam->NH)	{
				for (int n=0; n<NHNcat; n++)	{
					CloneNH(from,mode,n);
				}
			}
			else	{
				mMatrixArray[mode]->Copy(from->mMatrixArray[mode]);
			}
		}
	}
	else	{
		if (mParam->ModePoisson)	{
			ModeRateFactor[mode] = from->ModeRateFactor[mode];
			if (mParam->NH)	{
				for (int n=0; n<NHNcat; n++)	{
					CloneNH(from,mode,n);
				}
			}
			else	{
				for (int i=0; i<mParam->Nsite; i++)	{
					if (Mode[i] == mode)	{
						for (int j=0; j<mParam->ZipSize[i]; j++)	{
							ZipStationary[i][j] = from->ZipStationary[i][j];
						}
						if (mParam->HeteroMode == Covarion)	{
							// Copy all the site-specific covarion matrices
							for (int k=0; k<Ncov; k++)	{
								mCovZipMatrixArray[i][k]->Copy(from->mCovZipMatrixArray[i][k]);
							}
						}
					}
				}
			}
		}
		else	{
			if (mParam->Qmode)	{
				for (int j=0; j<Nrr; j++)	{
					RR[mode][j] = from->RR[mode][j];
				}
			}
			if (mParam->ZipGTR == 4)	{
				for (int i=0; i<mParam->Nsite; i++)	{
					if (Mode[i] == mode)	{
						/*
						for (int j=0; j<mParam->ZipSize[i]; j++)	{
							ZipStationary[i][j] = from->ZipStationary[i][j];
						}
						int nrr = mParam->ZipSize[i] * (mParam->ZipSize[i] - 1) / 2;
						for (int j=0; j<nrr; j++)	{
							SiteZipRR[i][j] = from->SiteZipRR[i][j];
						}
						*/
						mZipMatrixArray[i]->Copy(from->mZipMatrixArray[i]);
					}
				}
			}
			else if (mParam->ZipGTR >= 2)	{
				if (mParam->NH)	{
					for (int n=0; n<NHNcat; n++)	{
						CloneNH(from,mode,n);
					}
				}
				else	{
					for (int site=0; site<mParam->Nsite; site++)	{
						if (Mode[site] == mode)	{
							for (int i=0; i<SiteZipSize[site]; i++)	{
								ZipStationary[site][i] = from->ZipStationary[site][i];
							}
						}
					}
				}
				for (int site=0; site<mParam->Nsite; site++)	{
					if (Mode[site] == mode)	{
						SiteZipSize[site] = from->SiteZipSize[site];
						for (int i=0; i<Nstate; i++)	{
							SiteAASupport[site][i] = from->SiteAASupport[site][i];
						}
						for (int i=0; i<SiteZipSize[site]; i++)	{
							SiteAAZip[site][i] = from->SiteAAZip[site][i];
							SiteZipAA[site][i] = from->SiteZipAA[site][i];
						}
						int nrr = SiteZipSize[site] * (SiteZipSize[site] - 1) / 2;
						for (int i=0; i<nrr; i++)	{
							SiteZipRR[site][i] = from->SiteZipRR[site][i];
						}
					}
				}
			}
			else	{
				if (mParam->NH)	{
					for (int n=0; n<NHNcat; n++)	{
						CloneNH(from,mode,n);
					}
				}
				else	{
					for (int i=0; i<ModeZipSize[mode]; i++)	{
						ModeZipStationary[mode][i] = from->ModeZipStationary[mode][i];
					}
				}
				ModeZipSize[mode] = from->ModeZipSize[mode];
				for (int i=0; i<Nstate; i++)	{
					ModeAASupport[mode][i] = from->ModeAASupport[mode][i];
				}
				for (int i=0; i<ModeZipSize[mode]; i++)	{
					ModeAAZip[mode][i] = from->ModeAAZip[mode][i];
					ModeZipAA[mode][i] = from->ModeZipAA[mode][i];
				}
				int nrr = ModeZipSize[mode] * (ModeZipSize[mode] - 1) / 2;
				for (int i=0; i<nrr; i++)	{
					ModeZipRR[mode][i] = from->ModeZipRR[mode][i];
				}
			}
		}
	}
}


// ---------------------------------------------------------------------------
//		 CloneHetero
// ---------------------------------------------------------------------------

void
PhyloBayes::CloneHetero(PhyloBayes* from)	{

	xi0 = from->xi0;
	invprob = from->invprob;
}


	
// ---------------------------------------------------------------------------
//		 CloneRateMode()
// ---------------------------------------------------------------------------


void PhyloBayes::CloneRateMode(PhyloBayes* from)	{

	// what about Pconst0 ?

	for (int i=0; i<mParam->Nsite; i++)	{
		RateMode[i] = from->RateMode[i];
		ConstantStatus[i] = from->ConstantStatus[i];
	}
	NRateMode = from->NRateMode;
	for (int i=0; i<NRateMode; i++)	{
		ModeRate[i] = from->ModeRate[i];
		RateSiteNumber[i] = from->RateSiteNumber[i];
		RateModeWeight[i] = from->RateModeWeight[i];
	}
	RateAlpha = from->RateAlpha;
	for (int i=0; i<mParam->Nsite; i++)	{
		unirate[i] = from->unirate[i];
	}
}


// ---------------------------------------------------------------------------
//		 CloneRef
// ---------------------------------------------------------------------------

void PhyloBayes::CloneRef(PhyloBayes* from)	{

	for (int k=0; k<Nstate; k++)	{
		RefStationary[k] = from->RefStationary[k];
	} 
	for (int j=0; j<Nrr; j++)	{
		RefRR[j] = from->RefRR[j];
	}
	RefRateFactor = from->RefRateFactor;
	if (mParam->RefFastCompute)	{
		for (int i=0; i<mParam->Nsite; i++)	{
			for (int j=0; j<mParam->ZipSize[i]; j++)	{
				RefZipStationary[i][j] = from->RefZipStationary[i][j];
			}
		}
	}
	else	{
		mSubMatrix->Copy(from->mSubMatrix);
	}
}

// ---------------------------------------------------------------------------
//		 CloneGene(PhyloBayes* from)
// ---------------------------------------------------------------------------

void PhyloBayes::CloneGene(PhyloBayes* from)	{

	for (int i=0; i<mParam->Ngene; i++)	{
		for (int j=0; j<mParam->Nnode; j++)	{
			GeneBL[i][j] = from->GeneBL[i][j];
			GeneRigidBL[i][j] = from->GeneRigidBL[i][j];
			GeneBranch[i][j] = from->GeneBranch[i][j];
		}
	}
	if (!mParam->NormalApprox)	{
		for (int i=0; i<mParam->Ngene; i++)	{
			if (mParam->GeneRateMode)	{
				GeneRate[i] = from->GeneRate[i];
			}
			if (mParam->GeneGammaMode)	{
				GeneGamma[i] = from->GeneGamma[i];
			}
			if (mParam->GeneStationaryMode)	{
				CloneGeneMatrix(from, i);
			}
		}
	}
}


// ---------------------------------------------------------------------------
//		 CloneGeneMatrix(PhyloBayes* from, int gene)
// ---------------------------------------------------------------------------


void PhyloBayes::CloneGeneMatrix(PhyloBayes* from, int gene)	{

	if (! mParam->GeneStationaryMode)	{
		cerr << "error in PhyloBayes::CloneGeneMatrix()\n";
		exit(1);
	}	
	for (int j=0; j<Nstate; j++)	{
		GeneStationary[gene][j] = from->GeneStationary[gene][j];
	}
	if (! mParam->RefFastCompute)	{
		mGeneMatrixArray[gene]->Copy(from->mGeneMatrixArray[gene]);
	}
	else	{
		GeneRateFactor[gene] = from->GeneRateFactor[gene];
		for (int i=0; i<mParam->GeneSize[gene]; i++)	{
			int site = mParam->GeneFirstSite[gene] + i;
			for (int j=0; j<mParam->ZipSize[site]; j++)	{
				GeneZipStationary[site][j] = from->GeneZipStationary[site][j];
			}
		}
	}
}


// ---------------------------------------------------------------------------
//		 CloneState
// ---------------------------------------------------------------------------

void PhyloBayes::CloneState(PhyloBayes* from)	{

	for (int j=0; j<mParam->Nnode; j++)	{
		for (int i=0; i<mParam->Nsite; i++)	{
			State[j][i] = from->State[j][i];
		}
	}
}


// ---------------------------------------------------------------------------
//		 CloneSub
// ---------------------------------------------------------------------------

void
PhyloBayes::CloneSub(PhyloBayes* from)	{

	if (Nmode != from->Nmode)	{
		cerr << "error in CloneSub()\n";
		exit(1);
	}

	// copy nsub, ntotalsub, modetotal and modegrandtotal
	for (int i=0; i<mParam->Nsite; i++)	{
		TotalSub[i] = from->TotalSub[i];
		for (int k=0; k<Nstate; k++)	{
			Nsub[i][k] = from->Nsub[i][k];
		}
	}
	for (int mode = 0; mode < Nmode; mode++)	{
		ModeGrandTotal[mode] = from->ModeGrandTotal[mode];
		for (int k=0; k<Nstate; k++)	{
			ModeTotal[mode][k] = from->ModeTotal[mode][k];
		}
	}
}

// ---------------------------------------------------------------------------
//		 CloneLogProbs
// ---------------------------------------------------------------------------

void
PhyloBayes::CloneLogProbs(PhyloBayes* from)	{

	mLogPrior = from->mLogPrior;
	mLogSampling = from->mLogSampling;
	mLogPosterior = from->mLogPosterior;

	if (mParam->NormalApprox)	{
		for (int k=0; k<mParam->Ngene; k++)	{
			mFlexGeneLogSampling[k] = from->mFlexGeneLogSampling[k];
			mRigidGeneLogSampling[k] = from->mRigidGeneLogSampling[k];
			mGeneLogSampling[k] = from->mGeneLogSampling[k];
		}
		mSeparateNormalLogSampling = from->mSeparateNormalLogSampling;
		mConcatenateNormalLogSampling = from->mConcatenateNormalLogSampling;
	}
	else	{
		for (int i=0; i<mParam->Nsite; i++)	{
			mSiteLogSampling[i] = from->mSiteLogSampling[i];
		}

		if (mParam->SumOverModes)	{
			CloneSumMode(from);
		}
		if (mParam->SumOverRateModes)	{
			CloneSumRateMode(from);
		}
	}
}

void
PhyloBayes::CloneLogProbs(PhyloBayes* from, int i)	{

	mLogPrior = from->mLogPrior;
	mLogSampling = from->mLogSampling;
	mLogPosterior = from->mLogPosterior;
	mSiteLogSampling[i] = from->mSiteLogSampling[i];
}

void PhyloBayes::CloneSumMode(PhyloBayes* from)	{

	for (int i=0; i<mParam->Nsite; i++)	{
		mSiteLogSampling0[i] = from->mSiteLogSampling0[i];
		for (int mode=0; mode < Nmode; mode++)	{
			mModeSiteLogSampling[mode][i] = from->mModeSiteLogSampling[mode][i];
		}
	}
}

void PhyloBayes::CloneSumRateMode(PhyloBayes* from)	{

	for (int i=0; i<mParam->Nsite; i++)	{
		mSiteLogSampling0[i] = from->mSiteLogSampling0[i];
		for (int mode=0; mode < NRateMode; mode++)	{
			mRateModeSiteLogSampling[mode][i] = from->mRateModeSiteLogSampling[mode][i];
		}
	}
}

// ---------------------------------------------------------------------------
//		 CloneClock
// ---------------------------------------------------------------------------

void
PhyloBayes::CloneClock(PhyloBayes* from)	{

	Chi = from->Chi;
	Chi2 = from->Chi2;
	MulTheta = from->VarTheta;
	MulSigma = from->VarSigma;
	MeanTheta = from->VarTheta;
	MeanSigma = from->VarSigma;
	VarTheta = from->VarTheta;
	VarSigma = from->VarSigma;
	Theta = from->Theta;
	Sigma = from->Sigma;
	Scale = from->Scale;
	VarScale = from->VarScale;

	if (mParam->ContNchar)	{
		for (int k=0; k<mParam->ContNchar+1; k++)	{
			for (int i=0; i<mParam->Nnode; i++)	{
				ContRho[k][i] = from->ContRho[k][i];
			}
			ContRhoCenter[k] = from->ContRhoCenter[k];
		}
	}
	else	{
		for (int i=0; i<mParam->Nnode; i++)	{
			Rho[i] = from->Rho[i];
		}
	}

	Mu = from->Mu;

	for (int k=0; k<mParam->Ngene; k++)	{
		for (int i=0; i<mParam->Nnode; i++)	{
			GeneRho[k][i] = from->GeneRho[k][i];
		}
		GeneTheta[k] = from->GeneTheta[k];
		GeneSigma[k] = from->GeneSigma[k];
		GeneMu[k] = from->GeneMu[k];
	}
}		

int PhyloBayes::blfree(int j)	{

	if (j == root->label)	{
		return 0;
	}
	if ((! mParam->ActivateClock) && (! mParam->NH) && (! mParam->midpointrooting))	{
		if (mParam->FixTopo)	{
			if (! tree[j].branchLength)	{
				if ((j == root->left->label) || (j == root->right->label))	{
					return 0;
				}
				else	{
					cerr << "error : null branch length for non root nor root son\n";
					exit(1);
				}
			}
			return 1;
		}
		else	{
			if (j == root->left->label)	{
				/*
				if (tree[j].branchLength)	{
					cerr << "error : non null bl for left branch\n";
				}
				*/
				return 0;
			}
			if (j == root->right->label)	{
				/*
				if (! tree[j].branchLength)	{
					cerr << "error : null bl for right branch\n";
				}
				*/
				return 1;
			}
		}
	}
	return 1;
}

int PhyloBayes::nodefree(int j)	{

	if ((! mParam->ActivateClock) && (! mParam->NH) && (! mParam->midpointrooting))	{
		if ((j == root->label) && (!BL[root->left->label]))	{
			return 0;
		}
		return 1;
	}
	return 1;
}


void PhyloBayes::checkbl()	{

	if (BL[root->label])	{
		cerr << "error : non zera branch length for root : " << BL[root->label] << '\n';
		exit(1);
	}
	if ((! mParam->ActivateClock) && (! mParam->NH))	{
		if (BL[root->left->label])	{
			cerr << "error : non zera branch length for root left : " << BL[root->left->label] << '\n';
			exit(1);
		}
	}
}

double PhyloBayes::LengthNormFactor()	{

	if (mParam->ZipGTR)	{
		return 1;
	}
	if (mParam->NormalApprox)	{
		return 1;
	}
	double norm = 0;
	/*
	UpdateBaseFields();
	if (mParam->ModePoisson && (!mParam->HeteroMode == Covarion))	{
		for (int i=0; i<mParam->Nsite; i++)	{
			double temp = 0;
			for (int j=0; j<Nstate-1 ; j++)	{
				for (int k=j+1; k<Nstate; k++)	{
					temp += BaseStationary[i][j] * BaseStationary[i][k];
				}
			}
			norm += 2 * temp;
		}
		norm /= mParam->Nsite;
	}
	else	{
		for (int i=0; i<mParam->Nsite; i++)	{
			double temp = 0;
			SubMatrix* mat = BaseMatrix[i];
			double* matstat = mat->GetStationaries();
			int Nstate = mat->Nstate;
			for (int j=0; j<Nstate-1 ; j++)	{
				for (int k=j+1; k<Nstate; k++)	{
					temp += matstat[j] * (*mat)(j,k);
				}
			}
			norm += 2 * temp;
		}
		norm /= mParam->Nsite;
	}
	*/

	UpdateSiteNumber();
	if (mParam->ModePoisson)	{
		for (int mode=0; mode<Nmode; mode++)	{
			double temp = 0;
			for (int j=0; j<Nstate-1 ; j++)	{
				for (int k=j+1; k<Nstate; k++)	{
					temp += Stationary[mode][j] * Stationary[mode][k];
				}
			}
			norm += (SiteNumber[mode] + 1) * 2 * temp;
		}
		norm /= (mParam->Nsite + Nmode);
	}
	/*
	else if (mParam->ZipGTR == 4)	{
		for (int i=0; i<mParam->Nsite; i++)	{
			double temp = 0;
			int nstate = mParam->ZipSize[i];
			for (int j=0; j<nstate-1 ; j++)	{
				for (int k=j+1; k<nstate; k++)	{
					temp += ZipStationary[i][j] * (*mZipMatrixArray[i])(j,k);
				}
			}
			norm += 2 * temp;
		}
		norm /= mParam->Nsite;
	}
	*/
	else	{
		for (int mode=0; mode<Nmode; mode++)	{
			double* rr = 0;
			if (mParam->Qmode)	{
				rr = RR[mode];
			}
			else	{
				rr = ModeRR;
			}
			double temp = 0;
			if ((mParam->ZipGTR == 4) || (mParam->HeteroMode == Covarion))	{
				for (int j=0; j<Nstate-1 ; j++)	{
					for (int k=j+1; k<Nstate; k++)	{
						temp += Stationary[mode][j] * Stationary[mode][k] * rr[(2*Nstate - j - 1) * j / 2 + k - j - 1];
					}
				}
			}
			else	{
				SubMatrix* mat = mMatrixArray[mode];
				double* matstat = mMatrixArray[mode]->GetStationaries();
				for (int j=0; j<Nstate-1 ; j++)	{
					for (int k=j+1; k<Nstate; k++)	{
						temp += matstat[j] * (*mat)(j,k);
					}
				}
			}
			norm += (SiteNumber[mode] + 1) * 2 * temp;
		}
		norm /= (mParam->Nsite + Nmode);
	}
	norm *= GetMeanRate() * (1 - pconst);
	return norm;
}

void PhyloBayes::NormaliseLengths()	{
	double norm = LengthNormFactor();
	for (int j=0; j<mParam->Nnode; j++)	{
		BL[j] *= norm;
	}
}

void PhyloBayes::DenormaliseLengths()	{
	double norm = LengthNormFactor();
	for (int j=0; j<mParam->Nnode; j++)	{
		BL[j] /= norm;
	}
}

void PhyloBayes::RootBipartition(Bipartition& bp)	{
	FillBipartition(bp,root->left);	
	bp.Modulo();
}

void PhyloBayes::FillBipartition(Bipartition& bp, Node* node)	{

	if (node->isLeaf())	{
		bp.mArray[node->label] = 1;
	}
	else	{
		FillBipartition(bp,node->left);
		FillBipartition(bp,node->right);
	}
}

// ***************************************************************************
// ***************************************************************************
//		 Updating
// ***************************************************************************
// ***************************************************************************


// ---------------------------------------------------------------------------
//		 Update()
// ---------------------------------------------------------------------------


void
PhyloBayes::Update()	{

	DCM = 0;

	if (!mParam->ActivateClock)	{
		if (root->left->branchLength != 0)	{
			// cerr << "error: unrooted tree with uncorrect pattern of branch lengths\n";
			// root->left->branchLength = 0;
		}
	}
	
	if (mParam->NormalApprox == No)	{
		// Integral = 0;
		UpdateTree();
		UpdateRateMode();
		UpdateMode();
		UpdateMBL();
		UpdateRef();
		UpdateGene();
		UpdateBaseFields();
	}
	if (mParam->ActivateClock)	{
		SetTau();
		if (mParam->Nrho)	{
			UpdateDrho();
		}
		if (!NDir)	{
			RegisterDirichletTimePrior();
		}
	}
 	UpdateLogProbs();


}



// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//		 UpdateTree()
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------


void PhyloBayes::UpdateTree()	{

	UpdateTotalLength();
	if (mParam->ActivateClock)	{
		if (!NDir)	{
			RegisterDirichletTimePrior();
		}
	}
}

double PhyloBayes::UpdateTotalLength()	{
	TotalLength = 0;
	for (int i=0; i<mParam->Nnode; i++)	{
		TotalLength += BL[i];
	}
	if (mParam->MBL == 2)	{
		for (int i=0; i<mParam->Nnode; i++)	{
			if (blfree(i))	{
				UniBLMul[i] = TotalLength;
			}
			else	{
				UniBLMul[i] = 0;
			}
		}
	}
	return TotalLength;
}

// ---------------------------------------------------------------------------
//		 UpdateMBL()
// ---------------------------------------------------------------------------

void PhyloBayes::UpdateMBL()	{

	if (mParam->MBL)	{
		for (int i=0; i<NMixBL; i++)	{
			MixBLSiteNumber[i] = 0;
		}
		for (int i=0; i<mParam->Nsite; i++)	{
			MixBLSiteNumber[MixBLMode[i]]++;
		}
	}
}


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//		 UpdateMode()
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void PhyloBayes::UpdateSiteNumber()	{

	// restore inversion map between modes and sites
	for (int i=0; i<Nmode; i++)	{
		SiteNumber[i] = 0;
	}
	for (int i=0; i<mParam->Nsite; i++)	{
		SiteNumber[Mode[i]]++;
	}
}

void PhyloBayes::UpdateMode()	{
	
	UpdateSiteNumber();
	if (mParam->NH)	{
		if (mParam->NHPrior)	{
			ComputeNHStatDistorter();
			if (mParam->NHPrior == 2)	{
				ComputeNHCov();
			}
		}
		ComputeNHStat();
	}
	if ((mParam->ZipGTR != 4) && (mParam->ZipGTR >= 2) && (mParam->ZipGTRDP))	{
		for (int i=0; i<mParam->Nsite; i++)	{
			for (int k=0; k<Nstate; k++)	{
				SiteAASupport[i][k] = ModeAASupport[ZipMode[i]][k];
			}
		}
	}
	for (int mode=0; mode<Nmode; mode++)	{
		UpdateMode(mode);
	}
	if ((mParam->ZipGTR != 4) && (mParam->ModeFastCompute) && (!mParam->ModePoisson))	{
		UpdateZipOccupationNumber();
		UpdateZipGTRData();
	}
}

void PhyloBayes::UpdateMode(int mode)	{

	if (mParam->ModeFastCompute)	{
		if (mParam->ModePoisson)	{
			if (mParam->NH)	{
				for (int n=0; n<NHNcat; n++)	{
					UpdateNH(mode,n);
				}
			}
			else	{
				ComputeRateFactor(mode);
				for (int i=0; i<mParam->Nsite; i++)	{
					if (Mode[i] == mode)	{
						UpdateZip(i);
					}
				}
			}
		}
		else	{
			if (mParam->ZipGTR == 4)	{
				ComputeRateFactor(mode);
				for (int i=0; i<mParam->Nsite; i++)	{
					if (Mode[i] == mode)	{
						UpdateZip(i);
					}
				}
			}
			else if (mParam->ZipGTR >= 2)	{
				for (int i=0; i<mParam->Nsite; i++)	{
					if (Mode[i] == mode)	{	
						SiteZip(i);
						RecreateSiteMatrix(i);
						if (mParam->NH)	{
							for (int n=0; n<NHNcat; n++)	{
								mNHMatrixArray[i][n]->ComputeArray();
							}
						}
						else	{
							mZipMatrixArray[i]->ComputeArray();
						}
					}
				}
			}
			else	{
				ModeZip(mode);
				RecreateModeMatrix(mode);
				if (mParam->NH)	{
					for (int n=0; n<NHNcat; n++)	{
						mNHMatrixArray[mode][n]->ComputeArray();
					}
				}
				else	{
					mMatrixArray[mode]->ComputeArray();
				}
			}
		}
	}
	else	{
		if (mParam->HeteroMode == Covarion)	{
			if (mParam->ExternalCov)	{
				mCovMatrixArray[mode][0]->ComputeArray();
			}
			else	{
				for (int k=0; k<NRateMode; k++)	{
					mCovMatrixArray[mode][k]->ComputeArray();
				}
			}
		}
		else	{
			mMatrixArray[mode]->ComputeArray();
		}
	}	
}


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//		 UpdateRateMode()
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------


void PhyloBayes::UpdateRateMode()	{

	if (mParam->GammaNcat)	{
		DiscCateRate();
	}

	for (int i=0; i<NRateMode; i++)	{
		RateSiteNumber[i] = 0;
	}

	for (int i=0; i<mParam->Nsite; i++)	{
		// if (! ConstantStatus[i])	{
			RateSiteNumber[RateMode[i]]++;
		// }
	}
	MakeRates();
}

// ---------------------------------------------------------------------------
//		 MakeRates()
// ---------------------------------------------------------------------------


void PhyloBayes::MakeRates()	{
	for (int i=0; i<mParam->Nsite; i++)	{
		unirate[i] = ModeRate[RateMode[i]];
	}
	if ((mParam->HeteroMode == Covarion) && (! mParam->ExternalCov))	{
		UpdateMode();
	}
}

// ---------------------------------------------------------------------------
//		 UpdateBaseFields()
// ---------------------------------------------------------------------------

void PhyloBayes::UpdateZipGTRData()	{

	if (mParam->ZipGTR == 4)	{
		cerr << "error  in update zip gtr data\n";
		exit(1);
	}
	if (mParam->ZipGTR >= 2)	{
		for (int i=0; i<mParam->Nsite; i++)	{
			for (int j=0; j<mParam->Ntaxa; j++)	{
				if (Data[j][i] == unknown)	{
					ZipGTRData[j][i] = unknown;
				}
				else	{
					ZipGTRData[j][i] = SiteAAZip[i][Data[j][i]];
				}
			}
		}
	}
	else	{
		for (int i=0; i<mParam->Nsite; i++)	{
			int mode = Mode[i];
			for (int j=0; j<mParam->Ntaxa; j++)	{
				if (Data[j][i] == unknown)	{
					ZipGTRData[j][i] = unknown;
				}
				else	{
					ZipGTRData[j][i] = ModeAAZip[mode][Data[j][i]];
				}
			}
		}
	}
}

// ---------------------------------------------------------------------------
//		 UpdateBaseFields()
// ---------------------------------------------------------------------------

void
PhyloBayes::UpdateBaseFields()	{
	if (mParam->ModeFastCompute)	{
		if (mParam->ModePoisson || (mParam->ZipGTR == 4))	{
			BaseData = ZipData;
		}
		else	{
			UpdateZipGTRData();
			BaseData = ZipGTRData;
		}
	}
	else	{
		BaseData = Data;
	}
	for (int i=0; i<mParam->Nsite; i++)	{
		UpdateBaseFields(i);	
	}
}

void
PhyloBayes::UpdateBaseFields(int i)	{
	
	if (mParam->SeparateModelSwitch == 0)	{
		if (mParam->BLModelSwitch == 0)	{
			BaseBL[i] = RigidBL;
		}
		else	{
			BaseBL[i] = BL;
		}
		// BaseBL[i] = BL;
		if (mParam->MBL == 2)	{
			BaseBL[i] = UniBLMul;
			BaseBLMul[i] = MixBL[MixBLMode[i]];
		}
		else if (mParam->MBL == 1)	{
			BaseBLMul[i] = MixBL[MixBLMode[i]];
		}
		else	{
			BaseBLMul[i] = UniBLMul;
		}
		/*
		if (mParam->MBL == 2)	{
			BaseGeneRate[i] = &TotalLength;
		}
		else	{
			BaseGeneRate[i] = &UniGeneRate;
		}
		*/
		BaseGeneRate[i] = &UniGeneRate;
		BaseRefStationary[i] = RefStationary;
		BaseRefMatrix[i] = mSubMatrix;
		BaseRefRateFactor[i] = &RefRateFactor;
	}
	else	{
		if (mParam->ActivateClock)	{
			if (mParam->BLModelSwitch != 1)	{
				cerr << "error in SiteLogSamplingSeparate : works only under flexible models (BLModelSwitch = 0)\n";
				exit(1);
			}
			/*
			if (mParam->GeneBLMultiplier)	{
				cerr << "gene multiplier not yet implemented under clock\n";
				exit(1);
			}
			*/
		}
		int gene = mParam->Gene[i];
		if (mParam->GeneBLMode)	{
			if (mParam->GeneBLMultiplier)	{
				if (mParam->BLModelSwitch == 0)	{
					cerr << "error: gene BL multipliers only with BLModelSwitch == 1\n";
					exit(1);
				}
				if (mParam->ActivateClock)	{
					BaseBL[i] = RigidBL;
				}
				else	{
					BaseBL[i] = BL;
				}
				BaseBLMul[i] = GeneBL[gene];
			}
			else	{
				if (mParam->BLModelSwitch == 0)	{
					BaseBL[i] = GeneRigidBL[gene];
					BaseBLMul[i] = UniBLMul;
				}
				else	{
					BaseBL[i] = GeneBL[gene];
					BaseBLMul[i] = UniBLMul;
				}
			}
			/*
			if (mParam->GeneBLMultiplier)	{
				BaseBL[i] = BL;
				BaseBLMul[i] = GeneBL[gene];
			}
			else	{
				BaseBL[i] = GeneBL[gene];
				BaseBLMul[i] = UniBLMul;
			}
			*/
		}
		else	{
			if (mParam->BLModelSwitch == 0)	{
				BaseBL[i] = RigidBL;
			}
			else	{
				BaseBL[i] = BL;
			}
			/*
			BaseBL[i] = BL;
			BaseBLMul[i] = UniBLMul;
			*/
		}
		if (mParam->GeneRateMode)	{
			BaseGeneRate[i] = &GeneRate[gene];
		}
		else	{
			BaseGeneRate[i] = &UniGeneRate;
		}
		if (mParam->GeneStationaryMode)	{
			BaseRefStationary[i] = GeneStationary[gene];
			BaseRefMatrix[i] = mGeneMatrixArray[gene];
			BaseRefRateFactor[i] = &GeneRateFactor[gene];
		}
		else	{
			BaseRefStationary[i] = RefStationary;
			BaseRefMatrix[i] = mSubMatrix;
			BaseRefRateFactor[i] = &RefRateFactor;
		}
	}

	if (mParam->RASModelSwitch == 0)	{
		BasePconst = &pconst;

		if ((mParam->HeteroMode == Covarion) && (!mParam->ExternalCov))	{
			BaseRate[i] = &One;
		}
		else	{
			if (ConstantStatus[i])	{
				BaseRate[i] = &Zero;
			}
			else	{
				BaseRate[i] = &unirate[i];
			}
		}
	}
	else	{
		BaseRate[i] = &rate[i];
		BasePconst = &pconst;
	}

	if (mParam->SUBModelSwitch == 0)	{
		BaseStationary[i] = BaseRefStationary[i];
		BaseMatrix[i] = BaseRefMatrix[i];
		BaseRateFactor[i] = BaseRefRateFactor[i];
		BasePoisson = mParam->RefPoisson;
	}
	else	{
		int mode = Mode[i];
		if (mParam->HeteroMode == Covarion)	{
			BasePoisson = mParam->ModePoisson;
			int index = mParam->ExternalCov? 0 : RateMode[i];
			if (index <0)	{
				cerr << "negative index\n";
				exit(1);
			}
			if (index >=NRateMode)	{
				cerr << "rate mode index overflow\n";
				cerr << index << '\n';
				exit(1);
			}
			if (BasePoisson)	{
				BaseMatrix[i] = mCovZipMatrixArray[i][index];
			}
			else	{
				BaseMatrix[i] = mCovMatrixArray[mode][index];
			}
			BaseStationary[i] = BaseMatrix[i]->GetStationaries();
		}
		else	{
			BasePoisson = mParam->ModePoisson;
			if (BasePoisson)	{
				BaseMatrix[i] = 0;
				if (mParam->NH)	{
					BaseStationary[i] = NHZipStationary[i][NHcat[root->label]];
				}
				else	{
					BaseStationary[i] = ZipStationary[i];
				}
			}
			else	{
				if (mParam->NH)	{
					if (mParam->ModeFastCompute)	{
						if (mParam->ZipGTR >= 2)	{
							BaseMatrix[i] = mNHMatrixArray[i][NHcat[root->label]];
						}
						else	{
							BaseMatrix[i] = mNHMatrixArray[mode][NHcat[root->label]];
						}
					}
					else	{
						BaseMatrix[i] = mNHMatrixArray[mode][NHcat[root->label]];
					}
				}
				else	{
					if (mParam->ModeFastCompute)	{
						if (mParam->ZipGTR >= 2)	{
							BaseMatrix[i] = mZipMatrixArray[i];
						}
						else	{
							BaseMatrix[i] = mMatrixArray[mode];
						}
					}
					else	{
						BaseMatrix[i] = mMatrixArray[mode];
					}
				}
				BaseStationary[i]= BaseMatrix[i]->GetStationaries();
			}
		}
		BaseRateFactor[i] = &ModeRateFactor[mode];
	}
}




// ---------------------------------------------------------------------------
//		 UpdateZip()
// ---------------------------------------------------------------------------

void
PhyloBayes::UpdateZip()	{

	for (int i=0; i<mParam->Nsite; i++)	{
		UpdateZip(i);
	}
}

// ---------------------------------------------------------------------------
//		 UpdateZip(int i)
// ---------------------------------------------------------------------------


void
PhyloBayes::UpdateZip(int i)	{

	int mode = Mode[i];

	if (mParam->ModeFastCompute)	{

		for (int j=0; j<mParam->OrbitSize[i]; j++)	{
			ZipStationary[i][j] =  Stationary[mode][mParam->Indices[i][j]];
		}
		if (mParam->ZipSize[i]>mParam->OrbitSize[i])	{
			double total1 = 0;
			for (int j=0; j<Nstate; j++)	{
				if (! mParam->Orbit[i][j])	{
					double temp = Stationary[mode][j];
					total1 += temp;
				}
			}
			ZipStationary[i][mParam->OrbitSize[i]] = total1;
		}

		if (mParam->ZipGTR == 4)	{
			// update rr

			double* stat = Stationary[mode];
			double* rr = 0;
			if (mParam->Qmode)	{
				rr = RR[mode];
			}
			else	{
				rr = ModeRR;
			}
			int zipsize = mParam->ZipSize[i];
			
			double num[zipsize][zipsize];
			double denom[zipsize];
			for (int a=0; a<zipsize; a++)	{
				denom[a] = 0;
				for (int b=0; b<zipsize; b++)	{
					num[a][b] = 0;
				}
			}
	
			for (int j=0; j<Nstate; j++)	{
				int a = mParam->ZipIndices[i][j];
				denom[a] += stat[j];		
				for (int k=0; k<Nstate; k++)	{
					int b = mParam->ZipIndices[i][k];
					num[a][b] += stat[j] * stat[k] * rr[Random::rrindex(j,k,Nstate)];
				}
			}

			for (int a=0; a<zipsize; a++)	{
				for (int b=a+1; b<zipsize; b++)	{
					SiteZipRR[i][Random::rrindex(a,b,zipsize)] = num[a][b] / denom[a]/denom[b];
					/*
					if (SiteZipRR[i][Random::rrindex(a,b,zipsize)] < 1e-6)	{
						cerr << "too small zip rr\n";
						cerr << SiteZipRR[i][Random::rrindex(a,b,zipsize)] << '\n';
						cerr << a << '\t' << b << '\t' << zipsize << '\n';
						cerr << mParam->Indices[i][a] << '\t' << mParam->Indices[i][b] << '\t' << zipsize << '\n';
			int u = 0;
			for (int j=0; j<Nstate; j++)	{
				for (int k=j+1; k<Nstate; k++)	{
					cerr << rr[u] << '\t';
					u++;
				}
			}
			cerr << '\n';
						cerr << mode << '\n';
						
						exit(1);
					}
					*/
					
				}
			}
			if (RecomputeMatrix)	{
				mZipMatrixArray[i]->ComputeArray();
			}
		}
		if (mParam->HeteroMode == Covarion)	{
			// Covarion :
			// if rates have changed,
			// and if the Tuffley ans Steel covarion is used 
			// one needs to recompute all the matrices
			if (mParam->ExternalCov)	{
				mCovZipMatrixArray[i][0]->ComputeArray();
			}
			else	{
				for (int k=0; k<NRateMode; k++)	{
					mCovZipMatrixArray[i][k]->ComputeArray();
				}
			}
		}
	}
}

void PhyloBayes::UpdateZipNH(int i, int n)	{

	int mode = Mode[i];
	for (int j=0; j<mParam->OrbitSize[i]; j++)	{
		NHZipStationary[i][n][j] =  NHStationary[mode][n][mParam->Indices[i][j]];
	}
	if (mParam->ZipSize[i]>mParam->OrbitSize[i])	{
		double total1 = 0;
		for (int j=0; j<Nstate; j++)	{
			if (! mParam->Orbit[i][j])	{
				double temp = NHStationary[mode][n][j];
				total1 += temp;
			}
		}
		NHZipStationary[i][n][mParam->OrbitSize[i]] = total1;
	}
}

void PhyloBayes::UpdateZipNH(int i)	{

	for (int n=0; n<NHNcat; n++)	{
		UpdateZipNH(i,n);
	}
}


// ---------------------------------------------------------------------------
//		 ComputeRateFactor()
// ---------------------------------------------------------------------------

void
PhyloBayes::ComputeRateFactor()	{

	for (int i=0; i<GetModeNumber(); i++)	{
		ComputeRateFactor(i);
	}
}


// ---------------------------------------------------------------------------
//		 ComputeRateFactor(int i)
// ---------------------------------------------------------------------------

void
	PhyloBayes::ComputeRateFactor(int i)	{

	if (! mParam->ModeFastCompute)	{
		cerr << "what are you doing here ?? (PhyloBayes::ComputeRateFactor())\n";
		exit(1);
	}

	if (mParam->Normalise)	{
		double temp = 0;
		for (int j=0; j<Nstate-1 ; j++)	{
			for (int k=j+1; k<Nstate; k++)	{
				temp += Stationary[i][j] * Stationary[i][k];
			}
		}
		ModeRateFactor[i] = 1.0 / 2/ temp;
	}
	else	{
		ModeRateFactor[i] = 1.0;
	}
}


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//		 UpdateRef()
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void
PhyloBayes::UpdateRef()	{

	if (! mParam->RefFastCompute)	{
		mSubMatrix->ComputeArray();
	}
	else	{
		UpdateRefZip();
		ComputeRefRateFactor();
	}
}

// ---------------------------------------------------------------------------
//		 UpdateRefZipForAllSites()
// ---------------------------------------------------------------------------

void
PhyloBayes::UpdateRefZip()	{

	for (int i=0; i<mParam->Nsite; i++)	{
		UpdateRefZip(i);
	}
}


// ---------------------------------------------------------------------------
//		 UpdateRefZip(int i)
// ---------------------------------------------------------------------------

void
PhyloBayes::UpdateRefZip(int i)	{

	if (! mParam->RefFastCompute)	{
		cerr << "error : should not be in PhyloBayes::UpdateRefZip\n";
		exit(1);
	}

	for (int j=0; j<mParam->OrbitSize[i]; j++)	{
		RefZipStationary[i][j] =  RefStationary[mParam->Indices[i][j]];
	}

	if (mParam->ZipSize[i]>mParam->OrbitSize[i])	{
		double total1 = 0;
		for (int j=0; j<Nstate; j++)	{
			if (! mParam->Orbit[i][j])	{
				double temp = RefStationary[j];
				total1 += temp;
			}
		}
		RefZipStationary[i][mParam->OrbitSize[i]] = total1;
	}
}


// ---------------------------------------------------------------------------
//		 ComputeRefRateFactor()
// ---------------------------------------------------------------------------

void
	PhyloBayes::ComputeRefRateFactor()	{

	if (! mParam->RefFastCompute)	{
		cerr << "what are you doing here ?? (PhyloBayes::RefComputeRateFactor())\n";
		exit(1);
	}

	if (mParam->Normalise)	{
		double temp = 0;
		for (int j=0; j<Nstate-1 ; j++)	{
			for (int k=j+1; k<Nstate; k++)	{
				temp += RefStationary[j] * RefStationary[k];
			}
		}
		RefRateFactor = 1.0 / 2/ temp;
	}
	else	{
		RefRateFactor = 1.0;
	}
}


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//		 UpdateGene()
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void PhyloBayes::UpdateGene()	{

	// UpdateGeneBranch();
	if (mParam->GeneStationaryMode)	{
		if (! mParam->RefFastCompute)	{
			for (int i=0; i<mParam->Ngene; i++)	{
				mGeneMatrixArray[i]->ComputeArray();
			}
		}
		else	{
			UpdateGeneZipForAllSites();
			ComputeGeneRateFactorForAllGenes();
		}
	}
}

void PhyloBayes::UpdateGeneBranch()	{

	for (int gene=0; gene<mParam->Ngene; gene++)	{
		for (int j=0; j<mParam->Nnode; j++)	{
			GeneBranch[gene][j] = 1;
		}
		UpdateGeneBranch(gene,root);
		GeneBranch[gene][root->label] = 0;
	}
}

int PhyloBayes::UpdateGeneBranch(int gene, Node* node)	{

	int returnValue = 0;
	if (node->isLeaf())	{
		returnValue = mParam->GeneTaxa[gene][node->label];
	}
	else	{
		returnValue |= UpdateGeneBranch(gene,node->left);
		returnValue |= UpdateGeneBranch(gene,node->right);
	}
	GeneBranch[gene][node->label] = returnValue;
	return returnValue;
}

// ---------------------------------------------------------------------------
//		 UpdateGene(int gene)
// ---------------------------------------------------------------------------

void PhyloBayes::UpdateGene(int gene)	{

	if (! mParam->RefFastCompute)	{
		mGeneMatrixArray[gene]->ComputeArray();
	}
	else	{
		for (int i=0; i<mParam->GeneSize[gene]; i++)	{
			UpdateGeneZip(mParam->GeneFirstSite[gene] + i);
		}	
		ComputeGeneRateFactor(gene);
	}

}


// ---------------------------------------------------------------------------
//		 ComputeGeneRateFactorForAllGenes()
// ---------------------------------------------------------------------------

void
PhyloBayes::ComputeGeneRateFactorForAllGenes()	{

	for (int i=0; i<mParam->Ngene; i++)	{
		ComputeGeneRateFactor(i);
	}
}



// ---------------------------------------------------------------------------
//		 ComputeGeneRateFactor(int i)
// ---------------------------------------------------------------------------

void
	PhyloBayes::ComputeGeneRateFactor(int i)	{

	if (mParam->Normalise)	{
		double temp = 0;
		for (int j=0; j<Nstate-1 ; j++)	{
			for (int k=j+1; k<Nstate; k++)	{
				temp += GeneStationary[i][j] * GeneStationary[i][k];
			}
		}
		GeneRateFactor[i] = 1.0 / 2/ temp;
	}
	else	{
		GeneRateFactor[i] = 1.0;
	}
}


// ---------------------------------------------------------------------------
//		 UpdateGeneZipForAllSites()
// ---------------------------------------------------------------------------

void
PhyloBayes::UpdateGeneZipForAllSites()	{

	for (int i=0; i<mParam->Nsite; i++)	{
		UpdateGeneZip(i);
	}
}

// ---------------------------------------------------------------------------
//		 UpdateGeneZip(int i)
// ---------------------------------------------------------------------------

void
PhyloBayes::UpdateGeneZip(int i)	{

	int gene = mParam->Gene[i];

	for (int j=0; j<mParam->OrbitSize[i]; j++)	{
		GeneZipStationary[i][j] =  GeneStationary[gene][mParam->Indices[i][j]];
	}

	if (mParam->ZipSize[i]>mParam->OrbitSize[i])	{
		double total1 = 0;
		for (int j=0; j<Nstate; j++)	{
			if (! mParam->Orbit[i][j])	{
				double temp = GeneStationary[gene][j];
				total1 += temp;
			}
		}
		GeneZipStationary[i][mParam->OrbitSize[i]] = total1;
	}
}



// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//		 UpdateLogProbs()
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void PhyloBayes::UpdateLogProbs()	{

	if (!mParam->NormalApprox)	{
		for (int i=0; i<mParam->Nnode; i++)	{
			tbloffset[i] = 0;
		}
	}

	if (!mParam->NormalApprox && (mParam->ZipSub == 2))	{
		for (int i=0; i<mParam->Nsite; i++)	{
			UpdateMissingFlag(root->label,i);
		}	
	}

	if (Beta == 0)	{
		mLogSampling = 0;
		mConcatenateNormalLogSampling = 0;
		mSeparateNormalLogSampling = 0;
		mFlexGeneLogSampling = 0;
		mRigidGeneLogSampling = 0;
	}

	mLogPrior = -1;
	mLogSampling = -1;
	mLogPosterior = -1;
	logSampling();
	logPrior();
	logPosterior();
}

// ***************************************************************************
// ***************************************************************************
//		 Accessors
// ***************************************************************************
// ***************************************************************************


// ---------------------------------------------------------------------------
//		 GetRefRR(int i, int j) mode i, state j
// ---------------------------------------------------------------------------

double&
	PhyloBayes::GetRefRR(int i, int j)	{

		return (i<j) ?
			RefRR[(2 * Nstate - i - 1) * i / 2 + j - i - 1] :
			RefRR[(2 * Nstate - j - 1) * j / 2 + i - j - 1] ;
	}



// ---------------------------------------------------------------------------
//		 GetRR(int mode, int i, int j) mode i, state j
// ---------------------------------------------------------------------------

double&
	PhyloBayes::GetRR(int mode, int i, int j)	{

		return (i<j) ?
			RR[mode][(2 * Nstate - i - 1) * i / 2 + j - i - 1] :
			RR[mode][(2 * Nstate - j - 1) * j / 2 + i - j - 1] ;
	}


// ---------------------------------------------------------------------------
//		 GetModeRR(int i, int j) mode i, state j
// ---------------------------------------------------------------------------

double&
	PhyloBayes::GetModeRR(int i, int j)	{

		return (i<j) ?
			ModeRR[(2 * Nstate - i - 1) * i / 2 + j - i - 1] :
			ModeRR[(2 * Nstate - j - 1) * j / 2 + i - j - 1] ;
	}


// ---------------------------------------------------------------------------
//		 GetStationary(int i, int j) mode i, state j
// ---------------------------------------------------------------------------

double&
	PhyloBayes::GetStationary(int i, int j)	{

		return Stationary[i][j];

	}


// ---------------------------------------------------------------------------
//		 GetSiteStationary(int i, int j) mode i, state j
// ---------------------------------------------------------------------------

double
	PhyloBayes::GetSiteStationary(int i, int j)	{

		return GetStationary(Mode[i],j);
	}

// ---------------------------------------------------------------------------
//		 SetBeta()
// ---------------------------------------------------------------------------

void
PhyloBayes::SetBeta(double inBeta)	{

	Beta = inBeta;
	mLogPosterior = mLogPrior + Beta * mLogSampling;
}


// ---------------------------------------------------------------------------
//		 SetData()
// ---------------------------------------------------------------------------

void
PhyloBayes::SetData(int chain)	{

	if (chain == -1)	{
		Data = mParam->Data;
		ZipData = mParam->ZipData;
	}
	else	{
		cerr << "error : cannot set data in PhyloBayes::SetData\n";
		exit(1);
	}

}

// ---------------------------------------------------------------------------
//		 SwapMatrices()
// ---------------------------------------------------------------------------

void
PhyloBayes::SwapMatrices(PhyloBayes* from)	{

	if (mParam->ModeFastCompute)	{
		cerr << "error in swap matrices : mode fast compute\n";
		exit(1);
	}
	for (int i=0; i<mParam->NmodeMax; i++)	{
		SubMatrix* temp = mMatrixArray[i];
		mMatrixArray[i] = from->mMatrixArray[i];
		from->mMatrixArray[i] = temp;
	}
}


// ***************************************************************************
// ***************************************************************************
//		 Marginals for monitoring
// ***************************************************************************
// ***************************************************************************


// ---------------------------------------------------------------------------
//		 GetMeanBranchStat()
// ---------------------------------------------------------------------------

void PhyloBayes::GetMeanBranchStat(double* avstat, int j)	{

	for (int k=0; k<Nstate; k++)	{
		avstat[k] = 0;
	}
	for (int i=0; i<Nmode; i++)	{
		double* stat = NHStationary[i][NHcat[j]];
		double weight = SiteNumber[i];
		for (int k=0; k<Nstate; k++)	{
			avstat[k] += weight * stat[k];
		}
	}
	for (int k=0; k<Nstate; k++)	{
		avstat[k] /= mParam->Nsite;
	}
}

// ---------------------------------------------------------------------------
//		 GetMeanPopEffSize
// ---------------------------------------------------------------------------

double PhyloBayes::GetRootPopEff()	{
	return PopEffSize[NHcat[root->label]];
}

double PhyloBayes::GetMeanPopEffSize()	{

	double total = 0;
	for (int j=0; j<NHNcat; j++)	{
		total += PopEffSize[j];
	}
	return total / NHNcat;
}

double PhyloBayes::GetVarPopEffSize()	{

	double mean = 0;
	double var = 0;
	for (int j=0; j<NHNcat; j++)	{
		mean += PopEffSize[j];
		var += PopEffSize[j] * PopEffSize[j];
	}
	mean /= NHNcat;
	var /= NHNcat;
	var -= mean * mean;
	return var;
}

// ---------------------------------------------------------------------------
//		 GetMeanContChar
// ---------------------------------------------------------------------------

double PhyloBayes::GetMeanContChar(int k)	{

	double total = 0;
	for (int j=0; j<mParam->Nnode; j++)	{
		total += ContRho[k][j];
	}
	return total/mParam->Nnode;
}

double PhyloBayes::GetMeanLeafContChar(int k)	{

	double total = 0;
	for (int j=0; j<mParam->Ntaxa; j++)	{
		total += ContRho[k][j];
	}
	return total/mParam->Ntaxa;
}


// ---------------------------------------------------------------------------
//		 GetMeanOrbitSize
// ---------------------------------------------------------------------------

double PhyloBayes::GetMeanOrbitSize()	{

	double total = 0;
	if (mParam->ZipGTR == 4)	{
		cerr << "errpr in get mean orbit size\n";
		exit(1);
	}
	if (mParam->ZipGTR >= 2)	{
		for (int i=0; i<mParam->Nsite; i++)	{
			total += SiteZipSize[i];
		}
	}
	else	{
		for (int i=0; i<Nmode; i++)	{
			total += SiteNumber[i] * ModeZipSize[i];
		}
	}
	return total / mParam->Nsite;
}

double PhyloBayes::GetMeanConstantSiteOrbitSize()	{

	double total = 0;
	int nconst = 0;
	if ((mParam->ZipGTR != 4) && (mParam->ZipGTR >= 2))	{
		for (int i=0; i<mParam->Nsite; i++)	{
			if (mParam->OrbitSize[i] == 1)	{
				total += SiteZipSize[i];
				nconst ++;
			}
		}
	}
	return total / mParam->Nsite;
}


// ---------------------------------------------------------------------------
//		 GetSpeedRatio
// ---------------------------------------------------------------------------

double PhyloBayes::GetSpeedRatio()	{

	double total = 0;
	if ((mParam->ZipGTR != 4) && (mParam->ZipGTR >= 2))	{
		for (int i=0; i<mParam->Nsite; i++)	{
			double temp = ((double) SiteZipSize[i]) / Nstate;
			total +=  temp * temp;
		}
	}
	else	{
		for (int i=0; i<Nmode; i++)	{
			double temp = ((double) ModeZipSize[i]) / Nstate;
			total += SiteNumber[i] * temp * temp;
		}
	}
	return total / mParam->Nsite;
}

double PhyloBayes::GetOptimalSpeedRatio()	{

	double total = 0;
	for (int mode=0; mode<Nmode; mode++)	{
		int aa[Nstate];
		for (int k=0; k<Nstate; k++)	{
			aa[k] = 0;
		}
		for (int i=0; i<mParam->Nsite; i++)	{
			if (Mode[i] == mode)	{
				for (int j=0; j<mParam->Ntaxa; j++)	{
					if (mParam->Data[j][i] != unknown)	{
						aa[mParam->Data[j][i]] = 1;
					}
				}
			}
		}
		double size = 0;
		for (int k=0; k<Nstate; k++)	{
			size += aa[k];
		}
		total += SiteNumber[mode] * ((double) size * size) / Nstate / Nstate;
	}
	return total / mParam->Nsite;
}

// ---------------------------------------------------------------------------
//		 GetModeAffEntropy()
// ---------------------------------------------------------------------------

double PhyloBayes::GetModeAffEntropy()	{

	UpdateSiteNumber();
	double total = 0;
	for (int i=0; i<Nmode; i++)	{
		double temp = ((double) SiteNumber[i]) / mParam->Nsite;
		if (temp > 1e-6)	{
			total -= temp*log(temp);
		}
	}
	return total;
}

// ---------------------------------------------------------------------------
//		 GetEffNMixBL()
// ---------------------------------------------------------------------------

double PhyloBayes::GetEffNMixBL()	{

	UpdateMBL();
	double total = 0;
	for (int i=0; i<NMixBL; i++)	{
		double temp = ((double) MixBLSiteNumber[i]) / mParam->Nsite;
		total -= temp * log(temp);
	}
	return exp(total);
}

// ---------------------------------------------------------------------------
//		 GetEffNmode()
// ---------------------------------------------------------------------------

int PhyloBayes::GetNoccupied()	{

	UpdateSiteNumber();
	int total = 0;
	for (int i=0; i<Nmode; i++)	{
		if (SiteNumber[i])	{
			total++;
		}
	}
	return total;
}

double PhyloBayes::GetEffNmode()	{

	UpdateSiteNumber();
	double total = 0;
	for (int i=0; i<Nmode; i++)	{
		double temp = ((double) SiteNumber[i]) / mParam->Nsite;
		total -= temp * log(temp);
	}
	return exp(total);
}
		
// ---------------------------------------------------------------------------
//		 GetMeanGeneBL()
// ---------------------------------------------------------------------------

double
	PhyloBayes::GetMeanGeneBL()	{

	double total = 0;
	int n=0;
	for (int gene=0; gene<mParam->Ngene; gene++)	{
		for (int j=0; j<mParam->Nnode; j++)	{
			if (blfree(j))	{
				total += GeneBL[gene][j];
				n++;
			}
		}
	}
	return total / n;
}

// ---------------------------------------------------------------------------
//		 GetMeanMixBL()
// ---------------------------------------------------------------------------

double
	PhyloBayes::GetMeanMixBL()	{

	double total = 0;
	int n=0;
	for (int i=0; i<NMixBL; i++)	{
		for (int j=0; j<mParam->Nnode; j++)	{
			if (blfree(j))	{
				total += MixBL[i][j];
				n++;
			}
		}
	}
	return total / NMixBL;
}

// ---------------------------------------------------------------------------
//		 GetRateEntropy()
// ---------------------------------------------------------------------------

double
	PhyloBayes::GetRateEntropy()	{

		double total = 0;
		for (int i=0; i<mParam->Nsite; i++)	{
			if (mParam->RASModelSwitch)	{
				total += rate[i];
			}
			else	{
				total += ModeRate[i];
			}
		}

		double temp = 0;
		for (int i=0; i<mParam->Nsite; i++)	{
			double t = 0;
			if (mParam->RASModelSwitch)	{
				t = rate[i] / total;
			}
			else	{
				t = ModeRate[i] / total;
			}
			temp += ( t < mParam->TooSmall) ? 0 : -t * log (t);
		}
		return temp;
	}

// ---------------------------------------------------------------------------
//		 GetMeanModeRate()
// ---------------------------------------------------------------------------

double
	PhyloBayes::GetMeanModeRate()	{
		double total = 0;
		int weight = 0;
		for (int i=0; i<NRateMode; i++)	{
			total += ModeRate[i] * RateSiteNumber[i];
			weight += RateSiteNumber[i];
		}
		if (weight != mParam->Nsite)	{
			cerr << "error : total site number is not equal to Nsite\n";
		}
		return total / mParam->Nsite;
	}

// ---------------------------------------------------------------------------
//		 GetMeanRate()
// ---------------------------------------------------------------------------

double
	PhyloBayes::GetMeanRate()	{

		double temp = 0;
		if (mParam->SumOverRateModes)	{
			temp = 0;
			double weight = 0;
			for (int i=0; i<NRateMode; i++)	{
				temp += RateModeWeight[i] * ModeRate[i];
				weight += RateModeWeight[i];
			}
			temp /= weight;
		}
		else	{
			if (mParam->RASModelSwitch == 0)	{
				MakeRates();
			}
			for (int i=0; i<mParam->Nsite; i++)	{
				temp += (mParam->RASModelSwitch == 1) ? rate[i] : unirate[i];
			}
			temp /= mParam->Nsite;
		}
		return temp;
	}


// ---------------------------------------------------------------------------
//		 GetMeanNHVar()
// ---------------------------------------------------------------------------

double PhyloBayes::GetMeanNHVar()	{

	if (mParam->NH == 2){
		if (mParam->NHWishart)	{
			logNHCovProb();
			return NHTrace / ContNstate;
		}
	}
	double total = 0;
	for (int k=0; k<ContNstate; k++)	{
		total += NHVar[k];
	}
	return total / ContNstate;
}

double PhyloBayes::GetMeanNHContVar()	{

	double total = 0;
	for (int k=0; k<mParam->ContNchar+1; k++)	{
		total += NHContVar[k];
	}
	return total / (mParam->ContNchar+1);
}

double PhyloBayes::GetMeanNHStatCenter()	{

	double total = 0;
	for (int k=0; k<ContNstate; k++)	{
		total += NHStatCenter[k];
	}
	return total / ContNstate;
}

double PhyloBayes::GetVarNHStatCenter()	{

	double total = 0;
	for (int k=0; k<ContNstate; k++)	{
		total += NHStatCenter[k] * NHStatCenter[k];
	}
	return total / ContNstate;
}

double PhyloBayes::GetNHDim()	{
	if (mParam->NHPrior == 2)	{
		if (mParam->NHWishart)	{
			logNHCovProb();
		}
		return NHDim;
	}
	double total = 0;
	for (int k=0; k<ContNstate; k++)	{
		total += NHVar[k];
	}
	double entropy = 0;
	for (int k=0; k<ContNstate; k++)	{
		double tmp = NHVar[k] / total;
		if (tmp > 1e-8)	{
			entropy -= tmp * log(tmp);
		}
	}
	return exp(entropy);
}

// ---------------------------------------------------------------------------
//		 GetNHStatEntropy()
// ---------------------------------------------------------------------------

double PhyloBayes::GetNHStatEntropy()	{

	double total = 0;
	double totalweight = 0;
	for (int n=0; n<NHNcat; n++)	{
		double temp = 0;
		for (int k=0; k<Nstate; k++)	{
			double tmp = NHStatDistorter[n][k];
			if (tmp > 1e-8)	{
				temp -= tmp * log(tmp);
			}
		}
		total += NHWeight[n] * temp;
		totalweight += NHWeight[n];
	}
	return total / totalweight;
}

// ---------------------------------------------------------------------------
//		 GetNHWeightEntropy()
// ---------------------------------------------------------------------------

double PhyloBayes::GetNHWeightEntropy()	{
	
	double totalweight = 0;
	for (int n=0; n<NHNcat; n++)	{
		totalweight += NHWeight[n];
	}
	double total = 0;
	for (int n=0; n<NHNcat; n++)	{
		double tmp = NHWeight[n] / totalweight;
		if (tmp > 1e-8)	{
			total -= tmp * log(tmp);
		}
	}
	return total;
}

int PhyloBayes::GetNHActiveCatNumber()	{

	int sitenumber[NHNcat];
	for (int n=0; n<NHNcat; n++)	{
		sitenumber[n] = 0;
	}
	for (int j=0; j<mParam->Nnode; j++)	{
		sitenumber[NHcat[j]]++;
	}
	int N = 0;
	for (int n=0; n<NHNcat; n++)	{
		if (sitenumber[n]) N++;
	}
	return N;
}


// ---------------------------------------------------------------------------
//		 GetStationaryEntropy()
// ---------------------------------------------------------------------------

double PhyloBayes::GetFiniteTimeEntropy(double time, int mode)	{


	UpdateMode(mode);

	double meanent = 0;

	double** temp = new double*[Nstate];
	for (int i=0; i<Nstate; i++)	{
		temp[i] = new double[Nstate];
	}

	if (mParam->ModeFastCompute)	{
		
		for (int k=0; k<Nstate; k++)	{
			for (int l=0; l<Nstate; l++)	{
				temp[k][l] = (1 - exp(-time)) * Stationary[mode][l];
				if (k == l)	{
					temp[k][l] += exp(-time);
				}
			}
		}
	}
	else	{

		mMatrixArray[mode]->ComputeExponential(time,temp);
		
	}

	for (int k=0; k<Nstate; k++)	{
		double ent = 0;
		for (int l=0; l<Nstate; l++)	{
			ent -= temp[l][k] * log(temp[l][k]);
		}
		meanent += Stationary[mode][k] * ent;
	}

	for (int i=0; i<Nstate; i++)	{
		delete[] temp[i];
	}
	delete[] temp;

	return meanent;
}

double PhyloBayes::GetFiniteTimeEntropy(double time)	{

	double mean = 0;
	for (int mode=0; mode<Nmode; mode++)	{
		mean += GetFiniteTimeEntropy(time,mode);
	}
	mean /= Nmode;
	return mean;
}


double PhyloBayes::GetStatEnt()	{
	return mParam->SUBModelSwitch ? GetStationaryEntropy() : GetRefStationaryEntropy();
}

double PhyloBayes::GetStatEnt(int site)	{
	
	double temp = 0;
	double* stat = mParam->SUBModelSwitch ? Stationary[Mode[site]] : RefStationary;
	for (int j=0; j<Nstate; j++)	{
		double t = stat[j];
		temp += (t <mParam->TooSmall) ? 0 : -t*log(t);
	}
	return temp;
}

double
	PhyloBayes::GetStationaryEntropy()	{

	if (!mParam->SumOverModes)	{
		UpdateSiteNumber();
		double temp = 0;
		for (int i=0; i<Nmode; i++)	{
			for (int j=0; j<Nstate; j++)	{
				double t = Stationary[i][j];
				temp += SiteNumber[i] * ( ( t < mParam->TooSmall) ? 0 : -t * log (t));
			}
		}
		return temp / mParam->Nsite;
	}
	else	{

		double total = 0;
		for (int i=0; i<Nmode; i++)	{
			double temp = 0;
			for (int j=0; j<Nstate; j++)	{
				double t = Stationary[i][j];
				temp +=  ( ( t < mParam->TooSmall) ? 0 : -t * log (t));
			}
			total += ModeWeight[i] * temp;
		}
		return total;
	}
	return 0;
		
}


double PhyloBayes::GetSiteStationaryEntropy(int site)	{

	double total = 0;
	int mode = Mode[site];
	for (int j=0; j<Nstate; j++)	{
		double t = Stationary[mode][j];
		total += ( ( t < mParam->TooSmall) ? 0 : -t * log (t));
	}
	return total;
}

double PhyloBayes::GetStatCenterEntropy()	{
	
	double temp = 0;
	for (int j=0; j<Nstate; j++)	{
		double t = ModeStatCenter[j];
		temp += (t <mParam->TooSmall) ? 0 : -t*log(t);
	}
	return temp;
}



// ---------------------------------------------------------------------------
//		 GetRefStationaryEntropy()
// ---------------------------------------------------------------------------

double
	PhyloBayes::GetRefStationaryEntropy()	{

		double temp = 0;
		for (int j=0; j<Nstate; j++)	{
			double t = RefStationary[j];
			temp += ( ( t < mParam->TooSmall) ? 0 : -t * log (t));
		}
		return temp;
	}



// ---------------------------------------------------------------------------
//		 GetGeneStationaryEntropy()
// ---------------------------------------------------------------------------

double
	PhyloBayes::GetGeneStationaryEntropy()	{

		double temp = 0;
		for (int i=0; i<mParam->Ngene; i++)	{
			for (int j=0; j<Nstate; j++)	{
				double t = GeneStationary[i][j];
				temp += mParam->GeneSize[i] * ( ( t < mParam->TooSmall) ? 0 : -t * log (t));
			}
		}
		return temp / mParam->Nsite;
	}




// ---------------------------------------------------------------------------
//		 GetRRACGT_Entropy()
// ---------------------------------------------------------------------------

double PhyloBayes::GetRRACGT_Entropy()	{

	int Nrr = Nnuc * (Nnuc-1) / 2;
	double norm = 0;
	for (int i=0; i<Nrr; i++)	{
		norm += RefRRACGT[i];
	}
	double total = 0;
	for (int i=0; i<Nrr; i++)	{
		double tmp = RefRRACGT[i] / norm;
		if (tmp > 1e-6)	{
			total -= tmp * log(tmp);
		}
	}
	return total;
}	

double PhyloBayes::GetMeanRRACGT()	{

	int Nrr = Nnuc * (Nnuc-1) / 2;
	double norm = 0;
	for (int i=0; i<Nrr; i++)	{
		norm += RefRRACGT[i];
	}
	return norm / Nrr;
}	



// ---------------------------------------------------------------------------
//		 GetRREntropy()
// ---------------------------------------------------------------------------

double
	PhyloBayes::GetRREntropy()	{

		if (mParam->SUBModelSwitch == 0)	{
			return GetRefRREntropy();
		}

		double temp = 0;
		if (mParam->Qmode)	{
			for (int mode=0; mode<Nmode; mode++)	{
				double total = 0;
				for (int i=0; i<Nstate; i++)	{
					for (int j=i+1; j<Nstate; j++)	{
						total += RR[mode][(2 * Nstate - i - 1) * i / 2 + j - i - 1] ;
					}
				}
				for (int i=0; i<Nstate; i++)	{
					for (int j=i+1; j<Nstate; j++)	{
						double t = RR[mode][(2 * Nstate - i - 1) * i / 2 + j - i - 1]  / total;
						temp +=  ( t < mParam->TooSmall) ? 0 : -t * log (t);
					}
				}
			}
			temp /= Nmode;
		}
		else	{
			double total = 0;
			for (int i=0; i<Nstate; i++)	{
				for (int j=i+1; j<Nstate; j++)	{
					total += ModeRR[(2 * Nstate - i - 1) * i / 2 + j - i - 1] ;
				}
			}
			for (int i=0; i<Nstate; i++)	{
				for (int j=i+1; j<Nstate; j++)	{
					double t = ModeRR[(2 * Nstate - i - 1) * i / 2 + j - i - 1]  / total;
					temp +=  ( t < mParam->TooSmall) ? 0 : -t * log (t);
				}
			}
		}
		return temp;
	}


// ---------------------------------------------------------------------------
//		 GetMeanRR()
// ---------------------------------------------------------------------------

double
	PhyloBayes::GetMeanRR()	{

		double total = 0;
		if (mParam->Qmode)	{
			for (int mode=0; mode<Nmode; mode++)	{
				for (int i=0; i<Nstate; i++)	{
					for (int j=i+1; j<Nstate; j++)	{
						total += RR[mode][(2 * Nstate - i - 1) * i / 2 + j - i - 1] ;
					}
				}
			}
			total /= Nmode * Nstate * (Nstate - 1) / 2;
		}
		else	{
			for (int i=0; i<Nstate; i++)	{
				for (int j=i+1; j<Nstate; j++)	{
					total += ModeRR[(2 * Nstate - i - 1) * i / 2 + j - i - 1] ;
				}
			}
			total /= Nstate * (Nstate - 1) / 2;
		}
		return total;
	}

double
	PhyloBayes::GetMeanRefRR()	{

		double total = 0;
		for (int i=0; i<Nstate; i++)	{
			for (int j=i+1; j<Nstate; j++)	{
				total += RefRR[(2 * Nstate - i - 1) * i / 2 + j - i - 1] ;
			}
		}
		total /= Nstate * (Nstate - 1) / 2;
		return total;
	}

// ---------------------------------------------------------------------------
//		 GetRefRREntropy()
// ---------------------------------------------------------------------------

double
	PhyloBayes::GetRefRREntropy()	{

		double total = 0;
		for (int i=0; i<Nstate; i++)	{
			for (int j=i+1; j<Nstate; j++)	{
				total += RefRR[(2 * Nstate - i - 1) * i / 2 + j - i - 1] ;
			}
		}
		double temp = 0;
		for (int i=0; i<Nstate; i++)	{
			for (int j=i+1; j<Nstate; j++)	{
				double t = RefRR[(2 * Nstate - i - 1) * i / 2 + j - i - 1]  / total;
				temp +=  ( t < mParam->TooSmall) ? 0 : -t * log (t);
			}
		}
		return temp;
	}

// ---------------------------------------------------------------------------
//		 GetLength()
// ---------------------------------------------------------------------------

double PhyloBayes::GetLength()	{
	double total = 0;
	for (int j=0; j<mParam->Nnode; j++)	{
		total += BL[j];
	}	
	return total;
}

double PhyloBayes::GetRigidLength()	{
	double total = 0;
	for (int j=0; j<mParam->Nnode; j++)	{
		total += RigidBL[j];
	}	
	return total;
}

double PhyloBayes::GetGeneLength()	{
	double total = 0;
	if (mParam->Ngene)	{
		for (int k=0; k<mParam->Ngene; k++)	{
			total += GetGeneLength(k);
		}
		return total / mParam->Ngene;
	}
	return 0;
}

double PhyloBayes::GetGeneLength(int gene)	{
	double total = 0;
	for (int j=0; j<mParam->Nnode; j++)	{
		total += GeneBL[gene][j];
	}	
	return total;
}

double PhyloBayes::GetGeneRigidLength()	{
	double total = 0;
	if (mParam->Ngene)	{
		for (int k=0; k<mParam->Ngene; k++)	{
			total += GetGeneRigidLength(k);
		}
		return total / mParam->Ngene;
	}
	return 0;
}

double PhyloBayes::GetGeneRigidLength(int k)	{
	double total = 0;
	for (int j=0; j<mParam->Nnode; j++)	{
		total += GeneRigidBL[k][j];
	}	
	return total;
}

// ---------------------------------------------------------------------------
//		 GetMeanSigma()
// ---------------------------------------------------------------------------

double PhyloBayes::GetMeanSigma()	{
	double total = 0;
	if (mParam->SeparateRhoPrior)	{
		for (int k=0; k<mParam->Ngene; k++)	{
			total += GeneSigma[k];
		}
		total /= mParam->Ngene;
	}
	else	{
		total = Sigma;
	}
	return total;
}

// ---------------------------------------------------------------------------
//		 GetMeanTheta()
// ---------------------------------------------------------------------------

double PhyloBayes::GetMeanTheta()	{
	double total = 0;
	if (mParam->SeparateRhoPrior)	{
		for (int k=0; k<mParam->Ngene; k++)	{
			total += GeneTheta[k];
		}
		total /= mParam->Ngene;
	}
	else	{
		total = Theta;
	}
	return total;
}

// ---------------------------------------------------------------------------
//		 HomoplasyPerSite()
// ---------------------------------------------------------------------------

double PhyloBayes::HomoplasyPerSite()	{

	int** nsub = BasePoisson ? TrueSub : Nsub;
	double total = 0;
	for (int i=0; i<mParam->Nsite; i++)	{
		for (int k=0; k<Nstate; k++)	{
			if (nsub[i][k])	{
				total += nsub[i][k] - 1;
			}
		}
	}
	return total / mParam->Nsite;
}

double PhyloBayes::NTrueSubPerSite(double* tmp)	{

	int** nsub = BasePoisson ? TrueSub : Nsub;
	double total = 0;
	for (int i=0; i<mParam->Nsite; i++)	{
		double tot = 0;
		for (int k=0; k<Nstate; k++)	{
			tot += nsub[i][k];
		}
		total += tot;
		if (tmp)	{
			tmp[i] = tot;
		}
	}
	return total / mParam->Nsite;
}



// ---------------------------------------------------------------------------
//		 GetModeWeightEntropy()
// ---------------------------------------------------------------------------

double PhyloBayes::GetModeWeightEntropy()	{
	double total = 0;
	double totalweight = 0;
	for (int mode=0; mode<Nmode; mode++)	{
		if (ModeWeight[mode] > 1e-6)	{
			total -= ModeWeight[mode] * log(ModeWeight[mode]);
		}
		totalweight += ModeWeight[mode];
	}
	return total;
}

double PhyloBayes::GetCATAlphaDerivative()	{

	double d = Nmode;
	for (int i=0; i<mParam->Nsite; i++)	{
		d -= alpha / (alpha + i);
	}
	d *= log(mParam->AlphaMax) - log(mParam->AlphaMin);
	return d;
}

/*
void PhyloBayes::MakeCurrentProfileLogo(string filename)	{

	int pageNumber = 1;
	int lineNumber = 0;
	ostringstream appel;
	appel <<  "cp " << auxpath + modeheader << " " << filename;
	system(appel.str().c_str());

	ofstream logo_os(filename.c_str(), IOS_APPEND);
	logo_os << "%%Page: " << 1 << ' ' << 1 << '\n';
	logo_os << "startpage\n";
	logo_os << "startline\n";
	int Nstate = mParam->Nstate;
	double array[Nstate];
	
	for (int i=0; i<Nmode; i++)	{
		for (int k=0; k<Nstate; k++)	{

			logo_os << "gsave\n";

			array[k] = Stationary[i][k] * MaxHeight ;
			if (array[k] > logoThreshold)	{
				logo_os << array[k] << " (" << AminoAcids[k] << ") numchar\n";
			}

			logo_os << "grestore\n";
			logo_os << "shift\n";
		}
		
		logo_os << "endline\n";

		logo_os << "0 0 0 setrgbcolor\n";
		logo_os << "240 30 moveto\n";
		
		ostringstream s;
		s << '(' << SiteNumber[i] << ") show\n" ;
		logo_os << s.str();
		
		lineNumber ++;
		if (lineNumber == LinePerPage)	{
			lineNumber = 0;
			pageNumber++;
			logo_os << "endpage\n";

			logo_os << "%%Page: " << pageNumber << ' ' << pageNumber << '\n';
			logo_os << "startpage\n";
			logo_os << "startline\n";
		}
		else		{
			logo_os << " 0 " << SpaceBetweenLines << " translate\n";
			logo_os << "startline\n";
		}
	}
	logo_os << "endline\n";
	logo_os << "endpage\n";

	logo_os << "%%Trailer\n";
	logo_os << "%%Pages: " << pageNumber << "\n";

	logo_os.close();

}
*/
void PhyloBayes::SetRateMode(int i, int ratemode)	{

	// site i is set to rate-category number "ratemode"

	if (ratemode<0)	{
		cerr << "error in SetRateMode: mode out of range : " << ratemode << '\n';
		exit(1);
	}

	RateMode[i] =ratemode;
	unirate[i] = ModeRate[ratemode];

	// maintain Base Fields
	if ((mParam->HeteroMode == Covarion) && (!mParam->ExternalCov))	{
		BaseRate[i] = &One;
		if (BasePoisson)	{
			BaseMatrix[i] = mCovZipMatrixArray[i][RateMode[i]];
		}
		else	{
			BaseMatrix[i] = mCovMatrixArray[Mode[i]][RateMode[i]];
		}
		BaseStationary[i] = BaseMatrix[i]->GetStationaries();
	}
	else	{
		if (ConstantStatus[i])	{
			BaseRate[i] = &Zero;
		}
		else	{
			BaseRate[i] = &unirate[i];	
		}
	}
}

void PhyloBayes::SetMode(int i, int mode)	{

	// site i is set to profile-category number "mode"

	if (mode<0) {
		cerr << "error in SetMode: mode out of range : " << mode << '\n';
		exit(1);
	} 
	Mode[i] = mode;
	if (mParam->ModeFastCompute)	{
		if ((mParam->ZipGTR == 4) || (mParam->ZipGTR == 0))	{
			UpdateZip(i);
		}
	}
	// maintain Base Fields
	if (mParam->HeteroMode == Covarion)	{
		BasePoisson = mParam->ModePoisson;
		if (BasePoisson)	{
			if (mParam->ExternalCov)	{
				BaseMatrix[i] = mCovZipMatrixArray[i][0];
			}
			else	{
				BaseMatrix[i] = mCovZipMatrixArray[i][RateMode[i]];
			}
		}
		else	{
			if (mParam->ExternalCov)	{
				BaseMatrix[i] = mCovMatrixArray[mode][0];
			}
			else	{
				BaseMatrix[i] = mCovMatrixArray[mode][RateMode[i]];
			}
		}
		BaseStationary[i] = BaseMatrix[i]->GetStationaries();
	}
	else	{
		BasePoisson = mParam->ModePoisson;
		if (BasePoisson)	{
			BaseMatrix[i] = 0;
			if (mParam->NH)	{
				if (mParam->NH)	{
					for (int n=0; n<NHNcat; n++)	{
						UpdateZipNH(i,n);
					}
				}
				BaseStationary[i] = NHZipStationary[i][NHcat[root->label]];
			}
			else	{
				BaseStationary[i] = ZipStationary[i];
			}
		}
		else	{
			if (mParam->NH)	{
				if (mParam->ModeFastCompute)	{
					if (mParam->ZipGTR >= 2)	{
						BaseMatrix[i] = mNHMatrixArray[i][NHcat[root->label]];
					}
					else	{
						BaseMatrix[i] = mNHMatrixArray[mode][NHcat[root->label]];
					}
				}
				else	{
					BaseMatrix[i] = mNHMatrixArray[mode][NHcat[root->label]];
				}
			}
			else	{
				if (mParam->ModeFastCompute)	{
					if (mParam->ZipGTR >= 2)	{
						BaseMatrix[i] = mZipMatrixArray[i];
					}
					else	{
						BaseMatrix[i] = mMatrixArray[mode];
					}
				}
				else	{
					BaseMatrix[i] = mMatrixArray[mode];
				}
			}
			BaseStationary[i]= BaseMatrix[i]->GetStationaries();
		}
	}
	BaseRateFactor[i] = &ModeRateFactor[mode];
}
	
// ***************************************************************************
// ***************************************************************************
//		 Manipulating Tree structure
// ***************************************************************************
// ***************************************************************************


// ---------------------------------------------------------------------------
//		 Phylip
// ---------------------------------------------------------------------------

void	PhyloBayes::Phylip(ostream& os, int withLengths, int withProbs, int withSpeciesNames, int withInternalLabels)	 {

	Tree tree(this);
	tree.Phylip(os, withLengths, withProbs, withSpeciesNames, withInternalLabels);

}


// ---------------------------------------------------------------------------
//		¥ ToMrBayes
// ---------------------------------------------------------------------------

void	PhyloBayes::ToMrBayes(ostream& os) {

	Tree tree(this);
	tree.ToMrBayes(os);

}


// ---------------------------------------------------------------------------
//		 SetTree
// ---------------------------------------------------------------------------

void
	PhyloBayes::ReadPhylipFromStream(istream& is)	{

	SetTree(is);
}


void
	PhyloBayes::SetTree(istream& is)	{

	// get the tree from file, through the 'Tree' object
	Tree tree;
	tree.ReadFromStream(is,1);
	SetTree(tree);
}

// ---------------------------------------------------------------------------
//		 SetTree
// ---------------------------------------------------------------------------

void PhyloBayes::SetTree(Tree& inTree)	{

	for (int i=0; i<mParam->Nnode;i++)	{
		BL[i] = 0;
	}
	if (mParam->ActivateClock)	{

		/*
		inTree.WriteToStream(cout);
		cout << '\n';
		*/
		inTree.Dichotomise();

		if (! inTree.RegisterWithData(mParam->SpeciesNames, mParam->Ntaxa) )	{
			cerr << "could not register tree from this stream with species from data file " << mParam->DataFileSpec << '\n';
			exit(1);
		}

		if (mParam->OutGroup)	{
			TaxaParameters taxaparam(mParam);
			Bipartition bp(&taxaparam);
			for (int j=0; j<mParam->OutNtaxa; j++)	{
				int k=0;
				while ((k<mParam->Ntaxa) && (mParam->SpeciesNames[k] != mParam->OutGroup[j])) k++;
				if (k == mParam->Ntaxa)	{
					cerr << "error when rerooting : did not find : " << mParam->OutGroup[j] << '\n';
					exit(1);
				}
				bp.mArray[k] = 1;
			}
			inTree.SetParameters(&taxaparam);
			inTree.RootAt(bp);
		}


		// transfer this tree's structure to the present "PhyloBayes" object
		TraversePolyNode(inTree.GetRoot());
		
		if (! mParam->KeepTimes)	{
			for (int j=0; j<mParam->Nnode; j++)	{
				if ((j!= root->label) && (tree[j].branchLength == 0))	{
					tree[j].branchLength = 0.001;
				}
			}
		}
		SetBL();
		// Clockify();
		SetCalibrations();
		if (! mParam->KeepTimes)	{
			MakeClockTree();
		}
		else	{
			SetTimes();
		}

		IsInitBranchLength = Yes;
	}
	else	{
		if (mParam->RandomInitLength)	{
			// erase branch lengths
			for (int i=0; i<mParam->Nnode;i++)	{
					tree[i].branchLength = 0;
			}
		}

		// dichotomise the tree
		inTree.Dichotomise();
		
		// register with data : get the relation between labels and/or species names on the tree
		// and the taxa of the data matrix
		if (! inTree.RegisterWithData(mParam->SpeciesNames, mParam->Ntaxa) )	{
			cerr << "could not register tree from this stream with species from data file " << mParam->DataFileSpec << '\n';
			exit(1);
		}

		if (mParam->OutGroup)	{
			TaxaParameters taxaparam(mParam);
			Bipartition bp(&taxaparam);
			for (int j=0; j<mParam->OutNtaxa; j++)	{
				int k=0;
				while ((k<mParam->Ntaxa) && (mParam->SpeciesNames[k] != mParam->OutGroup[j])) k++;
				if (k == mParam->Ntaxa)	{
					cerr << "error when rerooting : did not find : " << mParam->OutGroup[j] << '\n';
					exit(1);
				}
				bp.mArray[k] = 1;
			}
			inTree.SetParameters(&taxaparam);
			inTree.RootAt(bp);
		}

		// transfer this tree's structure to the present "PhyloBayes" object
		TraversePolyNode(inTree.GetRoot());
		if (mParam->ActivateEM && (! mParam->RandomInitLength))	{
			root->right->branchLength += root->left->branchLength;
			root->left->branchLength = 0;
			for (int i=0; i<mParam->Nnode;i++)	{
				if (blfree(i))	{
					if (tree[i].branchLength <= 0)	{
						cerr << "WARNING : zero branch lengths\n";
						cerr << i << '\t' << tree[i].branchLength << '\n';
						/*
						if (&tree[i] == root->left)	{
							cerr << "root left\n";
						}
						if (&tree[i] == root->right)	{
							cerr << "root right\n";
						}
						if (tree[i].isLeaf())	{
							cerr << mParam->SpeciesNames[i] << '\n';
						}
						*/
					}
				}
			}
			SetBL();
		}
		else	{
			DrawLengths();
		}
		if (mParam->MBL == 2)	{
			SetMBL();
		}
	}
}

void PhyloBayes::SetCalibrations()	{

	if (mParam->NCalib)	{
		cerr << "calibrations\n";
		if (! mParam->isCalibrated)	{
			mParam->isCalibrated = new int[mParam->Nnode];
		}
		for (int j=0; j<mParam->Nnode; j++)	{
			mParam->isCalibrated[j] = 0;
		}
		for (int i=0; i<mParam->NCalib; i++)	{
			int k1 = 0;
			while ((k1<mParam->Ntaxa) && (mParam->CalibTaxon1[i] != mParam->SpeciesNames[k1])) k1++;
			if (k1 == mParam->Ntaxa)	{
				cerr << "error in calibrations: " << mParam->CalibTaxon1[i] << " not among taxa of the dataset\n";
				exit(1);
			}
			int k2 = 0;
			while ((k2<mParam->Ntaxa) && (mParam->CalibTaxon2[i] != mParam->SpeciesNames[k2])) k2++;
			if (k2 == mParam->Ntaxa)	{
				cerr << "error in calibrations: " << mParam->CalibTaxon2[i] << " not among taxa of the dataset\n";
				exit(1);
			}
			int ind1[mParam->Ntaxa-1];
			int ind2[mParam->Ntaxa-1];
			for (int j=0; j<mParam->Ntaxa-1; j++)	{
				ind1[j] = -1;
				ind2[j] = -1;
			}
			int node1 = tree[k1].label;
			int l1= 0;
			while (!tree[node1].isRoot())	{
				node1 = tree[node1].up->label;
				ind1[l1] = node1;
				l1++;
			}
			int node2 = tree[k2].label;
			int l2= 0;
			while (!tree[node2].isRoot())	{
				node2 = tree[node2].up->label;
				ind2[l2] = node2;
				l2++;
			}
			int m1 = 0;
			int m2 = l2;
			while ((m1 < l1) && (m2 == l2))	{
				m2 = 0;
				while ((m2 < l2) && (ind1[m1] != ind2[m2]))	{
					m2++;
				}
				if (m2 == l2)	{
					m1 ++;
				}
			}
			if (m1 == l1)	{
				cerr << "error in set calibrations\n";
				cerr << m1 << '\t' << l1 << '\n';
				exit(1);
			}		
			int found = ind1[m1];
			if ((found != ind2[m2]) || (found == -1))	{
				cerr << "error in set calibrations\n";
				cerr << "found node : " << found << '\n';
				cerr << ind1[m1] << '\t' << ind2[m2] << '\n';
				exit(1);
			}
			cerr << mParam->CalibTaxon1[i] << '\t' << mParam->CalibTaxon2[i] << '\t' << found << '\n';
			mParam->CalibIndex[i] = found;	
			mParam->isCalibrated[found] = 1;

			if ((mParam->CalibUpper[i] != -1) && (mParam->CalibUpper[i] < 0))	{
				cerr << "error in calibrations : does not recognise upper date : " << mParam->CalibUpper[i] << '\n'; 
				cerr << "dates should be positive numbers\n";
				cerr << "except for missing dates (-1)\n";
				exit(1);
			}
			if ((mParam->CalibLower[i] != -1) && (mParam->CalibLower[i] < 0))	{
				cerr << "error in calibrations : does not recognise lower date : " << mParam->CalibLower[i] << '\n'; 
				cerr << "dates should be positive numbers\n";
				cerr << "except for missing dates (-1)\n";
				exit(1);
			}
			if ((mParam->CalibUpper[i] != -1) && (mParam->CalibLower[i] != -1) && (mParam->CalibUpper[i] < mParam->CalibLower[i]))	{
				cerr << "error in calibrations : inconsistent dates " << mParam->CalibUpper[i] << " and " << mParam->CalibLower[i] << '\n';
				cerr << "upper date should be larger than lower date\n";
				exit(1);
			}
		}
	}

	if (mParam->SoftBounds)	{
		// check that the same node is not being calibrated separately
		for (int i=0; i<mParam->NCalib; i++)	{
			for (int j=i+1; j<mParam->NCalib; j++)	{
				if (mParam->CalibIndex[i] == mParam->CalibIndex[j])	{
					cerr << "error : calibration number " << i+1 << " and " << j+1 << " refer to the same node\n";
					cerr << mParam->CalibTaxon1[i] << '\t' << mParam->CalibTaxon2[i] << '\t' << mParam->CalibUpper[i] << '\t' << mParam->CalibLower[i] << '\n';
					cerr << mParam->CalibTaxon1[j] << '\t' << mParam->CalibTaxon2[j] << '\t' << mParam->CalibUpper[j] << '\t' << mParam->CalibLower[j] << '\n';
					cerr << "these 2 calibrations should be merged into one single command\n";
					exit(1);
				}
			}
		}
	}
}

void PhyloBayes::SwapBL()	{
	double* temp = BL;
	BL = RigidBL;
	RigidBL = temp;
	for (int k=0; k<mParam->Ngene; k++)	{
		double* temp = GeneBL[k];
		GeneBL[k] = GeneRigidBL[k];
		GeneRigidBL[k] = temp;
	}
}

void PhyloBayes::SetBL()	{
	for (int i=0; i<mParam->Nnode; i++)	{
		BL[i] = tree[i].branchLength;
	}
}

void PhyloBayes::ReverseSetBL()	{
	for (int i=0; i<mParam->Nnode; i++)	{
		tree[i].branchLength = BL[i];
	}
}

		

// ---------------------------------------------------------------------------
//		 TraversePolyNode
// ---------------------------------------------------------------------------

void PhyloBayes::TraversePolyNode(PolyNode* inNode)	{

	
	if ((inNode->GetLabel() < 0) || (inNode->GetLabel() > mParam->Nnode))	{
		cerr << "error in traverse polynode: tree structure not recognised\n";
		cerr << inNode->GetLabel() << '\n';
		if (inNode->IsRoot())	{
			cerr << "is root\n";
		}
		else if (inNode->up->IsRoot())	{
			cerr << "is son of root\n";
		}
		exit(1);
	}
	tree[inNode->GetLabel()].label = inNode->GetLabel();
	tree[inNode->GetLabel()].branchLength = inNode->GetBranchLength();

	if (inNode->IsRoot())	{
		tree[inNode->GetLabel()].up = 0;
		root = &tree[inNode->GetLabel()];
	}
	else	{
		tree[inNode->GetLabel()].up = &tree[inNode->Up()->GetLabel()];
	}
	if (inNode->IsLeaf())	{
		tree[inNode->GetLabel()].left = 0;
		tree[inNode->GetLabel()].right = 0;
	}
	else	{
		if ((inNode->Down()->GetLabel() < 0) || (inNode->Down()->GetLabel() > mParam->Nnode))	{
			cerr << "error in traverse polynode: tree structure not recognised\n";
			cerr << "down\n";
			cerr << inNode->GetLabel() << '\t' << inNode->Down()->GetLabel() << '\n';
			if (inNode->IsRoot())	{
				cerr << "is root\n";
			}
			else if (inNode->up->IsRoot())	{
				cerr << "is son of root\n";
			}
			exit(1);
		}
		if ((inNode->Down()->Next()->GetLabel() < 0) || (inNode->Down()->Next()->GetLabel() > mParam->Nnode))	{
			cerr << "error in traverse polynode: tree structure not recognised\n";
			cerr << "down next\n";

			cerr << inNode->GetLabel() << '\t' << inNode->Down()->Next()->GetLabel() << '\n';
			if (inNode->IsRoot())	{
				cerr << "is root\n";
			}
			else if (inNode->up->IsRoot())	{
				cerr << "is son of root\n";
			}
			exit(1);
			exit(1);
		}
		tree[inNode->GetLabel()].left = &tree[inNode->Down()->GetLabel()];
		TraversePolyNode(inNode->Down());

		tree[inNode->GetLabel()].right = &tree[inNode->Down()->Next()->GetLabel()];
		TraversePolyNode(inNode->Down()->Next());
	}

}


// ---------------------------------------------------------------------------
//		 Pluck
// ---------------------------------------------------------------------------

void PhyloBayes::GetLeafSet(Node* node, int* taxa)	{

	if (! node->isLeaf())	{
		GetLeafSet(node->left,taxa);
		GetLeafSet(node->right,taxa);
	}
	else	{
		if ((node->label < 0) || (node->label >= mParam->Ntaxa))	{
			cerr << "error in PhyloBayes::GetLeafSet: label: " << node->label << '\n';
		}
		taxa[node->label] = 1;
	}
}

void PhyloBayes::RegisterSubTree(Node* node, int range, int* map, int* reversemap, int& index)	{

	if (! node->isLeaf())	{
		RegisterSubTree(node->left, range, map, reversemap, index);
		RegisterSubTree(node->right, range, map, reversemap, index);
		reversemap[node->label] = index,
		map[index] = node->label;
		index++;
	}
}

void PhyloBayes::Pluck(Node* node, int range, int* map, int* reversemap, int& index)	{

	if (! node->isLeaf())	{
		Pluck(node->left, range, map, reversemap, index);
		Pluck(node->right, range, map, reversemap, index);
	}
	if ((index >= range) || (index < 0))	{
		cerr << "error in PhyloBayes::Pluck: index out of range : " << index << '\n';
		exit(1);
	}
	if ((node->label < 0) || (node->label >= mParam->Nnode-1))	{
		cerr << "error in PhyloBayes::Pluck: node label out of range : " << node->label << '\n';
		exit(1);
	}
	reversemap[node->label] = index,
	map[index] = node->label;
	index++;
}

// ---------------------------------------------------------------------------
//		 SetClockTree
// ---------------------------------------------------------------------------

int PhyloBayes::GetSubTreeSize(Node* node)	{
	if (node->isLeaf())	{
		return 1;
	}
	else{
		return GetSubTreeSize(node->left) + GetSubTreeSize(node->right);
	}
	return 0;
}

void PhyloBayes::SetClockTree()	{

	TraversePolyNode(mParam->TemplateTree->GetRoot());
	// Clockify();
	cerr << "set cal\n";
	SetCalibrations();
	cerr << "ok\n";
	if (! mParam->KeepTimes)	{
		MakeClockTree();
	}
}

void PhyloBayes::MakeClockTree()	{

	double* lowerlimit = new double[mParam->Nnode];
	double* upperlimit = new double[mParam->Nnode];
	double* bklowerlimit = new double[mParam->Nnode];
	double* bkupperlimit = new double[mParam->Nnode];
	for (int j=0; j<mParam->Nnode; j++)	{
		lowerlimit[j] = -1;
		upperlimit[j] = -1;
	}
	double maxupper = -1;
	for (int i=0; i<mParam->NCalib; i++)	{
		int index = mParam->CalibIndex[i];
		if ((index < 0)	|| (index >=mParam->Nnode))	{
			cerr << "in make clock tree : out of bound\n";
			exit(1);
		}
		lowerlimit[index] = mParam->CalibLower[i];
		upperlimit[index] = mParam->CalibUpper[i];
		if (maxupper == -1)	{
			maxupper = mParam->CalibUpper[i];
		}
		else if (mParam->CalibUpper[i] != -1)	{
			if (maxupper < mParam->CalibUpper[i])	{
				maxupper = mParam->CalibUpper[i];
			}
		}
	}
	if (maxupper == -1)	{
		maxupper = 100;
	}
	else	{
		maxupper *=2;
	}
	maxupper = -1;
	PropagateAgeConstraints(root->label,lowerlimit,upperlimit,maxupper);
	for (int j=0; j<mParam->Nnode; j++)	{
		bklowerlimit[j] = lowerlimit[j];
		bkupperlimit[j] = upperlimit[j];
	}
	double bkmaxupper = maxupper;
	double* age = new double[mParam->Nnode];
	int maxntrial = 100;
	int k = 0;
	while (!DrawNodeAge(root->label,lowerlimit,upperlimit,age))	{
		k++;
		for (int j=0; j<mParam->Nnode; j++)	{
			lowerlimit[j] = bklowerlimit[j];
			upperlimit[j] = bkupperlimit[j];
		}
		maxupper = bkmaxupper;
		if (k == maxntrial)	{
			cerr << "repeated error when drawing node ages: please check that imposed calibrations do not imply 0 branch lengths\n";
			exit(1);
		}
	}

	cerr << root->label << '\t';
	cerr << mParam->Nnode << '\n';
	double norm = age[root->label];

	for (int j=0; j<mParam->Nnode; j++)	{
		tree[j].height = age[j] / norm;
	}
	if (mParam->NCalib)	{
		Scale = norm;
	}
	for (int j=0; j<mParam->Nnode; j++)	{
		if (!tree[j].isRoot())	{
			if (tree[j].height == ((AmphiNode*) tree[j].up)->height)	{
				cerr << "error in make clock tree\n";
				cerr << j << '\t' << age[j] << '\t' << age[tree[j].up->label] << '\n';
				exit(1);
			}
		}
	}
	AgeToLength(root->label);

	for (int i=0; i<mParam->Nnode;i++)	{
		if (tree[i].isRoot())	{
			BL[i] = 0;
		}
		else	{
			BL[i] = -MeanLength * log(1 - Random::Uniform());
		}
	}
	for (int k=0; k<mParam->Ngene; k++)	{
		for (int i=0; i<mParam->Nnode;i++)	{
			if (tree[i].isRoot())	{
				GeneBL[k][i] = 0;
			}
			else	{
				GeneBL[k][i] = -MeanLength * log(1 - Random::Uniform());
			}
		}
	}
	delete[] age;
	delete[] lowerlimit;
	delete[] upperlimit;
	delete[] bklowerlimit;
	delete[] bkupperlimit;
	
}

void PhyloBayes::PropagateAgeConstraints(int label, double* lowerlimit, double* upperlimit, double maxupper)	{

	if (! tree[label].isRoot())	{
		int uplabel = tree[label].up->label;
		if (upperlimit[label] == -1)	{
			upperlimit[label] = upperlimit[uplabel];
		}
		else if (upperlimit[uplabel] != -1)	{
			if (upperlimit[label] > upperlimit[uplabel])	{
				upperlimit[label] = upperlimit[uplabel];
			}
		}
	}
	if (! tree[label].isLeaf())	{
		PropagateAgeConstraints(tree[label].left->label,lowerlimit,upperlimit,maxupper);
		PropagateAgeConstraints(tree[label].right->label,lowerlimit,upperlimit,maxupper);
		
		int leftlabel = tree[label].left->label;
		int rightlabel = tree[label].right->label;
		double limit = lowerlimit[label];
		if (limit == -1)	{
			limit = lowerlimit[leftlabel];
		}
		else if (lowerlimit[leftlabel] != -1)	{
			if (limit < lowerlimit[leftlabel])	{
				limit = lowerlimit[leftlabel];
			}
		}
		if (limit == -1)	{
			limit = lowerlimit[rightlabel];
		}
		else if (lowerlimit[rightlabel] != -1)	{
			if (limit < lowerlimit[rightlabel])	{
				limit = lowerlimit[rightlabel];
			}
		}
		lowerlimit[label] = limit;
	}
	else	{
		lowerlimit[label] = 0;
	}
	/*
	if (upperlimit[label] == -1)	{
		upperlimit[label] = maxupper;
	}
	*/
}

/*
int PhyloBayes::DrawNodeAge(int label, double* lowerlimit, double* upperlimit, double* age)	{

	int succeed = 1;
	if (! tree[label].isLeaf())	{
		if (upperlimit[label] == -1)	{
			cerr << "error : upper limit\n";
			exit(1);
		}
		if (lowerlimit[label] == -1)	{
			cerr << "error : lower limit\n";
			exit(1);
		}
		if (upperlimit[label] == lowerlimit[label])	{
			cerr << "error : upper limit is equal to lower limit\n";
			exit(1);
		}
		if (upperlimit[label] < lowerlimit[label])	{
			cerr << "error : upper limit is less than lower limit\n";
			exit(1);
		}
		if (! tree[label].isRoot())	{
			double upage = age[tree[label].up->label];
			if (upperlimit[label] > upage)	{
				if (lowerlimit[label] > upage)	{
					cerr << "error : lowerlimit older than parent node\n";
					exit(1);
				}
				upperlimit[label] = upage;
			}
		}
		age[label] = Random::Uniform() * (upperlimit[label] - lowerlimit[label]) + lowerlimit[label];	
		if (age[label] == upperlimit[label])	{
			cerr << "age is equal to upperlimit\n";
			exit(1);
		}
		if (age[label] == lowerlimit[label])	{
			cerr << "age is equal to lowerlimit\n";
			exit(1);
		}
		if (! tree[label].isRoot())	{
			if (fabs(age[label] - age[tree[label].up->label]) < 1e-4)	{
				succeed = 0;
			}
		}

		succeed &= DrawNodeAge(tree[label].left->label,lowerlimit,upperlimit,age);
		succeed &= DrawNodeAge(tree[label].right->label,lowerlimit,upperlimit,age);
	}
	else	{
		age[label] = 0;
	}
	return succeed;
}
*/

int PhyloBayes::DrawNodeAge(int label, double* lowerlimit, double* upperlimit, double* age)	{

	int succeed = 1;
	if ((label < 0) || (label >= mParam->Nnode))	{
		cerr << "error in draw node age\n";
		exit(1);
	}
	if (! tree[label].isLeaf())	{
		succeed &= DrawNodeAge(tree[label].left->label,lowerlimit,upperlimit,age);
		succeed &= DrawNodeAge(tree[label].right->label,lowerlimit,upperlimit,age);

		if (lowerlimit[label] == -1)	{
			cerr << "error : lower limit\n";
			exit(1);
		}
		/*
		if (upperlimit[label] == -1)	{
			cerr << "error : upper limit\n";
			exit(1);
		}
		*/
		if (upperlimit[label] == lowerlimit[label])	{
			cerr << "error : upper limit is equal to lower limit\n";
			exit(1);
		}
		if ((upperlimit[label] != -1) && (upperlimit[label] < lowerlimit[label]))	{
			cerr << "error : upper limit is less than lower limit\n";
			exit(1);
		}
		
		double leftage = age[tree[label].left->label];
		double rightage = age[tree[label].right->label];
		double lowerage = leftage;
		if (lowerage < rightage)	{
			lowerage = rightage;
		}
		if (lowerlimit[label]<lowerage)	{
			if ((upperlimit[label] != -1) && (upperlimit[label] < lowerage))	{
				cerr << "error : upperlimit older than child node\n";
				exit(1);
			}
			lowerlimit[label] = lowerage;
		}
		if (upperlimit[label] == -1)	{
			age[label] = lowerlimit[label] + 10;
		}
		else	{
			double t = Random::Uniform();
			// age[label] = t * (upperlimit[label] - lowerlimit[label]) + lowerlimit[label];	
			age[label] = 0.5 * t * (upperlimit[label] - lowerlimit[label]) + lowerlimit[label];	
		}
		if (age[label] == upperlimit[label])	{
			cerr.precision(20);
			cerr << "age is equal to upperlimit\n";
			cerr << upperlimit[label] << '\t' << lowerlimit[label] << '\t' << age[label] << '\n';
			exit(1);
		}
		if (age[label] == lowerlimit[label])	{
			cerr.precision(20);
			cerr << "age is equal to lowerlimit\n";
			cerr << upperlimit[label] << '\t' << lowerlimit[label] << '\t' << age[label] << '\n';
			exit(1);
		}
		if (! tree[label].isLeaf())	{
			if (fabs(age[tree[label].left->label] - age[label]) < 1e-4)	{
				succeed = 0;
			}
			if (fabs(age[tree[label].right->label] - age[label]) < 1e-4)	{
				succeed = 0;
			}
		}
		/*
		if (! tree[label].isRoot())	{
			if (fabs(age[label] - age[tree[label].up->label]) < 1e-4)	{
				succeed = 0;
			}
		}
		*/
	}
	else	{
		age[label] = 0;
	}
	return succeed;
}

void PhyloBayes::AgeToLength(int label)	{
	
	if (tree[label].isRoot())	{
		tree[label].branchLength = 0;
	}
	else	{
		tree[label].branchLength = ((AmphiNode*) (tree[label].up))->height - ((AmphiNode) (tree[label])).height;
		if (tree[label].branchLength == 0)	{
			cerr << "in AgeToLength : 0 bl\n";
			cerr << "node : " << label << '\n';
			if (tree[label].isLeaf())	{
				cerr << " leaf node\n";
			}
			else	{
				cerr << " internal node\n";
			}
			exit(1);
		}
	}
	if (!tree[label].isLeaf())	{
		AgeToLength(tree[label].left->label);
		AgeToLength(tree[label].right->label);
	}
}

void PhyloBayes::Clockify()	{

	// clockify the tree
    // deprecated
	cerr << "clockify\n";
	exit(1);
	
	for (int j=0; j<mParam->Nnode; j++)	{
		if (!tree[j].isRoot())	{
			if (! mParam->KeepTimes)	{
				if (tree[j].branchLength < 0.001)	{
					tree[j].branchLength = 0.001;
				}
			}
		}
	}
	AmphiNode* first = root->order();
	AmphiNode* current = first;
	double maxheight = first->height;
	do	{
		if (current->height > maxheight)	{
			maxheight = current->height;
		}
		current = current->next->next;
	} while (current != first->next);
	current = first;
	do	{
		current->branchLength += maxheight - current->height;
		current = current->next->next;
	} while (current != first->next);
	current = first;
	do	{
		current->branchLength /= maxheight;
		current = current->next;
	}	
 	while (current != first);
	first = root->order();
	current = first;
	do	{
		current = current->next;
	}	
 	while (current != first);
	
	root->branchLength = 0;

	for (int i=0; i<mParam->Nnode;i++)	{
		if (tree[i].isRoot())	{
			BL[i] = 0;
		}
		else	{
			BL[i] = -MeanLength * log(1 - Random::Uniform());
		}
	}
	for (int k=0; k<mParam->Ngene; k++)	{
		for (int i=0; i<mParam->Nnode;i++)	{
			if (tree[i].isRoot())	{
				GeneBL[k][i] = 0;
			}
			else	{
				GeneBL[k][i] = -MeanLength * log(1 - Random::Uniform());
			}
		}
	}
}	

// ---------------------------------------------------------------------------
//		 FindRoot
// ---------------------------------------------------------------------------

AmphiNode* PhyloBayes::FindRoot()	{

	int i=0;
	while ( (i < mParam->Nnode) && (! tree[i].isRoot()) )	{
		i++;
	}
	if ( i == mParam->Nnode)	{
		cerr << "error : could not find root\n";
		throw;
	}
	return &tree[i];
}


// ***************************************************************************
// ***************************************************************************
//		 Streams and Initialisations
// ***************************************************************************
// ***************************************************************************


// ---------------------------------------------------------------------------
//		 ostream 
// ---------------------------------------------------------------------------

ostream& operator<<(ostream& os, PhyloBayes& state)	{

	// for both pruning-based and normal approx

	os << state.Version << '\t';
	os << state.SubVersion << '\t';
	os << state.Beta << '\n';

	// tree
	AmphiNode* p =state.tree;
	for (int i=0; i< state.mParam->Nnode; i++)	{
		os << ((p->up) ? p->up->label : -1) << '\n';
		os << ((p->left) ? p->left->label : -1) << '\n';
		os << ((p->right) ? p->right->label  : -1) << '\n';
		os << p->branchLength << '\n';
		p++;
	}
	os << state.root->label << '\n';
	for (int j=0; j<state.mParam->Nnode; j++)	{
		os << state.BL[j] << '\t';
		os << state.RigidBL[j] << '\n';
	}
	if (state.mParam->MBL)	{
		os << state.MixBLAlpha << '\n';
		os << state.NMixBL << '\n';
		for (int i=0; i<state.NMixBL; i++)	{
			for (int j=0; j<state.mParam->Nnode; j++)	{
				os << state.MixBL[i][j] << '\t';
			}
			os << '\n';
		}
		for (int i=0; i<state.mParam->Nsite; i++)	{
			os << state.MixBLMode[i] << '\t';
		}
		os << '\n';
	}

	if (!state.mParam->SaveAll)	{
		os << state.MeanLength << '\n';
		os << state.gamma << '\n';
		os << state.pconst << '\n';
		os << state.alpha << '\n';
		os << state.Nmode << '\n';
		os << state.xi0 << '\n';
		os << state.invprob << '\n';
		os << state.Sigma << '\n';
		os << state.Theta << '\n';
		os << state.Mu << '\n';
		os << state.Scale << '\n';
		os << state.mLogPrior << '\n';
		os << state.mLogSampling << '\n';
	}
	else	{
			
		for (int i=0; i<state.mParam->Ngene; i++)	{
			for (int j=0; j<state.mParam->Nnode; j++)	{
				os << state.GeneBL[i][j] << '\t';
				os << state.GeneRigidBL[i][j] << '\n';
			}
		}

		os << state.MeanLength << '\n';
		os << state.VarLength << '\n';
		os << state.LengthGamma << '\n';

		if (state.mParam->NormalApprox == No)	{
			// rate
			for (int i=0; i<state.mParam->Nsite; i++)	{
				os << state.rate[i] << '\t';
			}
			os << '\n';
			os << state.gamma << '\n';
			os << state.pconst << '\n';

			// gene
			for (int i=0; i<state.mParam->Ngene; i++)	{
				os << state.GeneGamma[i] << '\t';
				os << state.GeneRate[i] << '\t';
				for (int j=0; j<state.Nstate; j++)	{
					os << state.GeneStationary[i][j] << '\t';
				}
				os << '\n';
			}

			// ref
			for (int i=0; i<state.Nstate; i++)	{
				os << state.RefStationary[i] << '\t';
			}
			os << '\n';
			for (int i=0; i<state.Nstate*(state.Nstate-1)/2; i++)	{
				os << state.RefRR[i] << '\t';
			}
			os << '\n';

			// mode
			os << state.alpha << '\n';
			os << state.Nmode << '\n';
			for (int i=0; i<state.Nmode; i++)	{
				if (state.mParam->SaveStat)	{
					for (int j=0; j<state.Nstate; j++)	{
						os << state.Stationary[i][j] << '\t';
					}
				}
				os << '\n';
				if (state.mParam->Qmode)	{
					for (int j=0; j<state.Nstate*(state.Nstate-1)/2; j++)	{
						os << state.RR[i][j] << '\t';
					}
				}			
				os << '\n';
				os << state.ModeWeight[i] << '\n';
			}
			if (! state.mParam->SaveStat)	{
				os << state.GetStatEnt() << '\n';
			}
			for (int i=0; i<state.Nstate*(state.Nstate-1)/2; i++)	{
				os << state.ModeRR[i] << '\t';
			}
			os << '\n';

			for (int i=0; i<state.Nstate; i++)	{
				os << state.ModeStatCenter[i] << '\t';
			}
			os << '\n';
			os << state.ModeStatAlpha << '\n';

			for (int i=0; i<state.mParam->Nsite; i++)	{
				os << state.Mode[i] << '\t';
			}
			os << '\n';

			// ratemode
			os << state.NRateMode << '\n';
			for (int i=0; i<state.NRateMode; i++)	{
				os << state.ModeRate[i] << '\t';
				os << state.RateModeWeight[i] << '\n';
			}
			for (int i=0; i<state.mParam->Nsite; i++)	{
				os << state.RateMode[i] << '\t';
			}
			os << '\n';
			os << state.RateAlpha << '\n';

			// log probs
			if (! isnan(state.mLogPrior))	{
				os << state.mLogPrior << '\n';
			}
			else 	{
				if (state.mParam->Beta0 != 0)	{
					cerr << "error : log prior is nan\n";
				}
				os << 0 << '\n';
			}
			os << state.mLogSampling << '\n';
			if ( state.mParam->SavePartialLogLikelihoods)	{
				for (int i=0; i<state.mParam->Nsite; i++)	{
					os << state.mSiteLogSampling[i] << '\t';
				}
				os << '\n';
			}
		
			// hetero
			os << state.xi0 << '\n';
			os << state.invprob << '\n';
		}
		else	{	// only for normal approx

			if (! isnan(state.mLogPrior))	{
				os << state.mLogPrior << '\n';
			}
			else 	{
				if (state.mParam->Beta0 != 0)	{
					cerr << "error : log prior is nan\n";
				}
				os << 0 << '\n';
			}
			os << state.mLogSampling << '\n';

		}
		if (state.mParam->ActivateClock)	{
			os << state.Sigma << '\n';
			os << state.Theta << '\n';
			os << state.Chi << '\n';
			os << state.Chi2 << '\n';
			/*
			os << state.MulSigma << '\n';
			os << state.MulTheta << '\n';
			*/
			os << state.MeanSigma << '\n';
			os << state.MeanTheta << '\n';
			os << state.VarSigma << '\n';
			os << state.VarTheta << '\n';
			os << state.Mu << '\n';
			os << state.Scale << '\n';
			os << state.VarScale << '\n';
			if (state.mParam->ContNchar)	{
				for (int k=0; k<state.mParam->ContNchar+1; k++)	{
					for (int j=0; j<state.mParam->Nnode; j++)	{
						os << state.ContRho[k][j] << '\t';
					}
					os << state.ContRhoCenter[k] << '\n';
				}
			}
			else	{
				for (int j=0; j<state.mParam->Nnode; j++)	{
					os << state.Rho[j] << '\t';
				}
			}
			os << '\n';
			for (int k=0; k<state.mParam->Ngene; k++)	{
				for (int j=0; j<state.mParam->Nnode; j++)	{
					os << state.GeneRho[k][j] << '\t';
				}
				os << '\n';
				os << state.GeneSigma[k] << '\n';
				os << state.GeneTheta[k] << '\n';
				os << state.GeneMu[k] << '\n';
			}
		}
	}

	if ((state.mParam->ZipGTR) && (state.mParam->ZipGTR != 4))	{
		os << state.ZipAlpha << '\n';
		os << state.ZipRate << '\n';
		os << state.NZipMode << '\n';
		for (int i=0; i<state.mParam->Nsite; i++)	{
			os <<  state.ZipMode[i] << '\t';
		}
		os << '\n';
		for (int k=0; k<state.Nstate; k++)	{
			os << state.ModeZipProb[k] << '\t';
			os << state.ZipA0[k] << '\t';
			os << state.ZipA1[k] << '\n';
		}
		if (state.mParam->ZipGTR == 1)	{
			for (int i=0; i<state.Nmode; i++)	{
				for (int k=0; k<state.Nstate; k++)	{
					os <<  state.ModeAASupport[i][k] << '\t';
				}
				os << '\n';
			}
		}
		else if ((state.mParam->ZipGTR >= 2) && (state.mParam->ZipGTRDP))	{
			for (int i=0; i<state.NZipMode; i++)	{
				for (int k=0; k<state.Nstate; k++)	{
					os << state.ModeAASupport[i][k] << '\t';
				}
				os << '\n';
			}
		}
		else	{
			for (int i=0; i<state.mParam->Nsite; i++)	{
				for (int k=0; k<state.Nstate; k++)	{
					os << state.SiteAASupport[i][k] << '\t';
				}
				os << '\n';
			}
		}
	}


	if (state.mParam->NH)	{

		os << state.NHNcat << '\n';
		for (int j=0; j<state.mParam->Nnode; j++)	{
			os << state.NHcat[j] << '\t';
		}
		os << '\n';
		for (int n=0; n<state.NHNcat; n++)	{
			os << state.NHWeight[n] << '\n';
			for (int k=0; k<state.mParam->Nstate; k++)	{
				os << state.NHStatDistorter[n][k] << '\t';
			}
			if (state.mParam->NHPrior)	{
				for (int k=0; k<state.ContNstate; k++)	{
					os << state.logNHStatDistorter[n][k] << '\t';
				}
			}
			os << '\n';
		}
		os << state.NHpswitch << '\n';
		os << state.NHStatAlpha << '\n';
		for (int k=0; k<state.ContNstate; k++)	{
			os << state.NHStatCenter[k] << '\t';
		}
		os << '\n';
		if (state.mParam->NHPrior)	{
			for (int k=0; k<state.mParam->ContNchar+1; k++)	{
				os << state.NHContVar[k] << '\t';
			}
			os << '\n';
			for (int k=0; k<state.ContNstate; k++)	{
				os << state.NHVar[k] << '\t';
			}
			os << '\n';
			if (state.mParam->NHPrior == 2)	{
				for (int k=0; k<state.ContNrr; k++)	{
					os << state.NHCovIndex[k] << '\t';
				}
				os << '\n';
			}
		}	
	}

	if (state.mParam->Qmode)	{
		os << state.RRAlpha << '\n';
	}

	os << state.GWf << '\n';
	if (state.mParam->PopEff)	{
		os << state.PopAlpha << '\n';
		os << state.PopBeta << '\n';
		for (int j=0; j<state.NHNcat; j++)	{
			os << state.PopEffSize[j] << '\n';
		}
	}
	if (state.mParam->MutMode == 2)	{
		os << state.Kappa << '\n';
		for (int i=0; i<Nnuc; i++)	{
			os << state.RefACGT[i] << '\t';
		}
		os << '\n';
	}
	if (state.mParam->MutMode == 3)	{
		for (int i=0; i<Nnuc; i++)	{
			os << state.RefACGT[i] << '\t';
		}
		os << '\n';
		for (int i=0; i<Nnuc * (Nnuc-1) / 2; i++)	{
			os << state.RefRRACGT[i] << '\t';
		}
		os << '\n';
	}

	return os;
}


// ---------------------------------------------------------------------------
//		 istream 
// ---------------------------------------------------------------------------

istream& operator>>(istream& is, PhyloBayes& state)	{

	is >> state.Version;
	is >> state.SubVersion;
	
	if (state.Version < 2)	{
		cerr << "error: version 1 is not supported by current version of phylobayes\n";
		cerr << "sorry\n";
		exit(1);
	}

	if ((state.Version == 3) && (state.SubVersion == 0))	{
		is >> state.Beta;

		// tree
		AmphiNode* p =state.tree;
		for (int i=0; i< state.mParam->Nnode; i++)	{
			int s;
			is >> s;
			p->up = (s == -1) ? 0 : &state.tree[s];
			is >> s;
			p->left = (s == -1) ? 0 : &state.tree[s];
			is >> s;
			p->right = (s== -1) ? 0 : &state.tree[s];
			is >> p->branchLength;
			p++;
		}
		int rootlabel;
		is >> rootlabel;
		state.root = & state.tree[rootlabel];

		for (int j=0; j<state.mParam->Nnode; j++)	{
			is >> state.BL[j];
			if (state.BL[j] < 1e-8)	{
				state.BL[j] = 1e-8;
			}
			is >> state.RigidBL[j];
			if (state.RigidBL[j] < 1e-8)	{
				state.RigidBL[j] = 1e-8;
			}
		}
		if (state.mParam->MBL)	{
			is >> state.MixBLAlpha;
			is >> state.NMixBL;
			for (int i=0; i<state.NMixBL; i++)	{
				for (int j=0; j<state.mParam->Nnode; j++)	{
					is >>state.MixBL[i][j];
				}
				state.MixBLSiteNumber[i] = 0;
			}
			for (int i=0; i<state.mParam->Nsite; i++)	{
				is >> state.MixBLMode[i];
			}
		}
		if (!state.mParam->SaveAll)	{
			is >> state.MeanLength;
			is >> state.gamma;
			is >> state.pconst;
			is >> state.alpha;
			is >> state.Nmode;
			is >> state.xi0;
			is >> state.invprob;
			is >> state.Sigma;
			is >> state.Theta;
			is >> state.Mu;
			is >> state.Scale;
			is >> state.mLogPrior;
			is >> state.mLogSampling;
		}
		else	{
			for (int i=0; i<state.mParam->Ngene; i++)	{
				for (int j=0; j<state.mParam->Nnode; j++)	{
					is >> state.GeneBL[i][j];
					is >> state.GeneRigidBL[i][j];
				}
			}

			is >> state.MeanLength;
			is >> state.VarLength;
			is >> state.LengthGamma;

			if (state.mParam->NormalApprox == No)	{

				// rate
				for (int i=0; i<state.mParam->Nsite; i++)	{
					is >> state.rate[i];
				}
				is >> state.gamma;
				is >> state.pconst;

				// gene
				for (int i=0; i<state.mParam->Ngene; i++)	{
					is >> state.GeneGamma[i];
					is >> state.GeneRate[i];
					for (int j=0; j<state.Nstate; j++)	{
						is >> state.GeneStationary[i][j];
					}
				}

				// ref
				for (int i=0; i<state.Nstate; i++)	{
					is >> state.RefStationary[i];
				}
				for (int i=0; i<state.Nstate*(state.Nstate-1)/2; i++)	{
					is >> state.RefRR[i];
				}

				// mode
				is >> state.alpha;
				is >> state.Nmode;
				for (int i=0; i<state.Nmode; i++)	{
					if (state.mParam->SaveStat)	{
						double total = 0;
						for (int j=0; j<state.Nstate; j++)	{
							is >> state.Stationary[i][j];
							total += state.Stationary[i][j];
						}
						for (int j=0; j<state.Nstate; j++)	{
							state.Stationary[i][j] /= total;
						}
					}
					if (state.mParam->Qmode)	{
						for (int j=0; j<state.Nstate*(state.Nstate-1)/2; j++)	{
							is >> state.RR[i][j];
						}
					}			
					is >> state.ModeWeight[i];
				}
				if (! state.mParam->SaveStat)	{
					is >> state.mStationaryEntropy;
				}

				for (int i=0; i<state.Nstate*(state.Nstate-1)/2; i++)	{
					is >> state.ModeRR[i];
				}

				for (int i=0; i<state.Nstate; i++)	{
					is >> state.ModeStatCenter[i];
				}
				is >> state.ModeStatAlpha;

				for (int i=0; i<state.mParam->Nsite; i++)	{
					is >> state.Mode[i];
				}
				// ratemode
				is >> state.NRateMode;
				for (int i=0; i<state.NRateMode; i++)	{
					is >> state.ModeRate[i];
					if (state.ModeRate[i] == 0) {
						state.ModeRate[i] = 1e-8;
					}
					is >> state.RateModeWeight[i];
				}
				for (int i=0; i<state.mParam->Nsite; i++)	{
					is >> state.RateMode[i];
				}
				is >> state.RateAlpha;

				// log probs
				is >> state.mLogPrior;
				is >> state.mLogSampling;
				state.mLogPosterior = state.mLogPrior + state.Beta * state.mLogSampling;

				if ( state.mParam->SavePartialLogLikelihoods)	{
					double temp = 0;
					for (int i=0; i<state.mParam->Nsite; i++)	{
						is >> state.mSiteLogSampling[i];
						temp += state.mSiteLogSampling[i];
					}
					if (temp != state.mLogSampling)	{
						cerr << "warning : partial likelihoods non concordant with total logsampling advertised in this file\n";
					}
				}
			
				// hetero
				is >> state.xi0;
				is >> state.invprob;

			}
			else	{ // normal approx

				is >> state.mLogPrior;
				is >> state.mLogSampling;
				state.mLogPosterior = state.mLogPrior + state.Beta *  state.mLogSampling;

			}
			if (state.mParam->ActivateClock)	{
				is >> state.Sigma;
				is >> state.Theta;
				is >> state.MulSigma;
				is >> state.MulTheta;
				is >> state.MeanSigma;
				is >> state.MeanTheta;
				is >> state.VarSigma;
				is >> state.VarTheta;
				is >> state.Mu;
				is >> state.Scale;
				is >> state.VarScale;
				if (state.mParam->ContNchar)	{
					for (int k=0; k<state.mParam->ContNchar+1; k++)	{
						for (int j=0; j<state.mParam->Nnode; j++)	{
							is >> state.ContRho[k][j];
						}
						is >> state.ContRhoCenter[k];
					}
				}
				else	{
					for (int j=0; j<state.mParam->Nnode; j++)	{
						is >> state.Rho[j];
					}
				}
				for (int k=0; k<state.mParam->Ngene; k++)	{
					for (int j=0; j<state.mParam->Nnode; j++)	{
						is >> state.GeneRho[k][j];
					}
					is >> state.GeneSigma[k];
					is >> state.GeneTheta[k];
					is >> state.GeneMu[k];
				}
			}
		}

		if (state.mParam->ZipGTR)	{
			is >> state.ZipAlpha;
			is >> state.ZipRate;
			is >> state.NZipMode;
			for (int i=0; i<state.mParam->Nsite; i++)	{
				is >> state.ZipMode[i];
			}
			for (int k=0; k<state.Nstate; k++)	{
				is >> state.ModeZipProb[k];
				is >> state.ZipA0[k];
				is >> state.ZipA1[k];
			}
			if (state.mParam->ZipGTR == 1)	{
				for (int i=0; i<state.Nmode; i++)	{
					for (int k=0; k<state.Nstate; k++)	{
						is >> state.ModeAASupport[i][k];
					}
				}
			}
			else if ((state.mParam->ZipGTR >= 2) && (state.mParam->ZipGTRDP))	{
				for (int i=0; i<state.NZipMode; i++)	{
					for (int k=0; k<state.Nstate; k++)	{
						is >> state.ModeAASupport[i][k];
					}
				}
			}
			else	{
				for (int i=0; i<state.mParam->Nsite; i++)	{
					for (int k=0; k<state.Nstate; k++)	{
						is >> state.SiteAASupport[i][k];
					}
				}
			}
		}

		if (state.mParam->NH)	{
			is >> state.NHNcat;
			for (int j=0; j<state.mParam->Nnode; j++)	{
				is >> state.NHcat[j];
			}
			for (int n=0; n<state.NHNcat; n++)	{
				is >> state.NHWeight[n];
				for (int k=0; k<state.mParam->Nstate; k++)	{
					is >>state.NHStatDistorter[n][k];
				}
				if (state.mParam->NHPrior)	{
					for (int k=0; k<state.ContNstate; k++)	{
						is >>state.logNHStatDistorter[n][k];
					}
				}
			}
			is >> state.NHpswitch;
			is >> state.NHStatAlpha;
			for (int k=0; k<state.ContNstate; k++)	{
				is >> state.NHStatCenter[k];
			}
			if (state.mParam->NHPrior)	{
				for (int k=0; k<state.mParam->ContNchar+1; k++)	{
					is >> state.NHContVar[k];
				}
				for (int k=0; k<state.ContNstate; k++)	{
					is >> state.NHVar[k];
				}
				if (state.mParam->NHPrior == 2)	{
					for (int k=0; k<state.ContNrr; k++)	{
						is >> state.NHCovIndex[k];
					}
				}
			}	
		}
		if (state.mParam->Qmode)	{
			is >> state.RRAlpha;
		}

		is >> state.GWf;

		if (state.mParam->PopEff)	{
			is >> state.PopAlpha;
			for (int j=0; j<state.NHNcat; j++)	{
				is >> state.PopEffSize[j];
			}
		}
		if (state.mParam->MutMode == 2)	{
			is >> state.Kappa;
			for (int i=0; i<Nnuc; i++)	{
				is >> state.RefACGT[i];
			}
		}
		if (state.mParam->MutMode == 3)	{
			for (int i=0; i<Nnuc; i++)	{
				is >> state.RefACGT[i];
			}
			for (int i=0; i<Nnuc * (Nnuc-1) / 2; i++)	{
				is >> state.RefRRACGT[i];
			}
		}

		return is;
	}
	else if ((state.Version == 3) && (state.SubVersion == 1))	{

		is >> state.Beta;

		// tree
		AmphiNode* p =state.tree;
		for (int i=0; i< state.mParam->Nnode; i++)	{
			int s;
			is >> s;
			p->up = (s == -1) ? 0 : &state.tree[s];
			is >> s;
			p->left = (s == -1) ? 0 : &state.tree[s];
			is >> s;
			p->right = (s== -1) ? 0 : &state.tree[s];
			is >> p->branchLength;
			p++;
		}
		int rootlabel;
		is >> rootlabel;
		state.root = & state.tree[rootlabel];

		for (int j=0; j<state.mParam->Nnode; j++)	{
			is >> state.BL[j];
			is >> state.RigidBL[j];
		}
		if (state.mParam->MBL)	{
			is >> state.MixBLAlpha;
			is >> state.NMixBL;
			for (int i=0; i<state.NMixBL; i++)	{
				for (int j=0; j<state.mParam->Nnode; j++)	{
					is >>state.MixBL[i][j];
				}
				state.MixBLSiteNumber[i] = 0;
			}
			for (int i=0; i<state.mParam->Nsite; i++)	{
				is >> state.MixBLMode[i];
			}
		}
		if (!state.mParam->SaveAll)	{
			is >> state.MeanLength;
			is >> state.gamma;
			is >> state.pconst;
			is >> state.alpha;
			is >> state.Nmode;
			is >> state.xi0;
			is >> state.invprob;
			is >> state.Sigma;
			is >> state.Theta;
			is >> state.Mu;
			is >> state.Scale;
			is >> state.mLogPrior;
			is >> state.mLogSampling;
		}
		else	{
			for (int i=0; i<state.mParam->Ngene; i++)	{
				for (int j=0; j<state.mParam->Nnode; j++)	{
					is >> state.GeneBL[i][j];
					is >> state.GeneRigidBL[i][j];
				}
			}

			is >> state.MeanLength;
			is >> state.VarLength;
			is >> state.LengthGamma;

			if (state.mParam->NormalApprox == No)	{

				// rate
				for (int i=0; i<state.mParam->Nsite; i++)	{
					is >> state.rate[i];
				}
				is >> state.gamma;
				is >> state.pconst;

				// gene
				for (int i=0; i<state.mParam->Ngene; i++)	{
					is >> state.GeneGamma[i];
					is >> state.GeneRate[i];
					for (int j=0; j<state.Nstate; j++)	{
						is >> state.GeneStationary[i][j];
					}
				}

				// ref
				for (int i=0; i<state.Nstate; i++)	{
					is >> state.RefStationary[i];
				}
				for (int i=0; i<state.Nstate*(state.Nstate-1)/2; i++)	{
					is >> state.RefRR[i];
				}

				// mode
				is >> state.alpha;
				is >> state.Nmode;
				for (int i=0; i<state.Nmode; i++)	{
					if (state.mParam->SaveStat)	{
						double total = 0;
						for (int j=0; j<state.Nstate; j++)	{
							is >> state.Stationary[i][j];
							total += state.Stationary[i][j];
						}
						for (int j=0; j<state.Nstate; j++)	{
							state.Stationary[i][j] /= total;
						}
					}
					if (state.mParam->Qmode)	{
						for (int j=0; j<state.Nstate*(state.Nstate-1)/2; j++)	{
							is >> state.RR[i][j];
						}
					}			
					is >> state.ModeWeight[i];
				}
				if (! state.mParam->SaveStat)	{
					is >> state.mStationaryEntropy;
				}

				for (int i=0; i<state.Nstate*(state.Nstate-1)/2; i++)	{
					is >> state.ModeRR[i];
				}

				for (int i=0; i<state.Nstate; i++)	{
					is >> state.ModeStatCenter[i];
				}
				is >> state.ModeStatAlpha;

				for (int i=0; i<state.mParam->Nsite; i++)	{
					is >> state.Mode[i];
				}
				// ratemode
				is >> state.NRateMode;
				for (int i=0; i<state.NRateMode; i++)	{
					is >> state.ModeRate[i];
					is >> state.RateModeWeight[i];
				}
				for (int i=0; i<state.mParam->Nsite; i++)	{
					is >> state.RateMode[i];
				}
				is >> state.RateAlpha;

				// log probs
				is >> state.mLogPrior;
				is >> state.mLogSampling;
				state.mLogPosterior = state.mLogPrior + state.Beta * state.mLogSampling;

				if ( state.mParam->SavePartialLogLikelihoods)	{
					double temp = 0;
					for (int i=0; i<state.mParam->Nsite; i++)	{
						is >> state.mSiteLogSampling[i];
						temp += state.mSiteLogSampling[i];
					}
					if (temp != state.mLogSampling)	{
						cerr << "warning : partial likelihoods non concordant with total logsampling advertised in this file\n";
					}
				}
			
				// hetero
				is >> state.xi0;
				is >> state.invprob;

			}
			else	{ // normal approx

				is >> state.mLogPrior;
				is >> state.mLogSampling;
				state.mLogPosterior = state.mLogPrior + state.Beta *  state.mLogSampling;

			}
			if (state.mParam->ActivateClock)	{
				is >> state.Sigma;
				is >> state.Theta;
				is >> state.Chi;
				is >> state.Chi2;
				/*	
				is >> state.MulSigma;
				is >> state.MulTheta;
				*/
				is >> state.MeanSigma;
				is >> state.MeanTheta;
				is >> state.VarSigma;
				is >> state.VarTheta;
				is >> state.Mu;
				is >> state.Scale;
				is >> state.VarScale;
				if (state.mParam->ContNchar)	{
					for (int k=0; k<state.mParam->ContNchar+1; k++)	{
						for (int j=0; j<state.mParam->Nnode; j++)	{
							is >> state.ContRho[k][j];
						}
						is >> state.ContRhoCenter[k];
					}
				}
				else	{
					for (int j=0; j<state.mParam->Nnode; j++)	{
						is >> state.Rho[j];
					}
				}
				for (int k=0; k<state.mParam->Ngene; k++)	{
					for (int j=0; j<state.mParam->Nnode; j++)	{
						is >> state.GeneRho[k][j];
					}
					is >> state.GeneSigma[k];
					is >> state.GeneTheta[k];
					is >> state.GeneMu[k];
				}
			}
		}

		if ((state.mParam->ZipGTR) && (state.mParam->ZipGTR != 4))	{
			is >> state.ZipAlpha;
			is >> state.ZipRate;
			is >> state.NZipMode;
			for (int i=0; i<state.mParam->Nsite; i++)	{
				is >> state.ZipMode[i];
			}
			for (int k=0; k<state.Nstate; k++)	{
				is >> state.ModeZipProb[k];
				is >> state.ZipA0[k];
				is >> state.ZipA1[k];
			}
			if (state.mParam->ZipGTR == 1)	{
				for (int i=0; i<state.Nmode; i++)	{
					for (int k=0; k<state.Nstate; k++)	{
						is >> state.ModeAASupport[i][k];
					}
				}
			}
			else if ((state.mParam->ZipGTR >= 2) && (state.mParam->ZipGTRDP))	{
				for (int i=0; i<state.NZipMode; i++)	{
					for (int k=0; k<state.Nstate; k++)	{
						is >> state.ModeAASupport[i][k];
					}
				}
			}
			else	{
				for (int i=0; i<state.mParam->Nsite; i++)	{
					for (int k=0; k<state.Nstate; k++)	{
						is >> state.SiteAASupport[i][k];
					}
				}
			}
		}

		if (state.mParam->NH)	{
			is >> state.NHNcat;
			for (int j=0; j<state.mParam->Nnode; j++)	{
				is >> state.NHcat[j];
			}
			for (int n=0; n<state.NHNcat; n++)	{
				is >> state.NHWeight[n];
				for (int k=0; k<state.mParam->Nstate; k++)	{
					is >>state.NHStatDistorter[n][k];
				}
				if (state.mParam->NHPrior)	{
					for (int k=0; k<state.ContNstate; k++)	{
						is >>state.logNHStatDistorter[n][k];
					}
				}
			}
			is >> state.NHpswitch;
			is >> state.NHStatAlpha;
			for (int k=0; k<state.ContNstate; k++)	{
				is >> state.NHStatCenter[k];
			}
			if (state.mParam->NHPrior)	{
				for (int k=0; k<state.mParam->ContNchar+1; k++)	{
					is >> state.NHContVar[k];
				}
				for (int k=0; k<state.ContNstate; k++)	{
					is >> state.NHVar[k];
				}
				if (state.mParam->NHPrior == 2)	{
					for (int k=0; k<state.ContNrr; k++)	{
						is >> state.NHCovIndex[k];
					}
				}
			}	
		}
		if (state.mParam->Qmode)	{
			is >> state.RRAlpha;
		}

		is >> state.GWf;

		if (state.mParam->PopEff)	{
			is >> state.PopAlpha;
			for (int j=0; j<state.NHNcat; j++)	{
				is >> state.PopEffSize[j];
			}
		}
		if (state.mParam->MutMode == 2)	{
			is >> state.Kappa;
			for (int i=0; i<Nnuc; i++)	{
				is >> state.RefACGT[i];
			}
		}
		if (state.mParam->MutMode == 3)	{
			for (int i=0; i<Nnuc; i++)	{
				is >> state.RefACGT[i];
			}
			for (int i=0; i<Nnuc * (Nnuc-1) / 2; i++)	{
				is >> state.RefRRACGT[i];
			}
		}
	}
	// CURRENT VERSION
	else if (state.Version >= 4)	{

		is >> state.Beta;

		// tree
		AmphiNode* p =state.tree;
		for (int i=0; i< state.mParam->Nnode; i++)	{
			int s;
			is >> s;
			p->up = (s == -1) ? 0 : &state.tree[s];
			is >> s;
			p->left = (s == -1) ? 0 : &state.tree[s];
			is >> s;
			p->right = (s== -1) ? 0 : &state.tree[s];
			is >> p->branchLength;
			p++;
		}
		int rootlabel;
		is >> rootlabel;
		state.root = & state.tree[rootlabel];

		for (int j=0; j<state.mParam->Nnode; j++)	{
			is >> state.BL[j];
			is >> state.RigidBL[j];
		}
		if (state.mParam->MBL)	{
			is >> state.MixBLAlpha;
			is >> state.NMixBL;
			for (int i=0; i<state.NMixBL; i++)	{
				for (int j=0; j<state.mParam->Nnode; j++)	{
					is >>state.MixBL[i][j];
				}
				state.MixBLSiteNumber[i] = 0;
			}
			for (int i=0; i<state.mParam->Nsite; i++)	{
				is >> state.MixBLMode[i];
			}
		}
		if (!state.mParam->SaveAll)	{
			is >> state.MeanLength;
			is >> state.gamma;
			is >> state.pconst;
			is >> state.alpha;
			is >> state.Nmode;
			is >> state.xi0;
			is >> state.invprob;
			is >> state.Sigma;
			is >> state.Theta;
			is >> state.Mu;
			is >> state.Scale;
			is >> state.mLogPrior;
			is >> state.mLogSampling;
		}
		else	{
			for (int i=0; i<state.mParam->Ngene; i++)	{
				for (int j=0; j<state.mParam->Nnode; j++)	{
					is >> state.GeneBL[i][j];
					is >> state.GeneRigidBL[i][j];
				}
			}

			is >> state.MeanLength;
			is >> state.VarLength;
			is >> state.LengthGamma;

			if (state.mParam->NormalApprox == No)	{

				// rate
				for (int i=0; i<state.mParam->Nsite; i++)	{
					is >> state.rate[i];
				}
				is >> state.gamma;
				is >> state.pconst;

				// gene
				for (int i=0; i<state.mParam->Ngene; i++)	{
					is >> state.GeneGamma[i];
					is >> state.GeneRate[i];
					for (int j=0; j<state.Nstate; j++)	{
						is >> state.GeneStationary[i][j];
					}
				}

				// ref
				for (int i=0; i<state.Nstate; i++)	{
					is >> state.RefStationary[i];
				}
				for (int i=0; i<state.Nstate*(state.Nstate-1)/2; i++)	{
					is >> state.RefRR[i];
				}

				// mode
				is >> state.alpha;
				is >> state.Nmode;
				for (int i=0; i<state.Nmode; i++)	{
					if (state.mParam->SaveStat)	{
						double total = 0;
						for (int j=0; j<state.Nstate; j++)	{
							is >> state.Stationary[i][j];
							total += state.Stationary[i][j];
						}
						for (int j=0; j<state.Nstate; j++)	{
							state.Stationary[i][j] /= total;
						}
					}
					if (state.mParam->Qmode)	{
						for (int j=0; j<state.Nstate*(state.Nstate-1)/2; j++)	{
							is >> state.RR[i][j];
						}
					}			
					is >> state.ModeWeight[i];
				}
				if (! state.mParam->SaveStat)	{
					is >> state.mStationaryEntropy;
				}

				for (int i=0; i<state.Nstate*(state.Nstate-1)/2; i++)	{
					is >> state.ModeRR[i];
				}

				for (int i=0; i<state.Nstate; i++)	{
					is >> state.ModeStatCenter[i];
				}
				is >> state.ModeStatAlpha;

				for (int i=0; i<state.mParam->Nsite; i++)	{
					is >> state.Mode[i];
				}
				// ratemode
				is >> state.NRateMode;
				for (int i=0; i<state.NRateMode; i++)	{
					is >> state.ModeRate[i];
					is >> state.RateModeWeight[i];
				}
				for (int i=0; i<state.mParam->Nsite; i++)	{
					is >> state.RateMode[i];
				}
				is >> state.RateAlpha;

				// log probs
				is >> state.mLogPrior;
				is >> state.mLogSampling;
				state.mLogPosterior = state.mLogPrior + state.Beta * state.mLogSampling;

				if ( state.mParam->SavePartialLogLikelihoods)	{
					double temp = 0;
					for (int i=0; i<state.mParam->Nsite; i++)	{
						is >> state.mSiteLogSampling[i];
						temp += state.mSiteLogSampling[i];
					}
					if (temp != state.mLogSampling)	{
						cerr << "warning : partial likelihoods non concordant with total logsampling advertised in this file\n";
					}
				}
			
				// hetero
				is >> state.xi0;
				is >> state.invprob;

			}
			else	{ // normal approx

				is >> state.mLogPrior;
				is >> state.mLogSampling;
				state.mLogPosterior = state.mLogPrior + state.Beta *  state.mLogSampling;

			}
			if (state.mParam->ActivateClock)	{
				is >> state.Sigma;
				is >> state.Theta;
				is >> state.Chi;
				is >> state.Chi2;
				/*	
				is >> state.MulSigma;
				is >> state.MulTheta;
				*/
				is >> state.MeanSigma;
				is >> state.MeanTheta;
				is >> state.VarSigma;
				is >> state.VarTheta;
				is >> state.Mu;
				is >> state.Scale;
				is >> state.VarScale;
				if (state.mParam->ContNchar)	{
					for (int k=0; k<state.mParam->ContNchar+1; k++)	{
						for (int j=0; j<state.mParam->Nnode; j++)	{
							is >> state.ContRho[k][j];
						}
						is >> state.ContRhoCenter[k];
					}
				}
				else	{
					for (int j=0; j<state.mParam->Nnode; j++)	{
						is >> state.Rho[j];
					}
				}
				for (int k=0; k<state.mParam->Ngene; k++)	{
					for (int j=0; j<state.mParam->Nnode; j++)	{
						is >> state.GeneRho[k][j];
					}
					is >> state.GeneSigma[k];
					is >> state.GeneTheta[k];
					is >> state.GeneMu[k];
				}
			}
		}

		if ((state.mParam->ZipGTR) && (state.mParam->ZipGTR != 4))	{
			is >> state.ZipAlpha;
			is >> state.ZipRate;
			is >> state.NZipMode;
			for (int i=0; i<state.mParam->Nsite; i++)	{
				is >> state.ZipMode[i];
			}
			for (int k=0; k<state.Nstate; k++)	{
				is >> state.ModeZipProb[k];
				is >> state.ZipA0[k];
				is >> state.ZipA1[k];
			}
			if (state.mParam->ZipGTR == 1)	{
				for (int i=0; i<state.Nmode; i++)	{
					for (int k=0; k<state.Nstate; k++)	{
						is >> state.ModeAASupport[i][k];
					}
				}
			}
			else if ((state.mParam->ZipGTR >= 2) && (state.mParam->ZipGTRDP))	{
				for (int i=0; i<state.NZipMode; i++)	{
					for (int k=0; k<state.Nstate; k++)	{
						is >> state.ModeAASupport[i][k];
					}
				}
			}
			else	{
				for (int i=0; i<state.mParam->Nsite; i++)	{
					for (int k=0; k<state.Nstate; k++)	{
						is >> state.SiteAASupport[i][k];
					}
				}
			}
		}

		if (state.mParam->NH)	{
			is >> state.NHNcat;
			for (int j=0; j<state.mParam->Nnode; j++)	{
				is >> state.NHcat[j];
			}
			for (int n=0; n<state.NHNcat; n++)	{
				is >> state.NHWeight[n];
				for (int k=0; k<state.mParam->Nstate; k++)	{
					is >>state.NHStatDistorter[n][k];
				}
				if (state.mParam->NHPrior)	{
					for (int k=0; k<state.ContNstate; k++)	{
						is >>state.logNHStatDistorter[n][k];
					}
				}
			}
			is >> state.NHpswitch;
			is >> state.NHStatAlpha;
			for (int k=0; k<state.ContNstate; k++)	{
				is >> state.NHStatCenter[k];
			}
			if (state.mParam->NHPrior)	{
				for (int k=0; k<state.mParam->ContNchar+1; k++)	{
					is >> state.NHContVar[k];
				}
				for (int k=0; k<state.ContNstate; k++)	{
					is >> state.NHVar[k];
				}
				if (state.mParam->NHPrior == 2)	{
					for (int k=0; k<state.ContNrr; k++)	{
						is >> state.NHCovIndex[k];
					}
				}
			}	
		}
		if (state.mParam->Qmode)	{
			is >> state.RRAlpha;
		}

		is >> state.GWf;

		if (state.mParam->PopEff)	{
			is >> state.PopAlpha;
			is >> state.PopBeta;
			for (int j=0; j<state.NHNcat; j++)	{
				is >> state.PopEffSize[j];
			}
		}
		if (state.mParam->MutMode == 2)	{
			is >> state.Kappa;
			for (int i=0; i<Nnuc; i++)	{
				is >> state.RefACGT[i];
			}
		}
		if (state.mParam->MutMode == 3)	{
			for (int i=0; i<Nnuc; i++)	{
				is >> state.RefACGT[i];
			}
			for (int i=0; i<Nnuc * (Nnuc-1) / 2; i++)	{
				is >> state.RefRRACGT[i];
			}
		}
	}
	else if ((state.Version == 3) && (state.SubVersion == 3))	{

		is >> state.Beta;

		// tree
		AmphiNode* p =state.tree;
		for (int i=0; i< state.mParam->Nnode; i++)	{
			int s;
			is >> s;
			p->up = (s == -1) ? 0 : &state.tree[s];
			is >> s;
			p->left = (s == -1) ? 0 : &state.tree[s];
			is >> s;
			p->right = (s== -1) ? 0 : &state.tree[s];
			is >> p->branchLength;
			p++;
		}
		int rootlabel;
		is >> rootlabel;
		state.root = & state.tree[rootlabel];

		for (int j=0; j<state.mParam->Nnode; j++)	{
			is >> state.BL[j];
			is >> state.RigidBL[j];
		}
		if (state.mParam->MBL)	{
			is >> state.MixBLAlpha;
			is >> state.NMixBL;
			for (int i=0; i<state.NMixBL; i++)	{
				for (int j=0; j<state.mParam->Nnode; j++)	{
					is >>state.MixBL[i][j];
				}
				state.MixBLSiteNumber[i] = 0;
			}
			for (int i=0; i<state.mParam->Nsite; i++)	{
				is >> state.MixBLMode[i];
			}
		}
		if (!state.mParam->SaveAll)	{
			is >> state.MeanLength;
			is >> state.gamma;
			is >> state.pconst;
			is >> state.alpha;
			is >> state.Nmode;
			is >> state.xi0;
			is >> state.invprob;
			is >> state.Sigma;
			is >> state.Theta;
			is >> state.Mu;
			is >> state.Scale;
			is >> state.mLogPrior;
			is >> state.mLogSampling;
		}
		else	{
			for (int i=0; i<state.mParam->Ngene; i++)	{
				for (int j=0; j<state.mParam->Nnode; j++)	{
					is >> state.GeneBL[i][j];
					is >> state.GeneRigidBL[i][j];
				}
			}

			is >> state.MeanLength;
			is >> state.VarLength;
			is >> state.LengthGamma;

			if (state.mParam->NormalApprox == No)	{

				// rate
				for (int i=0; i<state.mParam->Nsite; i++)	{
					is >> state.rate[i];
				}
				is >> state.gamma;
				is >> state.pconst;

				// gene
				for (int i=0; i<state.mParam->Ngene; i++)	{
					is >> state.GeneGamma[i];
					is >> state.GeneRate[i];
					for (int j=0; j<state.Nstate; j++)	{
						is >> state.GeneStationary[i][j];
					}
				}

				// ref
				for (int i=0; i<state.Nstate; i++)	{
					is >> state.RefStationary[i];
				}
				for (int i=0; i<state.Nstate*(state.Nstate-1)/2; i++)	{
					is >> state.RefRR[i];
				}

				// mode
				is >> state.alpha;
				is >> state.Nmode;
				for (int i=0; i<state.Nmode; i++)	{
					if (state.mParam->SaveStat)	{
						double total = 0;
						for (int j=0; j<state.Nstate; j++)	{
							is >> state.Stationary[i][j];
							total += state.Stationary[i][j];
						}
						for (int j=0; j<state.Nstate; j++)	{
							state.Stationary[i][j] /= total;
						}
					}
					if (state.mParam->Qmode)	{
						for (int j=0; j<state.Nstate*(state.Nstate-1)/2; j++)	{
							is >> state.RR[i][j];
						}
					}			
					is >> state.ModeWeight[i];
				}
				if (! state.mParam->SaveStat)	{
					is >> state.mStationaryEntropy;
				}

				for (int i=0; i<state.Nstate*(state.Nstate-1)/2; i++)	{
					is >> state.ModeRR[i];
				}

				for (int i=0; i<state.Nstate; i++)	{
					is >> state.ModeStatCenter[i];
				}
				is >> state.ModeStatAlpha;

				for (int i=0; i<state.mParam->Nsite; i++)	{
					is >> state.Mode[i];
				}
				// ratemode
				is >> state.NRateMode;
				for (int i=0; i<state.NRateMode; i++)	{
					is >> state.ModeRate[i];
					is >> state.RateModeWeight[i];
				}
				for (int i=0; i<state.mParam->Nsite; i++)	{
					is >> state.RateMode[i];
				}
				is >> state.RateAlpha;

				// log probs
				is >> state.mLogPrior;
				is >> state.mLogSampling;
				state.mLogPosterior = state.mLogPrior + state.Beta * state.mLogSampling;

				if ( state.mParam->SavePartialLogLikelihoods)	{
					double temp = 0;
					for (int i=0; i<state.mParam->Nsite; i++)	{
						is >> state.mSiteLogSampling[i];
						temp += state.mSiteLogSampling[i];
					}
					if (temp != state.mLogSampling)	{
						cerr << "warning : partial likelihoods non concordant with total logsampling advertised in this file\n";
					}
				}
			
				// hetero
				is >> state.xi0;
				is >> state.invprob;

			}
			else	{ // normal approx

				is >> state.mLogPrior;
				is >> state.mLogSampling;
				state.mLogPosterior = state.mLogPrior + state.Beta *  state.mLogSampling;

			}
			if (state.mParam->ActivateClock)	{
				is >> state.Sigma;
				is >> state.Theta;
				is >> state.Chi;
				is >> state.Chi2;
				/*	
				is >> state.MulSigma;
				is >> state.MulTheta;
				*/
				is >> state.MeanSigma;
				is >> state.MeanTheta;
				is >> state.VarSigma;
				is >> state.VarTheta;
				is >> state.Mu;
				is >> state.Scale;
				is >> state.VarScale;
				if (state.mParam->ContNchar)	{
					for (int k=0; k<state.mParam->ContNchar+1; k++)	{
						for (int j=0; j<state.mParam->Nnode; j++)	{
							is >> state.ContRho[k][j];
						}
						is >> state.ContRhoCenter[k];
					}
				}
				else	{
					for (int j=0; j<state.mParam->Nnode; j++)	{
						is >> state.Rho[j];
					}
				}
				for (int k=0; k<state.mParam->Ngene; k++)	{
					for (int j=0; j<state.mParam->Nnode; j++)	{
						is >> state.GeneRho[k][j];
					}
					is >> state.GeneSigma[k];
					is >> state.GeneTheta[k];
					is >> state.GeneMu[k];
				}
			}
		}

		if ((state.mParam->ZipGTR) && (state.mParam->ZipGTR != 4))	{
			is >> state.ZipAlpha;
			is >> state.ZipRate;
			is >> state.NZipMode;
			for (int i=0; i<state.mParam->Nsite; i++)	{
				is >> state.ZipMode[i];
			}
			for (int k=0; k<state.Nstate; k++)	{
				is >> state.ModeZipProb[k];
				is >> state.ZipA0[k];
				is >> state.ZipA1[k];
			}
			if (state.mParam->ZipGTR == 1)	{
				for (int i=0; i<state.Nmode; i++)	{
					for (int k=0; k<state.Nstate; k++)	{
						is >> state.ModeAASupport[i][k];
					}
				}
			}
			else if ((state.mParam->ZipGTR >= 2) && (state.mParam->ZipGTRDP))	{
				for (int i=0; i<state.NZipMode; i++)	{
					for (int k=0; k<state.Nstate; k++)	{
						is >> state.ModeAASupport[i][k];
					}
				}
			}
			else	{
				for (int i=0; i<state.mParam->Nsite; i++)	{
					for (int k=0; k<state.Nstate; k++)	{
						is >> state.SiteAASupport[i][k];
					}
				}
			}
		}

		if (state.mParam->NH)	{
			is >> state.NHNcat;
			for (int j=0; j<state.mParam->Nnode; j++)	{
				is >> state.NHcat[j];
			}
			for (int n=0; n<state.NHNcat; n++)	{
				is >> state.NHWeight[n];
				for (int k=0; k<state.mParam->Nstate; k++)	{
					is >>state.NHStatDistorter[n][k];
				}
				if (state.mParam->NHPrior)	{
					for (int k=0; k<state.ContNstate; k++)	{
						is >>state.logNHStatDistorter[n][k];
					}
				}
			}
			is >> state.NHpswitch;
			is >> state.NHStatAlpha;
			for (int k=0; k<state.ContNstate; k++)	{
				is >> state.NHStatCenter[k];
			}
			if (state.mParam->NHPrior)	{
				for (int k=0; k<state.mParam->ContNchar+1; k++)	{
					is >> state.NHContVar[k];
				}
				for (int k=0; k<state.ContNstate; k++)	{
					is >> state.NHVar[k];
				}
				if (state.mParam->NHPrior == 2)	{
					for (int k=0; k<state.ContNrr; k++)	{
						is >> state.NHCovIndex[k];
					}
				}
			}	
		}
		if (state.mParam->Qmode)	{
			is >> state.RRAlpha;
		}

		is >> state.GWf;

		if (state.mParam->PopEff)	{
			is >> state.PopAlpha;
			is >> state.PopBeta;
			for (int j=0; j<state.NHNcat; j++)	{
				is >> state.PopEffSize[j];
			}
		}
		if (state.mParam->MutMode == 2)	{
			is >> state.Kappa;
			for (int i=0; i<Nnuc; i++)	{
				is >> state.RefACGT[i];
			}
		}
		if (state.mParam->MutMode == 3)	{
			for (int i=0; i<Nnuc; i++)	{
				is >> state.RefACGT[i];
			}
			for (int i=0; i<Nnuc * (Nnuc-1) / 2; i++)	{
				is >> state.RefRRACGT[i];
			}
		}
	}
	else	{

		is >> state.Beta;

		// tree
		AmphiNode* p =state.tree;
		for (int i=0; i< state.mParam->Nnode; i++)	{
			int s;
			is >> s;
			p->up = (s == -1) ? 0 : &state.tree[s];
			is >> s;
			p->left = (s == -1) ? 0 : &state.tree[s];
			is >> s;
			p->right = (s== -1) ? 0 : &state.tree[s];
			is >> p->branchLength;
			p++;
		}
		int rootlabel;
		is >> rootlabel;
		state.root = & state.tree[rootlabel];

		for (int j=0; j<state.mParam->Nnode; j++)	{
			is >> state.BL[j];
			is >> state.RigidBL[j];
		}
		if (state.mParam->MBL)	{
			is >> state.MixBLAlpha;
			is >> state.NMixBL;
			for (int i=0; i<state.NMixBL; i++)	{
				for (int j=0; j<state.mParam->Nnode; j++)	{
					is >>state.MixBL[i][j];
				}
				state.MixBLSiteNumber[i] = 0;
			}
			for (int i=0; i<state.mParam->Nsite; i++)	{
				is >> state.MixBLMode[i];
			}
		}
		if (!state.mParam->SaveAll)	{
			is >> state.MeanLength;
			is >> state.gamma;
			is >> state.pconst;
			is >> state.alpha;
			is >> state.Nmode;
			is >> state.xi0;
			is >> state.invprob;
			is >> state.Sigma;
			is >> state.Theta;
			is >> state.Mu;
			is >> state.Scale;
			is >> state.mLogPrior;
			is >> state.mLogSampling;
		}
		else	{
			for (int i=0; i<state.mParam->Ngene; i++)	{
				for (int j=0; j<state.mParam->Nnode; j++)	{
					is >> state.GeneBL[i][j];
					is >> state.GeneRigidBL[i][j];
				}
			}

			is >> state.MeanLength;
			is >> state.VarLength;
			is >> state.LengthGamma;

			if (state.mParam->NormalApprox == No)	{

				// rate
				for (int i=0; i<state.mParam->Nsite; i++)	{
					is >> state.rate[i];
				}
				is >> state.gamma;
				is >> state.pconst;

				// gene
				for (int i=0; i<state.mParam->Ngene; i++)	{
					is >> state.GeneGamma[i];
					is >> state.GeneRate[i];
					for (int j=0; j<state.Nstate; j++)	{
						is >> state.GeneStationary[i][j];
					}
				}

				// ref
				for (int i=0; i<state.Nstate; i++)	{
					is >> state.RefStationary[i];
				}
				for (int i=0; i<state.Nstate*(state.Nstate-1)/2; i++)	{
					is >> state.RefRR[i];
				}

				// mode
				is >> state.alpha;
				is >> state.Nmode;
				for (int i=0; i<state.Nmode; i++)	{
					if (state.mParam->SaveStat)	{
						double total = 0;
						for (int j=0; j<state.Nstate; j++)	{
							is >> state.Stationary[i][j];
							total += state.Stationary[i][j];
						}
						for (int j=0; j<state.Nstate; j++)	{
							state.Stationary[i][j] /= total;
						}
					}
					if (state.mParam->Qmode)	{
						for (int j=0; j<state.Nstate*(state.Nstate-1)/2; j++)	{
							is >> state.RR[i][j];
						}
					}			
					is >> state.ModeWeight[i];
				}
				if (! state.mParam->SaveStat)	{
					is >> state.mStationaryEntropy;
				}

				for (int i=0; i<state.Nstate*(state.Nstate-1)/2; i++)	{
					is >> state.ModeRR[i];
				}

				for (int i=0; i<state.Nstate; i++)	{
					is >> state.ModeStatCenter[i];
				}
				is >> state.ModeStatAlpha;

				for (int i=0; i<state.mParam->Nsite; i++)	{
					is >> state.Mode[i];
				}
				// ratemode
				is >> state.NRateMode;
				for (int i=0; i<state.NRateMode; i++)	{
					is >> state.ModeRate[i];
					is >> state.RateModeWeight[i];
				}
				for (int i=0; i<state.mParam->Nsite; i++)	{
					is >> state.RateMode[i];
				}
				is >> state.RateAlpha;

				// log probs
				is >> state.mLogPrior;
				is >> state.mLogSampling;
				state.mLogPosterior = state.mLogPrior + state.Beta * state.mLogSampling;

				if ( state.mParam->SavePartialLogLikelihoods)	{
					double temp = 0;
					for (int i=0; i<state.mParam->Nsite; i++)	{
						is >> state.mSiteLogSampling[i];
						temp += state.mSiteLogSampling[i];
					}
					if (temp != state.mLogSampling)	{
						cerr << "warning : partial likelihoods non concordant with total logsampling advertised in this file\n";
					}
				}
			
				// hetero
				is >> state.xi0;
				is >> state.invprob;

			}
			else	{ // normal approx

				is >> state.mLogPrior;
				is >> state.mLogSampling;
				state.mLogPosterior = state.mLogPrior + state.Beta *  state.mLogSampling;

			}
			if (state.mParam->ActivateClock)	{
				is >> state.Sigma;
				is >> state.Theta;
				is >> state.Chi;
				is >> state.Chi2;
				/*
				is >> state.MulSigma;
				is >> state.MulTheta;
				*/
				is >> state.MeanSigma;
				is >> state.MeanTheta;
				is >> state.VarSigma;
				is >> state.VarTheta;
				is >> state.Mu;
				is >> state.Scale;
				is >> state.VarScale;
				if (state.mParam->ContNchar)	{
					for (int k=0; k<state.mParam->ContNchar+1; k++)	{
						for (int j=0; j<state.mParam->Nnode; j++)	{
							is >> state.ContRho[k][j];
						}
						is >> state.ContRhoCenter[k];
					}
				}
				else	{
					for (int j=0; j<state.mParam->Nnode; j++)	{
						is >> state.Rho[j];
					}
				}
				for (int k=0; k<state.mParam->Ngene; k++)	{
					for (int j=0; j<state.mParam->Nnode; j++)	{
						is >> state.GeneRho[k][j];
					}
					is >> state.GeneSigma[k];
					is >> state.GeneTheta[k];
					is >> state.GeneMu[k];
				}
			}
		}

		if ((state.mParam->ZipGTR) && (state.mParam->ZipGTR != 4))	{
			is >> state.ZipAlpha;
			is >> state.ZipRate;
			is >> state.NZipMode;
			for (int i=0; i<state.mParam->Nsite; i++)	{
				is >> state.ZipMode[i];
			}
			for (int k=0; k<state.Nstate; k++)	{
				is >> state.ModeZipProb[k];
				is >> state.ZipA0[k];
				is >> state.ZipA1[k];
			}
			if (state.mParam->ZipGTR == 1)	{
				for (int i=0; i<state.Nmode; i++)	{
					for (int k=0; k<state.Nstate; k++)	{
						is >> state.ModeAASupport[i][k];
					}
				}
			}
			else if ((state.mParam->ZipGTR >= 2) && (state.mParam->ZipGTRDP))	{
				for (int i=0; i<state.NZipMode; i++)	{
					for (int k=0; k<state.Nstate; k++)	{
						is >> state.ModeAASupport[i][k];
					}
				}
			}
			else	{
				for (int i=0; i<state.mParam->Nsite; i++)	{
					for (int k=0; k<state.Nstate; k++)	{
						is >> state.SiteAASupport[i][k];
					}
				}
			}
		}

		int cont;
		is >> cont;
		if (cont)	{
			if (state.mParam->NH)	{
				is >> state.NHNcat;
				for (int j=0; j<state.mParam->Nnode; j++)	{
					is >> state.NHcat[j];
				}
				for (int n=0; n<state.NHNcat; n++)	{
					is >> state.NHWeight[n];
					for (int k=0; k<state.mParam->Nstate; k++)	{
						is >>state.NHStatDistorter[n][k];
					}
					if (state.mParam->NHPrior)	{
						for (int k=0; k<state.ContNstate; k++)	{
							is >>state.logNHStatDistorter[n][k];
						}
					}
				}
				is >> state.NHpswitch;
				is >> state.NHStatAlpha;
				for (int k=0; k<state.ContNstate; k++)	{
					is >> state.NHStatCenter[k];
				}
				if (state.mParam->NHPrior)	{
					for (int k=0; k<state.mParam->ContNchar+1; k++)	{
						is >> state.NHContVar[k];
					}
					for (int k=0; k<state.ContNstate; k++)	{
						is >> state.NHVar[k];
					}
					if (state.mParam->NHPrior == 2)	{
						for (int k=0; k<state.ContNrr; k++)	{
							is >> state.NHCovIndex[k];
						}
					}
				}	
			}
			if (state.mParam->Qmode)	{
				is >> state.RRAlpha;
			}
			int cont;
			is >> cont;
			if (cont)	{
				is >> state.GWf;
				int cont;
				is >> cont;
				if (cont)	{
					cerr << "error in PhyloBayes::ifstream\n";
					exit(1);
				}
			}
		}
	}
	return is;
}


void PhyloBayes::WriteSubTotals(ofstream& os)	{
	for (int i=0; i<mParam->Nsite; i++)	{
		os << TotalSub[i] << '\t';
		for (int k=0; k<mParam->Nstate; k++)	{
			os << Nsub[i][k] << '\t';
		}
		os << '\n';
	}
}

void PhyloBayes::WriteRenormSubTotals(ofstream& os)	{
	for (int i=0; i<mParam->Nsite; i++)	{
		os << TotalSub[i] << '\t';
		for (int k=0; k<mParam->Nstate; k++)	{
			if (TotalSub[i])	{
				os << ((double) Nsub[i][k]) / TotalSub[i] << '\t';
			}
			else	{
				os << 0 << '\t';
			}
		}
		os << '\n';
	}
}

// ---------------------------------------------------------------------------
//		 InitFromStream
// ---------------------------------------------------------------------------

InitMode ReadMode(string s)	{

	if (s == "Uniform")	{
		return Uniform;
	}
	else if (s == "Random")	{
		return Random;
	}
	else if (s == "Empirical")	{
		return Empirical;
	}
	else if (s == "Fixed")	{
		return Fixed;
	}
	else if (s == "Poisson")	{
		return Poisson;
	}
	else if (s == "JTT")	{
		return JTT;
	}
	else if (s == "WAG")	{
		return WAG;
	}
	else if (s == "mtREV")	{
		return mtREV;
	}
	else if (s == "mtZOA")	{
		return mtZOA;
	}
	else if (s == "mtART")	{
		return mtART;
	}
	else if (s == "Fixed")	{
		return Fixed;
	}
	else return Error;
}

void PhyloBayes::InitFromStream(istream& is)	{

	try	{

		int through = 0;

		while( (! is.eof()) && (! through) )	{

			string key;
			is >> key;

			if (verbose)	{
				cerr << ":::" << key << '\n';
			}
			if (key == "Tree")	{
				ReadTree(is);
			}
			else if (key == "BranchLengths")	{
				ReadBranchLengths(is);
			}
			else if (key == "Rates")	{
				ReadRates(is);
			}
			else if (key == "Gamma")	{
				ReadGamma(is);
			}
			else if (key == "MeanLength")	{
				string s;
				is >> s;
				if (s == "Random")	{
					MeanLength= Random::sExpo();
					cerr << "MeanLength : " << MeanLength << '\n';
				}
				else	{
					MeanLength = Double(s);
					IsInitMeanLength = Yes;
				}
			}
			else if (key == "Pconst")	{
				ReadPconst(is);
			}
			else if (key == "Xi")	{
				ReadXi(is);
			}
			else if (key == "RefRR")	{
				ReadRefRR(is);
			}
			else if (key == "RefStationaries")	{
				ReadRefStationaries(is);
			}
			else if (key == "ModeRR")	{
				ReadModeRR(is);
			}
			else if (key == "ModeStationaries")	{
				ReadModeStationaries(is);
			}
			else if (key == "ModeStatCenter")	{
				ReadModeStatCenter(is);
			}
			else if (key == "ModeStatAlpha")	{
				ReadModeStatAlpha(is);
			}
			else if (key == "ModeAffiliations")	{
				ReadModeAffiliations(is);
			}
			else if (key == "Alpha")	{
				ReadAlpha(is);
			}
			else if (key == "RateAlpha")	{
				ReadRateAlpha(is);
			}
			else if(key == "Nmode")	{
				ReadNmode(is);
			}
			else if(key == "NRateMode")	{
				ReadNRateMode(is);
			}
			else if (key == "//")	{
				through = 1;
			}
			else	{
				cerr << "error : does not recognise : " << key << '\n';
				exit(1);

				throw(0);
			}
		}
	}

	catch(...)	{
		cerr << "cannot read phylobayes init state\n";
		cerr.flush();
		exit(1);
	}
	Initialise();
}



// ---------------------------------------------------------------------------
//		 SetBranchLengthFactor(double factor)
// ---------------------------------------------------------------------------

void PhyloBayes::SetBranchLengthFactor(double factor){

	for (int i=1; i<mParam->Nnode;i++)	{
		tree[i].branchLength *= factor;
	}
}


// ---------------------------------------------------------------------------
//		 ReadTree()
// ---------------------------------------------------------------------------

void PhyloBayes::ReadTree(istream& is)	{

	try	{

		if (SkipSeparators(is) == '(')	{	// read a tree in phylip format
			ReadPhylipFromStream(is);

			IsInitTree = Yes;
			/*
			cerr << "length : " << GetLength() << '\n';
			if (GetLength())	{
				IsInitBranchLength = Yes;
			}
			*/
		}
		else	{
			string temp;
			is >> temp;
			if (temp == "Random")	{
				InitTree = Random;
			}
			else	{
				throw;
			}
		}
	}
	catch(...)	{
		throw;
	}
}


// ---------------------------------------------------------------------------
//		 ReadBranchLengths()
// ---------------------------------------------------------------------------

void PhyloBayes::ReadBranchLengths(istream& is)	{

	try	{
		string temp;
		is >> temp;
		if (temp == "Random")	{
			InitBranchLength = Random;
		}
		else if (temp == "Uniform")	{
			InitBranchLength = Uniform;
			is >> InitInternalLength;
		}
		else	{
			throw;
		}
	}
	catch(...)	{
		cerr << "error in ReadBranchLengths\n";
		throw;
	}
}


// ---------------------------------------------------------------------------
//		 ReadRates()
// ---------------------------------------------------------------------------

void PhyloBayes::ReadRates(istream& is)	{

	try	{
		string temp;
		is >> temp;
		if (IsFloat(temp))	{
			rate[0] = Double(temp);
			for (int i=1; i<mParam->Nsite; i++)	{
				is >> rate[i];
			}
			IsInitRate = Yes;
		}
		else	{
			InitRate = ReadMode(temp);
		}
	}
	catch(...)	{
		throw;
	}
}


// ---------------------------------------------------------------------------
//		 ReadGamma()
// ---------------------------------------------------------------------------

void PhyloBayes::ReadGamma(istream& is)	{

	try	{

		string s;
		is >> s;

		if (s == "Random")	{
			gamma = mParam->GammaMax * Random::Uniform();
			IsInitGamma = Yes;
		}
		else if (IsFloat(s))	{
			gamma = Double(s);
			IsInitGamma = Yes;
		}
		else {
			InitGamma = ReadMode(s);
		}
	}
	catch(...)	{
		cerr << "error in ReadGamma()\n";
		throw;
	}
}

// ---------------------------------------------------------------------------
//		 ReadXi()
// ---------------------------------------------------------------------------

void PhyloBayes::ReadXi(istream& is)	{

	try	{

		string s;
		is >> s;

		if (IsFloat(s))	{
			xi0 = Double(s);
			IsInitXi = Yes;
		}
		else {
			InitXi = ReadMode(s);
		}
	}
	catch(...)	{
		cerr << "error in ReadXi()\n";
		throw;
	}
}


// ---------------------------------------------------------------------------
//		 ReadPconst()
// ---------------------------------------------------------------------------

void PhyloBayes::ReadPconst(istream& is)	{

	try	{

		string s;
		is >> s;

		if (IsFloat(s))	{
			pconst = Double(s);
			IsInitPconst = Yes;
		}
		else {
			InitPconst = ReadMode(s);
		}
	}
	catch(...)	{
		cerr << "error in ReadPconst()\n";
		throw;
	}
}


// ---------------------------------------------------------------------------
//		 ReadStationaries()
// ---------------------------------------------------------------------------

void PhyloBayes::ReadStationaries(double* stat, istream& is)	{

	try	{
		double total = 0;
		for (int i=0; i<Nstate; i++)	{
			is >> stat[i];
			total += stat[i];
		}
		for (int i=0; i<Nstate; i++)	{
			stat[i] /= total;
		}
	}
	catch(...)	{
		cerr << "error in ReadStationaries(double* istream&) \n";
		throw;
	}
}


// ---------------------------------------------------------------------------
//		 ReadRefStationaries()
// ---------------------------------------------------------------------------

void PhyloBayes::ReadRefStationaries(istream& is)	{

	try	{
		if ( IsDigit(SkipSeparators(is)) )	{
			ReadStationaries(RefStationary, is);
			IsInitRefStat = Yes;
		}
		else	{
			
			string temp;
			is >> temp;
			InitRefStat = ReadMode(temp);
		}
	}

	catch(...)	{
		cerr << "error in ReadRefStationaries\n";
		throw;
	}

}



// ---------------------------------------------------------------------------
//		 ReadRefRR()
// ---------------------------------------------------------------------------

void PhyloBayes::ReadRefRR(istream& is)	{

	try	{
		if ( IsDigit(SkipSeparators(is)) )	{
			cerr << "error : cannot read relative rates number by number\n";
			throw;
		}
		else	{

			string temp;
			is >> temp;
			InitRefRR = ReadMode(temp);
		}
	}

	catch(...)	{
		cerr << "error in ReadRefRR \n";
		throw;
	}

}


// ---------------------------------------------------------------------------
//		 ReadModeStationaries()
// ---------------------------------------------------------------------------

void PhyloBayes::ReadModeStationaries(istream& is)	{

	try	{
		if ( IsDigit(SkipSeparators(is)) )	{

			// should read them one by one
			for (int i=0; i<Nmode; i++)	{
				int temp;
				is >> temp;
				ReadStationaries(Stationary[i], is);
			}
			IsInitModeStat = Yes;
		}
		else	{

			string temp;
			is >> temp;
			InitModeStat = ReadMode(temp);
		}
	}

	catch(...)	{
		cerr << "error in ReadModeStationaries\n";
		throw;
	}

}


// ---------------------------------------------------------------------------
//		 ReadModeStatCenter()
// ---------------------------------------------------------------------------

void PhyloBayes::ReadModeStatCenter(istream& is)	{

	try	{
		if ( IsDigit(SkipSeparators(is)) )	{
			ReadStationaries(ModeStatCenter, is);
			IsInitModeStatCenter = Yes;
		}
		else	{

			string temp;
			is >> temp;
			InitModeStatCenter = ReadMode(temp);
		}
	}

	catch(...)	{
		cerr << "error in ReadModeStatCenter\n";
		throw;
	}

}

// ---------------------------------------------------------------------------
//		 ReadModeStatAlpha()
// ---------------------------------------------------------------------------

void PhyloBayes::ReadModeStatAlpha(istream& is)	{

	try	{
		string temp;
		is >> temp;
		ModeStatAlpha = Double(temp);
		IsInitModeStatAlpha = Yes;
	}

	catch(...)	{
		cerr << "error in ReadModeStatAlpha \n";
		throw;
	}

}


// ---------------------------------------------------------------------------
//		 ReadModeRR()
// ---------------------------------------------------------------------------

void PhyloBayes::ReadModeRR(istream& is)	{

	try	{
		if ( IsDigit(SkipSeparators(is)) )	{
			cerr << "error : cannot read relative rates number by number\n";
			throw;
		}
		else	{

			string temp;
			is >> temp;
			InitModeRR = ReadMode(temp);
		}
	}

	catch(...)	{
		cerr << "error in ReadModeRR \n";
		throw;
	}

}


// ---------------------------------------------------------------------------
//		 ReadModeAffiliations()
// ---------------------------------------------------------------------------


void PhyloBayes::ReadModeAffiliations(istream& is)	{

	try	{

		string temp;
		is >> temp;
		if (IsInt(temp))	{
			// then, read all the mode affiliations, one by one

			if ( (Nmode < 0) || (Nmode > mParam->Nsite) )	{
				cerr << "in reading mode affiliations : bad number of modes\n";
				throw;
			}
			Mode[0] = Int(temp);

			for (int i=1; i< mParam->Nsite; i++)	{
				is >> Mode[i];
			}
			IsInitModeAffiliation = Yes;
		}

		else	{
			// throw;
		}
	}

	catch(...)	{
		throw;
	}
}


// ---------------------------------------------------------------------------
//		 ReadAlpha()
// ---------------------------------------------------------------------------

void PhyloBayes::ReadAlpha(istream& is)	{

	try	{

		string s;
		is >> s;

		if (IsFloat(s))	{
			alpha = Double(s);
			IsInitAlpha = Yes;
		}
		else {
			InitAlpha = ReadMode(s);
		}
	}
	catch(...)	{
		cerr << "error in SetAlpha()\n";
		throw;
	}
}


// ---------------------------------------------------------------------------
//		 ReadRateAlpha()
// ---------------------------------------------------------------------------

void PhyloBayes::ReadRateAlpha(istream& is)	{

	try	{

		string s;
		is >> s;

		if (IsFloat(s))	{
			RateAlpha = Double(s);
			IsInitRateAlpha = Yes;
		}
		else {
			InitRateAlpha = ReadMode(s);
		}
	}
	catch(...)	{
		cerr << "error in SetRateAlpha()\n";
		throw;
	}
}


// ---------------------------------------------------------------------------
//		 ReadNmode()
// ---------------------------------------------------------------------------

void PhyloBayes::ReadNmode(istream& is)	{

	try	{

		string s;
		is >> s;

		if (IsInt(s))	{
			Nmode = Int(s);
			if ( (Nmode <=0 ) || (Nmode > mParam->Nsite) )	{
				cerr << "error in ReadNmode : bad number of modes " << Nmode << '\n';
				throw;
			}
			IsInitNmode = Yes;
		}
		else if ( s == "Nsite" )	{
			Nmode = mParam->Nsite;
			for (int i=0; i<mParam->Nsite; i++)	{
				Mode[i] = i;
			}
			IsInitNmode = Yes;
		}
		else {
			InitNmode = ReadMode(s);
		}
	}
	catch(...)	{
		cerr << "error in SetNmode()\n";
		throw;
	}
}


// ---------------------------------------------------------------------------
//		 ReadNRateMode()
// ---------------------------------------------------------------------------

void PhyloBayes::ReadNRateMode(istream& is)	{

	try	{

		string s;
		is >> s;

		if (IsInt(s))	{
			NRateMode = Int(s);
			if ( (NRateMode < 1) && (NRateMode > mParam->Nsite) && (NRateMode > mParam->NRateModeMax) )	{
				cerr << "error in ReadNmode : bad number of modes " << NRateMode << '\n';
				throw;
			}
			IsInitNRateMode = Yes;
		
			if (NRateMode == 1)	{	
				for (int i=0; i<mParam->Nsite; i++)	{
					RateMode[i] = 0;
				}
			}
			else	{
				for (int i=0; i<mParam->Nsite; i++)	{
					RateMode[i] = (int) (NRateMode * Random::Uniform());
				}
			}
		}
		else if ( s == "Nsite" )	{
			NRateMode = mParam->Nsite;
			for (int i=0; i<mParam->Nsite; i++)	{
				RateMode[i] = i;
			}
			IsInitNRateMode = Yes;
		}
		else {
			throw;
		}
	}
	catch(...)	{
		cerr << "error in SetNRateMode()\n";
		throw;
	}
}



// ---------------------------------------------------------------------------
//		 Initialise()
// ---------------------------------------------------------------------------

void PhyloBayes::Reinitialise()	{

	if (!mParam->FixTopo)	{
		IsInitTree = No;
	}
	IsInitBranchLength = No;
	if (! mParam->FixMeanLength)	{
		IsInitMeanLength = No;
	}

	if (!mParam->FixRate)	{
		IsInitRate = No;
	}
	if (!mParam->FixGamma)	{
		IsInitGamma = No;
	}
	IsInitRefRR = No;
	IsInitRefStat = No;

	if (!mParam->FixRR)	{
		IsInitModeRR = No;
	}
	if (!mParam->FixStat)	{
		IsInitModeStat = No;
	}
	if (!mParam->FixStatCenter)	{
		IsInitModeStatCenter = No;
		IsInitModeStatAlpha = No;
	}
	IsInitModeAffiliation = No;
	IsInitModeWeight = No;
	if (! mParam->FixNmode)	{
		IsInitNmode = No;
	}
	IsInitAlpha = No;
	IsInitNRateMode = No;
	IsInitRateAlpha = No;
	
	InitTree = Random;
	InitBranchLength = Random;
	InitRate = Random;
	InitGamma = Random;

	InitRefRR = LG;
	// InitRefRR = Uniform;
	InitRefStat = Random;

	InitModeRR = Random;
	InitModeStat = Random;
	InitModeStatCenter = Random;
	InitModeStatAlpha = Random;

	InitModeAffiliation = Random;
	InitNmode = Random;
	InitAlpha = Random;
	InitNRateMode = Random;
	InitRateAlpha = Random;

	Initialise();

}

void PhyloBayes::Initialise()	{

		if (! IsInitTree)	{
			SetTree(InitTree);
			IsInitTree = Yes;
		}

		if (! IsInitMeanLength)	{
			MeanLength = Random::sExpo();
			cerr << "MeanLength : " << MeanLength << '\n';
			IsInitMeanLength = Yes;
		}

		if (! IsInitBranchLength)	{
			SetBranchLengths(InitBranchLength);
			IsInitBranchLength = Yes;
		}

		if (! IsInitGamma)	{
			SetGamma(InitGamma);
			IsInitGamma = Yes;
		}

		if (! IsInitXi)	{
			SetXi(InitXi);
			IsInitXi = Yes;
		}

		if (! IsInitPconst)	{
			SetPconst(InitPconst);
			IsInitPconst = Yes;
		}

		if (! IsInitAlpha)	{
			SetAlpha(InitAlpha);
			IsInitAlpha = Yes;
		}

		if (! IsInitRateAlpha)	{
			SetRateAlpha(InitRateAlpha);
			IsInitRateAlpha = Yes;
		}

		if (! IsInitNH)	{
			SetNH();
			IsInitNH = Yes;
		}

		SetMode();
		SetRateMode();

		if (! IsInitRate)	{
			SetRates(InitRate);
			IsInitRate = Yes;
		}

		if (! IsInitRefRR)	{
			SetRelativeRates(RefRR, InitRefRR);
			IsInitRefRR = Yes;
		}
		
		if (! IsInitRefStat)	{
			SetStationaries(RefStationary, InitRefStat);
			IsInitRefStat = Yes;
		}
		
		if (! IsInitModeRR)	{
			SetRelativeRates(ModeRR, InitModeRR);
			if (mParam->Qmode)	{
				if (InitModeStat == CatFix)	{
					if (Nmode != mParam->Ncat)	{
						cerr << "error in init mode stat: number of modes should be " << mParam->Ncat << '\n';
						exit(1);
					}
					for (int mode=0; mode<mParam->Ncat; mode++)	{
						for (int k=0; k<Nrr; k++)	{
							RR[mode][k] = mParam->statrr[mode][k];
						}
					}
				}
				else	{
					for (int i=0; i<Nmode; i++)	{
						SetRelativeRates(RR[i], InitModeRR);
					}	
				}
			}	
			IsInitModeRR = Yes;
		}

		if (! IsInitModeStatCenter)	{
			SetStationaries(ModeStatCenter, InitModeStatCenter);
			IsInitModeStatCenter = Yes;
		}

		if (! IsInitModeStatAlpha)	{
			SetModeStatAlpha(InitModeStatAlpha);
			IsInitModeStatAlpha = Yes;
		}

		if (! IsInitModeStat)	{
			if (mParam->StatFlag)	{
				if (Nmode != mParam->Ncat)	{
					cerr << "error : Nmode and Ncat differ\n";
					exit(1);
				}
				for (int i=0; i<Nmode; i++)	{
					double total = 0;
					for (int j=0; j<Nstate; j++)	{
						Stationary[i][j] = mParam->StatFlag[i][j] + 0.01;
						total += Stationary[i][j];
					}
					for (int j=0; j<Nstate; j++)	{
						Stationary[i][j] /= total;
					}
				}
			}
			else if (InitModeStat == Random)	{
				for (int i=0; i<Nmode; i++)	{
					MakeRandomStationaries(Stationary[i], mParam->ModeStatPrior, ModeStatCenter, ModeStatAlpha);
				}
			}
			else if (InitModeStat == SiteEmpirical)	{

				for (int mode=0; mode<Nmode; mode++)	{
					for (int k=0; k<Nstate; k++)	{
						Stationary[mode][k] = 1;
					}
					for (int i=0; i<mParam->Nsite; i++)	{
						if (Mode[i] == mode)	{
							for (int k=0; k<Nstate; k++)	{
								Stationary[mode][k] += mParam->SiteEmpiricalCount[i][k];
							}
						}
					}
					double total = 0;
					for (int k=0; k<Nstate; k++)	{
						total += Stationary[mode][k];
					}
					for (int k=0; k<Nstate; k++)	{
						Stationary[mode][k] /= total;
					//	cerr << Stationary[mode][k] << '\t';
					}
					// cerr << '\n';
				}
			}
			else if (InitModeStat == CatFix)	{
				if (Nmode != mParam->Ncat)	{
					cerr << "error in init mode stat: number of modes should be " << mParam->Ncat << '\n';
					exit(1);
				}
				for (int mode=0; mode<mParam->Ncat; mode++)	{
					double total = 0;
					for (int k=0; k<Nstate; k++)	{
						if (mParam->statfix[mode][k] == -1)	{
							Stationary[mode][k] = mParam->EmpiricalFreq[k];
						}
						else	{	
							Stationary[mode][k] = mParam->statfix[mode][k];
						}
						total += Stationary[mode][k];
					}
					for (int k=0; k<Nstate; k++)	{
						Stationary[mode][k] /= total;
					}
					ModeWeight[mode] = mParam->weight[mode];
				}
			}
			else	{
				for (int i=0; i<Nmode; i++)	{
					SetStationaries(Stationary[i], InitModeStat);
				}
			}

			IsInitModeStat = Yes;
		}

		SetGene();

		if (verbose)	{
			cerr << "done\n";
			cerr.flush();
		}
	}


// ---------------------------------------------------------------------------
//		 SetStationaries()
// ---------------------------------------------------------------------------

void PhyloBayes::SetStationaries(double* stat, InitMode inMode)	{

	try	{

		// make a set of random stationaries
		if (inMode == Random)	{
			double total = 0;
			for (int i=0; i<Nstate; i++)	{
				stat[i] = Random::sGamma(1);
				total += stat[i];
			}
			for (int i=0; i<Nstate; i++)	{
				stat[i] /= total;
			}
		}

		// make a set of uniform stationaries
		else if (inMode == Uniform)	{
			for (int i=0; i<Nstate; i++)	{
				stat[i] = 1.0 / Nstate;
			}
		}

		else if (inMode == Empirical)	{
			for (int i=0; i<Nstate; i++)	{
				stat[i] = mParam->EmpiricalFreq[i];
			}
		}

		else if (inMode == JTT)	{
			for (int i=0; i<Nstate; i++)	{
				stat[i] = JTT_Stat[i];
			}
		}

		else if (inMode == WAG)	{
			for (int i=0; i<Nstate; i++)	{
				stat[i] = WAG_Stat[i];
			}
		}

		else if (inMode == mtREV)	{
			for (int i=0; i<Nstate; i++)	{
				stat[i] = mtREV_Stat[i];
			}
		}

		else if (inMode == mtZOA)	{
			for (int i=0; i<Nstate; i++)	{
				stat[i] = mtZOA_Stat[i];
			}
		}

		else if (inMode == mtART)	{
			for (int i=0; i<Nstate; i++)	{
				stat[i] = mtART_Stat[i];
			}
		}

		else if (inMode == LG)	{
			for (int i=0; i<Nstate; i++)	{
				stat[i] = LG_Stat[i];
			}
		}

		else	{
			cout << "could not find how to initialise\n";
			cout << inMode;
			throw;
		}
	}
	catch(...)	{
		throw;
	}
}


// ---------------------------------------------------------------------------
//		 SetGene()
// ---------------------------------------------------------------------------

void
PhyloBayes::SetGene()	{

	for (int i=0; i<mParam->Ngene; i++)	{
		GeneGamma[i] = 1;
		GeneRate[i] = 1;
		if (mParam->GeneBLMultiplier)	{
			for (int j=0; j<mParam->Nnode; j++)	{
				if (blfree(j))	{
					GeneBL[i][j] = 1;
				}
				else	{
					GeneBL[i][j] = 0;
				}
			}
		}
		else	{
			for (int j=0; j<mParam->Nnode; j++)	{
				GeneBL[i][j] = tree[j].branchLength;
			}
		}
		for (int j=0; j<Nstate; j++)	{
			GeneStationary[i][j] = RefStationary[j];
		}
	}
}

// ---------------------------------------------------------------------------
//		 SetNH()
// ---------------------------------------------------------------------------

void PhyloBayes::SetNH() {
	
	if (mParam->NH)	{
	if ((mParam->NH == 1) || (mParam->NH == 3))	{
		for (int j=0; j<mParam->Nnode; j++)	{
			NHcat[j] = (int) (NHNcat * Random::Uniform());
			if (j == root->label)	{
				NHcat[j] = 0;
			}
		}
	}
	else { // yang and robert's style
		for (int j=0; j<mParam->Nnode; j++)	{
			NHcat[j] = j;
			if (j == root->label)	{
				NHcat[j] = 0;
			}
		}
		NHcat[0] = root->label;
	}

	for (int n=0; n<NHNcat; n++)	{
		// MakeRandomStationaries(NHStatDistorter[n],MultiGamma,NHStatCenter,NHStatAlpha);
		for (int k=0; k<Nstate; k++)	{
			NHStatDistorter[n][k] = 1.0/Nstate;
		}
		if (mParam->NHPrior)	{
			for (int k=0; k<ContNstate; k++)	{
				logNHStatDistorter[n][k] = 0.1 * (Random::Uniform()-0.5);
			}
			if (mParam->NHClampFirst)	{
				logNHStatDistorter[n][0] = 0;
			}
			if (mParam->ContNchar)	{
				for (int j=0; j<mParam->Ntaxa; j++)	{
					for (int i=0; i<mParam->ContNchar; i++)	{
						if (! mParam->ContMissingData[j][i])	{
							int n = NHcat[j];
							logNHStatDistorter[n][Nstate + i] = mParam->ContData[j][i];
						}
					}
				}
			}
		}
		NHWeight[n] = 1.0;
	}
	for (int k=0; k<Nstate; k++)	{
		NHStatDistorter[0][k] = 1.0 / Nstate;
	}
	if (mParam->NHPrior)	{
		for (int k=0; k<ContNstate; k++)	{
			logNHStatDistorter[0][k] = 0;
			NHVar[k] = 1;
			NHStatCenter[k] = 0;
		}
		if (mParam->NHPrior == 2)	{
			for (int k=0; k<ContNrr; k++)	{
				NHCovIndex[k] = 0;
			}
		}
	}
	NHpswitch = 1;
	IsInitNH = Yes;
	}
}


// ---------------------------------------------------------------------------
//		 SetTree()
// ---------------------------------------------------------------------------

void PhyloBayes::SetTree(InitMode inMode)	{

	switch (inMode)	{

		case Random :
			MakeRandomTree();
			break;

		default :
			cerr << "error in SetTree\n";
			throw;
	}
}


// ---------------------------------------------------------------------------
//		 MakeRandomTree()
// ---------------------------------------------------------------------------

void
	PhyloBayes::MakeRandomTree( )	{

		for (int i=0; i<mParam->Nnode;i++)	{
			tree[i].label = i;
			tree[i].up = 0;
			tree[i].left = 0;
			tree[i].right = 0;
			tree[i].branchLength = 0;
		}
		
		for (int i=0; i<mParam->Ntaxa; i++)	{
			tree[i].height = 1;
		}

		for (int i=mParam->Ntaxa; i<mParam->Nnode; i++)	{
			tree[i].left = &tree[i];
			tree[i].right = &tree[i];
			tree[i].height = Random::Uniform();
		}


		if (mParam->OutNtaxa)	{
			int* taxflag = new int[mParam->Ntaxa];
			for (int i=0; i<mParam->Ntaxa; i++)	{
				taxflag[i] = 0;
				for (int j=0; j<mParam->OutNtaxa; j++)	{
					if (mParam->SpeciesNames[i] == mParam->OutGroup[j])	{
						taxflag[i] = 1;
					}
				}
			}
			int k = mParam->Ntaxa;
			int j = 0;
			int i = 0;
			int firstout = 0;
			int lastout = 0;
			while (i < mParam->OutNtaxa)	{
				while ((j<mParam->Ntaxa) && (!taxflag[j]))	{
					j++;
				}
				if (j == mParam->Ntaxa)	{
					cerr << "error in make random tree: overflow\n";
					cerr << "outgroup\n";
					exit(1);
				}

				// cerr << mParam->SpeciesNames[j] << '\t' << k << '\n';
				if (i)	{
					tree[j].prev = &tree[k-1];
					tree[k-1].next = &tree[j];
				}
				else	{
					firstout = j;
				}

				i++;

				if (i<mParam->OutNtaxa)	{
					tree[j].next = &tree[k];
					tree[k].prev = &tree[j];
					k++;
				}
				else	{
					lastout = j;
				}
				j++;

			}

			tree[lastout].next = &tree[firstout];
			tree[firstout].prev = &tree[lastout];
			// cerr << "outtree\n";
			AmphiNode* outroot = tree[firstout].makeTree();
			// cerr << "ok\n";

			j = 0;
			i = 0;
			int firstin = 0;
			int lastin = 0;

			while (i < mParam->Ntaxa - mParam->OutNtaxa)	{
				while ((j<mParam->Ntaxa) && (taxflag[j]))	{
					j++;
				}
				if (j == mParam->Ntaxa)	{
					cerr << "error in make random tree: overflow\n";
					cerr << "ingroup\n";
					exit(1);
				}

				if (i)	{
					tree[j].prev = &tree[k-1];
					tree[k-1].next = &tree[j];
				}
				else	{
					firstin = j;
				}

				i++;
				if (i<mParam->Ntaxa-mParam->OutNtaxa)	{
					tree[j].next = &tree[k];
					tree[k].prev = &tree[j];
					k++;
				}
				else	{
					lastin = j;
				}
				j++;
			}

			tree[lastin].next = &tree[firstin];
			tree[firstin].prev = &tree[lastin];
			// cerr << "intree\n";
			AmphiNode* inroot = tree[firstin].makeTree();
			// cerr << "ok\n";

			if (k != mParam->Nnode-1)	{
				cerr << "error in make random tree\n";
				cerr << "with outgroup\n";
				cerr << k << '\t' << mParam->Nnode << '\n';
				exit(1);
			}
			/*
			tree[lastin].next = &tree[k];
			tree[k].prev = &tree[lastin];
		
			tree[k].next = &tree[firstout];
			tree[firstout].prev = &tree[k];
		
			tree[lastout].next = &tree[firstin];
			tree[firstin].prev = &tree[lastout];
			*/
			tree[k].up = 0;
			tree[k].left = outroot;
			outroot->up = &tree[k];
			tree[k].right = inroot;
			inroot->up = &tree[k];
			root = &tree[k];

			/*
			tree[lastout].next = &tree[k];
			tree[k].prev = &tree[lastout];
		
			tree[k].next = &tree[firstin];
			tree[firstin].prev = &tree[k];
		
			tree[lastin].next = &tree[firstout];
			tree[firstout].prev = &tree[lastin];
			*/

			// tree[k].height = 0;

			// cerr << "maketree\n";
			// root = tree[firstin].makeTree();
			// root = tree[firstout].makeTree();

			// cerr << "set bl\n";
			SetBL();
			// cerr << "output tree\n";
			Phylip(cerr,1,0,1,0);
			delete[] taxflag;


		}
		else	{
			for (int i=mParam->Ntaxa; i<mParam->Nnode; i++)	{
				tree[i].prev = &tree[i-mParam->Ntaxa];
				tree[i].next = &tree[i-mParam->Ntaxa+1];
				tree[i-mParam->Ntaxa].next = &tree[i];
				tree[i-mParam->Ntaxa+1].prev= &tree[i];
			}
			tree[0].prev = &tree[mParam->Ntaxa-1];
			tree[mParam->Ntaxa-1].next = &tree[0];

			// cerr << "maketree\n";
			root = tree[0].makeTree();
			SetBL();
			Phylip(cerr,1,0,1,0);
		}
}

/*
void
	PhyloBayes::MakeRandomTree( )	{

		if (mParam->Ntaxa < 4)	{
			cerr << "error in PhyloBayes::MakeRandomTree: Ntaxa : " << mParam->Ntaxa << '\n';
			exit(1);
		}
		
		for (int i=0; i<mParam->Nnode;i++)	{
			tree[i].label = i;
			tree[i].up = 0;
			tree[i].left = 0;
			tree[i].right = 0;
			tree[i].branchLength = 0;
		}
		
		root = &tree[mParam->Nnode-1];
		int current = mParam->Ntaxa;
		root->left = &tree[0];
		tree[0].up = root;
		root->right = &tree[current];
		tree[current].up = root;
		tree[current].left = &tree[1];
		tree[1].up = &tree[current];
		tree[current].right = &tree[2];
		tree[2].up = &tree[current];
		current++;
		
		for (int i=3; i<mParam->Ntaxa; i++)	{
			// choose a node	
			int choose = (int) (Random::Uniform() * (2*i-2));
			if (choose >= i)	{
				choose += mParam->Ntaxa - i;
			}
			if (choose >= current)	{
				cerr << "error in PhyloBayes::MakeRandomTree : overflow\n";
				exit(1);
			}
			tree[current].right = &tree[i];
			tree[i].up = &tree[current];
			tree[current].left = &tree[choose];
			tree[current].up = tree[choose].up;
			if (tree[choose].up->left == &tree[choose])	{
				tree[choose].up->left = &tree[current];
			}
			else	{
				tree[choose].up->right = &tree[current];
			}
			tree[choose].up = &tree[current];
			current++;
		}
		if (current != mParam->Nnode-1)	{
			cerr << "error in PhyloBayes::MakeRandomTree\n";
			exit(1);
		}
		
		if (! root->isRoot())	{
			cerr << " PhyloBayes::PhyloBayes()	: root is not root\n";
		}

		for (int i=0; i<mParam->Ntaxa; i++)	{
			if (! tree[i].isLeaf()) 	{
				cerr << i << " is not leaf\n";
			}
		}
}
*/

// ---------------------------------------------------------------------------
//		 SetBranchLengths()
// ---------------------------------------------------------------------------

void PhyloBayes::SetBranchLengths(InitMode inMode)	{

	switch (inMode)	{

		case Random :
			// cerr << "draw random lengths\n";
			DrawLengths();
			break;

		case Uniform :
			MakeUniformBranchLengths();
			break;

		default :
			cerr << "error in SetBranchLengths\n";
			throw;
	}
	SetBL();
	SetMBL();

}

// ---------------------------------------------------------------------------
//		 SetMBL()
// ---------------------------------------------------------------------------

void PhyloBayes::SetMBL()	{

	UpdateTotalLength();
	if (mParam->MBL)	{
		NMixBL = 1;
		MixBLAlpha = 1;
		for (int i=0; i<mParam->Nsite; i++)	{
			MixBLMode[i] = 0;
		}
		/*
		for (int j=0; j<mParam->Nnode; j++)	{
			if (blfree(j))	{
				MixBL[0][j] =  1.0 / (mParam->Nnode-2);
			}
			else	{
				MixBL[0][j] =  1;
			}
		}
		*/

		/*
		NMixBL = mParam->Nsite;
		MixBLAlpha = 1;
		for (int i=0; i<mParam->Nsite; i++)	{
			MixBLMode[i] = i;
		}
		*/

		/*
		for (int i=0; i<mParam->Nsite; i++)	{
			for (int j=0; j<mParam->Nnode; j++)	{
				if (blfree(j))	{
					MixBL[i][j] =  1;
				}
				else	{
					MixBL[0][j] =  1;
				}
			}
		}
		*/
		for (int i=0; i<NMixBL; i++)	{
			DrawMixBL(i);
		}
	}
}

// ---------------------------------------------------------------------------
//		 MakeUniformBranchLengths()
// ---------------------------------------------------------------------------

void
	PhyloBayes::MakeUniformBranchLengths()	{

		for (int i=0; i<mParam->Nnode;i++)	{
			if (blfree(i))	{
				tree[i].branchLength = InitInternalLength;
			}
			else	{
				tree[i].branchLength = 0;
			}
		}
		SetBL();
	}

// ---------------------------------------------------------------------------
//		 SetRates()
// ---------------------------------------------------------------------------


void
	PhyloBayes::SetRates(InitMode inMode)	{

		switch (inMode)	{

			case Uniform:

				for (int i=0; i<mParam->Nsite; i++)	{
					rate[i] = 1.0;
				}
				break;

			case Random:
				DrawRates();
				break;

			default:
				cerr << "error : set rates\n";
				throw;
		}

	}

// ---------------------------------------------------------------------------
//		 SetGamma()
// ---------------------------------------------------------------------------


void
	PhyloBayes::SetGamma(InitMode inMode)	{

		if (inMode == Fixed)	{
			gamma = DefaultInitGammaValue;
		}
		else if (inMode == Random)	{
			if (mParam->GammaPrior == Exponential)	{
				gamma = Random::sExpo();
			}
			else if (mParam->GammaPrior == PowerLaw)	{
				double x = log(mParam->GammaMin) + (log(mParam->GammaMax) - log(mParam->GammaMin)) * Random::Uniform();
				gamma = exp(x);
			}
			else if (mParam->GammaPrior == Cauchy)	{
				double x = Random::sExpo();
				double y = Random::sExpo();
				gamma = x / y;
			}
			else	{
				cerr << "error: prior on alpha parameter not recognized\n";
				exit(1);
			}
		}
		else	{
			cerr << "error : does not know how to initialise gamma\n";
			throw;
		}
	}

// ---------------------------------------------------------------------------
//		 SetXi()
// ---------------------------------------------------------------------------


void
	PhyloBayes::SetXi(InitMode inMode)	{

		if (inMode == Fixed)	{
			xi0 = DefaultInitXiValue;
		}
		else if (inMode == Random)	{
			xi0 = Random::sExpo();
		}
		else	{
			cerr << "error : does not know how to initialise xi\n";
			throw;
		}
	}


// ---------------------------------------------------------------------------
//		 SetPconst()
// ---------------------------------------------------------------------------


void
	PhyloBayes::SetPconst(InitMode inMode)	{
		if (inMode == Fixed)	{
			pconst = DefaultInitPconstValue;
		}
		else if (inMode == Random)	{
			pconst = Random::Uniform();
		}
		else	{
			cerr << "error : do not know how to initialise Pconst\n";
			throw;
		}
		cerr << "pconst : " << pconst << '\n';
	}


// ---------------------------------------------------------------------------
//		 SetRelativeRates()
// ---------------------------------------------------------------------------


void
	PhyloBayes::SetRelativeRates(double* rr, InitMode inMode)	{

		if (! rr)	{
			cerr << "error in PhyloBayes::SetRelativeRates : 0 pointer !\n";
			throw;
		}

		if (inMode == WAG)	{
			double total = 0;
			for (int i=0; i<Nrr; i++)	{
				rr[i]= WAG_RR[i];
				total += rr[i];
			}
			for (int i=0; i<Nrr; i++)	{
			 	rr[i] /= total / Nrr;
			}
		}

		else if (inMode == WAGHSSP)	{
			double total = 0;
			for (int i=0; i<Nrr; i++)	{
				rr[i]= WAGHSSP_RR[i];
				total += rr[i];
			}
			for (int i=0; i<Nrr; i++)	{
			 	rr[i] /= total / Nrr;
			}
		}

		else if (inMode == LG)	{
			double total = 0;
			for (int i=0; i<Nrr; i++)	{
				rr[i]= LG_RR[i];
				total += rr[i];
			}
			for (int i=0; i<Nrr; i++)	{
			 	rr[i] /= total / Nrr;
			}
		}

		else if (inMode == JTT)	{
			double total = 0;
			for (int i=0; i<Nrr; i++)	{
				rr[i]= JTT_RR[i];
				total += rr[i];
			}
			for (int i=0; i<Nrr; i++)	{
				rr[i] /= total / Nrr;
			}
		}
		else if (inMode == mtREV)	{
			double total = 0;
			for (int i=0; i<Nrr; i++)	{
				rr[i]= mtREV_RR[i];
				total += rr[i];
			}
			for (int i=0; i<Nrr; i++)	{
				rr[i] /= total / Nrr;
			}
		}
		else if (inMode == mtZOA)	{
			double total = 0;
			for (int i=0; i<Nrr; i++)	{
				rr[i]= mtZOA_RR[i];
				total += rr[i];
			}
			for (int i=0; i<Nrr; i++)	{
				rr[i] /= total / Nrr;
			}
		}
		else if (inMode == mtART)	{
			double total = 0;
			for (int i=0; i<Nrr; i++)	{
				rr[i]= mtART_RR[i];
				total += rr[i];
			}
			for (int i=0; i<Nrr; i++)	{
				rr[i] /= total / Nrr;
			}
		}
		else if (inMode == Poisson)	{
			for (int i=0; i<Nrr; i++)	{
				rr[i]= 1.0;
			}
		}
		else if (inMode == Random)	{
			for (int i=0; i<Nrr; i++)	{
				rr[i]= Random::sGamma(1);
			}
		}
		else if (inMode == Uniform)	{
			for (int i=0; i<Nrr; i++)	{
				rr[i]= 1;
			}
		}
		else if (inMode == Cstom)	{
			if (! mParam->CustomRR)	{
				cerr << "error : Custom Rel Rates not defined\n";
				exit(1);
			}	
			double total = 0;
			for (int i=0; i<Nrr; i++)	{
				rr[i]= mParam->CustomRR[i];
				total += rr[i];
			}
			for (int i=0; i<Nrr; i++)	{
				rr[i] /= total / Nrr;
			}
		}	
		else	{
			cerr << "error in PhyloBayes:SetRelativeRates : dont know how to initialise\n";
			throw;
		}
	}


// ---------------------------------------------------------------------------
//		 SetModeStatAlpha()
// ---------------------------------------------------------------------------


void
	PhyloBayes::SetModeStatAlpha(InitMode inMode)	{

		if (inMode == Fixed)	{
			ModeStatAlpha = DefaultInitModeStatAlphaValue;
		}
		else if (inMode == Random)	{
			ModeStatAlpha = 0;
			while (ModeStatAlpha < ModeStatAlphaMin)	{
				ModeStatAlpha = Nstate * Random::sExpo();
			}
		}
		else	{
			cerr << "error : does not know how to initialise mode stat alpha\n";
			throw;
		}
	}

// ---------------------------------------------------------------------------
//		 SetModeAffiliations()
// ---------------------------------------------------------------------------


void PhyloBayes::SetModeAffiliations()	{

	if (Nmode == 1)	{
		for (int i=0; i<mParam->Nsite; i++)	{
			Mode[i] = 0;
		}
	}
	else if (Nmode == mParam->Nsite)	{
		for (int i=0; i<mParam->Nsite; i++)	{
			Mode[i] = i;
		}
	}
	else	{
		for (int i=0; i<mParam->Nsite; i++)	{
			Mode[i] = (int) (Nmode * Random::Uniform());
		}
	}
}

// ---------------------------------------------------------------------------
//		 SetAlpha()
// ---------------------------------------------------------------------------

void PhyloBayes::SetAlpha(InitMode inMode)	{

	switch(inMode)	{

		case Fixed :
			alpha = DefaultInitAlphaValue;
			break;

		case Random :
			if (mParam->AlphaPrior == Exponential)	{
				// cerr << "random alpha\n";
				alpha = 10 * Random::sExpo();
			}
			else	{
				alpha = mParam->AlphaMax * Random::Uniform();
			}
			break;

		default :
			cerr << "error in SetAlpha\n";
			throw;
	}
}

// ---------------------------------------------------------------------------
//		 SetRateAlpha()
// ---------------------------------------------------------------------------

void PhyloBayes::SetRateAlpha(InitMode inMode)	{

	switch(inMode)	{

		case Fixed :
			RateAlpha = DefaultInitAlphaValue;
			break;

		case Random :
			RateAlpha = Random::sExpo();
			break;

		default :
			cerr << "error in SetRateAlpha\n";
			throw;
	}
}

double PhyloBayes::CheckStat()	{

	double max = 0;
	for (int i=0; i<Nmode; i++)	{
		double total = 0;
		for (int k=0; k<Nstate; k++)	{
			total += Stationary[i][k];
		}
		double tmp = fabs(total-1);
		if (max < tmp)	{
			max = tmp;
		}
	}
	return max;
}

// ---------------------------------------------------------------------------
//		 SetNmode()
// ---------------------------------------------------------------------------

void PhyloBayes::SetNmode(InitMode inMode)	{

	if (inMode == Fixed)	{
		Nmode = DefaultInitNmodeValue;
	}
	else	{
		cerr << "error in SetNmode\n";
		throw;
	}
}


void PhyloBayes::SetMode()	{

	if (mParam->StartDP)	{
		if (mParam->StartDP == 1)	{
			cerr << "nmode = 1\n";
			Nmode = 1;
			for (int i=0; i<mParam->Nsite; i++)	{
				Mode[i] = 0;
			}
		}
		else if (mParam->StartDP == 2)	{
			cerr << "nmode = nsite\n";
			Nmode = mParam->Nsite;
			for (int i=0; i<mParam->Nsite; i++)	{
				Mode[i] = i;
			}
		}
		mParam->EmpiricalDP = 0;
		mParam->IncrementalDP = 0;
		IsInitNmode = Yes;
		IsInitModeAffiliation = Yes;
	}

	if (mParam->EmpiricalDP == 1)	{
		cerr << "draw empirical mode\n";
		DrawEmpiricalDPMode();
		IsInitNmode = Yes;
		IsInitModeAffiliation = Yes;
		IsInitModeWeight = Yes;
	}
	if ((mParam->ModeFastCompute) && (! mParam->ModePoisson) && (mParam->ZipGTR != 4))	{
		Nmode = 1;
		IsInitNmode = Yes;
		for (int k=0; k<Nstate; k++)	{
			ModeZipProb[k] = 0.5;
			Stationary[0][k] = Random::Uniform();
		}
		ZipAlpha = 1;
		for (int k=0; k<Nstate; k++)	{
			ZipA0[k] = 1;
			ZipA1[k] = 1;
		}
		NZipMode = 1;
		for (int i=0; i<mParam->Nsite; i++)	{
			ZipMode[i] = 0;
		}
		/*
		NZipMode = mParam->Nsite;
		for (int i=0; i<mParam->Nsite; i++)	{
			ZipMode[i] = i;
		}
		*/
		if ((mParam->ZipGTR != 4) && (mParam->ZipGTR >= 2))	{
			ZipRate = 1;
			if (mParam->ZipGTRDP)	{
				NZipMode = mParam->Nsite;
				for (int i=0; i<mParam->Nsite; i++)	{
					ZipMode[i] = i;
				}
				/*
				for (int k=0; k<Nstate; k++)	{
					ModeAASupport[0][k] = 1;
				}
				*/
				
				for (int i=0; i<mParam->Nsite; i++)	{
					for (int k=0; k<Nstate; k++)	{
						if (mParam->Orbit[i][k])	{
							SiteAASupport[i][k] = 1;
							ModeAASupport[i][k] = 1;
						}
						else	{
							SiteAASupport[i][k] = 0;
							ModeAASupport[i][k] = 0;
						}
						// SiteAASupport[i][k] = 1;
					}
					/*
					for (int k=0; k<Nstate; k++)	{
						SiteAASupport[i][k] = 1;
					}
					*/
				}
			}
			else	{
				for (int i=0; i<mParam->Nsite; i++)	{
					for (int k=0; k<Nstate; k++)	{
						if (mParam->Orbit[i][k])	{
							SiteAASupport[i][k] = 1;
						}
						else	{
							SiteAASupport[i][k] = 0;
						}
						// SiteAASupport[i][k] = 1;
					}
				}
			}
		}
		else	{
			for (int k=0; k<Nstate; k++)	{
				ModeAASupport[0][k] = 1;
			}

			Nmode = mParam->Nsite;
			for (int i=0; i<mParam->Nsite; i++)	{
				Mode[i] = i;
				for (int k=0; k<Nstate; k++)	{
					if (mParam->Orbit[i][k])	{
						ModeAASupport[i][k] = 1;
						Stationary[i][k] = Random::Uniform();
					}
					else	{
						ModeAASupport[i][k] = 0;
						Stationary[i][k] = 0;
					}
					// ModeAASupport[i][k] = 1;
				}
			}
		}
	}
	if (IsInitNmode)	{
		if (Nmode == 1)	{
			for (int i=0; i<mParam->Nsite; i++)	{
				Mode[i] = 0;
			}
		}
		else if (Nmode == mParam->Nsite)	{
			for (int i=0; i<mParam->Nsite; i++)	{
				Mode[i] = i;
			}
		}
		else	{
			// SetModeAffiliations();
			if (!IsInitModeWeight)	{
				// fraw from random distribution
				double total = 0;
				for (int i=0; i<Nmode; i++)	{
					ModeWeight[i] = Random::sGamma(1);
					total += ModeWeight[i];
				}
				for (int i=0; i<Nmode; i++)	{
					ModeWeight[i] /= total;
				}
				IsInitModeWeight = Yes;
			}
			if (! IsInitModeAffiliation)	{
				// draw from multinomial, given the weights
				for (int i=0; i<mParam->Nsite; i++)	{
					double q = Random::Uniform();
					double total = ModeWeight[0];
					int k = 0;
					while ((k<Nmode) && (q>total))	{
						k++;
						if (k<Nmode)	{
							total += ModeWeight[k];
						}
					}
					if (k==Nmode)	{
						cerr << "error in mode affiliations\n";
						exit(1);
					}
					Mode[i] = k;
				}
				IsInitModeAffiliation = Yes;
			}
		}
	}
	else	{
	
		// assumes Dirichlet process
		cerr << "draw dirichlet process\n";
		cerr.flush();
		DrawDPMode();
		IsInitNmode = Yes;
		IsInitModeAffiliation = Yes;
		IsInitModeWeight = Yes;
	}
}
				
void PhyloBayes::SetRateMode()	{

	if (mParam->GammaNcat)	{
		NRateMode = mParam->GammaNcat;
		for (int k=0; k<NRateMode; k++)	{
			RateModeWeight[k] = 1.0/NRateMode;
		}
		for (int i=0; i<mParam->Nsite; i++)	{
			RateMode[i] = (int) (Random::Uniform() * NRateMode);
		}
		IsInitNRateMode = Yes;
		IsInitRateModeWeight = Yes;
		IsInitRateModeAffiliation = Yes;
	}
	else if (IsInitNRateMode)	{
		if (!IsInitRateModeWeight)	{
			// fraw from random distribution
		}
		if (! IsInitRateModeAffiliation)	{
			// draw from multinomial, given the weights
		}
	}
	else	{
		/*
		cerr << "draw dp rate \n";
		// assumes Dirichlet process
		DrawDPRateMode();
		for (int k=0; k<NRateMode; k++)	{
			ModeRate[k] = Random::Gamma(gamma, gamma);
		}
		IsInitNRateMode = Yes;
		IsInitRateModeAffiliation = Yes;
		IsInitRateModeWeight = Yes;
		*/

		NRateMode = 1;
		for (int i=0; i<mParam->Nsite; i++)	{
			RateMode[i] = 0;
		}
		/*
		NRateMode = mParam->Nsite;
		for (int i=0; i<mParam->Nsite; i++)	{
			RateMode[i] = i;
		}
		*/
	}
}
				
// ---------------------------------------------------------------------------
//		 SetNRateMode()
// ---------------------------------------------------------------------------

void PhyloBayes::SetNRateMode(InitMode inMode)	{

	if (mParam->GammaNcat)	{
		NRateMode = mParam->GammaNcat;
	}
	else	{
		NRateMode = 1;
	}

	if (NRateMode == 1)	{	
		for (int i=0; i<mParam->Nsite; i++)	{
			RateMode[i] = 0;
		}
	}
	else if (NRateMode == mParam->Nsite)	{
		for (int i=0; i<mParam->Nsite; i++)	{
			RateMode[i] = i;
		}
	}
	else	{
		for (int i=0; i<mParam->Nsite; i++)	{
			RateMode[i] = (int) (NRateMode * Random::Uniform());
		}
	}
	IsInitNRateMode = Yes;
}

// ---------------------------------------------------------------------------
//		 SimulateData()
// ---------------------------------------------------------------------------

void PhyloBayes::ResampleModeWeight()	{

	// UpdateSiteNumber();
	double total = 0;	
	for (int k=0; k<Nmode; k++)	{
		ModeWeight[k] = Random::sGamma(1 + SiteNumber[k]);
		// ModeWeight[k] = 1 + SiteNumber[k];
		total += ModeWeight[k];
	}
	for (int k=0; k<Nmode; k++)	{
		ModeWeight[k] /= total;
	}
}


void PhyloBayes::DrawSiteVariables(int priormode, int priorratemode)	{

	UpdateBaseFields();
	if (Nmode > 1)	{
		if (priormode)	{
			if (!mParam->SumOverModes)	{
				ResampleModeWeight();
			}
			double total = 0;
			for (int k=0; k<Nmode; k++)	{
				total += ModeWeight[k];
			}
			for (int i=0; i<mParam->Nsite; i++)	{
				double U = Random::Uniform() * total;
				double t = ModeWeight[0];
				int index = 0;
				while ((index < Nmode) && (U > t))	{
					index++;
					if (index == Nmode)	{
						cerr << "error in simu data: resampling mode affiliations\n";
						exit(1);
					}
					t+= ModeWeight[index];
				}
				SetMode(i,index);
			}
		}
		else	{
			if (mParam->SumOverModes)	{
				UpdateModeSiteLogSampling();
				SumModeSiteLogSampling();
			}
		}
	}
	UpdateBaseFields();
	if (priorratemode)	{
		DrawRates();
	}
	else	{
		if (mParam->SumOverRateModes)	{
			UpdateRateModeSiteLogSampling();
			SumRateModeSiteLogSampling();
		}
	}
	UpdateBaseFields();

}


void PhyloBayes::RestoreData()	{

	if (! BKData)	{
		cerr << "error in restore data : no backup\n";
		exit(1);
	}
	for (int j=0; j<mParam->Ntaxa; j++)	{
		for (int i=0; i<mParam->Nsite; i++)	{
			Data[j][i] = BKData[j][i];
			if (mParam->ModeFastCompute)	{
				ZipData[j][i] = BKZipData[j][i];
			}
		}
	}
	if ((mParam->ZipGTR) && (mParam->ZipGTR != 4))	{
		UpdateZipGTRData();
	}
}

void PhyloBayes::SimulateData(int priormode, int priorratemode, int priorrootstate, int redraw, ostream* os, int* sitemask)	{

	if (redraw)	{
		DrawSiteVariables(priormode, priorratemode);
	}

	if ((! priorrootstate) && redraw)	{
		ResampleSub();
	}
	if (! BKData)	{
		BKData = new int*[mParam->Ntaxa];
		for (int j=0; j<mParam->Ntaxa; j++)	{
			BKData[j] = new int[mParam->Nsite];
		}
		if (mParam->ModeFastCompute)	{
			BKZipData = new int*[mParam->Ntaxa];
			for (int j=0; j<mParam->Ntaxa; j++)	{
				BKZipData[j] = new int[mParam->Nsite];
			}
		}
	}

	for (int j=0; j<mParam->Ntaxa; j++)	{
		for (int i=0; i<mParam->Nsite; i++)	{
			BKData[j][i] = Data[j][i];
			Data[j][i] = unknown;
			if (mParam->ModeFastCompute)	{
				BKZipData[j][i] = ZipData[j][i];
				ZipData[j][i] = unknown;
				if ((mParam->ZipGTR) && (mParam->ZipGTR != 4))	{
					ZipGTRData[j][i] = unknown;
				}
			}
		}
	}

	/*
	if (mParam->ZipGTR)	{
		UpdateZipGTRData();
	}
	*/
	for (int i=0; i<mParam->Nsite; i++)	{
		ConstantStatus[i] = DrawConstantStatus(i);
	}

	if (! priorrootstate)	{
		mParam->ResampleStateAtRoot = No;
	}

	ResampleSub(os,sitemask);

	if (! priorrootstate)	{
		mParam->ResampleStateAtRoot = Yes;
	}

	for (int j=0; j<mParam->Ntaxa; j++)	{
		for (int i=0; i<mParam->Nsite; i++)	{
			if (BKData[j][i] == unknown)	{
				Data[j][i] = unknown;
			}
			else	{
				if (mParam->ModeFastCompute)	{
					Data[j][i] = TrueState[j][i];
				}
				else	{
					Data[j][i] = State[j][i];
				}
			}
		}
	}
}

