#include "phylo.h"

extern int SamCompare(int nchain, int burnin, int stop, string* ChainName, double& disc, double& overlap, double& effsize, string outname);

// make a new chain

Chain::Chain(MCParameters* param, string Name)	{

	mParam = param;
	ChainName = Name;
	mParam->Reset();
	MakeFiles();
}

Chain::Chain(string InitFile, string Name)	{

	ChainName = Name;
	mParam = new MCParameters();
	mParam->InitFromFile(InitFile);
	mParam->Reset();
	MakeFiles();

}

// read an already existing chain

Chain::Chain(MCParameters* param, string Name, int a)	{

	mParam = param;
	if (mParam->HowManySaved)	{
		mParam->IncrementalDP = 0;
		mParam->IncrementalRateDP = 0;
	}
	ChainName = Name;
	if (! ifstream((ChainName+".param").c_str()))	{
		cerr << "error : non existing chain\n";
		cerr << ChainName << '\n';
		exit(1);
	}
	ifstream Param_is((ChainName + ".param").c_str());
	Param_is >> *mParam;
	mParam->Update();
}

Chain::Chain(string Name, int a, string inPath)	{

	ChainName = Name;
	mParam = new MCParameters();
	mParam->SetPath(inPath);
	ifstream Param_is((ChainName+".param").c_str());
	if (! Param_is)	{
		cerr << "?? non existing chain : " << ChainName << '\n';
		exit(1);
	}
	Param_is >> *mParam;
	if (mParam->HowManySaved)	{
		mParam->IncrementalDP = 0;
		mParam->IncrementalRateDP = 0;
	}
	mParam->Update();
}

Chain::~Chain()	{
	// delete mParam;
}


void Chain::MakeFiles()	{

	int nchain = 1;
	if (mParam->SaveAllChains)	{
		nchain = mParam->Nchain;
	}

	if (mParam->MSMode != None)	{
		ofstream os((ChainName + ".thermo").c_str());
	}
	if (mParam->SaveAll)	{
		if ((mParam->MSMode == None) && (nchain > 1))	{
			for (int chain = 0; chain < nchain; chain++)	{
				ostringstream name;
				name << ChainName << '.' << chain;
				ofstream os((name.str() + ".chain").c_str());
			}
		}
		else	{
			ofstream os((ChainName + ".chain").c_str());
			// os << *mParam->GetCurrentState(chain);
		}	

		if (mParam->Nchain > 1)	 {
			ofstream SwapHisto_os((ChainName + ".swaphisto").c_str());
			SwapHisto_os << mParam->HowManySaved << '\t';
			for (int i=0; i<mParam->Nchain; i++)	{
				int j=0;
				while ((j<mParam->Nchain) && (mParam->SwapPermut[j] != i))	{
					j++;
				}
				if (j == mParam->Nchain)	{
					cerr << "error in swap monitoring\n";
					exit(1);
				}
				SwapHisto_os << j << '\t';
				// SwapHisto_os << mParam->SwapPermut[i] << '\t';
			}
			SwapHisto_os << '\n';
		}
	}

	for (int chain=0; chain<mParam->Nchain; chain++)	{
		ostringstream s;
		if (mParam->Nchain == 1)	{
			s << ChainName << ".trace";
		}
		else	{
			s << ChainName << "." << chain+1 << ".trace";
		}
		ofstream Trace_os(s.str().c_str());
		Trace_os << "#cycle\t#treegen\ttime\tloglik\tlength";
		if (mParam->ActivateClock)	{
			Trace_os << "\tsigma";
			if ((mParam->ClockModel == CIRrigid) || (mParam->ClockModel == CIRrigid2) || (mParam->FlexClockModel == CIRflex)) 	{
				Trace_os << "\ttheta";
			}
			if (mParam->ClockModel == AutoRegressive)	{
				Trace_os << "\ttheta";
			}
			Trace_os << "\tmu";
			Trace_os << "\tmeanrates";
			// Trace_os << "\tvarrates";
			if (mParam->TimePrior == Dirich)	{
				Trace_os << "\tdir_alpha";
			}
			if (mParam->TimePrior == BD)	{
				Trace_os << "\tp1\tp2";
			}
			if (mParam->NCalib)	{
				Trace_os << "\tscale";
			}
		}
		if (mParam->NormalApprox == No)	{
			Trace_os << "\talpha";
			if (mParam->MutMode)	{
				if (mParam->MutMode == 2)	{
					Trace_os << "\tkappa\tA\tC\tG\tT";
				}
				else if (mParam->MutMode == 3)	{
					Trace_os << "\trrent\tmeanrr\tA\tC\tG\tT";
				}
			}
			if (! mParam->FixPconst)	{
				Trace_os << "\tpinv";
			}
			if (mParam->NRateModeMax == mParam->Nsite + 5)	{
				Trace_os << "\tnratecat\tmeanrate";
			}
			if (mParam->GeneBLMultiplier)	{
				Trace_os << "\tgenegamma";
				// Trace_os << "\tgenemeanbl";
			}
			if (mParam->MBL)	{
				Trace_os << "\tmblncat\tmbleffncat";
			}
			if (mParam->HeteroMode == Covarion)	{
				Trace_os << "\txi\tinvprob";
			}
			Trace_os << "\tnmode\tstat";
			Trace_os << "\tstatalpha";
			if (!mParam->FixRR)	{
				Trace_os << "\trrent";
				Trace_os << "\tmeanrr";
				if (mParam->RRPrior == GammaDistributed)	{
					Trace_os << "\trralpha";
				}
				if (mParam->Qmode)	{
					Trace_os << "\trefrrent";
				}
				if (mParam->MutMode == 1)	{
					Trace_os << "\trefstat";
				}
			}
			Trace_os << "\tkappa\tallocent";
			if ((mParam->NH == 1) || (mParam->NH == 3))	{
				Trace_os << "\tnhn_occ";
			}
			if (mParam->NH == 1)	{
				Trace_os  << "\tpsw";
			}
			if (mParam->NH && (! mParam->PopEff || mParam->NHPop))	{
				Trace_os << "\tnhent\tnhconc";
			}
			if (mParam->PopEff)	{
				Trace_os << "\tmeanbeta\tstdev\talpha_param\tbeta_param";
			}
			if (mParam->ZipGTR == 4)	{
				Trace_os << "\torbsize";
			}
			else if (mParam->ZipGTR && (mParam->ZipGTR != 3))	{
				Trace_os << "\torbsize\tminorbsize";
			}
			if (mParam->PopEff)	{
				Trace_os << "\tmeanpopeff\tstdev";
			}
			// Trace_os << "\tlogprior";
		}
		Trace_os << '\n';
	}

	ofstream Run_os((ChainName + ".run").c_str());
	ofstream TreeList_os((ChainName + ".treelist").c_str());

	ofstream Param_os((ChainName + ".param").c_str());
	Param_os.precision(30);
	Param_os << *mParam;
}

// ---------------------------------------------------------------------------
//		 RunningStatus
// ---------------------------------------------------------------------------

int Chain::RunningStatus()	{
	ifstream ris((ChainName + ".run").c_str());
	int i;
	ris >> i;
	ris.close();
	return i;
}

// ---------------------------------------------------------------------------
//		 CheckDiff
// ---------------------------------------------------------------------------

int Chain::CheckDiff()	{

	int burninfactor = mParam->BurnInDone; // in fact, burnin factor
	if (mParam->Nchain == 1) return 1;
	if ((mParam->Nchain > 1) && (mParam->HowManySaved >= burninfactor*mParam->BurnIn) && !(mParam->HowManySaved % mParam->BurnIn))	{
		int bpdiff = 0;
		string* name = new string[mParam->Nchain];
		for (int chain = 0; chain<mParam->Nchain; chain++)	{
			ostringstream s;
			s << ChainName << "." << chain+1;
			name[chain] = s.str();
		}

		if (mParam->FixTopo)	{
			bpdiff = 0;
		}
		else	{
			double diff = BPCompare(name, mParam->Nchain, mParam->HowManySaved / burninfactor, 10, -1, 0, 0, 1, ChainName, mParam->FinalBeta, 0.5);
			ofstream os((ChainName + ".diff").c_str());
			os << "max diff :  " << diff << "\n";
			bpdiff = (diff > mParam->FinalBeta);
		}

		int contdiff = 0;
		double disc = 0;
		double overlap = 0;
		double effsize = 0;
		for (int chain = 0; chain<mParam->Nchain; chain++)	{
			name[chain] = name[chain] + ".trace";
		}
		// cerr << burninfactor << '\n';
		SamCompare(mParam->Nchain,mParam->HowManySaved/burninfactor,mParam->HowManySaved-1,name,disc,overlap,effsize,ChainName + ".contdiff");
	
		// double shift = 1 - overlap / 100;
		contdiff = (effsize < mParam->BetaStep) || (disc  > mParam->FinalBeta);
		// contdiff = (effsize < mParam->HowManySaved/burninfactor) || (shift  > mParam->FinalBeta);
		delete[] name;

		return contdiff || bpdiff;
	}



	/*
	if ((mParam->Nchain > 1) && (mParam->HowManySaved >= 2*mParam->BurnIn) && !(mParam->HowManySaved % mParam->BurnIn))	{
		string* name = new string[mParam->Nchain];
		for (int chain = 0; chain<mParam->Nchain; chain++)	{
			ostringstream s;
			s << ChainName << "." << chain+1;
			name[chain] = s.str();
		}
		double diff = BPCompare(name, mParam->Nchain, mParam->BurnIn, 10, -1, 0, 0, 1, ChainName, mParam->FinalBeta, 0.5);
		ofstream os((ChainName + ".diff").c_str());
		os << "max diff :  " << diff << "\n";
		return (diff > mParam->FinalBeta);
	}
	*/
	return 1;
}


// ---------------------------------------------------------------------------
//		 SavePoint
// ---------------------------------------------------------------------------

void Chain::SavePoint(string ext)	{

	if (ext == "")	{
		ext = ".chain";
	}

	if (mParam->SaveAll)	{
		ofstream Chain_os((ChainName + ext).c_str(), APPEND);
		Chain_os << *mParam->GetCurrentState();
		Chain_os.close();
	}

	if (mParam->Nchain == 1)	{
		ofstream TreeList_os((ChainName + ".treelist").c_str(), APPEND);
		if (! mParam->NormalApprox)	{
			mParam->GetCurrentState()->NormaliseLengths();
		}
		mParam->GetCurrentState()->Phylip(TreeList_os,1,0,1,0);
		if (! mParam->NormalApprox)	{
			mParam->GetCurrentState()->DenormaliseLengths();
		}
		TreeList_os.close();
	}
	else	{
		for (int chain=0; chain<mParam->Nchain; chain++)	{
			ostringstream s;
			s << ChainName << "." << chain+1 << ".treelist";
			ofstream TreeList_os(s.str().c_str(), APPEND);
			if (! mParam->NormalApprox)	{
				mParam->GetCurrentState(chain)->NormaliseLengths();
			}
			mParam->GetCurrentState(chain)->Phylip(TreeList_os,1,0,1,0);
			if (! mParam->NormalApprox)	{
				mParam->GetCurrentState(chain)->DenormaliseLengths();
			}
			TreeList_os.close();
		}
	}

	mParam->HowManySaved++;
	ofstream Param_os((ChainName + ".param").c_str());
	ostringstream dummy;
	dummy << 0;

	Param_os.precision(30);
	Param_os << *mParam;
	Param_os.close();
}

// ---------------------------------------------------------------------------
//		 Start
// ---------------------------------------------------------------------------


int Chain::Start()	{
	ifstream ris((ChainName + ".run").c_str());
	int check = 0;
	ris >> check;
	if (check != -1)	{
		ofstream ros((ChainName + ".run").c_str());
		ros << 1 << '\n';
		ros.close();
		return Run();
	}
	return 0;
}


// ---------------------------------------------------------------------------
//		 ThermoIntegrate
// ---------------------------------------------------------------------------


double Chain::Annealing(int burnin, int every, int size)	{

	double returnValue = 0;
	double t01 = ThermoIntegrate(burnin, every, size, 1, Thermo, 1);
	double t10 = ThermoIntegrate(burnin, every, size, -1, Thermo, 1);
	ofstream thermo_os((ChainName + ".thermo").c_str());
	thermo_os << t01 << '\n';
	thermo_os << t10 << '\n';
	returnValue = 0.5 * (t01 + t10);
	return returnValue;
}



// ---------------------------------------------------------------------------
//		 HomoProb
// ---------------------------------------------------------------------------


void Chain::HomoProb(int nsub)	{

	double homoprob[nsub];
	for (int k=0; k<nsub; k++)	{
		homoprob[k] = 0;
	}

	PhyloBayes* pb = mParam->GetCurrentState();
	int Nstate = mParam->Nstate;
	int Nsite = mParam->Nsite;

	double probtrans[Nstate][Nstate];

	for (int site=0; site<Nsite; site++)	{

		// first make the transition probability matrix
		for (int i=0; i<Nstate; i++)	{
			double total = 0;
			for (int j=0; j<Nstate; j++)	{
				if (j==i)	{
					probtrans[i][j] = 0;
				}
				else 	{
					if (mParam->ModeFastCompute)	{ // Poisson	
						probtrans[i][j] = pb->Stationary[pb->Mode[site]][j];
					}
					else	{
						probtrans[i][j] = pb->mMatrixArray[pb->Mode[site]]->operator()(i,j);
					}
					total += probtrans[i][j];
				}
			}
			for (int j=0; j<Nstate; j++)	{
				probtrans[i][j] /= total;
			}
		}

		// iterate over one, two ... substitutions
		double p[Nstate];
		double q[Nstate];
		double* stat = pb->Stationary[pb->Mode[site]];

		for (int i=0; i<Nstate; i++)	{
			for (int j=0; j<Nstate; j++)	{
				p[j] = 0;
			}
			p[i] = 1;
			
			for (int k=0; k<nsub; k++)	{
				for (int l=0; l<Nstate; l++)	{
					double total = 0;
					for (int m=0; m<Nstate; m++)	{
						total += p[m] * probtrans[m][l];
					}
					q[l] = total;
				}
				for (int l=0; l<Nstate; l++)	{
					p[l] = q[l];
				}
				homoprob[k] += stat[i] * p[i];		
			}
		}	
	}
	for (int k=0; k<nsub; k++)	{
		homoprob[k] /= Nsite;
		cerr << k+1 << '\t' << homoprob[k] << '\n';
	}
}

// ---------------------------------------------------------------------------
//		 MeanSiteStat
// ---------------------------------------------------------------------------


/*
void Chain::MeanSiteStat(string ext, string outext)	{

	int style = 1;
	int Nstate = mParam->Nstate;
	int Nsite = mParam->Nsite;

	double stat[Nsite][Nstate];
	double rate[Nsite];
	for (int i=0; i<Nsite; i++)	{
		rate[i] = 0;
		for (int k=0; k<Nstate; k++)	{
			stat[i][k] = 0;
		}
	}
	
	ifstream is((ChainName + ext).c_str());
	int size;
	is >> size;
	cerr << "size : " << size << '\n';
	cerr.flush();

	PhyloBayes* pb = mParam->GetCurrentState();

	for (int m=0; m<size; m++)	{
		is >> *pb;
		for (int i=0; i<Nsite; i++)	{
			rate[i] += pb->rate[i];
			for (int k=0; k<Nstate; k++)	{
				stat[i][k] += pb->Stationary[pb->Mode[i]][k];
			}
		}
	}
	
	for (int i=0; i<Nsite; i++)	{	
		rate[i] /= size;
		for (int k=0; k<Nstate; k++)	{
			stat[i][k] /= size;
		}
	}

	string logoheader = auxpath + seqheader;
	string appel = "cp " + logoheader + " " + ChainName + ".logo";
	system(appel.c_str());
	ofstream os((ChainName + ".logo").c_str(), IOS_APPEND);


	int LinePerPage = 2;
	int SitePerLine = 50;
	double MaxHeight = 3;
	double logoThreshold = 0.001;
	int numbering = 1;
	double SpaceBetweenLines = -130;
	double EntropyEpsilon = 1e-10;

	double entropy[Nsite];
	double MaxEnt = log(Nstate);
	double max = 0;
	double min = 0;

	if(style == 1)	{ // heights proportional to stat information content

		min = MaxEnt;

		for (int i=0; i<Nsite; i++)	{
			double h = 0;
			for (int j=0; j<Nstate; j++)	{
				double temp = stat[i][j];
				h += (temp<EntropyEpsilon) ? 0 : -temp * log(temp);
			}
			entropy[i] = MaxEnt -  h;

			if (max < entropy[i])	{
				max = entropy[i];
			}
			if (min > entropy[i])	{
				min= entropy[i];
			}
		}
	}

	cerr << "max : " << max << " / " <<  MaxEnt << '\n';

	if (style == 2)	{ // heights all equal

		for (int i=0; i<Nsite; i++)	{
			entropy[i] = 1;
			max = 1;
			min = 0;
		}
	}

	if (style == 3)	{ // heights proportional to rates

		min = 1;
		max = 0;
		
		for (int i=0; i<Nsite; i++)	{

			entropy[i] = - log(rate[i]);
			
			if (max < entropy[i])	{
				max = entropy[i];
			}
			if (min > entropy[i])	{
				min = entropy[i];
			}
		}
	}

	int pageNumber = 1;
	int lineNumber = 0;
	int siteNumber = 0;

	os << "%%Page: " << pageNumber << ' ' << pageNumber << '\n';
	os << "startpage\n";
	os << "startline\n";

	double array[Nstate];
	int letters[Nstate];

	for (int i=0; i<Nsite; i++)	{
		entropy[i] = MaxHeight * (entropy[i] - min) /( max - min) + 0.3;

		for (int k=0; k<Nstate; k++)	{
			array[k] = stat[i][k] *entropy[i] ;
			letters[k] = k;
		}

		for (int j=0; j<Nstate-1; j++)	{
			for (int k=Nstate-1; k>j; k--)	{
				if (array[k] < array[k-1])	{
					double temp = array[k];
					array[k] = array[k-1];
					array[k-1] = temp;

					int intemp = letters[k];
					letters[k] = letters[k-1];
					letters[k-1] = intemp;
				}
			}
		}

		if (numbering)	{
			os << "numbering {(" << i << ") makenumber} if\n";
		}
		os << "gsave\n";

		for (int k=0; k<Nstate; k++)	{
			double height = (array[k] > 0) ? array[k] : -array[k];
			if (height > logoThreshold)	{
				os << height << " (" << AminoAcids[letters[k]] << ") numchar\n";
			}
		}

		os << "grestore\n";
		os << "shift\n";

		siteNumber++;
		if (siteNumber == SitePerLine)	{
			os << "endline\n";
			os << " 0 " << SpaceBetweenLines << " translate\n";

			siteNumber = 0;
			lineNumber ++;

			if (lineNumber == LinePerPage)	{
				lineNumber = 0;
				pageNumber++;
				os << "endpage\n";
				os << "%%Page: " << pageNumber << " " << pageNumber << "\n";
				os << "startpage\n";
			}
			os << "startline\n";
		}
	}

	os << "endline\n";
	os << "endpage\n";

	os << "%%Trailer\n";
	os << "%%Pages: " << pageNumber << "\n";

}
*/
// ---------------------------------------------------------------------------
//		 PosteriorMean
// ---------------------------------------------------------------------------


double Chain::PosteriorMean(PhyloBayes* pb, int burnin, int every, int size, string ext)	{

	if (ext != "")	{
		ofstream os((ChainName + ext).c_str());
		os << size << '\n';
	}

	int Nstate = mParam->Nstate;
	int Nrr = Nstate * (Nstate - 1) / 2;

	PhyloBayes* next = mParam->GetNextState();
	pb->Clone(next);

	pb->gamma = 0;
	pb->MeanLength = 0;
	pb->ModeStatAlpha = 0;
	for (int k=0; k<mParam->Nstate; k++)	{
		pb->ModeStatCenter[k] = 0;
	}
	
	for (int j=0; j<mParam->Nnode; j++)	{
		pb->BL[j] = 0;
	}
	for (int i=0; i<Nrr; i++)	{
		pb->ModeRR[i] = 0;
	}
	
	if (pb->Nmode == 1)	{
		for (int k=0; k<mParam->Nstate; k++)	{
			pb->Stationary[0][k] = 0;
		}
	}
	pb->mLogSampling = 0;
	
	cerr << "burnin \n";
	cerr.flush();

	int n = 0;
	while ((n<burnin*every)  && RunningStatus())	{
		cerr << n << '\n';
		cerr.flush();
		mParam->MoveAll();
		n++;
	}

	cerr << "first round: posterior mean\n";
	cerr.flush();

	n = 0;
	while ((n<size)  && RunningStatus())	{

		next->UpdateLogProbs();
		if (ext != "")	{
			SavePoint(ext);
		}
			
		pb->mLogSampling += next->mLogSampling;
		pb->gamma  += next->gamma;
		pb->MeanLength += next->MeanLength;
		pb->ModeStatAlpha += next->ModeStatAlpha;
		for (int k=0; k<mParam->Nstate; k++)	{
			pb->ModeStatCenter[k] += next->ModeStatCenter[k];
		}
		
		for (int j=0; j<mParam->Nnode; j++)	{
			pb->BL[j] += next->BL[j];
		}
		for (int i=0; i<Nrr; i++)	{
			pb->ModeRR[i] += next->ModeRR[i];
		}
		
		if (pb->Nmode == 1)	{
			for (int k=0; k<mParam->Nstate; k++)	{
				pb->Stationary[0][k] += next->Stationary[0][k];
			}
		}
		for (int j=0; j<every; j++)	{
			mParam->MoveAll();
		}
		n++;
	}

	pb->mLogSampling /= size;
	pb->gamma  /= size;
	pb->MeanLength /= size;
	pb->ModeStatAlpha /= size;
	for (int k=0; k<mParam->Nstate; k++)	{
		pb->ModeStatCenter[k] /= size;
	}
	
	for (int j=0; j<mParam->Nnode; j++)	{
		pb->BL[j] /= size;
	}
	double total = 0;
	for (int i=0; i<Nrr; i++)	{
		pb->ModeRR[i] /= size;
		total += pb->ModeRR[i];
	}
	total /= Nrr;
	for (int i=0; i<Nrr; i++)	{
		pb->ModeRR[i] /= total;
	}
	for (int j=0; j<mParam->Nnode; j++)	{
		pb->BL[j] *= total;
	}
	
	if (pb->Nmode == 1)	{
		for (int k=0; k<mParam->Nstate; k++)	{
			pb->Stationary[0][k] /= size;
		}
	}
	return pb->mLogSampling;
}

	
// ---------------------------------------------------------------------------
//		 PosteriorSample
// ---------------------------------------------------------------------------


void Chain::PosteriorSample(int burnin, int every, int size, string ext)	{

	if (ext != "")	{
		ofstream os((ChainName + ext).c_str());
		os << size << '\n';
	}


	// burnin
	int i=0; 
	while ((i<burnin*every) && RunningStatus())	{	
		mParam->MoveSiteVariables();
		i++;
	}

	// sampling
	int n = 0;
	while ((n<size) && RunningStatus())	{
	
		// next->UpdateLogProbs();
		if (ext != "")	{
			SavePoint(ext);
		}

		if (n != (size-1))	{
			for (int k=0; k<every; k++)	{
				mParam->MoveSiteVariables();
			}
		}
		n++;
	}
}

// ---------------------------------------------------------------------------
//		 fixed point
// ---------------------------------------------------------------------------


double Chain::FixedPointStep(int burnin, int every, int size, string ext)	{

	if (ext != "")	{
		ofstream os((ChainName + ext).c_str());
		os << size << '\n';
	}

	int Nnode = mParam->Nnode;
	double v[Nnode];
	double r[Nnode];
	for (int j=0; j<Nnode; j++)	{
		v[j] = 0;
		r[j] = 0;
	}

	PhyloBayes* pb = mParam->GetNextState();

	// burnin
	int i=0; 
	while ((i<burnin*mParam->SaveEvery) && RunningStatus())	{	
		mParam->Move();
		MonitorCurrentPoint();
		i++;
	}

	// sampling
	int n = 0;
	while ((n<size) && RunningStatus())	{
	
		pb->ResampleSub();
		for (int j=0; j<Nnode; j++)	{
			if (! pb->tree[j].isRoot())	{
				v[j] += pb->BranchTotalSub[j];
				r[j] += pb->BranchEffLength[j];
			}
		}
		if (n != (size-1))	{
			for (int k=0; k<mParam->SaveEvery; k++)	{
				mParam->Move();
				MonitorCurrentPoint();
			}
		}
		n++;
	}

	double lengthdiff = 0;
	for (int j=0; j<Nnode; j++)	{
		if ((j != pb->root->label) && (j != pb->root->left->label)) {
			v[j] /= size;
			r[j] /= size;
			lengthdiff += fabs(v[j] / pb->BL[j] - r[j]);
			pb->BL[j] = v[j] / r[j];
			if (pb->BL[j] < mParam->LengthMin)	{
				pb->BL[j] = mParam->LengthMin;
			}
		}
	}
	mParam->GetCurrentState()->Clone(mParam->GetNextState());
	mParam->Update();

	double difflog = lengthdiff;
	return difflog;
}


// ---------------------------------------------------------------------------
//		 ThermoIntegrate
// ---------------------------------------------------------------------------


double Chain::ThermoIntegrate(int burnin, int every, int size, int orientation, ModelSwitchMode mode, int moveall)	{
		
	mParam->MSMode = Thermo;

	double initbeta = 0;
	double betastep = 1.0 / size;
	if (orientation == -1)	{
		initbeta = 1;
		betastep = -1.0 / size;
	}

	string ext;
	if (orientation == 1)	{
		ext = ".01.thermo";
	}
	else	{
		ext = ".10.thermo";
	}

	ofstream Thermo_os((ChainName + ext).c_str());

	// initialise beta or what is useful
	switch (mode)	{

		case Thermo:
		case SA:
		mParam->Beta0 = initbeta;
		break;

		case SUB:	
		mParam->SUBModelSwitch = initbeta;
		break;

		default:
		cerr << "????\n";

	}
	mParam->Update();

	// burnin
	int n=0; 
	while ((n<burnin) && RunningStatus())	{	
		if (moveall)	{
			mParam->MoveAll(0);
		}
		else	{
			mParam->MoveSiteVariables(0);
		}
		n++;
	}
	
	double total = 0;

	n = 0;
	while ((n<= size) && RunningStatus())	{	

		// save log samplings
		double temp = 0;
		double potential = 0;
		switch(mode)	{

			case SA:
			case Thermo:
			potential = mParam->GetNextState()->logSampling();
			break;

			case SUB:
			temp = mParam->SUBModelSwitch;
			mParam->SUBModelSwitch = 1;
			mParam->GetNextState()->UpdateLogProbs();
			potential = mParam->GetNextState()->logSampling();
			mParam->SUBModelSwitch = 0;
			mParam->GetNextState()->UpdateLogProbs();
			potential -= mParam->GetNextState()->logSampling();
			mParam->SUBModelSwitch = temp;
			mParam->GetNextState()->UpdateLogProbs();
			break;

			default:
			cerr << "???\n";
			
		}

		if ((n==0) || (n==size))	{
			total += 0.5 * potential;
		}
		else	{
			total += potential;
		}

		ofstream Thermo_os((ChainName + ext).c_str(), APPEND);
		Thermo_os << mParam->Beta0 << '\t' << potential << '\n';
		Thermo_os.close();

		switch (mode)	{

			case Thermo:
			mParam->Beta0 += betastep;
			break;

			case SA:
			mParam->Beta0 *= betastep;
			break;

			case RAS:
			mParam->RASModelSwitch += betastep;
			break;

			case SUB:	
			mParam->SUBModelSwitch += betastep;
			break;

			default:
			cerr << "????\n";

		}
		mParam->Update();
		int i = 0;
		while ((i<= every) && RunningStatus())	{	
			if (moveall)	{
				mParam->MoveAll(0);
			}
			else	{
				mParam->MoveSiteVariables(0);
			}
			i++;
		}
		n++;
	}
	
	total /= size;
	return total;
}


// ---------------------------------------------------------------------------
//		 MonitorCurrentPoint
// ---------------------------------------------------------------------------

void Chain::MonitorCurrentPoint()	{

	double totaltime = 0;
	for (int k=0; k<mParam->TopoMoveTypeNumber; k++)	{
		for (int j=0; j<mParam->SuperMoveTypeNumber[k]; j++)	{
			for (int i=0; i< mParam->MoveTypeNumber[k][j]; i++)	{
				for (int chain = 0 ; chain < mParam->Nchain ; chain++)	{
					totaltime += mParam->TimeArray[k][j][i][chain];
				}
			}
		}
	}

	for (int chain=0; chain<mParam->Nchain; chain++)	{

		PhyloBayes* pb = mParam->GetCurrentState(chain);

		ostringstream s;
		s << ChainName;
		if (mParam->Nchain > 1)	{
			s << "." << chain+1;
		}
		string Name = s.str();
		
		ofstream currenttree_os((Name + ".currenttree").c_str());
		if (mParam->ActivateClock)	{
			if (! mParam->NormalApprox)	{
				pb->SwapBL();
				pb->NormaliseLengths();
				pb->Phylip(currenttree_os,1,0,1,0);
				pb->DenormaliseLengths();
				pb->SwapBL();
			}
			double bl[mParam->Nnode];
			for (int j=0; j<mParam->Nnode; j++)	{
				bl[j] = pb->BL[j];
			}
			pb->SetBL();
			pb->Phylip(currenttree_os,1,0,1,0);
			for (int j=0; j<mParam->Nnode; j++)	{
				pb->BL[j] = bl[j];
			}
		}
		else	{
			pb->NormaliseLengths();
			pb->Phylip(currenttree_os,1,0,1,0);
			pb->DenormaliseLengths();
		}

		ofstream os((Name + ".trace").c_str(), APPEND);

		os << mParam->HowManySaved << '\t';
		os << mParam->Ngen << '\t';
		/*
		if (mParam->FixTopo || mParam->ActivateClock)	{
			os << mParam->HowManySaved << '\t';
		}
		else	{
			os << mParam->Ngen << '\t';
		}
		*/
		os << totaltime << '\t' << -pb->logSampling();
		if (mParam->ActivateClock)	{
			if (mParam->SeparateModelSwitch == 1)	{
				if (mParam->BLModelSwitch == 1)	 {
					os << '\t' << pb->GetGeneLength();
				}
				else	{
					os << '\t' << pb->GetGeneRigidLength();
				}
			}
			else	{
				if (mParam->BLModelSwitch == 1)	 {
					os << '\t' << pb->GetLength() * pb->LengthNormFactor();
				}
				else	{
					os << '\t' << pb->GetRigidLength() * pb->LengthNormFactor();
				}
			}
			os << '\t' << pb->Sigma;
			if ((mParam->ClockModel == CIRrigid) || (mParam->ClockModel == CIRrigid2) || (mParam->FlexClockModel == CIRflex)) 	{
				os << '\t' << pb->Theta;
			}
			if (mParam->ClockModel == AutoRegressive)	{
				os << '\t' << pb->Theta;
			}
			os << '\t' << pb->Mu;
			os << '\t' << pb->MeanRho();
			// os << '\t' << pb->VarRho();
			if (mParam->TimePrior == Dirich)	{
				os << '\t' << pb->Chi;
			}
			else if (mParam->TimePrior == BD)	{
				os << '\t' << pb->Chi << '\t' << pb->Chi2;
				/*
				pb->RefreshBDLogGG();
				os << '\t' << -pb->logAbsBDUncalibratedPrior();
				*/
			}
			if (mParam->NCalib)	{
				os << '\t' << pb->Scale;
			}
			/*
			pb->RefreshBDLogGG();
			os << '\t' << pb->logAbsBDPrior();
			*/
		}
		else	{
			// os << '\t' << pb->GetLength();
			os << '\t' << pb->GetLength() * pb->LengthNormFactor();
		}

		if (mParam->NormalApprox == No)	{
			os << '\t' << pb->gamma;
			if (mParam->MutMode)	{
				if (mParam->MutMode == 2)	{
					os << '\t' << pb->Kappa;
					for (int i=0; i<Nnuc; i++)	{
						os << '\t' << pb->RefACGT[i];
					}
				}
				else if (mParam->MutMode == 3)	{
					os << '\t' << pb->GetRRACGT_Entropy();
					os << '\t' << pb->GetMeanRRACGT();
					for (int i=0; i<Nnuc; i++)	{
						os << '\t' << pb->RefACGT[i];
					}
				}
			}
			if (! mParam->FixPconst)	{
				os << '\t' << pb->pconst;
			}
			if (mParam->NRateModeMax == mParam->Nsite + 5)	{
				os << '\t' << pb->NRateMode;
				os << '\t' << pb->GetMeanRate();
			}
			if (mParam->GeneBLMultiplier)	{
				if (mParam->LengthGammaPrior)	{
					os << '\t' << 1.0 / pb->LengthGamma;
				}
				else	{
					os << '\t' << pb->LengthGamma;
				}
				// os << '\t' << pb->GetMeanGeneBL();
			}
			if (mParam->MBL)	{
				os << '\t' << pb->LengthGamma;
				os << '\t' << pb->NMixBL;
				os << '\t' << pb->GetEffNMixBL();
			}
			if (mParam->HeteroMode == Covarion)	{
				os << '\t' << pb->xi0 << '\t' << pb->invprob;
			}
			if (mParam->FixNmode && (pb->Nmode > 1))	{
				os << '\t' << pb->GetNoccupied();
			}
			else	{
				os << '\t' << pb->Nmode;
			}
			os << '\t' << pb->GetStationaryEntropy();
			os << '\t' << pb->ModeStatAlpha;
			// os << '\t' << pb->Nmode << '\t' << pb->GetEffNmode() << '\t' << pb->GetStationaryEntropy();
			if (mParam->MSMode == CATAlpha)	{
				os << '\t' << pb->alpha;
			}
			if (!mParam->FixRR)	{
				os << '\t' << pb->GetRREntropy() ;
				os << '\t' << pb->GetMeanRR() ;
				if (mParam->RRPrior == GammaDistributed)	{
					os << '\t' << pb->RRAlpha;
				}
				if (mParam->Qmode)	{
					os << '\t' << pb->GetRefRREntropy();
				}
				if (mParam->MutMode==1)	{
					os << '\t' << pb->GetRefStationaryEntropy();
				}
			}
			os << '\t' << pb->alpha << '\t' << pb->GetModeAffEntropy();
			if (mParam->NH)	{
				if ((mParam->NH == 1) || (mParam->NH == 3))	{
					os << '\t' << pb->GetNHActiveCatNumber();
				}
				if (mParam->NH == 1)	{
					os  << '\t' << pb->NHpswitch;
				}
				if (! mParam->PopEff || mParam->NHPop)	{
					os << '\t' << pb->GetNHStatEntropy();
					os << '\t' << pb->NHStatAlpha;
				}
				if (mParam->PopEff)	{
					os << '\t' << pb->GetMeanPopEffSize() << '\t' << sqrt(pb->GetVarPopEffSize()) << '\t' << pb->GetRootPopEff() << '\t' << pb->PopAlpha << '\t' << pb->PopBeta;
				}
			}
		}
		else	{
		}
		if (mParam->ZipGTR == 4)	{
			os << '\t' << mParam->MeanOrbitSize;
		}
		else if (mParam->ZipGTR && (mParam->ZipGTR != 3))	{
			os << '\t' << pb->GetMeanOrbitSize() << '\t' << mParam->MeanOrbitSize;
		}
		/*
		if  (!mParam->ModeFastCompute)	{
			os << '\t' << mParam->meannpow;
		}
		*/
		// os << '\t' << -pb->logPrior();
		os << '\n';
		os.close();
	}

	ofstream Monitor_os((ChainName + ".monitor").c_str());
	int chain = 0;
	for (int k=0; k< mParam->TopoMoveTypeNumber; k++)	{
		for (int j=0; j<mParam->SuperMoveTypeNumber[k]; j++)	{
			for (int i=0; i< mParam->MoveTypeNumber[k][j]; i++)	{
				Monitor_os << '\t';
				string temp = MoveTypeName[mParam->MoveTypeArray[k][j][i]];
				int n = 30 - temp.length();
				for (int l=0; l<n; l++)	{
					temp += " ";
				}
				Monitor_os << temp << '\t';
				Monitor_os << ((double) ((int) mParam->TimeArray[k][j][i][chain] * 100)) / 100 << '\t';
				if (mParam->NCallArray[k][j][i][chain])	{
					Monitor_os << 100* mParam->SuccessArray[k][j][i][chain] / mParam->NCallArray[k][j][i][chain] << '\t';
				}
				else	{
					Monitor_os << '-' << '\t';
				}
				Monitor_os << '\n';
			}
		}
	}
	
	Monitor_os << '\n';
	/*
	for (int chain=1; chain<mParam->Nchain; chain++)	{
		Monitor_os << mParam->Beta[chain-1] << mParam->Beta[chain] << '\t' << ((double) mParam->AcceptedSwapArray[chain]) /  ((double) mParam->TrialSwapArray[chain]) << '\n';
	}
	*/

	Monitor_os << '\n';
	Monitor_os << '\t' << "stat inf count : " << mParam->StatInfCount << '\n';
	Monitor_os << '\t' << "inf prob count : " << mParam->LogProbInfCount << '\n';
	Monitor_os << '\t' << "sub overflow count : " << mParam->SubOverflowCount << '\n';
	Monitor_os << '\t' << "sub max count : " << mParam->ObservedUniSubNmax << '\n';
	Monitor_os << '\n';
	Monitor_os.close();
	
}


// ---------------------------------------------------------------------------
//		 Run
// ---------------------------------------------------------------------------

int Chain::Run()	{


	if (RunningStatus() && ((mParam->StopAfter == -1) || (mParam->HowManySaved < mParam->StopAfter)) && CheckDiff())	{
	while(RunningStatus() && ((mParam->StopAfter == -1) || (mParam->HowManySaved < mParam->StopAfter)))	{

		// Make one cycle
		for (int i=0; i<mParam->SaveEvery; i++)	{
			mParam->Move();
		}

		if (mParam->MSMode != None) {
			// --------------------------
			// thermodynamic integration
			// --------------------------

			// save current point

			SavePoint();

			// save log samplings
			double temp = 0;
			double tmp = 0;
			double c1 = 0;
			double c0 = 0;
			switch(mParam->MSMode)	{

				case Thermo:
				break;

				case SUB:
				temp = mParam->SUBModelSwitch;
				mParam->SUBModelSwitch = 1;
				mParam->GetNextState()->UpdateLogProbs();
				c1 = mParam->GetNextState()->logSampling();
				mParam->SUBModelSwitch = 0;
				mParam->GetNextState()->UpdateLogProbs();
				c0 = mParam->GetNextState()->logSampling();
				mParam->SUBModelSwitch = temp;
				mParam->GetNextState()->UpdateLogProbs();
				break;

				case RAS:
				temp = mParam->RASModelSwitch;
				mParam->RASModelSwitch = 1;
				mParam->GetNextState()->UpdateLogProbs();
				c1 = mParam->GetNextState()->logSampling();
				mParam->RASModelSwitch = 0;
				mParam->GetNextState()->UpdateLogProbs();
				c0 = mParam->GetNextState()->logSampling();
				mParam->RASModelSwitch = temp;
				mParam->GetNextState()->UpdateLogProbs();
				break;

				case BL:
				temp = mParam->BLModelSwitch;
				mParam->BLModelSwitch = 0;
				mParam->GetNextState()->UpdateLogProbs();
				c1 = mParam->GetNextState()->logPosterior();
				mParam->BLModelSwitch = 1;
				mParam->GetNextState()->UpdateLogProbs();
				c0 = mParam->GetNextState()->logPosterior();
				mParam->BLModelSwitch = temp;
				mParam->GetNextState()->UpdateLogProbs();
				break;

				case ClockPrior:
				temp = mParam->ClockPriorModelSwitch;
				mParam->ClockPriorModelSwitch = 1;
				c1 = mParam->GetNextState()->logClockPrior();
				mParam->ClockPriorModelSwitch = 0;
				c0 = mParam->GetNextState()->logClockPrior();
				mParam->ClockPriorModelSwitch = temp;
				break;

				case Autocorrel:
				temp = mParam->AutoCorrelModelSwitch;
				tmp = mParam->GetNextState()->AutoCorrelLogDeriv();
				mParam->AutoCorrelModelSwitch = 1;
				c1 = mParam->GetNextState()->logClockPrior() + tmp;
				mParam->AutoCorrelModelSwitch = 0;
				c0 = mParam->GetNextState()->logClockPrior();
				mParam->AutoCorrelModelSwitch = temp;
				break;

				case FlexClock:
				c1 = mParam->GetNextState()->logBLClockPrior();
				c0 = mParam->GetNextState()->logBLDeconstrainedPrior();
				break;

				case CATAlpha:
				break;

				case LengthMixGamma:
				c1 = mParam->GetNextState()->logGeneBLDerivative1();
				c0 = mParam->GetNextState()->logGeneBLDerivative0();
				break;

				case RateGamma:
				c1 = mParam->GetNextState()->logRateGammaDerivative1();
				c0 = mParam->GetNextState()->logRateGammaDerivative0();
				break;

				default:
				cerr << "non recognised thermo mode\n";
				exit(1);
				break;
			}


			if ((mParam->BurnInDone == 0) || (mParam->BurnInDone == 2))	{

				if (mParam->HowManySaved == (mParam->BurnIn))	{
					mParam->BurnInDone++;
				}
				if (mParam->HowManySaved > (mParam->BurnIn))	{
					cerr << "error in thermodynamic integration : did not catch when mParam->HowManySaved == BurnIn\n";
					cerr.flush();
					throw;
				}

			}

			else if((mParam->BurnInDone == 1) || (mParam->BurnInDone == 3))	{
	
				if ((mParam->MSMode == CATAlpha) || (mParam->MSMode == LengthMixGamma) || (mParam->MSMode == RateGamma))	{
					double currentms = 0;
					if (mParam->MSMode == CATAlpha) 	{
						currentms = mParam->GetCurrentState()->alpha;
					}
					else if (mParam->MSMode == LengthMixGamma)	{
						currentms = mParam->GetCurrentState()->LengthGamma;
					}
					else if (mParam->MSMode == RateGamma)	{
						currentms = mParam->GetCurrentState()->gamma;
					}
					ofstream Thermo_os((ChainName+ ".thermo").c_str(), IOS_APPEND);
					Thermo_os.precision(30);
					Thermo_os << currentms << '\t' << c1 << '\t' << c0 << '\n';
					Thermo_os.close();
					if (mParam->MSMode == CATAlpha) 	{
						mParam->GetCurrentState()->alpha *= mParam->BetaStep;
						mParam->GetNextState()->alpha *= mParam->BetaStep;
					}
					else if (mParam->MSMode == LengthMixGamma)	{
						mParam->GetCurrentState()->LengthGamma *= mParam->BetaStep;
						mParam->GetNextState()->LengthGamma *= mParam->BetaStep;
					}
					else if (mParam->MSMode == RateGamma)	{
						mParam->GetCurrentState()->gamma *= mParam->BetaStep;
						mParam->GetNextState()->gamma *= mParam->BetaStep;
					}
				}
				else	{
					if (mParam->Tempering)	{
						mParam->Zeta += mParam->BetaStep;
						mParam->GetCurrentState()->SetBeta(mParam->Zeta);
						mParam->GetNextState()->SetBeta(mParam->Zeta);
					}
					else	{	
						int integration  = ((int) fabs((mParam->FinalBeta - mParam->InitBeta) / mParam->BetaStep)) + 1;

						double currentms = mParam->ArcTh(mParam->Zeta);
						ofstream Thermo_os((ChainName+ ".thermo").c_str(), IOS_APPEND);
						Thermo_os.precision(30);
						Thermo_os << currentms << '\t' << c1 << '\t' << c0 << '\n';
						Thermo_os.close();

						double lastms = 0;
						double nextms = 0;
						temp = c0 - c1;
						if (mParam->BFcount == 0)	{
							nextms = mParam->ArcTh(mParam->Zeta+mParam->BetaStep);
							if (mParam->BurnInDone == 1)	{
								mParam->logBF += 0.5*fabs(nextms-currentms)*temp;
								mParam->logBFmax += fabs(nextms-currentms)*temp;
							}
							else	{
								mParam->logBFrev += 0.5*fabs(nextms-currentms)*temp;
								mParam->logBFmaxrev += fabs(nextms-currentms)*temp;
							}
						}
						else if (mParam->BFcount==integration)	{
							lastms = mParam->ArcTh(mParam->Zeta-mParam->BetaStep);
							if (mParam->BurnInDone == 1)	{
								mParam->logBF += 0.5*fabs(currentms-lastms)*temp;
								mParam->logBFmin += fabs(currentms-lastms)*temp;
							}
							else	{
								mParam->logBFrev += 0.5*fabs(currentms-lastms)*temp;
								mParam->logBFminrev += fabs(currentms-lastms)*temp;
							}
						}
						else	{
							lastms = mParam->ArcTh(mParam->Zeta-mParam->BetaStep);
							nextms = mParam->ArcTh(mParam->Zeta+mParam->BetaStep);
							if (mParam->BurnInDone == 1)	{
								mParam->logBF += 0.5 *(fabs(nextms-lastms))*temp;
								mParam->logBFmax += fabs(nextms-currentms)*temp;
								mParam->logBFmin += fabs(currentms-lastms)*temp;
							}
							else	{
								mParam->logBFrev += 0.5 *(fabs(nextms-lastms))*temp;
								mParam->logBFmaxrev += fabs(nextms-currentms)*temp;
								mParam->logBFminrev += fabs(currentms-lastms)*temp;
							}
						}

						mParam->Zeta += mParam->BetaStep;
						mParam->BFcount ++;
						if ( ((mParam->InitBeta - mParam->FinalBeta) * (mParam->Zeta - mParam->FinalBeta)) < 0)	{
							mParam->Zeta = mParam->FinalBeta;
						}
						if (mParam->BFcount > integration)	{
							mParam->BurnInDone++;
							mParam->BetaStep *= -1;
							double temp = mParam->InitBeta;
							mParam->InitBeta = mParam->FinalBeta;
							mParam->FinalBeta = temp;
							mParam->BFcount = 0;
							mParam->HowManySaved = 0;
						}
						double newbeta = mParam->ArcTh(mParam->Zeta);

						switch (mParam->MSMode)	{
				
							case RAS:
							mParam->RASModelSwitch = newbeta;
							break;

							case ClockPrior:
							mParam->ClockPriorModelSwitch = newbeta;
							break;

							case BL:
							mParam->BLModelSwitch = newbeta;
							break;

							case FlexClock:
							mParam->FlexClockModelSwitch = newbeta;
							break;

							case SUB:	
							mParam->SUBModelSwitch = newbeta;
							break;

							case Autocorrel:
							mParam->AutoCorrelModelSwitch = newbeta;
							break;

							default:
							break;
						}
						mParam->Update();
					}
				}
			}
			else	{

				//check for end of run
				if (mParam->BurnInDone == 4) 	{

					ofstream Run_os((ChainName + ".run").c_str());
					Run_os << 0 << '\n';
					Run_os.close();

					if (! mParam->Tempering)	{
					if (mParam->MSMode == ClockPrior)	{
						mParam->Update();
						double offset= mParam->GetCurrentState()->logUniNorm();
						mParam->logBFmin += offset;
						mParam->logBFmax += offset;
						mParam->logBF += offset;
						mParam->logBFminrev += offset;
						mParam->logBFmaxrev += offset;
						mParam->logBFrev += offset;
					}

					ofstream bf_os((ChainName + ".log").c_str(),APPEND);
					
					double minfor = mParam->logBFmin < mParam->logBFmax ? mParam->logBFmin : mParam->logBFmax;
					double maxfor = mParam->logBFmin > mParam->logBFmax ? mParam->logBFmin : mParam->logBFmax;
					double minrev = mParam->logBFminrev < mParam->logBFmaxrev ? mParam->logBFminrev : mParam->logBFmaxrev;
					double maxrev = mParam->logBFminrev > mParam->logBFmaxrev ? mParam->logBFminrev : mParam->logBFmaxrev;
					double min = minfor > minrev ? minrev : minfor;
					double max = maxfor < maxrev ? maxrev : maxfor;
					bf_os << "logBF in [" << min << " : " << max << "]\n";
					bf_os << "discrete error : " << ((maxfor-minfor) > (maxrev-minrev) ? maxfor-minfor:maxrev-minrev) << "\n";
					bf_os << "thermic lag    : " << fabs(mParam->logBF - mParam->logBFrev) << "\n";
					bf_os.close();
					cout << "logBF = " << mParam->logBF << "   [ " << min << " : " << max << "]\n";
					}
				}
			}

			// save parameter
			ofstream Param_os((ChainName + ".param").c_str());
			Param_os.precision(30);
			Param_os << *mParam;
			Param_os.close();
		}

		else	{

			// --------------------
			// plain MC elongation
			// --------------------

			// Save chain
			SavePoint();

			// update running status
			if ( (mParam->StopAfter != -1) && (mParam->HowManySaved >= mParam->StopAfter) )	{
				ofstream Run_os((ChainName + ".run").c_str());
				Run_os << 0 << '\n';
				Run_os.close();
			}

			if (! CheckDiff())	{
				ofstream Run_os((ChainName + ".run").c_str());
				Run_os << 0 << '\n';
				Run_os.close();
			}
		}
		// monitor chain
		MonitorCurrentPoint();	
	}
	}
	cerr << '\n';
	cerr << '\n' << ChainName << " stopped at : " << mParam->HowManySaved << " points\n";
	return 1 - RunningStatus();
}

void Chain::AutoCorrel(int burnin)	{

		ifstream Chain_is((ChainName + ".0.chain").c_str());
		int size = mParam->HowManySaved - burnin;	
		double logsamp[size];
		double length[size];
		PhyloBayes PB(mParam);

		cerr << "loading\n";
		cerr.flush();

		for (int i=0; i<burnin; i++)	{
			Chain_is >> PB;
		}

		double meanlogsamp = 0;
		double meanlength = 0;
		for (int i=0; i<size; i++)	{
			Chain_is >> PB;
			logsamp[i] = PB.logSampling();
			length[i] = PB.GetLength();
			meanlogsamp += logsamp[i];
			meanlength += length[i];
		}
		meanlogsamp /= size;
		meanlength /= size;

		cerr << "computing autocorrel\n";
		cerr.flush();
		
		double logsampac[size/2];
		double lengthac[size/2];
		for (int i=0; i<size/2; i++)	{
			double total = 0;
			for (int j=0; j<size-i; j++)	{
				total += (length[j] - meanlength) * (length[j+i] - meanlength);
			}
			total /= (size-i);
			lengthac[i] = total;

			total = 0;
			for (int j=0; j<size-i; j++)	{
				total += (logsamp[j] - meanlogsamp) * (logsamp[j+i] - meanlogsamp);
			}
			total /= (size-i);
			logsampac[i] = total;
		}

		cerr << "ouput\n";
		cerr.flush();

		ofstream os((ChainName + ".autocorrel").c_str());
		for (int i=0; i<size/2; i++)	{
			os << i << '\t' << logsampac[i] << '\t' << lengthac[i] << '\n';
		}

		cerr << "done\n";
		cerr.flush();
}



// ---------------------------------------------------------------------------
//		 ReadThermo
// ---------------------------------------------------------------------------

double Chain::ReadThermo(double min, double max, int thick)	{
	
	ifstream Thermo_is((ChainName + ".thermo").c_str());
	int integration = ((int) (fabs(mParam->InitBeta - mParam->FinalBeta) / fabs(mParam->BetaStep))) + 1;
	double maxerror = 0;
	double zeta = mParam->InitBeta;
	double totalmin = 0;
	double totalmax = 0;
	double total = 0;
	double totalerror = 0;

	for (int i=0; i<= integration; i++)	{
		double modelswitch;
		double check1 , check0;
		Thermo_is >> modelswitch >> check1 >> check0;
		double temp = check0 - check1;

		double currentms = mParam->ArcTh(zeta);
		double lastms = 0;
		double nextms = 0;
		if (i == 0)	{
			nextms = mParam->ArcTh(zeta+mParam->BetaStep);
			total += 0.5*fabs(nextms-currentms)*temp;
			totalerror += 0.5*fabs(nextms-currentms)*temp*mParam->ArcW(zeta);
			// totalerror -= temp*mParam->ArcH(zeta)*mParam->ArcH(zeta);
			totalmax += fabs(nextms-currentms)*temp;
		}
		else if (i==integration)	{
			lastms = mParam->ArcTh(zeta-mParam->BetaStep);
			total += 0.5*fabs(currentms-lastms)*temp;
			totalerror += 0.5*fabs(currentms-lastms)*temp*mParam->ArcW(zeta);
			// totalerror += temp*mParam->ArcH(zeta)*mParam->ArcH(zeta);
			totalmin += fabs(currentms-lastms)*temp;
		}
		else	{
			lastms = mParam->ArcTh(zeta-mParam->BetaStep);
			nextms = mParam->ArcTh(zeta+mParam->BetaStep);
			total += 0.5 *(fabs(nextms-lastms))*temp;
			totalerror += 0.5 *(fabs(nextms-lastms))*temp*mParam->ArcW(zeta);
			totalmax += fabs(nextms-currentms)*temp;
			totalmin += fabs(currentms-lastms)*temp;
		}
		zeta += mParam->BetaStep;
	}
	if ((mParam->MSMode == ClockPrior) && (mParam->TimePrior == BD))	{
		mParam->Update();
		total += mParam->GetCurrentState()->logUniNorm();
		totalmin += mParam->GetCurrentState()->logUniNorm();
		totalmax += mParam->GetCurrentState()->logUniNorm();
	}

	totalerror /= mParam->ArcZ() *mParam->ArcZ();
	cerr << "error1 : " << totalerror << '\n';
	if (maxerror < totalerror)	{
		maxerror = totalerror;
	}

	min = totalmin;
	if (min > totalmax)	{
		min = totalmax;
	}
	max = totalmin;
	if (max < totalmax)	{
		max = totalmax;
	}

	zeta = mParam->FinalBeta;
	totalmin = 0;
	totalmax = 0;
	total = 0;
	totalerror = 0;

	for (int i=0; i<= integration; i++)	{
		double modelswitch;
		double check1 , check0;
		Thermo_is >> modelswitch >> check1 >> check0;
		double temp = check0 - check1;

		double currentms = mParam->ArcTh(zeta);
		double lastms = 0;
		double nextms = 0;
		if (i == 0)	{
			nextms = mParam->ArcTh(zeta+mParam->BetaStep);
			total += 0.5*fabs(nextms-currentms)*temp;
			totalerror += 0.5*fabs(nextms-currentms)*temp*mParam->ArcW(zeta);
			// totalerror -= temp*mParam->ArcH(zeta)*mParam->ArcH(zeta);
			totalmax += fabs(nextms-currentms)*temp;
		}
		else if (i==integration)	{
			lastms = mParam->ArcTh(zeta-mParam->BetaStep);
			total += 0.5*fabs(currentms-lastms)*temp;
			totalerror += 0.5*fabs(currentms-lastms)*temp*mParam->ArcW(zeta);
			// totalerror += temp*mParam->ArcH(zeta)*mParam->ArcH(zeta);
			totalmin += fabs(currentms-lastms)*temp;
		}
		else	{
			lastms = mParam->ArcTh(zeta-mParam->BetaStep);
			nextms = mParam->ArcTh(zeta+mParam->BetaStep);
			total += 0.5 *(fabs(nextms-lastms))*temp;
			totalerror += 0.5 *(fabs(nextms-lastms))*temp*mParam->ArcW(zeta);
			totalmax += fabs(nextms-currentms)*temp;
			totalmin += fabs(currentms-lastms)*temp;
		}
		zeta -= mParam->BetaStep;
	}
	if ((mParam->MSMode == ClockPrior) && (mParam->TimePrior == BD))	{
		mParam->Update();
		total += mParam->GetCurrentState()->logUniNorm();
		totalmin += mParam->GetCurrentState()->logUniNorm();
		totalmax += mParam->GetCurrentState()->logUniNorm();
	}
	totalerror /= mParam->ArcZ() * mParam->ArcZ();

	cerr << "error1 : " << totalerror << '\n';
	if (maxerror < totalerror)	{
		maxerror = totalerror;
	}

	if (min > totalmin)	{
		min = totalmin;
	}
	if (min > totalmax)	{
		min = totalmax;
	}

	if (max < totalmin)	{
		max = totalmin;
	}
	if (max < totalmax)	{
		max = totalmax;
	}

	double DiscreteError = fabs(max - min);
	double SamplingError = 1.645 * sqrt(maxerror * 10 / integration);

	ofstream os((ChainName + ".thermo.summary").c_str());
	os << ChainName << '\t' << min - SamplingError << '\t' << max + SamplingError << '\n';
	os << "discrete error : " << DiscreteError << '\n';
	os << "sampling error : " << SamplingError << '\n';
	os.close();
	return 0.5 * (min + max) ;
}			

