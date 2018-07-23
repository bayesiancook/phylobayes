#include "phylo.h"
#include <list>
#include "TexTab.h"

// ---------------------------------------------------------------------------
//		 Sample
// ---------------------------------------------------------------------------


Sample::Sample(string inSampleName)	{

	SampleName = inSampleName;
	extract = 0;
	eof = 0;
	Chain_is = 0;
	Sample_is = new ifstream((SampleName+".sample").c_str());
	if (! ifstream((SampleName+".sample").c_str()))	{
		cerr << "error : non existing sample\n";
		exit(1);
	}
	mParam = new MCParameters();
	*Sample_is >> *mParam >> Burnin >> Every >> mSize >> Nchain;
	mParam->Update();
	currentPB = mParam->GetNextState();
	// currentPB = new PhyloBayes(mParam);
	skip = 0;
	Nchain = mParam->Nchain;
	currentChain = 0;
	currentIndex = 0;
	int size = GetSize();
	if (!size)	{
		cerr << "empty sample\n";
		return;
	}
}

Sample::Sample(string inChainName, int burnin, int every, int until, string Path)	{

	ChainName = inChainName;
	SampleName = inChainName + "_sample";
	eof = 0;
	extract = 1;
	Chain_is = 0;
	Sample_is = 0;
	if (! ifstream((ChainName+".param").c_str()))	{
		cerr << "error : non existing chain\n";
		cerr << ChainName << '\n';
		exit(1);
	}
	ifstream Param_is((ChainName + ".param").c_str());
	mParam = new MCParameters();
	mParam->SetPath(Path);
	Param_is >> *mParam;
	mParam->Update();
	currentPB = mParam->GetNextState();
	// currentPB = new PhyloBayes(mParam);
	skip = 0;
	Nchain = mParam->Nchain;
	if (burnin == -1)	{
		burnin = mParam->HowManySaved / 5;
	}
	if (until == -1)	{
		until = mParam->HowManySaved;
	}
	if ( (burnin > mParam->HowManySaved) || (until > mParam->HowManySaved) || (burnin > until) )	{
		cerr << "error in Sample::Sample : extraction boundaries are not correct\n";
		cerr << burnin << '\t' << every << '\t' << until << '\n';
		cerr << mParam->HowManySaved << '\n';
		exit(1);
	}
	Burnin = burnin;
	Every = every;
	mSize =  (until - Burnin) / Every;

	currentChain = 0;
	Chain_is = 0;
	OpenChainStream();
	int size = GetSize();
	if (!size)	{
		cerr << "empty sample\n";
		return;
	}

}

// ---------------------------------------------------------------------------
//		 ~Sample
// ---------------------------------------------------------------------------

Sample::~Sample()	{

	delete Sample_is;
	delete Chain_is;
	// delete currentPB;
	delete mParam;
}

// ---------------------------------------------------------------------------
//		 OpenChainStream
// ---------------------------------------------------------------------------

void Sample::OpenChainStream()	{

	delete Chain_is;
	ostringstream name;
	name << ChainName;

	Chain_is = new ifstream((name.str() + ".chain").c_str());
	if (! *Chain_is)	{
		cerr << "error in Sample::OpenChainStream\n";
		ifstream is ((name.str() + ".trace").c_str());
		if (is)	{
			cerr << "cannot find .chain file\n";
			cerr << "you may run your chain using -s option\n";
		}
		exit(1);
	}
	for (int i=0; i<Burnin; i++)	{
		*Chain_is >> *currentPB;
	}
	currentIndex = 0;
}

// ---------------------------------------------------------------------------
//		 Reset
// ---------------------------------------------------------------------------

void Sample::Reset()	{
	currentChain = 0;
	eof = 0;
	OpenChainStream();
}

// ---------------------------------------------------------------------------
//		 ToFile
// ---------------------------------------------------------------------------

void Sample::ToFile(string Name)	{

	if (Name == "")	{
		Name = SampleName;
	}
	ofstream os((Name + ".sample").c_str());
	os << *mParam;
	os << Burnin << '\t' << Every << '\t' << mSize << '\t' << Nchain << '\n';

	PhyloBayes* PB;
	for (int i=0; i<mSize; i++)	{
		PB = GetNextPB();
		os << *PB;
	}
	os.close();
	cerr << "sample written onto file : " << Name << "\n";
	cerr.flush();
	Reset();
}
	

// ---------------------------------------------------------------------------
//		 GetNextPB
// ---------------------------------------------------------------------------

PhyloBayes* Sample::GetNextPB()	{

	if (eof)	{
		cerr << "error in Sample::GetNextPB() : reading past end of file\n";
		exit(1);
	}
	if (! extract)	{
		*Sample_is >> *currentPB;
		currentIndex ++;
		if (currentIndex == mSize)	{
			currentIndex = 0;
			currentChain++;
			if (currentChain == Nchain)	{
				eof = 1;
			}

		}
	}
	else	{
		if (skip)	{
			for (int i=0; i<Every-1; i++)	{
				*Chain_is >> *currentPB;
			}
			skip = 0;
		}
		*Chain_is >> *currentPB;
		currentIndex++;
		if (currentIndex == mSize)	{
			currentChain++;
			if (currentChain == Nchain)	{
				eof = 1;
			}
			else	{
				OpenChainStream();
			}
		}
		else	{
			skip = 1;
		}
	}
	return currentPB;
}

// ---------------------------------------------------------------------------
//		 FiniteTimeEntropy
// ---------------------------------------------------------------------------

void Sample::ReadFiniteTimeEntropy(double timemax, int N)	{

	int size = GetSize();

	double* array = new double[N];
	for (int n=0; n<N; n++)	{
		array[n] = 0;
	}

	for (int i=0; i<size; i++)	{
		cerr << '.';
		PhyloBayes* pb = GetNextPB();
		pb->Update();
		for (int n=0; n<N; n++)	{
			double time = (n+1) * timemax / N;
			array[n] += pb->GetFiniteTimeEntropy(time);
		}
	}
	cerr << '\n';

	ofstream os((SampleName + ".finitetimeentropy").c_str());
	for (int n=0; n<N; n++)	{
		array[n] /= size;
		double time = (n+1) * timemax / N;
		os << time << '\t' << array[n] << '\t' << exp(array[n]) << '\n';
	}

	cerr << '\n';
	cerr << "finite-time entropies in " << SampleName << ".finitetimeentropy\n";
	cerr << '\n';
}	

// ---------------------------------------------------------------------------
//		 SiteLogLikelihood()
// ---------------------------------------------------------------------------

void Sample::ReadSiteLogLikelihood()	{

	int size = GetSize();
	ofstream os((SampleName + ".sitelogl").c_str());
	for (int i=0; i<size; i++)	{
		cerr << '.';
		PhyloBayes* pb = GetNextPB();
		pb->Update();
		for (int j=0; j<mParam->Nsite; j++)	{
			os << -pb->mSiteLogSampling[j] << '\t';
		}
		os << '\n';
	}
	cerr << '\n';
	cerr << "site log likelihoods in " << SampleName << ".sitelogl\n";
	cerr << '\n';
}	

// ---------------------------------------------------------------------------
//		 RelativeRates
// ---------------------------------------------------------------------------

void Sample::ReadSummedLogLikelihood()	{

	int size = GetSize();
	double meanlogl = 0;
	for (int i=0; i<size; i++)	{
		cerr << '.';
		PhyloBayes* pb = GetNextPB();
		pb->Update();
		meanlogl += pb->GetSumModeLogSampling();
	}
	cerr << '\n';
	meanlogl /= size;
	ofstream os((SampleName + ".lnl").c_str());
	os << "mean lnL : " << meanlogl << '\n';
	cerr << '\n';
	cerr << "mean lnL : " << meanlogl << '\n';
	cerr << '\n';
}	

// ---------------------------------------------------------------------------
//		 RelativeRates
// ---------------------------------------------------------------------------

void Sample::MeanRR()	{

	int size = GetSize();
	int Nstate = mParam->Nstate;
	cerr << "Nstate " << Nstate << '\n';
	int Nrr = Nstate * (Nstate - 1) / 2;
	double meanrr[Nrr];
	for (int i=0; i<Nrr; i++)	{
		meanrr[i] = 0;
	}
	for (int i=0; i<size; i++)	{
		PhyloBayes* pb = GetNextPB();
		double total = 0;
		for (int k=0; k<Nrr; k++)	{
			total += pb->ModeRR[k];
		}
		total /= Nrr;
		for (int k=0; k<Nrr; k++)	{
			meanrr[k] += pb->ModeRR[k] / total;
		}
	}
	for (int k=0; k<Nrr; k++)	{
		meanrr[k] /= size;
	}
	/*
	double total = 0;
	for (int i=0; i<Nrr; i++)	{
		total += meanrr[i];
	}
	*/
	ofstream os((SampleName + ".meanrr").c_str());
	for (int k=0; k<Nstate; k++)	{
		os << mParam->Alphabet[k] << ' ';
	}
	os << '\n';
	os << '\n';
	/*
	for (int i=0; i<Nrr; i++)	{
		meanrr[i] /= total;
		os << meanrr[i] << '\n';
	}
	*/
	for (int i=0; i<Nstate; i++)	{
		for (int j=i+1; j<Nstate; j++)	{
			// meanrr[(2 * Nstate - i - 1) * i / 2 + j - i - 1] /= total;
			// os << meanrr[(2 * Nstate - i - 1) * i / 2 + j - i - 1];
			// os << '\n';
			os << mParam->Alphabet[i] << '\t' << mParam->Alphabet[j] << '\t' << meanrr[(2 * Nstate - i - 1) * i / 2 + j - i - 1]  << '\n';
		}
	}
	cerr << "mean relative rates in " << SampleName << ".meanrr\n";
}	

// ---------------------------------------------------------------------------
//		CATGTRPoint 
// ---------------------------------------------------------------------------

void Sample::ReadCATGTRPoint()	{

	int size = GetSize();
	int Nstate = mParam->Nstate;

	for (int i=0; i<size; i++)	{
		PhyloBayes* pb = GetNextPB();
		ostringstream s;
		s << SampleName << "_" << i;
		ofstream rros((s.str()+ ".rr").c_str());
		rros << Nstate * (Nstate -1) / 2 << '\t';
		for (int k=0; k<Nstate * (Nstate -1)/2; k++)	{
			rros << pb->ModeRR[k] << '\t';
		}
		rros << '\n';

		ofstream tos ((s.str() + ".tree").c_str());
		pb->Phylip(tos,1,0,1,0);

		ofstream os ((s.str() + ".mix").c_str());
		int Nmode = pb->Nmode;
		os << Nmode << '\n';
		double* weight = new double[Nmode];
		pb->UpdateSiteNumber();
		double total = 0;
		for (int j=0; j<Nmode; j++)	{
			weight[j] = Random::sGamma(1 + pb->SiteNumber[j]);
			total += weight[j];
		}
		os << Nmode << '\t';
		for (int j=0; j<Nmode; j++)	{
			weight[j] /= total;
			os << weight[j] << '\t';
		}
		os << '\n';
		
		for (int j=0; j<pb->Nmode; j++)	{
			os << Nstate << '\t';
			for (int k=0; k<Nstate; k++)	{
				os << pb->Stationary[j][k] << '\t';
			}
			os << '\n';
		}
		
		ofstream aos((s.str() + ".alpha").c_str());
		aos << pb->GetGamma() << '\n';

	}
}	

// ---------------------------------------------------------------------------
//		 CV
// ---------------------------------------------------------------------------

double Sample::CV(string testfile)	{

	if (mParam->NRateModeMax == mParam->Nsite + 5)	{
		cerr << "error in readcv : cannot compute cross validation score under the -ratecat model\n";
		cerr << "should use discrete gamma distribution (-dgam)\n";
		exit(1);
	}

	if (mParam->RASModelSwitch == 1)	{
		cerr << "error in readcv : cannot compute cross validation score under the continuous gamma model (-cgam)\n";
		cerr << "should use discrete gamma distribution (-dgam)\n";
		exit(1);
	}

	MCParameters* testParam = new MCParameters;
	testParam->ReadDataFromFile(testfile);
	if (mParam->Recoding)	{
		cerr << "recoding\n";
		testParam->Recoding = mParam->Recoding;
		testParam->RecodingFile = mParam->RecodingFile;
		testParam->LoadRecoding();
		cerr << "ok\n";
	}

	testParam->RegisterWithData();
	testParam->SimpleSampling = No;
	mParam->SimpleSampling = No;

	int size = GetSize();
	double logsamp[size];
	double min = 0;
	for (int i=0; i<size; i++)	{
		cerr.flush();
		PhyloBayes* pb = GetNextPB();
		cerr << '.';
		if (pb->Nmode > 1)	{
			mParam->SumOverModes = Yes;
			if (!pb->mModeSiteLogSampling)	{
				pb->mModeSiteLogSampling = new double*[mParam->NmodeMax];
				for (int i=0; i<mParam->NmodeMax; i++)	{
					pb->mModeSiteLogSampling[i] = new double[mParam->Nsite];
				}
			}
			pb->UpdateSiteNumber();
			pb->ResampleModeWeight();
		}
		mParam->SwapData(testParam);
		pb->SetData(-1);
		pb->Update();
		logsamp[i] = pb->mLogSampling;
		if ((!i) || (min > logsamp[i]))	{
			min = logsamp[i];
		}
		mParam->SwapData(testParam);
		pb->SetData(-1);
	}
	cerr << "\n";
	cerr.flush();
	double total = 0;
	for (int i=0; i<size; i++)	{
		total += exp(min - logsamp[i]);
	}
	return min - log(total);
}



// ---------------------------------------------------------------------------
//		 Check
// ---------------------------------------------------------------------------

void Sample::Check(PhyloBayes* truepb, int verbose, int topo, int ncat, int* truehisto, int* tothisto, double& palpha, double& pgamma, double& plength, double& pstatent, double& prrent)	{

	int size = mSize;
	double truealpha = truepb->alpha;
	double truegamma = truepb->gamma;
	double truelength = truepb->GetLength();
	int truenmode = truepb->Nmode;
	double truestatent = truepb->GetStationaryEntropy();
	double truerrent = truepb->GetRREntropy();

	int alpha= 0;
	int gamma= 0;
	int length= 0;
	int statent = 0;
	int rrent = 0;
	int nmode = 0;

	double meanalpha = 0;
	double meangamma = 0;
	double meanlength = 0;
	double meannmode = 0;
	double meanstatent = 0;
	double meanrrent = 0;
	double varalpha = 0;
	double vargamma = 0;
	double varlength = 0;
	double varnmode = 0;
	double varstatent = 0;
	double varrrent = 0;

	TaxaParameters* taxaparam = 0;
	TreeList* mTreeList = 0;
	if (topo)	{
		taxaparam = new TaxaParameters(mParam);
		mTreeList = new TreeList(taxaparam,size);
	}

	for (int i=0; i<size; i++)	{
		PhyloBayes* pb = GetNextPB();
		double temp = pb->alpha;
		meanalpha += temp;
		varalpha += temp*temp;
		if (temp < truealpha) alpha++;
		temp = pb->gamma;
		meangamma += temp;
		vargamma += temp*temp;
		if (temp < truegamma) gamma++;
		temp = pb->GetLength();
		meanlength += temp;
		varlength += temp * temp;
		if (temp < truelength) length++;
		temp = pb->Nmode;
		meannmode += temp;
		varnmode += temp * temp;
		if (temp < truenmode) nmode++;
		temp = pb->GetStationaryEntropy();
		meanstatent += temp;
		varstatent += temp * temp;
		if (temp < truestatent)	statent++;
		temp = pb->GetRREntropy();
		meanrrent += temp;
		varrrent += temp * temp;
		if (temp < truerrent) rrent++;
		if (topo)	{
			mTreeList->mTreeArray[i] = new Tree(pb);
		}
	}
	meanalpha /= size;
	varalpha /= size;
	varalpha -= meanalpha * meanalpha;
	meangamma /= size;
	vargamma /= size;
	vargamma -= meangamma * meangamma;
	meanlength /= size;
	varlength /= size;
	varlength -= meanlength * meanlength;
	meannmode /= size;
	varnmode /= size;
	varnmode -= meannmode * meannmode;
	meanstatent /= size;
	varstatent /= size;
	varstatent -= meanstatent * meanstatent;
	meanrrent /= size;
	varrrent /= size;
	varrrent -= meanrrent * meanrrent;

	if (topo)	{
		Tree truetree(truepb);
		truetree.Trichotomise();
		BipartitionList truebplist(&truetree);
		cerr << "true ok\n";
		BipartitionList bplist(mTreeList);	

		for (int i=0; i<truebplist.GetSize(); i++)	{
			Bipartition bp = truebplist.GetBipartition(i);
			int cat = 0;
			int k=0;
			while ((k<bplist.GetSize()) && (bp != bplist.GetBipartition(k))) k++;
			if (bplist.GetProb(k) != 1)	{
			if (k < GetSize())	{
				cat = (int) (bplist.GetProb(k) * ncat);
				if (cat == ncat) cat--;
			}
			truehisto[cat]++;
			}
		}
		for (int i=0; i<bplist.GetSize(); i++)	{
			if (bplist.GetProb(i) != 1)	{
			int cat = (int) (ncat * bplist.GetProb(i));
			if (cat == ncat) cat--;
			tothisto[cat] ++;
			}
		}
			
	}

	if (verbose)	{
		cout << "true\t" << truealpha << '\t' << truenmode << '\t' << truegamma << '\t' << truelength << '\t' << truestatent << '\t' << truerrent << '\n';	
		cout << "mean\t" << meanalpha << '\t' << meannmode << '\t' << meangamma << '\t' << meanlength << '\t' << meanstatent << '\t' << meanrrent << '\n';	
		cout << "stderr\t" << sqrt(varalpha) << '\t' << sqrt(varnmode) << '\t' << sqrt(vargamma) << '\t' << sqrt(varlength) << '\t' << sqrt(varstatent) << '\t' << sqrt(varrrent) << '\n';	
	}
	cout << size << '\t' << alpha << '\t' << nmode << '\t' << gamma << '\t' << length << '\t' << statent << '\t' << rrent << '\n';

	cout.flush();

	// palpha = 100 * ((double) alpha) / size;
	palpha = 100 * ((double) nmode) / size;
	pgamma = 100 * ((double) gamma) / size;
	plength = 100 * ((double) length) / size;
	pstatent = 100 * ((double) statent) / size;
	prrent = 100 * ((double) rrent) / size;

	delete mTreeList;
}

// ---------------------------------------------------------------------------
//		 Read
// ---------------------------------------------------------------------------

double gammadensity(double rate, double gamma)	{
	double temp = Random::logGamma(gamma) + gamma*(rate - log(gamma)) + (1 - gamma) * log(rate);
	if (temp < 10)	{
		return exp(-temp);
	}
	return 0;
}

int Sample::ReadBench()	{

	int Nsite = mParam->Nsite;
	int size = mSize;
	double meanstatent[Nsite];
	for (int i=0; i<Nsite; i++)	{
		meanstatent[i] = 0;
	}
	list<double> lengthlist;
	list<double> gammalist;
	list<double> ncomplist;
	list<double> statentlist;

	for (int i=0; i<size; i++)	{

		cerr << '.';

		PhyloBayes* pb = GetNextPB();

		lengthlist.push_back(pb->GetLength() * pb->LengthNormFactor());
		gammalist.push_back(pb->gamma);
		ncomplist.push_back(pb->Nmode);
		statentlist.push_back(pb->GetStationaryEntropy());

		for (int l=0; l<Nsite; l++)	{
			meanstatent[l] += pb->GetSiteStationaryEntropy(l);
		}
	}
	cerr << '\n';

	ofstream stos((SampleName + ".statent").c_str());
	for (int i=0; i<Nsite; i++)	{
		meanstatent[i] /= size;
		stos << meanstatent[i] << '\n';
	}

	ofstream os((SampleName + ".tex").c_str());
	os << textabentry(lengthlist,true,false,true) << '\n';
	os << textabentry(gammalist,true,false,true) << '\n';
	os << textabentry(ncomplist,true,false,true) << '\n';
	os << textabentry(statentlist,true,false,true) << '\n';

	return 1;

}

void Sample::ReadModeProfiles()	{

	int nmode = mParam->GetCurrentState()->Nmode;
	int Nstate = mParam->Nstate;
	int size = mSize;

	double** modestat = new double*[nmode];
	for (int j=0; j<nmode; j++)	{
		modestat[j] = new double[Nstate];
		for (int k=0; k<Nstate; k++)	{
			modestat[j][k] = 0;
		}
	}
	double* modeweight = new double[nmode];
	for (int j=0; j<nmode; j++)	{
		modeweight[j] = 0;
	}

	for (int i=0; i<size; i++)	{

		cerr << '.';

		PhyloBayes* pb = GetNextPB();
		PhyloBayes& PB = *pb;
		PB.UpdateSiteNumber();

		if (pb->Nmode != nmode)	{
			cerr << "error in readmodeprofiles: number of categories of the mixture is not fixed\n";
			exit(1);
		}
		for (int j=0; j<nmode; j++)	{
			modeweight[j] += pb->SiteNumber[j];
			for (int k=0; k<Nstate; k++)	{
				modestat[j][k] += pb->Stationary[j][k];
			}
		}
	}
	cerr << '\n';

	ofstream os((SampleName + ".modeprofiles").c_str());
	os << nmode << '\t' << Nstate << '\n';
	for (int j=0; j<nmode; j++)	{
		modeweight[j] /= mParam->Nsite * size;
		os << modeweight[j];
		for (int k=0; k<Nstate; k++)	{
			modestat[j][k] /= size;
			os << '\t' << modestat[j][k];
		}
		os << '\n';
	}

	cerr << "mode profiles in " << SampleName << ".modeprofiles\n";
	cerr << '\n';
}

int Sample::Read(int rates, int modes, int sitestat, double cutoff, int ncat, int ps)	{

	int cons = 1;
	int treelist = 1;

	int Nsite = mParam->Nsite;
	int Nnode = mParam->Nnode;
	int Nstate = mParam->Nstate;
	int size = mSize;

	double Stationary[Nsite][Nstate];
	for (int i=0; i<Nsite; i++)	{
		for (int j=0; j<Nstate; j++)	{
			Stationary[i][j] = 0;
		}
	}

	double AbsRate[Nsite];
	double AbsRateError[Nsite];
	for (int i=0; i<Nsite; i++)	{
		AbsRate[i] = 0;
		AbsRateError[i] = 0;
	}
	double minrate = 1000;
	double maxrate = 0;

	double Nmode = 0;
	double NmodeError = 0;
	double Alpha = 0;
	double AlphaError = 0;

	double LogSampling = 0;
	double LogSamplingError = 0;
	double TotalLength = 0;
	double TotalLengthError = 0;
	double Gamma = 0;
	double GammaError = 0;
	double pconst = 0;
	double pconstError = 0;

	double Length[Nnode];
	double LengthError[Nnode];
	for (int i=0; i<Nnode; i++)	{
		Length[i] = 0;
		LengthError[i] = 0;
	}

	int* ModeNumber = new int[size];
	int  ModeOccupancyHisto[Nsite];
	int  modehisto[Nsite];
	for (int i=0; i<Nsite; i++)	{
		ModeOccupancyHisto[i] = 0;
		modehisto[i] = 0;
	}

	double ModeStatCenter[Nstate];
	for (int k=0; k<Nstate; k++)	{
		ModeStatCenter[k] = 0;
	}

	double ModeStatAlpha = 0;
	double ModeStatAlphaError = 0;

	double MeanLength = 0;
	double MeanLengthError = 0;
	
	TaxaParameters* taxaparam = new TaxaParameters(mParam);
	TreeList* mTreeList = new TreeList(taxaparam,size);

	// Reading and writing
	for (int i=0; i<size; i++)	{

		cerr << '.';

		PhyloBayes* pb = GetNextPB();
		PhyloBayes& PB = *pb;
		PB.UpdateSiteNumber();
		pb->NormaliseLengths();
		pb->ReverseSetBL();

		double temp;
		int tmp;
		
		temp = PB.ModeStatAlpha;
		ModeStatAlpha += temp;
		ModeStatAlphaError += temp * temp;
		for (int k=0; k<Nstate; k++)	{
			ModeStatCenter[k] += PB.ModeStatCenter[k];
		}
		tmp = PB.GetModeNumber();
		ModeNumber[i] = tmp;
		Nmode += tmp;
		NmodeError += tmp*tmp;
		modehisto[tmp] ++;

		for (int j=0; j<PB.GetModeNumber(); j++)	{
			ModeOccupancyHisto[PB.GetSiteNumber(j)]++;
		}

		temp = PB.GetAlpha();
		Alpha += temp;
		AlphaError += temp*temp;

		temp = PB.logSampling();
		LogSampling += temp;
		LogSamplingError += temp*temp;

		// rates
		temp = PB.GetGamma();
		Gamma += temp;
		GammaError += temp * temp;

		temp = PB.GetPconst();
		pconst += temp;
		pconstError += temp * temp;

		double length = PB.GetLength();
		TotalLength += length;
		TotalLengthError += length*length;

		for (int j=0; j<Nnode; j++)	{
			double temp = PB.GetBranchLength(j);
			Length[j] += temp;
			LengthError[j] += temp * temp;
		}

		temp = PB.MeanLength;
		MeanLength += temp;
		MeanLengthError += temp * temp;

		if (rates)	{
			if (mParam->ActivateSumOverRateModes)	{
				PB.Update();
				PB.UpdateRateModeSiteLogSampling();
				PB.ComputeMeanPosteriorRates();
			}

			for (int l=0; l<Nsite; l++)	{
				temp = PB.rate[l] * length;
				AbsRate[l] += temp;
				AbsRateError[l] += temp*temp;

				if (maxrate < temp)	{
					maxrate = temp;
				}
				if (minrate > temp)	{
					minrate = temp;
				}
			}
		}

		if (modes || sitestat)	{
			/*
			int threshold = 1;
			int Nmode = PB.GetModeNumber();

			int already[Nmode];
			for (int j=0; j<Nmode; j++)	{
				already[j] = 0;
			}

			int j= 0;
			while (j<Nmode)	{

				// look for the most populated mode
				int mode = 0;
				int sitenumber = 0;
				for (int k = 0; k<Nmode; k++)	{
					if (! already[k])	{
						if (sitenumber < PB.GetSiteNumber(k))	{
							sitenumber = PB.GetSiteNumber(k);
							mode = k;
						}
					}
				}

				already[mode] = 1;

				if (sitenumber >= threshold)	{
					j ++;
				}
				else	{
					j = Nmode;
				}

			}
			*/

			// site stationaries
			for (int j=0; j<Nsite; j++)	{
				for (int k=0; k<mParam->Nstate; k++)	{
					Stationary[j][k] += PB.GetSiteStationary(j,k);
				}
			}
		}

		if (cons || treelist)	{
			// PB.root->rootAt(&PB.tree[0]);
			mTreeList->mTreeArray[i] = new Tree(&PB);
		}

	}
	cerr << '\n';

	// computing averages and standard errors
	for (int i=0; i<Nsite; i++)	{
		AbsRate[i] /= size;
		AbsRateError[i] /= size;
		AbsRateError[i] -= AbsRate[i] * AbsRate[i];

		for (int j=0; j<mParam->Nstate; j++)	{
			Stationary[i][j] /= size;
		}
	}

	MeanLength /= size;
	MeanLengthError /= size;
	MeanLengthError -= MeanLength * MeanLength;

	ModeStatAlpha /= size;
	ModeStatAlphaError /= size;
	ModeStatAlphaError -= ModeStatAlpha * ModeStatAlpha;
	ofstream diros((SampleName + ".dirweight").c_str());
	for (int k=0; k<Nstate; k++)	{
		ModeStatCenter[k] /= size;
		diros << ModeStatAlpha * ModeStatCenter[k] << '\n';
	}

	Nmode /= size;
	NmodeError /= size;
	NmodeError -= Nmode * Nmode;
	Alpha /= size;
	AlphaError /= size;
	AlphaError -= Alpha * Alpha;
	Gamma /= size;
	GammaError /= size;
	GammaError -= Gamma * Gamma;
	pconst /= size;
	pconstError /= size;
	pconstError -= pconst * pconst;

	TotalLength /= size;
	TotalLengthError /= size;
	TotalLengthError -= TotalLength * TotalLength;
	LogSampling /= size;
	LogSamplingError /= size;
	LogSamplingError -= LogSampling * LogSampling;

	for (int i=0; i<Nnode; i++)	{
		Length[i] /= size;
		LengthError[i] /= size;
		LengthError[i] -= Length[i] * Length[i] ;
	}

	// writing averages and standard errors
	ofstream Average_os((SampleName + ".general").c_str());
	Average_os << '\n';
	Average_os << "                     mean\tstandard error" << '\n' << '\n';
	Average_os << "-log likelihood      " << LogSampling << '\t' << sqrt(LogSamplingError) << '\n';
	Average_os << "Length               " << TotalLength << '\t' << sqrt(TotalLengthError) << '\n';
	Average_os << "rate alpha param     " << Gamma << '\t' << sqrt(GammaError) << '\n';
	Average_os << "mode number          " << Nmode << '\t' << sqrt(NmodeError) << '\n';
	Average_os << "mode alpha param     " << Alpha << '\t' << sqrt(AlphaError) << '\n';
	Average_os << "mean length     " << MeanLength << '\t' << sqrt(MeanLengthError) << '\n';
	Average_os << "mode stat alpha     " << ModeStatAlpha << '\t' << sqrt(ModeStatAlphaError) << '\n';
	Average_os << '\n';

	//writing mean rates
	if (rates)	{
		ofstream Rate_os((SampleName + ".rate").c_str());
		Rate_os << "site\trate\tstd error\n\n";
		for (int i=0; i<Nsite; i++)	{
			Rate_os << i+1 << '\t' << AbsRate[i] << '\t' << sqrt(AbsRateError[i]) << '\n';
		}
		cerr << "site specific rates in \t" << SampleName << ".rate\n";

		/*
		Reset();
		int histo[ncat];
		for (int cat=0; cat<ncat; cat++)	{
			histo[cat] = 0;
		}

		for (int i=0; i<size; i++)	{
			PhyloBayes* pb = GetNextPB();
			// double length = pb->GetLength();
			pb->MakeRates();
			double total = 0;
			for (int j=0; j<mParam->Nsite; j++)	{
				double temp = pb->unirate[j];
				total += temp;
			}
			total /= mParam->Nsite;
			for (int j=0; j<mParam->Nsite; j++)	{
				double temp = pb->rate[j] / total;
				int g = (int) (ncat * (temp - minrate) / (maxrate - minrate));
				if (g==ncat)	{
					g--;
				}
				histo[g]++;
			}
		}
		ofstream os((SampleName+ ".ratehisto").c_str());
		for (int cat=0; cat<ncat; cat++)	{
			os << minrate + cat * (maxrate-minrate) / ncat << '\t' << ((double) histo[cat]) / size / mParam->Nsite *ncat / (maxrate - minrate) << '\t' << gammadensity(minrate + cat*(maxrate-minrate),Alpha) << '\n';
		}
		cerr << "rate histogram in  \t" << SampleName << ".ratehisto\n";
		*/
	}

	if (modes)	{
		// writing down mode histogram
		ofstream ModeHisto_os((SampleName +  ".modehisto").c_str());
		int i = Nsite - 1;
		while (! modehisto[i])	{
			i--;
		}
		for (int k=0; k<=i; k++)	{
			ModeHisto_os << k << '\t' << modehisto[k] << '\n';
		}
		ModeHisto_os << i << '\t' << 0 << '\n';
		cerr << "mode histogram in  \t" << SampleName << ".modehisto\n";

		ofstream mos((SampleName + ".modenumber").c_str());
		for (int i=0; i<size; i++)	{
			mos << ModeNumber[i] << '\n';
		}
		cerr << "mode number distribution in \t" << SampleName << ".modenumber\n";
	}

	//writing site stats
	if (sitestat)	{
		ofstream vos((SampleName + ".varsitestat").c_str());
		ofstream SiteStat_os((SampleName + ".sitestat").c_str());
		/*
		SiteStat_os << "site" << '\t';
		for (int j=0; j<mParam->Nstate; j++)	{
			SiteStat_os << mParam->Alphabet[j] << '\t' ;
		}
		SiteStat_os << '\n';
		*/
		SiteStat_os << Nsite << '\t' << mParam->Nstate << '\n';
		for (int i=0; i<Nsite; i++)	{
			SiteStat_os << i+1;
			double var = 0;
			for (int j=0; j<mParam->Nstate; j++)	{
				SiteStat_os << '\t' << Stationary[i][j];
				var += Stationary[i][j] * Stationary[i][j];
			}
			SiteStat_os << '\n';
			var /= mParam->Nstate;
			vos << i+1 << '\t' << var << '\n';
		}
		cerr << "site stationaries in \t" << SampleName << ".sitestat\n";
	}
	if (sitestat)	{
		double freq[Nsite][Nstate];
		for (int i=0; i<Nsite; i++)	{
			int tot = 0;
			for (int j=0; j<Nstate; j++)	{
				tot += mParam->SiteEmpiricalCount[i][j];
			}
			for (int j=0; j<Nstate; j++)	{
				freq[i][j] = ((double) mParam->SiteEmpiricalCount[i][j]) / tot ;
			}
		}
		ofstream SiteStat_os((SampleName + ".empsitestat").c_str());
		SiteStat_os << "site" << '\t';
		for (int j=0; j<mParam->Nstate; j++)	{
			SiteStat_os << mParam->Alphabet[j] << '\t' ;
		}
		SiteStat_os << '\n';
		for (int i=0; i<Nsite; i++)	{
			SiteStat_os << i+1;
			for (int j=0; j<mParam->Nstate; j++)	{
				SiteStat_os << '\t' << freq[i][j];
			}
			SiteStat_os << '\n';
		}
		cerr << "empirical stationaries in \t" << SampleName << ".empsitestat\n";
	}

	if (treelist)	{
		mTreeList->SetParameters();
		ofstream TreeList_os( (SampleName + ".treelist").c_str());
		mTreeList->WriteToStream(TreeList_os);
		cerr << "tree list in\t\t" << SampleName << ".treelist\n";
	}

	if (cons)	{

		mTreeList->SetParameters();
		BipartitionList* blist = new BipartitionList(mTreeList);
		ofstream blist_os((SampleName + ".bplist").c_str());
		blist->WriteToStream(blist_os);
		blist_os.close();
	
		Consensus* cons = new Consensus(blist,cutoff);
		ofstream cons_os((SampleName + ".con.tre").c_str());
		cons->Trichotomise();
		cons->Phylip(cons_os, 1, 1, 1, 0);

		cerr << "bipartition list in\t" << SampleName << ".bplist\n";
		cerr << "consensus in\t\t" << SampleName << ".con.tre\n";
		if (ps)	{
			cons->ToPS(SampleName + ".con.tre", 12, 20, 1, 1, 1, 0);
			cerr << "postscript in\t\t" << SampleName << ".con.tre.ps\n";
		}
		delete cons;
	}
	cerr << "summary statistics in\t" << SampleName << ".general\n";
	cerr << '\n';
	delete mTreeList;
	delete taxaparam;

	return 1;

}

// ---------------------------------------------------------------------------
//		 ClusterModes
// ---------------------------------------------------------------------------

int Sample::ClusterModes(int SizeThreshold, double ClusteringRange, int ps)	{

	int NClusterMax = 1000;

	int Nstate = mParam->Nstate;
	int size = GetSize();

	// ---------------------------------------
	// creating clustering related arrays :
	// ---------------------------------------
	
	// NCluster : number of identified clusters
	// clusterTable(i,j) : frequency of amino acid j in cluster i
	// ClusterWeight[i] : weight of cluster i (mean number of sites affiliated to modes aggregated to this cluster)
	// ClusterSize[i] : size of cluster i (number of modes aggregated to this cluster)
	
	double  clusterTable[NClusterMax][Nstate];
	for (int i=0; i<NClusterMax; i++)	{
		for (int j=0; j<Nstate; j++)	{
			clusterTable[i][j] = 0;
		}
	}
	int NCluster = 0;
	
	int* ClusterWeight = new int[NClusterMax];
	for (int i=0; i<NClusterMax; i++)	{
		ClusterWeight[i] = 0;
	}

	int* ClusterSize = new int[NClusterMax];
	for (int i=0; i<NClusterMax; i++)	{
		ClusterSize[i] = 0;
	}

	int** IsAffiliated = new int*[NClusterMax];
	for (int i=0; i<NClusterMax; i++)	{
		IsAffiliated[i] = new int[size];
		for (int j=0; j<size; j++)	{
			IsAffiliated[i][j] = 0;
		}
	}
	int ClusterStability[NClusterMax];

	double** SiteAff = new double*[mParam->Nsite];
	for (int i=0; i<mParam->Nsite; i++)	{
		SiteAff[i] = new double[NClusterMax];
		for (int j=0; j<NClusterMax; j++)	{
			SiteAff[i][j] = 0;
		}
	}

	// -------------------------
	// reading sample file
	// -------------------------
	
	for (int i=0; i<size; i++)	{
		cerr << '.';
		cerr.flush();
		PhyloBayes* pb = GetNextPB();
		PhyloBayes& PB = *pb;
		PB.UpdateSiteNumber();
	
		for (int j=0; j<PB.GetModeNumber(); j++)	{
			
			double* stat = PB.Stationary[j];
			int sitenumber = PB.SiteNumber[j];

			int affiliated = -1;

			for (int k=0; k<NCluster; k++)	{
				// compare current mode to kth cluster
				// if close enough :
				// 	if not yet affiliated, then affiliated to cluster k
				// 	else affiliation cluster and cluster k are merged
				// 		and cluster k is flagged for elimination

				double kl = 0;
				double quad = 0;
				if (ClusterSize[k] == 0)	{
					cerr << "error while aggregating clusters : found an empty cluster\n";
					exit(1);
				}

				for (int l=0; l<Nstate; l++)	{
					kl += clusterTable[k][l]/ClusterWeight[k] * (log(clusterTable[k][l])/ClusterWeight[k] - log(stat[l]));
					double temp = clusterTable[k][l]/ClusterWeight[k] - stat[l];
					quad += temp * temp;
				}

				if (kl < 0)	{
				//	cerr << "error : negative KL divergence\n";
				//	exit(1);
				}

				if (quad < 0)	{
					cerr << "error : negative quadratic distance\n";
					exit(1);
				}

				if (quad < ClusteringRange)	{
					
					if (affiliated == -1)	{
						for (int l=0; l<Nstate; l++)	{
							clusterTable[k][l] += sitenumber * stat[l];
						}
						ClusterSize[k] ++;
						ClusterWeight[k] += sitenumber;
						affiliated = k;
						IsAffiliated[k][i] = 1;
						for (int site=0; site<mParam->Nsite; site++)	{
							if (PB.Mode[site] == j)	{
								SiteAff[site][k]++;
							}
						}
					}
					else	{
						// merging clusters
					
						for (int l=0; l<Nstate; l++)	{
							clusterTable[affiliated][l] += clusterTable[k][l];
							clusterTable[k][l] = 0;
						}
						ClusterSize[affiliated] += ClusterSize[k];
						for (int ii=0; ii<i; ii++)	{	
							IsAffiliated[affiliated][ii] |= IsAffiliated[k][ii];
							IsAffiliated[k][ii] = 0;
						}
						IsAffiliated[affiliated][i] = 1;
						IsAffiliated[k][i] = 0;
						ClusterSize[k] = 0;
						ClusterWeight[affiliated] += ClusterWeight[k];
						ClusterWeight[k] = 0;
						for (int site=0; site<mParam->Nsite; site++)	{
							SiteAff[site][affiliated] += SiteAff[site][k];
							SiteAff[site][k] = 0;
						}
					}
				}

			}

			// look whether should stack down clusters

			// scan for first empty cluster : k;
			// scan for next non empty cluster : l;
			// while l < NCluster
			// else
			//	tranfer l to k
			// and increment l until a new non empty cluster is found or l == NCluster
			
			int k=0;
			while ( (k<NCluster) && (ClusterSize[k]))	{
				k++;
			}
			
			int l = k+1;
			while (l < NCluster)	{
			
				while ( (l<NCluster) && (! ClusterSize[l]))	{
					l++;
				}
				if (l<NCluster)	{
					// tranfer l to k;

					ClusterSize[k] = ClusterSize[l];
					for (int m=0; m<Nstate; m++)	{
						clusterTable[k][m] = clusterTable[l][m];
						clusterTable[l][m] = 0;
					}
					for (int ii=0; ii<=i; ii++)	{	
						IsAffiliated[k][ii] = IsAffiliated[l][ii];
						IsAffiliated[l][ii] = 0;
					}
					for (int site=0; site<mParam->Nsite; site++)	{
						SiteAff[site][k] = SiteAff[site][l];
						SiteAff[site][l] = 0;
					}
					ClusterSize[l] = 0;
					ClusterWeight[k] = ClusterWeight[l];
					ClusterWeight[l] = 0;
					k++;
					l++;
				}
			}

			NCluster = k;
			if (NCluster > NClusterMax)	{
				cerr << "error : NCluster overflow\n";
				exit(1);
			}
			
			// if not affiliated, create a new cluster

			if (affiliated == -1)	{
				ClusterSize[NCluster] = 1;
				for (int l=0; l<Nstate; l++)	{
					clusterTable[NCluster][l] = sitenumber * stat[l];
				}
				ClusterWeight[NCluster] = sitenumber;
				for (int site=0; site<mParam->Nsite; site++)	{
					if (PB.Mode[site] == j)	{
						SiteAff[site][NCluster]=1;
					}
				}
				NCluster++;
				if (NCluster >= NClusterMax)	{
					cerr << "error : NCluster overflow\n";
					exit(1);
				}

			}

		}
	}
	cerr << '\n';

	// ---------------------------------------
	// computing cluster related statistics
	// ---------------------------------------

	// compute cluster profiles
	
	for (int k=0; k<NCluster; k++)	{

		if (! ClusterSize[k])	{
			cerr << "error : void cluster\n";
			exit(1);
		}

		for (int m=0; m<Nstate; m++)	{
			clusterTable[k][m] /= ClusterWeight[k];
		}
	}

	int totalweight = 0;
	for (int l=0; l<NCluster; l++)	{
		totalweight += ClusterWeight[l];
	}

	// cluster stability
	for (int l=0; l<NCluster; l++)	{
		ClusterStability[l] = 0;
	}
	for (int k=0; k<NCluster; k++)	{
		for (int l=0; l<size; l++)	{
			ClusterStability[k] += IsAffiliated[k][l];
		}
	}
	// ---------
	// output
	// ---------
	
	// prepare file for output

	ofstream Cluster_os( (SampleName + ".clustermodes").c_str() );
	ofstream ClusterStability_os( (SampleName + ".clusterstability").c_str() );

	// write down cluster	
	// and make logos
	
	double totalWeight = 0;

	Cluster_os << "profile\tweight";
	for (int k=0; k<Nstate; k++)	{
		Cluster_os << '\t' << mParam->Alphabet[k];
	}
	Cluster_os << '\n' << '\n';	

	int* clusterlist = new int[NCluster];
	int FinalNCluster = 0;

	double garbageclusterstat[Nstate];
	double garbageclusterweight = 0;
	for (int k=0; k<Nstate; k++)	{
		garbageclusterstat[k] = 0;
	}

	for (int i=0; i<NCluster; i++)	{

		int max = 0;
		int cluster = 0;
		for (int j=0; j<NCluster; j++)	{
			if (max < ClusterWeight[j])	{
				max = ClusterWeight[j];
				cluster = j;
			}
		}


		clusterlist[i] = cluster;

		double weight = ((double) ClusterWeight[cluster]) / size;

		totalWeight += weight;
		ClusterWeight[cluster] =0;

		if (weight >= SizeThreshold)	{

			FinalNCluster++;

			double max = 0;
			double min = 1;

			for (int k=0; k<Nstate; k++)	{
				if (min > clusterTable[cluster][k])	{
					min = clusterTable[cluster][k];
				}
				if (max < clusterTable[cluster][k])	{
					max = clusterTable[cluster][k];
				}
			}

			for (int l=0; l<Nstate; l++)	{
				double importance = clusterTable[cluster][l]/max;
				if (importance > 0.6)	{
					Cluster_os << AminoAcids[l];
				}
				else if (importance > 0.30)	{
					Cluster_os << aminoacids[l];
				}
			}

			for (int l=0; l<Nstate; l++)	{
				double importance = clusterTable[cluster][l]/max;
				if (importance > 0.6)	{
					ClusterStability_os << AminoAcids[l];
				}
				else if (importance > 0.30)	{
					ClusterStability_os << aminoacids[l];
				}
			}

			ClusterStability_os << '\t' << ((int) ( ((double) (ClusterStability[cluster]) / size) * 100)) << '\n';

			Cluster_os << '\t' << weight ;
			for (int l=0; l<Nstate; l++)	{
				Cluster_os << '\t' << clusterTable[cluster][l];
			}
			Cluster_os << '\n';
		}

		else 	{

			for (int k=0; k<Nstate; k++)	{
				garbageclusterstat[k] += weight * clusterTable[cluster][k];
			}
			garbageclusterweight += weight;
		
		}
	}

	FinalNCluster++;
	for (int k=0; k<Nstate; k++)	{
		garbageclusterstat[k] /= garbageclusterweight;
	}

	Cluster_os << '\n';

	double max = 0;
	double min = 1;

	for (int k=0; k<Nstate; k++)	{
		if (min > garbageclusterstat[k])	{
			min = garbageclusterstat[k];
		}
		if (max < garbageclusterstat[k])	{
			max = garbageclusterstat[k];
		}
	}

	for (int l=0; l<Nstate; l++)	{
		double importance = garbageclusterstat[l]/max;
		if (importance > 0.6)	{
			Cluster_os << AminoAcids[l];
		}
		else if (importance > 0.30)	{
			Cluster_os << aminoacids[l];
		}
	}
	Cluster_os << '\t' << garbageclusterweight ;
	for (int l=0; l<Nstate; l++)	{
		Cluster_os << '\t' << garbageclusterstat[l];
	}
	Cluster_os << '\n';

	cout << "number of clusters : " << FinalNCluster << '\n';

	Cluster_os << '\n';
	Cluster_os << "total number of clusters (not including the garbage-cluster) : " << FinalNCluster - 1 << '\n';
	Cluster_os << "total weight of garbage cluster : " << ((int) (garbageclusterweight / totalWeight * 100)) << " \%" << '\n';
	// site affiliations
	ofstream SiteAff_os( (SampleName + ".siteaff").c_str() );
	SiteAff_os << "site\taffiliation frequencies\n";

	double check[FinalNCluster];
	for (int i=0; i<FinalNCluster; i++)	{
		check[i] = 0;
	}
	for (int i=0; i<mParam->Nsite; i++)	{
		SiteAff_os << i;
		double total = 0;
		for (int j=0; j<NCluster; j++)	{
			total += SiteAff[i][clusterlist[j]];
		}
		double wgarb = 0;
		for (int j=0; j<FinalNCluster-1; j++)	{
			double tmp = SiteAff[i][clusterlist[j]] / total;
			SiteAff_os << '\t' << tmp;
			check[j] += tmp;
			wgarb += tmp;
		}
		check[FinalNCluster-1] += 1 - wgarb;
		SiteAff_os << '\t' << 1 - wgarb << '\n';
	}

	return 1;
}

// ---------------------------------------------------------------------------
//		 GrepModes
// ---------------------------------------------------------------------------

void Sample::GrepModes()	{

	if (! GetSize())	{
		cerr << "error in Sample::GrepModes: sample of size 0\n";
		exit(1);
	}
	ofstream os((SampleName + ".profiles").c_str());
	PhyloBayes* pb = GetNextPB();
	os << pb->Nmode << '\t' << mParam->Nstate << '\n';
	cerr << pb->Nmode << " distinct profiles\n";
	for (int k=0; k<pb->Nmode; k++)	{
		int max = 0;
		int imax = 0;
		for (int l=0; l<pb->Nmode; l++)	{
			if (max < pb->SiteNumber[l])	{
				max = pb->SiteNumber[l];
				imax = l;
			}
		}
		os << ((double) pb->SiteNumber[imax]) / mParam->Nsite;
		for (int j=0; j<mParam->Nstate; j++)	{
			os << '\t' << pb->Stationary[imax][j];
		}
		os << '\n';
		pb->SiteNumber[imax] = 0;
	}
	cerr << "profiles of first point of the sample into " << SampleName << ".profiles\n";
	cerr << "\n";
}



// ---------------------------------------------------------------------------
//		 ComputeStatePostProb()
// ---------------------------------------------------------------------------


void Sample::ComputeStatePostProb(double*** StatePostProb, int singlefreq)	{

	/*
	int singlefreq = 0;
	int doublefreq = 0;
	int triplefreq = 0;
	*/
	int allinc = 0;

	// ofstream os((SampleName + ".ancstatepostprob").c_str());
	// double*** StatePostProb = new double**[mParam->Nnode];
	if (! StatePostProb)	{
		StatePostProb = new double**[mParam->Nnode];
		for (int i=0; i<mParam->Nnode; i++)	{
			StatePostProb[i] = new double*[mParam->Nsite];
			for (int j=0; j<mParam->Nsite; j++)	{
				StatePostProb[i][j] = new double[mParam->Nstate];
				for (int k=0; k<mParam->Nstate; k++)	{
					StatePostProb[i][j][k] = 0;
				}
			}
		}
	}
	else	{	
		for (int i=0; i<mParam->Nnode; i++)	{
			for (int j=0; j<mParam->Nsite; j++)	{
				for (int k=0; k<mParam->Nstate; k++)	{
					StatePostProb[i][j][k] = 0;
				}
			}
		}
	}
	for (int i=0; i<GetSize(); i++)	{
		cerr << '.';
		PhyloBayes* pb = GetNextPB();
		if (mParam->midpointrooting)	{
			pb->RootAtMidPoint();
		}
		pb->ComputePartialStatePostProb();
		for (int i=0; i<mParam->Nnode; i++)	{
			for (int j=0; j<mParam->Nsite; j++)	{
				for (int k=0; k<mParam->Nstate; k++)	{
					StatePostProb[i][j][k] += pb->StatePostProb[i][j][k];
				}
			}
		}
	}
	cerr << '\n';

	/*
	os << "tax1\ttax2";
	for (int j=0; j<mParam->Nsite; j++)	{
		os << '\t';
		for (int k=0; k<mParam->Nstate; k++)	{
			os << '\t' << mParam->Alphabet[k];
		}
	}
	os << '\n';
	*/
	for (int i=0; i<mParam->Nnode; i++)	{
		for (int j=0; j<mParam->Nsite; j++)	{
			for (int k=0; k<mParam->Nstate; k++)	{
				StatePostProb[i][j][k] /= GetSize();
			}
		}
	}

	for (int i=0; i<mParam->Nnode; i++)	{
		if (allinc ||  mParam->GetCurrentState()->nodefree(i))	{
			ostringstream s;
			s << SampleName << '_' << i << '_' << mParam->GetCurrentState()->GetLeftNodeName(i) << '_' << mParam->GetCurrentState()->GetRightNodeName(i);
			ofstream nos((s.str() + ".ancstatepostprob").c_str());
			nos << mParam->Nsite << '\t' << mParam->Nstate;
			for (int k=0; k<mParam->Nstate; k++)	{
				nos << '\t' << mParam->Alphabet[k];
			}
			nos << '\n';
			// os << mParam->GetCurrentState()->GetLeftNodeName(i) << '\t' << mParam->GetCurrentState()->GetRightNodeName(i);
			for (int j=0; j<mParam->Nsite; j++)	{
				// os << '\t';
				nos << j+1 << '\t';
				for (int k=0; k<mParam->Nstate; k++)	{
					// os << '\t' << StatePostProb[i][j][k];
					nos << '\t' << StatePostProb[i][j][k];
				}
				nos << '\n';
			}
			// os << '\n';
		}
	}

	Tree* tree = new Tree(mParam->GetCurrentState());
	// tree->Dichotomise();
	tree->ChangeNames(mParam->SpeciesNames, mParam->SpeciesNames, mParam->Ntaxa);
	tree->mRoot->SortLeavesAlphabetical();

	ofstream label_os((SampleName + ".labels").c_str());
	tree->Phylip(label_os, 0, 0, 1, 1);
	label_os.close();

	cerr << "ancestral sequence posterior probabilities in " << SampleName << "_<nodelabel>_<taxon1>_<taxon2>.ancstatepostprob\n";
	cerr << "node labels in " << SampleName << ".labels\n";

	int nentry = 0;
	if (! allinc)	{
		for (int i=0; i<mParam->Nnode; i++)	{
			if (mParam->GetCurrentState()->nodefree(i))	{
				nentry++;
			}
		}
	}
	else	{
		nentry = mParam->Nnode;
	}

	if (singlefreq)	{
		ofstream sos((SampleName + ".singlecomp").c_str());
		sos << nentry << '\t' << mParam->Nstate << '\n';
		for (int i=0; i<mParam->Nnode; i++)	{
			if (allinc || mParam->GetCurrentState()->nodefree(i))	{
				sos << mParam->GetCurrentState()->GetLeftNodeName(i) << '\t' << mParam->GetCurrentState()->GetRightNodeName(i);
				for (int k=0; k<mParam->Nstate; k++)	{
					double tot = 0;
					for (int j=0; j<mParam->Nsite; j++)	{
						tot += StatePostProb[i][j][k];
					}
					sos << '\t' << tot/mParam->Nsite;
				}
				sos << '\n';
			}
		}

		if (mParam->Nstate == 4)	{
			ofstream gcos((SampleName + ".gc").c_str());
			gcos << mParam->Nnode - 1 << '\t' << 1 << '\n';
			for (int i=0; i<mParam->Nnode; i++)	{
				// if (! mParam->GetCurrentState()->tree[i].isRoot())	{
				// if (blfreeonly && mParam->GetCurrentState()->blfree(i))	{
				if (allinc || mParam->GetCurrentState()->nodefree(i))	{
					gcos << mParam->GetCurrentState()->GetLeftNodeName(i) << '\t' << mParam->GetCurrentState()->GetRightNodeName(i);
					double tot = 0;
					for (int j=0; j<mParam->Nsite; j++)	{
						tot += StatePostProb[i][j][1] + StatePostProb[i][j][2];
					}
					tot /= mParam->Nsite;
					double x = log(tot / (1-tot));
					gcos << '\t' << x;
					gcos << '\n';
				}
			}
		}
	}

	/*
	if (doublefreq)	{
		ofstream dos((SampleName + ".doublecomp").c_str());
		dos << nentry << '\t' << mParam->Nstate * mParam->Nstate << '\n';
		for (int i=0; i<mParam->Nnode; i++)	{
			if (allinc || mParam->GetCurrentState()->nodefree(i))	{
				dos << mParam->GetCurrentState()->GetLeftNodeName(i) << '\t' << mParam->GetCurrentState()->GetRightNodeName(i);
				for (int k=0; k<mParam->Nstate; k++)	{
					for (int l=0; l<mParam->Nstate; l++)	{
						double tot = 0;
						for (int j=0; j<mParam->Nsite-1; j++)	{
							tot += StatePostProb[i][j][k] * StatePostProb[i][j+1][l];
						}
						dos << '\t' << tot/(mParam->Nsite - 1);
					}
				}
				dos << '\n';
			}
		}
	}

	if (triplefreq)	{
		ofstream tos((SampleName + ".triplecomp").c_str());
		tos << nentry << '\t' << mParam->Nstate * mParam->Nstate * mParam->Nstate << '\n';
		for (int i=0; i<mParam->Nnode; i++)	{
			if (allinc || mParam->GetCurrentState()->nodefree(i))	{
				tos << mParam->GetCurrentState()->GetLeftNodeName(i) << '\t' << mParam->GetCurrentState()->GetRightNodeName(i);
				for (int k=0; k<mParam->Nstate; k++)	{
					for (int l=0; l<mParam->Nstate; l++)	{
						for (int m=0; m<mParam->Nstate; m++)	{
							double tot = 0;
							for (int j=0; j<mParam->Nsite-2; j++)	{
								tot += StatePostProb[i][j][k] * StatePostProb[i][j+1][l] * StatePostProb[i][j+2][m];
							}
							tos << '\t' << tot/(mParam->Nsite - 2);
						}
					}
				}
				tos << '\n';
			}
		}
	}
	*/
}

// ---------------------------------------------------------------------------
//		 SampleSub()
// ---------------------------------------------------------------------------


void Sample::SampleSub(double*** StatePostProb, int* mask)	{

	int allinc = 0;

	if (! StatePostProb)	{
		StatePostProb = new double**[mParam->Nnode];
		for (int i=0; i<mParam->Nnode; i++)	{
			StatePostProb[i] = new double*[mParam->Nsite];
			for (int j=0; j<mParam->Nsite; j++)	{
				StatePostProb[i][j] = new double[mParam->Nstate];
				for (int k=0; k<mParam->Nstate; k++)	{
					StatePostProb[i][j][k] = 0;
				}
			}
		}
	}
	else	{	
		for (int i=0; i<mParam->Nnode; i++)	{
			for (int j=0; j<mParam->Nsite; j++)	{
				for (int k=0; k<mParam->Nstate; k++)	{
					StatePostProb[i][j][k] = 0;
				}
			}
		}
	}

	double** compo = new double*[mParam->Nnode];
	for (int i=0; i<mParam->Nnode; i++)	{
		compo[i] = new double[mParam->Nstate];
		for (int k=0; k<mParam->Nstate; k++)	{
			compo[i][k] = 0;
		}
	}

	double** instantcompo = new double*[mParam->Nnode];
	for (int i=0; i<mParam->Nnode; i++)	{
		instantcompo[i] = new double[mParam->Nstate];
		for (int k=0; k<mParam->Nstate; k++)	{
			instantcompo[i][k] = 0;
		}
	}

	int nfree = 0;
	for (int i=0; i<mParam->Nnode; i++)	{
		if (allinc ||  mParam->GetCurrentState()->nodefree(i))	{
			nfree++;
		}
	}

	for (int i=0; i<GetSize(); i++)	{
		cerr << '.';

		for (int j=0; j<mParam->Nnode; j++)	{
			for (int k=0; k<mParam->Nstate; k++)	{
				instantcompo[j][k] = 0;
			}
		}

		PhyloBayes* pb = GetNextPB();
		pb->Update();
		if (!i)	{
			pb->CreateTrueSub();
		}
		pb->DrawSiteVariables(0,0);
		
		pb->ResampleSub();
		for (int ii=0; ii<mParam->Nnode; ii++)	{
			for (int j=0; j<mParam->Nsite; j++)	{
				if (mask[j])	{
					if (mParam->ModeFastCompute)	{
						StatePostProb[ii][j][pb->TrueState[ii][j]] ++ ;
						compo[ii][pb->TrueState[ii][j]] ++ ;
						instantcompo[ii][pb->TrueState[ii][j]] ++ ;
					}
					else	{
						StatePostProb[ii][j][pb->State[ii][j]] ++ ;
						compo[ii][pb->State[ii][j]] ++ ;
						instantcompo[ii][pb->State[ii][j]] ++ ;
					}
				}
			}
		}

		for (int j=0; j<mParam->Nnode; j++)	{
			double tot = 0;
			for (int k=0; k<mParam->Nstate; k++)	{
				tot += instantcompo[j][k];
			}
			for (int k=0; k<mParam->Nstate; k++)	{
				instantcompo[j][k] /= tot;
			}
		}
		
		ostringstream s;
		s << SampleName << "_" << i;
		ofstream compos((s.str() + ".comp").c_str());
		compos << nfree << '\t' << mParam->Nstate << '\n';

		for (int j=0; j<mParam->Nnode; j++)	{
			if (allinc ||  mParam->GetCurrentState()->nodefree(j))	{
				compos << mParam->GetCurrentState()->GetLeftNodeName(j) << '\t' << mParam->GetCurrentState()->GetRightNodeName(j);
				for (int k=0; k<mParam->Nstate; k++)	{
					compos << '\t' << instantcompo[j][k];
				} 
				compos << '\n';
			}
		}
	}
	cerr << '\n';

	for (int i=0; i<mParam->Nnode; i++)	{
		for (int j=0; j<mParam->Nsite; j++)	{
			for (int k=0; k<mParam->Nstate; k++)	{
				StatePostProb[i][j][k] /= GetSize();
			}
		}
		double tot = 0;
		for (int k=0; k<mParam->Nstate; k++)	{
			tot += compo[i][k];
		}
		for (int k=0; k<mParam->Nstate; k++)	{
			compo[i][k] /= tot;
		}

	}

	/*
	os << "tax1\ttax2";
	for (int j=0; j<mParam->Nsite; j++)	{
		os << '\t';
		for (int k=0; k<mParam->Nstate; k++)	{
			os << '\t' << mParam->Alphabet[k];
		}
	}
	os << '\n';
	*/

	ofstream compos((SampleName + ".comp").c_str());
	compos << nfree << '\t' << mParam->Nstate << '\n';

	ofstream hos((SampleName + ".entropy").c_str());
	hos << nfree << '\t' << mParam->Nstate << '\n';


	for (int i=0; i<mParam->Nnode; i++)	{
		if (allinc ||  mParam->GetCurrentState()->nodefree(i))	{
			ostringstream s;
			s << SampleName << '_' << i << '_' << mParam->GetCurrentState()->GetLeftNodeName(i) << '_' << mParam->GetCurrentState()->GetRightNodeName(i);
			// s << SampleName << '_' << mParam->GetCurrentState()->GetLeftNodeName(i) << '_' << mParam->GetCurrentState()->GetRightNodeName(i);
			ofstream nos((s.str() + ".ancstatepostprob").c_str());
			nos << mParam->Nsite << '\t' << mParam->Nstate;
			for (int k=0; k<mParam->Nstate; k++)	{
				nos << '\t' << mParam->Alphabet[k];
			}
			nos << '\n';
			// os << mParam->GetCurrentState()->GetLeftNodeName(i) << '\t' << mParam->GetCurrentState()->GetRightNodeName(i);
			double meanh = 0;
			for (int j=0; j<mParam->Nsite; j++)	{
				// os << '\t';
				nos << j+1 << '\t';
				double h = 0;
				for (int k=0; k<mParam->Nstate; k++)	{
					// os << '\t' << StatePostProb[i][j][k];
					nos << '\t' << StatePostProb[i][j][k];
					if (StatePostProb[i][j][k] > 1e-8)	{
						h -= StatePostProb[i][j][k] * log(StatePostProb[i][j][k]);
					}
				}
				nos << '\n';
				meanh += h;
			}
			// os << '\n';
			meanh /= mParam->Nsite;
			hos << mParam->GetCurrentState()->GetLeftNodeName(i) << '\t' << mParam->GetCurrentState()->GetRightNodeName(i);
			hos << '\t' << meanh << '\t' << exp(meanh) << '\n';

			// global compositions
			compos << mParam->GetCurrentState()->GetLeftNodeName(i) << '\t' << mParam->GetCurrentState()->GetRightNodeName(i);
			for (int k=0; k<mParam->Nstate; k++)	{
				compos << '\t' << compo[i][k];
			} 
			compos << '\n';
		}
	}

	Tree* tree = new Tree(mParam->GetCurrentState());
	// tree->Dichotomise();
	tree->ChangeNames(mParam->SpeciesNames, mParam->SpeciesNames, mParam->Ntaxa);
	tree->mRoot->SortLeavesAlphabetical();

	ofstream label_os((SampleName + ".labels").c_str());
	tree->Phylip(label_os, 0, 0, 1, 1);
	label_os.close();

	cerr << "ancestral sequence posterior probabilities in " << SampleName << "_<nodelabel>_<taxon1>_<taxon2>.ancstatepostprob\n";
	cerr << "node labels in " << SampleName << ".labels\n";


}

void Sample::ConstantSiteCorrection(double rho)	{

	ofstream os((SampleName + ".cst").c_str());

	double* logweight = new double[GetSize()];
	double* weight = new double[GetSize()];

	double meanpv = 0;

	double max = 0;
	for (int i=0; i<GetSize(); i++)	{
		cerr << '.';
		PhyloBayes* pb = GetNextPB();
		pb->Update();

		double tmp = 0;
		double correct = pb->ConstantSiteCorrection(rho, tmp);
		meanpv += tmp;
		logweight[i] = correct;
		if ((!i) || (max < correct))	{
			max = correct;
		}
	}
	double sum = 0;
	double sum2 = 0;
	for (int i=0; i<GetSize(); i++)	{
		weight[i] = exp(logweight[i] - max);
		sum += weight[i];
		sum2 += weight[i] * weight[i];
		os << weight[i] << '\n';
	}

	cerr << "effective sample size : " << sum * sum / sum2  << '\t' << GetSize() << '\n';
	meanpv /= GetSize();
	cerr << "mean prob of being variable across sites: " << meanpv << '\n';
}

void Sample::SampleSub(int nrep, int ppred, int priormode, int priorratemode, int priorrootstate, int* sitemask)	{

	ofstream os((SampleName + ".sub").c_str());
	for (int i=0; i<GetSize(); i++)	{
		cerr << '.';
		PhyloBayes* pb = GetNextPB();
		pb->Update();
		if (!i)	{
			pb->CreateTrueSub();
		}
		pb->DrawSiteVariables(priormode, priorratemode);
		
		os << "point " << i+1 << '\n';
		for (int rep=0; rep<nrep; rep++)	{
			os << "rep " << rep+1 << '\n';
			if (ppred)	{
				os << "posterior\n";
			}
			pb->ResampleSub(&os, sitemask);
			if (ppred)	{
				os << "posterior predictive\n";
				pb->SimulateData(priormode, priorratemode, priorrootstate,1,&os,sitemask);
				pb->RestoreData();
				os << '\n';
			}
		}
		os << '\n';
	}
	cerr << '\n';
	
	cerr << "substitution mappings in " << SampleName << ".sub\n";
}
/*
void Sample::SampleSubWS(int priormode, int priorratemode, int priorrootstate)	{

	ofstream os((SampleName + ".sub").c_str());
	PhyloBayes* pb = GetNextPB();
	pb->Update();
	// pb->CreateTrueSub();
	pb->DrawSiteVariables(priormode, priorratemode);
	pb->ResampleSub();
	for (int i=0; i<mParam->Nnode; i++)	{
		Node* node = &(mParam->GetCurrentState()->tree[i]);
		if ((node != mParam->GetCurrentState()->root) && (node != mParam->GetCurrentState()->root->left) && (node != mParam->GetCurrentState()->root->right))	{
			int totw = 0;
			int tots = 0;
			for (int j=0; j<mParam->Nsite; j++)	{
				totw += pb->NHNsub[i][j][0];
				tots += pb->NHNsub[i][j][1];
			}
			os << mParam->GetCurrentState()->GetLeftNodeName(i) << '\t' << mParam->GetCurrentState()->GetRightNodeName(i) << '\t' << totw << '\t' << tots << '\n';
		}
	}
	cerr << "substitution mappings in " << SampleName << ".sub\n";
}
*/

void Sample::SampleSubWS(int priormode, int priorratemode, int priorrootstate)	{

	ofstream os((SampleName + ".sub").c_str());
	PhyloBayes* pb = GetNextPB();
	pb->Update();
	// pb->CreateTrueSub();
	pb->DrawSiteVariables(priormode, priorratemode);
	pb->ResampleSub();
	for (int i=0; i<mParam->Nnode; i++)	{
		Node* node = &(mParam->GetCurrentState()->tree[i]);
		if ((node != mParam->GetCurrentState()->root) && (node != mParam->GetCurrentState()->root->left) && (node != mParam->GetCurrentState()->root->right))	{
			int Nstate = mParam->Nstate;
			int tot[Nstate];
			for (int k=0; k<Nstate; k++)	{
				tot[k] = 0;
			}
			for (int j=0; j<mParam->Nsite; j++)	{
				for (int k=0; k<Nstate; k++)	{
					tot[k] += pb->NHNsub[i][j][k];
				}
			}
			os << mParam->GetCurrentState()->GetLeftNodeName(i) << '\t' << mParam->GetCurrentState()->GetRightNodeName(i);
			for (int k=0; k<Nstate; k++)	{
				os << '\t' << tot[k];
			}
			os << '\n';
		}
	}
	cerr << "substitution mappings in " << SampleName << ".sub\n";
}

// ---------------------------------------------------------------------------
//		 Coaff()
// ---------------------------------------------------------------------------

void Sample::Coaff(int type, int step)	{

	int Nsite = mParam->Nsite;
	int N = Nsite / step;
	double coaff[N][N];
	for (int i=0; i<N; i++)	{
		for (int j=i+1; j<N; j++)	{
			coaff[i][j] = 0;
		}
	}
	for (int i=0; i<GetSize(); i++)	{
		PhyloBayes* pb = GetNextPB();
		int* mode = 0;
		if (type == 1)	{
			mode = pb->Mode;
		}
		else if (type == 2)	{
			mode = pb->RateMode;
		}
		else if (type == 3)	{
			mode = pb->MixBLMode;
		}
		for (int i=0; i<N; i++)	{
			for (int j=i+1; j<N; j++)	{
				if (mode[step*i] == mode[step*j])	{
					coaff[i][j] ++;
				}
			}
		}
	}

	ofstream os((SampleName + ".coaff").c_str());
	for (int i=0; i<N; i++)	{
		for (int j=i+1; j<N; j++)	{
			coaff[i][j] /= GetSize();
			os << i << '\t' << j << '\t' << coaff[i][j] << '\n';
		}
	}
}

// ---------------------------------------------------------------------------
//		 ReadContCharacters()
// ---------------------------------------------------------------------------

void Sample::ReadContCharacters()	{

	int Nnode = mParam->Nnode;
	int Nchar = mParam->ContNchar;
	cerr << "number of continuous characters: " << Nchar << '\n';
	cerr.flush();

	double mean[Nnode][Nchar];
	double var[Nnode][Nchar];
	for (int j=0; j<Nnode; j++)	{
		for (int i=0; i<Nchar; i++)	{
			mean[j][i] = 0;
			var[j][i] = 0;
		}
	}
	for (int i=0; i<GetSize(); i++)	{
		PhyloBayes* pb = GetNextPB();
		for (int j=0; j<Nnode; j++)	{
			for (int i=0; i<Nchar; i++)	{
				int k = pb->NHcat[j];
				double tmp = pb->logNHStatDistorter[k][mParam->Nstate + i];
				mean[j][i] += tmp;
				var[j][i] += tmp * tmp;
			}
		}
	}
	ofstream os((SampleName + ".contchar").c_str());
	for (int j=0; j<Nnode; j++)	{
		os << j << '\t';
		if (j<mParam->Ntaxa)	{
			os << mParam->SpeciesNames[j] << '\t';
		}
		else	{
			os << "internal\t";
		}
		for (int i=0; i<Nchar; i++)	{
			mean[j][i] /= GetSize();
			var[j][i] /= GetSize();
			var[j][i] -= mean[j][i] * mean[j][i];
			os << mean[j][i] << '\t';
			// os << mean[j][i] << " +/- " << sqrt(var[j][i]) << '\t';
		}
		os << '\n';
	}
}

// ---------------------------------------------------------------------------
//		 TrueLength(int nrep)
// ---------------------------------------------------------------------------


void Sample::TrueLength(int nrep, int sitebranch)	{
		
	double meanbl[mParam->Nnode];
	double varbl[mParam->Nnode];

	double meanppbl[mParam->Nnode];
	double varppbl[mParam->Nnode];

	double meanlength = 0;
	double varlength = 0;

	double meanpplength = 0;
	double varpplength = 0;

	double meanl = 0;
	double varl = 0;

	double meanpbl[mParam->Nnode];
	double varpbl[mParam->Nnode];

	for (int j=0; j<mParam->Nnode; j++)	{
		meanbl[j] = 0;
		varbl[j] = 0;
		meanpbl[j] = 0;
		varpbl[j] = 0;
		meanppbl[j] = 0;
		varppbl[j] = 0;
	}

	double meansite[mParam->Nsite];
	double varsite[mParam->Nsite];
	for (int i=0; i<mParam->Nsite; i++)	{
		meansite[i] = 0;
		varsite[i] = 0;
	}

	double** meansitebranch = 0;
	double** varsitebranch = 0;
	if (sitebranch)	{
		meansitebranch = new double*[mParam->Nsite];
		varsitebranch = new double*[mParam->Nsite];
		for (int i=0; i<mParam->Nsite; i++)	{
			meansitebranch[i] = new double[mParam->Nnode];
			varsitebranch[i] = new double[mParam->Nnode];
			for (int j=0; j<mParam->Nnode; j++)	{
				meansitebranch[i][j] = 0;
				varsitebranch[i][j] = 0;
			}
		}
	}

	ofstream tos((SampleName + ".alltrees.efflength").c_str());

	for (int i=0; i<GetSize(); i++)	{
		cerr << ".";
		cerr.flush();

		PhyloBayes* pb = GetNextPB();
		pb->Update();

		if (!i)	{
			pb->CreateTrueSub();
		}
		if (mParam->ActivateClock && (mParam->BLModelSwitch == 0))	{
			pb->SwapBL();
		}
		double norm = pb->LengthNormFactor();
		double l = pb->GetLength() * norm;
		meanl += l;
		varl += l * l;

		for (int j=0; j<mParam->Nnode; j++)	{
			double tmp = pb->BL[j] * norm;
			meanpbl[j] += tmp;
			varpbl[j] += tmp * tmp;
		}

		if (mParam->ActivateClock && (mParam->BLModelSwitch == 0))	{
			pb->SwapBL();
		}

		if (mParam->SumOverRateModes)	{
			pb->UpdateRateModeSiteLogSampling();
			pb->SumRateModeSiteLogSampling();
		}

		for (int rep=0; rep<nrep; rep++)	{
			pb->ResampleSub();
			pb->UpdateTotalTrueSub();
			if (sitebranch)	{
				for (int i=0; i<mParam->Nsite; i++)	{
					for (int j=0; j<mParam->Nnode; j++)	{
						meansitebranch[i][j] += pb->BranchSiteTotalTrueSub[j][i];
						varsitebranch[i][j] += pb->BranchSiteTotalTrueSub[j][i] * pb->BranchSiteTotalTrueSub[j][i];
					}
				}
			}
			for (int j=0; j<mParam->Nnode; j++)	{
				pb->BL[j] = ((double) pb->BranchTotalTrueSub[j]) / mParam->Nsite;
			}
			pb->ReverseSetBL();
			pb->Phylip(tos,1,0,1,0);

			for (int i=0; i<mParam->Nsite; i++)	{
				meansite[i] += pb->TotalTrueSub[i];
				varsite[i] += pb->TotalTrueSub[i] * pb->TotalTrueSub[i];
			}

			double total = 0;
			for (int j=0; j<mParam->Nnode; j++)	{
				double temp = ((double) pb->BranchTotalTrueSub[j]) / mParam->Nsite;
				meanbl[j] += temp;
				varbl[j] += temp * temp;
				total += temp;
			}
			meanlength += total;
			varlength += total * total;
		}
	}
	cerr << '\n';
	tos.close();
	int size = GetSize() * nrep;

	meanl /= GetSize();
	varl /= GetSize();
	varl -= meanl * meanl;

	meanlength /= size;
	varlength /= size;
	varlength -= meanlength * meanlength;

	meanpplength /= size;
	varpplength /= size;
	varpplength -= meanpplength * meanpplength;

	for (int j=0; j<mParam->Nnode; j++)	{
		meanbl[j] /= size;
		varbl[j] /= size;
		varbl[j] -= meanbl[j] * meanbl[j];
		meanpbl[j] /= size;
		varpbl[j] /= size;
		varpbl[j] -= meanpbl[j] * meanpbl[j];
		meanppbl[j] /= size;
		varppbl[j] /= size;
		varppbl[j] -= meanppbl[j] * meanppbl[j];
	}

	PhyloBayes* pb = mParam->GetCurrentState();
	// PhyloBayes* pb = mParam->GetNextState();
	if (sitebranch)	{
		ofstream os((SampleName + ".sitebranch").c_str());
		for (int i=0; i<mParam->Nsite; i++)	{
			os << i << '\t';
			for (int j=0; j<mParam->Nnode; j++)	{
				meansitebranch[i][j] /= size;
				varsitebranch[i][j] /= size;
				varsitebranch[i][j] -= meansitebranch[i][j] * meansitebranch[i][j];
				// os << '\t' << meansitebranch[i][j];
				pb->BL[j] = meansitebranch[i][j];
			}
			pb->ReverseSetBL();
			pb->Phylip(os,1,0,1,0);
		}
		os.close();
		for (int i=0; i<mParam->Nsite; i++)	{
			delete[] meansitebranch[i];
			delete[] varsitebranch[i];
		}
		delete[] meansitebranch;
		delete[] varsitebranch;
	}
	for (int i=0; i<mParam->Nsite; i++)	{
		meansite[i] /= size;
		varsite[i] /= size;
		varsite[i] -= meansite[i] * meansite[i];
	}
	ofstream ros((SampleName + ".effabsrate").c_str());
	for (int i=0; i<mParam->Nsite; i++)	{
		ros << i << '\t' << meansite[i] << '\t' << sqrt(varsite[i]) << '\n';
	}

	ofstream pos((SampleName + ".param.tre").c_str());
	ofstream subos((SampleName + ".subs.tre").c_str());
	ofstream os((SampleName + ".efflength").c_str());

	for (int j=0; j<mParam->Nnode; j++)	{
		pb->BL[j] = meanpbl[j];
	}
	pb->ReverseSetBL();
	pb->Phylip(pos,1,0,1,0);
	
	for (int j=0; j<mParam->Nnode; j++)	{
		pb->BL[j] = meanbl[j];
	}
	pb->ReverseSetBL();
	pb->Phylip(subos,1,0,1,0);

	/*
	for (int j=0; j<mParam->Nnode; j++)	{
		pb->BL[j] = meanppbl[j];
	}
	pb->ReverseSetBL();
	pb->Phylip(ppsubos,1,0,1,0);
	*/

	os << "#\tparam\tsubs\tstderr\n";
	os << '\n';
	for (int j=0; j<mParam->Nnode; j++)	{
		os << j << '\t' << meanpbl[j] << '\t' << meanbl[j] << '\t' << sqrt(varppbl[j]) << '\n';
	}
	os << '\n';
	os << "parametric length   : " << meanl << " +/- " << sqrt(varl) << '\n';
	// os << "ppred subst  length : " << meanpplength << " +/- " << sqrt(varpplength) << '\n';
	os << "substitution length : " << meanlength << " +/- " << sqrt(varlength) << '\n';
	os << '\n';

	cout << '\n';
	cout << "parametric length   : " << meanl << " +/- " << sqrt(varl) << '\n';
	// cout << "ppred subst  length : " << meanpplength << " +/- " << sqrt(varpplength) << '\n';
	cout << "substitution length : " << meanlength << " +/- " << sqrt(varlength) << '\n';
	cout << '\n';

	for (int j=0; j<mParam->Nnode; j++)	{
		pb->tree[j].branchLength = j;
		pb->BL[j] = j;
	}
	ofstream label_os((SampleName + ".labels").c_str());
	pb->Phylip(label_os, 1, 0, 1, 0);
	label_os.close();

	cerr << "\n";
	cerr << "tree with parametric branch lengths in                 : " << SampleName << ".params.tre\n";
	cerr << "tree with substitutional branch lengths in             : " << SampleName << ".subs.tre\n";
	cerr << "labelled tree in                                       : " << SampleName << ".labels\n";
	cerr << "\n";
	cerr << "branch-tabulated lengths in                            : " << SampleName << ".efflength\n";
	cerr << "site-tabulated mean total number of substitutions in   : " << SampleName << ".effabsrate\n";
	cerr << "\n";

}


// ---------------------------------------------------------------------------
//		 SiteCount()
// ---------------------------------------------------------------------------


void Sample::SiteCount(int nrep)	{

	for (int i=0; i<GetSize(); i++)	{

		cerr << ".";
		cerr.flush();

		PhyloBayes* pb = GetNextPB();
		pb->Update();
		
		for (int rep=0; rep<nrep; rep++)	{
			ostringstream s;
			if (nrep > 1)	{
				s << SampleName << "_" << i << "_" << rep << ".sitecount";
			}
			else	{
				s << SampleName << "_" << i << ".ali";
			}
			string name = s.str();

			pb->ResampleSub();
			pb->WriteSiteCount(name);
		}
	}

	cerr << '\n';
	cerr << "site counts in " << SampleName << "_n.sitecount\n"; 
	cerr << '\n';
}

// ---------------------------------------------------------------------------
//		 PosteriorPredictive()
// ---------------------------------------------------------------------------


void Sample::PosteriorPredictive(int nrep, int priormode, int priorratemode, int priorrootstate)	{

	for (int i=0; i<GetSize(); i++)	{

		cerr << ".";
		cerr.flush();

		PhyloBayes* pb = GetNextPB();
		pb->Update();
		
		for (int rep=0; rep<nrep; rep++)	{
			ostringstream s;
			if (nrep > 1)	{
				s << SampleName << "_" << i << "_" << rep << ".ali";
			}
			else	{
				s << SampleName << "_" << i << ".ali";
			}
			string name = s.str();

			pb->SimulateData(priormode, priorratemode, priorrootstate);
			mParam->WriteDataToFile(name);
			pb->RestoreData();
		}
	}

	cerr << '\n';
	cerr << "replicates in " << SampleName << "_n.ali\n"; 
	cerr << '\n';
}


// ---------------------------------------------------------------------------
//		ClockRate()
// ---------------------------------------------------------------------------


void Sample::ClockRateAutoCorrel(int nrep)	{

	double meanobs = 0;
	double varobs = 0;
	double meanpred = 0;
	double varpred = 0;
	double pp = 0;
	int size = GetSize();
	for (int i=0; i<GetSize(); i++)	{

		cerr << ".";
		cerr.flush();

		PhyloBayes* pb = GetNextPB();
		
		for (int rep=0; rep<nrep; rep++)	{
			pb->SetTimes();
			pb->SetTau();
			double obs = pb->GetAutoCorrelation();
			meanobs += obs;
			varobs += obs * obs;
			pb->ResampleClockRate(pb->root->label);
			if (mParam->BLModelSwitch == 1)	{
				pb->ResampleClockLength(pb->root->label);
			}
			pb->SetTimes();
			pb->SetTau();
			double pred = pb->GetAutoCorrelation();
			meanpred += pred;
			varpred += pred * pred;
			if (obs < pred)	{
				pp++;
			}
		}
	}
	cerr << '\n';
	cerr << '\n';

	meanobs /= size;
	varobs /= size;
	varobs -= meanobs * meanobs;
	meanpred /= size;
	varpred /= size;
	varpred -= meanpred * meanpred;

	ofstream os((SampleName + "autocorrel").c_str());
	os << "observed   autocorrelation : " << meanobs << "\t+/-\t" << sqrt(varobs) << '\n';
	os << "predictive autocorrelation : " << meanpred << "\t+/-\t" << sqrt(varpred) << '\n';
	os << "p-value : " << pp / size << '\n';

}

void Sample::ClockRate(int nrep)	{

	ofstream postos((SampleName + ".postlength").c_str());
	ofstream predos((SampleName + ".predlength").c_str());
	double meanobsl = 0;
	double varobsl = 0;
	double meanpredl = 0;
	double varpredl = 0;
	double pp = 0;
	int size = GetSize();
	for (int i=0; i<GetSize(); i++)	{

		cerr << ".";
		cerr.flush();

		PhyloBayes* pb = GetNextPB();
		pb->Update();
		
		for (int rep=0; rep<nrep; rep++)	{
			pb->SetTimes();
			pb->SetTau();
			if (mParam->BLModelSwitch == 0)	{
				pb->SwapBL();
			}
			double obsl = pb->GetLength();
			meanobsl += obsl;
			varobsl += obsl * obsl;
			Tree tree(pb);
			tree.Dichotomise();
			tree.ChangeNames(mParam->SpeciesNames, mParam->SpeciesNames, mParam->Ntaxa);
			tree.mRoot->SortLeavesAlphabetical();
			tree.Phylip(postos, 1, 0, 1, 0);
			if (mParam->BLModelSwitch == 0)	{
				pb->SwapBL();
			}
			if (mParam->BLModelSwitch == 0)	{
				pb->ResampleClockRate(pb->root->label);
				pb->SetTau();
				pb->SwapBL();
			}
			else	{
				pb->ResampleClockRate(pb->root->label);
				pb->ResampleClockLength(pb->root->label);
			}
			double predl = pb->GetLength();
			meanpredl += predl;
			varpredl += predl * predl;
			if (obsl > predl)	{
				pp++;
			}
			Tree tree2(pb);
			tree2.Dichotomise();
			tree2.ChangeNames(mParam->SpeciesNames, mParam->SpeciesNames, mParam->Ntaxa);
			tree2.mRoot->SortLeavesAlphabetical();
			tree2.Phylip(predos, 1, 0, 1, 0);
			if (mParam->BLModelSwitch == 0)	{
				pb->SwapBL();
			}
		}
	}
	cerr << '\n';
	cerr << '\n';

	meanobsl /= size;
	varobsl /= size;
	varobsl -= meanobsl * meanobsl;
	meanpredl /= size;
	varpredl /= size;
	varpredl -= meanpredl * meanpredl;

	cout << "observed   total length  : " << meanobsl << "\t+/-\t" << sqrt(varobsl) << '\n';
	cout << "predictive total length  : " << meanpredl << "\t+/-\t" << sqrt(varpredl) << '\n';
	cout << "p-value : " << pp / size << '\n';

	cerr << '\n';
	cerr << "posterior  lengths in " << SampleName << ".postlength\n"; 
	cerr << "predictive lengths in " << SampleName << ".predlength\n"; 
	cerr << '\n';
}

// ---------------------------------------------------------------------------
//		 Posterior predictive LRT 
// ---------------------------------------------------------------------------

void Sample::postpredDoubleCount(int nrep, int ncat, int priormode, int priorratemode, int priorrootstate)	{

	int size = GetSize();
	double meanobsdisp = 0;
	double varobsdisp = 0;
	double meandisp = 0;
	double vardisp = 0;
	double pvdisp = 0;

	int Nstate = mParam->Nstate;
	int Nrr = Nstate * (Nstate-1) / 2;
	double** mat = new double*[mParam->Nstate];
	for (int k=0; k<Nstate; k++)	{
		mat[k] = new double[mParam->Nstate];
	}

	double meanobsrr[Nrr];
	double meanpredrr[Nrr];
	double varobsrr[Nrr];
	double varpredrr[Nrr];
	double* rr = new double[Nrr];
	for (int m=0; m<Nrr; m++)	{
		meanobsrr[m] = 0;
		meanpredrr[m] = 0;
	}

	for (int i=0; i<size; i++)	{
		cerr << '.';
		cerr.flush();
		PhyloBayes* pb = GetNextPB();
		pb->Update();
		pb->countmultiplesub = 1;
		pb->CreateMultipleSub();
		
		for (int rep=0; rep<nrep; rep++)	{
			pb->DrawSiteVariables(0,0);
			pb->ResampleSub();

			double obs = pb->CountMatrix(mat,rr);
			for (int m=0; m<Nrr; m++)	{
				meanobsrr[m] += rr[m];
				varobsrr[m] += rr[m] * rr[m];
			}
			
			pb->SimulateData(priormode, priorratemode, priorrootstate,0);

			double pred = pb->CountMatrix(mat,rr);
			for (int m=0; m<Nrr; m++)	{
				meanpredrr[m] += rr[m];
				varpredrr[m] += rr[m] * rr[m];
			}

			meanobsdisp += obs;
			varobsdisp += obs*obs;
			meandisp += pred;
			vardisp += pred*pred;
			if (pred > obs)	{
				pvdisp++;
			}
			pb->RestoreData();
		}
		pb->DeleteMultipleSub();
	}
	cerr << '\n';
	cerr.flush();

	meandisp /= size;
	vardisp /= size;
	vardisp -= meandisp * meandisp;
	meanobsdisp /= size;
	varobsdisp /= size;
	varobsdisp -= meanobsdisp * meanobsdisp;

	pvdisp /= size;

	ofstream sos((SampleName + ".doubledisp").c_str());

	sos << "observed double disp   : " << meanobsdisp << " +/- " << sqrt(varobsdisp) << '\n';
	sos << "predictive             : " << meandisp << " +/- " << sqrt(vardisp) << '\n';
	sos << "pp value               : " << pvdisp << '\t' << "(" << size * nrep << ")" << '\n';
	sos.close();

	cout << '\n';
	system (("cat " + SampleName + ".doubledisp").c_str());

	ofstream oos((SampleName + ".obsrr").c_str());
	ofstream pos((SampleName + ".predrr").c_str());
	ofstream opos((SampleName + ".obspredrr").c_str());
	for (int m=0; m<Nrr; m++)	{
		meanobsrr[m] /= size;
		varobsrr[m] /= size;
		varobsrr[m] -= meanobsrr[m] * meanobsrr[m];
		meanpredrr[m] /= size;
		varpredrr[m] /= size;
		varpredrr[m] -= meanpredrr[m] * meanpredrr[m];
		oos << meanobsrr[m]/Nrr << '\t' << sqrt(varobsrr[m])/Nrr << '\n';
		pos << meanpredrr[m]/Nrr << '\t' << sqrt(varpredrr[m])/Nrr << '\n';
		opos << meanobsrr[m]/Nrr << '\t' << sqrt(varobsrr[m])/Nrr << '\t' << meanpredrr[m]/Nrr << '\t' << sqrt(varpredrr[m])/Nrr << '\n';
	}

	for (int k=0; k<Nstate; k++)	{
		delete[] mat[k];
	}
	delete[] mat;
	delete[] rr;
}

	

void Sample::postpredLRT(int nrep, int ncat, int priormode, int priorratemode, int priorrootstate)	{

	int size = GetSize();

	double* obsrate = new double[nrep * size];
	double meanobsrate = 0;
	double varobsrate = 0;
	double* predrate = new double[nrep * size];
	double meanpredrate = 0;
	double varpredrate = 0;
	double pvrate = 0;
	double meanratediff = 0;
	double varratediff = 0;

	double* obsprofile = new double[nrep * size];
	double meanobsprofile = 0;
	double varobsprofile = 0;
	double* predprofile = new double[nrep * size];
	double meanpredprofile = 0;
	double varpredprofile = 0;
	double pvprofile = 0;
	double meanprofilediff = 0;
	double varprofilediff = 0;

	double* obsdouble = new double[nrep * size];
	double meanobsdouble = 0;
	double varobsdouble = 0;
	double* preddouble = new double[nrep * size];
	double meanpreddouble = 0;
	double varpreddouble = 0;
	double pvdouble = 0;
	double meandoublediff = 0;
	double vardoublediff = 0;

	double obsdiversity = mParam->MeanDiversity();
	double ppmeandiversity = 0;
	double ppvardiversity = 0;
	double pvmeandiversity = 0;

	double obsinv = mParam->PropInvariant();
	double meaninv = 0;
	double varinv = 0;
	double pvinv = 0;

	for (int i=0; i<size; i++)	{
		cerr << '.';
		cerr.flush();
		PhyloBayes* pb = GetNextPB();
		pb->Update();
		if (!i)	{
			pb->CreateTrueSub();
		}
		pb->countmultiplesub = 1;
		pb->CreateMultipleSub();
		
		for (int rep=0; rep<nrep; rep++)	{
			double mean = 0;
			double max = 0;
			pb->DrawSiteVariables(0,0);
			pb->ResampleSub();

			pb->RateDeltaLog(mean,max);
			obsrate[i*nrep + rep] = mean;
			meanobsrate += mean;
			varobsrate += mean * mean;
			
			pb->ProfileDeltaLog(mean,max);
			// pb->OptimiseSiteProfile(mean,max);
			obsprofile[i*nrep + rep] = mean;
			meanobsprofile += mean;
			varobsprofile += mean * mean;
			
			// mean = pb->DoubleBias();
			pb->DoubleDeltaLog(mean,max);
			obsdouble[i*nrep + rep] = mean;
			meanobsdouble += mean;
			varobsdouble += mean * mean;
			
			pb->SimulateData(priormode, priorratemode, priorrootstate,0);

			double diversity = mParam->MeanDiversity();
			ppmeandiversity += diversity;
			ppvardiversity += diversity*diversity;
			if (diversity < obsdiversity)	{
				pvmeandiversity ++;
			}

			double inv = mParam->PropInvariant();
			meaninv += inv;
			varinv += inv*inv;
			if (inv > obsinv)	{
				pvinv ++;
			}

			pb->RateDeltaLog(mean,max);
			predrate[i*nrep + rep] = mean;
			meanpredrate += mean;
			varpredrate += mean * mean;
			
			pb->ProfileDeltaLog(mean,max);
			predprofile[i*nrep + rep] = mean;
			meanpredprofile += mean;
			varpredprofile += mean * mean;
			
			// mean = pb->DoubleBias();
			pb->DoubleDeltaLog(mean,max);
			preddouble[i*nrep + rep] = mean;
			meanpreddouble += mean;
			varpreddouble += mean * mean;

			pb->RestoreData();

			double tmp = obsrate[i*nrep + rep] - predrate[i*nrep + rep];
			meanratediff += tmp;
			varratediff += tmp * tmp;

			tmp = obsprofile[i*nrep + rep] - predprofile[i*nrep + rep];
			meanprofilediff += tmp;
			varprofilediff += tmp * tmp;

			tmp = obsdouble[i*nrep + rep] - preddouble[i*nrep + rep];
			meandoublediff += tmp;
			vardoublediff += tmp * tmp;

			if (predrate[i*nrep + rep] > obsrate[i*nrep + rep])	{
				pvrate ++;
			}
			if (predprofile[i*nrep + rep] > obsprofile[i*nrep + rep])	{
				pvprofile ++;
			}
			if (preddouble[i*nrep + rep] > obsdouble[i*nrep + rep])	{
				pvdouble ++;
			}
		}
		pb->DeleteMultipleSub();
	}
	cerr << '\n';
	cerr.flush();

	ppmeandiversity /= size * nrep;
	ppvardiversity /= size * nrep;
	ppvardiversity -= ppmeandiversity * ppmeandiversity;
	double zdiversity = (ppmeandiversity - obsdiversity) / sqrt(ppvardiversity);

	meaninv /= size*nrep;
	varinv /= size*nrep;
	varinv -= meaninv*meaninv;
	pvinv /= size*nrep;

	meanratediff /= size*nrep;
	varratediff /= size*nrep;
	varratediff -= meanratediff * meanratediff;

	meanprofilediff /= size*nrep;
	varprofilediff /= size*nrep;
	varprofilediff -= meanprofilediff * meanprofilediff;

	meandoublediff /= size*nrep;
	vardoublediff /= size*nrep;
	vardoublediff -= meandoublediff * meandoublediff;

	meanobsrate /= size*nrep;
	varobsrate /= size*nrep;
	varobsrate -= meanobsrate * meanobsrate;

	meanpredrate /= size*nrep;
	varpredrate /= size*nrep;
	varpredrate -= meanpredrate * meanpredrate;

	pvrate /= size*nrep;

	meanobsprofile /= size*nrep;
	varobsprofile /= size*nrep;
	varobsprofile -= meanobsprofile * meanobsprofile;

	meanpredprofile /= size*nrep;
	varpredprofile /= size*nrep;
	varpredprofile -= meanpredprofile * meanpredprofile;

	pvprofile /= size*nrep;

	meanobsdouble /= size*nrep;
	varobsdouble /= size*nrep;
	varobsdouble -= meanobsdouble * meanobsdouble;

	meanpreddouble /= size*nrep;
	varpreddouble /= size*nrep;
	varpreddouble -= meanpreddouble * meanpreddouble;

	pvdouble /= size*nrep;

	ofstream sos((SampleName + ".pplrt").c_str());

	sos << "observed diversity  : " << obsdiversity << '\n';
	sos << "posterior predictive: " << ppmeandiversity << '\n';
	sos << "z-score             : " << zdiversity << '\n';
	sos << "pp value            : " << ((double) pvmeandiversity) / size / nrep << '\t' << "(" << pvmeandiversity << "/" << size * nrep << ")" << '\n';
	sos << '\n';

	sos << "observed prop inv     : " << obsinv << '\n';
	sos << "posterior predictive  : " << meaninv << " +/- " << sqrt(varinv) << '\n';
	sos << "z-score               : " << (obsinv - meaninv) / sqrt(varinv) << '\n';
	sos << "pvalue                : " << pvinv << '\n';
	sos << '\n';
	
	sos << "observed rate lrt      : " << meanobsrate << " +/- " << sqrt(varobsrate) << '\n';
	sos << "predictive rate lrt    : " << meanpredrate << " +/- " << sqrt(varpredrate) << '\n';
	sos << "z-score                : " << meanratediff / sqrt(varratediff) << '\n';
	sos << "pp value               : " << pvrate << '\t' << "(" << size * nrep << ")" << '\n';
	sos << '\n';

	sos << "observed profile lrt   : " << meanobsprofile << " +/- " << sqrt(varobsprofile) << '\n';
	sos << "predictive profile lrt : " << meanpredprofile << " +/- " << sqrt(varpredprofile) << '\n';
	sos << "z-score                : " << meanprofilediff / sqrt(varprofilediff) << '\n';
	sos << "pp value               : " << pvprofile << '\t' << "(" << size * nrep << ")" << '\n';
	sos << '\n';

	sos << "observed double bias   : " << meanobsdouble << " +/- " << sqrt(varobsdouble) << '\n';
	sos << "predictive double bias : " << meanpreddouble << " +/- " << sqrt(varpreddouble) << '\n';
	sos << "z-score                : " << meandoublediff / sqrt(vardoublediff) << '\n';
	sos << "pp value               : " << pvdouble << '\t' << "(" << size * nrep << ")" << '\n';
	sos << '\n';

	sos.close();

	cout << '\n';
	system (("cat " + SampleName + ".pplrt").c_str());

	MakeHisto(obsrate,nrep*size,SampleName + ".obsrate.histo",ncat);
	MakeHisto(predrate,nrep*size,SampleName + ".predrate.histo",ncat);

	MakeHisto(obsprofile,nrep*size,SampleName + ".obsprofile.histo",ncat);
	MakeHisto(predprofile,nrep*size,SampleName + ".predprofile.histo",ncat);

	MakeHisto(obsdouble,nrep*size,SampleName + ".obsdouble.histo",ncat);
	MakeHisto(preddouble,nrep*size,SampleName + ".preddouble.histo",ncat);

}

void Sample::DiscrepancyTest(int nrep, int ncat, int priormode, int priorratemode, int priorrootstate)	{

	mParam->UniSubOnly = Yes;
	int size = GetSize();

	double* obsrate = new double[nrep * size];
	double meanobsrate = 0;
	double varobsrate = 0;
	double* predrate = new double[nrep * size];
	double meanpredrate = 0;
	double varpredrate = 0;
	double pvrate = 0;

	double* obsprofile = new double[nrep * size];
	double meanobsprofile = 0;
	double varobsprofile = 0;
	double* predprofile = new double[nrep * size];
	double meanpredprofile = 0;
	double varpredprofile = 0;
	double pvprofile = 0;

	double* obsdouble = new double[nrep * size];
	double meanobsdouble = 0;
	double varobsdouble = 0;
	double* preddouble = new double[nrep * size];
	double meanpreddouble = 0;
	double varpreddouble = 0;
	double pvdouble = 0;

	double obsdiversity = mParam->MeanDiversity();
	double ppmeandiversity = 0;
	double ppvardiversity = 0;
	double pvmeandiversity = 0;

	ofstream xyos((SampleName + ".ppdiscxy").c_str());

	for (int i=0; i<size; i++)	{
		cerr << '.';
		cerr.flush();
		PhyloBayes* pb = GetNextPB();
		pb->Update();
		if (!i)	{
			pb->CreateTrueSub();
		}
		pb->countmultiplesub = 1;
		pb->CreateMultipleSub();
		
		for (int rep=0; rep<nrep; rep++)	{
			double mean = 0;
			double max = 0;
			pb->DrawSiteVariables(0,0);
			pb->ResampleSub();

			pb->RateDiscrepancy(mean,max);
			xyos << mean << '\t';
			obsrate[i*nrep + rep] = mean;
			meanobsrate += mean;
			varobsrate += mean * mean;
			
			pb->ProfileDiscrepancy(mean,max);
			xyos << mean << '\t';
			// pb->OptimiseSiteProfile(mean,max);
			obsprofile[i*nrep + rep] = mean;
			meanobsprofile += mean;
			varobsprofile += mean * mean;
			
			mean = pb->DoubleDiscrepancy();
			xyos << mean << '\t';
			obsdouble[i*nrep + rep] = mean;
			meanobsdouble += mean;
			varobsdouble += mean * mean;
			
			pb->SimulateData(priormode, priorratemode, priorrootstate,0);

			double diversity = mParam->MeanDiversity();
			ppmeandiversity += diversity;
			ppvardiversity += diversity*diversity;
			if (diversity < obsdiversity)	{
				pvmeandiversity ++;
			}

			pb->RateDiscrepancy(mean,max);
			xyos << mean << '\t';
			predrate[i*nrep + rep] = mean;
			meanpredrate += mean;
			varpredrate += mean * mean;
			
			pb->ProfileDiscrepancy(mean,max);
			xyos << mean << '\t';
			predprofile[i*nrep + rep] = mean;
			meanpredprofile += mean;
			varpredprofile += mean * mean;
			
			mean = pb->DoubleDiscrepancy();
			xyos << mean << '\n';
			preddouble[i*nrep + rep] = mean;
			meanpreddouble += mean;
			varpreddouble += mean * mean;

			pb->RestoreData();

			if (predrate[i*nrep + rep] > obsrate[i*nrep + rep])	{
				pvrate ++;
			}
			if (predprofile[i*nrep + rep] > obsprofile[i*nrep + rep])	{
				pvprofile ++;
			}
			if (preddouble[i*nrep + rep] > obsdouble[i*nrep + rep])	{
				pvdouble ++;
			}
		}
		
		pb->DeleteMultipleSub();
	}
	cerr << '\n';
	cerr.flush();

	ppmeandiversity /= size * nrep;
	ppvardiversity /= size * nrep;
	ppvardiversity -= ppmeandiversity * ppmeandiversity;
	double zdiversity = (ppmeandiversity - obsdiversity) / sqrt(ppvardiversity);

	meanobsrate /= size*nrep;
	varobsrate /= size*nrep;
	varobsrate -= meanobsrate * meanobsrate;

	meanpredrate /= size*nrep;
	varpredrate /= size*nrep;
	varpredrate -= meanpredrate * meanpredrate;

	pvrate /= size*nrep;

	meanobsprofile /= size*nrep;
	varobsprofile /= size*nrep;
	varobsprofile -= meanobsprofile * meanobsprofile;

	meanpredprofile /= size*nrep;
	varpredprofile /= size*nrep;
	varpredprofile -= meanpredprofile * meanpredprofile;

	pvprofile /= size*nrep;

	meanobsdouble /= size*nrep;
	varobsdouble /= size*nrep;
	varobsdouble -= meanobsdouble * meanobsdouble;

	meanpreddouble /= size*nrep;
	varpreddouble /= size*nrep;
	varpreddouble -= meanpreddouble * meanpreddouble;

	pvdouble /= size*nrep;

	ofstream sos((SampleName + ".ppdisc").c_str());

	sos << "observed diversity  : " << obsdiversity << '\n';
	sos << "posterior predictive: " << ppmeandiversity << '\n';
	sos << "z-score             : " << zdiversity << '\n';
	sos << "pp value            : " << ((double) pvmeandiversity) / size / nrep << '\t' << "(" << pvmeandiversity << "/" << size * nrep << ")" << '\n';
	sos << '\n';
	sos << "observed rate disc      : " << meanobsrate << " +/- " << sqrt(varobsrate) << '\n';
	sos << "predictive rate disc    : " << meanpredrate << " +/- " << sqrt(varpredrate) << '\n';
	sos << "pp value               : " << pvrate << '\t' << "(" << size * nrep << ")" << '\n';
	sos << '\n';
	sos << "observed profile dics   : " << meanobsprofile << " +/- " << sqrt(varobsprofile) << '\n';
	sos << "predictive profile dics : " << meanpredprofile << " +/- " << sqrt(varpredprofile) << '\n';
	sos << "pp value               : " << pvprofile << '\t' << "(" << size * nrep << ")" << '\n';
	sos << '\n';
	sos << "observed double disc   : " << meanobsdouble << " +/- " << sqrt(varobsdouble) << '\n';
	sos << "predictive double disc : " << meanpreddouble << " +/- " << sqrt(varpreddouble) << '\n';
	sos << "pp value               : " << pvdouble << '\t' << "(" << size * nrep << ")" << '\n';
	sos << '\n';
	sos.close();

	cout << '\n';
	system (("cat " + SampleName + ".ppdisc").c_str());

	/*
	MakeHisto(obsrate,nrep*size,SampleName + ".obsratedics.histo",ncat);
	MakeHisto(predrate,nrep*size,SampleName + ".predratedisc.histo",ncat);

	MakeHisto(obsprofile,nrep*size,SampleName + ".obsprofiledisc.histo",ncat);
	MakeHisto(predprofile,nrep*size,SampleName + ".predprofiledisc.histo",ncat);

	MakeHisto(obsdouble,nrep*size,SampleName + ".obsdoubledisc.histo",ncat);
	MakeHisto(preddouble,nrep*size,SampleName + ".preddoubledisc.histo",ncat);
	*/

}

// ---------------------------------------------------------------------------
//		 Homoplasy test / Stochastic Mapping version
// ---------------------------------------------------------------------------

void Sample::HomoplasyStochastic(int nrep, int ncat, int priormode, int priorratemode, int priorrootstate, int redraw, double cutoff)	{

	int size = GetSize();

	double* obs = new double[nrep * size];
	double meanobs = 0;
	double varobs = 0;
	
	double* pred = new double[nrep * size];
	double meanpred = 0;
	double varpred = 0;

	double pvsat = 0;

	double obsdiv = mParam->MeanDiversity();
	double meandiv = 0;
	double vardiv = 0;
	double pvdiv = 0;
	double obsinv = mParam->PropInvariant();
	double meaninv = 0;
	double varinv = 0;
	double pvinv = 0;

	double meanobsnsub = 0;
	double varobsnsub = 0;
	double meanprednsub = 0;
	double varprednsub = 0;
	double ppnsub = 0;

	double* meansitesat = new double[mParam->Nsite];
	double* tmpsitesub = new double[mParam->Nsite];
	double* sitediv = new double[mParam->Nsite];

	for (int i=0; i<mParam->Nsite; i++)	{
		meansitesat[i] = 0;
	}

	mParam->MeanDiversity(0,sitediv);

	ofstream xyos((SampleName + ".homoxy").c_str());
	for (int i=0; i<size; i++)	{
		cerr << '.';
		cerr.flush();
		PhyloBayes* pb = GetNextPB();
		pb->Update();
		if (!i)	{
			pb->CreateTrueSub();
		}

		for (int rep=0; rep<nrep; rep++)	{
			pb->ResampleSub();
			pb->UpdateTotalTrueSub();
			double tmpobs = pb->HomoplasyPerSite();
			double tmpobsnsub = pb->NTrueSubPerSite(tmpsitesub);
			for (int j=0; j<mParam->Nsite; j++)	{
				meansitesat[j] += tmpsitesub[j] - sitediv[j] + 1;
			}
			obs[i*nrep + rep] = tmpobs;
			meanobs += tmpobs;
			varobs += tmpobs * tmpobs;
			meanobsnsub += tmpobsnsub;
			varobsnsub += tmpobsnsub * tmpobsnsub;

			pb->SimulateData(priormode, priorratemode, priorrootstate,redraw);
			pb->UpdateTotalTrueSub();
			double diversity = mParam->MeanDiversity();
			meandiv += diversity;
			vardiv += diversity * diversity;
			if (diversity < obsdiv)	{
				pvdiv++;
			}
			double inv = mParam->PropInvariant();
			meaninv += inv;
			varinv += inv*inv;
			if (inv > obsinv)	{
				pvinv ++;
			}

			double tmppred = pb->HomoplasyPerSite();
			double tmpprednsub = pb->NTrueSubPerSite();
			pred[i*nrep + rep] = tmppred;
			meanpred += tmppred;
			varpred += tmppred * tmppred;
			meanprednsub += tmpprednsub;
			varprednsub += tmpprednsub * tmpprednsub;

			if (tmppred > tmpobs)	{
				pvsat++;
			}
			if (tmpprednsub > tmpobsnsub)	{
				ppnsub++;
			}
			pb->RestoreData();
			xyos << tmpobs << '\t' << tmppred << '\n';
		}
	}
	cerr << '\n';
	cerr << '\n';
	cerr.flush();

	ofstream satos((SampleName + ".sitehomoplasy").c_str());
	ofstream maskos((SampleName + ".sitehomoplasymask").c_str());
	satos << mParam->Nsite << '\n';
	maskos << mParam->Nsite << '\n';
	int totnsite = 0;
	for (int i=0; i<mParam->Nsite; i++)	{
		meansitesat[i] /= size;
		satos << meansitesat[i] << '\t';
		if ((sitediv[i] > 1) && (meansitesat[i] < cutoff))	{
			totnsite++;
			maskos << 1 << '\t';
		}
		else	{
			maskos << 0 << '\t';
		}
	}
	satos << '\n';
	maskos << '\n';
	cerr << "site specific saturation in " << SampleName << ".sitehomoplasy\n";
	cerr << "corresponding masks ( > " << cutoff << ") in " << SampleName << ".sitehomoplasymask\n";
	cerr << "number of sites included : " << totnsite  << '\n';
	cerr << '\n';

	meanobs /= size*nrep;
	varobs /= size*nrep;
	varobs -= meanobs * meanobs;

	meanpred /= size*nrep;
	varpred /= size*nrep;
	varpred -= meanpred * meanpred;

	pvsat /= size*nrep;

	meanobsnsub /= size*nrep;
	varobsnsub /= size*nrep;
	varobsnsub -= meanobsnsub * meanobsnsub;

	meanprednsub /= size*nrep;
	varprednsub /= size*nrep;
	varprednsub -= meanprednsub * meanprednsub;

	ppnsub /= size*nrep;

	meandiv /= size*nrep;
	vardiv /= size*nrep;
	vardiv -= meandiv*meandiv;
	pvdiv /= size*nrep;

	meaninv /= size*nrep;
	varinv /= size*nrep;
	varinv -= meaninv*meaninv;
	pvinv /= size*nrep;

	ofstream sos((SampleName + ".homoplasy").c_str());

	sos << "observed subs number  : " << meanobsnsub << " +/- " << sqrt(varobsnsub) << '\n';
	sos << "posterior predictive  : " << meanprednsub << " +/- " << sqrt(varprednsub) << '\n';
	sos << "pvalue                : " << ppnsub << '\n';
	sos << '\n';
	
	sos << "observed homoplasy    : " << meanobs << " +/- " << sqrt(varobs) << '\n';
	sos << "posterior predictive  : " << meanpred << " +/- " << sqrt(varpred) << '\n';
	sos << "pvalue                : " << pvsat << '\n';
	sos << '\n';
	
	sos << "observed diversity    : " << obsdiv << '\n';
	sos << "posterior predictive  : " << meandiv << " +/- " << sqrt(vardiv) << '\n';
	sos << "z-score               : " << (meandiv - obsdiv) / sqrt(vardiv) << '\n';
	sos << "pvalue                : " << pvdiv << '\n';
	sos << '\n';
	
	sos << "observed prop inv     : " << obsinv << '\n';
	sos << "posterior predictive  : " << meaninv << " +/- " << sqrt(varinv) << '\n';
	sos << "z-score               : " << (meaninv - obsinv) / sqrt(varinv) << '\n';
	sos << "pvalue                : " << pvinv << '\n';
	sos << '\n';
	
	sos.close();

	MakeHisto(obs,nrep*size,SampleName + ".obssat.histo",ncat);
	MakeHisto(pred,nrep*size,SampleName + ".predsat.histo",ncat);

	system(("cat " + SampleName + ".homoplasy").c_str());

	cout << '\n';
	cout << "histograms in " << SampleName << ".obssat.histo and " << SampleName << ".predsat.histo\n";
	cout << '\n';

	delete[] obs;
	delete[] pred;
}


// ---------------------------------------------------------------------------
//		 Gene effects
// ---------------------------------------------------------------------------

void Sample::GenePostPred(string partition, int nrep, int ncat, int priormode, int priorratemode, int priorrootstate, int redraw)	{

	int size = GetSize();

	ifstream is(partition.c_str());
	int Ngene;
	is >> Ngene;
	int* GeneSize = new int[Ngene];
	for (int gene=0; gene<Ngene; gene++)	{
		is >> GeneSize[gene];
	}

	int* GeneTotalSub = new int[Ngene];
	int** GeneBranchTotalSub = new int*[Ngene];
	for (int gene=0; gene<Ngene; gene++)	{
		GeneBranchTotalSub[gene] = new int[mParam->Nnode];
	}

	double meanobsgene = 0;
	double varobsgene = 0;
	double meanpredgene = 0;
	double varpredgene = 0;
	double ppgene = 0;

	double meanobsgenebl = 0;
	double varobsgenebl = 0;
	double meanpredgenebl = 0;
	double varpredgenebl = 0;
	double ppgenebl = 0;

	double meanobsmaxgenebl = 0;
	double varobsmaxgenebl = 0;
	double meanpredmaxgenebl = 0;
	double varpredmaxgenebl = 0;
	double ppmaxgenebl = 0;

	double ppbl[mParam->Nnode];
	for (int j=0; j<mParam->Nnode; j++)	{
		ppbl[j] = 0;
	}

	PhyloBayes* pb = currentPB;
	pb->CreateTrueSub();

	for (int i=0; i<size; i++)	{
		cerr << '.';
		cerr.flush();
		GetNextPB();
		pb->Update();

		for (int rep=0; rep<nrep; rep++)	{
			pb->ResampleSub();
			pb->UpdateTotalTrueSub();
			double tmpobsgene = 0;
			double tmpobsgenebl = 0;
			double tmpobsmaxgenebl = 0;
			double tmpobsppbl[mParam->Nnode];
			pb->GetGeneVar(tmpobsgene,tmpobsgenebl,tmpobsmaxgenebl,tmpobsppbl,Ngene,GeneSize,GeneTotalSub,GeneBranchTotalSub);
			meanobsgene += tmpobsgene;
			varobsgene += tmpobsgene * tmpobsgene;
			meanobsgenebl += tmpobsgenebl;
			varobsgenebl += tmpobsgenebl * tmpobsgenebl;
			meanobsmaxgenebl += tmpobsmaxgenebl;
			varobsmaxgenebl += tmpobsmaxgenebl * tmpobsmaxgenebl;
			
			pb->SimulateData(priormode, priorratemode, priorrootstate,redraw);
			pb->UpdateTotalTrueSub();
			double tmppredgene = 0;
			double tmppredgenebl = 0;
			double tmppredmaxgenebl = 0;
			double tmppredppbl[mParam->Nnode];
			pb->GetGeneVar(tmppredgene,tmppredgenebl,tmppredmaxgenebl,tmppredppbl,Ngene,GeneSize,GeneTotalSub,GeneBranchTotalSub);
			meanpredgene += tmppredgene;
			varpredgene += tmppredgene * tmppredgene;
			meanpredgenebl += tmppredgenebl;
			varpredgenebl += tmppredgenebl * tmppredgenebl;
			meanpredmaxgenebl += tmppredmaxgenebl;
			varpredmaxgenebl += tmppredmaxgenebl * tmppredmaxgenebl;
			if (tmppredgene > tmpobsgene)	{
				ppgene ++;
			}
			if (tmppredgenebl > tmpobsgenebl)	{
				ppgenebl++;
			}
			if (tmppredmaxgenebl > tmpobsmaxgenebl)	{
				ppmaxgenebl++;
			}
			for (int j=0; j<mParam->Nnode; j++)	{
				if (tmppredppbl[j] > tmpobsppbl[j])	{
					ppbl[j]++;
				}
			}
			pb->RestoreData();
		}
	}
	cerr << '\n';
	cerr << '\n';
	cerr.flush();

	meanobsgene /= size*nrep;
	varobsgene /= size*nrep;
	varobsgene -= meanobsgene * meanobsgene;

	meanobsgenebl /= size*nrep;
	varobsgenebl /= size*nrep;
	varobsgenebl -= meanobsgenebl * meanobsgenebl;

	meanobsmaxgenebl /= size*nrep;
	varobsmaxgenebl /= size*nrep;
	varobsmaxgenebl -= meanobsmaxgenebl * meanobsmaxgenebl;

	meanpredgene /= size*nrep;
	varpredgene /= size*nrep;
	varpredgene -= meanpredgene * meanpredgene;

	meanpredgenebl /= size*nrep;
	varpredgenebl /= size*nrep;
	varpredgenebl -= meanpredgenebl * meanpredgenebl;

	meanpredmaxgenebl /= size*nrep;
	varpredmaxgenebl /= size*nrep;
	varpredmaxgenebl -= meanpredmaxgenebl * meanpredmaxgenebl;

	ppgene /= size*nrep;
	ppgenebl /= size*nrep;
	ppmaxgenebl /= size*nrep;

	for (int j=0; j<mParam->Nnode; j++)	{
		ppbl[j] /= size*nrep;
	}

	ofstream sos((SampleName + ".genepp").c_str());

	sos << "observed among gene variance : " << meanobsgene << " +/- " << sqrt(varobsgene) << '\n';
	sos << "posterior predictive         : " << meanpredgene << " +/- " << sqrt(varpredgene) << '\n';
	sos << "pvalue                       : " << ppgene << '\n';
	sos << '\n';
	
	sos << "mean observed among gene bl variance : " << meanobsgenebl << " +/- " << sqrt(varobsgenebl) << '\n';
	sos << "posterior predictive                 : " << meanpredgenebl << " +/- " << sqrt(varpredgenebl) << '\n';
	sos << "pvalue                               : " << ppgenebl << '\n';
	sos << '\n';

	/*
	sos << "max observed among gene bl variance : " << meanobsmaxgenebl << " +/- " << sqrt(varobsmaxgenebl) << '\n';
	sos << "posterior predictive                : " << meanpredmaxgenebl << " +/- " << sqrt(varpredmaxgenebl) << '\n';
	sos << "pvalue                              : " << ppmaxgenebl << '\n';
	sos << '\n';
	*/

	if (mParam->FixTopo)	{
		for (int j=0; j<mParam->Nnode; j++)	{
			if (pb->blfree(j))	{
				sos << pb->GetLeftNodeName(j) << '\t' << pb->GetRightNodeName(j) << '\t' << ppbl[j] << '\n';
			}
		}
	}

	sos.close();

	system(("cat " + SampleName + ".genepp").c_str());
}


// ---------------------------------------------------------------------------
//		 Homoplasy test / Maximum Parsimony version
// ---------------------------------------------------------------------------

void Sample::Homoplasy(int nrep, int priormode, int priorratemode, int priorrootstate)	{

	int Nstate = mParam->Nstate;

	int size = GetSize();

	double* obs = new double[Nstate];
	int* tmp = new int[Nstate];
	double* pred = new double[Nstate];
	for (int k=0; k<Nstate; k++)	{
		pred[k] = 0;
	}

	double obsinv = mParam->PropInvariant();
	double meaninv = 0;
	double varinv = 0;
	double pvinv = 0;

	double meandiversity = mParam->MeanDiversity(tmp);
	double ppmeandiversity = 0;
	double ppvardiversity = 0;
	double pvmeandiversity = 0;
	double meanhomoplasy = mParam->GetCurrentState()->ParsimonyScore() / mParam->Nsite - meandiversity +1;
	double ppmeanhomoplasy = 0;
	double ppvarhomoplasy = 0;
	double pvmeanhomoplasy = 0;
	double meanpars = mParam->GetCurrentState()->ParsimonyScore() / mParam->Nsite;
	double ppmeanpars = 0;
	double ppvarpars = 0;
	double pvmeanpars = 0;

	for (int k=0; k<Nstate; k++)	{
		obs[k] = (double) tmp[k];
	}

	for (int i=0; i<size; i++)	{
		cerr << '.';
		cerr.flush();
		PhyloBayes* pb = GetNextPB();
		pb->Update();
		for (int rep=0; rep<nrep; rep++)	{
			pb->SimulateData(priormode, priorratemode, priorrootstate);
			double diversity = mParam->MeanDiversity(tmp);
			double pars = ((double) pb->ParsimonyScore()) / mParam->Nsite;
			double homoplasy = pars - diversity + 1;
			for (int k=0; k<Nstate; k++)	{
				pred[k] += tmp[k];
			}
			ppmeandiversity += diversity;
			ppvardiversity += diversity*diversity;
			if (diversity < meandiversity)	{
				pvmeandiversity ++;
			}
			ppmeanhomoplasy += homoplasy;
			ppvarhomoplasy += homoplasy * homoplasy;
			if (homoplasy > meanhomoplasy)	{
				pvmeanhomoplasy ++;
			}
			ppmeanpars += pars;
			ppvarpars += pars * pars;
			if (pars < meanpars)	{
				pvmeanpars ++;
			}

			double inv = mParam->PropInvariant();
			meaninv += inv;
			varinv += inv*inv;
			if (inv > obsinv)	{
				pvinv ++;
			}

			pb->RestoreData();
		}
	}
	cerr << '\n';
	cerr << '\n';
	cerr.flush();

	meaninv /= size*nrep;
	varinv /= size*nrep;
	varinv -= meaninv*meaninv;
	pvinv /= size*nrep;

	ppmeandiversity /= size * nrep;
	ppvardiversity /= size * nrep;
	ppvardiversity -= ppmeandiversity * ppmeandiversity;
	double zdiversity = (ppmeandiversity - meandiversity) / sqrt(ppvardiversity);

	ppmeanhomoplasy /= size * nrep;
	ppvarhomoplasy /= size * nrep;
	ppvarhomoplasy -= ppmeanhomoplasy * ppmeanhomoplasy;
	double zhomoplasy = (ppmeanhomoplasy - meanhomoplasy) / sqrt(ppvarhomoplasy);

	ppmeanpars /= size * nrep;
	ppvarpars /= size * nrep;
	ppvarpars -= ppmeanpars * ppmeanpars;
	double zpars = (ppmeanpars - meanpars) / sqrt(ppvarpars);

	for (int k=0; k<Nstate; k++)	{
		pred[k] /= size * nrep;
	}

	ofstream sos((SampleName + ".parshomoplasy").c_str());

	sos << "observed prop inv     : " << obsinv << '\n';
	sos << "posterior predictive  : " << meaninv << " +/- " << sqrt(varinv) << '\n';
	sos << "z-score               : " << (meaninv - obsinv) / sqrt(varinv) << '\n';
	sos << "pvalue                : " << pvinv << '\n';
	sos << '\n';
	
	sos << "observed parsimony  : " << meanpars << '\n';
	sos << "posterior predictive: " << ppmeanpars << '\n';
	sos << "z-score             : " << zpars << '\n';
	sos << "pp value            : " << ((double) pvmeanpars) / size / nrep << '\t' << "(" << pvmeanpars << "/" << size * nrep << ")" << '\n';
	sos << '\n';

	sos << "observed diversity  : " << meandiversity << '\n';
	sos << "posterior predictive: " << ppmeandiversity << '\n';
	sos << "z-score             : " << zdiversity << '\n';
	sos << "pp value            : " << ((double) pvmeandiversity) / size / nrep << '\t' << "(" << pvmeandiversity << "/" << size * nrep << ")" << '\n';
	sos << '\n';

	sos << "observed homoplasy  : " << meanhomoplasy << '\n';
	sos << "posterior predictive: " << ppmeanhomoplasy << '\n';
	sos << "z-score             : " << zhomoplasy << '\n';
	sos << "pp value            : " << ((double) pvmeanhomoplasy) / size / nrep << '\t' << "(" << pvmeanhomoplasy << "/" << size * nrep << ")" << '\n';
	sos << '\n';

	sos.close();

	ofstream dos((SampleName + ".divhisto").c_str());
	for (int k=0; k<Nstate; k++)	{
		dos << k << '\t' << obs[k] << '\t' << pred[k] << '\n';
	}

	string cat = "cat "  + SampleName + ".parshomoplasy";
	system(cat.c_str());

	delete[] tmp;
	delete[] obs;
	delete[] pred;

}


// ---------------------------------------------------------------------------
//		 Homogeneity Test
// ---------------------------------------------------------------------------

double Sample::HomogeneityMeasure(int** data, double* dev, double& mean)	{
	
	int Nstate = mParam->Nstate;
	int Ntaxa = mParam->Ntaxa;
	int Nsite = mParam->Nsite;

	double refp[Nstate];
	double p[Ntaxa][Nstate];

	for (int k=0; k<Nstate; k++)	{
		refp[k] = 0;
	}
	for (int j=0; j<Ntaxa; j++)	{
		for (int k=0; k<Nstate; k++)	{
			p[j][k] = 0;
		}
	}

	// observed
	for (int j=0; j<Ntaxa; j++)	{
		int count = 0;
		for (int i=0; i<Nsite; i++)	{
			if (data[j][i] != unknown)	{
				count++;
				p[j][data[j][i]] += 1;
			}
		}
		double total = 0;
		for (int k=0; k<Nstate; k++)	{
			p[j][k] /= count;
			total += p[j][k];
			refp[k] += p[j][k];
		}
		if (fabs(total-1) > 1e-3)	{
			cerr << "numerical error in homogeneity test\n";
			exit(1);
		}
	}
	double total = 0;
	for (int k=0; k<Nstate; k++)	{
		refp[k] /= Ntaxa;
		total += refp[k];
	}
	if (fabs(total-1) > 1e-3)	{
		cerr << "numerical error in homogeneity test\n";
		exit(1);
	}

	double max = 0;
	mean = 0;
	for (int j=0; j<Ntaxa; j++)	{
		dev[j] = 0;
		for (int k=0; k<Nstate; k++)	{
			dev[j] += (refp[k] - p[j][k]) * (refp[k] - p[j][k]);
		}
		if (max < dev[j])	{
			max = dev[j];
		}
		mean += dev[j];
	}
	mean /= Ntaxa;
	return max;
}

double Sample::ArginineSkew(int** data, double* dev)	{
	
	int Nstate = mParam->Nstate;
	int Ntaxa = mParam->Ntaxa;
	int Nsite = mParam->Nsite;

	if  (Nstate != 4)	{
		cerr << "error : arginine test only on nucleotide data\n";
		exit(1);
	}
	if (Nsite % 3)	{
		cerr << "error : number of positions should be multiple of 3 (codon alignment)\n";
		exit(1);
	}

	// observed
	double skew[Ntaxa];
	double meanskew = 0;
	for (int j=0; j<Ntaxa; j++)	{
		double cgn = 0;
		double agr = 0;
		for (int i=0; i<Nsite/3; i++)	{
			if ((data[j][3*i] == 1) && (data[j][3*i+1] == 2))	{
				cgn++;
			}
			else if ((data[j][3*i] == 0) && (data[j][3*i+1] == 2) && ((data[j][3*i+2] == 0) || (data[j][3*i+2] == 2)))	{
				agr++;
			}
		}
		skew[j] = (cgn - agr) / (cgn + agr);
		meanskew += skew[j];
	}		
	meanskew /= Ntaxa;

	double sum2 = 0;

	for (int j=0; j<Ntaxa; j++)	{
		double tmp = (skew[j]-meanskew) * (skew[j]-meanskew) / meanskew;
		// double tmp = (skew[j]-meanskew) * (skew[j]-meanskew) / meanskew;
		dev[j] = tmp;
		sum2 += tmp;
	}
	return sum2 / Ntaxa;
}

double Sample::LeucineSkew(int** data, double* dev)	{
	
	int Nstate = mParam->Nstate;
	int Ntaxa = mParam->Ntaxa;
	int Nsite = mParam->Nsite;

	if  (Nstate != 4)	{
		cerr << "error : leucine test only on nucleotide data\n";
		exit(1);
	}
	if (Nsite % 3)	{
		cerr << "error : number of positions should be multiple of 3 (codon alignment)\n";
		exit(1);
	}

	// observed
	double skew[Ntaxa];
	double meanskew = 0;
	for (int j=0; j<Ntaxa; j++)	{
		double ctn = 0;
		double ttr = 0;
		for (int i=0; i<Nsite/3; i++)	{
			if ((data[j][3*i] == 1) && (data[j][3*i+1] == 3))	{
				ctn++;
			}
			else if ((data[j][3*i] == 3) && (data[j][3*i+1] == 3) && ((data[j][3*i+2] == 0) || (data[j][3*i+2] == 2)))	{
				ttr++;
			}
		}
		skew[j] = (ctn - ttr) / (ctn + ttr);
		meanskew += skew[j];
	}		
	meanskew /= Ntaxa;

	double sum2 = 0;

	for (int j=0; j<Ntaxa; j++)	{
		double tmp = (skew[j]-meanskew) * (skew[j]-meanskew) / meanskew;
		// double tmp = (skew[j]-meanskew) * (skew[j]-meanskew) / meanskew;
		dev[j] = tmp;
		sum2 += tmp;
	}
	return sum2 / Ntaxa;
}

double Sample::SerineSkew(int** data, double* dev)	{
	
	int Nstate = mParam->Nstate;
	int Ntaxa = mParam->Ntaxa;
	int Nsite = mParam->Nsite;

	if  (Nstate != 4)	{
		cerr << "error : serine test only on nucleotide data\n";
		exit(1);
	}
	if (Nsite % 3)	{
		cerr << "error : number of positions should be multiple of 3 (codon alignment)\n";
		exit(1);
	}

	// observed
	double skew[Ntaxa];
	double meanskew = 0;
	for (int j=0; j<Ntaxa; j++)	{
		double tcn = 0;
		double agy = 0;
		for (int i=0; i<Nsite/3; i++)	{
			if ((data[j][3*i] == 3) && (data[j][3*i+1] == 1))	{
				tcn++;
			}
			else if ((data[j][3*i] == 0) && (data[j][3*i+1] == 2) && ((data[j][3*i+2] == 1) || (data[j][3*i+2] == 3)))	{
				agy++;
			}
		}
		skew[j] = (tcn - agy) / (tcn + agy);
		meanskew += skew[j];
	}		
	meanskew /= Ntaxa;

	double sum2 = 0;

	for (int j=0; j<Ntaxa; j++)	{
		double tmp = (skew[j]-meanskew) * (skew[j]-meanskew) / meanskew;
		// double tmp = (skew[j]-meanskew) * (skew[j]-meanskew) / meanskew;
		dev[j] = tmp;
		sum2 += tmp;
	}
	return sum2 / Ntaxa;
}

double Sample::CodonHomogeneityMeasure(int** data, double* dev, double& mean)	{
	
	int Nstate = mParam->Nstate;
	int Ntaxa = mParam->Ntaxa;
	int Nsite = mParam->Nsite;
	int Naa = 20;

	if  (Nstate != 4)	{
		cerr << "error : serine test only on nucleotide data\n";
		exit(1);
	}
	if (Nsite % 3)	{
		cerr << "error : number of positions should be multiple of 3 (codon alignment)\n";
		exit(1);
	}


	double refp[Naa];
	double p[Ntaxa][Naa];

	for (int k=0; k<Naa; k++)	{
		refp[k] = 0;
	}
	for (int j=0; j<Ntaxa; j++)	{
		for (int k=0; k<Naa; k++)	{
			p[j][k] = 0;
		}
	}

	for (int j=0; j<Ntaxa; j++)	{
		int count = 0;
		for (int i=0; i<Nsite/3; i++)	{
			int i1 = data[j][3*i];
			int i2 = data[j][3*i+1];
			int i3 = data[j][3*i+2];
			if ((i1 != unknown) && (i2 != unknown) && (i3 != unknown))	{
				count++;
				int k = 0;
				while ((k<64) && (codonpos[0][k] != i1) && (codonpos[1][k] != i2) && (codonpos[2][k] != i3))	{
					k++;
				}
				if (k == 64)	{
					cerr << "error in codon test\n";
					exit(1);
				}
				p[j][k] ++;
			}
		}
		double total = 0;
		for (int k=0; k<Naa; k++)	{
			p[j][k] /= count;
			total += p[j][k];
			refp[k] += p[j][k];
		}
		if (fabs(total-1) > 1e-3)	{
			cerr << "numerical error in homogeneity test\n";
			exit(1);
		}
	}
	double total = 0;
	for (int k=0; k<Naa; k++)	{
		refp[k] /= Ntaxa;
		total += refp[k];
	}
	if (fabs(total-1) > 1e-3)	{
		cerr << "numerical error in homogeneity test\n";
		exit(1);
	}

	double max = 0;
	mean = 0;
	for (int j=0; j<Ntaxa; j++)	{
		dev[j] = 0;
		for (int k=0; k<Naa; k++)	{
			dev[j] += (refp[k] - p[j][k]) * (refp[k] - p[j][k]);
		}
		if (max < dev[j])	{
			max = dev[j];
		}
		mean += dev[j];
	}
	mean /= Ntaxa;
	return max;
}

void Sample::HomogeneityTest(int nrep, double threshold, int priormode, int priorratemode, int priorrootstate, int skew)	{

	int size = GetSize();
	int Ntaxa = mParam->Ntaxa;

	double meanobs = 0;
	double* obs = new double[Ntaxa];
	double maxobs = 0;
	if (skew == 1)	{
		maxobs = SerineSkew(mParam->Data,obs);
	}
	else if (skew == 2)	{
		maxobs = LeucineSkew(mParam->Data,obs);
	}
	else if (skew == 3)	{
		maxobs = ArginineSkew(mParam->Data,obs);
	}
	else if (skew == 4)	{
		maxobs = CodonHomogeneityMeasure(mParam->Data,obs,meanobs);
	}
	else	{
		maxobs = HomogeneityMeasure(mParam->Data,obs,meanobs);	
	}
	double* pp = new double[Ntaxa];

	double* val = new double[size * nrep];
	int l= 0;
	
	double mean[Ntaxa];
	double var[Ntaxa];
	double pv[Ntaxa];
	double maxpv = 0;
	double maxmean = 0;
	double maxvar = 0;
	double meanpv = 0;
	double meanmean = 0;
	double meanvar = 0;
	for (int j=0; j<Ntaxa; j++)	{
		mean[j] = 0;
		var[j] = 0;
		pv[j] = 0;
	}

	for (int i=0; i<GetSize(); i++)	{

		cerr << ".";
		cerr.flush();

		PhyloBayes* pb = GetNextPB();
		pb->Update();
		if (pb->InitModeStat == CatFix)	{
			pb->UpdateSumMode(mParam->GetCurrentState());
			pb->ResampleModeAff(mParam->GetCurrentState());
		}
		for (int rep=0; rep<nrep; rep++)	{
			pb->SimulateData(priormode, priorratemode, priorrootstate);
			double meanpp = 0;
			double maxpp = 0;
			if (skew == 1)	{
				maxpp = SerineSkew(mParam->Data,pp);
			}
			else if (skew == 2)	{
				maxpp = LeucineSkew(mParam->Data,pp);
			}
			else if (skew == 3)	{
				maxpp = ArginineSkew(mParam->Data,pp);
			}
			else if (skew == 4)	{
				maxpp = CodonHomogeneityMeasure(mParam->Data, pp, meanpp);
			}
			else	{
				maxpp = HomogeneityMeasure(mParam->Data, pp, meanpp);
			}
			val[l++] = maxpp;
			if (maxobs <= maxpp)	{
				maxpv ++;
			}
			maxmean += maxpp;
			maxvar += maxpp * maxpp;
			if (meanobs <= meanpp)	{
				meanpv ++;
			}
			meanmean += meanpp;
			meanvar += meanpp * meanpp;
			for (int j=0; j<Ntaxa; j++)	{
				if (obs[j] <= pp[j])	{
					pv[j] ++;
				}
				mean[j] += pp[j];
				var[j] += pp[j] * pp[j];
			}
			pb->RestoreData();
		}
	}

	cerr << '\n' << '\n';
	maxmean /= size*nrep;
	maxvar /= size*nrep;
	maxvar -= maxmean * maxmean;
	maxpv /= size*nrep;
	
	meanmean /= size*nrep;
	meanvar /= size*nrep;
	meanvar -= meanmean * meanmean;
	meanpv /= size*nrep;
	
	for (int j=0; j<Ntaxa; j++)	{
		mean[j] /= size*nrep;
		var[j] /= size*nrep;
		var[j] -= mean[j] * mean[j];
		pv[j] /= size*nrep;
	}

	int maxsize = 0;
	for (int j=0; j<Ntaxa; j++)	{
		if (maxsize < ((int) mParam->SpeciesNames[j].length()))	{
			maxsize = (int) (mParam->SpeciesNames[j].length());
		}
	}

	// MakeHisto(val,nrep*size,SampleName + ".ppcomphisto",10);
	ofstream os((SampleName + ".ppht").c_str());
	os << "homogeneity test\n";
	os << '\n';
	os << "   taxon";
	for (int i=0; i<maxsize ; i++)	{
		os << ' ';
	}
	os << "p-value\tz-score\n";
	os << '\n';
	int* taxa = new int[mParam->Ntaxa];
	int kept = 0;
	for (int j=0; j<Ntaxa; j++)	{
		if (threshold == -1)	{
			if (pv[j] < 0.05)	{
				os << " * ";
			}
			else	{
				os << "   ";
			}
		}
		else	{
			double z = (obs[j] - mean[j]) / sqrt(var[j]);
			if (z > threshold)	{
				os << " * ";
				taxa[j] = 0;
			}
			else	{
				os << "   ";
				taxa[j] = 1;
				kept++;
			}
		}
		os << mParam->SpeciesNames[j];
		for (int i=0; i<maxsize + 5 -((int) (mParam->SpeciesNames[j].length())); i++)	{
		 	os << ' ';
		}	
		os << ((double) ((int) (1000 * pv[j])))/1000 << '\t' << ((double) ((int) (1000*(obs[j] - mean[j]) / sqrt(var[j]))))/1000 << '\n';
	}
	os << '\n';
	os << "global test:\n";	
	os << '\n';
	os << "max \n";
	os << "observed   : " << maxobs << '\n';
	os << "mean pred  : " << maxmean << '\n';
	os << "p-value    : " << maxpv << '\n';
      	os << "z-score    : " << (maxobs - maxmean) / sqrt(maxvar) << '\n';
	os << '\n';	

	if (! skew)	{
	os << "mean \n";
	os << "observed   : " << meanobs << '\n';
	os << "mean pred  : " << meanmean << '\n';
	os << "p-value    : " << meanpv << '\n';
      	os << "z-score    : " << (meanobs - meanmean) / sqrt(meanvar) << '\n';
	os << '\n';	
	}
	os.close();
	
	system(("cat " + SampleName + ".ppht").c_str());

	if (threshold != -1)	{
		ostringstream s;
		s << SampleName << "_z" << threshold << ".ali";
		mParam->WriteDataToFile(s.str(),taxa);
		cerr << mParam->Ntaxa - kept;
		if (mParam->Ntaxa -kept > 1)	{
			cerr << " taxa were removed\n";
		}
		else	{
			cerr << " taxon was removed\n";
		}
		cerr << "resulting dataset in " << s.str() << '\n';
		cerr << '\n';
	}	
	delete[] taxa;
	delete[] val;
}


// ---------------------------------------------------------------------------
//		 Diversity Test
// ---------------------------------------------------------------------------

void Sample::DifferentialDiversity(string taxonmask, int nrep, int ncat, int priormode, int priorratemode, int priorrootstate)	{

	int* mask = new int[mParam->Ntaxa];
	for (int j=0; j<mParam->Ntaxa; j++)	{
		mask[j] = 0;
	}
	ifstream is(taxonmask.c_str());
	int N;
	is >> N;
	for (int k=0; k<N; k++)	{
		string name;
		is >> name;
		int j = 0;
		while ((j<mParam->Ntaxa) && (mParam->SpeciesNames[j] != name))	{
			j++;
		}
		if (j == mParam->Ntaxa)	{
			cerr << "error in differential diversity test; does not recognise : " << name << '\n';
			exit(1);
		}
		mask[j] = 1;
		cerr << mParam->SpeciesNames[j] << '\n';
	}

	int size = GetSize();
	double meandiversity = mParam->MeanDifferentialDiversity(mask);
	double ppmeandiversity = 0;
	double ppvardiversity = 0;
	double pvmeandiversity = 0;

	for (int i=0; i<size; i++)	{
		cerr << '.';
		cerr.flush();
		PhyloBayes* pb = GetNextPB();
		pb->Update();
		for (int rep=0; rep<nrep; rep++)	{
			pb->SimulateData(priormode, priorratemode, priorrootstate);
			double diversity = mParam->MeanDifferentialDiversity(mask);
			ppmeandiversity += diversity;
			ppvardiversity += diversity*diversity;
			if (diversity < meandiversity)	{
				pvmeandiversity ++;
			}
			pb->RestoreData();
		}
	}

	cerr << '\n';
	cerr.flush();

	ppmeandiversity /= size * nrep;
	ppvardiversity /= size * nrep;
	ppvardiversity -= ppmeandiversity * ppmeandiversity;
	double z = (ppmeandiversity - meandiversity) / sqrt(ppvardiversity);

	ofstream sos((SampleName + ".diffdiversity").c_str());
	sos << "observed diff diversity  : " << meandiversity << '\n';
	sos << "posterior predictive: " << ppmeandiversity << '\n';
	sos << "z-score             : " << z << '\n';
	sos << "pp value            : " << ((double) pvmeandiversity) / size / nrep << '\t' << "(" << pvmeandiversity << "/" << size * nrep << ")" << '\n';
	sos << '\n';

	cout << "observed diff diversity  : " << meandiversity << '\n';
	cout << "posterior predictive: " << ppmeandiversity << '\n';
	cout << "z-score             : " << z << '\n';
	cout << "pp value            : " << ((double) pvmeandiversity) / size / nrep << '\t' << "(" << pvmeandiversity << "/" << size * nrep << ")" << '\n';
	cout << '\n';

}


void Sample::Diversity(int nrep, int ncat, int priormode, int priorratemode, int priorrootstate)	{

	int Nstate = mParam->Nstate;
	int size = GetSize();

	double* obs = new double[Nstate];
	int* tmp = new int[Nstate];
	double* pred = new double[Nstate];
	for (int k=0; k<Nstate; k++)	{
		pred[k] = 0;
	}
	double meandiversity = mParam->MeanDiversity(tmp);
	double ppmeandiversity = 0;
	double ppvardiversity = 0;
	double pvmeandiversity = 0;

	for (int k=0; k<Nstate; k++)	{
		obs[k] = (double) tmp[k];
	}

	for (int i=0; i<size; i++)	{
		cerr << '.';
		cerr.flush();
		PhyloBayes* pb = GetNextPB();
		pb->Update();
		for (int rep=0; rep<nrep; rep++)	{
			pb->SimulateData(priormode, priorratemode, priorrootstate);
			double diversity = mParam->MeanDiversity(tmp);
			for (int k=0; k<Nstate; k++)	{
				pred[k] += tmp[k];
			}
			ppmeandiversity += diversity;
			ppvardiversity += diversity*diversity;
			if (diversity < meandiversity)	{
				pvmeandiversity ++;
			}
			pb->RestoreData();
		}
	}
	cerr << '\n';
	cerr.flush();

	ppmeandiversity /= size * nrep;
	ppvardiversity /= size * nrep;
	ppvardiversity -= ppmeandiversity * ppmeandiversity;
	double z = (ppmeandiversity - meandiversity) / sqrt(ppvardiversity);

	for (int k=0; k<Nstate; k++)	{
		pred[k] /= size * nrep;
	}

	ofstream sos((SampleName + ".diversity").c_str());
	sos << "observed diversity  : " << meandiversity << '\n';
	sos << "posterior predictive: " << ppmeandiversity << " +/- " << sqrt(ppvardiversity) << '\n';
	sos << "z-score             : " << z << '\n';
	sos << "pp value            : " << ((double) pvmeandiversity) / size / nrep << '\t' << "(" << pvmeandiversity << "/" << size * nrep << ")" << '\n';
	sos << '\n';

	ofstream dos((SampleName + ".divhisto").c_str());
	for (int k=0; k<Nstate; k++)	{
		dos << k << '\t' << obs[k] << '\t' << pred[k] << '\n';
	}
	cout << "observed diversity  : " << meandiversity << '\n';
	cout << "posterior predictive: " << ppmeandiversity << " +/- " << sqrt(ppvardiversity) << '\n';
	cout << "z-score             : " << z << '\n';
	cout << "pp value            : " << ((double) pvmeandiversity) / size / nrep << '\t' << "(" << pvmeandiversity << "/" << size * nrep << ")" << '\n';
	cout << '\n';

	delete[] tmp;
	delete[] obs;
	delete[] pred;
}

// ---------------------------------------------------------------------------
//		 Reconstituting evolution of continuous characters
// ---------------------------------------------------------------------------

void Sample::ReadCont(double meanlogT, int ps)	{

	int size = GetSize();
	// double* inf95 = new double[mParam->Nnode];
	// double* sup95 = new double[mParam->Nnode];
	double meanT[mParam->Nnode];
	double varT[mParam->Nnode];
	double meanGC[mParam->Nnode];
	double varGC[mParam->Nnode];
	for (int j=0; j<mParam->Nnode; j++)	{
		meanT[j] = 0;
		varT[j] = 0;
		meanGC[j] = 0;
		varGC[j] = 0;
	}
	
	for (int i=0; i<size; i++)	{
		PhyloBayes* pb = GetNextPB();
		for (int j=0; j<mParam->Nnode; j++)	{
			double T = exp(pb->logNHStatDistorter[pb->NHcat[j]][mParam->Nstate] + meanlogT) - 273;
			if (meanlogT >= 10)	{
				T = pb->logNHStatDistorter[pb->NHcat[j]][mParam->Nstate] + meanlogT;
			}
			double GC = 100.0 / (1 + exp(-pb->logNHStatDistorter[pb->NHcat[j]][mParam->Nstate+1]));
			meanT[j] += T;
			varT[j] += T * T;
			meanGC[j] += GC;
			varGC[j] += GC * GC;
		}
	}
	
	for (int j=0; j<mParam->Nnode; j++)	{
		meanT[j] /= size;
		varT[j] /= size;
		varT[j] -= meanT[j] * meanT[j];
		if (fabs(varT[j]) < 1e-6)	{
			varT[j] = 0;
		}
		meanGC[j] /= size;
		varGC[j] /= size;
		varGC[j] -= meanGC[j] * meanGC[j];
		if (fabs(varGC[j]) < 1e-6)	{
			varGC[j] = 0;
		}
	}
	PhyloBayes* pb = mParam->GetCurrentState();

	ofstream os((SampleName + ".tgc").c_str());
	os << "#node\tname\tmeanT\tstderr\tmeanGC\tstderr\n";
	for (int j=0; j<mParam->Nnode; j++)	{
		os << j << '\t';
		if (pb->tree[j].isLeaf())	{
			os << mParam->SpeciesNames[j] << '\t';
		}
		else if (pb->tree[j].isRoot())	{
			os << "Root\t";
		}
		else	{
			os << "Internal\t";
		}

		os << meanT[j] << '\t' << sqrt(varT[j]) << '\t';
		os << meanGC[j] << '\t' << sqrt(varGC[j]) << '\n';
	}

	Tree* tree = new Tree(pb);
	tree->Dichotomise();
	tree->ChangeNames(mParam->SpeciesNames, mParam->SpeciesNames, mParam->Ntaxa);
	tree->mRoot->SortLeavesAlphabetical();

	ofstream label_os((SampleName + ".labels").c_str());
	tree->Phylip(label_os, 1, 0, 1, 1);
	label_os.close();
	if (ps)	{
		tree->ToPS(SampleName + ".T", 10 ,24,meanT, varT);
		tree->ToPS(SampleName + ".GC", 10 ,24,meanGC, varGC);
	}
	cerr << "labelled tree in : " << SampleName << ".labels\n";
	if (ps)	{
		cerr << "postscript in    : " << SampleName << ".T.ps and in " << SampleName << ".GC.ps\n";
	}
	cerr << '\n';
}	


void Sample::ReadCont2(int ps)	{

	int size = GetSize();
	double meanRate[mParam->Nnode];
	double varRate[mParam->Nnode];
	double meanCont[mParam->ContNchar][mParam->Nnode];
	double varCont[mParam->ContNchar][mParam->Nnode];
	for (int j=0; j<mParam->Nnode; j++)	{
		meanRate[j] = 0;
		varRate[j] = 0;
		for (int k=0; k<mParam->ContNchar; k++)	{
			meanCont[k][j] = 0;
			varCont[k][j] = 0;
		}
	}
	
	for (int i=0; i<size; i++)	{
		PhyloBayes* pb = GetNextPB();
		for (int j=0; j<mParam->Nnode; j++)	{
			for (int k=0; k<mParam->ContNchar; k++)	{
				double cont = pb->ContRho[k+1][j];
				// double cont = exp(pb->ContRho[k+1][j]);
				meanCont[k][j] += cont;
				varCont[k][j] += cont * cont;
			}
			double rate = exp(pb->Rho[j]);
			meanRate[j] += rate;
			varRate[j] += rate * rate;
		}
	}
	
	for (int j=0; j<mParam->Nnode; j++)	{
		meanRate[j] /= size;
		varRate[j] /= size;
		varRate[j] -= meanRate[j] * meanRate[j];
		if (fabs(varRate[j]) < 1e-6)	{
			varRate[j] = 0;
		}
		
		for (int k=0; k<mParam->ContNchar; k++)	{
			meanCont[k][j] /= size;
			varCont[k][j] /= size;
			varCont[k][j] -= meanCont[k][j] * meanCont[k][j];
			if (fabs(varCont[k][j]) < 1e-6)	{
				varCont[k][j] = 0;
			}
		}
	}
	PhyloBayes* pb = mParam->GetCurrentState();

	ofstream os((SampleName + ".longevity").c_str());
	os << "#node\tname\tmeanRate\tstderr\tmeanLongevityt\tstderr\n";
	for (int j=0; j<mParam->Nnode; j++)	{
		os << j << '\t';
		if (pb->tree[j].isLeaf())	{
			os << mParam->SpeciesNames[j] << '\t';
		}
		else if (pb->tree[j].isRoot())	{
			os << "Root\t";
		}
		else	{
			os << "Internal\t";
		}

		os << meanRate[j] << '\t' << sqrt(varRate[j]) << '\t';
		for (int k=0; k<mParam->ContNchar; k++)	{
			os << meanCont[k][j] << '\t' << sqrt(varCont[k][j]) << '\t';
		}
		os << '\n';
	}

	Tree* tree = new Tree(pb);
	tree->Dichotomise();
	tree->ChangeNames(mParam->SpeciesNames, mParam->SpeciesNames, mParam->Ntaxa);
	tree->mRoot->SortLeavesAlphabetical();

	ofstream label_os((SampleName + ".labels").c_str());
	tree->Phylip(label_os, 1, 0, 1, 1);
	label_os.close();
	if (ps)	{
		tree->ToPS(SampleName + ".rate", 10 ,24,meanRate, varRate);
		for (int k=0; k<mParam->ContNchar; k++)	{
			ostringstream s;
			s << ".cont" << k;
			tree->ToPS(SampleName + s.str(), 10 ,24,meanCont[k], varCont[k]);
		}
	}
	cerr << "labelled tree in : " << SampleName << ".labels\n";
	cerr << '\n';
}	

void Sample::ReadPopEff(int ps)	{

	int size = GetSize();
	double meanRate[mParam->Nnode];
	double varRate[mParam->Nnode];
	for (int j=0; j<mParam->Nnode; j++)	{
		meanRate[j] = 0;
		varRate[j] = 0;
	}
	
	for (int i=0; i<size; i++)	{
		PhyloBayes* pb = GetNextPB();
		for (int j=0; j<mParam->Nnode; j++)	{
			double rate = pb->PopEffSize[pb->NHcat[j]];
			meanRate[j] += rate;
			varRate[j] += rate * rate;
		}
	}
	
	for (int j=0; j<mParam->Nnode; j++)	{
		meanRate[j] /= size;
		varRate[j] /= size;
		varRate[j] -= meanRate[j] * meanRate[j];
		if (fabs(varRate[j]) < 1e-6)	{
			varRate[j] = 0;
		}
	}
	PhyloBayes* pb = mParam->GetCurrentState();

	ofstream os((SampleName + ".popeffsize").c_str());
	os << "#node\tname\tmeanRate\tstderr\n";
	for (int j=0; j<mParam->Nnode; j++)	{
		os << j << '\t';
		if (pb->tree[j].isLeaf())	{
			os << mParam->SpeciesNames[j] << '\t';
		}
		else if (pb->tree[j].isRoot())	{
			os << "Root\t";
		}
		else	{
			os << "Internal\t";
		}

		os << meanRate[j] << '\t' << sqrt(varRate[j]) << '\n';
	}

	Tree* tree = new Tree(pb);
	tree->Dichotomise();
	tree->ChangeNames(mParam->SpeciesNames, mParam->SpeciesNames, mParam->Ntaxa);
	tree->mRoot->SortLeavesAlphabetical();

	ofstream tree_os((SampleName + ".popeff.tre").c_str());
	tree->ChronoPhylip(tree_os, 1, 1, varRate, meanRate, true); // withleaf = true
	tree_os.close();

	ofstream label_os((SampleName + ".labels").c_str());
	tree->Phylip(label_os, 1, 0, 1, 1);
	label_os.close();
	if (ps)	{
		for (int j=0; j<mParam->Nnode; j++)	{
			meanRate[j] *= 100;
			varRate[j] *= 100;
		}
		tree->ToPS(SampleName + ".popeff", 10 ,24,meanRate, varRate);
		cerr << "tree in : " << SampleName << ".popeff.ps\n";
		cerr << '\n';
	}
}	

void Sample::ReadBranchFreqs()	{

	if (!mParam->NH)	{
		cerr << "error in Sample::ReadBranchFreqs: only for branch-heterogeneous models\n";
		exit(1);
	}
	int size = GetSize();

	double meanfreq[mParam->Nnode][mParam->Nstate];
	for (int j=0; j<mParam->Nnode; j++)	{
		for (int k=0; k<mParam->Nstate; k++)	{
			meanfreq[j][k] = 0;
		}
	}
	
	list<double>** freq = new list<double>*[mParam->Nnode];
	for (int j=0; j<mParam->Nnode; j++)	{
		freq[j] = new list<double>[mParam->Nstate];
	}

	double* tmpstat = new double[mParam->Nstate];
	for (int i=0; i<size; i++)	{

		PhyloBayes* pb = GetNextPB();
		pb->UpdateSiteNumber();
		pb->UpdateNH();

		for (int j=0; j<mParam->Nnode; j++)	{

			pb->GetMeanBranchStat(tmpstat,j);

			for (int k=0; k<mParam->Nstate; k++)	{
				meanfreq[j][k] += tmpstat[k];
				freq[j][k].push_front(tmpstat[k]);
			}
		}
	}
	delete[] tmpstat;

	for (int j=0; j<mParam->Nnode; j++)	{
		for (int k=0; k<mParam->Nstate; k++)	{
			meanfreq[j][k] /= size;
		}
	}

	PhyloBayes* pb = mParam->GetCurrentState();

	ofstream os((SampleName + ".branchfreqs").c_str());
	os << "#node\tstatus\ttax1\ttax2";
	for (int k=0; k<mParam->Nstate; k++)	{
		os << "\tpi_" << k;
		os << " mean : median (inf,sup)";
	}
	os << '\n';

	for (int j=0; j<mParam->Nnode; j++)	{
		os << j << '\t';
		if (pb->tree[j].isLeaf())	{
			os << "Leaf"; // mParam->SpeciesNames[j];
		}
		else if (pb->tree[j].isRoot())	{
			os << "Root";
		}
		else	{
			os << "Intr";
		}

		os << '\t' << mParam->GetCurrentState()->GetLeftNodeName(j) << '\t' << mParam->GetCurrentState()->GetRightNodeName(j);

		/*
		for (int k=0; k<mParam->Nstate; k++)	{
			os << '\t' << meanfreq[j][k];
		}
		*/

		double median[mParam->Nstate];
		double inf95[mParam->Nstate];
		double sup95[mParam->Nstate];

		double alpha = 0.05;

		int take = (int) (alpha * size / 2);
		int mediansize = size / 2;

		for (int k=0; k<mParam->Nstate; k++)	{
			freq[j][k].sort();
			list<double>::const_iterator i = freq[j][k].begin();
			for (int l=0; l<take; l++)	{
				i++;
			}
			inf95[k] = *i;
			for (int l=take; l<mediansize; l++)	{
				i++;
			}
			median[k] = *i;
			for (int l=mediansize; l<size - take; l++)	{
				i++;
			}
			sup95[k] = *i;
			os << '\t' << meanfreq[j][k] << " : " << median[k] << " (" << inf95[k] << "," << sup95[k] << ")";
		}
		os << '\n';
	}

	Tree* tree = new Tree(pb);
	tree->Dichotomise();
	tree->ChangeNames(mParam->SpeciesNames, mParam->SpeciesNames, mParam->Ntaxa);
	tree->mRoot->SortLeavesAlphabetical();

	ofstream label_os((SampleName + ".labels").c_str());
	tree->Phylip(label_os, 1, 0, 1, 1);
	label_os.close();
}	



void Sample::Dating(int ps, int verbose, double alpha)	{

	int size = GetSize();
	list<double>* date = new list<double>[mParam->Nnode];
	double* inf95 = new double[mParam->Nnode];
	double* sup95 = new double[mParam->Nnode];
	double* offinf95 = new double[mParam->Nnode];
	double* offsup95 = new double[mParam->Nnode];
	double mean[mParam->Nnode];
	double var[mParam->Nnode];
	double meanrate[mParam->Nnode];
	double varrate[mParam->Nnode];
	double meaninstantrate[mParam->Nnode];
	double varinstantrate[mParam->Nnode];
	for (int j=0; j<mParam->Nnode; j++)	{
		mean[j] = 0;
		var[j] = 0;
	}
	for (int j=0; j<mParam->Nnode; j++)	{
		meanrate[j] = 0;
		varrate[j] = 0;
	}
	for (int j=0; j<mParam->Nnode; j++)	{
		meaninstantrate[j] = 0;
		varinstantrate[j] = 0;
	}
	
	double meanl[mParam->Nnode];
	double varl[mParam->Nnode];
	for (int j=0; j<mParam->Nnode; j++)	{
		meanl[j] = 0;
		varl[j] = 0;
	}
	
	double meanchi = 0;
	double meanchi2 = 0;
	double varchi = 0;
	double varchi2 = 0;

	double meanrescaledsigma = 0;
	double varrescaledsigma = 0;
	double meanrescaledtheta = 0;
	double varrescaledtheta = 0;
	double meanrescalednu = 0;
	double varrescalednu = 0;


	double* calunderflow = 0;
	double* caloverflow = 0;
	if (mParam->NCalib)	{
		calunderflow = new double[mParam->NCalib];
		caloverflow = new double[mParam->NCalib];
		for (int cal=0; cal<mParam->NCalib; cal++)	{
			calunderflow[cal] = 0;
			caloverflow[cal] = 0;
		}
	}

	if (verbose)	{
		ofstream dd_os((SampleName + ".datedist").c_str());
	}
	double meanheight = 0;
	for (int i=0; i<size; i++)	{

		cerr << '.';

		PhyloBayes* pb = GetNextPB();

		meanchi += pb->Chi;
		varchi += pb->Chi * pb->Chi;
		meanchi2 += pb->Chi2;
		varchi2 += pb->Chi2 * pb->Chi2;

		pb->SetTimes();

		double tmp = pb->Sigma / pb->Scale;
		meanrescaledsigma += tmp;
		varrescaledsigma += tmp * tmp;

		tmp = pb->Theta / pb->Scale;
		meanrescaledtheta += tmp;
		varrescaledtheta += tmp * tmp;

		tmp = pb->Mu / pb->Scale;
		meanrescalednu += tmp;
		varrescalednu += tmp * tmp;


		// calibrations
		for (int cal=0; cal<mParam->NCalib; cal++)	{
			int index = mParam->CalibIndex[cal];
			double t = pb->Scale* pb->tree[index].height;
			double t_l = mParam->CalibLower[cal];
			double t_u = mParam->CalibUpper[cal];

			if (t_u != -1)	{
				if (t > t_u)	{
					caloverflow[cal] ++;
				}
			}
		
			if (t_l != -1)	{
				if (t < t_l)	{
					calunderflow[cal] ++;
				}
			}
		}

		if ((! mParam->NormalApprox) && (mParam->BLModelSwitch == 0))	{
			pb->SwapBL();
			pb->NormaliseLengths();
			pb->SwapBL();
		}
		
		for (int j=0; j<mParam->Nnode; j++)	{
			double temp = pb->tree[j].height * pb->Scale;
			double duration = temp;
			if (pb->tree[j].isRoot())	{
				duration = 0;
			}
			else	{
				duration -= ((AmphiNode*) (pb->tree[j].up))->height * pb->Scale;
			}
			date[j].push_front(temp);
			mean[j] += temp;
			var[j] += temp * temp;
			double temp2 = 0;
			if (mParam->BLModelSwitch == 0)	{
				temp2 = pb->RigidBL[j];
			}
			else {
				temp2 = pb->BL[j];
			}
			meanl[j] += temp2;
			varl[j] += temp2 * temp2;
			double rate = 0;
			if (duration < 0)	{
				rate = - temp2 / duration;
			}
			meanrate[j] += rate;
			varrate[j] += rate*rate;
			double instantrate = pb->Rho[j] / pb->Scale * pb->Mu;
			meaninstantrate[j] += instantrate;
			varinstantrate[j] += instantrate * instantrate;
		}
		meanheight += pb->Scale;
		if (verbose)	{
			for (int j=0; j<mParam->Nnode; j++)	{
				pb->tree[j].branchLength *= pb->Scale;
			}
			pb->SetBL();
			Tree tree(pb);
			tree.Dichotomise();
			tree.ChangeNames(mParam->SpeciesNames, mParam->SpeciesNames, mParam->Ntaxa);
			tree.mRoot->SortLeavesAlphabetical();

			ofstream dd_os((SampleName + ".datedist").c_str(),IOS_APPEND);
			// tree.Phylip(dd_os, 1, 0, 1, 0);
			tree.ChronoPhylip(dd_os,1,1,0,0);
			dd_os.close();
			for (int j=0; j<mParam->Nnode; j++)	{
				pb->tree[j].branchLength /= pb->Scale;
			}
		}

	}

	cerr << '\n';

	for (int cal=0; cal<mParam->NCalib; cal++)	{
		calunderflow[cal] /= size;
		caloverflow[cal] /= size;
	}

	meanchi /= size;
	varchi /= size;
	varchi -= meanchi * meanchi;
	if (varchi < 0)	{
		varchi = 0;
	}
	meanchi2 /= size;
	varchi2 /= size;
	varchi2 -= meanchi2 * meanchi2;
	if (varchi2 < 0)	{
		varchi2 = 0;
	}

	meanrescaledsigma /= size;
	varrescaledsigma /= size;
	varrescaledsigma -= meanrescaledsigma * meanrescaledsigma;
	if (varrescaledsigma < 0)	{
		varrescaledsigma = 0;
	}
	meanrescalednu /= size;
	varrescalednu /= size;
	varrescalednu -= meanrescalednu * meanrescalednu;
	if (varrescalednu < 0)	{
		varrescalednu = 0;
	}
	meanrescaledtheta /= size;
	varrescaledtheta /= size;
	varrescaledtheta -= meanrescaledtheta * meanrescaledtheta;
	if (varrescaledtheta < 0)	{
		varrescaledtheta = 0;
	}

	for (int j=0; j<mParam->Nnode; j++)	{
		mean[j] /= size;
		var[j] /= size;
		var[j] -= mean[j] * mean[j];
		meanl[j] /= size;
		varl[j] /= size;
		varl[j] -= meanl[j] * meanl[j];
		meanrate[j] /= size;
		varrate[j] /= size;
		varrate[j] -= meanrate[j] * meanrate[j];
		meaninstantrate[j] /= size;
		varinstantrate[j] /= size;
		varinstantrate[j] -= meaninstantrate[j] * meaninstantrate[j];
	}
	meanheight /= size;
	PhyloBayes* pb = mParam->GetCurrentState();

	int take = (int) (alpha * size / 2);
	// int take = size / 40;
	for (int j=0; j<mParam->Nnode; j++)	{
		date[j].sort();
		list<double>::const_iterator i = date[j].begin();
		for (int k=0; k<take; k++)	{
			i++;
		}
		inf95[j] = *i;
		for (int k=0; k<=size - 2*take; k++)	{
			i++;
		}
		sup95[j] = *i;
	}
		
	if (mParam->NCalib)	{
		ofstream os((SampleName + ".calib").c_str());
		os << "upper\tlower\t\t\%overflow\t\%underflow\n\n";
		for (int cal=0; cal<mParam->NCalib; cal++)	{
			os << mParam->CalibUpper[cal] << '\t' << mParam->CalibLower[cal] << '\t' << '\t';
			if (mParam->CalibUpper[cal] != -1)	{
				os << ((double) ((int) (caloverflow[cal] * 10000))) / 100 << '\t';
			}
			else	{
				os << '-' << '\t';
			}
			os << '\t';
			if (mParam->CalibLower[cal] != -1)	{
				os << ((double) ((int) (calunderflow[cal] * 10000))) / 100 << '\t';
			}
			else	{
				os << '-' << '\t';
			}
			os << '\n';
		}
	}

	ofstream os((SampleName + ".dates").c_str());
	os << "#node\t(meandate\tstderr\tinf" << ((int) (100 * (1-alpha))) << "\tsup" << ((int) (100 * (1-alpha))) << "\tinstant rate\tstderr\taverage rate\tstderr\n";
	for (int j=0; j<mParam->Nnode; j++)	{
	// for (int j=mParam->Ntaxa; j<mParam->Nnode; j++)	{
		if ((!mParam->NCalib) && (pb->tree[j].isRoot()))	{
			var[j] = 0;
		}
		if (pb->tree[j].isLeaf())	{
			os << mParam->SpeciesNames[j] << '\t' << 0 << '\t' << 0 << '\t' << 0 << '\t' << 0 << '\t';
		}
		else	{
			os << j << '\t' << mean[j] << '\t' << sqrt(var[j]) << '\t' << inf95[j] << '\t' << sup95[j] << '\t';
		}
		// os << j << '\t' << mean[j] << '\t' << sqrt(var[j]) << '\t' << meanl[j] << '\t' << sqrt(varl[j]) << '\n';
		os << meaninstantrate[j] << '\t' << sqrt(varinstantrate[j]) << '\t';
		os << meanrate[j] << '\t' << sqrt(varrate[j]) << '\n';
	}

	for (int j=0; j<mParam->Nnode; j++)	{
		pb->tree[j].height = meanheight - mean[j];
		// pb->tree[j].label = (int) (sqrt(var[j]));
		pb->RigidBL[j] = meanl[j];
		meanl[j] *= 100;
	}
	pb->root->computeLengths();
	// pb->root->renormLengths(0.001);
	pb->SetBL();
	Tree* tree = new Tree(pb);
	tree->Dichotomise();
	tree->ChangeNames(mParam->SpeciesNames, mParam->SpeciesNames, mParam->Ntaxa);
	tree->mRoot->SortLeavesAlphabetical();

	ofstream label_os((SampleName + ".labels").c_str());
	tree->Phylip(label_os, 1, 0, 1, 1);
	label_os.close();
	ofstream cons_os((SampleName + ".chronogram").c_str());
	tree->ChronoPhylip(cons_os, 1, 1, inf95, sup95);
	cons_os.close();

	for (int j=0; j<mParam->Nnode; j++)	{
		meanrate[j] *= 1000;
	}
	ofstream rate_os((SampleName + ".ratetree").c_str());
	tree->ChronoPhylip(rate_os, 1, 1, 0, meanrate,1);
	rate_os.close();
	
	for (int j=0; j<mParam->Nnode; j++)	{
		offinf95[j] = inf95[j] - mean[j];
		offsup95[j] = sup95[j] - mean[j];
	}

	if (ps)	{
		tree->ChronoPS(SampleName + ".chronogram", 15 ,34,offinf95, offsup95);
	}
	cons_os.close();
	if (ps)	{
		// tree->ToPS(SampleName + ".chronogram", 15 ,34,mean, var);
	}
	ofstream tree_os((SampleName + ".tre").c_str());
	pb->SwapBL();
	Tree* tree2 = new Tree(pb);
	tree2->Dichotomise();
	tree2->ChangeNames(mParam->SpeciesNames, mParam->SpeciesNames, mParam->Ntaxa);
	tree2->mRoot->SortLeavesAlphabetical();
	tree2->Phylip(tree_os, 1, 0, 1, 0);
	tree_os.close();
	if (ps)	{
		tree2->ToPS(SampleName + ".tre", 15, 34, 1, 0, 1, 0);
	}

	int rootlabel = pb->root->label;
	ofstream hyperos((SampleName + ".hyperparam").c_str());
	hyperos << "Age of the root                 : "  << mean[rootlabel] << " +/- " << sqrt(var[rootlabel]) << '\n';
	hyperos << '\n';
	hyperos << "Substitution rate per site (nu) : " << meanrescalednu << " +/- " << sqrt(varrescalednu) << '\n';
	cout << "Age of the root                 : "  << mean[rootlabel] << " +/- " << sqrt(var[rootlabel]) << '\n';
	cout << '\n';
	cout << "Substitution rate per site (nu) : " << meanrescalednu << " +/- " << sqrt(varrescalednu) << '\n';
	if ((mParam->ClockModel == KT) || (mParam->ClockModel == UndebiasedKT))	{
		hyperos << "sigma^2                         : " << meanrescaledsigma << " +/- " << sqrt(varrescaledsigma) << '\n';
		cout << "sigma^2                         : " << meanrescaledsigma << " +/- " << sqrt(varrescaledsigma) << '\n';
	}
	else if (mParam->ClockModel == CIRrigid)	{
		hyperos << "sigma^2                         : " << meanrescaledsigma << " +/- " << sqrt(varrescaledsigma) << '\n';
		hyperos << "theta                           : " << meanrescaledtheta << " +/- " << sqrt(varrescaledtheta) << '\n';
		cout << "sigma^2                         : " << meanrescaledsigma << " +/- " << sqrt(varrescaledsigma) << '\n';
		cout << "theta                           : " << meanrescaledtheta << " +/- " << sqrt(varrescaledtheta) << '\n';
	}
	hyperos << '\n';
	cout << '\n';


	if (mParam->TimePrior == BD)	{
		meanchi /= mParam->MeanScale;
		varchi /= mParam->MeanScale * mParam->MeanScale;
		meanchi2 /= mParam->MeanScale;
		varchi2 /= mParam->MeanScale * mParam->MeanScale;
	}

	if (mParam->TimePrior)	{
		if (mParam->TimePrior == Dirich)	{
			cout << "Dirichlet prior hyperparameter\n";
			cout << meanchi << " +/- " << sqrt(varchi) << '\n';
		}
		if (mParam->TimePrior == BD)	{
			cout << '\n';
			cout << "Birth death prior hyperparameters\n";
			cout << "p1 : " << meanchi << " +/- " << sqrt(varchi) << '\n';
			cout << "p2 : " << meanchi2 << " +/- " << sqrt(varchi2) << '\n';
			cout << '\n';

		}
		if (mParam->TimePrior == Dirich)	{
			hyperos << "Dirichlet prior hyperparameter\n";
			hyperos << meanchi << " +/- " << sqrt(varchi) << '\n';
		}
		if (mParam->TimePrior == BD)	{
			hyperos << "Birth death prior hyperparameters\n";
			hyperos << "p1 : " << meanchi << " +/- " << sqrt(varchi) << '\n';
			hyperos << "p2 : " << meanchi2 << " +/- " << sqrt(varchi2) << '\n';
		}
	}
			
	cerr << "posterior means  : " << SampleName << ".hyperparam\n";
	
	cerr << "labelled tree in : " << SampleName << ".labels\n";
	cerr << "chronogram in    : " << SampleName << ".chronogram\n";
	if (ps)	{
		cerr << "postscript in    : " << SampleName << ".chronogram.ps\n";
	}
	cerr << "efflength tree in: " << SampleName << ".tre\n";
	if (ps)	{
		cerr << "postscript in    : " << SampleName << ".tre.ps\n";
	}
	if (mParam->NCalib)	{
		cerr << "info about over- or underflow of calibrated nodes in " << SampleName << ".calib\n";
	}
	cerr << '\n';
	delete[] date;

}	


