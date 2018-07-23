
#include "phylo.h"

double** Stationary;

int main(int argc, char* argv[])	{

	// initialise random 
	// Random::Random();

	string SampleName;
	string ChainName;

	int burnin = -1;
	int every = 1;
	int until = -1;
	int nrep = 1;

	int ppred = 0;
	int priormode = 0;
	int priorratemode = 0;
	int priorrootstate = 1;

	int midpointrooting = 0;

	int imin = -1;
	int imax = -1;
	// read arguments

	string mask = "None";

	// int singlefreq = 0;

	int statepostprob = 1;

	int ws = 0;

	try	{

		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];

			if (s == "-nrep")	{
				i++;
				if (i == argc) throw(0);
				s = argv[i];
				if (! IsInt(s))	{
					throw(0);
				}
				nrep = atoi(argv[i]);
			}
			else if (s == "-cat")	{
				i++;
				s = argv[i];
				if (s == "prior")	{
					priormode = 1;
				}
				else if ((s == "post") || (s == "posterior"))	{
					priormode = 0;
				}
				else	{
					cerr << "error after -cat : should say either \"prior\", or \"post\"\n";
					cerr << argv[i] << '\n';
					exit(1);
				}
			}
			else if (s == "-rate")	{
				i++;
				s = argv[i];
				if (s == "prior")	{
					priorratemode = 1;
				}
				else if ((s == "post") || (s == "posterior"))	{
					priorratemode = 0;
				}
				else	{
					cerr << "error after -rate : should say either \"prior\", or \"post\"\n";
					cerr << argv[i] << '\n';
					exit(1);
				}
			}
			else if (s == "-root")	{
				i++;
				s = argv[i];
				if (s == "prior")	{
					priorrootstate = 1;
				}
				else if ((s == "post") || (s == "posterior"))	{
					priorrootstate = 0;
				}
				else	{
					cerr << "error after -root : should say either \"prior\", or \"post\"\n";
					exit(1);
				}
			}
			else if (s == "-midpointrooting")	{
				midpointrooting = 1;
			}
			else if (s == "-anc")	{
				statepostprob = 1;
			}
            /*
			else if (s == "-comp")	{
				singlefreq = 1;
			}
            */
			else if (s == "-mask")	{
				i++;
				mask = argv[i];
			}
			else if (s == "-sub")	{
				statepostprob = 0;
			}
			else if (s == "-ws")	{
				statepostprob = 0;
				ws = 1;
			}
			else if (s == "-ppred")	{
				ppred = 1;
			}
			else if (s == "-i")	{
				i++;
				imin = atoi(argv[i]);
				imin--;
				i++;
				imax = atoi(argv[i]);
			}
			else if ( (s == "-x") || (s == "-extract") )	{
				i++;
				if (i == argc) throw(0);
				s = argv[i];
				if (! IsInt(s))	{
					throw(0);
				}
				burnin = atoi(argv[i]);
				i++;
				if (i == argc) throw(0);
				s = argv[i];
				if (IsInt(s))	{
					every = atoi(argv[i]);
					i++;
					if (i == argc) throw(0);
					s = argv[i];
					if (IsInt(s))	{
						until = atoi(argv[i]);
					}
					else	{
						i--;
					}
				}
				else	{
					i--;
				}
			}
			else	{
				if (i != (argc -1))	{
					throw(0);
				}
				ChainName = argv[i];
			}
			i++;
		}
		if ((SampleName == "") && (ChainName == ""))	{
			throw(0);
		}
	}
	catch(...)	{
		cerr << "ancestral [-x <burnin> <every> <until>] [-sub -ppred -nrep <nrep> -i <imin> <imax>] <chainname>\n";
		exit(1);
	}

	if (SampleName == "")	{
		SampleName = ChainName + "_sample";
	}

	try	{

		Sample* sample = new Sample(ChainName,burnin,every,until);

		MCParameters* mParam = sample->GetParameters();
		if (midpointrooting)	{
			mParam->midpointrooting = 1;
		}
		int Nsite = mParam->Nsite;
		int Ntaxa = mParam->Ntaxa;

		if (! mParam->FixTopo)	{
			cerr << "error : ancestral works only under fixed topology\n";
			exit(1);
		}

		cout << '\n';
		cout << "Nsite : " << Nsite << '\n';
		cout << "Ntaxa : " << Ntaxa << '\n';
		cout << sample->GetSize() << " points to read\n";
		cout << '\n';
		cout.flush();

		int* sitemask = new int[mParam->Nsite];
		if (mask != "None")	{
			ifstream is(mask.c_str());
			int nsite;
			is >> nsite;
			if (nsite != mParam->Nsite)	{
				cerr << "error when reading mask: non matching number of sites : " << nsite << '\n';
				exit(1);
			}

			for (int i=0; i<nsite; i++)	{
				is >> sitemask[i];
			}
		}
		else	{
			for (int i=0; i<mParam->Nsite; i++)	{
				sitemask[i] = 0;
			}
			if (imin == -1)	{
				imin = 0;
			}
			if (imax == -1)	{
				imax = mParam->Nsite;
			}
			for (int i=imin; i<imax; i++)	{
				sitemask[i] = 1;
			}
		}

		/*
		double*** StatePostProb1 = new double**[mParam->Nnode];
		for (int i=0; i<mParam->Nnode; i++)	{
			StatePostProb1[i] = new double*[mParam->Nsite];
			for (int j=0; j<mParam->Nsite; j++)	{
				StatePostProb1[i][j] = new double[mParam->Nstate];
			}
		}

		double*** StatePostProb2 = new double**[mParam->Nnode];
		for (int i=0; i<mParam->Nnode; i++)	{
			StatePostProb2[i] = new double*[mParam->Nsite];
			for (int j=0; j<mParam->Nsite; j++)	{
				StatePostProb2[i][j] = new double[mParam->Nstate];
			}
		}

		if (mParam->NH)	{
			sample->SampleSub(StatePostProb2);
		}
		sample->ComputeStatePostProb(StatePostProb1);
		sample->Reset();
		sample->SampleSub(nrep,ppred,priormode,priorratemode,priorrootstate,sitemask,StatePostProb2);

		ofstream os((SampleName + ".comp").c_str());
		for (int i=mParam->Nstate; i<mParam->Nnode; i++)	{
			for (int j=0; j<mParam->Nsite; j++)	{
				for (int k=0; k<mParam->Nstate; k++)	{
					if (fabs(StatePostProb1[i][j][k] - StatePostProb2[i][j][k]) > 0.001)	{
						os << StatePostProb1[i][j][k] << '\t' << StatePostProb2[i][j][k] << '\t';
						os << i << '\t' << j << '\t' <<  k << '\t';
						os << mParam->GetCurrentState()->GetLeftNodeName(i) << '\t' << mParam->GetCurrentState()->GetRightNodeName(i) << '\n';
					}
				}
			}
		}
		os.close();
		exit(1);
		*/

		mParam->ZipSub = 1;

		if (statepostprob)	{
			sample->SampleSub(0,sitemask);
			/*
			if (mParam->NH)	{
				sample->SampleSub();
			}
			else	{
				sample->ComputeStatePostProb(0,singlefreq);
			}
			*/
			exit(1);
		}
		else if (ws)	{
			sample->SampleSubWS(priormode,priorratemode,priorrootstate);
			exit(1);
		}
		sample->SampleSub(nrep, ppred, priormode, priorratemode, priorrootstate, sitemask);

	}
	catch(...)	{
		exit(1);
	}
}

